        program ppmosaic_multitile

! Construct a mosaic of PPMAP images covering one or more Hi-GAL tiles.
!
! Input parameters:
!	ntiles  =	number of Hi-GAL tiles for which to merge the
!	                PPMAP output
!
!       fieldname1, fieldname2, ... fieldname<ntiles>
!               =       names of fields corresponding to the tiles,
!                       e.g. l136 l138 etc
!
        use iso_c_binding
        use omp_lib
        use m_matmul_omp

        implicit real (a-h, o-z)
        implicit integer (i-n)

        real,      parameter   :: pc      = 3.08568e18  ! cm
        real,      parameter   :: mu      = 2.8         ! mean molecular weight
        real,      parameter   :: mH      = 1.6726e-24  ! mass of H atom [g]
        real,      parameter   :: Msun    = 1.9891e33   ! solar mass [g]
        real,      parameter   :: fwhmobj = 2.		! FWHM of individual 
                                                        ! object [pixels]
        real(8),   parameter   :: dtor    = 0.0174533   ! degrees to radians
        real(8),   parameter   :: pi      = 3.1415927

        real(4),   allocatable :: a(:,:), buffer(:), outcube(:,:,:)
        real(4),   allocatable :: Tgrid(:), cctab(:,:), cct(:)
        real(4),   allocatable :: wavelengths(:)
        real(4),   allocatable :: cover(:,:),bigcover(:,:),rms(:)
        real(4),   allocatable :: cdens(:,:),temp(:,:)
        real(4),   allocatable :: tvar(:,:),tskew(:,:),tkurt(:,:)
        real(4),   allocatable :: rho(:,:,:),drho(:,:,:)
        real(4),   allocatable :: rh(:,:,:),rht(:,:),rhv(:,:),rhoff(:,:,:)
        real(4),   allocatable :: rhot(:,:),rhotx(:,:)
        real(4),   allocatable :: drhot(:,:),drhotc(:,:)
        real(4),   allocatable :: rhsum(:,:),uncscale(:,:),rchisq(:)
        real(4),   allocatable :: ukern(:,:),idelt(:),jdelt(:)
        real(4),   allocatable :: backcorr(:,:),nrat(:)
        real(4),   allocatable :: deltaset(:,:,:),dback(:,:)
        real(4),   allocatable :: diffset(:,:,:,:),diffsetok(:),diffs(:)
        real(4),   allocatable :: nplaneset(:,:,:,:),nplane(:,:),nplaneok(:,:)
        real(4),   allocatable :: wx(:),wy(:),w(:,:),twt(:,:),kern(:,:)
        real(4),   allocatable :: uncback(:),uncpeak(:),ell(:),bll(:)
        integer(4),allocatable :: x0set(:),y0set(:)
        integer(4),allocatable :: mfields(:),lookuptile(:)
        integer(4),allocatable :: nvalues(:,:),noverset(:,:),noset(:,:,:)
        integer(4),allocatable :: ifield(:),jfield(:),rmask(:,:),cmask(:,:)
        integer(4),allocatable :: noplaneset(:,:,:,:),noplane(:,:)
        integer(4),allocatable :: noplaneok(:,:)
        integer(4),allocatable :: nodiffset(:,:,:,:)
        integer(4),allocatable :: wt(:,:),wto(:,:),overlap(:,:)
        integer(4),allocatable :: ok(:,:,:,:),counter(:)
        character(len=10),        allocatable :: fieldnames(:)

        character(len=8)       :: date, ctype1, ctype2
        character(len=10)      :: time
        character(len=5)       :: zone
        character(len=40)      :: fieldname, mosaicname, sntiles
        character(len=80)      :: line,imagefile,inname,outfile,removeblanks
        character(len=80)      :: rhofile,sigrhofile,cdensfile,tempfile
        character(len=80)      :: tvarfile,tskewfile,tkurtfile
        integer, dimension (8) :: values
        integer, dimension (1) :: npk
        integer    status,cwin,dx0,dy0,x0diff,y0diff,edge
        real(4)    Nexp, Nprior, kappa300
        logical    incomplete, notfinished, bgcorr, verbose

        bgcorr = .true.         ! Apply background correction
        verbose = .false.       ! Minimal narration

        call date_and_time(date,time,zone,values)
        write(*,'("Begin PPMOSAIC_MULTITILE on ",a8," at ",a10)') date,time

! Read command line.
        m = iargc()
        if (m < 2) then
            print *, &
              'Syntax: ppmosaic_multitile ntiles tile1, tile2, ... tile<ntiles>'
            stop
        endif
        
        call getarg(1,sntiles)
        read(sntiles,*) ntiles

        allocate (fieldnames(ntiles))
        allocate (mfields(ntiles))
        allocate (idelt(ntiles))
        allocate (jdelt(ntiles))

        do n = 1,ntiles
            call getarg(n+1,fieldname)
            fieldnames(n) = fieldname
        enddo

! Construct output name for mosaic.
        mosaicname = fieldnames(1)
        if (ntiles > 1) then
            do n = 2,ntiles
                line = mosaicname//'-'//fieldnames(n)
                mosaicname = removeblanks(line)
            enddo
        endif
        print *,' '
        print *,'Mosaic name: ',mosaicname
        print *,' '

        rnan = -1.e35
        hnan = rnan/2.

! First, establish the new field size and centre coordinates. Will need
! some information from the PREMAP input parameters files.
        allocate (ell(ntiles))
        allocate (bll(ntiles))
        elmin = 1.e35
        elmax = -1.e35
        bmin = 1.e35
        bmax = -1.e35
        bmean = 0.

      do n = 1,ntiles
        fieldname = fieldnames(n)
        line = removeblanks(fieldname//'_premap.inp')
        write (*,'("Reading input parameters in ",a)') line
        open (unit=2, form='formatted', file=line,status='old', &
            action='read')
        getbands = 0
        nbands = 0
        Nt = 0
        notfinished = .true.

        do while (notfinished)
            read(2,'(a)',IOSTAT=iostatus) line
            if (iostatus /= 0) then 
                notfinished = .false.
            else
	        k = index(line,'<gloncent>')
	        if (k /= 0) read(line(1:k-1),*) glon0
	        k = index(line,'<glatcent>')
	        if (k /= 0) read(line(1:k-1),*) glat0
	        k = index(line,'<fieldsize>')
	        if (k /= 0) read(line(1:k-1),*) fieldwid,fieldhgt
	        k = index(line,'<pixel>')
	        if (k /= 0) read(line(1:k-1),*) pixel
            endif
        enddo
        close(2)
        if (n==1) pix1 = pixel
        if (pixel /= pix1) then
            print *,'Pixel size must be the same for all tiles'
            stop
        endif
        e = glon0 - 0.5*fieldwid/cos(glat0*dtor)
        if (e < elmin) elmin = e
        e = glon0 + 0.5*fieldwid/cos(glat0*dtor)
        if (e > elmax) elmax = e
        b = glat0 - 0.5*fieldhgt
        if (b < bmin) bmin = b
        b = glat0 + 0.5*fieldhgt
        if (b > bmax) bmax = b
        bmean = bmean + glat0
        ell(n) = glon0 + 0.5*fieldwid/cos(glat0*dtor)
        bll(n) = glat0 - 0.5*fieldwid
      enddo

        bmean = bmean/ntiles
        elcent = (elmin + elmax)/2.
        bcent = (bmin + bmax)/2.
        nxbig = nint((elmax - elmin)*cos(bmean*dtor)/(pixel/3600.))
        nybig = nint((bmax - bmin)/(pixel/3600.))
        ixbig = nxbig/2 + 1
        iybig = nybig/2 + 1
        allocate (bigcover(nxbig,nybig))
        bigcover = 0.
        print *,'Output mosaic size is',nxbig,'  x',nybig
        print *,'Centre coordinates:',elcent,bcent,' deg'
        print *,' '

      do n = 1,ntiles

! Read the PPMAP input parameters file.
        fieldname = fieldnames(n)
        line = removeblanks(fieldname//'_ppmap.inp')
        write (*,'("Reading input parameters in ",a)') line(1:40)
        open (unit=2, form='formatted', file=line,status='old', &
            action='read')
        read(2,*) nbands,nflds,Nt,icc,maxiterations,distance,kappa300,beta,eta
        mfields(n) = nflds
        read(2,'(a)') line
        if (n==1) allocate (Tgrid(Nt))
        read(2,*) Tgrid
        if (n==1) allocate (cctab(Nt,nbands))
        cctab = 1.

        if (icc /= 0) then      ! Colour corrections to be applied
            allocate (cct(Nt))
            do i = 1,nbands
                read(2,*) cct
                do it = 1,Nt
                    cctab(it,i) = cct(it)
                enddo
            enddo
            deallocate(cct)
        endif

        if (n==1) allocate (wavelengths(nbands))
        beammax = -1.e35

! Read in the observational images.
        do i = 1,nbands
	    read(2,'(a)') imagefile
	    if (i==1) then
                call readheader(imagefile,nximage,nyimage,cp1,cp2,cd2,status)
                nbuffer = nximage*nyimage
                allocate (buffer(nbuffer))
                allocate (a(nximage,nyimage))
            endif
            call readimage_wcs(imagefile,a,nximage,nyimage,buffer, &
                ctype1,ctype2,crpix1,crpix2,crval1,crval2,cdelt1,cdelt2, &
                crota2,pixel,wl,sig,status)
            if (status > 0) then
                print *,'Error reading image:'
                print *,imagefile
                stop
            endif
            wavelengths(i) = wl
            beam = 2.*cdelt2*3600.
            if (beam > beammax) beammax = beam
            if (i==1) then
                xref = crpix1 - 1.
                yref = crpix2 - 1.
            endif
        enddo
        close(2)
        deallocate(a)
        print *,mfields(n),' subfields present in '//fieldname

! Read coverage map.
        allocate (cover(nximage,nyimage))
        line = removeblanks(fieldname//'_coverage.fits')
        call readimage_basic(line,cover,nximage,nyimage,buffer,status)
        if (status > 0) then
            print *,'Error reading coverage map'
            stop
        endif
        ill = ixbig - nint((ell(n) - elcent)*cos(bmean*dtor)/(pixel/3600.))
        jll = iybig - nint((bcent - bll(n))/(pixel/3600.))
        print *,'tile',n,'   lower left pixel =',ill,jll
        do j = 1,nyimage
        do i = 1,nximage
             if (cover(i,j)==1.) bigcover(ill+i-1, jll+j-1) = cover(i,j)
        enddo
        enddo
        idelt(n) = ill - 1
        jdelt(n) = jll - 1
        deallocate(buffer)
        deallocate(cover)
      enddo

        deallocate(ell)
        deallocate(bll)

! Write out new coverage map.
        outfile = removeblanks(mosaicname//'_coverage.fits')
        print *,'Writing out '//outfile(1:40)
        call writeimage_wcs(outfile,bigcover,nxbig,nybig,ctype1,ctype2, &
            ixbig+1.,iybig+1.,elcent,bcent,-pixel/3600.,pixel/3600., &
            crota2,pixel,0.,0.,status)
        allocate (cover(nxbig,nybig))
        cover = bigcover
        deallocate(bigcover)

! Set up arrays for mosaicing.
        nfields = sum(mfields)
        allocate (lookuptile(nfields))
        n = 0
        do itile = 1,ntiles
            do m = 1,mfields(itile)
                n = n+1
                lookuptile(n) = itile
            enddo
        enddo
        nximage = nxbig
        nyimage = nybig
        allocate (cdens(nximage,nyimage))
        allocate (temp(nximage,nyimage))
        allocate (tvar(nximage,nyimage))
        allocate (tskew(nximage,nyimage))
        allocate (tkurt(nximage,nyimage))
        allocate (rho(nximage,nyimage,Nt))
        allocate (drho(nximage,nyimage,Nt))
        allocate (rchisq(nfields))
        allocate (nrat(nfields))
        allocate (x0set(nfields))
        allocate (y0set(nfields))

      do itile = 1,ntiles
        fieldname = fieldnames(itile)
        incomplete = .true.
        do while (incomplete)
          nfound = 0
          do n = 1,mfields(itile)
            write(inname,'(a40,"_",i5.5)') fieldname,n
            line = removeblanks(inname//'_rho.fits')
            call readheader3d(line,nx,ny,nz,status)
            if (status /= 0) then
                print *,' '
                print *,'Terminating due to missing PPMAP output file:'
                print *,line
                stop
            endif
            nfound = nfound + 1
          enddo
          if (nfound==mfields(itile)) incomplete = .false.
        enddo
        write(*,'(a,i5)') 'PPMAP output files all present for tile ',itile
      enddo

        allocate (rh(nx,ny,nz))
        allocate (rht(nx,ny))
        allocate (rhoff(nx,ny,Nt))
        allocate (rhot(nx,ny))
        allocate (rhotx(nx,ny))
        nbuffer = nx*ny*Nt
        allocate (buffer(nbuffer))

        do n = 1,nfields
            itile = lookuptile(n)
            fieldname = fieldnames(itile)
            noff = 0
            if (itile > 1) noff = mfields(itile-1)
            write(inname,'(a40,"_",i5.5)') fieldname,n-noff
            line = removeblanks(inname//'_rho.fits')
            call readrho(line,rh,nx,ny,Nt,buffer,ix0,iy0, &
                et,rhoobj,Nprior,Nexp,rchi,status)
            if (status==0) then
                x0set(n) = ix0 + idelt(itile)
                y0set(n) = iy0 + jdelt(itile)
                rchisq(n) = rchi
            else
                print *,'Could not read '//line
                stop
            endif
        enddo

        dx0 = 10000
        dy0 = 10000
        do n = 2,nfields
            x0diff = abs(x0set(n) - x0set(n-1))
            y0diff = abs(y0set(n) - y0set(n-1))
            if (x0diff > 0 .and. x0diff < dx0) dx0 = x0diff
            if (y0diff > 0 .and. y0diff < dy0) dy0 = y0diff
        enddo

        edge = nint(beammax/pixel)
        if (nfields > 1) edge = min(min(edge, (nx-dx0)/2-1), (ny-dy0)/2-1)
        do i = 1,edge
        do j = 1,nyimage
            cover(i,j) = 0.
            cover(nximage-edge+i,j) = 0.
        enddo
        enddo
        do j = 1,edge
        do i = 1,nximage
            cover(i,j) = 0.
            cover(i,nyimage-edge+j) = 0.
        enddo
        enddo

! Find background correction for each of the NFIELDS tiles.
        maxpasses = 1000
        cwin = 100
        tol = 0.01
        allocate (backcorr(nfields,Nt))
        backcorr = 0.

      if(bgcorr) then
        print *,'Begin edge matching ...'
        nxfields = (maxval(x0set) - minval(x0set))/dx0 + 1
        nyfields = (maxval(y0set) - minval(y0set))/dy0 + 1
        allocate (nvalues(nxfields,nyfields))
        allocate (noverset(nfields,Nt))
        allocate (deltaset(nx*ny,nfields,Nt))
        allocate (noset(nx*ny,nfields,Nt))
        allocate (ifield(nfields))
        allocate (jfield(nfields))
        allocate (dback(nfields,Nt))
        ifield = (x0set - minval(x0set))/dx0
        jfield = (y0set - minval(y0set))/dy0

        nvalues = 0
        do n = 1,nfields
            nvalues(ifield(n)+1,jfield(n)+1) = n
        enddo

        allocate (rms(maxpasses))
        iter = 0
        nx3 = 3*nx
        ny3 = 3*ny
        allocate (rhv(nx3,ny3))
        allocate (rmask(nx3,ny3))
        allocate (overlap(nx3,ny3))
        allocate (wt(nx3,ny3))
        allocate (cmask(nx3,ny3))
        allocate (nplaneset(nx3,ny3,3,3))
        allocate (noplaneset(nx3,ny3,3,3))
        allocate (diffset(nx3,ny3,3,3))
        allocate (nodiffset(nx3,ny3,3,3))
        allocate (nplane(nx3,ny3))
        allocate (noplane(nx3,ny3))
        allocate (nplaneok(nx3,ny3))
        allocate (noplaneok(nx3,ny3))
        allocate (ok(nx3,ny3,3,3))
        allocate (wto(nx,ny))
        ipass = -1
        notfinished = .true.

        do while (ipass < maxpasses-1 .and. notfinished)
	  ipass = ipass + 1
	  dback = 0.
	  if (ipass==0) then
            do n = 1,nfields
                itile = lookuptile(n)
                fieldname = fieldnames(itile)
                noff = 0
                if (itile > 1) noff = mfields(itile-1)
                write(inname,'(a40,"_",i5.5)') fieldname,n-noff
                line = removeblanks(inname//'_rho.fits')
                call readrho(line,rh,nx,ny,Nt,buffer, &
                    ix0,iy0,et,rhoobj,Nprior,Nexp,rchi,status)
		do it = 1,Nt
		    rht = rh(1:nx,1:ny,it) + backcorr(n,it)
		    rhv = rnan
                    do j = edge+1,ny-edge
                    do i = edge+1,nx-edge 
                        if (rh(i,j,it) /= 0.) rhv(nx+i,ny+j) = rht(i,j)
                    enddo
                    enddo
	  	    rmask = 0
                    where(rhv > hnan) rmask = 1

! Collect neighboring pixels.
		    nplaneset = rnan
		    noplaneset = 0
		    wt = 0
		    do joff = -1,1
                    do ioff = -1,1
		      ifo = ifield(n) + ioff
		      jfo = jfield(n) + joff
		      no = nvalues(ifo+1,jfo+1)
		      if (no /= 0 .and. ifo >= 0 .and. ifo < nxfields .and. &
                         jfo >= 0 .and. jfo < nyfields .and. (ioff /= 0 .or. &
                         joff /= 0)) then 
                        itile = lookuptile(no)
                        fieldname = fieldnames(itile)
                        noff = 0
                        if (itile > 1) noff = mfields(itile-1)
                        write(inname,'(a40,"_",i5.5)') fieldname,no-noff
                        line = removeblanks(inname//'_rho.fits')
                        call readrho(line,rhoff,nx,ny,Nt,buffer, &
                            ix0,iy0,et,rhoobj,Nprior,Nexp,rchi,status)
			rhot = rhoff(1:nx,1:ny,it) + backcorr(no,it)
			rhotx = rnan
                        do j = edge+1,ny-edge
                        do i = edge+1,nx-edge
			    if (rhoff(i,j,it) /= 0.) rhotx(i,j) = rhot(i,j)
                        enddo
                        enddo
			wto = 0
                        where(rhotx > hnan) wto = 1
			ill = nx + x0set(no) - x0set(n)
			jll = ny + y0set(no) - y0set(n)
                        do j = 1,ny
                        do i = 1,nx
			    nplaneset(ill+i,jll+j,ioff+2,joff+2) = rhot(i,j)
			    noplaneset(ill+i,jll+j,ioff+2,joff+2) = no
			    wt(ill+i,jll+j) = wt(ill+i,jll+j) + wto(i,j)
                        enddo
                        enddo
	    	      endif
		    enddo
                    enddo
		    cmask = 0
                    where(wt /= 0.) cmask = 1
                    overlap = 0
                    where(rmask*cmask /= 0.) overlap = 1
                    nover = sum(overlap)
		    if (nover /= 0) then
 	  		diffset = rnan
			nodiffset = 0
	  		do joff = -1,1
                        do ioff = -1,1
	    		    if ((ioff /= 0 .or. joff /= 0) .and. &
				    abs(ioff) /= abs(joff)) then
	      		      nplaneok = rnan
	      		      noplaneok = 0
	      		      nplane = nplaneset(1:nx3,1:ny3,ioff+2,joff+2)
	      		      noplane = noplaneset(1:nx3,1:ny3,ioff+2,joff+2)
                              where(overlap==1) nplaneok = nplane - rhv
			      where(overlap==1) noplaneok = noplane
                              do j = 1,ny3
                              do i = 1,nx3
	      		          diffset(i,j,ioff+2,joff+2) = nplaneok(i,j)
			          nodiffset(i,j,ioff+2,joff+2) = noplaneok(i,j)
                              enddo
                              enddo
	    		    endif
	  		enddo
                        enddo
                        ok = 0
                        where(diffset > hnan) ok = 1
                        nok = sum(ok)
			noverset(n,it) = nok
			if (nok /= 0) then
                          allocate (diffsetok(nok))
                          iok = 0
	  		  do joff = -1,1 
                          do ioff = -1,1 
	    		    if ((ioff /= 0 .or. joff /= 0) .and. &
			      abs(ioff) /= abs(joff)) then 
                             do j3 = 1,ny3
                             do i3 = 1,nx3
                               if (diffset(i3,j3,ioff+2,joff+2) > hnan) then
                                 iok = iok + 1
                                 diffsetok(iok) = diffset(i3,j3,ioff+2,joff+2)
                                 deltaset(iok,n,it) = diffsetok(iok)
			         noset(iok,n,it)=nodiffset(i3,j3,ioff+2,joff+2)
                               endif
                             enddo
                             enddo
                            endif
                          enddo
                          enddo
                            
			  call trimavg(diffsetok,nok,dmean,dsig,2.)
	  		  dback(n,it) = dmean
                          deallocate(diffsetok)
			endif
		    else 
                        if (it==1 .and. verbose) &
                            print *,'Field',n,' has no overlap'
                    endif
      		enddo
	    enddo
	    rmsprev = sqrt(sum(dback**2))
	  else 
	    do n = 1,nfields
            do it = 1,Nt
		nover = noverset(n,it)
		if (nover /= 0) then
                    allocate (diffs(nover))
		    do i = 1,nover 
			no = noset(i,n,it)
			diffs(i) = deltaset(i,n,it) + &
			    backcorr(no,it) - backcorr(n,it)
		    enddo
		    call trimavg(diffs,nover,dmean,dsig,2.)
	  	    dback(n,it) = dmean
                    deallocate(diffs)
		endif
	    enddo
	    enddo
	  endif
          backcorr = backcorr + dback/2.
	  rms(ipass+1) = sqrt(sum(dback**2))
	  if (verbose) print *,'ipass,rms:',ipass,rms(ipass+1)
	  if (mod((ipass+1),cwin)==0) then 
	    rmsmean = sum(rms(ipass-cwin+2:ipass+1))/cwin
	    conv = (rmsprev - rmsmean)/rmsmean
	    rmsprev = rmsmean
	    iter = iter + 1
	    write(*,'("Iteration",i3,":   RMS background error =",f12.4,' &
                //'";   conv =",f8.2)') iter,rmsmean,conv
	    if (conv < tol) notfinished = .false.
	  endif
        enddo

        deallocate(rhv)
        deallocate(rmask)
        deallocate(overlap)
        deallocate(cmask)
        deallocate(nplaneset)
        deallocate(noplaneset)
        deallocate(diffset)
        deallocate(nodiffset)
        deallocate(nplane)
        deallocate(noplane)
        deallocate(nplaneok)
        deallocate(noplaneok)
        deallocate(ok)
        deallocate(nvalues)
        deallocate(noverset)
        deallocate(deltaset)
        deallocate(noset)
        deallocate(ifield)
        deallocate(jfield)
        deallocate(dback)
        deallocate(wto)
        deallocate(wt)

! Adjust background correction to zero median at all temperatures.
        do it = 1,Nt 
            call nmedian(backcorr(1:nfields,it),nfields,bcmed)
            do n = 1,nfields
                backcorr(n,it) = backcorr(n,it) - bcmed
            enddo
        enddo
      endif

        deallocate(x0set)
        deallocate(y0set)
        deallocate(rhot)
        deallocate(rht)
        deallocate(rhoff)
        deallocate(rhotx)

! Calculate mosaic. Weight the individual tiles with a 2D Hann window,
! obtained from the product of two orthogonal 1D Hann windows.
        print *,'Calculating mosaic'
        allocate (twt(nximage,nyimage))
        allocate (wx(nx))
        allocate (wy(ny))
        allocate (w(nx,ny))
        wx = 0.
        wy = 0.
        nh = nx - 2*edge
        do i = 1,nh
            wx(edge+i) = (sin(pi*(i-1)/(nh-1.)))**2   ! 1D Hann window (x axis)
        enddo
        nh = ny - 2*edge
        do j = 1,nh
            wy(edge+j) = (sin(pi*(j-1)/(nh-1.)))**2   ! 1D Hann window (y axis)
        enddo
        do i = 1,nx
        do j = 1,ny
            w(i,j) = wx(i)*wy(j)                      ! 2D version
        enddo
        enddo
        deallocate(wx)
        deallocate(wy)
        rho = 0.
        twt = 0.

        do n = 1,nfields
            itile = lookuptile(n)
            fieldname = fieldnames(itile)
            noff = 0
            if (itile > 1) noff = mfields(itile-1)
            write(inname,'(a40,"_",i5.5)') fieldname,n-noff
            line = removeblanks(inname//'_rho.fits')
            call readrho(line,rh,nx,ny,Nt,buffer,ix0,iy0, &
                et,rhoobj,Nprior,Nexp,rchi,status)
            ix0 = ix0 + idelt(itile)
            iy0 = iy0 + jdelt(itile)
	    nrat(n) = Nexp/Nprior

	    if(ix0>=0 .and. iy0>=0 .and. ix0+nx<nximage .and. iy0+ny<nyimage) &
              then
                do j = iy0+1,iy0+ny
                do i = ix0+1,ix0+nx
	            do it = 1,Nt
                        rho(i,j,it) = rho(i,j,it) + w(i-ix0,j-iy0) *  &
                        (rh(i-ix0,j-iy0,it) + backcorr(n,it))
                    enddo
	            twt(i,j) = twt(i,j) + w(i-ix0,j-iy0)
                enddo
                enddo
	    endif
        enddo

        deallocate(rh)
        deallocate(buffer)
        deallocate(w)
        deallocate(backcorr)
        deallocate(fieldnames)
        deallocate(mfields)
        deallocate(idelt)
        deallocate(jdelt)
        deallocate(lookuptile)
        if (verbose) print *,'Mosaic completed'

! Construct "single object" kernel.
        ik = nint(1.5*fwhmobj)
        nk = 2*ik+1
        allocate (kern(nk,nk))
        aln2 = log(2.)
        do j = 1,nk
        do i = 1,nk
            kern(i,j) = exp(-4.*(aln2/fwhmobj**2)*((i-ik-1)**2 + (j-ik-1)**2))
        enddo
        enddo
        kern = kern/sum(kern)

! Calculate image cube and convolve with single-object profile.
        allocate(rhot(nximage,nyimage))
        allocate(rhotx(nximage,nyimage))
        do it = 1,Nt
            do j = 1,nyimage
            do i = 1,nximage
                if (twt(i,j) /= 0.) rho(i,j,it) = rho(i,j,it)/twt(i,j)
            enddo
            enddo
            rhot = rho(1:nximage,1:nyimage,it)
            where(rhot < 0.) rhot = 0.
            call convol(rhot,kern,rhotx,nximage,nyimage,nk)
            do j = 1,nyimage
            do i = 1,nximage
                rho(i,j,it) = rhotx(i,j)
            enddo
            enddo
        enddo
        deallocate(rhot)
        deallocate(rhotx)

! Calculate uncertainties.
        print *,'Calculating uncertainties'
        npk = maxloc(nrat)
        npeak = npk(1)
        allocate (uncback(Nt))
        unback = 0.
        allocate (uncpeak(Nt))
        uncpeak = 0.
        write(inname,'(a40,"_",i5.5)') fieldname,npeak
        line = removeblanks(inname//'_rho.fits')
        call pperr(line,wavelengths,nbands,Tgrid,Nt,cctab,pixel, &
            kappa300,beta,uncback,uncpeak)
        deallocate(wavelengths)
        deallocate(cctab)
        allocate (rhsum(nximage,nyimage))
        do j = 1,nyimage
        do i = 1,nximage
            rsum = 0.
            do it = 1,Nt
                rsum = rsum + rho(i,j,it)
            enddo
            rhsum(i,j) = rsum
        enddo
        enddo
        allocate (ukern(nx,ny))
        ukern = 1.
        allocate (uncscale(nximage,nyimage))
        call convol(rhsum,ukern,uncscale,nximage,nyimage,nx)
        ihb = nx/2
        do i = 1,ihb
          do j = 1,nyimage
            uncscale(i,j) = uncscale(ihb+1,j)
            uncscale(nximage-i+1,j) = uncscale(nximage-ihb,j)
          enddo
          do j = 1,ihb
            uncscale(i,j) = uncscale(ihb+1,ihb+1)
            uncscale(i,nyimage-j+1) = uncscale(ihb+1,nyimage-ihb)
          enddo
        enddo
        do j = 1,ihb
          do i = 1,nximage
            uncscale(i,j) = uncscale(i,ihb+1)
            uncscale(i,nyimage-j+1) = uncscale(i,nyimage-ihb)
          enddo
          do i = 1,ihb
            uncscale(nximage-i+1,j) = uncscale(nximage-ihb,ihb+1)
            uncscale(nximage-i+1,nyimage-j+1)=uncscale(nximage-ihb,nyimage-ihb)
          enddo
        enddo

        uncm = maxval(uncscale)
        uncscale = sqrt(uncscale/uncm)
        do it = 1,Nt
            du = max(uncpeak(it)-uncback(it), 0.)
            do j = 1,nyimage
            do i = 1,nximage
	        drho(i,j,it) = uncback(it) + uncscale(i,j)*du
            enddo
            enddo
	    write(*,'("T =",f6.2," K:  bkgnd & peak uncs:",' &
                //'f10.4,",",f12.4,"  x 10^20 cm^-2")') Tgrid(it),uncback(it), &
                uncpeak(it)
        enddo
        deallocate(uncpeak)
        deallocate(uncscale)
        deallocate(ukern)

! Calculate integrated column density and mean temperature.
        cdens = 0.
        temp = 0.
        do it = 1,Nt
	    cdens = cdens + rho(1:nximage,1:nyimage,it)
	    temp = temp + rho(1:nximage,1:nyimage,it)*Tgrid(it)
        enddo
        where(cdens > 0.) temp = temp/cdens

! Calculate temperature variance, skew, and kurtosis.
        tvar = 0.
        tskew = 0.
        tkurt = 0.
        do it = 1,Nt
	    tvar = tvar + rho(1:nximage,1:nyimage,it)*(Tgrid(it)-temp)**2
	    tskew = tskew + rho(1:nximage,1:nyimage,it)*(Tgrid(it)-temp)**3
	    tkurt = tkurt + rho(1:nximage,1:nyimage,it)*(Tgrid(it)-temp)**4
        enddo
        where(cdens > 0.) tvar = tvar/cdens
        where(cdens > 0.) tskew = tskew/cdens
        where(cdens > 0.) tkurt = tkurt/cdens

        call nmedian(rchisq,nfields,rchisqmed)
        print *,'Median reduced chi squared   =',rchisqmed
        call nmedian(nrat,nfields,ratmed)
        deallocate(nrat)
        if (verbose) print *,'Median ratio of Npost/Nprior =',ratmed
        write(*,'(" Total mass                   =",1pe12.3," Msun")') &
            sum(rho)*(1.e20/Msun)*mu*mH*(pixel*(dtor/3600.)*distance*pc)**2
        rhofile = removeblanks(mosaicname//'_tdenscube.fits')
        sigrhofile = removeblanks(mosaicname//'_sigtdenscube.fits')
        cdensfile = removeblanks(mosaicname//'_cdens.fits')
        tempfile = removeblanks(mosaicname//'_temp.fits')
        tvarfile = removeblanks(mosaicname//'_tvar.fits')
        tskewfile = removeblanks(mosaicname//'_tskew.fits')
        tkurtfile = removeblanks(mosaicname//'_tkurt.fits')
        print *,'Writing out: ',cdensfile(1:60)
        print *,'             ',tempfile(1:60)
        print *,'             ',tvarfile(1:60)
        print *,'             ',tskewfile(1:60)
        print *,'             ',tkurtfile(1:60)
        print *,'             ',rhofile(1:60)
        print *,'             ',sigrhofile(1:60)

! Write out image cube of differential column density.
        allocate (rhot(nximage,nyimage))
        allocate (drhot(nximage,nyimage))
        do it = 1,Nt
	    rhot = rho(1:nximage,1:nyimage,it)
            where(twt==0. .or. cover==0.) rhot = 0.
	    drhot = drho(1:nximage,1:nyimage,it)
            where(twt==0. .or. cover==0.) drhot = 0.
            do j = 1,nyimage
            do i = 1,nximage
	        rho(i,j,it) = rhot(i,j)
	        drho(i,j,it) = drhot(i,j)
            enddo
            enddo
        enddo
        deallocate(rhot)
        deallocate(drhot)

        x0 = ixbig + 1.
        y0 = iybig + 1.
        glon0 = elcent
        glat0 = bcent
        pixd = pixel/3600.
        call writerho(rhofile,rho,nximage,nyimage,Nt,ctype1,ctype2, &
            x0,y0,glon0,glat0,-pixd,pixd,crota2,1,1,Tgrid(1),Tgrid(Nt),eta, &
            rhoobj,Nprior,ratmed*Nprior,1.,maxiterations,rchisqmed,status)
        call writerho(sigrhofile,drho,nximage,nyimage,Nt,ctype1,ctype2, &
            x0,y0,glon0,glat0,-pixd,pixd,crota2,1,1,Tgrid(1),Tgrid(Nt),eta, &
            rhoobj,Nprior,ratmed*Nprior,1.,maxiterations,rchisqmed,status)
        deallocate(Tgrid)
        deallocate(rho)
        deallocate(drho)

! Write out the maps of integrated properties.
        where(twt==0. .or. cover==0.) cdens = 0.
        where(twt==0. .or. cover==0. .or. temp<0.) temp = 0.
        where(twt==0. .or. cover==0. .or. temp<0.) tvar = 0.
        where(twt==0. .or. cover==0. .or. temp<0.) tskew = 0.
        where(twt==0. .or. cover==0. .or. temp<0.) tkurt = 0.
        call writeimage2d(cdensfile,cdens,nximage,nyimage,'10^20 cm^-2 ', &
            ctype1,ctype2,x0,y0,glon0,glat0,-pixd,pixd,crota2,status) 
        call writeimage2d(tempfile,temp,nximage,nyimage,'K           ', &
            ctype1,ctype2,x0,y0,glon0,glat0,-pixd,pixd,crota2,status) 
        call writeimage2d(tvarfile,tvar,nximage,nyimage,'K^2         ', &
            ctype1,ctype2,x0,y0,glon0,glat0,-pixd,pixd,crota2,status) 
        call writeimage2d(tskewfile,tskew,nximage,nyimage,'K^3         ', &
            ctype1,ctype2,x0,y0,glon0,glat0,-pixd,pixd,crota2,status) 
        call writeimage2d(tkurtfile,tkurt,nximage,nyimage,'K^4         ', &
            ctype1,ctype2,x0,y0,glon0,glat0,-pixd,pixd,crota2,status) 

        deallocate(cover)
        deallocate(cdens)
        deallocate(temp)
        deallocate(tvar)
        deallocate(tskew)
        deallocate(tkurt)
        deallocate(twt)

! Output a table of reduced chi squared values.
        outfile = removeblanks(fieldname//'_rchisq.txt')
        print *,'             ',outfile(1:60)
        open (unit=1,file=outfile,status='UNKNOWN')
        bin = 0.1
        nbins = 100
        allocate (counter(nfields))
        do i = 1,nbins
            bini = (i-1)*bin
            counter = 0
            where(rchisq >= bini-bin/2. .and. rchisq < bini+bin/2.) &
                counter = 1
            write(1,'(f12.2,i12)') bini,sum(counter)
        enddo
        close(1)
        deallocate(counter)
        deallocate(rchisq)

        call date_and_time(date,time,zone,values)
        write(*,'("End PPMOSAIC_MULTITILE on ",a8," at ",a10)') date,time
        stop

        end program ppmosaic_multitile
