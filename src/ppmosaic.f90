        subroutine ppmosaic(fieldname)

! Construct a mosaic of PPMAP images.
!
! Input parameters:
!	fieldname	=	name of field, e.g. 'l224'
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

        character(len=4),  allocatable :: bands(:)

        real(4),   allocatable :: a(:,:), buffer(:)
        real(4),   allocatable :: Tgrid(:), betagrid(:), cctab(:,:), cct(:)
        real(4),   allocatable :: wavelengths(:),beamsizes(:)
        real(4),   allocatable :: cover(:,:),rms(:)
        real(4),   allocatable :: rchisqind(:),chisqind(:),rcs(:)
        real(4),   allocatable :: cdens(:,:),temp(:,:),betamean(:,:)
        real(4),   allocatable :: tvar(:,:),tskew(:,:),tkurt(:,:),rcsmap(:,:)
        real(4),   allocatable :: rho(:,:,:,:),drho(:,:,:,:)
        real(4),   allocatable :: rh(:,:,:,:),rht(:,:),rhv(:,:),rhoff(:,:,:,:)
        real(4),   allocatable :: rhot(:,:),rhotx(:,:)
        real(4),   allocatable :: drhot(:,:),drhotc(:,:)
        real(4),   allocatable :: rhsum(:,:),uncscale(:,:)
        real(4),   allocatable :: ukern(:,:)
        real(4),   allocatable :: backcorr(:,:,:),nrat(:)
        real(4),   allocatable :: backarr(:,:),sba(:,:),bkern(:,:)
        real(4),   allocatable :: deltaset(:,:,:,:),dback(:,:,:)
        real(4),   allocatable :: diffset(:,:,:,:),diffsetok(:),diffs(:)
        real(4),   allocatable :: nplaneset(:,:,:,:),nplane(:,:),nplaneok(:,:)
        real(4),   allocatable :: wx(:),wy(:),w(:,:),twt(:,:),kern(:,:)
        real(4),   allocatable :: uncback(:,:),uncpeak(:,:),rchisq(:)
        integer(4),allocatable :: psfsizes(:),nobsind(:),ntotind(:)
        integer(4),allocatable :: x0set(:),y0set(:)
        integer(4),allocatable :: nvalues(:,:),noverset(:,:,:),noset(:,:,:,:)
        integer(4),allocatable :: ifield(:),jfield(:),rmask(:,:),cmask(:,:)
        integer(4),allocatable :: noplaneset(:,:,:,:),noplane(:,:)
        integer(4),allocatable :: noplaneok(:,:)
        integer(4),allocatable :: nodiffset(:,:,:,:)
        integer(4),allocatable :: wt(:,:),wto(:,:),overlap(:,:)
        integer(4),allocatable :: ok(:,:,:,:),counter(:)

        character(len=8)       :: date, ctype1, ctype2
        character(len=10)      :: time
        character(len=5)       :: zone
        character(len=20)      :: fieldname
        character(len=80)      :: line,imagefile,inname,outfile,removeblanks
        character(len=80)      :: rhofile,sigrhofile
        character(len=80)      :: cdensfile,tempfile,betafile
        character(len=80)      :: tvarfile,tskewfile,tkurtfile,rcsfile
        integer, dimension (8) :: values
        integer, dimension (1) :: npk
        integer    status,cwin,dx0,dy0,x0diff,y0diff,edge
        real(4)    Nexp, Nprior, kappa300
        logical    incomplete, notfinished, bgcorr, verbose

        bgcorr = .true.         ! Apply background correction
        verbose = .false.       ! Minimal narration

        call date_and_time(date,time,zone,values)
        write(*,'("Begin PPMOSAIC on ",a8," at ",a10)') date,time
        rnan = -1.e35
        hnan = rnan/2.
        aln2 = log(2.)
        fnan = 0
        fnan = 0/fnan

! Read the input parameters file.
        line = removeblanks(fieldname//'_ppmap.inp')
        write (*,'("Reading input parameters in ",a)') line(1:40)
        open (unit=2, form='formatted', file=line,status='old', &
            action='read')
        read(2,*) nbands,nfields,Nt,nbeta,icc,maxiterations,tdistance,kappa300,&
            eta,beta0,sigbeta,rchisqconv,trimlev,ihp
        distance = 1000.                ! scale to true distance at end
        allocate (betagrid(nbeta))
        read(2,*) betagrid
        allocate (psfsizes(nbands))
        read(2,*) psfsizes
        allocate (Tgrid(Nt))
        read(2,*) Tgrid
        allocate (cctab(Nt,nbands))
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

        deallocate (psfsizes)
        allocate (wavelengths(nbands))
        allocate (beamsizes(nbands))
        allocate (rchisqind(nbands))
        allocate (chisqind(nbands))
        allocate (nobsind(nbands))
        allocate (ntotind(nbands))
        allocate (bands(nbands))

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
            write(bands(i),'(i4.4)') nint(wl)
            beamsizes(i) = 2.*cdelt2*3600.
            if (i==1) then
                xref = crpix1 - 1.
                yref = crpix2 - 1.
            endif
        enddo
        close(2)
        deallocate(a)
        call nmedian(beamsizes,nbands,beammed)
        print *,'Tgrid:',Tgrid
        print *,'Betagrid:',betagrid
        print *,nfields,' fields in mosaic'

! Read coverage map.
        allocate (cover(nximage,nyimage))
        line = removeblanks(fieldname//'_coverage.fits')
        call readimage_wcs(line,cover,nximage,nyimage,buffer, &
                ctype1,ctype2,crpix1out,crpix2out,crval1out,crval2out, &
                cdelt1out,cdelt2out,crota2,pixel,wl,sig,status)
        if (status > 0) then
            print *,'Error reading coverage map'
            stop
        endif
        deallocate(buffer)

! Set up arrays for mosaicing.
        allocate (cdens(nximage,nyimage))
        allocate (temp(nximage,nyimage))
        allocate (betamean(nximage,nyimage))
        allocate (tvar(nximage,nyimage))
        allocate (tskew(nximage,nyimage))
        allocate (tkurt(nximage,nyimage))
        allocate (rcsmap(nximage,nyimage))
        allocate (rho(nximage,nyimage,Nt,nbeta))
        allocate (drho(nximage,nyimage,Nt,nbeta))
        allocate (rchisq(nfields))
        allocate (nrat(nfields))
        allocate (x0set(nfields))
        allocate (y0set(nfields))

        incomplete = .true.
        do while (incomplete)
          nfound = 0
          do n = 1,nfields
            write(inname,'(a40,"_",i5.5)') fieldname,n
            line = removeblanks(fieldname//'_results/'//inname//'_rho.fits')
            call readheader3d(line,nx,ny,nz,status)
            if (status==0) nfound = nfound + 1
          enddo
          if (nfound==nfields) incomplete = .false.
          call sleep(10)
        enddo
        write(*,'(a)') 'PPMAP output files all present'
        allocate (rh(nx,ny,nz,nbeta))
        allocate (rht(nx,ny))
        allocate (rhoff(nx,ny,Nt,nbeta))
        allocate (rhot(nx,ny))
        allocate (rhotx(nx,ny))
        allocate (rcs(nfields))
        rcs = 0.
        nbuffer = nx*ny*Nt*nbeta
        allocate (buffer(nbuffer))
        chisqind = 0.
        ntotind = 0
        nlimit = 0
        m = 0

        do n = 1,nfields
            write(inname,'(a40,"_",i5.5)') fieldname,n
            line = removeblanks(fieldname//'_results/'//inname//'_rho.fits')
            call readrho(line,rh,nx,ny,Nt,nbeta,buffer,ix0,iy0, &
                et,rhoobj,Nprior,Nexp,nsteps,rchi,status)
            if (status==0) then
                x0set(n) = ix0
                y0set(n) = iy0
                if (rchi > 0.) then
                    m = m+1
                    rchisq(m) = rchi
                endif
            else
                print *,'Could not read '//line
                stop
            endif
            if (nsteps==maxiterations) nlimit = nlimit + 1
            write(inname,'(a40,"_",i5.5)') fieldname,n
            line = removeblanks(fieldname//'_results/'//inname//'_rchisq.txt')
            open (unit=10,file=line,iostat=istat,status='UNKNOWN')
            if (istat==0) then
              read(10,*,iostat=istat) rchisqind
              if (istat==0) then
                read(10,*) rchisqind
                read(10,*) nobsind
                chisqind = chisqind + rchisqind*nobsind
                ntotind = ntotind + nobsind
                rcs(n) = sum(rchisqind*nobsind)/sum(nobsind)
              else
                close(10)
              endif
            endif
        enddo

        rchisqind = chisqind/ntotind
        rchisq0 = sum(chisqind)/sum(ntotind)
        mfields = m
        print *,'Fraction of subfields reaching iteration limit =', &
            float(nlimit)/nfields

        dx0 = 10000
        dy0 = 10000
        do n = 2,nfields
            x0diff = abs(x0set(n) - x0set(n-1))
            y0diff = abs(y0set(n) - y0set(n-1))
            if (x0diff > 0 .and. x0diff < dx0) dx0 = x0diff
            if (y0diff > 0 .and. y0diff < dy0) dy0 = y0diff
        enddo

        edge = nint(beammed/pixel)
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
        allocate (backcorr(nfields,Nt,nbeta))
        backcorr = 0.

      if(bgcorr .and. nfields >= 9) then
        print *,'Calculating background correction using a',trimlev, &
            'sigma trim level'
        nxfields = (maxval(x0set) - minval(x0set))/dx0 + 1
        nyfields = (maxval(y0set) - minval(y0set))/dy0 + 1
        allocate (backarr(nxfields,nyfields))
        allocate (sba(nxfields,nyfields))
        allocate (nvalues(nxfields,nyfields))
        allocate (noverset(nfields,Nt,nbeta))
        allocate (deltaset(9*nx*ny,nfields,Nt,nbeta))
        allocate (noset(9*nx*ny,nfields,Nt,nbeta))
        allocate (ifield(nfields))
        allocate (jfield(nfields))
        allocate (dback(nfields,Nt,nbeta))
        ifield = (x0set - minval(x0set))/dx0
        jfield = (y0set - minval(y0set))/dy0

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
                write(inname,'(a40,"_",i5.5)') fieldname,n
                line = removeblanks(fieldname//'_results/'//inname//'_rho.fits')
                call readrho(line,rh,nx,ny,Nt,nbeta,buffer, &
                    ix0,iy0,et,rhoobj,Nprior,Nexp,nsteps,rchi,status)
                do ib = 1,nbeta
		do it = 1,Nt
		    rht = rh(1:nx,1:ny,it,ib) + backcorr(n,it,ib)
		    rhv = rnan
                    do j = edge+1,ny-edge
                    do i = edge+1,nx-edge 
		        rhv(nx+i,ny+j) = rht(i,j)
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
		      if (ifo >= 0 .and. ifo < nxfields .and. jfo >= 0 .and. &
		         jfo < nyfields .and. (ioff /= 0 .or. joff /= 0)) then 
			no = nvalues(ifo+1,jfo+1)
                        write(inname,'(a40,"_",i5.5)') fieldname,no
                        line = removeblanks(fieldname//'_results/'//inname// &
                            '_rho.fits')
                        call readrho(line,rhoff,nx,ny,Nt,nbeta,buffer, &
                            ix0,iy0,et,rhoobj,Nprior,Nexp,nsteps,rchi,status)
			rhot = rhoff(1:nx,1:ny,it,ib) + backcorr(no,it,ib)
			rhotx = rnan
                        do j = edge+1,ny-edge
                        do i = edge+1,nx-edge
			    rhotx(i,j) = rhot(i,j)
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
			noverset(n,it,ib) = nok
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
                                 deltaset(iok,n,it,ib) = diffsetok(iok)
			         noset(iok,n,it,ib)=nodiffset(i3,j3,ioff+2,joff+2)
                               endif
                             enddo
                             enddo
                            endif
                          enddo
                          enddo
                            
			  call trimavg(diffsetok,nok,dmean,dsig,trimlev)
	  		  dback(n,it,ib) = dmean
                          deallocate(diffsetok)
			endif
		    else 
                        if (it==1) print *,'Field',n,' has no overlap'
                    endif
      		enddo
                enddo
	    enddo
	    rmsprev = sqrt(sum(dback**2))
	  else 
	    do n = 1,nfields
            do ib = 1,nbeta
            do it = 1,Nt
		nover = noverset(n,it,ib)
		if (nover /= 0) then
                    allocate (diffs(nover))
		    do i = 1,nover 
			no = noset(i,n,it,ib)
			diffs(i) = deltaset(i,n,it,ib) + &
			    backcorr(no,it,ib) - backcorr(n,it,ib)
		    enddo
		    call trimavg(diffs,nover,dmean,dsig,trimlev)
	  	    dback(n,it,ib) = dmean
                    deallocate(diffs)
		endif
	    enddo
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
        deallocate(dback)
        deallocate(wto)
        deallocate(wt)

! Remove any large scale trends from background correction and adjust to
! zero mean at all temperatures.
        iwb = 6                         ! FWHM of Gaussian smoothing function
        nbk = 2*iwb + 1
        allocate (bkern(nbk,nbk))
        bfwhmsq = iwb**2
        do j = 1,nbk
        do i = 1,nbk
            bkern(i,j) = exp(-4.*aln2*((i-iwb-1)**2 + (j-iwb-1)**2)/bfwhmsq)
        enddo
        enddo

        do ib = 1,nbeta
        do it = 1,Nt 
            backarr = 0.
            do n = 1,nfields
                backarr(ifield(n)+1,jfield(n)+1) = backcorr(n,it,ib)
            enddo
            do j = 1,nyfields
            do i = 1,nxfields
                ilo = i-iwb
                ihi = i+iwb
                jlo = j-iwb
                jhi = j+iwb
                bsum = 0.
                bwt = 0.
                do jj = jlo,jhi
                do ii = ilo,ihi
                    if (ii>0 .and. jj>0 .and. &
                      ii <= nxfields .and. jj <= nyfields) then
                        bsum = bsum + backarr(ii,jj)* &
                            bkern(iwb+1+ii-i, iwb+1+jj-j)
                        bwt = bwt + bkern(iwb+1+ii-i, iwb+1+jj-j)
                    endif
                enddo
                enddo
                sba(i,j) = bsum/bwt
            enddo
            enddo
            backarr = backarr - sba
            do n = 1,nfields
                backcorr(n,it,ib) = backarr(ifield(n)+1,jfield(n)+1)
            enddo
            bcmean = sum(backcorr(1:nfields,it,ib))/nfields
            do n = 1,nfields
                backcorr(n,it,ib) = backcorr(n,it,ib) - bcmean
            enddo
        enddo
        enddo
        deallocate(backarr)
        deallocate(sba)
        deallocate(bkern)
        deallocate(ifield)
        deallocate(jfield)
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
        rcsmap = 0.
        twt = 0.
        m = 0

        do n = 1,nfields
            write(inname,'(a40,"_",i5.5)') fieldname,n
            line = removeblanks(fieldname//'_results/'//inname//'_rho.fits')
            call readrho(line,rh,nx,ny,Nt,nbeta,buffer,ix0,iy0, &
                et,rhoobj,Nprior,Nexp,nsteps,rchi,status)
            if (rchi > 0.) then
                m = m+1
	        nrat(m) = Nexp/Nprior
            endif

	    if(ix0>=0 .and. iy0>=0 .and. ix0+nx<nximage .and. iy0+ny<nyimage) &
              then
                do j = iy0+1,iy0+ny
                do i = ix0+1,ix0+nx
                    do ib = 1,nbeta
	            do it = 1,Nt
                        rho(i,j,it,ib) = rho(i,j,it,ib) + w(i-ix0,j-iy0) *  &
                        (rh(i-ix0,j-iy0,it,ib) + backcorr(n,it,ib))
                    enddo
                    enddo
                    rcsmap(i,j) = rcsmap(i,j) + w(i-ix0,j-iy0)*rcs(n)
	            twt(i,j) = twt(i,j) + w(i-ix0,j-iy0)
                enddo
                enddo
	    endif
        enddo

        deallocate(rh)
        deallocate(buffer)
        deallocate(w)
        deallocate(backcorr)
        deallocate(rcs)
        if (verbose) print *,'Mosaic completed'

! Construct "single object" kernel.
        ik = nint(1.5*fwhmobj)
        nk = 2*ik+1
        allocate (kern(nk,nk))
        do j = 1,nk
        do i = 1,nk
            kern(i,j) = exp(-4.*(aln2/fwhmobj**2)*((i-ik-1)**2 + (j-ik-1)**2))
        enddo
        enddo
        kern = kern/sum(kern)

! Calculate image cube and convolve with single-object profile.
        allocate(rhot(nximage,nyimage))
        allocate(rhotx(nximage,nyimage))
        do ib = 1,nbeta
        do it = 1,Nt
            do j = 1,nyimage
            do i = 1,nximage
                if (twt(i,j) /= 0.) rho(i,j,it,ib) = rho(i,j,it,ib)/twt(i,j)
            enddo
            enddo
            rhot = rho(1:nximage,1:nyimage,it,ib)
            where(rhot < 0.) rhot = 0.
            call convol(rhot,kern,rhotx,nximage,nyimage,nk)
            do j = 1,nyimage
            do i = 1,nximage
                rho(i,j,it,ib) = rhotx(i,j)
            enddo
            enddo
        enddo
        enddo
        deallocate(rhot)
        deallocate(rhotx)
        do j = 1,nyimage
        do i = 1,nximage
            if (twt(i,j) /= 0.) rcsmap(i,j) = rcsmap(i,j)/twt(i,j)
        enddo
        enddo

! Calculate uncertainties.
        print *,'Calculating uncertainties'
        npk = maxloc(nrat)
        npeak = npk(1)
        allocate (uncback(Nt,nbeta))
        unback = 0.
        allocate (uncpeak(Nt,nbeta))
        uncpeak = 0.
        write(inname,'(a40,"_",i5.5)') fieldname,npeak
        line = removeblanks(fieldname//'_results/'//inname//'_rho.fits')
        call pperr(line,wavelengths,nbands,Tgrid,Nt,betagrid,nbeta,sigbeta, &
             cctab,pixel,kappa300,uncback,uncpeak)

        deallocate(bands)
        deallocate(cctab)
        allocate (rhsum(nximage,nyimage))
        do j = 1,nyimage
        do i = 1,nximage
            rsum = 0.
            do ib = 1,nbeta
            do it = 1,Nt
                rsum = rsum + rho(i,j,it,ib)
            enddo
            enddo
            rhsum(i,j) = rsum
        enddo
        enddo

! Construct convolution kernel.
        iuk = nx/2
        nuk = 2*iuk + 1
        allocate (ukern(nuk,nuk))
        ufwhmsq = (nuk/2.)**2
        do j = 1,nuk
        do i = 1,nuk
            ukern(i,j) = exp(-4.*aln2*((i-iuk-1)**2 + (j-iuk-1)**2)/ufwhmsq)
        enddo
        enddo

! Do the convolution.
        allocate (uncscale(nximage,nyimage))
        call convol(rhsum,ukern,uncscale,nximage,nyimage,nuk)

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

! Interpolate between uncback and uncpeak.
        where(uncscale < 0.) uncscale = 0.
        uncm = maxval(uncscale)
        uncscale = sqrt(uncscale/uncm)
        
        do ib = 1,nbeta
        do it = 1,Nt
	    write(*,'("T[K],beta:",f6.2,f5.2,": bkgnd & peak uncs:",' &
                //'f10.4,",",f12.4,"  x 10^20 cm^-2")') &
                Tgrid(it),betagrid(ib),uncback(it,ib),uncpeak(it,ib)
            do j = 1,nyimage
            do i = 1,nximage
                drho(i,j,it,ib) = &
                    max(uncpeak(it,ib)*uncscale(i,j), uncback(it,ib))
            enddo
            enddo
        enddo
        enddo
        deallocate(uncpeak)
        deallocate(uncscale)
        deallocate(ukern)

! Calculate integrated column density, mean temperature, and mean beta.
        cdens = 0.
        temp = 0.
        betamean = 0.
        do ib = 1,nbeta
        do it = 1,Nt
	    cdens = cdens + rho(1:nximage,1:nyimage,it,ib)
	    temp = temp + rho(1:nximage,1:nyimage,it,ib)*Tgrid(it)
	    betamean = betamean + rho(1:nximage,1:nyimage,it,ib)*betagrid(ib)
        enddo
        enddo
        where(cdens > 0.) temp = temp/cdens
        where(cdens > 0.) betamean = betamean/cdens

! Calculate temperature variance, skew, and kurtosis.
        tvar = 0.
        tskew = 0.
        tkurt = 0.
        do ib = 1,nbeta
        do it = 1,Nt
	    tvar = tvar + rho(1:nximage,1:nyimage,it,ib)*(Tgrid(it)-temp)**2
	    tskew = tskew + rho(1:nximage,1:nyimage,it,ib)*(Tgrid(it)-temp)**3
	    tkurt = tkurt + rho(1:nximage,1:nyimage,it,ib)*(Tgrid(it)-temp)**4
        enddo
        enddo
        where(cdens > 0.) tvar = tvar/cdens
        where(cdens > 0.) tskew = tskew/cdens
        where(cdens > 0.) tkurt = tkurt/cdens

! Output a metadata file.
        outfile = removeblanks(fieldname//'_metadata.txt')
        open (unit=1,file=outfile,status='UNKNOWN')
        write(1,'(a)') '\ Temperature grid:'
        write(1,'(a2,20f8.1)') '\ ',Tgrid
        write(1,'(a)') '\ Beta grid:'
        write(1,'(a2,20f8.1)') '\ ',betagrid
        write(1,'(a2)') '\ '
        write(1,'(a2,i8," fields in mosaic")') '\ ',nfields
        write(1,'(a2)') '\ '
        write(1,'(a)') '\ Reduced chi squared versus band:'
        write(1,'(a2,20f8.1)') '\ ',wavelengths
        write(1,'(a2,20f8.1)') '\ ',rchisqind
        write(1,'(a2,20i8)') '\ ',ntotind
        write(1,'(a2)') '\ '
        write(1,'("\ Global value of reduced chi squared =",f8.2)') rchisq0
        write(1,'(a2)') '\ '
        write(1,'(a)') '\ Histogram of tile-based reduced chi squared:'
        write(1,'(a2)') '\ '
        if (mfields > 1) then
            call nmedian(rchisq,mfields,rchisqmed)
        else
            rchisqmed = rchisq(1)
        endif
        binmax = max(3.*rchisqmed, 10.)
        nbins = max(min(100, nfields/6),2)
        bin = binmax/nbins
        allocate (counter(nfields))
        if (nfields > mfields) then
            do i = mfields+1, nfields
                rchisq(i) = -999.
            enddo
        endif
        ihistmax = -1
        do i = 1,nbins
            bini = (i-1)*bin
            counter = 0
            where(rchisq >= bini-bin/2. .and. rchisq < bini+bin/2.) &
                counter = 1
            ihist = sum(counter)
            if (ihist > ihistmax) then
                ihistmax = ihist
                rchisqmode = bini
            endif
            write(1,'(f12.2,i12)') bini,ihist
        enddo
        close(1)
        deallocate(counter)
        deallocate(wavelengths)
        deallocate(beamsizes)
        deallocate(chisqind)
        deallocate(rchisqind)
        deallocate(nobsind)
        deallocate(ntotind)

        if (mfields > 1) then
            call nmedian(nrat,mfields,ratmed)
            print *,'Median ratio of Npost/Nprior      =',ratmed
            print *,'Median reduced chi squared        =',rchisqmed
            if (rchisqmode /= 0. .and. mfields >= 100) &
                print *,'Mode of reduced chi squared       =',rchisqmode
        else
            print *,'Ratio of Npost/Nprior             =',nrat(1)
            print *,'Reduced chi squared               =',rchisq0
        endif
        deallocate(rchisq)
        deallocate(nrat)

        write(*,'(" Total mass                        =",1pe12.3," Msun")') &
            sum(rho)*(1.e20/Msun)*mu*mH*(pixel*(dtor/3600.)*tdistance*pc)**2 
        rhofile = removeblanks(fieldname//'_tdenscube.fits')
        sigrhofile = removeblanks(fieldname//'_sigtdenscube.fits')
        cdensfile = removeblanks(fieldname//'_cdens.fits')
        tempfile = removeblanks(fieldname//'_temp.fits')
        betafile = removeblanks(fieldname//'_beta.fits')
        tvarfile = removeblanks(fieldname//'_tvar.fits')
        tskewfile = removeblanks(fieldname//'_tskew.fits')
        tkurtfile = removeblanks(fieldname//'_tkurt.fits')
        rcsfile = removeblanks(fieldname//'_rchisq.fits')
        print *,'Writing out: ',cdensfile(1:60)
        print *,'             ',tempfile(1:60)
        print *,'             ',betafile(1:60)
        print *,'             ',tvarfile(1:60)
        print *,'             ',tskewfile(1:60)
        print *,'             ',tkurtfile(1:60)
        print *,'             ',rcsfile(1:60)
        print *,'             ',rhofile(1:60)
        print *,'             ',sigrhofile(1:60)
        print *,'             ',outfile(1:60)

! Write out image cube of differential column density.
        allocate (rhot(nximage,nyimage))
        allocate (drhot(nximage,nyimage))
        do ib = 1,nbeta
        do it = 1,Nt
	    rhot = rho(1:nximage,1:nyimage,it,ib)
            where(twt==0. .or. cover==0.) rhot = fnan
	    drhot = drho(1:nximage,1:nyimage,it,ib)
            where(twt==0. .or. cover==0.) drhot = fnan
            do j = 1,nyimage
            do i = 1,nximage
	        rho(i,j,it,ib) = rhot(i,j)
	        drho(i,j,it,ib) = drhot(i,j)
            enddo
            enddo
        enddo
        enddo
        deallocate(rhot)
        deallocate(drhot)
        glon0 = crval1
        glat0 = crval2
        pixd = pixel/3600.
        dkpc = tdistance/1000.
        call writerho(rhofile,rho,nximage,nyimage,Nt,nbeta,ctype1,ctype2, &
            crpix1out,crpix2out,crval1out,crval2out,cdelt1out,cdelt2out,crota2,&
            1,1,Tgrid(1),Tgrid(Nt),eta,rhoobj,Nprior,ratmed*Nprior,dkpc, &
            maxiterations,rchisq0,status)
        call writerho(sigrhofile,drho,nximage,nyimage,Nt,nbeta,ctype1,ctype2, &
            crpix1out,crpix2out,crval1out,crval2out,cdelt1out,cdelt2out,crota2,&
            1,1,Tgrid(1),Tgrid(Nt),eta,rhoobj,Nprior,ratmed*Nprior,dkpc, &
            maxiterations,rchisq0,status)
        deallocate(rho)
        deallocate(drho)

! Write out the maps of integrated properties.
        where(twt==0. .or. cover==0.) cdens = fnan
        where(twt==0. .or. cover==0. .or. temp<0.) betamean = fnan
        where(twt==0. .or. cover==0. .or. temp<0.) temp = fnan
        where(twt==0. .or. cover==0. .or. temp<0.) tvar = fnan
        where(twt==0. .or. cover==0. .or. temp<0.) tskew = fnan
        where(twt==0. .or. cover==0. .or. temp<0.) tkurt = fnan
        call writeimage2d(cdensfile,cdens,nximage,nyimage,'10^20 cm^-2 ', &
            ctype1,ctype2,crpix1out,crpix2out,crval1out,crval2out, &
            cdelt1out,cdelt2out,crota2,status) 
        call writeimage2d(betafile,betamean,nximage,nyimage,'K           ', &
            ctype1,ctype2,crpix1out,crpix2out,crval1out,crval2out, &
            cdelt1out,cdelt2out,crota2,status) 
        call writeimage2d(tempfile,temp,nximage,nyimage,'K           ', &
            ctype1,ctype2,crpix1out,crpix2out,crval1out,crval2out, &
            cdelt1out,cdelt2out,crota2,status) 
        call writeimage2d(tvarfile,tvar,nximage,nyimage,'K^2         ', &
            ctype1,ctype2,crpix1out,crpix2out,crval1out,crval2out, &
            cdelt1out,cdelt2out,crota2,status) 
        call writeimage2d(tskewfile,tskew,nximage,nyimage,'K^3         ', &
            ctype1,ctype2,crpix1out,crpix2out,crval1out,crval2out, &
            cdelt1out,cdelt2out,crota2,status) 
        call writeimage2d(tkurtfile,tkurt,nximage,nyimage,'K^4         ', &
            ctype1,ctype2,crpix1out,crpix2out,crval1out,crval2out, &
            cdelt1out,cdelt2out,crota2,status) 
        call writeimage2d(rcsfile,rcsmap,nximage,nyimage,'            ', &
            ctype1,ctype2,crpix1out,crpix2out,crval1out,crval2out, &
            cdelt1out,cdelt2out,crota2,status) 
        deallocate(cover)
        deallocate(cdens)
        deallocate(betamean)
        deallocate(temp)
        deallocate(tvar)
        deallocate(tskew)
        deallocate(tkurt)
        deallocate(rcsmap)
        deallocate(twt)
        deallocate(Tgrid)

        call date_and_time(date,time,zone,values)
        write(*,'("End PPMOSAIC on ",a8," at ",a10)') date,time
        return

        end subroutine ppmosaic
