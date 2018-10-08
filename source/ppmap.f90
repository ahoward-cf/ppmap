        program ppmap

! Mapping of extended dust structures on a set of predefined temperature planes
! using point process technique. All observed images must be input with arrays
! of the same size. They can have different pixel sizes (ideally corresponding
! to Nyquist sampling) but must be co-registered with respect to their 
! reference pixels. The PSF images must be square and all have the same 
! sampling interval, but not necessarily equal to the observational pixel size.

! Input parameters:
!	fieldname	=	name of field, e.g. 'l224'
!	isf1, isf2	=       range of subfield numbers. The subfields
!                       	are defined by an index contained within
!                       	the input parameters file.
!       mosaic          =       if present, calculate mosaic

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
        real,      parameter   :: reflambda= 300.	! opacity ref wavelen 
        real,      parameter   :: tmax = 1.		! termination value of t
        real(8),   parameter   :: dtor    = 0.0174533   ! degrees to radians

        character(len=80), allocatable :: obsimages(:)
        real(8),   allocatable :: grid(:,:)
        real(4),   allocatable :: a(:,:), amask(:,:)
        real(4),   allocatable :: buffer(:), outcube(:,:,:,:), rchisqind(:)
        real(4),   allocatable :: data(:), wldata(:), signu2(:), cover(:,:)
        real(4),   allocatable :: rho(:), rhox(:), rho0(:), rhosave(:)
        real(4),   allocatable :: B(:,:), Bsum(:), Brho(:), Bx(:,:)
        real(4),   allocatable :: phi1(:), rhogrid(:,:,:,:)
        real(4),   allocatable :: kern(:,:),Tgrid(:), cctab(:,:), cct(:)
        real(4),   allocatable :: refmodel(:,:)
        real(4),   allocatable :: wavelengths(:), betagrid(:), pixset(:)
        real(4),   allocatable :: psfset(:,:,:),psf(:,:)
        real(4),   allocatable :: obsset(:,:,:),obsmask(:,:,:),sigset(:,:,:)
        integer(4),allocatable :: psfsizes(:)
        integer(4),allocatable :: support(:), mskdata(:), emsk(:), nobsind(:)
        integer(4),allocatable :: xlo(:), xhi(:), ylo(:), yhi(:)
        character(len=8)       :: date, ctype1, ctype2
        character(len=10)      :: time
        character(len=5)       :: zone
        character(len=20)      :: fieldname,sf1,sf2
        character(len=80)      :: line,outfile,removeblanks
        integer, dimension (8) :: values

        real(8)    x,y,z,zx,zy,zxx,zyy,zxy,nref8
        real(8)    t_0, t_1, ra, dec, xpixr, ypixr
        real(4)    Nexp, Nprior, mag, lambda, kappa300
        integer(4) status, count
        logical  finished, ok, stepok, finetune, converged, diverged
        logical  makemaps, makemosaic

        call date_and_time(date,time,zone,values)
        write(*,'("Begin PPMAP on ",a8," at ",a10)') date,time

! Read command line.
        m = iargc()
        if (m < 3 .or. m > 4) then
            print *,'Syntax: ppmap fieldname isf1 isf2 [mosaic]'
            stop
        else if (m==4) then
            makemosaic = .true.
        else
            makemosaic = .false.
        endif
        
        call getarg(1,fieldname)
        call getarg(2,sf1)
        call getarg(3,sf2)
        read(sf1,*) isf1
        read(sf2,*) isf2

        if (isf1==0 .and. isf2==0) then
            makemaps = .false.
        else
            makemaps = .true.
        endif

        if (makemaps) then
            write(*,'("Will make maps for subfield numbers: ",i4," -",i4)') &
                isf1,isf2
            if (makemosaic) print *,'and then construct mosaic'
        else
            if (makemosaic) then
                print *,'(proceed directly to mosaicing step)'
            else
                print *,'No tasks given'
            endif
        endif

if (makemaps) then

! Read the input parameters file.
        line = removeblanks(fieldname//'_ppmap.inp')
        write (*,'("Reading input parameters in ",a)') line
        open (unit=2, form='formatted', file=line,status='old', &
            action='read')
        read(2,*)nbands,nsf,Nt,nbeta,icc,maxiterations,tdistance,kappa300,eta, &
            betamean,sigbeta,rchisqconv,trimlev,ihp

        distance = 1000.                ! scale results to true distance at end
        allocate (betagrid(nbeta))
        read(2,*) betagrid
        allocate (psfsizes(nbands))
        read(2,*) psfsizes
        allocate (Tgrid(Nt))
        read(2,*) Tgrid

        if (icc /= 0) then      ! Colour corrections to be applied
            allocate (cctab(Nt,nbands))
            allocate (cct(Nt))
            do i = 1,nbands
                read(2,*) cct
                do it = 1,Nt
                    cctab(it,i) = cct(it)
                enddo
            enddo
            deallocate(cct)
        endif

        allocate (obsimages(nbands))
        do i = 1,nbands
	    read(2,'(a)') line
	    obsimages(i) = line
	    if (i==1) then 
                call readheader(line,nx,ny,crpix1,crpix2,cdelt2,status)
                if (status > 0) then
                    print *,'Error reading header of '//line
                    stop
                endif
                call readheader('mask_'//line,nxm,nym,p1m,p2m,d2m,status)
                if (status > 0 .or. nxm /= nx .or. nym /= ny) then
                    print *,'Error reading header of mask_'//line
                    stop
                endif
                allocate (obsset(nbands,nx,ny))
                allocate (obsmask(nbands,nx,ny))
                allocate (sigset(nbands,nx,ny))
	    endif
        enddo
        allocate(xlo(nsf))
        allocate(xhi(nsf))
        allocate(ylo(nsf))
        allocate(yhi(nsf))

        do i = 1,nsf
            read(2,*) isn,xlo(i),xhi(i),ylo(i),yhi(i)
        enddo
        close(2)

        print *,'Opacity index grid:  '
        print *,betagrid
        if (nbeta > 1) then
            print *,'Parameters of beta prior:'
            print *,betamean,sigbeta
        endif
        print *,'Temperature grid [K]:'
        print *,Tgrid
        print *,' '
        if (icc /= 0) print *,'Colour corrections will be applied'
        if (ihp==1) print *, &
            'High pass filtering will be applied where necessary'

! Read PSFs.
        line = removeblanks(fieldname//'_psfset.fits')
        call readheader3d(line,nb,ncellsp,ncellspy,status)
        if (status > 0) then
            print *,'Could not read header for '//line
            stop
        endif
        allocate (psfset(nbands,ncellsp,ncellsp))
        nbuffer = nbands*ncellsp*ncellsp
        allocate (buffer(nbuffer))
        call readimage3d(line,psfset,nbands,ncellsp,ncellsp, &
            buffer,status)
        deallocate(buffer)

! Read in the observational images and bad pixel masks.
        nbuffer = nx*ny
        allocate (buffer(nbuffer))
        allocate (a(nx,ny))
        allocate (amask(nx,ny))
        allocate (pixset(nbands))
        allocate (wavelengths(nbands))
        do lam = 1,nbands
            call readimage_wcs(obsimages(lam),a,nx,ny,buffer,ctype1,ctype2, &
                crpix1,crpix2,crval1,crval2,cdelt1,cdelt2,crota2,pixel,wl,sig, &
                status)
            if (status > 0) then
                print *,'Error reading image:'
                print *,obsimages(lam)
                stop
            endif
            call readimage_wcs('mask_'//obsimages(lam),amask,nx,ny,buffer, &
                ctype1,ctype2,p1m,p2m,v1m,v2m,d1m,d2m,r2m,pm,wlm,sigm,status)
            if (status > 0) then
                print *,'Error reading image:'
                print *,'mask_'//obsimages(lam)
                stop
            endif
	    pixset(lam) = cdelt2*3600.
            crota2 = 0.
            wavelengths(lam) = wl
            print *,'Peak SNR at',wavelengths(lam),' microns =',maxval(a)/sig
            do j = 1,ny
            do i = 1,nx
                obsset(lam,i,j) = a(i,j)
                obsmask(lam,i,j) = amask(i,j)
                sigset(lam,i,j) = sig
            enddo
            enddo
        enddo

! Read coverage map.
        allocate (cover(nx,ny))
        line = removeblanks(fieldname//'_coverage.fits')
        call readimage_basic(line,cover,nx,ny,buffer,status)
        if (status > 0) then
            print *,'Error reading coverage map'
            stop
        endif
        deallocate(buffer)
        deallocate(a)
        deallocate(amask)

! Set up convolution kernel representing the profile of a single object.
        ikern = nint(fwhmobj)
        nkern = 2*ikern+1
        allocate (kern(nkern,nkern))
        fln2 = -4.*log(2.)

        do j = 1,nkern
        do i = 1,nkern
            arg = fln2*(float(i-ikern-1)**2 + float(j-ikern-1)**2)/ &
                fwhmobj**2
            kern(i,j) = exp(arg)
        enddo
        enddo

! At this point, sigset contains a set of Nyquist-sampled images in units
! of Jy/pixel. All images are nx by ny and all have the same reference
! pixel (nx/2, ny/2). The pixel size of ith image is pixset(i) arcsec.

! Set up the imaging grid.
        print *,'Image grid            =',nx,' x',ny,' pixels'
        print *,'Field of view         =',nx*pixel,' x',ny*pixel,' arcsec'
        print *,'Sampling  interval    =',pixel,' arcsec'

! Calculate a set of image cubes corresponding to subfield numbers from
! isf1 to isf2.

    do isf = isf1,isf2
        ilo = xlo(isf)
        ihi = xhi(isf)
        jlo = ylo(isf)
        jhi = yhi(isf)
        print *,' '
        print *,'Begin subfield',isf
        print *,'Portion of image plane:',ilo,ihi,jlo,jhi

! Set up imaging mask.
        iglo = max(ilo, 0)
        ighi = min(ihi, (nx-1))
        jglo = max(jlo, 0)
        jghi = min(jhi, (ny-1))

      if (any(cover(iglo+1:ighi+1,jglo+1:jghi+1) /= 0.)) then ! proceed
        npoints = (ighi-iglo+1)*(jghi-jglo+1)
        allocate (support(npoints))
        k = 0
        m = -1
        do j = 1,ny
        do i = 1,nx
            m = m+1
            if (i>iglo .and. i <= ighi+1 .and. j>jglo .and. j <= jghi+1) then
                k = k+1
                support(k) = m
            endif
        enddo
        enddo

! Count the observations.
        ix = nx/2
        iy = ny/2
        nobs = 0
        do i = 1,nbands
            mag = pixel/pixset(i)
            ilom = max(ix+nint((ilo-ix)*mag)-1, 0)
            ihim = min(ix+nint((ihi-ix)*mag)+1, nx-1)
            jlom = max(iy+nint((jlo-iy)*mag)-1, 0)
            jhim = min(iy+nint((jhi-iy)*mag)+1, ny-1)
            nobs = nobs + (ihim-ilom+1)*(jhim-jlom+1)
        enddo

! Calculate system matrix, B. Also, set up data and noise vectors.
        nstates = npoints*Nt*nbeta
        allocate (data(nobs))
        allocate (wldata(nobs))
        allocate (mskdata(nobs))
        allocate (emsk(nobs))
        allocate (signu2(nobs))
        allocate (rchisqind(nbands))
        allocate (nobsind(nbands))
        allocate (B(nobs,nstates))
        allocate (Bx(nobs,nstates))
        allocate (rho0(nstates))
        allocate (rho(nstates))
        allocate (rhosave(nstates))
        allocate (rhox(nstates))
        allocate (rhogrid(nx,ny,Nt,nbeta))
        allocate (Brho(nobs))
        allocate (Bsum(nstates))
        allocate (phi1(nstates))
        print *,'System matrix is',nobs,'  x',nstates
        mskdata = 0
        emsk = 0
        is = 0

        do lam = 1,nbands
	    mag = pixel/pixset(lam)
            ilom = max(ix+nint((ilo-ix)*mag)-1, 0)
            ihim = min(ix+nint((ihi-ix)*mag)+1, nx-1)
            jlom = max(iy+nint((jlo-iy)*mag)-1, 0)
            jhim = min(iy+nint((jhi-iy)*mag)+1, ny-1)
            margin = max((ihim-ilom+1)/4, (jhim-jlom)/4)
            count = 0
            do j = 1,ny
            do i = 1,nx
                if (i>ilom .and. i <= ihim+1 .and. &
                    j>jlom .and. j <= jhim+1) then
                        count = count+1
                        data(is+count) = obsset(lam,i,j)
                        signu2(is+count) = sigset(lam,i,j)**2
                        if (obsmask(lam,i,j) >= 0.5) then
                            signu2(is+count) = 10.*signu2(is+count)
                            mskdata(is+count) = 1
                        endif
                        if (i-ilom >= margin .and. ihim-i >= margin .and. &
                            j-jlom >= margin .and. jhim-j >= margin) &
                            emsk(is+count) = 1
                        wldata(is+count) = wavelengths(lam)
                endif
            enddo
            enddo

	    ncp = psfsizes(lam)
            allocate (psf(ncp,ncp))
            psf = psfset(lam,1:ncp,1:ncp)
	    nxsub = ihim - ilom + 1
	    nysub = jhim - jlom + 1
            nref = ighi - iglo + 1 + 2*ncp
            allocate (refmodel(nref,nref))
            allocate( grid(0:nref, 0:nref) )
            grid = 0.
	    n = 0

            do kb = 1,nbeta
	    do kt = 1,Nt
                call refmodelcalc(wavelengths(lam),reflambda,betagrid(kb), &
                    Tgrid(kt),pixel,fwhmobj,psf,ncp,refmodel,nref)
                if (icc /= 0) refmodel = refmodel/cctab(kt,lam)
                do jg = 1,nref
                do ig = 1,nref
                    grid(ig,jg) = refmodel(ig,jg)
                enddo
                enddo
                grid = grid/mag**2
                iref = nref/2 + 1
	        do np = 1,npoints
		    n = n + 1
		    i = mod(support(np), nx)
		    j = support(np)/nx
                    ioff = ix + iref - i
                    joff = iy + iref - j
                    rref = nref

! For each pixel (ii,jj) in the observational subgrid, determine the 
! response to a model source at position (i,j) in the imaging grid.
                    k = 0
                    nref8 = nref*1.d0
                    do jj = 1,nysub
                    do ii = 1,nxsub
                        k = k+1
                        x = min(max(ioff+(ilom+ii-ix-1)/mag, 1.), rref)
                        y = min(max(joff+(jlom+jj-iy-1)/mag, 1.), rref)
                        call intrp2(x, 0.d0, nref8, nref, y, 0.d0, &
                            nref8, nref, grid,z,zx,zy,zxx,zyy,zxy)
                        B(is+k,n) = z
                    enddo
                    enddo
	        enddo
	    enddo
            enddo
            deallocate(refmodel)
            deallocate(grid)
            deallocate (psf)
	    is = is + count
        enddo

! Set optical depth increment.
        dtau = 0.01
        B = B*dtau
        finished = .false.

! If necessary, set a default value for the a priori dilution, eta.
        if (eta==0.) then
            Npc = 4                     ! number of unknowns per component
            Nprior = nobs/float(Npc)    ! maximum number of components that
                                        ! the data could constrain
            eta = Nprior/nstates
        else
            Nprior = eta*nstates
        endif

! Calculate the expectation value of the density of points in state space.
        rho0 = eta                      ! a priori density

! Impose a Gaussian prior on beta if there are multiple possible values.
        if (nbeta > 1) then
	    n = 0
            do kb = 1,nbeta
	    do kt = 1,Nt
            do np = 1,npoints
 	        n = n + 1
                rho0(n) = eta * exp(-0.5*((betagrid(kb)-betamean)/sigbeta)**2)
            enddo
	    enddo
            enddo
        endif

        epsilon = 1.
        rchisq0 = 1.e35
        pixarea = (pixel*distance*pc*dtor/3600.)**2
        rhoobj = dtau * tau2mass(kappa300,pixel,distance) * sum(kern) * &
            1.e-20/(mu*mH*(pixarea/Msun)) ! column density of single object
        where (isnan(signu2)) signu2 = 0.
        if (minval(signu2)==0.) then
            finished = .true.
        else
            finished = .false.
        endif
        if (ihp==1) call hpcorr(data,wldata,mskdata,nobs)  ! Subtract mean from
                                                           ! ground-based images
        t_0 = omp_get_wtime()
        call date_and_time(date,time,zone,values)
        write(*,'("Begin integration on ",a8," at ",a10)') date,time

        do while(.not.finished)
            print *,'Integration step size =',epsilon
            Bsum = 0.

            do n = 1,nstates
                Bsum(n) = Bsum(n) + 0.5*sum(B(1:nobs,n)**2/signu2(1:nobs))
            enddo

            rho = rho0
            rhosave = rho0
            rchisq = 1.e35
            rchisqprev = 1.e35
            rchisqmin = 1.e35
            converged = .false.
            diverged = .false.
            it = -1
            ok = .true.
            stepok = .true.
            finetune = .true.
            epsum = 0.
            t = epsum/rchisq0

            do while (t < tmax .and. it < maxiterations .and. &
                .not.converged .and. .not.diverged .and. ok)
        	it = it+1
                Brho = matmul_omp(B,rho)
                if (ihp==1) call hpcorr(Brho,wldata,mskdata,nobs) 
                                                      ! Subtract mean from model
                                                      ! images corresponding to
                                                      ! ground-based data
	        rchisq = sum((data - Brho)**2/signu2)/nobs
	        if (it==0) rchisq0 = rchisq

	        if (mod(it,100)==0) then 
                    if (abs(rchisq - rchisqprev)/rchisq < rchisqconv) &
                        converged = .true.
	            if (rchisq > rchisqprev) diverged = .true.
	            rchisqprev = rchisq
	            if (.not.isnan(rchisq) .and. rchisq < 1.e6) then 
		        Nexp = sum(rho)*epsilon
		        t = epsum/rchisq0
                        t_1 = omp_get_wtime()
                        write(*, '("Step#",i6,"  Nexp =",f9.1,"  t =",'     &
                            //'f8.4,"  Rchisq =",f10.2,"  at",f10.2,"s")')  &
                            it,Nexp,t,rchisq,t_1-t_0
!                        print *,'          rho0(1) = ',rho0(1)
!                        print *,'Iterative rho(1)  = ',rho(1)

	            endif
                    if (rchisq < rchisqmin) then
	                epsilonx = epsilon
	                dtaux = dtau
	                Bx = B
	                rhox = rho
	                epsumx = epsum
	                itx = it
                        rchisqmin = rchisq
                    endif
	        endif
                if (rchisq <= rchisqmin) rhosave = rho

                phi1 = matmul_omp((data-Brho)/signu2,B)
                phi1 = Bsum - phi1
	        rho = rho*(1. - epsilon*phi1)
	        epsum = epsum + epsilon

	        if (mod(it,500)==0 .and. it >= 500 .and. it < maxiterations &
                    .and. stepok) then 
	            stepscale = 1.2
	            epsilon = epsilon*stepscale
	            dtau = dtau*stepscale
	            B = B*stepscale
	            rho = rho/stepscale
	            print *,'Integration step size increased to',epsilon
	        endif

                if (rchisq<0.8 .and. t>0.8*tmax .and. finetune .and. it>1) then
                    finescale = 0.1
	            epsilon = epsilon*finescale
	            dtau = dtau*finescale
	            B = B*finescale
	            rho = rho/finescale
	            print *,'At iteration',it,' decrease step size to',epsilon
                    finetune = .false.
                endif

	        if ((minval(rho) < -1.e-4 .or. sum(rho)==0  &
                  .or. isnan(sum(rho))) .and. finetune) then
	            if (it >= 500) then 
                      if (stepok) then
		        epsilon = epsilonx
		        dtau = dtaux
		        B = Bx
		        rho = rhox
		        epsum = epsumx
		        it = itx
		        print *,'Integration step size restored to ',epsilon
		        if (stepok) diverged = .false.
		        stepok = .false.
                      endif
	            else 
		        rhomin = minval(rho)
		        if (.not.isnan(rhomin)) &
		            print *,'Adjusting stepsize;  min(rho) =',rhomin
		        epsilon = epsilon/2.
		        dtau = dtau/2.
		        B = B/2.
		        ok = .false.
	            endif
	        endif

                if (isnan(sum(rho)) .and. .not.stepok) then
		    dtau = dtaux
		    B = Bx
		    rho = rhox
		    epsum = epsumx
		    it = itx
                    diverged = .true.
                    finished = .true.
                endif

! Write out the image cube of rho.
	        t = epsum/rchisq0
                if (it >= maxiterations .or. t >= tmax .or. converged &
                    .or. diverged) then
                    rho = rhosave
                    Brho = matmul_omp(B,rho)
                    if (ihp==1) call hpcorr(Brho,wldata,mskdata,nobs)
	            rchisq = sum((data - Brho)**2/signu2)/nobs
                    call rchisqindcalc(Brho,data,signu2,mskdata,emsk,wldata, &
                        nobs,wavelengths,rchisqind,nobsind,nbands)

                    rhogrid = 0.
	            n = 0
                    do kb = 1,nbeta
	            do kt = 1,Nt
                    do np = 1,npoints
 	                n = n + 1
	                i = mod(support(np),nx)
	                j = support(np)/nx
	                if (i /= iglo .and. i /= ighi .and.  &
                            j /= jglo .and. j /= jghi)       & 
                            rhogrid(i+1,j+1,kt,kb) = rho(n) 
                        if (mod(n-1,(npoints*Nt*nbeta)/5) == 0) then
                       !    print *,' ' 
                       !    print *,'n=',n
                       !    print *,'i=',i
                       !    print *,'j=',j
                       !    print *,'kt=',kt
                       !    print *,'kb=',kb
                           print *,'Dump Buffer.'      !'rho(n)=',rho(n)
                       !    print *,'rhogrid(i+1,j+1,kt,kb)='
                       !    print *,rhogrid(i+1,j+1,kt,kb)
                       !    print *,' '
                        endif
                    enddo
	            enddo
                    enddo
                    
                   ! print *,' '
                   ! print *,'rho(42) after copy to rhogrid = ',rho(42)
                   ! print *,'rhogrid(iglo+2,jglo+2,1,1) after copy = '
                   ! print *,rhogrid(iglo+2,jglo+2,1,1)
                   ! print *,' ' 

                    do kb = 1,nbeta
	                t2m = tau2mass(kappa300,pixel,distance) &
                            * dtau*sum(kern) * 1.e-20/(mu*mH*(pixarea/Msun))
                        do kt = 1,Nt
                        do j = 1,ny
                        do i = 1,nx
                            rhogrid(i,j,kt,kb) = t2m*rhogrid(i,j,kt,kb)
                        enddo
                        enddo
                        enddo
                    enddo
                     
                   ! print *,'t2m factor = ',t2m
                   ! print *,'rhogrid(iglo+2,jglo+2,1,1) after t2m = '
                   ! print *,rhogrid(iglo+2,jglo+2,1,1)			
                   ! print *,' '
                    write(*, '("Tile portion completed on ",a8," at ",a10)') &
                        date,time
                    print *,'Terminal value of t =',t
                    print *,'Number of steps     =',it
                    print *,'N(a priori)         =',Nprior
                    print *,'N(a posteriori)     =',Nexp
                    print *,'Reduced chi squared =',rchisq
                    if (diverged) print *,'Solution divergent'

! Write output image cube.
                   ! print *,'iglo = ',iglo
                   ! print *,'ighi = ',ighi
                   ! print *,'jglo = ',jglo
                   ! print *,'jgho = ',jghi
                    ncellsx = ighi-iglo+1
                    ncellsy = jghi-jglo+1
                    allocate (outcube(ncellsx,ncellsy,Nt,nbeta))
                    outcube=rhogrid(iglo+1:ighi+1,jglo+1:jghi+1,1:Nt,1:nbeta)
                    
                   ! print *,' '
                   ! print *,'outcube(2,2,1,1) = ',outcube(2,2,1,1)
                   ! print *,' '
                    
                    call date_and_time(date,time,zone,values)
                    write(outfile,'(a,"_results/",a,"_",i5.5,"_rho.fits")') &
                        fieldname,fieldname,isf
                    outfile = removeblanks(outfile)
                    print *,'Writing out image cube '//outfile
                    call writerho(outfile,outcube,ncellsx,ncellsy,Nt,nbeta, &
                        ctype1,ctype2,crpix1,crpix2,crval1,crval2,cdelt1, &
                        cdelt2,crota2,iglo,jglo,Tgrid(1),Tgrid(Nt),eta,rhoobj,&
                        Nprior,Nexp,tdistance/1000.,it,rchisq,status)
                    deallocate (outcube)
                    write(outfile,'(a,"_results/",a,"_",i5.5,"_rchisq.txt")') &
                        fieldname,fieldname,isf
                    outfile = removeblanks(outfile)
                    print *,'Writing out reduced chi squares file '//outfile
                    open (unit=10, form='formatted', file=outfile, &
                        status='unknown')
                    write(10,'(20f8.1)') wavelengths
                    where(nobsind==0) rchisqind = 0.
                    where(rchisqind > 99999.9) rchisqind = 99999.9
                    write(10,'(20f8.1)') rchisqind
                    write(10,'(20i8)') nobsind
                    close(10)
	        endif
            enddo

            if (ok) finished = .true.
            print *,'rho(1) after writing image cube = ',rho(1)
	enddo

        deallocate (data)
        deallocate (wldata)
        deallocate (mskdata)
        deallocate (emsk)
        deallocate (signu2)
        deallocate (rchisqind)
        deallocate (nobsind)
        deallocate (B)
        deallocate (Bx)
        deallocate (support)
        deallocate (rho0)
        deallocate (rho)
        deallocate (rhosave)
        deallocate (rhox)
        deallocate (rhogrid)
        deallocate (Brho)
        deallocate (Bsum)
        deallocate (phi1)
      else
        print *,'Subfield outside coverage region'
        ncellsx = ighi-iglo+1
        ncellsy = jghi-jglo+1
        allocate (outcube(ncellsx,ncellsy,Nt,nbeta))
        outcube = 0.
        write(outfile,'(a,"_results/",a,"_",i5.5,"_rho.fits")') &
            fieldname,fieldname,isf
        outfile = removeblanks(outfile)
        print *,'Writing out image cube '//outfile(1:50)
        cdelt1 = -pixel/3600.
        cdelt2 = pixel/3600.
        rhoobj = -999.
        Nprior = -999.
        Nexp = -999.
        t = -999.
        rchisq = -999.
        call writerho(outfile,outcube,ncellsx,ncellsy,Nt,nbeta,ctype1,ctype2, &
            crpix1,crpix2,crval1,crval2,cdelt1,cdelt2,crota2,iglo,jglo, &
            Tgrid(1),Tgrid(Nt),eta,rhoobj,Nprior,Nexp,tdistance/1000.,it,&
            rchisq,status)
        deallocate (outcube)
      endif
    enddo

        deallocate (obsset)
        deallocate (obsmask)
        deallocate (sigset)
        deallocate (psfset)
        deallocate (Tgrid)
        deallocate (wavelengths)
        deallocate (pixset)
        deallocate (obsimages)
        deallocate (psfsizes)
        deallocate (kern)
        if (icc /= 0) deallocate(cctab)
endif

! If appropriate, combine all images into a mosaic.
        if (makemosaic) call ppmosaic(fieldname)
        stop

        end program ppmap
