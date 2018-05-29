        subroutine pperr(estrho,wavelengths,nbands,Tgrid,Nt,betagrid,nbeta, &
            sigbeta,cctab,pixel,kappa300,uncback,unccent)

! Error calculation for PPMAP.
!
! Input parameters:
!       estrho  =       file name for estimated image cube of expectation
!                       occupation numbers
!	wavelengths =	list of wavelengths [microns]
!       nbands  =       number of bands
!	Tgrid	=	grid of possible temperatures [K]
!       Nt      =       number of possible temperatures
!       betagrid=       grid of possible values of opacity index, beta
!       nbeta   =       number of possible beta values
!       cctab   =       colour correction table
!	pixel	=	pixel size in image grid [arcsec]
!       kappa300=       opacity at 300 microns [cm^2/g]
!
! Output parameters:
!	uncback	=	set of sky background uncertainties, one per temp.
!	unccent =	uncertainties at the central spatial location

        implicit real (a-h, o-z)
        implicit integer (i-n)

! Set some constants.
        real, parameter  :: mu = 2.8	     ! mean molecular weight
        real, parameter  :: mH = 1.6726e-24  ! mass of hydrogen atom [g]
        real, parameter  :: fwhmobj = 2.     ! FWHM of individual object [pix]
        real, parameter  :: pi = 3.141593       

        real(8),   allocatable :: gamma(:,:), gamma0(:,:), ginv(:,:)
        real(4),   allocatable :: diag(:), diag0(:)
        real(4),   allocatable :: buffer(:),rho(:,:,:,:),kern(:,:)
        real(4),   allocatable :: beamsizes(:),psfimageset(:,:,:)
        real(4),   allocatable :: sigcent(:),nyqpix(:),B(:,:)
        real(4),   allocatable :: a(:,:),ac(:,:)
        integer(4),allocatable :: counter(:,:)

        character(len=4),  allocatable :: bands(:)
        character(len=80)      :: estrho,line,imagefile,outfile,removeblanks
        character(len=40)      :: fieldname
        character(len=8)       :: ctype1, ctype2
        real(4)    cctab(Nt,nbands)
        real(4)    nprior,npost,kappa300
        real(4)    Tgrid(Nt),betagrid(nbeta)
        real(4)    wavelengths(nbands),uncback(Nt,nbeta),unccent(Nt,nbeta)
        integer    status

! Get pixel size and bands.
        call readheader4d(estrho,nxr,nyr,nz,n4,status)
        if (status > 0) then
            print *,'Could not read header for '//estrho
            stop
        endif
        if (nz /= Nt) then
            print *,'Image cube has incorrect number of temperatures'
            stop
        endif
        if (n4 /= nbeta) then
            print *,'Image cube has incorrect number of beta values'
            stop
        endif
        allocate (rho(nxr,nyr,Nt,nbeta))
        nbuffer = nxr*nyr*Nt*nbeta
        allocate (buffer(nbuffer))
        call readrho(estrho,rho,nxr,nyr,Nt,nbeta,buffer,ix0,iy0,eta,dzeta, &
            nprior,npost,nsteps,rchisq,status)
        etapost = npost/(nxr*float(nyr)*Nt*nbeta)
        allocate (bands(nbands))
        do i = 1,nbands 
            write(bands(i),'(i4.4)') nint(wavelengths(i))
        enddo
        deallocate(buffer)

! Construct "single object" kernel.
        ik = nint(1.5*fwhmobj)
        nk = 2*ik+1
        allocate (kern(nk,nk))
        aln2 = alog(2.)
        do j = 1,nk
        do i = 1,nk
            kern(i,j) = exp(-4.*(aln2/fwhmobj**2)*((i-ik-1)**2 + (j-ik-1)**2))
        enddo
        enddo
        kern = kern/sum(kern)
        rsskern = sqrt(sum(kern**2))

! Read in the PSFs.
        allocate (beamsizes(nbands))
        do i = 1,nbands
            imagefile = 'psf_'//bands(i)//'.fits'
            call readheader(imagefile,mxp,myp,crpix1,crpix2,cdelt2,status)
            if (status > 0) then
                print *,'Could not read header for '//imagefile
                stop
            endif
	    if (mxp /= myp) then 
	        print *,'PSF images must be square'
	        stop
	    endif
            if (i==1) then
                nbuffer = mxp*myp
                allocate (buffer(nbuffer))
                allocate (a(mxp,myp))
                allocate (ac(mxp,myp))
                allocate (counter(mxp,myp))
	        ncellsp = mxp
	        icentp = nint(crpix2 - 1.)
	        allocate (psfimageset(nbands,ncellsp,ncellsp))
            endif
	    if (mxp /= ncellsp) then 
	        print *,'PSF images must all be the same size'
	        stop
            endif 
            call readimage_basic(imagefile,a,mxp,myp,buffer,status)
            counter = 0
            acut = maxval(a)/2.
            where(a >= acut) counter = 1
            nhi = sum(counter)
            pixp = cdelt2*3600.
	    beamsizes(i) = 2.*sqrt(float(nhi)/pi)*pixp

! Convolve with single-object profile.
            call convol(a,kern,ac,mxp,myp,nk)
            asum = sum(ac)
            do n = 1,ncellsp
            do m = 1,ncellsp
	        psfimageset(i,m,n) = ac(m,n)/asum
            enddo
            enddo
        enddo
        deallocate(kern)
        deallocate(buffer)
        deallocate(a)
        deallocate(ac)
        deallocate(counter)

! Read measurement noise values.
        k = index(estrho,'_')
        if (k==0) then
	    print *,'Obs. images not found, or file names have wrong format'
	    stop
        endif
        fieldname = estrho(1:k-1)
        allocate (sigcent(nbands))
        allocate (nyqpix(nbands))
        do i = 1,nbands
            imagefile = removeblanks(fieldname//'_'//bands(i)//'.fits')
            call readheader(imagefile,nximage,nyimage,crpix1,crpix2,cdelt2, &
                status)
            if (status > 0) then
                print *,'Could not read header for '//imagefile
                stop
            endif
            if (i==1) then
                allocate (buffer(nximage*nyimage))
                allocate (a(nximage,nyimage))
            endif
            call readimage_wcs(imagefile,a,nximage,nyimage,buffer, &
                ctype1,ctype2,crpix1,crpix2,crval1,crval2,cdelt1,cdelt2, &
                crota2,pix,wl,sigback,status)
	    sigcent(i) = sigback
	    nyqpix(i) = cdelt2*3600.
        enddo
        deallocate(buffer)
        deallocate(a)

! Calculate system matrix.
        span = min(maxval(beamsizes), 5.*minval(beamsizes)) 
        icent = max(nint(0.5*span/pixel), 1)
        ncells = 2*icent + 1
        nstate = ncells**2 * Nt*nbeta
        nobs = 0
        do iband = 1,nbands
	    isize = max(nint(0.5*span/nyqpix(iband)), 1)
	    nsize = 2*isize + 1
	    nobs = nobs + nsize**2
            do n = 1,ncellsp
            do m = 1,ncellsp
                psfimageset(iband,m,n) = &
                    psfimageset(iband,m,n)*(nyqpix(iband)/pixp)**2
            enddo
            enddo
        enddo

        allocate (B(nobs,nstate))
        k = 0
        do ib = 1,nbeta
        do it = 1,Nt
          pfloor = 1.e-5*planckfn(5000./Tgrid(it),Tgrid(it))
          do n = 1,ncells
          do m = 1,ncells
	    k = k+1
	    i = 0
	    do iband = 1,nbands
	        bnu = max(planckfn(wavelengths(iband),Tgrid(it)), pfloor) &
                    /cctab(it,iband)
	        isize = max(nint(0.5*span/nyqpix(iband)), 1)
	        nsize = 2*isize + 1
	        do jj = 1,nsize
                do ii = 1,nsize
		    i = i+1
		    xobs = (ii-isize-1)*nyqpix(iband)
		    yobs = (jj-isize-1)*nyqpix(iband)
		    x = (m-icent-1)*pixel
		    y = (n-icent-1)*pixel
		    ip = icentp + nint((xobs - x)/pixp)
		    jp = icentp + nint((yobs - y)/pixp)
		    B(i,k) = psfimageset(iband,ip,jp)*bnu &
		        * kappa300*(300./wavelengths(iband))**betagrid(ib) &
                        * mu * (1.e20*mH)*nyqpix(iband)**2 / sigcent(iband)
	        enddo
	        enddo
	    enddo
          enddo
          enddo
        enddo
        enddo

        allocate (gamma(nstate,nstate))
        allocate (gamma0(nstate,nstate))
        allocate (ginv(nstate,nstate))
        allocate (diag(nstate))
        allocate (diag0(nstate))

        do j = 1,nstate
        do jj = 1,nstate
            bsum = 0.
            do i = 1,nobs
                bsum = bsum + B(i,j)*B(i,jj)
            enddo
            gamma(j,jj) = bsum
        enddo
        enddo
        deallocate(B)

        gamma(4,4) = gamma(4,4) + 1./sigbeta**2
        etaprime = eta*dzeta**2
        etapostprime = etapost*dzeta**2
        do n = 1,nstate
            gamma0(n,n) = gamma(n,n) + 1./etaprime
            gamma(n,n) = gamma(n,n) + 1./etapostprime
            diag0(n) = gamma0(n,n)
            diag(n) = gamma(n,n)
        enddo
        call inversep(gamma0,ginv,nstate)
        gamma0 = ginv
        call inversep(gamma,ginv,nstate)
        gamma = ginv
        deallocate(ginv)

        k = 0
        do ib = 1,nbeta
        do it = 1,Nt
        do n = 1,ncells
        do m = 1,ncells
	    k = k+1
	    if (m==icent+1 .and. n==icent+1) then 
	        if (gamma0(k,k) > 0. .and. gamma0(k,k) < 100./diag0(k)) then
                    uncback(it,ib) = max(sqrt(gamma0(k,k)),sqrt(1./diag0(k)))/2.
	        else 
                    uncback(it,ib) = sqrt(1./diag0(k))/2.
                endif
	        if (gamma(k,k) > 0. .and. gamma(k,k) < 100./diag(k)) then
                    unccent(it,ib) = max(sqrt(gamma(k,k)), sqrt(1./diag(k)))
	        else 
                    unccent(it,ib) = sqrt(1./diag(k))
                endif
                uncback(it,ib) = min(uncback(it,ib),unccent(it,ib))
	    endif
        enddo
        enddo
        enddo
        enddo
        where(uncback > unccent) uncback = unccent
        uncback = uncback*rsskern
        unccent = unccent*rsskern

        deallocate(gamma0)
        deallocate(gamma)
        deallocate(diag0)
        deallocate(diag)
        deallocate(rho)
        deallocate(bands)
        deallocate(beamsizes)
        deallocate(psfimageset)
        deallocate(sigcent)
        deallocate(nyqpix)
        return

        end subroutine pperr
