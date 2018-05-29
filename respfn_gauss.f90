        subroutine respfn_gauss(lambda,reflambda,beta,nx,ny,obspix,x,y,T, &
	    nmod,pixel,fwhmobj,psf,ncp,refexists,refmodel,nref,rf)

! Calculate spatial response at wavelength lambda of a gaussian object of 
! size 'fwhmobj' [pixels], temperature T [K], and unity central optical depth 
! at reference wavelength, at (x,y) [arcsec] with respect to image center.  
! Output units: Jy/pixel.

        implicit real (a-h,o-z)
        implicit integer (i-n)

        real,   intent(in),                           &
        &       dimension(ncp,ncp)      ::   psf         ! point spread func

        real,   intent(inout),                        &
        &       dimension(nref,nref)    ::   refmodel    ! reference model

        real,   intent(out),                          &
        &       dimension(nx,ny)        ::   rf          ! output image

        real,   allocatable             ::   tau(:,:), image(:,:)
        real,   allocatable             ::   kernel(:,:), refshift(:,:)
        real lambda
        integer xc,yc
        logical refexists

        allocate (image(nref,nref))
        image = 0.
        iref = nref/2 + 1
        fwhm = fwhmobj*pixel
        coef = -4.*log(2.)*(pixel/fwhm)**2

        if (.not.refexists) then 
            allocate (tau(nref,nref))
            tau = 0.
            do j = 1,nref
            do i = 1,nref
                tau(i,j) = exp(coef*(float(i-iref)**2 + float(j-iref)**2))
            enddo
            enddo
            pfloor = 1.e-5*planckfn(5000./T,T)

            image = (tau*(reflambda/lambda)**beta * &
                    max(planckfn(lambda,T), pfloor))*pixel**2
            deallocate(tau)

! Convolve with PSF.
            allocate (kernel(ncp+1, ncp+1))
            icp = ncp/2 + 1
            kernel = 0.
            do j = 1,ncp
            do i = 1,ncp
                kernel(ncp-i+2, ncp-j+2) = psf(i,j) 
            enddo
            enddo
            kernel = kernel/sum(kernel)
            refmodel = 0.

            do j = icp, nref-icp
            do i = icp, nref-icp
                csum = 0.
                do jj = 1,ncp+1
                do ii = 1,ncp+1
                    csum = csum + kernel(ii,jj)*image(i+ii-icp, j+jj-icp)
                enddo
                enddo
                refmodel(i,j) = csum
            enddo
            enddo
            deallocate(kernel)
            refexists = .true.
        endif

! Refmodel is an array nref by nref with centre pixel (iref,iref).
! We need to interpolate at intervals of obspix/pixel (= 1/mag) onto
! an array nx by ny which is centred on (iref,iref).

! Shift the array by xc,yc.
        allocate (refshift(nref,nref))
        xc = nint(x/pixel)
        yc = nint(y/pixel)
        ilo = max(xc+1, 1)
        jlo = max(yc+1, 1)
        ihi = min(xc+nref, nref)
        jhi = min(yc+nref, nref)
        do j = jlo,jhi
        do i = ilo,ihi
            refshift(i-xc, j-yc) = refmodel(i,j)
        enddo
        enddo

! Resample to Nyquist-sampled pixel size.
        rmag = pixel/obspix
        call resample(refshift,nref,nref,image,nref,nref,rmag)
        deallocate(refshift)
        image = image/rmag**2

        if (nref >= nx .and. nref >= ny) then 
            rf = image(iref-nx/2:iref-nx/2+nx-1, iref-ny/2:iref-ny/2+ny-1)
        else 
            rf = 0.
            ilo = max(2-iref+nx/2, 1)
            jlo = max(2-iref+ny/2, 1)
            ihi = min(nref-iref+nx/2+1, nref)
            do j = 1,nref
            do i = 1,nref
                rf(i,j) = image(iref-nx/2+i-1, iref-ny/2+j-1)
            enddo
            enddo
        endif

        deallocate(image)
        return

        end subroutine respfn_gauss
