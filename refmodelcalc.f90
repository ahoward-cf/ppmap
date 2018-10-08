        subroutine refmodelcalc(lambda,reflambda,beta,T,pixel,fwhmobj, &
            psf,ncp,refmodel,nref)

! Calculate model image at wavelength lambda based on gaussian object of 
! size 'fwhmobj' [pixels], temperature T [K], and unity central optical depth 
! at reference wavelength. 
! Output refmodel is an array nref by nref with centre pixel (iref,iref).
! Output units: Jy/pixel.

        implicit real (a-h,o-z)
        implicit integer (i-n)

        real,   intent(in),                           &
        &       dimension(ncp,ncp)      ::   psf         ! point spread func

        real,   intent(out),                        &
        &       dimension(nref,nref)    ::   refmodel    ! reference model

        real,   allocatable             ::   tau(:,:), image(:,:)
        real,   allocatable             ::   kernel(:,:)
        real lambda
        real,      parameter   :: pi      = 3.141593

        allocate (image(nref,nref))
        image = 0.
        iref = nref/2 + 1
        fwhm = fwhmobj*pixel
        coef = -4.*log(2.)*(pixel/fwhm)**2

        allocate (tau(nref,nref))
        tau = 0.
        do j = 1,nref
        do i = 1,nref
            tau(i,j) = exp(coef*(float(i-iref)**2 + float(j-iref)**2))
        enddo
        enddo
        image = tau*(reflambda/lambda)**beta * planckfn(lambda,T)*pixel**2
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
        deallocate(image)
        return

        end subroutine refmodelcalc
