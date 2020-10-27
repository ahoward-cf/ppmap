    subroutine rchisqindcalc(Brho,data,signu2,mskdata,emsk,wldata,nobs, &
        wavelengths,rchisqind,nobsind,nbands)

! Calculate reduced chi squared for each individual wavelength.

    real Brho(nobs), data(nobs), signu2(nobs), wldata(nobs), &
        wavelengths(nbands),rchisqind(nbands)
    integer mskdata(nobs),emsk(nobs),nobs,nobsind(nbands),nbands,N,k
    integer,    allocatable :: counter(:), xcount(:)
    real,       allocatable :: x(:)

    allocate (counter(nobs))
    allocate (x(nobs))
    allocate (xcount(nobs))

    do k = 1,nbands
        counter = 0
        where (Brho > 0. .and. abs(wldata - wavelengths(k)) < 1.e-4) counter = 1
        x = 0.
        xcount = 0
        where (counter==1 .and. mskdata==0 .and. emsk==1) &
            x = (data - Brho)**2 / signu2
        where (counter==1 .and. mskdata==0 .and. emsk==1) xcount = 1
        nobsind(k) = sum(xcount)
        if (nobsind(k) /= 0) rchisqind(k) = sum(x)/nobsind(k)
    enddo

    deallocate (counter)
    deallocate (x)
    deallocate (xcount)
    return

    end subroutine rchisqindcalc
