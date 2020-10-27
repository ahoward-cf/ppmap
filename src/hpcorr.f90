    subroutine hpcorr(Brho,wldata,mskdata,nobs)

! Subtract mean value from model and observed images at wavelengths for which 
! the former had been subject to high pass filtering during the data reduction 
! stage. This normally means data from ground-based observatories. The ones
! currently allowed for are:
!
!       SABOCA 350 microns (designated as 351 to distinguish from SPIRE 350)
!       SCUBA2 450 microns
!       SCUBA2 850 microns
!       LABOCA 870 microns
!       Nika2 1150 microns
!       ALMA  1300 microns
!       Nika2 2000 microns
!       PdBI  3000 microns
!
! The subtraction is applied to the vector Brho, which contains all of the
! synthesized images for the current set of model parameters. It also gets
! applied to a data vector of the same dimensionality.

    real Brho(nobs), wldata(nobs), wlo(8), whi(8), bmean
    integer mskdata(nobs),nobs,N,nranges,k
    integer,    allocatable :: counter(:), xcount(:)
    real,       allocatable :: x(:)

    data nranges /8/                    ! number of separate high pass ranges
    data wlo/350.1, 440., 800., 860., 900., 1200., 1500., 2600./ 
                                        ! lower limits of wavelength ranges
    data whi/352,   460., 860., 900.,1200., 1400., 2500., 3500./ 
                                        ! upper limits of wavelength ranges

    allocate (counter(nobs))
    allocate (x(nobs))
    allocate (xcount(nobs))

    do k = 1,nranges
        counter = 0
        where (Brho > 0. .and. wldata > wlo(k) .and. wldata < whi(k)) counter=1
        x = 0.
        xcount = 0
        bmean = 0.
        where (counter==1 .and. mskdata==0) x = Brho
        where (counter==1 .and. mskdata==0) xcount = 1
        N = sum(xcount)
        if (N /= 0) bmean = sum(x)/N
        where (counter==1) Brho = Brho - bmean
    enddo

    deallocate (counter)
    deallocate (x)
    deallocate (xcount)
    return

    end subroutine hpcorr
