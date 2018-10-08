    function getnoise(a,nxo,nyo,fwhm)

! Estimate sky noise by iterative subtraction of smoothed background.

    implicit real (a-h,o-z)
    implicit integer (i,n)
    real(4),allocatable :: as(:,:),ac(:,:),atrim(:,:),wt(:,:),kern(:,:)
    logical(4),allocatable :: mask(:,:)
    real(4) a(nxo,nyo)
    logical converged

    isubsample = max(nint(fwhm/40.), 1)
    nx = nxo/isubsample
    ny = nyo/isubsample
    allocate (as(nx,ny))
    allocate (ac(nx,ny))
    allocate (atrim(nx,ny))
    allocate (wt(nx,ny))
    allocate (mask(nx,ny))

    if (isubsample==1) then
        as = a
    else
        do j = 1,ny
        do i = 1,nx
            as(i,j) = a(isubsample*i, isubsample*j)
        enddo
        enddo
    endif

! Construct Gaussian smoothing kernel.
    fwhms = fwhm/isubsample
    ik = nint(0.75*fwhms)
    nk = 2*ik+1
    allocate (kern(nk,nk))
    aln2 = log(2.)
    do j = 1,nk
    do i = 1,nk
        kern(i,j) = exp(-4.*(aln2/fwhms**2)*((i-ik-1)**2 + (j-ik-1)**2))
    enddo
    enddo
    kern = kern/sum(kern)

    wt = 1.
    where(as==0.) wt = 0.
    atrim = as*wt

! Begin iterative smoothing.
    converged = .false.
    sigprev = 1.e35
    maxit = 100
    tol = 0.001
    it = 0
    do while (.not.converged .and. it<maxit)
        it = it+1
        call convol(atrim,kern,ac,nx,ny,nk)
        sigma = sqrt(sum(wt*(atrim-ac)**2)/sum(wt))
        mask = .false.
        where(atrim-ac > 3.*sigma) mask = .true.
        where(mask) atrim = ac
        where(mask) wt = 0.
        if (abs(sigma-sigprev)/sigma < tol) converged = .true.
        sigprev = sigma
    enddo

    deallocate(as)
    deallocate(ac)
    deallocate(atrim)
    deallocate(wt)
    deallocate(mask)
    deallocate(kern)
    getnoise = sigma
    return

    end function getnoise
