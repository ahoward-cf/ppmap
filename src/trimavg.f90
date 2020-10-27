    subroutine trimavg(x,N,xmean,sigx,trim)

! Trimmed average.

    implicit real (a-h,o-z)
    implicit integer (i,n)
    real(4),allocatable :: wt(:)
    real(4) x(N)
    logical notfinished

    notfinished = .true.
    allocate (wt(N))
    wt = 1.

    do while (notfinished)
        xsum = sum(wt*x)
        xsumsq = sum(wt*x*x)
        twt = sum(wt)
	xmean = xsum/twt
	sigx = sqrt((1./(twt-1.))*(xsumsq - xsum*xsum/twt))
	where(abs(x-xmean) > trim*sigx) wt = 0.
        twtnew = sum(wt)
	if (twtnew==twt .or. twtnew <= 6) notfinished = .false.
    enddo

    deallocate(wt)

    return

    end subroutine trimavg
