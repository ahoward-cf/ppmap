    subroutine divider(nx,ny,ncells,ncellsbest,noverlap,nnodes,nsubx,nsuby)

! Divide up nx by ny array into subfields for processing by PPMAP.

    implicit real (a-h,o-z)
    implicit integer (i,n)

    real(4) minrem
    logical nofit

    nnodes = 16			! number of nodes available simultaneously
    minrem = float(nx)*ny
    ncellsbest = ncells
    nofit = .true.

    do nc = 2*noverlap,ncells,2 
	nsubxnom = (nx - noverlap)/(nc - noverlap)
	nsubynom = (ny - noverlap)/(nc - noverlap)
	nremx = nc - (nsubxnom*(nc-noverlap) + noverlap)
	nremy = nc - (nsubynom*(nc-noverlap) + noverlap)
	rem = float(nremx)*nremy
	if (rem < minrem) then 
	    ncellsbest = nc
	    nremxbest = nremx
	    nremybest = nremy
	    nsubx = nsubxnom
	    nsuby = nsubynom
	    minrem = rem
	    nofit = .false.
	endif
    enddo

    if (nofit) then 
	interval = ncells-noverlap
	nsubx = (nx-ncells)/interval + 1
	nsuby = (ny-ncells)/interval + 1
    endif

    print *,'Imaging array divided into',nsubx,' subfields in x, and', &
	nsuby,' in y,'
    print *,'with a subfield width of',ncellsbest,' pixels'
    print *,'Actual coverage:',0,(nsubx-1)*(ncellsbest-noverlap)+ncellsbest-1,&
        0,(nsuby-1)*(ncellsbest-noverlap)+ncellsbest-1
    print *,'Overlap =',noverlap

    end subroutine divider
