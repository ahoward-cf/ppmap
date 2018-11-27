    subroutine divider(nx,ny,ncells,ncellsbest,noverlap,nnodes, &
                nsubx,nsuby,nsubxmid,nsubymid,nxmid,nymid, &
                nsubxmidlo, nsubymidlo, ilostart, jlostart)

! Divide up nx by ny array into subfields for processing by PPMAP.

    implicit real (a-h,o-z)
    implicit integer (i,n)

    real(4) minrem
    logical nofit
    integer :: xstartprint, ystartprint

    nnodes = 16			! number of nodes available simultaneously
    minrem = float(nx)*ny
    ncellsbest = ncells
    nofit = .true.
    
    xstartprint = 0
    ystartprint = 0
    
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
    
    if (mod(nsubx,2) == 0) then
    nsubx = nsubx + 1
    endif
    if (mod(nsuby,2) == 0) then
    nsuby = nsuby + 1
    endif
    
    nsubxmid = ceiling(float(nsubx) / 2)
    nsubymid = ceiling(float(nsuby) / 2)
    
    nxmid = ceiling(float(nx) / 2)
    nymid = ceiling(float(ny) / 2)
    
    nsubxmidlo = nxmid - floor(float(ncellsbest) / 2)
    nsubymidlo = nymid - floor(float(ncellsbest) / 2)
    
    ilostart = nsubxmidlo - (floor(float(nsubx)/2) &
    * (ncellsbest - noverlap)) + 1
    jlostart = nsubymidlo - (floor(float(nsuby)/2) &
    * (ncellsbest - noverlap)) + 1
    
    if (ilostart < 0) then
        ilostart = ilostart + ncellsbest
        nsubx = nsubx - 2
    endif
    if (jlostart < 0) then
        jlostart = jlostart + ncellsbest
        nsuby = nsuby - 2
    endif
    

    print *,'Imaging array divided into',nsubx,' subfields in x, and', &
	nsuby,' in y,'
    print *,'with a subfield width of',ncellsbest,' pixels'
    print *,'Actual coverage:',ilostart,&
    (nsubx-1)*(ncellsbest-noverlap)+ncellsbest-1 + ilostart,&
    jlostart,(nsuby-1)*(ncellsbest-noverlap)+ncellsbest-1 + jlostart
    print *,'Overlap =',noverlap

    end subroutine divider
