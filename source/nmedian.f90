!*******************************************************
!* Given an array A of N numbers, returns their median *
!* value XMED. In this version, the input array A is   *
!* not destroyed.                                       *
!*******************************************************
SUBROUTINE NMEDIAN(a,N,XMED)
  real a(N)
  real(4),   allocatable :: X(:)
  allocate (x(N))
  X = a
  call hpsort(N,X)
  N2=N/2
  if (2*N2.eq.N) then
    XMED = 0.5*(X(N2)+X(N2+1))
  else
    XMED = X(N2+1)
  endif
  deallocate(X)
  return
END


!*****************************************************
!*  Sorts an array RA of length N in ascending order *
!*                by the Heapsort method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*          N     size of table RA                       *
!*      RA        table to be sorted                     *
!* OUTPUT:                                           *
!*          RA    table sorted in ascending order        *
!*                                                   *
!* NOTE: The Heapsort method is a N Log2 N routine,  *
!*       and can be used for very large arrays.      *
!*****************************************************         
SUBROUTINE HPSORT(N,RA)
  real RA(N)
  L=N/2+1
  IR=N
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
    L=L-1
    RRA=RA(L)
  else
    RRA=RA(IR)
    RA(IR)=RA(1)
    IR=IR-1
    if(IR.eq.1)then
      RA(1)=RRA
      return
    end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
  if(J < IR)then
    if(RA(J) < RA(J+1))  J=J+1
  end if
  if(RRA < RA(J))then
    RA(I)=RA(J)
    I=J; J=J+J
  else
    J=IR+1
  end if
  goto 20
  end if
  RA(I)=RRA
  goto 10
END
