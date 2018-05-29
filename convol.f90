   subroutine convol(a, b, c, idima, jdima, idimb)

!-------------------------------------------------------------
! Convolve array A with kernel B to yield output array C.
! Input matrix A is dimensioned idima x jdima
! Input matrix B is dimensioned idimb x idimb (idimb must be an odd number)
! Output matrix C is dimensioned idima x jdima
!-------------------------------------------------------------
    implicit real(a-h,o-z)
    implicit integer(i-n)
  
    real(4) a(idima,jdima)    ! matrix a
    real(4) b(idimb,idimb)    ! matrix b
    real(4) c(idima,jdima)    ! a convolved with b
   
    ihb = idimb/2
    c = sum(a)/(idima*jdima)
   
    do j = 1+ihb, jdima-ihb
    do i = 1+ihb, idima-ihb
        c(i,j) = sum(a(i-ihb:i+ihb, j-ihb:j+ihb)*b)
    enddo
    enddo

    return
   
end subroutine convol
