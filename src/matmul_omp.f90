! use OpenMP to parallelise matmul
module m_matmul_omp

   use iso_c_binding
   use omp_lib
   
   implicit none
   
   private
   public :: matmul_omp
   
   interface matmul_omp
      module procedure double_matmul_omp
      module procedure float_matmul_omp
      module procedure double_matvecmul_omp
      module procedure float_matvecmul_omp
      module procedure double_vecmatmul_omp
      module procedure float_vecmatmul_omp
   end interface
   
   contains
   
   ! matrix-matrix multiplication: C=AB
   function double_matmul_omp(a,b) result (c)
   
      ! argument declarations
      real(kind=c_double),intent(in) :: a(:,:)        ! matrix A
      real(kind=c_double),intent(in) :: b(:,:)        ! matrix B
      
      ! result declaration
      real(kind=c_double) :: c(size(a,1),size(b,2))   ! matrix C=AB

      ! variable declarations
      integer(kind=c_int) :: i_min                    ! lower bound
      integer(kind=c_int) :: i_max                    ! upper bound
      integer(kind=c_int) :: n_thread                 ! number of omp threads available
      integer(kind=c_int) :: thread_id                ! thread id
      integer(kind=c_int) :: i_chunk                  ! number of elements per chunk
      
      n_thread=min(omp_get_num_procs(),size(c,2))
      
      ! set chunk size
      i_chunk=ceiling(real(size(c,2))/real(n_thread))

      ! open a group of threads
      !$omp parallel num_threads(n_thread)&
      !$omp& default(private) shared(a,b,c,i_chunk)
      
         ! get thread ID
         thread_id=omp_get_thread_num()
         
         ! set array indicies
         i_min=thread_id*i_chunk+1
         i_max=min((thread_id+1)*i_chunk,size(c,2))
         
         ! perform matrix multiplication
         c(:,i_min:i_max)=&
            &matmul(a,b(:,i_min:i_max))
         
      !$omp end parallel
      
      return
   
   end function
   
   ! matrix-matrix multiplication: C=AB
   function float_matmul_omp(a,b) result (c)
   
      ! argument declarations
      real(kind=c_float),intent(in) :: a(:,:)         ! matrix A
      real(kind=c_float),intent(in) :: b(:,:)         ! matrix B
      
      ! result declaration
      real(kind=c_float) :: c(size(a,1),size(b,2))    ! matrix C=AB

      ! variable declarations
      integer(kind=c_int) :: i_min                    ! lower bound
      integer(kind=c_int) :: i_max                    ! upper bound
      integer(kind=c_int) :: n_thread                 ! number of omp threads available
      integer(kind=c_int) :: thread_id                ! thread id
      integer(kind=c_int) :: i_chunk                  ! number of elements per chunk
      
      n_thread=min(omp_get_num_procs(),size(c,2))
      
      ! set chunk size
      i_chunk=ceiling(real(size(c,2))/real(n_thread))

      ! open a group of threads
      !$omp parallel num_threads(n_thread)&
      !$omp& default(private) shared(a,b,c,i_chunk)
      
         ! get thread ID
         thread_id=omp_get_thread_num()
         
         ! set array indicies
         i_min=thread_id*i_chunk+1
         i_max=min((thread_id+1)*i_chunk,size(c,2))
         
         ! perform matrix multiplication
         c(:,i_min:i_max)=&
            &matmul(a,b(:,i_min:i_max))
         
      !$omp end parallel
      
      return
   
   end function
   
   ! matrix-vector multiplication: y=Ax
   function double_matvecmul_omp(a,x) result (y)
   
      ! argument declarations
      real(kind=c_double),intent(in) :: a(:,:)        ! matrix A
      real(kind=c_double),intent(in) :: x(:)          ! column vector x
      
      ! result declaration
      real(kind=c_double) :: y(size(a,1))             ! column vector y=Ax

      ! variable declarations
      integer(kind=c_int) :: i_min                    ! lower bound
      integer(kind=c_int) :: i_max                    ! upper bound
      integer(kind=c_int) :: n_thread                 ! number of omp threads available
      integer(kind=c_int) :: thread_id                ! thread id
      integer(kind=c_int) :: i_chunk                  ! number of elements per chunk
      
      n_thread=min(omp_get_num_procs(),size(y))
      
      ! set chunk size
      i_chunk=ceiling(real(size(y))/real(n_thread))
      

      ! open a group of threads
      !$omp parallel num_threads(n_thread)&
      !$omp& default(private) shared(a,x,y,i_chunk)
      
         ! get thread ID
         thread_id=omp_get_thread_num()
         
         i_min=thread_id*i_chunk+1
         i_max=min((thread_id+1)*i_chunk,size(y))
         
         ! perform matrix multiplication
         y(i_min:i_max)=matmul(a(i_min:i_max,:),x)
         
      !$omp end parallel
      
      return
   
   end function
   
   ! matrix-vector multiplication: y=Ax
   function float_matvecmul_omp(a,x) result (y)
   
      ! argument declarations
      real(kind=c_float),intent(in) :: a(:,:)         ! matrix A
      real(kind=c_float),intent(in) :: x(:)           ! column vector x
      
      ! result declaration
      real(kind=c_float) :: y(size(a,1))              ! column vector y=Ax

      ! variable declarations
      integer(kind=c_int) :: i_min                    ! lower bound
      integer(kind=c_int) :: i_max                    ! upper bound
      integer(kind=c_int) :: n_thread                 ! number of omp threads available
      integer(kind=c_int) :: thread_id                ! thread id
      integer(kind=c_int) :: i_chunk                  ! number of elements per chunk
      
      n_thread=min(omp_get_num_procs(),size(y))
      
      ! set chunk size
      i_chunk=ceiling(real(size(y))/real(n_thread))
      

      ! open a group of threads
      !$omp parallel num_threads(n_thread)&
      !$omp& default(private) shared(a,x,y,i_chunk)
      
         ! get thread ID
         thread_id=omp_get_thread_num()
         
         i_min=thread_id*i_chunk+1
         i_max=min((thread_id+1)*i_chunk,size(y))
         
         ! perform matrix multiplication
         y(i_min:i_max)=matmul(a(i_min:i_max,:),x)
         
      !$omp end parallel
      
      return
   
   end function
   
   ! vector-matrix multiplication: y=xB
   function double_vecmatmul_omp(x,b) result (y)
   
      ! argument declarations
      real(kind=c_double),intent(in) :: x(:)          ! row vector x
      real(kind=c_double),intent(in) :: b(:,:)        ! matrix B
      
      ! result declaration
      real(kind=c_double) :: y(size(b,2))             ! row vector y=xB

      ! variable declarations
      integer(kind=c_int) :: i_min                    ! lower bound
      integer(kind=c_int) :: i_max                    ! upper bound
      integer(kind=c_int) :: n_thread                 ! number of omp threads available
      integer(kind=c_int) :: thread_id                ! thread id
      integer(kind=c_int) :: i_chunk                  ! number of elements per chunk
      
      n_thread=min(omp_get_num_procs(),size(y))
      
      ! set chunk size
      i_chunk=ceiling(real(size(y))/real(n_thread))
      

      ! open a group of threads
      !$omp parallel num_threads(n_thread)&
      !$omp& default(private) shared(b,x,y,i_chunk)
      
         ! get thread ID
         thread_id=omp_get_thread_num()
         
         i_min=thread_id*i_chunk+1
         i_max=min((thread_id+1)*i_chunk,size(y))
         
         ! perform matrix multiplication
         y(i_min:i_max)=matmul(x,b(:,i_min:i_max))
         
      !$omp end parallel
      
      return
   
   end function
   
   ! vector-matrix multiplication: y=xB
   function float_vecmatmul_omp(x,b) result (y)
   
      ! argument declarations
      real(kind=c_float),intent(in) :: x(:)           ! row vector x
      real(kind=c_float),intent(in) :: b(:,:)         ! matrix B
      
      ! result declaration
      real(kind=c_float) :: y(size(b,2))              ! row vector y=xB

      ! variable declarations
      integer(kind=c_int) :: i_min                    ! lower bound
      integer(kind=c_int) :: i_max                    ! upper bound
      integer(kind=c_int) :: n_thread                 ! number of omp threads available
      integer(kind=c_int) :: thread_id                ! thread id
      integer(kind=c_int) :: i_chunk                  ! number of elements per chunk
      
      n_thread=min(omp_get_num_procs(),size(y))
      
      ! set chunk size
      i_chunk=ceiling(real(size(y))/real(n_thread))
      

      ! open a group of threads
      !$omp parallel num_threads(n_thread)&
      !$omp& default(private) shared(b,x,y,i_chunk)
      
         ! get thread ID
         thread_id=omp_get_thread_num()
         
         i_min=thread_id*i_chunk+1
         i_max=min((thread_id+1)*i_chunk,size(y))
         
         ! perform matrix multiplication
         y(i_min:i_max)=matmul(x,b(:,i_min:i_max))
         
      !$omp end parallel
      
      return
   
   end function

end module