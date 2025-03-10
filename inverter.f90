module matrix_utils
   use iso_c_binding, only: c_double
   implicit none
   integer, parameter :: dp = c_double
contains
   subroutine invert_matrix(matrix)
      implicit none
      real(dp), dimension(:,:), intent(inout) :: matrix
      integer :: N, LWORK, IERR
      integer, allocatable :: IPIV(:)
      real(dp), allocatable :: WORK(:)

      if (size(matrix,1) /= size(matrix,2)) then
         stop "Matrix is not square"
      end if
      N = size(matrix,1)
      
      allocate(IPIV(N), stat=IERR)
      if (IERR /= 0) stop "Failed to allocate IPIV"
      
      LWORK = N * N
      allocate(WORK(LWORK), stat=IERR)
      if (IERR /= 0) stop "Failed to allocate WORK"

      call dgetrf(N, N, matrix, N, IPIV, IERR)
      if (IERR /= 0) stop "Error in dgetrf: Matrix is singular"
      call dgetri(N, matrix, N, IPIV, WORK, LWORK, IERR)
      if (IERR /= 0) stop "Error in dgetri: Matrix is singular"

      deallocate(IPIV, WORK)
   end subroutine invert_matrix
end module matrix_utils
