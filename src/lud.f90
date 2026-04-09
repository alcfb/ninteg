module m_lud
    implicit none

    type, public :: t_lud
        real(8), allocatable :: A (:,:)
        integer :: rank
    contains
        procedure :: make => p_make
        procedure :: solv => p_solv
        procedure :: free => p_free
    end type t_lud

    private
    contains

    subroutine p_make (self, n)
        class (t_lud), intent (inout) :: self
        integer, intent (in) :: n
        allocate (self % A (n,n))
    end subroutine p_make

    subroutine p_solv (self, mat, x, n)
        class (t_lud), intent (inout) :: self
        integer, intent (in) :: n
        real(8), intent (in) :: mat (n,n)
        real(8), intent(out) :: x (n)
        real(8) :: w0 = 1.0001D0
        integer :: i, j

        x = 1.D0

        self % A = mat

        do i = 1, n
            if (self % A (i,i) == 0) then
                self % A (i,i) = 1.D0
            else
                self % A (i,i) = self % A (i,i) * w0
            endif
        enddo

        call luSolve (n, self % A, x)

    end subroutine p_solv

    subroutine p_free (self)
        class (t_lud), intent (inout) :: self
        deallocate (self % A)
    end subroutine p_free

end module m_lud