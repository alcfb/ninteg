module m_size

    use m_parameters
    implicit none

    type, public :: t_size
        real(8), allocatable :: dz(:), z0(:)
    contains
        procedure :: make => p_make
        procedure :: predict => p_predict
        procedure :: correct => p_correct
        procedure :: free => p_free
    end type t_size

    private
    contains

    subroutine p_make (this, rank)
        implicit none
        integer, intent (in) :: rank
        class (t_size), intent (inout) :: this
        allocate (this % dz (rank))
        allocate (this % z0 (rank))
    end subroutine p_make

    function ratio (q, qnew, dx) result (r)
        integer, intent (in) :: q, qnew
        real(8), intent (in) :: dx
        real(8) :: r, tau, dq, C1
        integer :: C0
        if (q < qnew) then
            Dq = dx * FACT(q) / (q+2)
            C0 = q+2
            C1 = 1.4D0
        elseif (q > qnew) then
            Dq = dx * FACT(q-1)
            C0 = q
            C1 = 1.3D0
        else
            Dq = dx * FACT(q) / (q+1)
            C0 = q+1
            C1 = 1.2D0
        endif
        r = 1.D0 / (C1 * (Dq**(1.D0/C0) + 1.D-6))
    end function ratio

    subroutine p_predict (this, n, q, x, zc, e0, e1, atol, rtol, rs, ru, rd)
        ! Compute next step size
        implicit none
        class (t_size), intent (inout) :: this
        integer, intent (in) :: q, n
        real(8), intent (in) :: zc (MAX_BDF_ORDER+1, n)
        real(8), intent (in) :: e0 (n), e1 (n), x(n)
        real(8), intent (in) :: atol (n), rtol (n)
        real(8), intent (inout) :: rs, rd, ru
        real(8) :: err

        ! for order q
        call erms (n, e0, x, atol, rtol, err)
        rs = ratio (q, q, err)

        ! for order q-1
        this % dz = zc (q+1,:)
        call erms (n, this % dz, x, atol, rtol, err)
        rd = ratio (q, q-1, err)

        ! for order q+1
        this % dz = e0 - e1
        call erms (n, this % dz, x, atol, rtol, err)
        ru = ratio (q, q+1, err)

    end subroutine p_predict

    subroutine p_correct (this, h, time, tend)
        implicit none
        class (t_size), intent (inout) :: this
        real(8), intent (in) :: time, tend
        real(8), intent (inout) :: h
        integer :: n
        ! Check how far we are from the end of time
        if (time + 20 * h < tend) return
        ! Compute the number of steps by the end
        n = int ((tend - time)/ h) + 1
        ! Correct to the end time
        h = (tend - time) / n
    end subroutine p_correct

    subroutine p_free (this)
        implicit none
        class (t_size), intent (inout) :: this
        deallocate (this % dz)
        deallocate (this % z0)
    end subroutine p_free

end module m_size