module m_bdf

    use m_parameters
    use m_nonlin
    implicit none

    type, public :: t_bdf
        integer :: status
        integer, private :: rank
        real (8), allocatable, private :: zp (:,:), b (:), x (:)
        type (t_nonlin) :: nonlin
    contains
        procedure :: make => p_make
        procedure :: next => p_next
        procedure :: free => p_free
    end type t_bdf

    private
    contains

    subroutine p_make (this, rank)

        implicit none
        integer, intent (in) :: rank
        class (t_bdf), intent (inout) :: this

        this % rank = rank

        allocate (this % zp (MAX_BDF_ORDER+1, rank))
        allocate (this % x (rank))
        allocate (this % b (rank))

        call this % nonlin % make (rank)
        call this % nonlin % deft

    end subroutine p_make

    subroutine p_next (this, n, q, t, h, zc, e)
        ! Solves the system of equations
        ! Nordsieck vector: ZC = (ZC(0),...,ZC(q))^T
        ! Predictor vector: ZP = PASCAL * ZC0
        ! Residual value:   dz = ZC(0) - ZP(0)
        !
        ! Polynomial passes q nodes in [t,t+h]:
        ! ZC(1) - ZP(1) = l1 * dz               (1)
        ! ZC(q) - ZP(q) = lq * dz, q > 1
        !
        ! Polynomial's first derivative
        ! ZC(1) = h * f(ZC(0))                  (2)
        ! 
        ! Resolve 0th component from Eqs. 1 and 2:
        ! ZC(0) : h * f(ZC(0)) - l1 * ZC(0) = ZP(1) - l1 * ZP(0)
        !
        ! Resolve other components:
        ! ZC(q) = ZP(q) + lq * (ZC(0) - ZP(0)), q > 0
        !

        implicit none
        class (t_bdf), intent (inout) :: this
        real(8), intent (inout) :: zc (MAX_BDF_ORDER+1, n)
        real(8), intent (out) :: e (n)
        real(8), intent (in) :: t, h
        integer, intent (in) :: q, n
        integer :: i, j
        real(8) :: alpha, c0

        ! Predictor ZP = PASCAL * ZC
        this % zp = 0
        do i = 1, q + 1
            do j = 1, q + 1
                if (PASCAL (i,j) > 0) then
                    this % zp (i,:) = this % zp (i,:) + PASCAL (i,j) * zc (j,:)
                endif
            enddo
        enddo

        ! Right hand side B = YP - l1 * h * YP'
        this % b = this % zp (1,:) - this % zp (2,:) / CONST_L (q,1)

        ! Next time step solution, x: x - l1 * h * f(x) = b
        alpha = h / CONST_L (q,1)
        this % nonlin % max_err = 0.1D0 * (q+1)/(q+2) / FACT(q) / CONST_L (q,q) * CONST_L (q,1)
        this % x = this % zp (1,:)
        call this % nonlin % solv (n, alpha, t + h, this % b, this % x)

        ! Corrector ZC = ZP + li/l1 * (X - YP)
        zc (1,:) = this % x
        this % b = zc (1,:) - this % zp (1,:)
        do i = 2, q + 1
            zc (i,:) = this % zp (i,:) + this % b * CONST_L (q,i-1)
        enddo

        e = zc(q+1,:) - this % zp (q+1,:)

        this % status = 0

    end subroutine p_next

    subroutine p_free (this)
        implicit none
        class (t_bdf), intent (inout) :: this
        deallocate (this % x)
        deallocate (this % b)
        deallocate (this % zp)
        call this % nonlin % free
    end subroutine p_free
end module m_bdf