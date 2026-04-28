module m_bdf

    use m_parameters
    use m_nonlin
    implicit none

    type, public :: t_bdf
        integer :: status
        integer, private :: rank
        integer, allocatable, private :: e2 (:)
        real(8), allocatable, private :: PL (:,:), r0(:), r1(:)
        real(8), allocatable, private :: zp (:,:), b (:), x (:)
        type (t_nonlin) :: nonlin
    contains
        procedure :: make => p_make
        procedure :: next => p_next
        procedure :: free => p_free
        procedure :: interpolate => p_interpolate
        procedure :: rescale => p_rescale
    end type t_bdf

    private
    contains

    subroutine p_make (this, rank)

        implicit none
        integer, intent (in) :: rank
        class (t_bdf), intent (inout) :: this
        integer :: i, j

        this % rank = rank

        allocate (this % PL (MAX_BDF_ORDER, MAX_BDF_ORDER+1))
        allocate (this % zp (MAX_BDF_ORDER+1, rank))
        allocate (this % x (rank))
        allocate (this % b (rank))
        allocate (this % e2 (MAX_BDF_ORDER+1))
        allocate (this % r0 (MAX_BDF_ORDER))
        allocate (this % r1 (MAX_BDF_ORDER))

        call this % nonlin % make (rank)
        call this % nonlin % deft

        do i = 0, MAX_BDF_ORDER
            this % e2 (i+1) = i
        enddo

        do i = 1, MAX_BDF_ORDER
            do j = 1, MAX_BDF_ORDER + 1
                this % PL (i,j) = sum (PASCAL (j,:) * LL (i,:))
            enddo
        enddo

        do i = 1, MAX_BDF_ORDER
            this % r0 (i) = sum (LL(i,:))
            this % r1 (i) = sum (LL(i,:) * this % e2)
        enddo

    end subroutine p_make



    subroutine p_next (this, n, q, t, h, zc, e)
        implicit none
        class (t_bdf), intent (inout) :: this
        real(8), intent (inout) :: zc (MAX_BDF_ORDER+1, n)
        real(8), intent (out) :: e (n)
        real(8), intent (in) :: t, h
        integer, intent (in) :: q, n
        integer :: i, j
        real(8) :: alpha, c0

        ! ------ ZP = PASCAL * ZC' ------

        this % zp = 0
        do i = 1, q + 1
            do j = 1, q + 1
                if (PASCAL (i,j) > 0) then
                    this % zp (i,:) = this % zp (i,:) + PASCAL (i,j) * zc (j,:)
                endif
            enddo
        enddo

        ! ------ RHS = ZP(1) - alpha * ZP(2) ------

        alpha = this % r0 (q) / this % r1 (q)

        this % b = this % zp (1,:) - this % zp (2,:) * alpha

        ! ------ In-Step Equation Solution Error ------

        this % nonlin % max_err = 0.1D0 * (q+1) / (q+2) / FACT(q) * this % PL (q,2)

        ! ------ X - alpha * h * F (t + h, X) = RHS => ZC(1) = X ------

        this % x = this % zp (1,:)

        call this % nonlin % solv (n, alpha * h, t + h, this % b, this % x)

        zc (1,:) = this % x

        ! ------ BETA = (ZC(1) - ZP(1)) / sum(L) ------

        this % b = (zc (1,:) - this % zp (1,:)) / this % r0 (q)

        ! ------ ZC(i) = ZP(i) + BETA * PASCAL * L ------

        do i = 2, q + 1
            zc (i,:) = this % zp (i,:) + this % b * this % PL (q,i)
        enddo

        ! ------ Local Truncation Error ------

        e = zc(q+1,:) - this % zp (q+1,:)

        this % status = 0

    end subroutine p_next


    subroutine p_interpolate (this, n, q, zc)

        implicit none
        class (t_bdf), intent (inout) :: this
        real(8), intent (inout) :: zc (MAX_BDF_ORDER+1, n)
        integer, intent (in) :: q, n
        integer :: i, j

        do i = 1, q + 1
            do j = i + 1, q + 1
                zc (i,:) = zc (i,:) + PASCAL (i,j) * zc (j,:)
            enddo
        enddo

    end subroutine p_interpolate



    subroutine p_rescale (this, n, z, h0, h1)
        implicit none
        class (t_bdf), intent (inout) :: this
        integer, intent (in) :: n
        real(8), intent (inout) :: z (MAX_BDF_ORDER+1, n)
        real(8), intent (in) :: h0, h1
        integer :: i
        real(8) :: ratio
        ratio = h1 / h0
        do i = 2, MAX_BDF_ORDER+1
            z (i,:) = ratio**(i-1) * z(i,:)
        enddo
    end subroutine p_rescale


    subroutine p_free (this)
        implicit none
        class (t_bdf), intent (inout) :: this
        deallocate (this % x)
        deallocate (this % b)
        deallocate (this % zp)
        deallocate (this % r0)
        deallocate (this % r1)
        deallocate (this % e2)
        deallocate (this % PL)
        call this % nonlin % free
    end subroutine p_free
end module m_bdf