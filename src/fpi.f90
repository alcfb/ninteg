module m_fpi

    use m_lud, only : t_lud

    implicit none

    type, public :: t_fpi
        real(8) :: alpha, err, max_err
        real(8), allocatable :: e (:,:), x (:,:), b (:), mat (:,:), w(:)
        real(8), allocatable :: atol (:), rtol (:)
        integer :: max_order = 5, imethod
        integer :: rank, i_iter, status=-10
        integer :: i_order, i_shift, max_iter
        type (t_lud) :: lud
    contains
        procedure :: make => p_make
        procedure :: default => p_default
        procedure :: init => p_init
        procedure :: update => p_update
        procedure :: free => p_free
    end type t_fpi

    private
    contains

    subroutine p_make (this, rank, order)
        class (t_fpi), intent (inout) :: this
        integer, intent (in) :: rank, order
        this % rank = rank
        this % max_order = order
        allocate (this % atol (rank))
        allocate (this % rtol (rank))
        allocate (this % mat (order, order))
        allocate (this % e (rank, order))
        allocate (this % x (rank, order))
        allocate (this % b (order))
        allocate (this % w (rank))
        call this % lud % make (order)
        this % status = -11
    end subroutine p_make

    subroutine p_default (this)
        class (t_fpi), intent (inout) :: this
        this % alpha = 0.5D0
        this % rtol = 1.D-6
        this % atol = 1.D-12
        this % max_err = 1.D0
        this % w = 0.01D0
        this % max_iter = 100
        this % status = -12
        this % imethod = 2
    end subroutine p_default

    subroutine p_init (this)
        class (t_fpi), intent (inout) :: this
        this % i_iter = 0
        this % i_order = 0
        this % i_shift = 0
        this % x = 0
        this % e = 0
        this % b = 0
        this % err = 1.D+10
        this % status = 1
    end subroutine p_init

    subroutine p_update (this, x, rank)
        class (t_fpi), intent (inout) :: this
        integer, intent (in) :: rank
        real(8), intent (inout) :: x (rank)
        if (this % imethod == 1) call update_picard (this, x, rank)
        if (this % imethod == 2) call update_anderson (this, x, rank)
    end subroutine p_update

    subroutine update_picard (this, x, rank)
        class (t_fpi), intent (inout) :: this
        integer, intent (in) :: rank
        real(8), intent (inout) :: x (rank)

        if (this % i_iter == 0) then
            this % x (:,1) = x
            this % err = 1.D+10
        endif

        if (this % i_iter >= 1) then
            this % e (:,1) = x - this % x (:,1)
            call erms (rank, this % e (:,1), this % x (:,1), this % atol, this % rtol, this % err)
            x = this % x (:,1) + this % alpha * this % e(:,1)
            this % x (:,1) = x
        endif

        this % i_iter = this % i_iter + 1
        if (this % err < this % max_err) this % status = 0
        if (this % i_iter > this % max_iter) this % status = -1

    end subroutine update_picard

    subroutine update_anderson (this, x, rank)
        class (t_fpi), intent (inout) :: this
        integer, intent (in) :: rank
        real(8), intent (inout) :: x (rank)
        real(8) :: total, x_norm, e_norm
        integer :: i, j, k, n, m, info

        if (this % i_iter == 0) then

            this % x (:,1) = x

            this % err = 1.D+10

        endif

        if (this % i_iter == 1) then

            this % e (:,1) = x - this % x (:,1)

            call erms (rank, this % e (:,1), this % x (:,1), this % atol, this % rtol, this % err)

            this % x (:,2) = x

            this % i_shift = mod (this % i_shift, this % max_order) + 1

            this % i_order = min (this % i_order + 1, this % max_order)

        endif

        if (this % i_iter > 1) then

            m = mod (this % i_shift, this % max_order) + 1

            n = min (this % i_order + 1, this % max_order)

            this % e (:,m) = x - this % x (:,m)

            this % w = 1.D0 / (this % atol + this % rtol * this % x (:,m))

            call erms (rank, this % e (:,m), this % x (:,m), this % atol, this % rtol, this % err)

            if (this % err > this % max_err) then

                this % mat = 0.D0
                do i = 1, this % max_order
                    do j = 1, this % max_order
                        this % mat (i,j) = sum (this % e (:,i) * this % e (:,j) * this % w**2)
                    enddo
                enddo

                call this % lud % solv (this % mat, this % b, this % max_order)

                this % b (1:n) = this % b (1:n) / maxval (abs(this % b (1:n)))
                total = sum (this % b (1:n))
                if (total .ne. 0) this % b (1:n) = this % b (1:n) / total

                x = 0
                do i = 1, n
                    x = x + this % b (i) * (this % x (:,i) + this % alpha * this % e (:,i))
                enddo

                this % i_shift = m
                this % i_order = n

            endif

            m = mod (m, this % max_order) + 1

            this % x (:,m) = x
        endif

        if (this % err < this % max_err) this % status = 0
        if (this % i_iter > this % max_iter) this % status = -1

        this % i_iter = this % i_iter + 1

    end subroutine update_anderson

    subroutine p_free (this)
        class (t_fpi), intent (inout) :: this
        deallocate (this % atol)
        deallocate (this % rtol)
        deallocate (this % mat)
        deallocate (this % x)
        deallocate (this % e)
        deallocate (this % b)
        deallocate (this % w)
        call this % lud % free
    end subroutine p_free

end module m_fpi