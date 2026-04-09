module m_nonlin

    use m_fpi, only : t_fpi
    use m_parameters, only : MAX_FPI_ORDER
    implicit none

    type, public :: t_nonlin
        integer :: status, i_mult_iter, i_solv_iter, rank, max_iter, imethod
        real (8), allocatable :: x1(:), x2(:)
        real (8) :: max_err
        procedure (arg_mult), pointer, nopass :: mult1 => NULL()
        procedure (arg_solv), pointer, nopass :: solv1 => NULL()
        type (t_fpi) :: fpi
    contains
        procedure :: make => p_make
        procedure :: deft => p_deft
        procedure :: solv => p_solv
        procedure :: free => p_free
        procedure :: stds => p_stds
    end type t_nonlin

    include "mult.h"
    include "solv.h"

    private
    contains

    subroutine p_make (this, rank)

        implicit none
        integer, intent (in) :: rank
        class (t_nonlin), intent (inout) :: this

        this % rank = rank

        allocate (this % x1(rank))
        allocate (this % x2(rank))

        call this % fpi % make (rank, MAX_FPI_ORDER)

        this % status = -1

    end subroutine p_make

    subroutine p_deft (this)
        implicit none
        class (t_nonlin), intent (inout) :: this
        call this % fpi % default
        this % status = -1
        this % max_iter = 25
        this % max_err = 1
        this % imethod = 22
        this % i_mult_iter = 0
        this % i_solv_iter = 0
    end subroutine p_deft

    subroutine p_solv (this, n, h, t, b, x)
        ! Selection of method
        ! im : func  iteration
        ! 11 : solv  -
        ! 12 : solv  picard
        ! 13 : solv  anderson
        ! 21 : mult  picard
        ! 22 : mult  anderson
        implicit none
        class (t_nonlin), intent (inout) :: this
        integer, intent (in) :: n
        real(8), intent (in) :: b (n), h, t
        real(8), intent (out) :: x (n)

        if ((this % imethod == 11).and.(.not. associated (this % solv1))) call error_message ("'solv' must be provided")
        if ((this % imethod == 12).and.(.not. associated (this % solv1))) call error_message ("'solv' must be provided")
        if ((this % imethod == 13).and.(.not. associated (this % solv1))) call error_message ("'solv' must be provided")
        if ((this % imethod == 21).and.(.not. associated (this % mult1))) call error_message ("'mult' must be provided")
        if ((this % imethod == 22).and.(.not. associated (this % mult1))) call error_message ("'mult' must be provided")

        this % i_mult_iter = 0
        this % i_solv_iter = 0

        if (this % imethod == 11) then
            call this % solv1 (n, h, t, b, x, this % max_err * this % fpi % rtol (1))
            this % i_solv_iter = this % i_solv_iter + 1
            this % status = 0
        endif

        if (this % imethod == 12) this % fpi % imethod = 1
        if (this % imethod == 21) this % fpi % imethod = 1
        if (this % imethod == 13) this % fpi % imethod = 2
        if (this % imethod == 22) this % fpi % imethod = 2

        if ((this % imethod == 12).or.(this % imethod == 13)) then
            call this % fpi % init
            this % fpi % max_err = this % max_err
            this % fpi % max_iter = this % max_iter
            do while (this % fpi % status == 1)
                call this % solv1 (n, h, t, b, x, this % max_err * this % fpi % rtol (1))
                this % i_solv_iter = this % i_solv_iter + 1
                call this % fpi % update (x, n)
            enddo
            if (this % fpi % status  < 0) this % status = -2
            if (this % fpi % status == 0) this % status =  0
        endif

        if ((this % imethod == 21).or.(this % imethod == 22)) then
            call this % fpi % init
            this % fpi % max_err = this % max_err
            this % fpi % max_iter = this % max_iter
            this % x1 = x
            call this % fpi % update (this % x1, n)
            do while (this % fpi % status == 1)
                call this % mult1 (n, t, this % x1, this % x2, this % max_err * this % fpi % rtol (1))
                this % i_mult_iter = this % i_mult_iter + 1
                this % x2 = this % x2 * h + b
                call this % fpi % update (this % x2, n)
                this % x1 = this % x2
            enddo
            if (this % fpi % status  < 0) this % status = -2
            if (this % fpi % status == 0) this % status =  0
            x = this % x1
        endif

    end subroutine p_solv

    subroutine p_free (this)
        implicit none
        class (t_nonlin), intent (inout) :: this
        deallocate (this % x1)
        deallocate (this % x2)
        call this % fpi % free
    end subroutine p_free

    subroutine p_stds (this, n, x)
        ! Selection of method
        ! im : func  iteration
        ! 21 : mult  picard
        ! 22 : mult  anderson
        implicit none
        class (t_nonlin), intent (inout) :: this
        integer, intent (in) :: n
        real(8), intent (inout) :: x (n)
        real(8) :: t = 0.D0

        if ((this % imethod .ne. 21).and.(this % imethod .ne. 22)) call error_message ("'stds' supports only methods 21 and 22")
        if (.not. associated (this % mult1)) call error_message ("'mult' must be provided for 'stds'")

        if (this % imethod == 21) this % fpi % imethod = 1
        if (this % imethod == 22) this % fpi % imethod = 2

        this % i_mult_iter = 0
        this % i_solv_iter = 0

        call this % fpi % init
        this % fpi % max_err = this % max_err
        this % fpi % max_iter = this % max_iter
        this % x1 = x
        call this % fpi % update (this % x1, n)
        do while (this % fpi % status == 1)
            call this % mult1 (n, t, this % x1, this % x2, this % max_err * this % fpi % rtol (1))
            this % i_mult_iter = this % i_mult_iter + 1
            this % x2 = this % x2 + this % x1
            call this % fpi % update (this % x2, n)
            this % x1 = this % x2
        enddo
        if (this % fpi % status  < 0) this % status = -2
        if (this % fpi % status == 0) this % status =  0
        x = this % x1
    end subroutine p_stds

end module m_nonlin