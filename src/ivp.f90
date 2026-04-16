module m_ivp

    use m_parameters
    use m_size
    use m_bdf
    use m_ctrl

    implicit none

    type, public :: t_ivp
        real(8), allocatable :: rwork (:)
        integer, allocatable :: iwork (:)
        integer :: status
        integer :: imethod
        integer :: rank, i_iter, i_step, i_step_rejects, i_iter_rejects, i_trials, n_calls
        real(8) :: t, t0, t1, lre
        real(8), allocatable :: x (:), atol(:), rtol(:)
        real(8), allocatable, private :: zc (:,:), z1 (:,:), dz(:), z0(:)
        real(8), allocatable, private :: e0(:), e1(:)
        real(8) :: h0, h1
        integer :: q0, q1
        type (t_bdf) :: bdf
        type (t_size) :: size
        type (t_ctrl) :: ctrl
        procedure (arg_mult), pointer, nopass :: mult => NULL()
        procedure (arg_solv), pointer, nopass :: solv => NULL()
    contains
        procedure :: make => p_make
        procedure :: deft => p_deft
        procedure :: init => p_init
        procedure :: next => p_next
        procedure :: free => p_free
    end type t_ivp

    include "mult.h"
    include "solv.h"

    type (t_ivp), pointer :: this_

    private
    contains

    subroutine p_make (this, rank)

        implicit none
        integer, intent (in) :: rank
        class (t_ivp), intent (inout), target :: this

        this % rank = rank

        allocate (this % x (rank))
        allocate (this % dz (rank))
        allocate (this % z0 (rank))
        allocate (this % e0 (rank))
        allocate (this % e1 (rank))
        allocate (this % zc (MAX_BDF_ORDER+1,rank))
        allocate (this % z1 (MAX_BDF_ORDER+1,rank))
        allocate (this % atol(rank))
        allocate (this % rtol(rank))

        call this % bdf % make (rank)
        call this % size % make (rank)
        call this % ctrl % make

        ! 0 - allocated memory
        ! 1 - initialized parameters
        ! 2 - computed state at the first point
        ! 3 - continue after successful step
        ! 4 - continue after rejected step
        ! 5 - computed state at the last point
        this % status = 0

        this_ => this

    end subroutine p_make


    subroutine p_deft (this)
        implicit none
        class (t_ivp), intent (inout) :: this
        if (.not. allocated (this % x)) call error_message ("IVP", "arrays are not allocated")
        this % t0 = 0.D0
        this % h0 = 1.D-4
        this % h1 = this % h0
        this % q0 = 1
        this % q1 = 1
        this % t1 = 1.D0
        this % x = 0.D0
        this % atol = 1.D-10
        this % rtol = 1.D-3
        this % status = 1
        this % imethod = 13 ! Anderson's scheme and solv
        this % ctrl % qmax = 5
        this % ctrl % hmin = 1.D-10
        this % ctrl % hmax = 1.D+10
    end subroutine p_deft


    subroutine p_init (this)

        implicit none
        class (t_ivp), intent (inout) :: this

        if (this % status < 1)  call error_message ("IVP", "'deft' must be called first")
        if (this % t0 > this % t1) call error_message ("IVP", "t0 is larger than t1")

        this % i_step = 0
        this % status = 2
        this % t = this % t0
        this % i_iter = 0
        this % i_step_rejects = 0
        this % i_iter_rejects = 0
        this % i_trials = 0
        this % n_calls = 0
        this % lre = 1

        call this % ctrl % init

    end subroutine p_init

    subroutine p_next (this)
        implicit none
        class (t_ivp), intent (inout) :: this
        integer :: q, i
        real(8) :: lre, h, ratio, ru, rs, rd

        ! Update pointers to the mult and solv subroutines
        this % bdf % nonlin % mult1 => this % mult
        this % bdf % nonlin % solv1 => this % solv
        ! Setup tolerance for inner iteration
        this % bdf % nonlin % fpi % atol = this % atol
        this % bdf % nonlin % fpi % rtol = this % rtol
        ! Setup iteration method
        this % bdf % nonlin % imethod = this % imethod

        if (this % status < 2) call error_message ("IVP","'init' must be called first")

        if (this % i_step == 0) then
            ! Always start with the first order
            this % q0 = 1
            this % zc = 0
            this % e0 = 0
            ! The 1st vector
            this % zc (1,:) = this % x
            this % zc (2,:) = 0.D0 ! <= initial guess
            ! The 2d vector
            !call this % mult (this % rank, this % t, this % x, this % zc (2,:))
            !this % n_calls = this % n_calls + 1
            !this % zc (2,:) = this % h0 * this % zc (2,:)
            ! Accept this step
            this % h1 = this % h0
            this % q1 = this % q0
            this % z1 = this % zc
            this % e1 = this % e0
            ! Local error
            this % lre = 1
            this % i_step = this % i_step + 1
            this % status = 2
        else
            q = this % q0
            ! Update ZC and ZP given q, t and h
            call this % bdf % next (this % rank, q, this % t, this % h0, this % zc, this % e0)
            this % n_calls = this % n_calls + this % bdf % nonlin % i_mult_iter
            this % n_calls = this % n_calls + this % bdf % nonlin % i_solv_iter
            ! Update local relative error
            call erms (this % rank, this % e0, this % x, this % atol, this % rtol, this % lre)
            this % lre = this % lre * FACT (q) / (q + 1)
            ! Check if reject is required
            if ((this % lre > 1).or.(this % bdf % nonlin % status .ne. 0)) then
                !write (*,*) 'LRE: ', this % lre
                !this % e0 = this % e0 / (this % x * this % rtol + this % atol)
                !write (*,*) maxval (this % e0 (1:144)), maxval (this % e0 (145:288))
                this % i_trials = this % i_trials + 1
                ! Recall state and setup new step size
                this % zc = this % z1
                this % q0 = this % q1
                this % e0 = this % e1
                this % h0 = this % h0 * 0.5D0
                ! Adjust ZC to the new step h
                call this % bdf % rescale (this % rank, this % zc, this % h1, this % h0)
                this % status = 4
            else
                ! Select time step
                call this % size % predict (this % rank, q, this % x, this % zc, this % e0, this % e1, this % atol, this % rtol, rs, ru, rd)
                call this % ctrl % select (this % q0, this % h0, rs, ru, rd, q, h)
                ! Adjust time step to the end of time
                !h = min (h, this % t1 - this % t0)
                ! call this % size % correct (h, this % t + this % h0, this % t1)
                ! Save the current time step and order
                this % z1 = this % zc
                this % e1 = this % e0
                this % h1 = this % h0
                this % q1 = this % q0
                ! Order increment
                if (q > this % q0) then
                    this % zc (q+1,:) = this % e0 / q
                endif
                ! Adjust ZC to the new step h
                if (this % ctrl % todo .ne. "o") then
                    call this % bdf % rescale (this % rank, this % zc, this % h0, h)
!                    this % e0 = this % e0 * (h / this % h0)**(this % q0+1)
!                    this % e1 = this % e1 * (h / this % h1)**(this % q1+1)
                endif
                ! Accept the new step size
                this % t = this % t + this % h1
                this % h0 = h
                this % q0 = q
                this % i_step = this % i_step + 1
                this % x = this % zc (1,:)
                call this % ctrl % accept
                ! Reset number of trials of rejected step
                if (this % i_trials > 0) then
                    this % i_iter_rejects = this % i_iter_rejects + this % i_trials
                    this % i_step_rejects = this % i_step_rejects + 1
                    this % i_trials = 0
                endif
                this % status = 3
            endif
        endif

        this % i_iter = this % i_iter + 1

        if (this % h0 < 0) call error_message ("IVP", "step size is negative")

        if (this % i_trials > 20) call error_message ("IVP", "too many trials at one time step")

        if (this % t >= this % t1) then

            ! The previous step h1 was too large
            ! Therefore, we must go back in time by h

            h = this % t - this % t1

            ! Rescale vector Zc to -h
            ! Shift Zc to the requested point t1
            ! Rescale Zc again to h0

            call this % bdf % rescale (this % rank, this % zc, this % h0, -h)
            call this % bdf % interpolate (this % rank, this % q0, this % zc)
            call this % bdf % rescale (this % rank, this % zc, -h, this % h0)

            ! Update the solution vector X

            this % x(:) = this % zc (1,:)

            ! Update the restart vector Z1

            this % h1 = this % h0 - h
            this % z1 = this % zc
            call this % bdf % rescale (this % rank, this % z1, -h, this % h1)

            this % t = this % t1

            this % status = 5
        endif

    end subroutine p_next

    subroutine p_free (this)
        implicit none
        class (t_ivp), intent (inout) :: this
        if (this % status < 0) call error_message ("IVP", "arrays were not allocated")
        deallocate (this % x)
        deallocate (this % z0)
        deallocate (this % dz)
        deallocate (this % e0)
        deallocate (this % e1)
        deallocate (this % zc)
        deallocate (this % z1)
        deallocate (this % atol)
        deallocate (this % rtol)
        call this % bdf % free
        call this % size % free
        call this % ctrl % free
        this % status = -1
    end subroutine p_free

end module m_ivp