module m_py_ninteg

    use iso_c_binding
    use m_ivp, only : t_ivp
    implicit none

    type (t_ivp) :: ivp

    logical :: is_allocated = .false.

    include "mult.h"
    include "solv.h"

    contains

    subroutine make (rank) bind(c, name='make')
        integer(c_int), intent(in) :: rank
        if (is_allocated) call clean
        call ivp % make (rank)
        call ivp % deft
        is_allocated = .true.
    end subroutine make

    subroutine init_t0 (t0) bind(c, name='init_t0')
        real(c_double), intent(in) :: t0
        ivp % t0 = t0
    end subroutine init_t0

    subroutine init_t1 (t1) bind(c, name='init_t1')
        real(c_double), intent(in) :: t1
        ivp % t1 = t1
    end subroutine init_t1

    subroutine init_h0 (h0) bind(c, name='init_h0')
        real(c_double), intent(in) :: h0
        ivp % h0 = h0
    end subroutine init_h0

    subroutine init_x0 (n, x) bind(c, name='init_x0')
        integer(c_int), intent(in) :: n
        real(c_double), intent(in) :: x(n)
        ivp % x(:) = x(:)
    end subroutine init_x0

    subroutine init_rtol (rtol) bind(c, name='init_rtol')
        real(c_double), intent(in) :: rtol
        ivp % rtol(:) = rtol
    end subroutine init_rtol

    subroutine init_atol (atol) bind(c, name='init_atol')
        real(c_double), intent(in) :: atol
        ivp % atol(:) = atol
    end subroutine init_atol

    subroutine init_mult (mult) bind(c, name='init_mult')
        procedure (arg_mult), pointer, intent(in) :: mult
        if (associated(mult)) ivp % mult => mult
    end subroutine init_mult

    subroutine init_solv (solv) bind(c, name='init_solv')
        procedure (arg_solv), pointer, intent(in) :: solv
        if (associated(solv)) ivp % solv => solv
    end subroutine init_solv

    subroutine init_im (im) bind(c, name='init_im')
        integer(c_int), intent(in) :: im
        ivp % imethod = im
    end subroutine init_im

    subroutine init_q_max (q) bind(c, name='init_q_max')
        integer(c_int), intent(in) :: q
        ivp % ctrl % qmax = q
    end subroutine init_q_max

    subroutine init_h_min (h) bind(c, name='init_h_min')
        real(c_double), intent(in) :: h
        ivp % ctrl % hmin = h
    end subroutine init_h_min

    subroutine init_h_max (h) bind(c, name='init_h_max')
        real(c_double), intent(in) :: h
        ivp % ctrl % hmax = h
    end subroutine init_h_max

    subroutine get_info (val) bind (c, name='get_info')

        real (c_double), intent(out) :: val (15)

        val ( 1) = real (ivp % i_iter, kind=c_double)
        val ( 2) = real (ivp % i_step, kind=c_double)
        val ( 3) = real (ivp % i_iter_rejects, kind=c_double)
        val ( 4) = real (ivp % i_step_rejects, kind=c_double)
        val ( 5) = real (ivp % i_trials, kind=c_double)
        val ( 6) = real (ivp % n_calls, kind=c_double)
        val ( 7) = real (ivp % q0, kind=c_double)
        val ( 8) = real (ivp % q1, kind=c_double)
        val ( 9) = ivp % t
        val (10) = ivp % h0
        val (11) = ivp % h1
        val (12) = ivp % lre
        val (13) = real (ivp % status, kind=c_double)
        val (14) = 0
        val (15) = 0

    end subroutine get_info


    subroutine get_t (t) bind(c, name='get_t')
        real (c_double), intent(out) :: t
        t = ivp % t
    end subroutine get_t

    subroutine get_x (n, x) bind(c, name='get_x')
        integer (c_int), intent(in) :: n
        real (c_double), intent(out) :: x(n)
        if (size (ivp % x) /= n) call error_message ('ninteg.f90', 'Size mismatch')
        x(:) = ivp % x (:)
    end subroutine get_x

    subroutine ivp_init () bind(c, name='ivp_init')
        call ivp % init
    end subroutine ivp_init

    subroutine ivp_next () bind (c, name='ivp_next')
        if (.not. is_allocated) call error_message ('ninteg.f90', "IVP not initialized. Call make() first.")
        call ivp % next
    end subroutine ivp_next

    subroutine clean () bind(c, name='clean')
        if (is_allocated) then
            call ivp % free
            is_allocated = .false.
        endif
    end subroutine clean

end module m_py_ninteg