module m_ctrl

    ! Step size and order control module
    ! Build in accordance with LSODE algorithm

    use m_parameters

    implicit none

    type, public :: t_ctrl
        integer :: h_hold, q_hold
        real(8) :: hmin, hmax
        integer :: qmax
        character*1 :: todo
    contains
        procedure :: make => p_make
        procedure :: init => p_init
        procedure :: select => p_select
        procedure :: accept => p_accept
        procedure :: free => p_free
    end type t_ctrl

    private
    contains

    subroutine p_make (this)
        implicit none
        class (t_ctrl), intent (inout) :: this
    end subroutine p_make

    subroutine p_init (this)
        implicit none
        class (t_ctrl), intent (inout) :: this
        this % q_hold = 0
        this % h_hold = 0
    end subroutine p_init

    subroutine p_select (this, q0, h0, rs, ru, rd, q, h)
        ! Select step size
        implicit none
        class (t_ctrl), intent (inout) :: this
        integer, intent (in) :: q0
        real(8), intent (in) :: h0, rs, ru, rd
        real(8), intent (out):: h
        integer, intent (out):: q

        if (q0 == MAX_BDF_ORDER) then
            if (rs >=rd) this % todo = "s"
            if (rs < rd) this % todo = "d"
        elseif (q0 == 1) then
            if (rs >=ru) this % todo = "s"
            if (rs < ru) this % todo = "u"
        else
            if ((rs >=ru).and.(rs >= rd)) this % todo = "s"
            if ((ru >=rs).and.(ru >= rd)) this % todo = "u"
            if ((rd >=ru).and.(rd >= rs)) this % todo = "d"
        endif

        if (this % q_hold > 0) this % todo = "s"
        if (this % h_hold > 0) then
            if (rs > 1) this % todo = "o"
        endif

        if (this % todo == "s") then
            h = rs * h0
            q = q0
        elseif (this % todo == "u") then
            h = ru * h0
            q = q0 + 1
        elseif (this % todo == "d") then
            h = rd * h0
            q = q0 - 1
        else
            h = h0
            q = q0
        endif

        h = min (h, 5*h0)

        h = min (h, this % hmax)
        h = max (h, this % hmin)
        q = min (q, this % qmax)

    end subroutine p_select

    subroutine p_accept (this)
        implicit none
        class (t_ctrl), intent (inout) :: this
        if ((this % todo == "d").or.(this % todo == "u")) then
            this % q_hold = 3
            this % h_hold = 3
        endif
        if (this % todo == "s") then
            this % q_hold = max (0, this % q_hold - 1)
            this % h_hold = 3
        endif
        if (this % todo == "o") then
            this % h_hold = max (0, this % h_hold - 1)
            this % q_hold = max (0, this % q_hold - 1)
        endif
    end subroutine p_accept

    subroutine p_free (this)
        implicit none
        class (t_ctrl), intent (inout) :: this
    end subroutine p_free

end module m_ctrl