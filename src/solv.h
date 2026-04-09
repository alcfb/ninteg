interface
    subroutine arg_solv (n, h, t, b, x, tol) bind(c, name='solv')
        integer, intent (in) :: n
        real(8), intent (in) :: b (n), h, t, tol
        real(8), intent (out) :: x (n)
    end subroutine arg_solv
end interface