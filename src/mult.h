interface
    subroutine arg_mult (n, t, u, x, tol) bind(c, name='mult')
        integer, intent (in) :: n
        real(8), intent (in) :: u (n), t, tol
        real(8), intent (out) :: x (n)
    end subroutine arg_mult
end interface