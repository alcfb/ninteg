subroutine erms (n, dx, x0, atol, rtol, res)
    implicit none
    integer, intent (in) :: n
    real (8), intent (in) :: dx(n), x0(n), atol(n), rtol(n)
    real (8), intent (out) :: res
    integer :: i
    res = 0.D0
    do i = 1, n
        res = res + (dx (i) / (rtol (i) * abs(x0 (i)) + atol (i)))**2
    enddo
    res = sqrt (res / n)
end subroutine erms
