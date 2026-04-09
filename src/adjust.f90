subroutine adjust_nordsieck (m, n, z, h0, h1)
    implicit none
    integer, intent (in) :: n, m
    real(8), intent (inout) :: z (m, n)
    real(8), intent (in) :: h0, h1
    integer :: i
    real(8) :: ratio
    ratio = h1/h0
    do i = 2, m
        z (i,:) = ratio**(i-1) * z(i,:)
    enddo
end subroutine adjust_nordsieck