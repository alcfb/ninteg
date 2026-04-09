subroutine LUsolve (n, A, b)

    implicit none
    integer :: i, j, k
    integer, intent (in) :: n
    real(8), intent (inout) :: A (n, n)
    real(8), intent (inout) :: b (n)
    real(8) :: sum

    do j = 1, n
        do i= 1, j
            do k = 1, i - 1
                A (i, j) = A (i, j) - A (i, k) * A (k, j)
            enddo
        end do

        do i = j + 1, n
            sum = A (i, j)
            do k = 1, j - 1
                sum = sum - A (i, k) * A (k, j)
            enddo
            A (i, j) = sum / A (j, j)
        end do
    end do

    do i = 1, n
        do k = 1, i - 1
            b(i) = b(i) - A (i, k) * b(k)
        end do
    end do

    do i = n, 1, -1
        sum = b (i)
        do k = i + 1, n
            sum = sum - A (i, k) * b (k)
        enddo
        b (i) = sum / A (i, i)
    end do

return
end subroutine LUsolve