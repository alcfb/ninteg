module m_parameters

    integer, parameter :: MAX_BDF_ORDER = 6

    integer, parameter :: MAX_FPI_ORDER = 5

    integer, parameter :: PASCAL (7,7) = reshape ((/ &
    1,1,1,1,1, 1, 1, &
    0,1,2,3,4, 5, 6, &
    0,0,1,3,6,10,15, &
    0,0,0,1,4,10,20, &
    0,0,0,0,1, 5,15, &
    0,0,0,0,0, 1, 6, &
    0,0,0,0,0, 0, 1/), (/7,7/), order=(/2,1/))

    real(8), parameter :: LL (6,7) = reshape ((/ &
    0.D0,  1.D0,  0.D0,   0.D0,  0.D0,  0.D0, 0.D0, &
    0.D0,  1.D0,  1.D0,   0.D0,  0.D0,  0.D0, 0.D0, &
    0.D0,  2.D0,  3.D0,   1.D0,  0.D0,  0.D0, 0.D0, &
    0.D0,  6.D0, 11.D0,   6.D0,  1.D0,  0.D0, 0.D0, &
    0.D0, 24.D0, 50.D0,  35.D0, 10.D0,  1.D0, 0.D0, &
    0.D0,120.D0,274.D0, 225.D0, 85.D0, 15.D0, 1.D0  &
     /), (/6,7/), order=(/2,1/))

    integer, parameter :: FACT (0:7) = (/1,1,2,6,24,120,720,5040/)

end module m_parameters