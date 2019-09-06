PROGRAM GENERATE_DELAUNAY_DATA
! Generate test data to be interpolated by Delaunay.
IMPLICIT NONE

! PROBLEM DIMENSIONS: 
! ADJUST VALUES D AND N TO GENERATE PROBLEMS OF VARYING SIZE.
!-----------------------------------------------------------!
INTEGER, PARAMETER :: D=20, N=500
!-----------------------------------------------------------!

! Get real precision.
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15, 307)
! Local variables.
REAL(KIND=DP) :: A(D), B(D) ! Temp arrays for storing output.
INTEGER :: I ! Loop index.
!INTEGER, ALLOCATABLE :: SEED(:)

! Initialize random seed.
CALL RANDOM_SEED()

! Store output in file deldata.txt below.
OPEN(D, FILE='deldata.txt')
WRITE(D, *) D, N, 1, 0
B(1:D) = 0.0_DP
! Begin generating data points at random.
DO I = 1, N 
   CALL RANDOM_NUMBER(A)
   WRITE(D, *) A(1:D)
   B(1:D) = B(1:D) + A(1:D)
END DO
! Get an interior point for interpolation.
B(1:D) = B(1:D) / REAL(N, KIND=DP)
WRITE(D, *) B(1:D)
CLOSE(D)

END PROGRAM GENERATE_DELAUNAY_DATA
