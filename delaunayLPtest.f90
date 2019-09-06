PROGRAM SAMPLE_DELAUNAY
! Locate the Delaunay simplex containing a point Q using linear programming.
! Consider a set of points PTS = {P_i} \in R^d, i = 1, ..., n.
! Let A be a (d+1) X (d+1) matrix whose rows are given by: A_i = [ p_i, 1 ];
! let B be a n-vector B_i = p_i \dot p_i; and
! let C be a (d+1)-vector C_i = -Q_i for i = 1, ..., d, C_{d+1} = -1.
!
! If the problem
!
! \max C^T X
! s.t. AX \leq B.
!
! has a unique solution, then the vertices of the simplex S \in DT(PTS) that
! contains Q are given by the solution basis for solving and the affine
! interpolation weights are given by the dual solution for the corresponding
! basis.
!
! In such a situation, since the dual solution is unique, it can be solved
! for via the asymetric dual
!
! \min B^T Y
! s.t. A^T Y = C
! Y \geq 0.
!
! If the solution is nonunique, the problem is significantly harder and cannot
! be solved via the above methodology. Therefore, this code is to be used for
! demonstrative purposes only.

USE DUALSIMP_MOD
IMPLICIT NONE

! Declare inputs and local data.
INTEGER :: D, N, M, DUMMY
REAL(KIND=R8), ALLOCATABLE :: A(:,:), B(:), C(:,:), X(:), Y(:), WEIGHTS(:,:)
REAL(KIND=R8) :: START, FINISH, EPS
INTEGER, ALLOCATABLE :: SIMPLEX(:,:), BASIS(:), IERR(:)
INTEGER :: I, J
CHARACTER(LEN=80) :: FILEPATH

! Open the file path $(filepath), and get the metadata from the
! first line (D, N, and M).
CALL GET_COMMAND_ARGUMENT(1, FILEPATH)
OPEN(1, FILE=TRIM(FILEPATH))
READ(1, *) D, N, M, DUMMY
IF(D .LE. 0 .OR. N .LE. 0 .OR. M .LE. 0) THEN
   WRITE(*,*) "Illegal input dimensions in input file, line 1."
END IF
! Allocate all necessarry arrays.
ALLOCATE(A(D+1,N), B(N), C(D+1,M), Y(N), X(D+1), WEIGHTS(D+1,M), &
   & SIMPLEX(D+1,M), BASIS(D+1), IERR(M), STAT=I)
IF(I .NE. 0) THEN
   WRITE(*,*) "Memory allocation error."
END IF
! Read the input data/training points into the matrix PTS(:,:).
DO I = 1, N
   READ(1, *) A(1:D, I)
   A(D+1,I) = 1.0_R8
   B(I) = DOT_PRODUCT(A(1:D,I),A(1:D,I))
END DO
A = -A
! Read the interpolation points into the matrix C(:,:).
DO I = 1, M
   READ(1, *) C(1:D, I)
   C(D+1,I) = 1.0_R8
END DO
CLOSE(1)
C = -C

EPS = SQRT(EPSILON(0.0_R8))

! Compute the interpolation results and time.
CALL CPU_TIME(START)
!DO I = 1, 20
DO I = 1, M
   CALL FEASIBLEBASIS(D+1, N, A, C(:,I), BASIS, IERR(I),EPS=EPS)
   IF (IERR(I) .EQ. 0) THEN
      CALL DUALSIMPLEX(D+1, N, A, B, C(:,I), BASIS, X, Y, &
         & IERR(I), OBASIS=SIMPLEX(:,I),EPS=EPS)
      IF (IERR(I) .EQ. 0) THEN
         FORALL (J = 1 : D+1) WEIGHTS(J,I) = Y(SIMPLEX(J,I))
      END IF
   END IF 
END DO
!END DO
CALL CPU_TIME(FINISH)
FINISH = (FINISH - START)! / 20.0_R8
! Print the timing data.
DO I = 1, M
   IF (IERR(I) .EQ. 2) THEN
      WRITE(*,*) 'Extrapolation at point ', I, '. No solution computed.'
   ELSE IF (IERR(I) .NE. 0) THEN
      WRITE(*,14) 'Error at point ', I, '. IERR = ', IERR(I)
   ELSE
      WRITE(*,10) 'Interpolation point: ', -C(1:D,I)
      WRITE(*,11) 'Simplex : ', SIMPLEX(:,I)
      WRITE(*,10) 'Weights: ', WEIGHTS(:,I)
   END IF
END DO

WRITE(*,15) M, ' points interpolated in ', FINISH, ' seconds.'

10 FORMAT(1X,A,/,(1X,5ES15.7))
11 FORMAT(1X,A,/,(10I7))
14 FORMAT(1X,A,I7,A,I2)
15 FORMAT(/,I7,A,ES16.8,A,/)

END PROGRAM SAMPLE_DELAUNAY
