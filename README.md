# DUALSIMPLEX code with Driver for Delaunay interpolation

The DUALSIMPLEX code is an open source Fortran 90 module for solving
linear programs of the form
``
max c^T x
``
``
s.t. A x <= b
``
where, c is a cost vector and A x <= b is a system of linear inequality
constraints.

DUALSIMPLEX solves this problem by applying the revised simplex method
of Dantzig to the asymetric dual
``
min b^T y
``
``
s.t. A^T y = c 
``
and
``
y >= 0.
``

To solve the dual problem, the revised simplex algorithm is applied. In
each iteration of the dual simplex algorithm, a complete LU factorization
is performed via DGETRF. Dantzig's rule is used to pivot, until either a
solution is found, or a pivot does not improve the objective. When a pivot
fails to improve the objective, Bland's rule is used until improvement
resumes, at which point Dantzig's pivoting strategy is resumed.

This strategy is most effective when A is dense and when M > N. For
efficient memory access patterns, the constraints are specified by
inputting A^T instead of A. For efficient linear algebra, LAPACK is
used for computing LU factorizations and performing triangular solve
operations, and BLAS is used for computing dot products and matrix-vector
multiplication.

While a complete LU factorization in each iteration produces maximum
numerical stability, a Woodbury update could alternatively be performed
if iteration speed is of greater concern. When M >> N, the cost of the
LU factorization is negligible compared to other computational costs,
and the difference in speed would be negligible.

## Delaunay interpolation

The problem of piecewise linear multivariate interpolation can be posed
as a linear program, where the constraint equations are determined by the
positions of the data points in the input (independent variable) space,
and the cost vector is determined by the interpolation point.

For input points P = { p1, ..., pn } and interpolation point q in d-dimensions, 
the constraint equations are given by
 - Ai = [ -(pi) 1 ], where Ai denotes the ith row of A;
 - bi = Ai dot Ai, where bi denotes the ith element of b and "dot" denotes
   the standard inner product; and
 - c = [-q 1]^T.

Then for the following linear program, a basis is primal feasibility 
if and only if it corresponds to the vertices of a Delaunay simplex,
and a basis is dual feasible if and only if it corresponds to the vertices
of a simplex containing q.
``
max c^T x
``
``
s.t. A x <= b
``

So, it follows that any basic solution to the above problem results in the
Delaunay simplex containing q, and the nonzero dual weights are exactly
the convex (barycentric) weights needed for solving the interpolation
problem.

This repository provides open source code for solving the dual formulation
of an arbitrary linear programming problem DUALSIMPLEX, as well as a driver
program for solving Delaunay interpolation problems, for points in
general position.
Note that if the input points are not in general position, then the basic
solution to the dual problem may not be primal feasible, and the above
methodology will in general fail.

## Contents

This repo contains the following files:
 - dualsimplex.f90 contains an open source Fortran 90 library for solving
   linear programs of the form: max c^T x, such that Ax <= b, via the dual
   formulation using the revised simplex algorithm.
   Note that dualsimplex.f90 uses a dense update and always performs a full
   Dantzig pivot, making it inefficient for most industrial applications.
 - delaunayLPtest.f90 tests the DUALSIMPLEX algorithm by performing Delaunay
   interpolation on points in general position.
 - generate\_data.f90 generates data points in general position along with
   an interpolation point in their convex hull, for usage with
   delaunayLPtest.f90
 - Makefile builds the project and executes the main programs.

## The DUALSIMPLEX module
dualsimplex.f90 contains a Fortran 90 module containing two subroutines:

 - DUALSIMPLEX for solving an LP when the initial basis is known (Phase II of
   the simplex algorithm),

 - FEASIBLEBASIS for finding an initial dual feasible basis when none is
   known (Phase I of the simplex algorithm) by solving the auxiliary
   problem using DUALSIMPLEX.


DUALSIMPLEX(N, M, AT, B, C, IBASIS, X, Y, IERR, EPS, IBUDGET, OBASIS)

On input:
 - N is the integer number of variables in the primal problem.
 - M is the integer number of constraints in the primal problem.
 - AT(N,M) is the transpose of the real valued constraint matrix A.
 - B(M) is the real valued vector of upper bounds for the constraints.
 - C(N) is the real valued cost vector for the objective function.
 - IBASIS(N) is an integer valued vector containing the indices (from AT)
   of an intitial basis that is dual feasible.

On output:
 - X(N) is a real valued vector, which contains the primal solution.
 - Y(M) is a real valued vector, which contains the dual solution.
 - IERR is an integer valued error flag. The error codes are listed in
   dualsimplex.f90, where DUALSIMPLEX is defined.

Optional arguments:
 - EPS contains the working precision for the problem. EPS must be a
   strictly positive real number, and by default EPS is the square-root of
   the unit roundoff.
 - IBUDGET contains the integer budget for the maximum number of pivots
   allowed. By default, IBUDGET=50,000.
 - When present, OBASIS(:) returns the integer indices of the final basis
   as listed in AT.

FEASIBLEBASIS (N, M, AT, C, BASIS, IERR, EPS, IBUDGET)

On input:
 - N is the integer number of variables in the primal problem.
 - M is the integer number of constraints in the primal problem.
 - AT(N,M) is the transpose of the real valued constraint matrix A.
 - C(N) is the real valued cost vector for the objective function.

On output:
 - BASIS(N) is an integer valued vector, which contains the indices of a
   dual feasible basis for AT.
 - IERR is an integer valued error flag. The error codes are listed in
   dualsimplex.f90, where FEASIBLEBASIS is defined.

Optional arguments:
 - EPS contains the working precision for the problem. EPS must be a
   strictly positive real number, and by default EPS is the square-root of
   the unit roundoff.
 - IBUDGET contains the integer budget for the maximum number of pivots
   allowed. By default, IBUDGET=50,000.

## Installation

### Prerequisites

DUALSIMPLEX requires both BLAS and LAPACK for efficient linear algebra.

### Building

To install, pull this repo and run
``
make -B
``

## Author

* ** Tyler Chang ** - *Primary author*

## Further reading

See Fukuda's FAQ in Polyhedral Computation
https://www.cs.mcgill.ca/~fukuda/soft/polyfaq/polyfaq.html
