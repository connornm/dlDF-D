#! /usr/bin/env python
# Note: Requires python3.3+ for use of .copy()

from copy import copy

def flin ( n, jsearch, l, f, x, nf, v, q0, q1, qd0, qd1, qa, qb, qc ):

#*****************************************************************************80
#
## FLIN is the function of one variable to be minimized by MINNY.
#
#  Discussion:
#
#    F(X) is a scalar function of a vector argument X.
#
#    A minimizer of F(X) is sought along a line or parabola.
#
#    This function has been modified, by removing the occurrence of a
#    common block, so that it looks more like a "normal" function that does
#    not rely so strongly on peculiarities of FORTRAN.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2016
#
#  Author:
#
#    Original FORTRAN77 version by Richard Brent.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Richard Brent,
#    Algorithms for Minimization with Derivatives,
#    Prentice Hall, 1973,
#    Reprinted by Dover, 2002.
#
#  Parameters:
#
#    Input, integer N, the number of variables.
#
#    Input, integer JSEARCH, indicates the kind of search.
#    If J is a legal column index, linear search in direction of V(*,JSEARCH).
#    Otherwise, then the search is parabolic, based on X, Q0 and Q1.
#
#    Input, real L, is the parameter determining the particular
#    point at which F is to be evaluated.
#    For a linear search, L is the step size.
#    For a quadratic search, L is a parameter which specifies
#    a point in the plane of X, Q0 and Q1.
#
#    Input, real F ( X, N ), the function to be minimized.
#
#    Input, real X(N), the base point of the search.
#
#    Input/output, integer NF, the function evaluation counter.
#
#    Input, real V(N,N), a matrix whose columns constitute
#    search directions.
#
#    Input, real Q0(N), Q1(N), two auxiliary points used to
#    determine the plane when a quadratic search is performed.
#
#    Input, real QD0, QD1, values needed to compute the
#    coefficients QA, QB, QC.
#
#    Input/output, real QA, QB, QC, coefficients used to combine
#    Q0, X, and A1 if a quadratic search is used.  (Yes, technically
#    these are input quantities as well.)
#
#    Output, real VALUE, the value of the function at the
#    minimizing point.
#
  import numpy as np

  t = np.zeros ( n )
#
#  The search is linear.
#
  if ( 0 <= jsearch ):

    t[0:n] = x[0:n] + l * v[0:n,jsearch]
#
#  The search is along a parabolic space curve.
#
  else:

    qa =                 l * ( l - qd1 ) /       ( qd0 + qd1 ) / qd0
    qb = - ( l + qd0 ) *     ( l - qd1 ) / qd1                 / qd0
    qc =   ( l + qd0 ) * l               / qd1 / ( qd0 + qd1 )

    
    # t[0:n] = qa * q0[0:n] + qb * x[0:n] + qc * q1[0:n]
    # Connor Dolan: Above line doesn't seem to work. Replaced with loop:
    for i in range(n):
      t[i] = qa*q0[i] + qb*x[i] + qc*q1[i]

#
#  The function evaluation counter NF is incremented.
#
  nf = nf + 1
#
#  Evaluate the function.
#
  value = f ( t, n )

  return value, nf, qa, qb, qc

def minfit ( n, tol, a ):

#*****************************************************************************80
#
## MINFIT computes the singular value decomposition of an N by N array.
#
#  Discussion:
#
#    This is an improved version of the EISPACK routine MINFIT
#    restricted to the case M = N and P = 0.
#
#    The singular values of the array A are returned in Q.  A is
#    overwritten with the orthogonal matrix V such that U * diag(Q) = A * V,
#    where U is another orthogonal matrix.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2016
#
#  Author:
#
#    Original FORTRAN77 version by Richard Brent.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Richard Brent,
#    Algorithms for Minimization with Derivatives,
#    Prentice Hall, 1973,
#    Reprinted by Dover, 2002.
#
#    James Wilkinson, Christian Reinsch,
#    Handbook for Automatic Computation,
#    Volume II, Linear Algebra, Part 2,
#    Springer Verlag, 1971.
#
#    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow, Yasuhiko Ikebe,
#    Virginia Klema, Cleve Moler,
#    Matrix Eigensystem Routines, EISPACK Guide,
#    Lecture Notes in Computer Science, Volume 6,
#    Springer Verlag, 1976,
#    ISBN13: 978-3540075462,
#    LC: QA193.M37.
#
#  Parameters:
#
#    Input, integer N, the order of the matrix A.
#
#    Input, real TOL, a tolerance which determines when a vector
#    (a column or part of a column of the matrix) may be considered
#    "essentially" equal to zero.
#
#    Input/output, real A(N,N).  On input, an N by N array whose
#    singular value decomposition is desired.  On output, the
#    SVD orthogonal matrix factor V.
#
#    Output, real Q(N), the singular values.
#
  import numpy as np
  from r8_epsilon import r8_epsilon
  from r8_hypot import r8_hypot

  from r8vec_print import r8vec_print

  kt_max = 30

  e = np.zeros ( n )
  q = np.zeros ( n )
#
#  Householder's reduction to bidiagonal form.
#
  if ( n == 1 ):
    q[0] = a[0,0]
    a[0,0] = 1.0
    return a, q

  g = 0.0
  x = 0.0

  for i in range ( 0, n ):

    e[i] = g
    l = i + 1

    s = 0.0
    for i1 in range ( i, n ):
      s = s + a[i1,i] ** 2

    g = 0.0

    if ( tol <= s ):

      f = a[i,i]

      g = np.sqrt ( s )
      if ( 0.0 <= f ):
        g = - g

      h = f * g - s
      a[i,i] = f - g

      for j in range ( l, n ):

        f = 0.0
        for i1 in range ( i, n ):
          f = f + a[i1,i] * a[i1,j]
        f = f / h

        for i1 in range ( i, n ):
          a[i1,j] = a[i1,j] + f * a[i1,i]

    q[i] = g

    s = 0.0
    for j1 in range ( l, n ):
      s = s + a[i,j1] ** 2

    g = 0.0

    if ( tol <= s ):

      if ( i < n - 1 ):
        f = a[i,i+1]

      g = np.sqrt ( s )
      if ( 0.0 <= f ):
        g = - g

      h = f * g - s

      if ( i < n - 1 ):

        a[i,i+1] = f - g

        for j1 in range ( l, n ):
          e[j1] = a[i,j1] / h

        for j in range ( l, n ):

          s = 0.0
          for j1 in range ( l, n ):
            s = s + a[j,j1] * a[i,j1]

          for j1 in range ( l, n ):
            a[j,j1] = a[j,j1] + s * e[j1]

    y = abs ( q[i] ) + abs ( e[i] )

    x = max ( x, y )
#
#  Accumulation of right-hand transformations.
#
  a[n-1,n-1] = 1.0
  g = e[n-1]
  l = n - 1

  for i in range ( n - 2, -1, -1 ):

    if ( g != 0.0 ):

      h = a[i,i+1] * g

      for i1 in range ( l, n ):
        a[i1,i] = a[i,i1] / h

      for j in range ( l, n ):

        s = 0.0
        for j1 in range ( l, n ):
          s = s + a[i,j1] * a[j1,j]

        for i1 in range ( l, n ):
          a[i1,j] = a[i1,j] + s * a[i1,i]

    for j1 in range ( l, n ):
      a[i,j1] = 0.0

    for i1 in range ( l, n ):
      a[i1,i] = 0.0

    a[i,i] = 1.0

    g = e[i]

    l = i
#
#  Diagonalization of the bidiagonal form.
#
  epsx = r8_epsilon ( ) * x

  for k in range ( n - 1, -1, -1 ):

    kt = 0

    while ( True ):

      kt = kt + 1

      if ( kt_max < kt ):
        e[k] = 0.0
        print ( '' )
        print ( 'MINFIT - Fatal error!' )
        print ( '  The QR algorithm failed to converge.' )
        exit ( 'MINFIT - Fatal error!' )

      skip = False

      for l2 in range ( k, -1, -1 ):

        l = l2

        if ( abs ( e[l] ) <= epsx ):
          skip = True
          break

        if ( 0 < l ):
          if ( abs ( q[l-1] ) <= epsx ):
            break
#
#  Cancellation of E(L) if 1 < L.
#
      if ( not skip ):

        c = 0.0
        s = 1.0

        for i in range ( l, k + 1 ):

          f = s * e[i]
          e[i] = c * e[i]

          if ( abs ( f ) <= epsx ):
            break

          g = q[i]
#
#  q(i) = h = sqrt(g*g + f*f).
#
          h = r8_hypot ( f, g )

          q[i] = h

          if ( h == 0.0 ):
            g = 1.0
            h = 1.0

          c =   g / h
          s = - f / h
#
#  Test for convergence for this index K.
#
      z = q[k]

      if ( l == k ):
        if ( z < 0.0 ):
          q[k] = - z
          for i1 in range ( 0, n ):
            a[i1,k] = - a[i1,k]
        break
#
#  Shift from bottom 2*2 minor.
#
      x = q[l]
      y = q[k-1]
      g = e[k-1]
      h = e[k]
      f = ( ( y - z ) * ( y + z ) + ( g - h ) * ( g + h ) ) / ( 2.0 * h * y )

      g = r8_hypot ( f, 1.0 )

      if ( f < 0.0 ):
        temp = f - g
      else:
        temp = f + g

      f = ( ( x - z ) * ( x + z ) + h * ( y / temp - h ) ) / x
#
#  Next QR transformation.
#
      c = 1.0
      s = 1.0

      for i in range ( l + 1, k + 1 ):

        g = e[i]
        y = q[i]
        h = s * g
        g = g * c

        z = r8_hypot ( f, h )

        e[i-1] = z

        if ( z == 0.0 ):
          f = 1.0
          z = 1.0

        c = f / z
        s = h / z
        f =   x * c + g * s
        g = - x * s + g * c
        h = y * s
        y = y * c

        for j in range ( 0, n ):
          x = a[j,i-1]
          z = a[j,i]
          a[j,i-1] = x * c + z * s
          a[j,i] = - x * s + z * c

        z = r8_hypot ( f, h )

        q[i-1] = z

        if ( z == 0.0 ):
          f = 1.0
          z = 1.0

        c = f / z
        s = h / z
        f =   c * g + s * y
        x = - s * g + c * y

      e[l] = 0.0
      e[k] = f
      q[k] = x

  return a, q

def minfit_test ( ):

#*****************************************************************************80
#
## MINFIT_TEST tests MINFIT, which is a sort of SVD computation.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform
  from r8_epsilon import r8_epsilon
  from r8mat_print import r8mat_print
  from r8vec_print import r8vec_print

  n = 5

  print ( '' )
  print ( 'MINFIT_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  MINFIT computes part of the SVD of a matrix A.' )
  print ( '    SVD: A = U * D * V\'' )
  print ( '  MINFIT is given A, and returns the diagonal D' )
  print ( '  and the orthogonal matrix V.' )

  a = np.zeros ( [ n, n ] )

  for i in range ( 0, n ):
    a[i,i] = 2.0

  for i in range ( 0, n - 1 ):
    a[i,i+1] = -1.0

  for i in range ( 1, n ):
    a[i,i-1] = -1.0

  r8mat_print ( n, n, a, '  The matrix A:' )
#
#  Numpy's EPS function is not easy to find!
#
  eps = r8_epsilon ( )

  tol = np.sqrt ( eps )

  a, d = minfit ( n, tol, a )

  r8mat_print ( n, n, a, '  The vector V:' )

  r8vec_print ( n, d, '  The singular values D:' )
#
#  Because A is positive definite symmetric, the "missing" matrix V = U.
#
  print ( '' )
  print ( '  Because A is positive definite symmetric,' )
  print ( '  we can reconstruct it as A = V * D * V\'' )

  a2 = np.zeros ( [ n, n ] )

  for i in range ( 0, n ):
    for j in range ( 0, n ):
      for k in range ( 0, n ):
        a2[i,j] = a2[i,j] + a[i,k] * d[k] * a[j,k]

  r8mat_print ( n, n, a2, '  The product A2 = V * D * V\'' )
#
#  Terminate.
#
  print ( '' )
  print ( 'MINFIT_TEST:' )
  print ( '  Normal end of execution.' )
  return

def minny ( n, jsearch, nits, d2, x1, f1, fk, f, x, t, h, v, q0, q1, \
  nl, nf, dmin, ldt, fx, qa, qb, qc, qd0, qd1 ):

#*****************************************************************************80
#
## MINNY minimizes a scalar function of N variables along a line.
#
#  Discussion:
#
#    MINNY minimizes F along the line from X in the direction V(*,J) unless
#    J is less than 1, when a quadratic search is made in the plane
#    defined by Q0, Q1 and X.
#
#    If FK = true, then F1 is FLIN(X1).  Otherwise X1 and F1 are ignored
#    on entry unless final FX is greater than F1.
#
#    This function was modified by removing the common blocks
#    and the use of labeled statements, 28 July 2016.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 August 2016
#
#  Author:
#
#    Original FORTRAN77 version by Richard Brent.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Richard Brent,
#    Algorithms for Minimization with Derivatives,
#    Prentice Hall, 1973,
#    Reprinted by Dover, 2002.
#
#  Parameters:
#
#    Input, integer N, the number of variables.
#
#    Input, integer JSEARCH, indicates the kind of search.
#    If JSEARCH is a legal column index, linear search in direction of V(*,J).
#    Otherwise, the search is parabolic, based on X, Q0 and Q1.
#
#    Input, integer NITS, the maximum number of times the interval
#    may be halved to retry the calculation.
#
#    Input/output, real D2, is either zero, or an approximation to
#    the value of (1/2) times the second derivative of F.
#
#    Input/output, real X1, on entry, an estimate of the
#    distance from X to the minimum along V(*,J), or, if J = 0, a curve.
#    On output, the distance between X and the minimizer that was found.
#
#    Input/output, real F1, ?
#
#    Input, logical FK if FK is TRUE, then on input F1 contains
#    the value FLIN(X1).
#
#    Input, real F ( X, N ), the function to be minimized.  
#
#    Input/output, real X(N), ?
#
#    Input, real T, ?
#
#    Input, real H, ?
#
#    Input, real V(N,N), a matrix whose columns are direction
#    vectors along which the function may be minimized.
#
#    Input, real Q0(N), an auxiliary point used to define
#    a curve through X.
#
#    Input, real Q1(N), an auxiliary point used to define
#    a curve through X.
#
#    Input/output, integer NL, the number of linear searches.
#
#    Input/output, integer NF, the number of function evaluations.
#
#    Input, real DMIN, an estimate for the smallest eigenvalue.
#
#    Input, real LDT, the length of the step.
#
#    Input/output, real FX, the value of F(X,N).
#
#    Input/output, real QA, QB, QC, ?
#
#    Input, real QD0, QD1, ?.
#
  import numpy as np
  from r8_epsilon import r8_epsilon

  machep = r8_epsilon ( )
  small = machep ** 2
  m2 = np.sqrt ( machep )
  m4 = np.sqrt ( m2 )
  sf1 = f1
  sx1 = x1
  k = 0
  xm = 0.0
  fm = fx
  f0 = fx
  dz = ( d2 < machep )
#
#  Find the step size.
#
  s = np.linalg.norm ( x )

  if ( dz ):
    temp = dmin
  else:
    temp = d2

  t2 = m4 * np.sqrt ( abs ( fx ) / temp + s * ldt ) + m2 * ldt
  s = m4 * s + t
  if ( dz and s < t2 ):
    t2 = s

  t2 = max ( t2, small )
  t2 = min ( t2, 0.01 * h )

  if ( fk and f1 <= fm ):
    xm = x1
    fm = f1

  if ( ( not fk ) or abs ( x1 ) < t2 ):

    if ( 0.0 <= x1 ):
      temp = 1.0
    else:
      temp = - 1.0

    x1 = temp * t2

    f1, nf, qa, qb, qc = flin ( n, jsearch, x1, f, x, nf, v, \
      q0, q1, qd0, qd1, qa, qb, qc )

  if ( f1 <= fm ):
    xm = x1
    fm = f1
#
#  Evaluate FLIN at another point and estimate the second derivative.
#
  while ( True ):

    if ( dz ):

      if ( f1 <= f0 ):
        x2 = 2.0 * x1
      else:
        x2 = - x1

      f2, nf, qa, qb, qc = flin ( n, jsearch, x2, f, x, nf, v, \
        q0, q1, qd0, qd1, qa, qb, qc )

      if ( f2 <= fm ):
        xm = x2
        fm = f2

      d2 = ( x2 * ( f1 - f0 ) - x1 * ( f2 - f0 ) ) \
        / ( ( x1 * x2 ) * ( x1 - x2 ) )
#
#  Estimate the first derivative at 0.
#
    d1 = ( f1 - f0 ) / x1 - x1 * d2
    dz = True
#
#  Predict the minimum.
#
    if ( d2 <= small ):

      if ( 0.0 <= d1 ):
        x2 = - h
      else:
        x2 = h

    else:

      x2 = ( - 0.5 * d1 ) / d2

    if ( h < abs ( x2 ) ):

      if ( x2 <= 0.0 ):
        x2 = - h
      else:
        x2 = h
#
#  Evaluate F at the predicted minimum.
#
    ok = True

    while ( True ):

      f2, nf, qa, qb, qc = flin ( n, jsearch, x2, f, x, nf, v, \
        q0, q1, qd0, qd1, qa, qb, qc )

      if ( nits <= k or f2 <= f0 ):
        break

      k = k + 1

      if ( f0 < f1 and 0.0 < x1 * x2 ):
        ok = False
        break

      x2 = 0.5 * x2

    if ( ok ):
      break
#
#  Increment the one-dimensional search counter.
#
  nl = nl + 1

  if ( fm < f2 ):
    x2 = xm
  else:
    fm = f2
#
#  Get a new estimate of the second derivative.
#
  if ( small < abs ( x2 * ( x2 - x1 ) ) ):
    d2 = ( x2 * ( f1 - f0 ) - x1 * ( fm - f0 ) ) / ( ( x1 * x2 ) * ( x1 - x2 ) )
  else:
    if ( 0 < k ):
      d2 = 0.0

  d2 = max ( d2, small )

  x1 = x2
  fx = fm

  if ( sf1 < fx ):
    fx = sf1
    x1 = sx1
#
#  Update X for linear but not parabolic search.
#
  if ( 0 <= jsearch ):
    x[0:n] = x[0:n] + x1 * v[0:n,jsearch]

  return d2, x1, f1, x, nl, nf, fx, qa, qb, qc

def praxis ( t0, h0, n, prin, x, f ):

#*****************************************************************************80
#
## PRAXIS seeks an N-dimensional minimizer X of a scalar function F(X).
#
#  Discussion:
#
#    PRAXIS returns the minimum of the function F(X,N) of N variables
#    using the principal axis method.  The gradient of the function is
#    not required.
#
#    The approximating quadratic form is
#
#      Q(x') = F(x,n) + (1/2) * (x'-x)' * A * (x'-x)
#
#    where X is the best estimate of the minimum and
#
#      A = inverse(V') * D * inverse(V)
#
#    V(*,*) is the matrix of search directions
#    D(*) is the array of second differences.
#
#    If F(X) has continuous second derivatives near X0, then A will tend
#    to the hessian of F at X0 as X approaches X0.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 August 2016
#
#  Author:
#
#    Original FORTRAN77 version by Richard Brent.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Richard Brent,
#    Algorithms for Minimization with Derivatives,
#    Prentice Hall, 1973,
#    Reprinted by Dover, 2002.
#
#  Parameters:
#
#    Input, real T0, is a tolerance.  PRAXIS attempts to return
#    praxis = f(x) such that if X0 is the true local minimum near X, then
#    norm ( x - x0 ) < T0 + sqrt ( EPSILON ( X ) ) * norm ( X ),
#    where EPSILON ( X ) is the machine precision for X.
#
#    Input, real H0, is the maximum step size.  H0 should be
#    set to about the maximum distance from the initial guess to the minimum.
#    If H0 is set too large or too small, the initial rate of
#    convergence may be slow.
#
#    Input, integer N, the number of variables.
#
#    Input, integer PRIN, controls printing intermediate results.
#    0, nothing is printed.
#    1, F is printed after every n+1 or n+2 linear minimizations.
#       final X is printed, but intermediate X is printed only
#       if N is at most 4.
#    2, the scale factors and the principal values of the approximating
#       quadratic form are also printed.
#    3, X is also printed after every few linear minimizations.
#    4, the principal vectors of the approximating quadratic form are
#       also printed.
#
#    Input/output, real X(N), is an array containing on entry a
#    guess of the point of minimum, on return the estimated point of minimum.
#
#    Input, real F ( X, N ), the function to be minimized.
#
#    Output, real PRAXIS, the function value at the minimizer.
#
#  Local parameters:
#
#    Local, real DMIN, an estimate for the smallest eigenvalue.
#
#    Local, real FX, the value of F(X,N).
#
#    Local, logical ILLC, is TRUE if the system is ill-conditioned.
#
#    Local, real LDT, the length of the step.
#
#    Local, integer NF, the number of function evaluations.
#
#    Local, integer NL, the number of linear searches.
#
  import numpy as np
  from r8_epsilon import r8_epsilon
  from r8mat_print import r8mat_print
  from r8vec_print import r8vec_print
#
#  Initialization.
#
  machep = r8_epsilon ( )
  small = machep * machep
  vsmall = small * small
  large = 1.0 / small
  vlarge = 1.0 / vsmall
  m2 = np.sqrt ( machep )
  m4 = np.sqrt ( m2 )
#
#  Heuristic numbers:
#
#  If the axes may be badly scaled (which is to be avoided if
#  possible), then set SCBD = 10.  Otherwise set SCBD = 1.
#
#  If the problem is known to be ill-conditioned, initialize ILLC = true.
#
#  KTM is the number of iterations without improvement before the
#  algorithm terminates.  KTM = 4 is very cautious usually KTM = 1
#  is satisfactory.
#
  scbd = 1.0
  illc = False
  ktm = 1

  if ( illc ):
    ldfac = 0.1
  else:
    ldfac = 0.01

  kt = 0
  nl = 0
  nf = 1
  fx = f ( x, n )
  qf1 = fx
  t = small + abs ( t0 )
  t2 = t
  dmin = small
  h = h0
  h = max ( h, 100.0 * t )
  ldt = h
#
#  The initial set of search directions V is the identity matrix.
#
  v = np.zeros ( [ n, n ] )
  for i in range ( 0, n ):
    v[i,i] = 1.0

  d = np.zeros ( n )
  y = np.zeros ( n )
  z = np.zeros ( n )
  qa = 0.0
  qb = 0.0
  qc = 0.0
  qd0 = 0.0
  qd1 = 0.0
  q0 = x.copy ( )
  q1 = x.copy ( )

  if ( 0 < prin ):
    print2 ( n, x, prin, fx, nf, nl )
#
#  The main loop starts here.
#
  while ( True ):

    sf = d[0]
    d[0] = 0.0
#
#  Minimize along the first direction V(*,1).
#
    jsearch = 0
    nits = 2
    d2 = d[0]
    s = 0.0
    value = fx
    fk = False

    d2, s, value, x, nl, nf, fx, qa, qb, qc = minny ( n, jsearch, nits, \
      d2, s, value, fk, f, x, t, h, v, q0, q1, nl, nf, dmin, ldt, \
      fx, qa, qb, qc, qd0, qd1 )

    d[0] = d2

    if ( s <= 0.0 ):
      for i1 in range ( 0, n ):
        v[i1,0] = - v[i1,0]

    if ( sf <= 0.9 * d[0] or d[0] <= 0.9 * sf ):
      d[1:n] = 0.0
#
#  The inner loop starts here.
#
    for k in range ( 2, n + 1 ):

      y = x.copy ( )
      sf = fx

      if ( 0 < kt ):
        illc = True

      while ( True ):

        kl = k
        df = 0.0
#
#  A random step follows, to avoid resolution valleys.
#
        if ( illc ):

          for j in range ( 0, n ):
            r = np.random.rand ( 1 )
            s = ( 0.1 * ldt + t2 * 10.0 ** kt ) * ( r - 0.5 )
            z[j] = s
            x[0:n] = x[0:n] + s * v[0:n,j]

          fx = f ( x, n )
          nf = nf + 1
#
#  Minimize along the "non-conjugate" directions V(*,K),...,V(*,N).
#
        for k2 in range ( k, n + 1 ):

          sl = fx

          jsearch = k2 - 1
          nits = 2
          d2 = d[k2-1]
          s = 0.0
          value = fx
          fk = False

          d2, s, value, x, nl, nf, fx, qa, qb, qc = minny ( n, jsearch, nits, \
            d2, s, value, fk, f, x, t, h, v, q0, q1, nl, nf, dmin, ldt, \
            fx, qa, qb, qc, qd0, qd1 )

          d[k2-1] = d2

          if ( illc ):
            s = d[k2-1] * ( ( s + z[k2-1] ) ** 2 )
          else:
            s = sl - fx

          if ( df <= s ):
            df = s
            kl = k2
#
#  If there was not much improvement on the first try, set
#  ILLC = true and start the inner loop again.
#
        if ( illc ):
          break

        if ( abs ( 100.0 * machep * fx ) <= df ):
          break

        illc = True

      if ( k == 2 and 1 < prin ):
        r8vec_print ( n, d, '  The second difference array:' )
#
#  Minimize along the "conjugate" directions V(*,1),...,V(*,K-1).
#
      for k2 in range ( 1, k ):

        jsearch = k2 - 1
        nits = 2
        d2 = d[k2-1]
        s = 0.0
        value = fx
        fk = False

        d2, s, value, x, nl, nf, fx, qa, qb, qc = minny ( n, jsearch, nits, \
          d2, s, value, fk, f, x, t, h, v, q0, q1, nl, nf, dmin, ldt, \
          fx, qa, qb, qc, qd0, qd1 )

        d[k2-1] = d2

      f1 = fx
      fx = sf

      for i in range ( 0, n ):
        temp = x[i]
        x[i] = y[i]
        y[i] = temp - y[i]

      lds = np.linalg.norm ( y )
#
#  Discard direction V(*,kl).
#
#  If no random step was taken, V(*,KL) is the "non-conjugate"
#  direction along which the greatest improvement was made.
#
      if ( small < lds ):

        for j in range ( kl - 1, k - 1, -1 ):
          v[0:n,j] = v[0:n,j-1]
          d[j] = d[j-1]

        d[k-1] = 0.0

        v[0:n,k-1] = y[0:n] / lds
#
#  Minimize along the new "conjugate" direction V(*,k), which is
#  the normalized vector:  (new x) - (old x).
#
        jsearch = k - 1
        nits = 4
        d2 = d[k-1]
        value = f1
        fk = True

        d2, lds, value, x, nl, nf, fx, qa, qb, qc = minny ( n, jsearch, nits, \
          d2, lds, value, fk, f, x, t, h, v, q0, q1, nl, nf, dmin, ldt, \
          fx, qa, qb, qc, qd0, qd1 )

        d[k-1] = d2

        if ( lds <= 0.0 ):
          lds = - lds
          v[0:n,k-1] = - v[0:n,k-1]

      ldt = ldfac * ldt
      ldt = max ( ldt, lds )

      if ( 0 < prin ):
        print2 ( n, x, prin, fx, nf, nl )

      t2 = m2 * np.linalg.norm ( x ) + t
#
#  See whether the length of the step taken since starting the
#  inner loop exceeds half the tolerance.
#
      if ( 0.5 * t2 < ldt ):
        kt = - 1

      kt = kt + 1

      if ( ktm < kt ):

        if ( 0 < prin ):
          r8vec_print ( n, x, '  X:' )

        value = fx

        return value, x
#
#  The inner loop ends here.
#
#  Try quadratic extrapolation in case we are in a curved valley.
#


    x, q0, q1, nl, nf, fx, qf1, qa, qb, qc, qd0, qd1 = quad ( n, f, x, t, h, v, q0, q1, nl, nf, dmin, ldt, fx, qf1, qa, qb, qc, qd0, qd1 )

    for j in range ( 0, n ):
      d[j] = 1.0 / np.sqrt ( d[j] )

    dn = max ( d )

    if ( 3 < prin ):
      r8mat_print ( n, n, v, '  The new direction vectors:' )

    for j in range ( 0, n ):
      v[0:n,j] = ( d[j] / dn ) * v[0:n,j]
#
#  Scale the axes to try to reduce the condition number.
#
    if ( 1.0 < scbd ):

      for i in range ( 0, n ):
        s = 0.0
        for j in range ( 0, n ):
          s = s + v[i,j] ** 2
        s = np.sqrt ( s )
        z[i] = max ( m4, s )

      s = min ( z )

      for i in range ( 0, n ):

        sl = s / z[i]
        z[i] = 1.0 / sl

        if ( scbd < z[i] ):
          sl = 1.0 / scbd
          z[i] = scbd

        v[i,0:n] = sl * v[i,0:n]
#
#  Calculate a new set of orthogonal directions before repeating
#  the main loop.
#
#  Transpose V for MINFIT:
#
    v = np.transpose ( v )
#
#  Call MINFIT to find the singular value decomposition of V.
#
#  This gives the principal values and principal directions of the
#  approximating quadratic form without squaring the condition number.
#
    v, d = minfit ( n, vsmall, v )
#
#  Unscale the axes.
#
    if ( 1.0 < scbd ):

      for i in range ( 0, n ):
        v[i,0:n] = z[i] * v[i,0:n]

      for j in range ( 0, n ):

        s = 0.0
        for i1 in range ( 0, n ):
          s = x + v[i1,j] ** 2
        s = sqrt ( s )

        d[j] = s * d[j]
        v[0:n,j] = v[0:n,j] / s

    for i in range ( 0, n ):

      dni = dn * d[i]

      if ( large < dni ):
        d[i] = vsmall
      elif ( dni < small ):
        d[i] = vlarge
      else:
        d[i] = 1.0 / dni ** 2
#
#  Sort the singular values and singular vectors.
#
    d, v = svsort ( n, d, v )
#
#  Determine the smallest eigenvalue.
#
    dmin = max ( d[n-1], small )
#
#  The ratio of the smallest to largest eigenvalue determines whether
#  the system is ill conditioned.
#
    if ( dmin < m2 * d[0] ):
      illc = True
    else:
      illc = False

    if ( 1 < prin ):

      if ( 1.0 < scbd ):
        r8vec_print ( n, z, '  The scale factors:' )

      r8vec_print ( n, d, '  Principal values of the quadratic form:' )

    if ( 3 < prin ):
      r8mat_print ( n, n, v, '  The principal axes:' )
#
#  The main loop ends here.
#
  if ( 0 < prin ):
    r8vec_print ( n, x, '  X:' )

  value = fx

  return value, x

def print2 ( n, x, prin, fx, nf, nl ):

#*****************************************************************************80
#
## PRINT2 prints certain data about the progress of the iteration.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    02 August 2016
#
#  Author:
#
#    Original FORTRAN77 version by Richard Brent.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Richard Brent,
#    Algorithms for Minimization with Derivatives,
#    Prentice Hall, 1973,
#    Reprinted by Dover, 2002.
#
#  Parameters:
#
#    Input, integer N, the number of variables.
#
#    Input, real X(N), the current estimate of the minimizer.
#
#    Input, integer PRIN, the user-specifed print level.
#    0, nothing is printed.
#    1, F is printed after every n+1 or n+2 linear minimizations.
#       final X is printed, but intermediate X is printed only
#       if N is at most 4.
#    2, the scale factors and the principal values of the approximating
#       quadratic form are also printed.
#    3, X is also printed after every few linear minimizations.
#    4, the principal vectors of the approximating quadratic form are
#       also printed.
#
#    Input, real FX, the smallest value of F(X) found so far.
#
#    Input, integer NF, the number of function evaluations.
#
#    Input, integer NL, the number of linear searches.
#
  from r8vec_print import r8vec_print

  print ( '' )
  print ( '  Linear searches      %d' % ( nl ) )
  print ( '  Function evaluations %d' % ( nf ) )
  print ( '  The function value FX = %g' % ( fx ) )

  if ( n <= 4 or 2 < prin ):
    r8vec_print ( n, x, '  X:' )

  return

def quad ( n, f, x, t, h, v, q0, q1, nl, nf, dmin, ldt, fx, qf1, qa, qb, \
  qc, qd0, qd1 ):

#*****************************************************************************80
#
## QUAD seeks to minimize the scalar function F along a particular curve.
#
#  Discussion:
#
#    The minimizer to be sought is required to lie on a curve defined
#    by Q0, Q1 and X.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2016
#
#  Author:
#
#    Original FORTRAN77 version by Richard Brent.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Richard Brent,
#    Algorithms for Minimization with Derivatives,
#    Prentice Hall, 1973,
#    Reprinted by Dover, 2002.
#
#  Parameters:
#
#    Input, integer N, the number of variables.
#
#    Input, real F ( x, n ), the function.
#
#    Input/output, real X(N), ?
#
#    Input, real T, ?
#
#    Input, rea H, ?
#
#    Input, real V(N,N), the matrix of search directions.
#
#    Input/output, real Q0(N), an auxiliary point used to define
#    a curve through X.
#
#    Input/output, real Q1(N), an auxiliary point used to define
#    a curve through X.
#
#    Input/output, integer NL, the number of linear searches.
#
#    Input/output, integer NF, the number of function evaluations.
#
#    Input, real DMIN, an estimate for the smallest eigenvalue.
#
#    Input, real LDT, the length of the step.
#
#    Input/output, real FX, the value of F(X,N).
#
#    Input/output, real QF1, QA, QB, QC, QD0, QD1, ?
#
  import numpy as np

  temp = fx
  fx   = qf1
  qf1  = temp

  temp = x
  x = q1
  q1 = temp
  print(x)
  print(q1)
  qd1 = np.linalg.norm ( np.array(x) - np.array(q1) )

  l = qd1
  s = 0.0

  if ( qd0 <= 0.0 or qd1 <= 0.0 or nl < 3 * n * n ):

    fx = qf1
    qa = 0.0
    qb = 0.0
    qc = 1.0

  else:

    jsearch = -1
    nits = 2
    value = qf1
    fk = True

    s, l, value, x, nl, nf, fx, qa, qb, qc = minny ( n, jsearch, nits, \
      s, l, value, fk, f, x, t, h, v, q0, q1, nl, nf, dmin, ldt, \
      fx, qa, qb, qc, qd0, qd1 )

    qa =                 l * ( l - qd1 )       / ( qd0 + qd1 ) / qd0
    qb = - ( l + qd0 )     * ( l - qd1 ) / qd1                 / qd0
    qc =   ( l + qd0 ) * l               / qd1 / ( qd0 + qd1 )

  qd0 = qd1

  xnew = np.zeros ( n )

#  xnew[0:n] = float(qa) * float(q0[0:n]) + float(qb) * float(x[0:n]) + float(qc) * float(q1[0:n])
#  Connor Dolan: Previous code above does not work (?). Replaced with a loop:
  for i in range(n):
    xnew[i] = float(qa)*float(q0[i]) + float(qb)*float(x[i]) + float(qc)*float(q1[i])

  q0[0:n] = x[0:n]
  x[0:n] = xnew[0:n]

  return x, q0, q1, nl, nf, fx, qf1, qa, qb, qc, qd0, qd1

def svsort ( n, d, v ):

#*****************************************************************************80
#
## SVSORT descending sorts singular values D and adjusts V.
#
#  Discussion:
#
#    A simple bubble sort is used on D.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2016
#
#  Author:
#
#    Original FORTRAN77 version by Richard Brent.
#    Python version by John Burkardt.
#
#  Reference:
#
#    Richard Brent,
#    Algorithms for Minimization with Derivatives,
#    Prentice Hall, 1973,
#    Reprinted by Dover, 2002.
#
#  Parameters:
#
#    Input, integer N, the length of D, and the order of V.
#
#    Input/output, real D(N), the vector to be sorted.
#    On output, the entries of D are in descending order.
#
#    Input/output, real V(N,N), an N by N array to be adjusted
#    as D is sorted.  In particular, if the value that was in D(I) on input is
#    moved to D(J) on output, then the input column V(*,I) is moved to
#    the output column V(*,J).
#
  for j in range ( 0, n - 1 ):

    j3 = j
    for j2 in range ( j + 1, n ):
      if ( d[j3] < d[j2] ):
        j3 = j2

    t     = d[j]
    d[j]  = d[j3]
    d[j3] = t

    for i in range ( 0, n ):
      t       = v[i,j]
      v[i,j]  = v[i,j3]
      v[i,j3] = t

  return d, v

def svsort_test ( ):

#*****************************************************************************80
#
## SVSORT_TEST tests SVSORT, which sorts singular value information.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    01 August 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform

  n = 5

  print ( '' )
  print ( 'SVSORT_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  SVSORT sorts a vector D, and the corresponding columns' )
  print ( '  of a matrix V.' )

  d = np.random.rand ( n )

  v = np.zeros ( [ n, n ] )

  for i in range ( 0, n ):
    for j in range ( 0, n ):
      v[i,j] = 10 * ( i + 1 ) + ( j + 1 )

  print ( '' )
  print ( '  First row = entries of D.' )
  print ( '  Corresponding columns of V below.' )
  print ( '' )
  for j in range ( 0, n ):
    print ( '%14.6g' % ( d[j] ) ),
  print ( '' )
  print ( '' )
  for i in range ( 0, n ):
    for j in range ( 0, n ):
      print ( '%14.6g' % ( v[i,j] ) ),
    print ( '' )

  d, v = svsort ( n, d, v )

  print ( '' )
  print ( '  After sorting D and rearranging V:' )
  print ( '' )
  for j in range ( 0, n ):
    print ( '%14.6g' % ( d[j] ) ),
  print ( '' )
  print ( '' )
  for i in range ( 0, n ):
    for j in range ( 0, n ):
      print ( '%14.6g' % ( v[i,j] ) ),
    print ( '' )
#
#  Terminate.
#
  print ( '' )
  print ( 'SVSORT_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  minfit_test ( )
  svsort_test ( )
  timestamp ( )
