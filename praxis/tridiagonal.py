#! /usr/bin/env python
#
def tridiagonal_f ( x, n ):

#*****************************************************************************80
#
## TRIDIAGONAL_F evaluates the tridiagonal function.
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
#    John Burkardt
#
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  value = x[0] ** 2 + 2.0 * sum ( x[1:n] ** 2 )

  for i in range ( 0, n - 1 ):
    value = value - 2.0 * x[i] * x[i+1]

  value = value - 2.0 * x[0]

  return value

def tridiagonal_test ( ):

#*****************************************************************************80
#
## TRIDIAGONAL_TEST calls PRAXIS for the Tridiagonal function.
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
#    John Burkardt
#
  import numpy as np
  import platform
  from praxis import praxis
  from r8vec_print import r8vec_print

  n = 4

  print ( '' )
  print ( 'TRIDIAGONAL_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Tridiagonal function.' )

  t0 = 0.00001
  h0 = 8.0
  prin = 0

  x = np.zeros ( n )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( tridiagonal_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, tridiagonal_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( tridiagonal_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'TRIDIAGONAL_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  tridiagonal_test ( )
  timestamp ( )

