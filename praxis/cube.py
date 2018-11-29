#! /usr/bin/env python
#
def cube_f ( x, n ):

#*****************************************************************************80
#
## CUBE_F evaluates the Cube function.
#
#  Discussion:
#
#    The function is the sum of the squares of two functions.
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
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  fx1 = 10.0 * ( x[1] - x[0] ** 3 )
  fx2 = 1.0 - x[0]

  value = fx1 ** 2 + fx2 ** 2

  return value

def cube_test ( ):

#*****************************************************************************80
#
## CUBE_TEST calls PRAXIS for the Cube function.
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
  from praxis import praxis
  from r8vec_print import r8vec_print

  n = 2

  print ( '' )
  print ( 'CUBE_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Cube function.' )

  t0 = 0.00001
  h0 = 1.0
  prin = 0

  x = np.array ( [ -1.2, -1.0 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( cube_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, cube_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( cube_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'CUBE_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  cube_test ( )
  timestamp ( )

