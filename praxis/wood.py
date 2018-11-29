#! /usr/bin/env python
#
def wood_f ( x, n ):

#*****************************************************************************80
#
## WOOD_F evaluates the Wood function.
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
  f1 = x[1] - x[0] ** 2
  f2 = 1.0 - x[0]
  f3 = x[3] - x[2] ** 2
  f4 = 1.0 - x[2]
  f5 = x[1] + x[3] - 2.0
  f6 = x[1] - x[3]

  value = \
      100.0 * f1 ** 2 \
    +         f2 ** 2 \
    +  90.0 * f3 ** 2 \
    +         f4 ** 2 \
    +  10.0 * f5 ** 2 \
    +   0.1 * f6 ** 2

  return value

def wood_test ( ):

#*****************************************************************************80
#
## WOOD_TEST calls PRAXIS for the Wood function.
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

  n = 4

  print ( '' )
  print ( 'WOOD_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Wood function.' )

  t0 = 0.00001
  h0 = 10.0
  prin = 0

  x = np.array ( [ -3.0, -1.0, -3.0, -1.0 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( wood_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, wood_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( wood_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'WOOD_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  wood_test ( )
  timestamp ( )

