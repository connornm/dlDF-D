#! /usr/bin/env python
#
def box_f ( x, n ):

#*****************************************************************************80
#
## BOX_F evaluates the Box function.
#
#  Discussion:
#
#    The function is formed by the sum of squares of 10 separate terms.
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
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  import numpy as np

  value = 0.0

  for i in range ( 1, 11 ):

    c = - i / 10.0

    fx = np.exp ( c * x[0] ) - np.exp ( c * x[1] ) \
      - x[2] * ( np.exp ( c ) - np.exp ( 10.0 * c ) )

    value = value + fx ** 2

  return value

def box_test ( ):

#*****************************************************************************80
#
## BOX_TEST calls PRAXIS for the Box function.
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
  from praxis import praxis
  from r8vec_print import r8vec_print

  n = 3

  print ( '' )
  print ( 'BOX_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Box function.' )

  t0 = 0.00001
  h0 = 20.0
  prin = 0

  x = np.array ( [ 0.0, 10.0, 20.0 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( box_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, box_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( box_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'BOX_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  box_test ( )
  timestamp ( )

