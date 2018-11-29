#! /usr/bin/env python
#
def helix_f ( x, n ):

#*****************************************************************************80
#
## HELIX_F evaluates the Helix function.
#
#  Discussion:
#
#    The function is the sum of the squares of three functions.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    04 August 2016
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

  r = np.linalg.norm ( x )

  if ( 0.0 <= x[0] ):
    theta = 0.5 * np.arctan2 ( x[1], x[0] ) / np.pi
  else:
    theta = 0.5 * ( np.arctan2 ( x[1], x[0] ) + np.pi ) / np.pi

  fx1 = 10.0 * ( x[2] - 10.0 * theta )
  fx2 = 10.0 * ( r - 1.0 )
  fx3 = x[2]

  value = fx1 ** 2 + fx2 ** 2 + fx3 ** 2

  return value

def helix_test ( ):

#*****************************************************************************80
#
## HELIX_TEST calls PRAXIS for the Fletcher-Powell Helix function.
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

  n = 3

  print ( '' )
  print ( 'HELIX_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Fletcher-Powell Helix function.' )

  t0 = 0.00001
  h0 = 1.0
  prin = 0

  x = np.array ( [ -1.0, 0.0, 0.0 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( helix_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, helix_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( helix_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'HELIX_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  helix_test ( )
  timestamp ( )

