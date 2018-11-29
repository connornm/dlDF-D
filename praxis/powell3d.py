#! /usr/bin/env python
#
def powell3d_f ( x, n ):

#*****************************************************************************80
#
## POWELL3D_F evaluates the Powell 3D function.
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
#  Reference:
#
#    M J D Powell,
#    An Efficient Method for Finding the Minimum of a Function of
#    Several Variables Without Calculating Derivatives,
#    Computer Journal,
#    Volume 7, Number 2, pages 155-162, 1964.
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

  value = 3.0 - 1.0 / ( 1.0 + ( x[0] - x[1] ) ** 2 ) \
    - np.sin ( 0.5 * np.pi * x[1] * x[2] ) \
    - np.exp ( - ( ( x[0] - 2.0 * x[1] + x[2] ) / x[1] ) ** 2 )

  return value

def powell3d_test ( ):

#*****************************************************************************80
#
## POWELL3D_TEST calls PRAXIS for the Powell 3D function.
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
  print ( 'POWELL3D_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Powell 3D function.' )

  t0 = 0.00001
  h0 = 1.0
  prin = 0

  x = np.array ( [ 0.0, 1.0, 2.0 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( powell3d_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, powell3d_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( powell3d_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'POWELL3D_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  powell3d_test ( )
  timestamp ( )

