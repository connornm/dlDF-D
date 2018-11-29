#! /usr/bin/env python
#
def rosenbrock_f ( x, n ):

#*****************************************************************************80
#
## ROSENBROCK_F evaluates the Rosenbrock function.
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
  value = 0.0

  for j in range ( 0, n ):
    if ( ( j % 2 ) == 0 ):
      value = value + ( 1.0 - x[j] ) ** 2
    else:
      value = value + 100.0 * ( x[j] - x[j-1] ** 2 ) ** 2

  return value

def rosenbrock_test ( ):

#*****************************************************************************80
#
## ROSENBROCK_TEST calls PRAXIS for the Rosenbrock function.
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

  n = 2

  print ( '' )
  print ( 'ROSENBROCK_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Rosenbrock function.' )

  t0 = 0.00001
  h0 = 1.0
  prin = 0

  x = np.array ( [ -1.2, 1.0 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( rosenbrock_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, rosenbrock_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( rosenbrock_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'ROSENBROCK_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  rosenbrock_test ( )
  timestamp ( )

