#! /usr/bin/env python
#
def singular_f ( x, n ):

#*****************************************************************************80
#
## SINGULAR_F evaluates the Powell Singular function.
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
  value = 0.0

  for j in range ( 0, n, 4 ):

    if ( j + 1 <= n - 1 ):
      xjp1 = x[j+1]
    else:
      xjp1 = 0.0

    if ( j + 2 <= n - 1 ):
      xjp2 = x[j+2]
    else:
      xjp2 = 0.0

    if ( j + 3 <= n - 1 ):
      xjp3 = x[j+3]
    else:
      xjp3 = 0.0

    f1 = x[j] + 10.0 * xjp1

    if ( j + 1 <= n - 1 ):
      f2 = xjp2 - xjp3
    else:
      f2 = 0.0

    if ( j + 2 <= n - 1 ):
      f3 = xjp1 - 2.0 * xjp2
    else:
      f3 = 0.0

    if ( j + 3 <= n - 1 ):
      f4 = x[j] - xjp3
    else:
      f4 = 0.0

    value = value \
      +        f1 ** 2 \
      +  5.0 * f2 ** 2 \
      +        f3 ** 4 \
      + 10.0 * f4 ** 4

  return value

def singular_test ( ):

#*****************************************************************************80
#
## SINGULAR_TEST calls PRAXIS for the Powell Singular function.
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
  print ( 'SINGULAR_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Powell Singular function.' )

  t0 = 0.00001
  h0 = 1.0
  prin = 0

  x = np.array ( [ 3.0, -1.0, 0.0, 1.0 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( singular_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, singular_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( singular_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'SINGULAR_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  singular_test ( )
  timestamp ( )

