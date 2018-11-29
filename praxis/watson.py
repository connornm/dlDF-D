#! /usr/bin/env python
#
def watson_f ( x, n ):

#*****************************************************************************80
#
## WATSON_F evaluates the Watson function.
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

  for i in range ( 0, 29 ):

    s1 = 0.0
    d = 1.0
    for j in range ( 1, n ):
      s1 = s1 + j * d * x[j]
      d = d * ( i + 1 ) / 29.0

    s2 = 0.0
    d = 1.0
    for j in range ( 0, n ):
      s2 = s2 + d * x[j]
      d = d * ( i + 1 ) / 29.0

    value = value + ( s1 - s2 * s2 - 1.0 ) ** 2

  value = value + x[0] ** 2 + ( x[1] - x[0] ** 2 - 1.0 ) ** 2
 
  return value

def watson_test ( ):

#*****************************************************************************80
#
## WATSON_TEST calls PRAXIS for the Watson function.
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

  n = 6

  print ( '' )
  print ( 'WATSON_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Watson function.' )

  t0 = 0.00001
  h0 = 1.0
  prin = 0

  x = np.zeros ( n )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( watson_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, watson_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( watson_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'WATSON_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  watson_test ( )
  timestamp ( )

