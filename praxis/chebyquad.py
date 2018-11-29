#! /usr/bin/env python
#
def chebyquad_f ( x, n ):

#*****************************************************************************80
#
## CHEBYQUAD_F evaluates the Chebyquad function.
#
#  Discussion:
#
#    The function is formed by the sum of squares of N separate terms.
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

  fvec = np.zeros ( n )

  for j in range ( 0, n ):

    t1 = 1.0;
    t2 = 2.0 * x[j] - 1.0
    t = 2.0 * t2

    for i in range ( 0, n ):
      fvec[i] = fvec[i] + t2
      th = t * t2 - t1
      t1 = t2
      t2 = th

  for i in range ( 0, n ):
    fvec[i] = fvec[i] / n
    if ( ( i % 2 ) == 1 ):
      fvec[i] = fvec[i] + 1.0 / ( i * ( i + 2 ) )
#
#  Compute F.
#
  value = 0.0
  for i in range ( 0, n ):
    value = value + fvec[i] ** 2

  return value

def chebyquad_test ( ):

#*****************************************************************************80
#
## CHEBYQUAD_TEST calls PRAXIS for the Chebyquad function.
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

  n = 8

  print ( '' )
  print ( 'CHEBYQUAD_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Chebyquad function.' )

  t0 = 0.00001
  h0 = 0.1
  prin = 0

  x = np.zeros ( n )

  for i in range ( 0, n ):
    x[i] = float ( i + 1 ) / float ( n + 1 )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( chebyquad_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, chebyquad_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( chebyquad_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'CHEBYQUAD_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  chebyquad_test ( )
  timestamp ( )

