#! /usr/bin/env python
#
def hilbert_f ( x, n ):

#*****************************************************************************80
#
## HILBERT_F evaluates the Hilbert function.
#
#  Discussion:
#
#    The function is a positive definite quadratic function of the form
#
#      f(x) = x' A x
#
#    where A is the Hilbert matrix, A(I,J) = 1/(I+J-1).
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

  for i in range ( 0, n ):
    for j in range ( 0, n ):
      value = value + x[i] * x[j] / ( i + j + 1 )

  return value

def hilbert_test ( ):

#*****************************************************************************80
#
## HILBERT_TEST calls PRAXIS for the Hilbert function.
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

  n = 10

  print ( '' )
  print ( 'HILBERT_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Hilbert function.' )

  t0 = 0.00001
  h0 = 10.0
  prin = 0

  x = np.ones ( n )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( hilbert_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, hilbert_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( hilbert_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'HILBERT_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  hilbert_test ( )
  timestamp ( )

