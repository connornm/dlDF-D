#! /usr/bin/env python
#
def beale_f ( x, n ):

#*****************************************************************************80
#
## BEALE_F evaluates the Beale function.
#
#  Discussion:
#
#    The function is the sum of the squares of three functions.
#
#    This function has a valley approaching the line X(2) = 1.
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
#  Reference:
#
#    E Beale,
#    On an Iterative Method for Finding a Local Minimum of a Function
#    of More than One Variable,
#    Technical Report 25, Statistical Techniques Research Group,
#    Princeton University, 1958.
#
#    Richard Brent,
#    Algorithms for Finding Zeros and Extrema of Functions Without
#    Calculating Derivatives,
#    Stanford University Technical Report STAN-CS-71-198.
#
#  Parameters:
#
#    Input, real X(N), the evaluation point.
#
#    Input, integer N, the number of variables.
#
#    Output, real VALUE, the function value.
#
  c1 = 1.5
  c2 = 2.25
  c3 = 2.625

  fx1 = c1 - x[0] * ( 1.0 - x[1]      )
  fx2 = c2 - x[0] * ( 1.0 - x[1] ** 2 )
  fx3 = c3 - x[0] * ( 1.0 - x[1] ** 3 )

  value = fx1 ** 2 + fx2 ** 2 + fx3 ** 2

  return value

def beale_test ( ):

#*****************************************************************************80
#
## BEALE_TEST calls PRAXIS for the Beale function.
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
  import numpy as np
  import platform
  from praxis import praxis
  from r8vec_print import r8vec_print
  
  n = 2

  print ( '' )
  print ( 'BEALE_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  The Beale function.' )

  t0 = 0.00001
  h0 = 0.25
  prin = 0

  x = np.array ( [ 0.1, 0.1 ] )

  r8vec_print ( n, x, '  Initial point:' )

  print ( '  Function value = %g' % ( beale_f ( x, n ) ) )

  pr, x = praxis ( t0, h0, n, prin, x, beale_f )

  r8vec_print ( n, x, '  Computed minimizer:' )

  print ( '  Function value = %g' % ( beale_f ( x, n ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'BEALE_TEST:' )
  print ( '  Normal end of execution.' )
  return
 
if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  beale_test ( )
  timestamp ( )

