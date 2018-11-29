#! /usr/bin/env python
#
def praxis_test ( ):

#*****************************************************************************80
#
## PRAXIS_TEST tests the PRAXIS library.
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
  import platform
  from beale        import beale_test
  from box          import box_test
  from chebyquad    import chebyquad_test
  from cube         import cube_test
  from helix        import helix_test
  from hilbert      import hilbert_test
  from praxis       import minfit_test
  from powell3d     import powell3d_test
  from rosenbrock   import rosenbrock_test
  from singular     import singular_test
  from praxis       import svsort_test
  from tridiagonal  import tridiagonal_test
  from watson       import watson_test
  from wood         import wood_test

  print ( '' )
  print ( 'PRAXIS_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  Test the PRAXIS library.' )
#
#  Minimization tests.
#
  beale_test ( )
  box_test ( )
  chebyquad_test ( )
  cube_test ( )
  helix_test ( )
  hilbert_test ( )
  powell3d_test ( )
  rosenbrock_test ( )
  singular_test ( )
  tridiagonal_test ( )
  watson_test ( )
  wood_test ( )
#
#  Utility tests.
#
  minfit_test ( )
  svsort_test ( )
#
#  Terminate.
#
  print ( '' )
  print ( 'PRAXIS_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  praxis_test ( )
  timestamp ( )
 

