// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// Test program for G4AnalyticalPolSolver class. 
//

#include "G4ios.hh"
#include "globals.hh"
#include "G4AnalyticalPolSolver.hh"





// Testing quadratic
//  Re=  1.0000000000 Im=  1.0000000000 
//  Re=  1.0000000000 Im= -1.0000000000 

// Testing cubic
//  Re=  1.0000000000 Im=  0.0000000000 
//  Re=  2.0000000000 Im=  0.0000000000 
//  Re=  3.0000000000 Im=  0.0000000000 

// Testing biquadratic
//  Re=  0.0000000000 Im=  1.4142135624 
//  Re=  0.0000000000 Im= -1.4142135624 
//  Re= -0.0000000000 Im=  1.4142135624 
//  Re= -0.0000000000 Im= -1.4142135624 




int main()
{
  G4int i, k, n;
  G4double p[5], r[3][5];
  G4AnalyticalPolSolver solver;


  for ( n = 2; n <= 4; n++ )
  {
    // Various test cases 

    if( n == 4 )   // roots: 1,2,3,4 
    { 
      p[0] =  1.; 
      p[1] = -10.; 
      p[2] =  35.; 
      p[3] = -50.; 
      p[4] =  24.; 

      p[0] =  1.; 
      p[1] =  0.; 
      p[2] =  4.; 
      p[3] =  0.; 
      p[4] =  4.; // roots: +-i sqrt(2) 
    }
    if( n == 3 )  // roots:  1,2,3 
    { 
      p[0] =  1.; 
      p[1] = -6.; 
      p[2] = 11.; 
      p[3] = -6.; 
    }
    if(n==2)  // roots : 1 +- i 
    { 
      p[0] =  1.; 
      p[1] = -2.; 
      p[2] =  2.;
    }

    if( n == 2 )
    {
      G4cout<<"Test QUADROOTS(p,r):"<<G4endl;    
      i = solver.QUADROOTS(p,r);
    }
    else if( n == 3 )
    {
      G4cout<<"Test CUBICROOTS(p,r):"<<G4endl;    
      i = solver.CUBICROOTS(p,r);
    }
    else if( n == 4 )
    {
      G4cout<<"Test BIQUADROOTS(p,r):"<<G4endl;    
      i = solver.BIQUADROOTS(p,r);
    }

    for( k = 1; k <= n; k++ )
    { 
      G4cout << r[1][k] << " " << r[2][k] <<" i" << G4endl;
    }
  }
  return 0;
}

