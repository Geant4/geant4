//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4MagHelicalStepper.cc 97598 2016-06-06 07:19:46Z gcosmo $
//
// --------------------------------------------------------------------

#include "G4MagHelicalStepper.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4LineSection.hh"
#include "G4Mag_EqRhs.hh"

// given a purely magnetic field a better approach than adding a straight line
// (as in the normal runge-kutta-methods) is to add helix segments to the
// current position


// Constant for determining unit conversion when using normal as integrand.
//
const G4double G4MagHelicalStepper::fUnitConstant = 0.299792458*(GeV/(tesla*m));


G4MagHelicalStepper::G4MagHelicalStepper(G4Mag_EqRhs *EqRhs)
   : G4MagIntegratorStepper(EqRhs, 6), // integrate over 6 variables only !!
                                       // position & velocity
     fPtrMagEqOfMot(EqRhs), fAngCurve(0.), frCurve(0.), frHelix(0.)
{
}

G4MagHelicalStepper::~G4MagHelicalStepper()
{
}

void
G4MagHelicalStepper::AdvanceHelix( const G4double  yIn[],
                                   G4ThreeVector   Bfld,    
                                   G4double  h,
                                   G4double  yHelix[],
                                   G4double  yHelix2[] )
{
  // const G4int    nvar = 6;
 
  // OLD  const G4double approc_limit = 0.05;
  // OLD  approc_limit = 0.05 gives max.error=x^5/5!=(0.05)^5/5!=2.6*e-9
  // NEW  approc_limit = 0.005 gives max.error=x^5/5!=2.6*e-14

  const G4double approc_limit = 0.005;
  G4ThreeVector  Bnorm, B_x_P, vperp, vpar;

  G4double B_d_P;
  G4double B_v_P;
  G4double Theta;
  G4double R_1;
  G4double R_Helix;
  G4double CosT2, SinT2, CosT, SinT;
  G4ThreeVector positionMove, endTangent;

  G4double Bmag = Bfld.mag();
  const G4double *pIn = yIn+3;
  G4ThreeVector initVelocity= G4ThreeVector( pIn[0], pIn[1], pIn[2]);
  G4double      velocityVal = initVelocity.mag();
  G4ThreeVector initTangent = (1.0/velocityVal) * initVelocity;
  
  R_1=GetInverseCurve(velocityVal,Bmag);

  // for too small magnetic fields there is no curvature
  // (include momentum here) FIXME

  if( (std::fabs(R_1) < 1e-10)||(Bmag<1e-12) )
  {
    LinearStep( yIn, h, yHelix );

    // Store and/or calculate parameters for chord distance

    SetAngCurve(1.);     
    SetCurve(h);
    SetRadHelix(0.);
  }
  else
  {
    Bnorm = (1.0/Bmag)*Bfld;

    // calculate the direction of the force
    
    B_x_P = Bnorm.cross(initTangent);

    // parallel and perp vectors

    B_d_P = Bnorm.dot(initTangent); // this is the fraction of P parallel to B

    vpar = B_d_P * Bnorm;       // the component parallel      to B
    vperp= initTangent - vpar;  // the component perpendicular to B
    
    B_v_P  = std::sqrt( 1 - B_d_P * B_d_P); // Fraction of P perp to B

    // calculate  the stepping angle
  
    Theta   = R_1 * h; // * B_v_P;

    // Trigonometrix
      
    if( std::fabs(Theta) > approc_limit )
    {
       SinT     = std::sin(Theta);
       CosT     = std::cos(Theta);
    }
    else
    {
      G4double Theta2 = Theta*Theta;
      G4double Theta3 = Theta2 * Theta;
      G4double Theta4 = Theta2 * Theta2;
      SinT     = Theta - 1.0/6.0 * Theta3;
      CosT     = 1 - 0.5 * Theta2 + 1.0/24.0 * Theta4;
    }

    // the actual "rotation"

    G4double R = 1.0 / R_1;

    positionMove  = R * ( SinT * vperp + (1-CosT) * B_x_P) + h * vpar;
    endTangent    = CosT * vperp + SinT * B_x_P + vpar;

    // Store the resulting position and tangent

    yHelix[0]   = yIn[0] + positionMove.x(); 
    yHelix[1]   = yIn[1] + positionMove.y(); 
    yHelix[2]   = yIn[2] + positionMove.z();
    yHelix[3] = velocityVal * endTangent.x();
    yHelix[4] = velocityVal * endTangent.y();
    yHelix[5] = velocityVal * endTangent.z();

    // Store 2*h step Helix if exist

    if(yHelix2)
    {
      SinT2     = 2.0 * SinT * CosT;
      CosT2     = 1.0 - 2.0 * SinT * SinT;
      endTangent    = (CosT2 * vperp + SinT2 * B_x_P + vpar);
      positionMove  = R * ( SinT2 * vperp + (1-CosT2) * B_x_P) + h*2 * vpar;
      
      yHelix2[0]   = yIn[0] + positionMove.x(); 
      yHelix2[1]   = yIn[1] + positionMove.y(); 
      yHelix2[2]   = yIn[2] + positionMove.z(); 
      yHelix2[3] = velocityVal * endTangent.x();
      yHelix2[4] = velocityVal * endTangent.y();
      yHelix2[5] = velocityVal * endTangent.z();
    }

    // Store and/or calculate parameters for chord distance

    G4double ptan=velocityVal*B_v_P;

    G4double particleCharge = fPtrMagEqOfMot->FCof() / (eplus*c_light); 
    R_Helix =std::abs( ptan/(fUnitConstant  * particleCharge*Bmag));
       
    SetAngCurve(std::abs(Theta));
    SetCurve(std::abs(R));
    SetRadHelix(R_Helix);
  }
}

//
//  Use the midpoint method to get an error estimate and correction
//  modified from G4ClassicalRK4: W.Wander <wwc@mit.edu> 12/09/97
//

void
G4MagHelicalStepper::Stepper( const G4double yInput[],
                              const G4double*,
                                    G4double hstep,
                                    G4double yOut[],
                                    G4double yErr[]  )
{  
   const G4int nvar = 6;

   G4int i;

   // correction for Richardson Extrapolation.
   // G4double  correction = 1. / ( (1 << IntegratorOrder()) -1 );
   
   G4double      yTemp[8], yIn[8] ;
   G4ThreeVector Bfld_initial, Bfld_midpoint;
   
   //  Saving yInput because yInput and yOut can be aliases for same array

   for(i=0;i<nvar;i++) { yIn[i]=yInput[i]; }

   G4double h = hstep * 0.5; 

   MagFieldEvaluate(yIn, Bfld_initial) ;      

   // Do two half steps

   DumbStepper(yIn,   Bfld_initial,  h, yTemp);
   MagFieldEvaluate(yTemp, Bfld_midpoint) ;     
   DumbStepper(yTemp, Bfld_midpoint, h, yOut); 

   // Do a full Step

   h = hstep ;
   DumbStepper(yIn, Bfld_initial, h, yTemp);

   // Error estimation

   for(i=0;i<nvar;i++)
   {
     yErr[i] = yOut[i] - yTemp[i] ;
   }
   
   return;
}

G4double
G4MagHelicalStepper::DistChord() const 
{
  // Check whether h/R >  pi  !!
  // Method DistLine is good only for <  pi

  G4double Ang=GetAngCurve();
  if(Ang<=pi)
  {
    return GetRadHelix()*(1-std::cos(0.5*Ang));
  }
  else
  {
    if(Ang<twopi)
    {
      return GetRadHelix()*(1+std::cos(0.5*(twopi-Ang)));
    }
    else  // return Diameter of projected circle
    {
      return 2*GetRadHelix();
    }
  }
}
