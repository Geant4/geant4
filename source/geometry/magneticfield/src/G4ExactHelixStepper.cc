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
// $Id: G4ExactHelixStepper.cc 66356 2012-12-18 09:02:32Z gcosmo $ 
//
//  Helix a-la-Explicity Euler: x_1 = x_0 + helix(h)
//   with helix(h) being a helix piece of length h
//   simplest approach for solving linear differential equations.
//  Take the current derivative and add it to the current position.
//
//  As the field is assumed constant, an error is not calculated.
// 
//  Author: J. Apostolakis, 28 Jan 2005
//     Implementation adapted from ExplicitEuler of W.Wander 
// -------------------------------------------------------------------

#include "G4ExactHelixStepper.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4LineSection.hh"

G4ExactHelixStepper::G4ExactHelixStepper(G4Mag_EqRhs *EqRhs)
  : G4MagHelicalStepper(EqRhs),
    fBfieldValue(DBL_MAX, DBL_MAX, DBL_MAX),
    fPtrMagEqOfMot(EqRhs)
{
  ;
}

G4ExactHelixStepper::~G4ExactHelixStepper() {} 

void
G4ExactHelixStepper::Stepper( const G4double yInput[],
                              const G4double*,
                                    G4double hstep,
                                    G4double yOut[],
                                    G4double yErr[]      )
{  
   const G4int nvar = 6;

   G4int i;
   G4ThreeVector Bfld_value;

   MagFieldEvaluate(yInput, Bfld_value);        
   AdvanceHelix(yInput, Bfld_value, hstep, yOut);

  // We are assuming a constant field: helix is exact
  //
  for(i=0;i<nvar;i++)
  {
    yErr[i] = 0.0 ;
  }

  fBfieldValue=Bfld_value;
}

void
G4ExactHelixStepper::DumbStepper( const G4double  yIn[],
                                        G4ThreeVector   Bfld,
                                        G4double  h,
                                        G4double  yOut[])
{
  // Assuming a constant field: solution is a helix

  AdvanceHelix(yIn, Bfld, h, yOut);

  G4Exception("G4ExactHelixStepper::DumbStepper",
              "GeomField0002", FatalException,
              "Should not be called. Stepper must do all the work." ); 
}  


// ---------------------------------------------------------------------------

G4double
G4ExactHelixStepper::DistChord() const 
{
  // Implementation : must check whether h/R >  pi  !!
  //   If( h/R <  pi)   DistChord=h/2*std::tan(Ang_curve/4)
  //   Else             DistChord=R_helix

  G4double distChord;
  G4double Ang_curve=GetAngCurve();

  if (Ang_curve<=pi)
  {
    distChord=GetRadHelix()*(1-std::cos(0.5*Ang_curve));
  }
  else if(Ang_curve<twopi)
  {
    distChord=GetRadHelix()*(1+std::cos(0.5*(twopi-Ang_curve)));
  }
  else
  {
    distChord=2.*GetRadHelix();  
  }

  return distChord;
}   

G4int
G4ExactHelixStepper::IntegratorOrder() const 
{
  return 1; 
}
