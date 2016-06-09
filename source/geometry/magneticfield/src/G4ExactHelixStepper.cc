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
// $Id: G4ExactHelixStepper.cc,v 1.4 2006/06/29 18:23:50 gunter Exp $ 
// GEANT4 tag $Name: geant4-08-02 $
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
#include "G4ThreeVector.hh"
#include "G4LineSection.hh"


G4ExactHelixStepper::G4ExactHelixStepper(G4Mag_EqRhs *EqRhs)
  : G4MagHelicalStepper(EqRhs),
    fBfieldValue(DBL_MAX, DBL_MAX, DBL_MAX), yInitialEHS(DBL_MAX), yFinalEHS(-DBL_MAX),
    fLastStepSize( DBL_MAX )
{
   const G4int nvar = 6 ;
   G4int i; 
   for(i=0;i<nvar;i++)  {
     fYInSav[i]= DBL_MAX;
   }
}

G4ExactHelixStepper::~G4ExactHelixStepper() {} 

void
G4ExactHelixStepper::Stepper( const G4double yInput[],
		              const G4double*,
		                    G4double hstep,
		                    G4double yOut[],
		                    G4double yErr[]      )
{  
   const G4int nvar = 6 ;

   G4int i;
   // G4double      yTemp[7], yIn[7] ;
   G4ThreeVector Bfld_value;

   for(i=0;i<nvar;i++)  {
      // yIn[i]=     yInput[i];
      fYInSav[i]= yInput[i];
   }

   MagFieldEvaluate(yInput, Bfld_value) ;        
   fBfieldValue= Bfld_value;  // Save it for chord if needed.
   fLastStepSize = hstep;     //   ditto

   // DumbStepper(yIn, Bfld_value, hstep, yTemp);
   AdvanceHelix(yInput, Bfld_value, hstep, yOut);

   // We are assuming a constant field: helix is exact.
   for(i=0;i<nvar;i++) {
     yErr[i] = 0.0 ;
   }

   yInitialEHS = G4ThreeVector( yInput[0],   yInput[1],   yInput[2]); 
   yFinalEHS   = G4ThreeVector( yOut[0],  yOut[1],  yOut[2]); 
}

void
G4ExactHelixStepper::DumbStepper( const G4double  yIn[],
				   G4ThreeVector   Bfld,
				   G4double  h,
				   G4double  yOut[])
{
  // Assuming a constant field: solution is a helix
  AdvanceHelix(yIn, Bfld, h, yOut);

  G4Exception("G4ExactHelixStepper::DumbStepper should not be called.",
	      "EHS:NoDumbStepper", FatalException, "Stepper must do all the work." ); 
}  

G4double
G4ExactHelixStepper::DistChord() const 
{
  //  Method below is good only for < 2 pi   ---> TO-DO
   const G4int nvar = 6 ;

  // Calculate yMidPoint
   G4double hHalf = fLastStepSize * 0.5; 
   G4ThreeVector Bfld_initial= fBfieldValue; 
   G4double      yStart[7], yMid[7] ;

   G4int i;
   for(i=0;i<nvar;i++) yStart[i]= fYInSav[i];

   // Do the half step
   // DumbStepper(yStart, Bfld_initial, hHalf, yMid);
   AdvanceHelix(yStart, Bfld_initial, hHalf, yMid);
   G4ThreeVector yMidPointEHS = G4ThreeVector( yMid[0],  yMid[1],  yMid[2]); 

  return G4LineSection::Distline( yMidPointEHS, yInitialEHS, yFinalEHS );
  // This is a class method that gives distance of Mid 
  //  from the Chord between the Initial and Final points.
}

G4int
G4ExactHelixStepper::IntegratorOrder() const 
{
  return 1; 
}
