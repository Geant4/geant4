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
// $Id: G4HelixExplicitEuler.cc 97598 2016-06-06 07:19:46Z gcosmo $
//
//
//  Helix Explicit Euler: x_1 = x_0 + helix(h)
//  with helix(h) being a helix piece of length h
//  most simple approach for solving linear differential equations.
//  Take the current derivative and add it to the current position.
//
//  W.Wander <wwc@mit.edu> 12/09/97 
// -------------------------------------------------------------------

#include "G4HelixExplicitEuler.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"


void G4HelixExplicitEuler::Stepper(  const G4double  yInput[7],
                               const G4double*,
                                     G4double Step,
                                     G4double yOut[7],
                                     G4double yErr[])

{

 //Estimation of the Stepping Angle

  G4ThreeVector Bfld;
  MagFieldEvaluate(yInput, Bfld); 
  
  const G4int nvar = 6 ;
  G4int i;
  G4double      yTemp[8], yIn[8] ;
  G4ThreeVector  Bfld_midpoint;
  //  Saving yInput because yInput and yOut can be aliases for same array
        for(i=0;i<nvar;i++) yIn[i]=yInput[i];
     
        G4double h = Step * 0.5;
 
     // Do full step and two half steps
        G4double yTemp2[7];
        AdvanceHelix(yIn,   Bfld,  h, yTemp2,yTemp);
        MagFieldEvaluate(yTemp2, Bfld_midpoint) ;     
        AdvanceHelix(yTemp2, Bfld_midpoint, h, yOut);
    
     // Error estimation
        for(i=0;i<nvar;i++) {
         yErr[i] = yOut[i] - yTemp[i] ;
       }
    
}

G4double G4HelixExplicitEuler::DistChord()   const 
{
  // Implementation : must check whether h/R > 2 pi  !!
  //   If( h/R <  pi) use G4LineSection::DistLine
  //   Else           DistChord=R_helix
  //
  G4double distChord;
  G4double Ang_curve=GetAngCurve();

      
	 if(Ang_curve<=pi){
	   distChord=GetRadHelix()*(1-std::cos(0.5*Ang_curve));
	 }
         else 
         if(Ang_curve<twopi){
           distChord=GetRadHelix()*(1+std::cos(0.5*(twopi-Ang_curve)));
         }
         else{
          distChord=2.*GetRadHelix();  
         }

  return distChord;
  
}
void
G4HelixExplicitEuler::DumbStepper( const G4double  yIn[],
				   G4ThreeVector   Bfld,
				   G4double        h,
				   G4double        yOut[])
{
    
       AdvanceHelix(yIn, Bfld, h, yOut);
               
}  
