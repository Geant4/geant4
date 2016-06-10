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
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:    G4QMDCollision.hh 
//
//      Author: Koi, Tatsumi (tkoi@slac.stanford.edu)       
// 
//      Creation date: 9 April 2007
// -----------------------------------------------------------------------------
//
// 081120 Add deltaT in signature of CalKinematicsOfBinaryCollisions

#ifndef G4QMDCollision_hh
#define G4QMDCollision_hh

#include "G4QMDSystem.hh"
#include "G4QMDMeanField.hh"

#include "G4Scatterer.hh"

class G4QMDCollision 
{
   public:
      G4QMDCollision();
      ~G4QMDCollision();

      void CalKinematicsOfBinaryCollisions( G4double );
      G4bool CalFinalStateOfTheBinaryCollision( G4int , G4int );
      G4bool CalFinalStateOfTheBinaryCollisionJQMD( G4double , G4double , G4ThreeVector , G4double , G4double , G4ThreeVector , G4double , G4int , G4int );
      //     CalFinalStateOfTheBinaryCollision ( sig , cutoff , pcm , prcm , srt, beta , gamma , i , j );

      void SetMeanField ( G4QMDMeanField* meanfield ){ theMeanField = meanfield; theSystem = meanfield->GetSystem(); }

      // Get, Set method of parameters; For expert only 
      void deltar( G4double x ){ fdeltar = x; };
      G4double deltar(){ return fdeltar; };
      void bcmax0( G4double x ){ fbcmax0 = x; };
      G4double bcmax0(){ return fbcmax0; };
      void bcmax1( G4double x ){ fbcmax1 = x; };
      G4double bcmax1(){ return fbcmax1; };
      void epse( G4double x ){ fepse = x; };
      G4double epse(){ return fepse; };

   private:
      //copy is unexpeced
      G4QMDCollision( const G4QMDCollision& ){;};
      const G4QMDCollision& operator= ( const G4QMDCollision& );

      G4QMDSystem* theSystem;
      G4QMDMeanField* theMeanField;

      G4double fdeltar;
      G4double fbcmax0 , fbcmax1;
      //G4double sig0 , sig1;
      G4double fepse; 

      G4Scatterer* theScatterer;
};

#endif
