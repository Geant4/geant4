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

      void CalKinematicsOfBinaryCollisions();
      G4bool CalFinalStateOfTheBinaryCollision( G4int , G4int );
      G4bool CalFinalStateOfTheBinaryCollisionJQMD( G4double , G4double , G4ThreeVector , G4double , G4double , G4ThreeVector , G4double , G4int , G4int );
      //     CalFinalStateOfTheBinaryCollision ( sig , cutoff , pcm , prcm , srt, beta , gamma , i , j );

      void SetMeanField ( G4QMDMeanField* meanfield ){ theMeanField = meanfield; theSystem = meanfield->GetSystem(); }


   private:
      G4QMDSystem* theSystem;
      G4QMDMeanField* theMeanField;

      G4double deltar;
      G4double bcmax0 , bcmax1;
      G4double sig0 , sig1;
      G4double epse; 

      G4Scatterer* theScatterer;
};

#endif
