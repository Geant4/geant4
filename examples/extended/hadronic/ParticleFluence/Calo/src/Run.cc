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
/// \file Run.cc
/// \brief Implementation of the Run class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run() : G4Run(), fNumEvents( 0 ),
             fPrimaryParticleId( 0 ), fPrimaryParticleEnergy( 0.0 ),
             fPrimaryParticleDirection( G4ThreeVector( 0.0, 0.0, 0.0 ) ),
             fAbsorberMaterialName( "" ), fActiveMaterialName( "" ),
             fCubicVolumeScoringUpDown( 1.0 ), fCubicVolumeScoringSide( 1.0 )
{
  fSteppingArray.fill( 0.0 );
  fTrackingArray1.fill( 0 );
  fTrackingArray2.fill( 0.0 );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::RecordEvent( const G4Event* anEvent ) {
  // This method is called automatically by the Geant4 kernel (not by the user!) at the end
  // of each event : in MT-mode, it is called only for the working thread that handled the event.
  G4int nEvt = anEvent->GetEventID();
  if ( nEvt % 10 == 0 ) G4cout << " Event#=" << nEvt << G4endl;
  G4Run::RecordEvent( anEvent );  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge( const G4Run* aRun ) {
  // This method is called automatically by the Geant4 kernel (not by the user!) only in the case
  // of multithreaded mode and only for working threads.
  const Run* localRun = static_cast< const Run* >( aRun );
  fPrimaryParticleId = localRun->GetPrimaryParticleId();
  fPrimaryParticleEnergy = localRun->GetPrimaryParticleEnergy();
  fPrimaryParticleDirection = localRun->GetPrimaryParticleDirection();
  fAbsorberMaterialName = localRun->GetAbsorberMaterialName();
  fActiveMaterialName = localRun->GetActiveMaterialName();
  fCubicVolumeScoringUpDown = localRun->GetCubicVolumeScoringUpDown();
  fCubicVolumeScoringSide = localRun->GetCubicVolumeScoringSide();
  fNumEvents += localRun->GetNumberOfEvent();
  for ( G4int i = 0; i < SteppingAction::fkNumberCombinations; ++i ) {
    fSteppingArray[i] += localRun->GetSteppingArray()[i];
  }
  for ( G4int i = 0; i < TrackingAction::fkNumberCombinations; ++i ) {
    fTrackingArray1[i] += localRun->GetTrackingArray1()[i];
    fTrackingArray2[i] += localRun->GetTrackingArray2()[i];
  }
  G4Run::Merge( aRun );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::PrintInfo() const {
  // This method is called by RunAction::EndOfRunAction. In MT-mode, only the master thread
  // calls it.
  const G4double floatingNumberOfEvents =
    std::max( 1.0, fNumEvents > 0 ? fNumEvents*1.0 : GetNumberOfEvent()*1.0 );
  // The fluence in the scoring volume is defined as sum of step lengths in that volume
  // divided by the volume of that scoring volume.
  const G4double conversionFactor = CLHEP::cm * CLHEP::cm;  // From  mm^-2  to  cm^-2
  const G4double factorUpDown =
    conversionFactor / ( fCubicVolumeScoringUpDown*floatingNumberOfEvents );
  const G4double factorSide =
    conversionFactor / ( fCubicVolumeScoringSide*floatingNumberOfEvents );
  G4cout << std::setprecision(6) << G4endl << G4endl
         << " ===============  Run::PrintInfo()  =============== \t RunID = " << GetRunID()
         << G4endl
         << " Primary particle PDG code = " << fPrimaryParticleId << G4endl
         << " Primary particle kinetic energy = " << fPrimaryParticleEnergy / CLHEP::GeV
         << " GeV" << G4endl
         << " Primary particle direction = " << fPrimaryParticleDirection << G4endl
         << " Absorber material = " << fAbsorberMaterialName << G4endl
         << " Active material   = " << fActiveMaterialName << G4endl
         << " Cubic-volume scoring up-down = " << fCubicVolumeScoringUpDown << " mm^3" << G4endl
         << " Cubic-volume scoring side    = " << fCubicVolumeScoringSide   << " mm^3" << G4endl
         << " Number of events = " << floatingNumberOfEvents << G4endl
         << " Conversion factor: fluence from mm^-2 to cm^-2 = " << conversionFactor << G4endl
         << " Particle fluence in unit of cm^-2 :" << G4endl;
  for ( G4int i = 0; i < SteppingAction::fkNumberScoringVolumes; ++i ) {
    G4double factor = ( i == 1 ? factorSide : factorUpDown );
    for ( G4int j = 0; j < SteppingAction::fkNumberKinematicRegions; ++j ) {
      for ( G4int k = 0; k < SteppingAction::fkNumberParticleTypes; ++k ) {
        G4int index = SteppingAction::GetIndex( i, j, k );
        //G4cout << "(i, j, k)=(" << i << ", " << j << ", " << k << ")  ->" << index;
        G4cout << "   case=" << std::setw(3) << index
               << "   " << std::setw(12) << SteppingAction::fkArrayScoringVolumeNames[i]
               << "   " << std::setw(12) << SteppingAction::fkArrayKinematicRegionNames[j]
               << "   " << std::setw(12) << SteppingAction::fkArrayParticleTypeNames[k]
               << "   " << std::setw( 8) << factor*fSteppingArray[index] << G4endl;
      }
    }
  }
  G4cout << " ------------------------------------------------------------- " << G4endl
         << " Extra information: particle production \t \t      <N>      <E_kin> <Sum_Ekin> [MeV]"
         << G4endl;
  const G4double normalization = 1.0 / floatingNumberOfEvents;
  for ( G4int i = 0; i < TrackingAction::fkNumberScoringVolumes; ++i ) {
    for ( G4int j = 0; j < TrackingAction::fkNumberKinematicRegions; ++j ) {
      for ( G4int k = 0; k < TrackingAction::fkNumberParticleTypes; ++k ) {
        G4int index = TrackingAction::GetIndex( i, j, k );
        //G4cout << "(i, j, k)=(" << i << ", " << j << ", " << k << ")  ->" << index;
        G4cout << "   case=" << std::setw(3) << index
               << "   " << std::setw(12) << TrackingAction::fkArrayScoringVolumeNames[i]
               << "   " << std::setw(12) << TrackingAction::fkArrayKinematicRegionNames[j]
               << "   " << std::setw(12) << TrackingAction::fkArrayParticleTypeNames[k]
               << "   " << std::setw( 8) << normalization * fTrackingArray1[index]
               << "   " << std::setw( 8) << ( fTrackingArray1[index] > 0 ?
                                              fTrackingArray2[index] / fTrackingArray1[index] :
                                              0.0 )
               << "   " << std::setw( 8) << normalization * fTrackingArray2[index]
               << G4endl;
      }
    }
  }  
  G4cout << " ============================================================= " << G4endl << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetSteppingArray( const std::array< G4double,
                            SteppingAction::fkNumberCombinations >& inputArray ) {
  for ( G4int i = 0; i < SteppingAction::fkNumberCombinations; ++i ) {
    fSteppingArray[i] = inputArray[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetTrackingArray1( const std::array< G4int,
                             TrackingAction::fkNumberCombinations >& inputArray ) {
  for ( G4int i = 0; i < TrackingAction::fkNumberCombinations; ++i ) {
    fTrackingArray1[i] = inputArray[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::SetTrackingArray2( const std::array< G4double,
                             TrackingAction::fkNumberCombinations >& inputArray ) {
  for ( G4int i = 0; i < TrackingAction::fkNumberCombinations; ++i ) {
    fTrackingArray2[i] = inputArray[i];
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
