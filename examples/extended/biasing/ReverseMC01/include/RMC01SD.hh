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
/// \file biasing/ReverseMC01/include/RMC01SD.hh
/// \brief Definition of the RMC01SD class
//
//
//////////////////////////////////////////////////////////////
//  Class Name:            RMC01SD
//        Author:               L. Desorgher
//        Organisation:         SpaceIT GmbH
//        Contract:        ESA contract 21435/08/NL/AT
//        Customer:             ESA/ESTEC
//////////////////////////////////////////////////////////////
// CHANGE HISTORY
//--------------
//      ChangeHistory:
//                 17-11-2009 creation by L. Desorgher
//
//-------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RMC01SD_h
#define RMC01SD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4ThreeVector.hh"

class G4Step;
class G4HCofThisEvent;
class G4StepPoint;
class G4LogicalVolume;
class G4VPhysicalVolume;

#include"G4THitsCollection.hh"
#include"G4ios.hh"
#include"RMC01DoubleWithWeightHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RMC01SD : public G4VSensitiveDetector
{
  public:
    RMC01SD(G4String name);
    virtual ~RMC01SD();
    virtual void Initialize(G4HCofThisEvent*HCE);
    virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    virtual void EndOfEvent(G4HCofThisEvent*HCE);
    virtual void Clear();
    virtual void DrawAll();
    virtual void PrintAll();

  private:
    G4double fTotalEventEdep;
    RMC01DoubleWithWeightHitsCollection* fEventEdepCollection;
    RMC01DoubleWithWeightHitsCollection* fProtonCurrentCollection;
    RMC01DoubleWithWeightHitsCollection* fGammaCurrentCollection;
    RMC01DoubleWithWeightHitsCollection* fElectronCurrentCollection;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

