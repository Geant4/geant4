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
// --------------------------------------------------------------------------------
//       MONTE CARLO SIMULATION OF REALISTIC GEOMETRY FROM MICROSCOPES IMAGES
//
// Authors and contributors:
// P. Barberet, S. Incerti, N. H. Tran, L. Morelli
//
// University of Bordeaux, CNRS, LP2i, UMR5797, Gradignan, France
//
// If you use this code, please cite the following publication:
// P. Barberet et al.,
// "Monte-Carlo dosimetry on a realistic cell monolayer
// geometry exposed to alpha particles."
// Ph. Barberet et al 2012 Phys. Med. Biol. 57 2189
// doi: 110.1088/0031-9155/57/8/2189
// --------------------------------------------------------------------------------

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class DetectorConstruction;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorMessenger: public G4UImessenger
{
  public:

    explicit DetectorMessenger(DetectorConstruction*);
    ~DetectorMessenger() override;

    void SetNewValue(G4UIcommand*, G4String) override;

  private:

    DetectorConstruction*      fDetector = nullptr;

    G4UIdirectory*             fPhantomDir = nullptr;
    G4UIdirectory*             fWorldDir = nullptr;

    G4UIcmdWithAString*        fNameCmd = nullptr;
    G4UIcmdWithAString*        fMatCmd = nullptr;

    G4UIcmdWithADoubleAndUnit* fDenRedCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fDenGreenCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fDenBlueCmd = nullptr;

    G4UIcmdWithADoubleAndUnit* fShiftXCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fShiftYCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fShiftZCmd = nullptr;

    G4UIcmdWithADoubleAndUnit* fMediumSizeXYCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fMediumSizeZCmd = nullptr;

    G4UIcmdWithADoubleAndUnit* fWorldSizeXYCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fWorldSizeZCmd = nullptr;
};

#endif
