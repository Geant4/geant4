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
/// \file radiobiology/include/DetectorMessenger.hh
/// \brief Definition of the RadioBio::DetectorMessenger class

#ifndef RadiobiologyDetectorMessenger_h
#define RadiobiologyDetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;

namespace RadioBio
{

// Forward declariation of other radiobiology classes
class DetectorConstruction;

class DetectorMessenger : public G4UImessenger
{
  public:
    DetectorMessenger(DetectorConstruction*);
    ~DetectorMessenger() override;

    void SetNewValue(G4UIcommand*, G4String) override;

  private:
    DetectorConstruction* fDetector = nullptr;

    // Geometry directory
    G4UIdirectory* fGeometryDir = nullptr;

    // To change world material
    G4UIcmdWithAString* fMaterCmd = nullptr;

    // To change world size
    G4UIcmdWithADoubleAndUnit* fSizeCmd = nullptr;

    // To change world size (giving a vector)
    G4UIcmdWith3VectorAndUnit* fSizeVectorCmd = nullptr;

    // To change world X size
    G4UIcmdWithADoubleAndUnit* fSizeXCmd = nullptr;

    // To change world Y size
    G4UIcmdWithADoubleAndUnit* fSizeYCmd = nullptr;

    // To change world Z size
    G4UIcmdWithADoubleAndUnit* fSizeZCmd = nullptr;
};
}  // namespace RadioBio

#endif  // DetectorMessenger_h
