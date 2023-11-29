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
// G4VUserDetectorConstruction
//
// Class description:
//
// This is the abstract base class for the user's mandatory initialization
// of the detector setup. It has only one pure virtual method Construct()
// which is invoked by G4RunManager when its Initialize() method is invoked.
// The Construct() method must return the G4VPhysicalVolume pointer which
// represents the world volume.

// Original author: M.Asai, 1999
// --------------------------------------------------------------------
#ifndef G4VUserDetectorConstruction_hh
#define G4VUserDetectorConstruction_hh 1

#include "globals.hh"

#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VUserParallelWorld;
class G4VSensitiveDetector;

class G4VUserDetectorConstruction
{
  public:
    G4VUserDetectorConstruction() = default;
    virtual ~G4VUserDetectorConstruction() = default;

    virtual G4VPhysicalVolume* Construct() = 0;

    // This method is used in multi-threaded applications to build
    // per-worker non-shared objects: SensitiveDetectors and Field managers.
    virtual void ConstructSDandField();

    virtual void CloneSD();
    virtual void CloneF();

    void RegisterParallelWorld(G4VUserParallelWorld*);

    G4int ConstructParallelGeometries();
    void ConstructParallelSD();

    G4int GetNumberOfParallelWorld() const;
    G4VUserParallelWorld* GetParallelWorld(G4int i) const;

  protected:
    void SetSensitiveDetector(const G4String& logVolName, G4VSensitiveDetector* aSD,
                              G4bool multi = false);
    void SetSensitiveDetector(G4LogicalVolume* logVol, G4VSensitiveDetector* aSD);

  private:
    std::vector<G4VUserParallelWorld*> parallelWorld;
};

#endif
