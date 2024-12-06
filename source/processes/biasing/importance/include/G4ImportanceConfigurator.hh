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
// Class G4ImportanceConfigurator
//
// Class description:
//
// This class builds and places the G4ImportanceProcess.
// If the object is deleted the process is removed from the 
// process list.
//
// Author: Michael Dressel, CERN
// --------------------------------------------------------------------
#ifndef G4ImportanceConfigurator_hh
#define G4ImportanceConfigurator_hh 1

#include "G4Types.hh"
#include "G4ProcessPlacer.hh"
#include "G4VSamplerConfigurator.hh"

class G4ImportanceProcess;
class G4VImportanceAlgorithm;
class G4VIStore;
class G4VPhysicalVolume;

class G4ImportanceConfigurator : public G4VSamplerConfigurator
{

 public:

  G4ImportanceConfigurator(const G4VPhysicalVolume* worldvolume, 
			   const G4String& particlename,
                                 G4VIStore& istore,
                           const G4VImportanceAlgorithm* ialg,
			         G4bool paraflag);

  G4ImportanceConfigurator(const G4String& worldvolumeName, 
			   const G4String& particlename,
                                 G4VIStore& istore,
                           const G4VImportanceAlgorithm* ialg,
			         G4bool paraflag);

  virtual ~G4ImportanceConfigurator();

  G4ImportanceConfigurator(const G4ImportanceConfigurator&) = delete;
  G4ImportanceConfigurator& operator=(const G4ImportanceConfigurator&) = delete;

  virtual void Configure(G4VSamplerConfigurator* preConf);
  virtual const G4VTrackTerminator* GetTrackTerminator() const;

  void SetWorldName(const G4String& Name);

 private:

  const G4VPhysicalVolume* fWorld = nullptr;
  G4String fWorldName;
  G4ProcessPlacer fPlacer;
  G4VIStore& fIStore;
  G4bool fDeleteIalg = false;
  const G4VImportanceAlgorithm* fIalgorithm = nullptr;
  G4ImportanceProcess* fImportanceProcess = nullptr;
  G4bool paraflag = false;
};

#endif
