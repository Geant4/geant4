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
// G4WeightWindowConfigurator
//
// Class description:
//
// Configuration of weight window processes.
//
// Author: Michael Dressel, CERN
// --------------------------------------------------------------------
#ifndef G4WeightWindowConfigurator_hh
#define G4WeightWindowConfigurator_hh 1

#include "G4Types.hh"
#include "G4ProcessPlacer.hh"
#include "G4VSamplerConfigurator.hh"
#include "G4PlaceOfAction.hh"

class G4VWeightWindowStore;
class G4VWeightWindowAlgorithm;
class G4WeightWindowProcess;
class G4VPhysicalVolume;

class G4WeightWindowConfigurator : public G4VSamplerConfigurator
{

 public:

  G4WeightWindowConfigurator(const G4VPhysicalVolume* worldvolume,
			     const G4String& particlename,
                                   G4VWeightWindowStore& wwstore,
                             const G4VWeightWindowAlgorithm* wwAlg,
                                   G4PlaceOfAction placeOfAction,
                                   G4bool paraflag);

  virtual ~G4WeightWindowConfigurator();

  G4WeightWindowConfigurator(const G4WeightWindowConfigurator &) = delete;
  G4WeightWindowConfigurator& operator=(const G4WeightWindowConfigurator &) = delete;

  virtual void Configure(G4VSamplerConfigurator* preConf);
  virtual const G4VTrackTerminator* GetTrackTerminator() const;

 private:

  const G4VPhysicalVolume* fWorld = nullptr;
  G4ProcessPlacer fPlacer;
  G4VWeightWindowStore& fWeightWindowStore;
  G4bool fDeleteWWalg = false;
  const G4VWeightWindowAlgorithm* fWWalgorithm = nullptr;
  G4WeightWindowProcess* fWeightWindowProcess = nullptr;
  G4PlaceOfAction fPlaceOfAction;
  G4bool paraflag = false;
};

#endif
