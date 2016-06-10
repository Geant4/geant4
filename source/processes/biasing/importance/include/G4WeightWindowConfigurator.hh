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
// $Id: G4WeightWindowConfigurator.hh 77477 2013-11-25 09:42:24Z gcosmo $
//
// ----------------------------------------------------------------------
// Class G4WeightWindowConfigurator
//
// Class description:
// Configuration of weight window processes.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4WeightWindowConfigurator_hh
#define G4WeightWindowConfigurator_hh G4WeightWindowConfigurator_hh

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

public:  // with description

  G4WeightWindowConfigurator(const G4VPhysicalVolume* worldvolume,
			     const G4String &particlename,
                              G4VWeightWindowStore &wwstore,
                              const G4VWeightWindowAlgorithm *wwAlg,
                              G4PlaceOfAction placeOfAction,
			     G4bool paraflag);

  virtual ~G4WeightWindowConfigurator();
  virtual void Configure(G4VSamplerConfigurator *preConf);
  virtual const G4VTrackTerminator *GetTrackTerminator() const;

private:

  G4WeightWindowConfigurator(const G4WeightWindowConfigurator &);
  G4WeightWindowConfigurator &
  operator=(const G4WeightWindowConfigurator &);

  const G4VPhysicalVolume* fWorld;
  G4ProcessPlacer fPlacer;
  G4VWeightWindowStore &fWeightWindowStore;
  G4bool fDeleteWWalg;
  const G4VWeightWindowAlgorithm *fWWalgorithm;
  G4WeightWindowProcess *fWeightWindowProcess;
  G4PlaceOfAction fPlaceOfAction;


  G4bool paraflag;

};

#endif
