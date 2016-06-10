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
// $Id: G4WeightCutOffConfigurator.hh 77477 2013-11-25 09:42:24Z gcosmo $
//
// ----------------------------------------------------------------------
// Class G4WeightCutOffConfigurator
//
// Class description:
// This class builds and places the G4WeightCutOffProcess.
// If the object is deleted the process is removed from the 
// process list.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4WeightCutOffConfigurator_hh
#define G4WeightCutOffConfigurator_hh G4WeightCutOffConfigurator_hh

#include "G4Types.hh"
#include "G4VSamplerConfigurator.hh"
#include "G4ProcessPlacer.hh"

class G4WeightCutOffProcess;
//class G4VGCellFinder;
class G4VIStore;
class G4VPhysicalVolume;

class G4WeightCutOffConfigurator : public G4VSamplerConfigurator
{

public:  // with description

  G4WeightCutOffConfigurator(const G4VPhysicalVolume* worldvolume,
				const G4String &particlename,
                             G4double wsurvival,
                             G4double wlimit,
                             G4double isource,
                             G4VIStore *istore,
                             //const G4VGCellFinder &aGCellFinder,
			     G4bool paraflag);

  virtual ~G4WeightCutOffConfigurator();
  virtual void Configure(G4VSamplerConfigurator *preConf);
  virtual const G4VTrackTerminator *GetTrackTerminator() const ;
  
private:

  G4WeightCutOffConfigurator(const G4WeightCutOffConfigurator&);
  G4WeightCutOffConfigurator &
  operator=(const G4WeightCutOffConfigurator&);
  const G4VPhysicalVolume* fWorld;
  G4ProcessPlacer fPlacer;
  G4WeightCutOffProcess *fWeightCutOffProcess;
  G4bool fPlaced;


  G4bool paraflag;

};

#endif
