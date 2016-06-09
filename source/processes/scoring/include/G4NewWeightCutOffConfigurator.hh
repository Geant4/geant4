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
// $Id: G4NewWeightCutOffConfigurator.hh,v 1.1 2006/11/20 10:02:05 ahoward Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// Class G4NewWeightCutOffConfigurator
//
// Class description:
// This class builds and places the G4NewWeightCutOffProcess.
// If the object is deleted the process is removed from the 
// process list.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4NewWeightCutOffConfigurator_hh
#define G4NewWeightCutOffConfigurator_hh G4NewWeightCutOffConfigurator_hh

#include "G4Types.hh"
#include "G4VSamplerConfigurator.hh"
#include "G4NewProcessPlacer.hh"

class G4NewWeightCutOffProcess;
class G4VGCellFinder;
class G4VIStore;
class G4VPhysicalVolume;

class G4NewWeightCutOffConfigurator : public G4VSamplerConfigurator
{

public:  // with description

  G4NewWeightCutOffConfigurator(G4VPhysicalVolume* worldvolume,
				const G4String &particlename,
                             G4double wsurvival,
                             G4double wlimit,
                             G4double isource,
                             G4VIStore *istore,
                             const G4VGCellFinder &aGCellFinder,G4bool paraflag);

  virtual ~G4NewWeightCutOffConfigurator();
  virtual void Configure(G4VSamplerConfigurator *preConf);
  virtual const G4VTrackTerminator *GetTrackTerminator() const ;
  
private:

  G4NewWeightCutOffConfigurator(const G4NewWeightCutOffConfigurator&);
  G4NewWeightCutOffConfigurator &
  operator=(const G4NewWeightCutOffConfigurator&);
  G4VPhysicalVolume* fWorld;
  G4NewProcessPlacer fPlacer;
  G4NewWeightCutOffProcess *fNewWeightCutOffProcess;
  G4bool fPlaced;


  G4bool paraflag;

};

#endif
