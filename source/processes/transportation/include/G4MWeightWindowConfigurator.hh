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
// $Id: G4MWeightWindowConfigurator.hh,v 1.6 2006/06/29 21:09:50 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// Class G4MWeightWindowConfigurator
//
// Class description:
// Configuration of weight window processes.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4MWeightWindowConfigurator_hh
#define G4MWeightWindowConfigurator_hh G4MWeightWindowConfigurator_hh

#include "G4Types.hh"
#include "G4ProcessPlacer.hh"
#include "G4VSamplerConfigurator.hh"
#include "G4PlaceOfAction.hh"

class G4VWeightWindowStore;
class G4VWeightWindowAlgorithm;
class G4MassWeightWindowProcess;

class G4MWeightWindowConfigurator : public G4VSamplerConfigurator
{

public:  // with description

  G4MWeightWindowConfigurator(const G4String &particlename,
                              G4VWeightWindowStore &wwstore,
                              const G4VWeightWindowAlgorithm *wwAlg,
                              G4PlaceOfAction placeOfAction);

  virtual ~G4MWeightWindowConfigurator();
  virtual void Configure(G4VSamplerConfigurator *preConf);
  virtual const G4VTrackTerminator *GetTrackTerminator() const;

private:

  G4MWeightWindowConfigurator(const G4MWeightWindowConfigurator &);
  G4MWeightWindowConfigurator &
  operator=(const G4MWeightWindowConfigurator &);

  G4ProcessPlacer fPlacer;
  G4VWeightWindowStore &fWeightWindowStore;
  G4bool fDeleteWWalg;
  const G4VWeightWindowAlgorithm *fWWalgorithm;
  G4MassWeightWindowProcess *fMassWeightWindowProcess;
  G4PlaceOfAction fPlaceOfAction;
};

#endif
