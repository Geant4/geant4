//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4MWeightWindowConfigurator.hh,v 1.5 2003/11/26 14:51:48 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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
