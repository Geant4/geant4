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
// $Id: G4PWeightWindowConfigurator.hh,v 1.5 2006/06/29 21:10:06 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// ----------------------------------------------------------------------
// Class G4PWeightWindowConfigurator
//
// Class description:
// Configures the weight window processes for the parallel geometry.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4PWeightWindowConfigurator_hh
#define G4PWeightWindowConfigurator_hh G4PWeightWindowConfigurator_hh

#include "G4Types.hh"
#include "G4ProcessPlacer.hh"
#include "G4VSamplerConfigurator.hh"
#include "G4PlaceOfAction.hh"
#include "G4WeightWindowExaminer.hh"

class G4VWeightWindowStore;
class G4VWeightWindowAlgorithm;
class G4ParallelWorld;
class G4ParallelWeightWindowProcess;

class G4PWeightWindowConfigurator : public G4VSamplerConfigurator
{

public:  // with description

  G4PWeightWindowConfigurator(const G4String &particlename,
                              G4ParallelWorld &parallelWorld,
                              G4VWeightWindowStore &wwstore,
                              const G4VWeightWindowAlgorithm *wwAlg,
                              G4PlaceOfAction placeOfAction);

  virtual ~G4PWeightWindowConfigurator();
  virtual void Configure(G4VSamplerConfigurator *preConf);
  virtual const G4VTrackTerminator *GetTrackTerminator() const;

private:

  G4PWeightWindowConfigurator(const G4PWeightWindowConfigurator &);
  G4PWeightWindowConfigurator &
  operator=(const G4PWeightWindowConfigurator &);

  G4ProcessPlacer fPlacer;
  G4ParallelWorld &fPWorld;
  G4bool fDeleteWWalg;
  const G4VWeightWindowAlgorithm *fWWalgorithm;
  G4WeightWindowExaminer fExaminer;
  G4VProcess *fParallelWWProcess;
  G4PlaceOfAction fPlaceOfAction;
  G4VTrackTerminator *fTrackTerminator;
};

#endif
