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
// $Id: G4PWeightWindowConfigurator.hh,v 1.1 2003-08-19 15:18:22 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4PWeightWindowConfigurator
//
// Class description:


// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4PWeightWindowConfigurator_hh
#define G4PWeightWindowConfigurator_hh G4PWeightWindowConfigurator_hh

#include "globals.hh"
#include "G4ProcessPlacer.hh"
#include "G4VSamplerConfigurator.hh"
#include "G4PlaceOfAction.hh"
#include "G4ParallelWeightWindowProcess.hh"
#include "G4WeightWindowExaminer.hh"

class G4VWeightWindowStore;
class G4VWeightWindowAlgorithm;
class G4ParallelWorld;

class G4PWeightWindowConfigurator : public G4VSamplerConfigurator{
public:
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
