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
// $Id: G4PImportanceConfigurator.cc,v 1.2 2002-10-16 16:27:00 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4PImportanceConfigurator
//

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#include "G4PImportanceConfigurator.hh"

#include "G4ParallelWorld.hh"
#include "G4ImportanceAlgorithm.hh"

G4PImportanceConfigurator::
G4PImportanceConfigurator(const G4String &particlename,
			  G4ParallelWorld &parallelWorld,
			  G4VIStore &istore,
			  const G4VImportanceAlgorithm *ialg) 
  :
  fPlacer(particlename),
  fPWorld(parallelWorld),
  fDeleteIalg( ( ! ialg) ),
  fIalgorithm(( (fDeleteIalg) ? 
		new G4ImportanceAlgorithm : ialg)),
  fExaminer(*fIalgorithm, 
	    parallelWorld.GetParallelStepper(),  
	    istore),
  fParallelImportanceProcess(0)
{
  if (!fIalgorithm) {
    G4std::G4Exception("ERROR:G4PImportanceConfigurator::G4PImportanceConfigurator: no fIalgorithm!");
  }
}


G4PImportanceConfigurator::~G4PImportanceConfigurator(){
  if (fParallelImportanceProcess) {
    fPlacer.RemoveProcess(fParallelImportanceProcess);
    delete fParallelImportanceProcess;
  }
  if (fDeleteIalg) {
    delete fIalgorithm;
  }
}

void G4PImportanceConfigurator::
Configure(G4VSamplerConfigurator *preConf){
  const G4VTrackTerminator *terminator = 0;
  if (preConf) {
    terminator = preConf->GetTrackTerminator();
  };

  fParallelImportanceProcess = 
    new G4ParallelImportanceProcess(fExaminer, 
				    fPWorld.GetGeoDriver(), 
				    fPWorld.GetParallelStepper(),
				    terminator);
  if (!fParallelImportanceProcess) {
    G4std::G4Exception("ERROR:G4PImportanceConfigurator::Configure: new failed to create G4ParallelImportanceProcess!");
  }
  
  fPlacer.AddProcessAsSecondDoIt(fParallelImportanceProcess); 
}

const G4VTrackTerminator *G4PImportanceConfigurator::
GetTrackTerminator() const {
  return fParallelImportanceProcess;
}

