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
// $Id: G4PImportanceConfigurator.cc,v 1.1 2002-10-10 13:25:31 dressel Exp $
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
{}


G4PImportanceConfigurator::~G4PImportanceConfigurator(){
  if (fParallelImportanceProcess) {
    fPlacer.RemoveProcess(fParallelImportanceProcess);
    delete fParallelImportanceProcess;
  }
  if (fDeleteIalg) delete fIalgorithm;
}

void G4PImportanceConfigurator::
Configure(G4VSamplerConfigurator *preConf){
  G4VTrackTerminator *terminator = 0;
  if (preConf) {
    terminator = preConf->GetTrackTerminator();
  };

  fParallelImportanceProcess = 
    new G4ParallelImportanceProcess(fExaminer, 
				    fPWorld.GetGeoDriver(), 
				    fPWorld.GetParallelStepper(),
				    terminator);
  
  fPlacer.AddProcessAsSecondDoIt(fParallelImportanceProcess); 
}

G4VTrackTerminator *G4PImportanceConfigurator::GetTrackTerminator(){
  return fParallelImportanceProcess;
}

