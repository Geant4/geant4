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
// $Id: G4MImportanceConfigurator.cc,v 1.3 2002-11-04 10:47:56 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4MImportanceConfigurator
//

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------


#include "G4MImportanceConfigurator.hh"

#include "G4MassImportanceProcess.hh"
#include "G4ProcessPlacer.hh"
#include "G4ImportanceAlgorithm.hh"

G4MImportanceConfigurator::
G4MImportanceConfigurator(const G4String &particlename,
			  G4VIStore &istore,
			  const G4VImportanceAlgorithm *ialg)
  :
  fPlacer(particlename),
  fIStore(istore),
  fDeleteIalg( ( ! ialg) ),
  fIalgorithm(( (fDeleteIalg) ? 
		new G4ImportanceAlgorithm : ialg)),
  fMassImportanceProcess(0)
{}

G4MImportanceConfigurator::
~G4MImportanceConfigurator(){
  if (fMassImportanceProcess) {
    fPlacer.RemoveProcess(fMassImportanceProcess);
    delete fMassImportanceProcess;
  }
  if (fDeleteIalg) {
    delete fIalgorithm;
  }
}
void  
G4MImportanceConfigurator::Configure(G4VSamplerConfigurator *preConf){
  const G4VTrackTerminator *terminator = 0;
  if (preConf) {
    terminator = preConf->GetTrackTerminator();
  };

  fMassImportanceProcess = 
    new G4MassImportanceProcess(*fIalgorithm, 
				fIStore, 
				terminator);
  if (!fMassImportanceProcess) {
    G4Exception("ERROR: G4MImportanceConfigurator::Configure: new failed to create  G4MassImportanceProcess!");
  }
  fPlacer.AddProcessAsSecondDoIt(fMassImportanceProcess);
}

const G4VTrackTerminator *G4MImportanceConfigurator::
GetTrackTerminator() const {
  return fMassImportanceProcess;
}
