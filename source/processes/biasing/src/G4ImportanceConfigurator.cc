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
// $Id$
//
// ----------------------------------------------------------------------
// Class G4ImportanceConfigurator
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#include "G4ImportanceConfigurator.hh"

#include "G4ImportanceProcess.hh"
#include "G4ProcessPlacer.hh"
#include "G4ImportanceAlgorithm.hh"

G4ImportanceConfigurator::
G4ImportanceConfigurator(G4VPhysicalVolume* worldvolume, 
			 const G4String &particlename,
                          G4VIStore &istore,
                          const G4VImportanceAlgorithm *ialg, G4bool para)
  : fWorld(worldvolume),
    fPlacer(particlename),
    fIStore(istore),
    fDeleteIalg( ( ! ialg) ),
    fIalgorithm(( (fDeleteIalg) ? 
                  new G4ImportanceAlgorithm : ialg)),
    fImportanceProcess(0),
    paraflag(para)
{
}

G4ImportanceConfigurator::~G4ImportanceConfigurator()
{
  if (fImportanceProcess)
  {
    fPlacer.RemoveProcess(fImportanceProcess);
    delete fImportanceProcess;
  }
  if (fDeleteIalg)
  {
    delete fIalgorithm;
  }
}

void  
G4ImportanceConfigurator::Configure(G4VSamplerConfigurator *preConf)
{
  G4cout << " entering importance configure, paraflag " << paraflag << G4endl;
  const G4VTrackTerminator *terminator = 0;
  if (preConf)
  {
    terminator = preConf->GetTrackTerminator();
  };

  G4cout << " creating importance process, paraflag is: " << paraflag << G4endl;
  fImportanceProcess = 
    new G4ImportanceProcess(*fIalgorithm, 
                                fIStore, 
                                terminator,"ImportanceProcess",paraflag);
  if (!fImportanceProcess)
  {
    G4Exception("G4ImportanceConfigurator::Configure()",
                "FatalError", FatalException,
                "Failed allocation of G4ImportanceProcess !");
  }
  
  if(paraflag) fImportanceProcess->SetParallelWorld(fWorld);

  fPlacer.AddProcessAsSecondDoIt(fImportanceProcess);
}

const G4VTrackTerminator *G4ImportanceConfigurator::
GetTrackTerminator() const
{
  return fImportanceProcess;
}
