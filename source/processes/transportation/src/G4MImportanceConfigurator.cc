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
// $Id: G4MImportanceConfigurator.cc,v 1.7 2006/11/13 16:13:50 japost Exp $
// GEANT4 tag $Name: geant4-08-02 $
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
  : fPlacer(particlename),
    fIStore(istore),
    fDeleteIalg( ( ! ialg) ),
    fIalgorithm(( (fDeleteIalg) ? 
                  new G4ImportanceAlgorithm : ialg)),
    fMassImportanceProcess(0)
{
}

G4MImportanceConfigurator::~G4MImportanceConfigurator()
{
  if (fMassImportanceProcess)
  {
    fPlacer.RemoveProcess(fMassImportanceProcess);
    delete fMassImportanceProcess;
  }
  if (fDeleteIalg)
  {
    delete fIalgorithm;
  }
}

void  
G4MImportanceConfigurator::Configure(G4VSamplerConfigurator *preConf)
{
  const G4VTrackTerminator *terminator = 0;
  if (preConf)
  {
    terminator = preConf->GetTrackTerminator();
  };

  fMassImportanceProcess = 
    new G4MassImportanceProcess(*fIalgorithm, 
                                fIStore, 
                                terminator);
  if (!fMassImportanceProcess)
  {
    G4Exception("G4MImportanceConfigurator::Configure()",
                "FatalError", FatalException,
                "Failed allocation of G4MassImportanceProcess !");
  }
  fPlacer.AddProcessAsSecondDoIt(fMassImportanceProcess);
}

const G4VTrackTerminator *G4MImportanceConfigurator::
GetTrackTerminator() const
{
  return fMassImportanceProcess;
}
