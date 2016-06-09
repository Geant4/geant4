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
// $Id: G4PImportanceConfigurator.cc,v 1.6.2.1 2006/06/29 21:11:22 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// ----------------------------------------------------------------------
// Class G4PImportanceConfigurator
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#include "G4PImportanceConfigurator.hh"

#include "G4ParallelWorld.hh"
#include "G4ImportanceAlgorithm.hh"
#include "G4ParallelImportanceProcess.hh"

G4PImportanceConfigurator::
G4PImportanceConfigurator(const G4String &particlename,
                                G4ParallelWorld &parallelWorld,
                                G4VIStore &istore,
                          const G4VImportanceAlgorithm *ialg) 
  : fPlacer(particlename),
    fPWorld(parallelWorld),
    fDeleteIalg( ( ! ialg) ),
    fIalgorithm(( (fDeleteIalg) ? new G4ImportanceAlgorithm : ialg)),
    fExaminer(*fIalgorithm, parallelWorld.GetParallelStepper(), istore),
    fParallelImportanceProcess(0)
{
  if (!fIalgorithm)
  {
    G4Exception("G4PImportanceConfigurator::G4PImportanceConfigurator()",
                "InvalidSetup", FatalException, "No algorithm specified !");
  }
}

G4PImportanceConfigurator::~G4PImportanceConfigurator()
{
  if (fParallelImportanceProcess)
  {
    fPlacer.RemoveProcess(fParallelImportanceProcess);
    delete fParallelImportanceProcess;
  }
  if (fDeleteIalg)
  {
    delete fIalgorithm;
  }
}

void
G4PImportanceConfigurator::Configure(G4VSamplerConfigurator *preConf)
{
  const G4VTrackTerminator *terminator = 0;
  if (preConf)
  {
    terminator = preConf->GetTrackTerminator();
  }

  fParallelImportanceProcess = 
    new G4ParallelImportanceProcess(fExaminer, 
                                    fPWorld.GetGeoDriver(), 
                                    fPWorld.GetParallelStepper(),
                                    terminator);
  if (!fParallelImportanceProcess)
  {
    G4Exception("G4PImportanceConfigurator::Configure()",
                "FatalError", FatalException,
                "Failed to allocate G4ParallelImportanceProcess !");
  }
  fPlacer.AddProcessAsSecondDoIt(fParallelImportanceProcess); 
}

const G4VTrackTerminator *
G4PImportanceConfigurator::GetTrackTerminator() const
{
  return fParallelImportanceProcess;
}

