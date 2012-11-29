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
// $Id: 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// --------------------------------------------------------------
//   New Physics scheme           8 Jan. 1997  H.Kurahige
// ------------------------------------------------------------

#include "G4IVContinuousDiscreteProcess.hh"
G4IVContinuousDiscreteProcess::G4IVContinuousDiscreteProcess()
  :G4VProcess("No Name Discrete Process"),
   valueGPILSelection(CandidateForSelection), 
   theNlambdaTable(0),theInverseNlambdaTable(0),
   BIGSTEP(1.e10)
{
  G4Exception("G4IVContinuousDiscreteProcess::G4IVContinuousDiscreteProcess","ProcMan102",
	      JustWarning,"Default constructor is called");
}

G4IVContinuousDiscreteProcess::G4IVContinuousDiscreteProcess(const G4String& aName , G4ProcessType aType)
  : G4VProcess(aName, aType),
    valueGPILSelection(CandidateForSelection), 
    theNlambdaTable(0),theInverseNlambdaTable(0),
    BIGSTEP(1.e10)
{
  enableAtRestDoIt = false;
}
G4IVContinuousDiscreteProcess::~G4IVContinuousDiscreteProcess()
{
}

G4IVContinuousDiscreteProcess::G4IVContinuousDiscreteProcess(G4IVContinuousDiscreteProcess& right)
  : G4VProcess(right),
    valueGPILSelection(right.valueGPILSelection),
    theNlambdaTable(0),theInverseNlambdaTable(0),
    BIGSTEP(right.BIGSTEP)
{
}


G4double
G4IVContinuousDiscreteProcess::PostStepGetPhysicalInteractionLength(
                              const G4Track&,
                              G4double,
                              G4ForceCondition*
                             )
{
  G4double value = DBL_MAX ;

  return value;
}







