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
// $Id: G4MassImportanceScoreManager.cc,v 1.4 2002-04-09 17:40:15 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4MassImportanceScoreManager.cc
//
// ----------------------------------------------------------------------

#include "G4MassImportanceScoreManager.hh"
#include "G4MassImportanceManager.hh"
#include "G4MassScoreManager.hh"


G4MassImportanceScoreManager::
G4MassImportanceScoreManager(G4VIStore &aIstore,
			     G4VPScorer &ascorer,
			     const G4String &particlename)
 : fMassImportanceManager(new G4MassImportanceManager(aIstore, particlename)),
   fMassScoreManager(new G4MassScoreManager(ascorer, particlename))
{}
  
G4MassImportanceScoreManager::
G4MassImportanceScoreManager(G4VIStore &aIstore,
			     G4VPScorer &ascorer,
			     const G4String &particlename,
			     const G4VImportanceAlgorithm &algorithm)
 : fMassImportanceManager(new G4MassImportanceManager(aIstore, 
                                                      particlename,
                                                      algorithm)),
   fMassScoreManager(new G4MassScoreManager(ascorer, particlename))
{}

G4MassImportanceScoreManager::~G4MassImportanceScoreManager()
{
  delete fMassScoreManager;
  delete fMassImportanceManager;
}

void G4MassImportanceScoreManager::Initialize()
{
  fMassScoreManager->Initialize();
  fMassImportanceManager->Initialize();
}
