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
// $Id: G4EMBuilder.cc,v 1.3 2005-11-29 16:54:48 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EMBuilder
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard
//
//----------------------------------------------------------------------------
//

#include "G4EMBuilder.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"

G4EMBuilder::G4EMBuilder(): wasActivated(false), synchOn(false),
  gammNucOn(false), theElectronSynch(0), thePositronSynch(0),
  theGNPhysics(0)
{
  theT = new G4EMTailorer(this);
}

G4EMBuilder::~G4EMBuilder()
{
  delete theT;
  if(theElectronSynch) delete theElectronSynch;
  if(thePositronSynch) delete thePositronSynch;
  if(theGNPhysics)     delete theGNPhysics;
}

void G4EMBuilder::Synch (G4String & newState)
{
  if(newState == "on") synchOn = true;
  else                 synchOn = false;
}

void G4EMBuilder::GammaNuclear (G4String & newState)
{
  if(newState == "on") gammNucOn = true;
  else                 gammNucOn = false;
}

void G4EMBuilder::Build()
{
  G4ProcessManager * pManager = 0;
  wasActivated = true;

  if(synchOn) {
    pManager = G4Electron::Electron()->GetProcessManager();
    theElectronSynch = new G4SynchrotronRadiation();
    pManager->AddDiscreteProcess(theElectronSynch);

    pManager = G4Positron::Positron()->GetProcessManager();
    thePositronSynch = new G4SynchrotronRadiation();
    pManager->AddDiscreteProcess(thePositronSynch);
  }
  if (gammNucOn) {
    theGNPhysics = new G4ElectroNuclearBuilder();
    theGNPhysics->Build();
  }
}
