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
// $Id: processTest.cc,v 1.5 2001-10-29 10:44:11 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
//      File name:     processTest
//
//      Author:        Maria Grazia Pia
// 
//      Creation date: 1 May 2001
//
//      Modifications: 
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/fstream"
#include "g4std/iomanip"

#include "G4VProcess.hh"
#include "G4Step.hh"
#include "G4Track.hh"

#include "G4ProcessTest.hh"
#include "G4TestSetup.hh"
#include "G4MaterialSetup.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4TestFactory.hh"

#include "G4TestUI.hh"

int main()
{
  // Setup

  G4MaterialSetup materialSetup;
  materialSetup.makeMaterials();

  G4TestUI ui;
  ui.configure();

  G4Material* material = ui.getSelectedMaterial();
  G4ParticleDefinition* definition = ui.getParticleDefinition();
  G4double eMin = ui.getMinEnergy();
  G4double eMax = ui.getMaxEnergy();
  G4TestSetup setup(material,definition,eMin,eMax);
  setup.makeGeometry();

  G4TestFactory testSelector;
  G4String type = ui.getProcessType();
  G4String category = ui.getProcessCategory();
  G4bool polarised = ui.getPolarisationSelection();
  G4ProcessTest* test = testSelector.createTestProcess(type,category,polarised);
  test->buildTables(category,polarised);
  // DoIt test

  G4int nIterations = ui.getNumberOfIterations();
  for (G4int iter=0; iter<nIterations; iter++)
    {
      G4cout << "---- Iteration " << iter << G4endl;
      G4Track* track = setup.makeTrack();
      G4Step* step = setup.makeStep();
      test->postStepTest(*track,*step);
    }

  cout << "End of the test" << G4endl;
}

