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
// $Id: processTest.cc,v 1.3 2001-10-15 17:36:13 pia Exp $
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



int main()
{
  // Setup

  G4int nIterations;
  G4cout << "How many interactions? " << G4endl;
  G4cin >> nIterations;
  if (nIterations <= 0) G4Exception("Wrong input");

  G4TestSetup setup;
  setup.init();

  // Process to be tested
  G4VProcess* process = setup.createTestProcess();

  G4ProcessTest test;

  // DoIt test
  for (G4int iter=0; iter<nIterations; iter++)
    {
      G4cout << "---- Iteration " << iter << G4endl;
      G4Track* track = setup.makeTrack();
      G4Step* step = setup.makeStep();
      test.postStepTest(process,*track,*step);
    }
  
  cout << "End of the test" << G4endl;
}

