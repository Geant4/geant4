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
// $Id: test33.cc,v 1.1 2002-10-29 15:43:03 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - test33
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "CLHEP/Random/Random.h"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4UIsession.hh"

#include "Tst33AppStarter.hh"

int main(int argc, char **argv)
{  

  HepRandom::setTheSeed(345354);

  Tst33AppStarter *appstarter = new Tst33AppStarter;

  if (argc!=2) {
    G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
    session = new G4UIterminal(new G4UItcsh);      
#else
    session = new G4UIterminal();
#endif    
    session->SessionStart();
  }
  else {
    G4UImanager::GetUIpointer()->
      ApplyCommand(G4String("/control/execute ") + G4String(argv[1]));
  }
  delete appstarter;
}
