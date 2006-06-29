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
// $Id: test33.cc,v 1.7 2006-06-29 21:59:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - test33
//
// --------------------------------------------------------------
// Comments:
// The main function may be given a macro to execute file as argument. 
// If no argument is given a session is started.
// --------------------------------------------------------------

#include "CLHEP/Random/Random.h"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4UIsession.hh"

#include "Tst33AppStarter.hh"

#include "geomdefs.hh"


// const G4double kCarTolerance = 1E-9*mm;
// const G4double kRadTolerance = 1E-9*mm;
// const G4double kAngTolerance = 1E-9*rad;


int main(int argc, char **argv)
{  

  CLHEP::HepRandom::setTheSeed(345354);

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
