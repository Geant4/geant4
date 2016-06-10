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
// $Id: G4GMocrenFile.cc 81943 2014-06-06 15:56:00Z gcosmo $
//
//
// GMocrenFile factory.
//
// Created:  Mar. 31, 2009  Akinori Kimura  
//
//

#include "G4GMocrenFile.hh"

#include "G4VSceneHandler.hh"
#include "G4GMocrenFileSceneHandler.hh"
#include "G4GMocrenFileViewer.hh"
#include "G4GMocrenMessenger.hh"


	//----- G4GMocrenFile, constructor
G4GMocrenFile::G4GMocrenFile ()
 : G4VGraphicsSystem ("gMocrenFile", "gMocrenFile",
		      "A gMocren file driver (ver.4)",
                      G4VGraphicsSystem::fileWriter),
		      //GMOCRENFILE_FEATURES, G4VGraphicsSystem::threeD),
   kViewer(NULL), kSceneHandler(NULL), kMessenger(new G4GMocrenMessenger()) {
  ;
}

	//----- G4GMocrenFile, destructor
G4GMocrenFile::~G4GMocrenFile () {
  delete kMessenger;
}


	//-----  G4GMocrenFile::CreateSceneHandler (const G4String& name) 
G4VSceneHandler* G4GMocrenFile::CreateSceneHandler (const G4String& name) 
{
  kSceneHandler = new G4GMocrenFileSceneHandler (*this, *kMessenger, name);
  return kSceneHandler;
}

	//-----  G4GMocrenFile::CreateViewer ()
G4VViewer* G4GMocrenFile::CreateViewer (G4VSceneHandler& scene,
					const G4String& name) 
{
  kViewer =  new G4GMocrenFileViewer((G4GMocrenFileSceneHandler&) scene,
				     *kMessenger, name);
  return kViewer;
}
