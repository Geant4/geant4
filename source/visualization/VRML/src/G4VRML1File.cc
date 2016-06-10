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
// $Id: G4VRML1File.cc 78838 2014-01-28 08:46:17Z gcosmo $
//
// G4VRML1File.cc
// Satoshi Tanaka & Yasuhide Sawada

#include <stdio.h> // sscanf
#include <stdlib.h> // getenv

#include "G4VSceneHandler.hh"

#include "G4VRML1File.hh"
#include "G4VRML1FileSceneHandler.hh"
#include "G4VRML1FileViewer.hh"

#include "G4FRClient.hh"


G4VRML1File::G4VRML1File() :
	G4VGraphicsSystem("VRML1FILE", "VRML1FILE", G4VGraphicsSystem::fileWriter)
{
}

G4VRML1File::~G4VRML1File()
{
}


G4VSceneHandler* G4VRML1File::CreateSceneHandler(const G4String& name) 
{
	G4VSceneHandler *p = NULL;

	p = new G4VRML1FileSceneHandler(*this, name);

	return p;
}

G4VViewer* G4VRML1File::CreateViewer(G4VSceneHandler& scene, const G4String& name)
{
	G4VViewer* pView = NULL;

	G4VRML1FileSceneHandler* pScene = (G4VRML1FileSceneHandler*)&scene;
	pView = new G4VRML1FileViewer(*pScene, name);

	return pView;
}
