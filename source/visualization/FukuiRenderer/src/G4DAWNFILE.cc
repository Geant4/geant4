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
// $Id: G4DAWNFILE.cc 78838 2014-01-28 08:46:17Z gcosmo $
//
// Satoshi TANAKA
// DAWNFILE factory.


//#define DEBUG_FR_SYSTEM


#include "G4DAWNFILE.hh"

#define __G_ANSI_C__

//#include "G4VisFeaturesOfDAWNFILE.hh"
#include "G4FRFeatures.hh" 
#include "G4VSceneHandler.hh"
#include "G4DAWNFILESceneHandler.hh"
#include "G4DAWNFILEViewer.hh"
#include "G4FRConst.hh"


	//----- G4DAWNFILE, constructor
G4DAWNFILE::G4DAWNFILE ():
  G4VGraphicsSystem ("DAWNFILE",
		     "DAWNFILE",
		     FR_DAWNFILE_FEATURES,
		     G4VGraphicsSystem::fileWriter)
{}

	//----- G4DAWNFILE, destructor
G4DAWNFILE::~G4DAWNFILE () 
{}


	//-----  G4DAWNFILE::CreateSceneHandler (const G4String& name) 
G4VSceneHandler* G4DAWNFILE::CreateSceneHandler (const G4String& name) 
{
	G4VSceneHandler* p = new G4DAWNFILESceneHandler (*this, name);
	return p;
}

	//-----  G4DAWNFILE::CreateViewer ()
G4VViewer* G4DAWNFILE::CreateViewer (G4VSceneHandler& scene, const G4String& name) 
{
       	G4VViewer* pView = 
	  new G4DAWNFILEViewer ((G4DAWNFILESceneHandler&) scene, name);
	return pView;
}
