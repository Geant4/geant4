// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DAWNFILE.cc,v 1.1 1999-01-07 16:14:37 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Satoshi TANAKA
// DAWNFILE factory.

//=================//
#ifdef G4VIS_BUILD_DAWNFILE_DRIVER
//=================//


#include "G4DAWNFILE.hh"

#define __G_ANSI_C__

//#include "G4VisFeaturesOfDAWNFILE.hh"
#include "G4FRFeatures.hh" 
#include "G4VScene.hh"
#include "G4DAWNFILEScene.hh"
#include "G4DAWNFILEView.hh"
#include "G4FRConst.hh"


	//----- G4DAWNFILE, constructor
G4DAWNFILE::G4DAWNFILE ():
  G4VGraphicsSystem ("DAWNFILE",
		     "DAWNFILE",
		     FR_DAWNFILE_FEATURES,
		     G4VGraphicsSystem::threeD)
{}

	//----- G4DAWNFILE, destructor
G4DAWNFILE::~G4DAWNFILE () 
{}


	//-----  G4DAWNFILE::CreateScene (const G4String& name) 
G4VScene* G4DAWNFILE::CreateScene (const G4String& name) 
{
	G4VScene* p = new G4DAWNFILEScene (*this, name);

	G4cout	<< G4DAWNFILEScene::GetSceneCount ()
		<< ' ' << fName << " scenes extanct." << endl;

	return p;
}

	//-----  G4DAWNFILE::CreateView ()
G4VView* G4DAWNFILE::CreateView (G4VScene& scene, const G4String& name) 
{
       	G4VView* pView = 
	  new G4DAWNFILEView ((G4DAWNFILEScene&) scene, name);
	return pView;
}


#endif // G4VIS_BUILD_DAWNFILE_DRIVER
