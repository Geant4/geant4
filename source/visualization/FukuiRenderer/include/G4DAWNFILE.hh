// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DAWNFILE.hh,v 1.4 1999-05-10 15:38:26 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Satoshi TANAKA
// DAWNFILE driver factory.

//=================//
#if defined (G4VIS_BUILD_DAWNFILE_DRIVER) || defined (G4VIS_USE_DAWNFILE)
//=================//

#ifndef G4DAWNFILE_HH
#define G4DAWNFILE_HH

#include "G4VGraphicsSystem.hh"

	//----- prototype
class G4VSceneHandler   ;

	//----------------------------//
	//----- class G4DAWNFILE -----// 
	//----------------------------//
class G4DAWNFILE: public G4VGraphicsSystem {

public:
  G4DAWNFILE ();
  virtual ~G4DAWNFILE ();
  G4VSceneHandler* CreateSceneHandler (const G4String& name = "");
  G4VViewer*  CreateViewer  (G4VSceneHandler&, const G4String& name = "");

private:


};

#endif
#endif 
