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
// $Id: G4VRML2.hh,v 1.7 2001-07-11 10:09:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML2.hh
// Satoshi Tanaka and Yasuhide Sawada

#if defined (G4VIS_BUILD_VRML_DRIVER) || defined (G4VIS_USE_VRML)

#ifndef G4VRML2_HH
#define G4VRML2_HH

#include "G4VGraphicsSystem.hh"
#include "FRClient.h"

class G4VSceneHandler;

#include "G4VRMLNetConfig.hh"
	//#define FR_VRML_DEFAULT_PORT	40801
	//#define FR_VRML_PORT_ENV	"FR_VRML_PORT"
	//#define FR_VRML_HOST_NAME_ENV	"FR_VRML_HOST_NAME"

class G4VRML2: public G4VGraphicsSystem {
public:
	G4VRML2(); 
	virtual ~G4VRML2();
	G4VSceneHandler* CreateSceneHandler(const G4String& name = "");
	G4VViewer*  CreateViewer(G4VSceneHandler&, const G4String& name = "");

public:
	const G4String& getHostName() { return fHostName; }
	G4int getPort() { return fPort; }

private:
	G4String	fHostName;
	G4int		fPort;

};

#endif //G4VRML2_HH
#endif //G4VIS_BUILD_VRML_DRIVER
