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
// $Id: G4VRML1File.hh,v 1.7 2001-07-11 10:09:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4VRML1File.hh
// Satoshi Tanaka & Yasuhide Sawada

#if defined (G4VIS_BUILD_VRMLFILE_DRIVER) || defined (G4VIS_USE_VRMLFILE)

#ifndef G4VRML1FILE_HH
#define G4VRML1FILE_HH

#include "G4VGraphicsSystem.hh"

class G4VSceneHandler;

class G4VRML1File: public G4VGraphicsSystem {
public:
	G4VRML1File(); 
	virtual ~G4VRML1File();
	G4VSceneHandler* CreateSceneHandler(const G4String& name = "");
	G4VViewer*  CreateViewer(G4VSceneHandler&, const G4String& name = "");

};

#endif //G4VRML1File_HH
#endif //G4VIS_BUILD_VRMLFILE_DRIVER

