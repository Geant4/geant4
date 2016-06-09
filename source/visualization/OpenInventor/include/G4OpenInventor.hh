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
// $Id: G4OpenInventor.hh,v 1.7 2004/04/08 09:39:37 gbarrand Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Guy Barrand 26 Mar 1998.
// OpenInventor graphics system factory.

#ifndef G4OPENINVENTOR_HH
#define G4OPENINVENTOR_HH

#if defined(G4VIS_BUILD_OI_DRIVER) || defined(G4VIS_USE_OI)

#include "G4VGraphicsSystem.hh"

class G4VInteractorManager;

class G4OpenInventor: public G4VGraphicsSystem {
public:
  G4OpenInventor(const G4String,const G4String,G4VGraphicsSystem::Functionality);
  virtual ~G4OpenInventor();
  void SetInteractorManager(G4VInteractorManager*);
  G4VInteractorManager* GetInteractorManager();
  G4VSceneHandler* CreateSceneHandler (const G4String& name);
  void InitNodes();
private:
  G4VInteractorManager* interactorManager;
};

#endif

#endif
