// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenInventor.hh,v 1.1 1999-01-07 16:15:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  G4VScene* CreateScene(const G4String& name = "");
  G4VView* CreateView(G4VScene&,const G4String& name = "");
  void SetInteractorManager(G4VInteractorManager*);
  G4VInteractorManager* GetInteractorManager();
  void InitHEPVis();
private:
  G4VInteractorManager* interactorManager;
};

#endif

#endif
