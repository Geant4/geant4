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

#ifndef G4VInteractiveSession_H
#define G4VInteractiveSession_H 1

#include "g4std/map"

#include "G4VInteractorManager.hh"

class G4UImessenger;

// Class description :
//
//  G4VInteractiveSession : a base class to isolate common things
// to various graphical "basic" G4UI sessions ; AddMenu, AddButton,
// a dictionnary for widgets.
//  The word "interactor" is for "piece of user interface" or
// "widget" (which means nothing).
//
// Class description - end :

class G4VInteractiveSession {
public:
  G4VInteractiveSession();
  virtual ~G4VInteractiveSession();
  virtual void AddMenu (const char*,const char*);
  virtual void AddButton (const char*,const char*,const char*);
  void AddInteractor(G4String,G4Interactor);
  G4Interactor GetInteractor(G4String);
private:
  G4UImessenger* messenger;
  typedef G4std::map<G4String,G4Interactor, G4std::less<G4String> > G4interactor_map;
  G4interactor_map interactors;
};

#endif

