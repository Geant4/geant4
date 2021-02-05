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

#ifndef G4VInteractiveSession_H
#define G4VInteractiveSession_H 1

#include <map>

#include "G4VInteractorManager.hh"

class G4UImessenger;

// Class description :
//
// G4VInteractiveSession is a base class to isolate common things
// to various graphical "basic" G4UI sessions ; AddMenu, AddButton,
// a dictionnary for widgets.
// The word "interactor" is for "piece of user interface" or "widget".

class G4VInteractiveSession
{
  public:  //with description

    G4VInteractiveSession();
    virtual ~G4VInteractiveSession();
    virtual void AddMenu (const char*,const char*);
    virtual void AddButton (const char*,const char*,const char*);
    virtual void AddIcon (const char*,const char*,const char*,const char*);
    virtual void DefaultIcons (bool);
    virtual void OutputStyle (const char*,const char*,const char*);
    void AddInteractor(G4String,G4Interactor);
    G4Interactor GetInteractor(G4String);

  private:

  G4UImessenger* messenger;
  typedef std::map<G4String,G4Interactor, std::less<G4String> > G4interactor_map;
  G4interactor_map interactors;
};

#endif

