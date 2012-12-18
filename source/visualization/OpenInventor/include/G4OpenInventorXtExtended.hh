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
// $Id: G4OpenInventorXt.hh,v 1.3 2006-06-29 21:20:38 gunter Exp $
//
// 
// Andrew Walkden  27th March 1996
// OpenInventor graphics system factory.
// Frederick Jones and TJR  October 2012
// Extended driver based on G4OpenInventorXt.hh
// Uses G4OpenInventorXtExaminerViewer.

#ifndef G4OPENINVENTORXTEXTENDED_HH
#define G4OPENINVENTORXTEXTENDED_HH

#include "G4OpenInventor.hh"

class G4OpenInventorXtExtended: public G4OpenInventor {
public:
  G4OpenInventorXtExtended ();
  virtual ~G4OpenInventorXtExtended ();
  G4VViewer* CreateViewer(G4VSceneHandler&,const G4String& name = "");
private:
  virtual void Initialize();
private:
  bool fInited;
};

#endif
