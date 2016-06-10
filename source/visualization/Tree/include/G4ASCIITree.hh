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
// $Id: G4ASCIITree.hh 77479 2013-11-25 10:01:22Z gcosmo $
//
// 
// John Allison  5th April 2001
// A graphics system to dump geometry hierarchy to standard output as
//   ASCII stream.

#ifndef G4ASCIITREE_HH
#define G4ASCIITREE_HH

#include "G4VTree.hh"

class G4ASCIITreeMessenger;

class G4ASCIITree: public G4VTree {
public:
  G4ASCIITree ();
  virtual ~G4ASCIITree ();
  G4VSceneHandler* CreateSceneHandler (const G4String& name = "");
  G4VViewer*  CreateViewer  (G4VSceneHandler&, const G4String& name = "");
  G4int    GetVerbosity() const {return fVerbosity;}
  G4String GetOutFileName () const {return fOutFileName;}
  void SetVerbosity    (G4int verbosity) {fVerbosity = verbosity;}
  void SetOutFileName (const G4String& name)  {fOutFileName = name;}
protected:
  G4int fVerbosity;
  G4ASCIITreeMessenger* fpMessenger;
  G4String fOutFileName;
};

#endif
