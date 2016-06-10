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
// $Id: G4VGraphicsSystem.hh 91686 2015-07-31 09:40:08Z gcosmo $
//
// 
// John Allison  27th March 1996
//
// Class description
//
// Abstract interface class for graphics systems.

#ifndef G4VGRAPHICSSYSTEM_HH
#define G4VGRAPHICSSYSTEM_HH

#include "globals.hh"

#include <vector>

class G4VSceneHandler;
class G4VViewer;

class G4VGraphicsSystem {

public: // With description

  enum Functionality {
    noFunctionality
    ,nonEuclidian       // e.g., tree representation of geometry hierarchy.
    ,twoD               // Simple 2D, e.g., X (no stored structures).
    ,twoDStore          // 2D with stored structures.
    ,threeD             // Passive 3D (with stored structures).
    ,threeDInteractive  // 3D with "pick" functionality.
    ,virtualReality     // Virtual Reality functionality.
    ,fileWriter         // File writer
  };

  G4VGraphicsSystem (const G4String& name,
		     Functionality f);

  G4VGraphicsSystem (const G4String& name,
		     const G4String& nickname,
		     Functionality f);

  G4VGraphicsSystem (const G4String& name,
		     const G4String& nickname,
		     const G4String& description,
		     Functionality f);

  virtual ~G4VGraphicsSystem ();

  virtual G4VSceneHandler* CreateSceneHandler (const G4String& name) = 0;

  virtual G4VViewer* CreateViewer (G4VSceneHandler&, const G4String& name) = 0;

  // Access functions.
  const G4String& GetName                   () const {return fName;}
  const G4String& GetNickname               () const {return fNicknames[0];}
  const std::vector<G4String>& GetNicknames () const {return fNicknames;}
  const G4String& GetDescription            () const {return fDescription;}
  Functionality   GetFunctionality          () const {return fFunctionality;}
  virtual G4bool  IsUISessionCompatible     () const;
  void AddNickname (const G4String& nickname) {fNicknames.push_back(nickname);}

protected:
  G4String fName;
  std::vector<G4String> fNicknames;
  G4String fDescription;
  Functionality  fFunctionality;
};

std::ostream& operator << (std::ostream& os, const G4VGraphicsSystem& gs);

#endif
