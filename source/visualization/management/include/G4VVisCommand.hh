// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VVisCommand.hh,v 1.8 2001-02-04 20:26:19 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// Base class for visualization commands - John Allison  9th August 1998
// It is really a messenger - we have one command per messenger.

#ifndef G4VVISCOMMAND_HH
#define G4VVISCOMMAND_HH

#include "G4UImessenger.hh"
#include "G4ThreeVector.hh"
#include "g4std/vector"

class G4VisManager;
class G4UIcommand;
class G4UIcmdWithAString;

class G4VVisCommand: public G4UImessenger {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VVisCommand ();
  virtual ~G4VVisCommand ();
  static void SetVisManager (G4VisManager*);
  static G4double G4VVisCommand::ValueOf(G4String unitName);
  static G4String ConvertToString(G4bool blValue);
  static G4String ConvertToString(G4double x, G4double y,
				  const char * unitName);
  static G4String ConvertToString(const G4ThreeVector& vec);
  static G4bool        GetNewBoolValue(const G4String& paramString);
  static G4double      GetNewDoubleValue(G4String paramString);
  static G4ThreeVector GetNew3VectorValue(G4String paramString);
  static void          GetNewDoublePairValue(G4String paramString,
					     G4double& xval,
					     G4double& yval);
protected:
  static G4VisManager* fpVisManager;
  static  G4std::vector<G4UIcommand*> sceneNameCommands;
  typedef G4std::vector<G4UIcommand*>::iterator sceneNameCommandsIterator; 
  static  G4std::vector<G4UIcommand*> sceneHandlerNameCommands;
  typedef G4std::vector<G4UIcommand*>::iterator
    sceneHandlerNameCommandsIterator;
  static  G4std::vector<G4UIcommand*> viewerNameCommands;
  typedef G4std::vector<G4UIcommand*>::iterator viewerNameCommandsIterator; 
};

#include "G4VVisCommand.icc"

#endif
