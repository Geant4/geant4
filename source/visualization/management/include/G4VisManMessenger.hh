// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisManMessenger.hh,v 1.6 2001-02-23 15:43:19 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// GEANT4 Visualization Manager Messenger - John Allison 22nd July 1996.

#ifndef G4VISMANMESSENGER_HH
#define G4VISMANMESSENGER_HH

#include "G4UImessenger.hh"
#include "globals.hh"
#include "g4std//vector"

class G4VisManager;
class G4UIcommand;

class G4VisManMessenger: public G4UImessenger {
public:
  G4VisManMessenger (G4VisManager* pVMan);
  virtual ~G4VisManMessenger ();
  void SetNewValue (G4UIcommand* command, G4String newValues);
  G4String GetCurrentValue (G4UIcommand* command);
private:
  G4VisManMessenger (const G4VisManMessenger&);
  G4VisManMessenger& operator = (const G4VisManMessenger&);
  void AddCommandCamera      ();
  //  void AddCommandClear       ();
  //  void AddCommandCopy        ();
  //  void AddCommandCreateScene ();
  void AddCommandCreateView  ();
  void AddCommandDraw        ();
  void AddCommandLights      ();
  //  void AddCommandPrint       ();
  //  void AddCommandRefresh     ();
  void AddCommandSet         ();
  //  void AddCommandShow        ();
  void AddCommandExpert      ();
  void DoCommandCamera      (const G4String& commandPath, G4String& newValues);
  //  void DoCommandClear       (const G4String& commandPath, G4String& newValues);
  //  void DoCommandCopy        (const G4String& commandPath, G4String& newValues);
  //  void DoCommandCreateScene (const G4String& commandPath, G4String& newValues);
  void DoCommandCreateView  (const G4String& commandPath, G4String& newValues);
  void DoCommandDraw        (const G4String& commandPath, G4String& newValues);
  void DoCommandLights      (const G4String& commandPath, G4String& newValues);
  //  void DoCommandPrint       (const G4String& commandPath, G4String& newValues);
  //  void DoCommandRefresh     (const G4String& commandPath, G4String& newValues);
  void DoCommandSet         (const G4String& commandPath, G4String& newValues);
  //  void DoCommandShow        (const G4String& commandPath, G4String& newValues);
  void DoCommandExpert      (const G4String& commandPath, G4String& newValues);
  G4bool ViewValid ();
  void RotateViewpointAboutUpVectorBy (G4double dbeta);
  G4VisManager* fpVMan;
  G4std::vector<G4UIcommand*> fCommandList;
};

#endif
