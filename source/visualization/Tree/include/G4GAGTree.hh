// Satoshi Tanaka  31 May 2001
// A graphics system to dump geometry hierarchy to GAG.
//

#ifndef G4GAGTREE_HH
#define G4GAGTREE_HH

#include "G4VTree.hh"

class G4GAGTreeMessenger;

class G4GAGTree: public G4VTree {
public:
  G4GAGTree ();
  virtual ~G4GAGTree ();
  G4VSceneHandler* CreateSceneHandler (const G4String& name = "");
  G4VViewer*  CreateViewer  (G4VSceneHandler&, const G4String& name = "");
  void SetVerbosity(G4int verbosity) {fVerbosity = verbosity;}
  G4int GetVerbosity() const {return fVerbosity;}
protected:
  G4int fVerbosity;
  G4GAGTreeMessenger* fpMessenger;
};

#endif
