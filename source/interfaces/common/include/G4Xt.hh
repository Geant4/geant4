// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Xt.hh,v 1.2 1999-04-13 01:26:31 yhajime Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  To unify X11 event treatment between 
// G4/interfaces Xt sessions and G4/visualizations Xt drivers.
// G.Barrand

#ifndef G4XT_HH
#define G4XT_HH

#if defined(G4INTY_BUILD_XT) || defined(G4INTY_USE_XT)

#include <X11/Intrinsic.h>

#include "G4VInteractorManager.hh"

class G4Xt : public G4VInteractorManager {
public:
  static G4Xt* getInstance                 ();
  static G4Xt* getInstance                 (int,char**,char*);
  void         PutStringInResourceDatabase (char*);
  G4bool       Inited                      ();
  void*        GetEvent                    ();
  void         FlushAndWaitExecution       ();
  virtual     ~G4Xt                        ();                     
private:
  G4Xt (int,char**,char*);                     
  static G4Xt* instance; // Pointer to single instance.
};

#endif //HAS_XT

#endif
