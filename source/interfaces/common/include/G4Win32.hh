// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Win32.hh,v 1.1 1999-01-07 16:09:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  To unify Windows message treatment between 
// G4/interfaces Windows sessions and G4/visualizations Windows drivers.
// G.Barrand

#ifndef G4WIN32_HH
#define G4WIN32_HH

#if defined(G4INTY_BUILD_WIN32) || defined(G4INTY_USE_WIN32)

#include <windows.h>

#include "G4VInteractorManager.hh"

class G4Win32 : public G4VInteractorManager {
public:
  static G4Win32*  getInstance           ();
  static G4Win32*  getInstance           (HINSTANCE,HINSTANCE,LPSTR,int);
  G4bool           Inited                ();
  void*            GetEvent              ();
  void             FlushAndWaitExecution ();
  static G4bool    dispatchWin32Event    (void*);
  void             getWinMainArguments   (HINSTANCE*,HINSTANCE*,LPSTR*,int*);
  virtual         ~G4Win32               ();                     
private:
  G4Win32 (HINSTANCE,HINSTANCE,LPSTR,int);                     
  static G4Win32* instance; // Pointer to single instance.
};

#endif //HAS_WIN32

#endif
