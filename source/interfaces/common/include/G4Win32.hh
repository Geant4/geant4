// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Win32.hh,v 1.5 1999-12-15 14:50:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  To unify Windows message treatment between 
// G4/interfaces Windows sessions and G4/visualizations Windows drivers.
// G.Barrand

#ifndef G4WIN32_HH
#define G4WIN32_HH

#if defined(G4INTY_BUILD_WIN32) || defined(G4INTY_USE_WIN32)

#include <windows.h>
#include <windowsx.h>

#include "G4VInteractorManager.hh"

// Class description :
//
//  G4Win32 : a singleton to handle GUI sessions and visualization 
// drivers built over Windows. It permits to have one Windows main 
// loop for the whole application. 
//
// Class description - end :

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
