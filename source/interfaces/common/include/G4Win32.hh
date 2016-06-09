//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4Win32.hh,v 1.7 2004/04/08 08:16:34 gbarrand Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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
  G4bool           Inited                ();
  void*            GetEvent              ();
  void             FlushAndWaitExecution ();
  static G4bool    dispatchWin32Event    (void*);
  virtual         ~G4Win32               ();                     
private:
  G4Win32();
  static G4Win32* instance; // Pointer to single instance.
};

#endif

#endif
