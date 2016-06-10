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
// $Id: G4Win32.hh 66892 2013-01-17 10:57:59Z gunter $
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
  static G4bool    DispatchWin32Event    (void*);
  virtual         ~G4Win32               ();                     
private:
  G4Win32();
  static G4Win32* instance; // Pointer to single instance.
};

#endif

#endif
