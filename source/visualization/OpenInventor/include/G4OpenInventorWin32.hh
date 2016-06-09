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
// $Id: G4OpenInventorWin32.hh,v 1.6 2004/04/08 09:39:38 gbarrand Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// OpenInventor graphics system factory.

#if defined (G4VIS_BUILD_OIWIN32_DRIVER) || defined (G4VIS_USE_OIWIN32)

#ifndef G4OPENINVENTORWIN32_HH
#define G4OPENINVENTORWIN32_HH

// For backward compatibilty.

#include "G4OpenInventorWin.hh"

class G4OpenInventorWin32: public G4OpenInventorWin {
public:
  G4OpenInventorWin32 () {}
  virtual ~G4OpenInventorWin32 () {}
};

#endif

#endif
