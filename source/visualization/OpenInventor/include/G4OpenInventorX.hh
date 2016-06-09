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
// $Id: G4OpenInventorX.hh,v 1.6 2004/04/08 09:39:38 gbarrand Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 
// Andrew Walkden  27th March 1996
// OpenInventor graphics system factory.

#if defined (G4VIS_BUILD_OIX_DRIVER) || defined (G4VIS_USE_OIX)

#ifndef G4OPENINVENTORX_HH
#define G4OPENINVENTORX_HH

// For backward compatibilty.

#include "G4OpenInventorXt.hh"

class G4OpenInventorX: public G4OpenInventorXt {
public:
  G4OpenInventorX () {}
  virtual ~G4OpenInventorX () {}
};

#endif

#endif
