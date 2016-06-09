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
// $Id: G4NURBStube.hh,v 1.7 2003/04/03 15:31:06 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// 
// Olivier Crumeyrolle  12 September 1996

// Tube builder prototype
// OC 090796

// Class Description:
// Tube builder prototype for NURBS.
// See documentation in graphics_reps/doc for details.

#ifndef __C_G4NURBStube__
#define __C_G4NURBStube__ 1 

#include "G4NURBS.hh"

class  G4NURBStube : public G4NURBS
{
  public:
    G4NURBStube(G4double RMIN, G4double RMAX, G4double DZ);
    virtual G4Visible&  operator = (const G4Visible& right);
    virtual G4VVisPrim& operator = (const G4VVisPrim& right);
    const char*  Whoami() const;
};
#endif
// end of __C_G4NURBStube__
