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
// $Id: G4CSGSolid.hh,v 1.5 2001-07-11 09:59:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//  
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4CSGSolid
//
// Class description:
//
//   An abstract class for Constructed Solids.
//   Used primarily to simplify inheritance tree.
//
//   Member functions:
//
//   As inherited from G4VSolid 

// History:
// 27.03.98 J.Apostolakis   Created first version.
// --------------------------------------------------------------------

#ifndef G4CSGSOLID_HH
#define G4CSGSOLID_HH

#include "G4VSolid.hh"

class G4CSGSolid : public G4VSolid
{
  public:

    G4CSGSolid(const G4String& pName);

    virtual ~G4CSGSolid();
};

#endif
