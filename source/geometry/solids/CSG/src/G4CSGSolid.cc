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
// $Id: G4CSGSolid.cc,v 1.4 2002-10-28 11:43:05 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4CSGSolid.hh"

// Constructor
//  - Base class constructor 

G4CSGSolid::G4CSGSolid(const G4String& name) :
   G4VSolid(name)
{
}

G4CSGSolid::~G4CSGSolid() 
{
}

G4std::ostream& G4CSGSolid::StreamInfo(G4std::ostream& os) const
{
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: " << GetEntityType() << "\n"
     << " Parameters: \n"
     << "   NOT available !\n"
     << "-----------------------------------------------------------\n";

  return os;
}
