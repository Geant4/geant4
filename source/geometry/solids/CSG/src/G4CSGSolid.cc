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
// $Id: G4CSGSolid.cc,v 1.11 2005/11/09 15:03:09 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// --------------------------------------------------------------------

#include "G4CSGSolid.hh"

#include "G4Polyhedron.hh"

//////////////////////////////////////////////////////////////////////////
//
// Constructor
//  - Base class constructor 

G4CSGSolid::G4CSGSolid(const G4String& name) :
  G4VSolid(name), fCubicVolume(0.), fpPolyhedron(0)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4CSGSolid::G4CSGSolid( __void__& a )
  : G4VSolid(a), fCubicVolume(0.), fpPolyhedron(0)
{
}

G4CSGSolid::~G4CSGSolid() 
{
  delete fpPolyhedron;
}

std::ostream& G4CSGSolid::StreamInfo(std::ostream& os) const
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

G4Polyhedron* G4CSGSolid::GetPolyhedron () const
{
  if (!fpPolyhedron ||
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
    {
      delete fpPolyhedron;
      fpPolyhedron = CreatePolyhedron();
    }
  return fpPolyhedron;
}
