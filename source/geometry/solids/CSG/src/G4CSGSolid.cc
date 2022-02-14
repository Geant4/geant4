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
// Implementation of G4CSGSolid
//
// 27.03.98 J.Apostolakis - First version.
// --------------------------------------------------------------------

#include <cmath>

#include "G4CSGSolid.hh"
#include "G4QuickRand.hh"
#include "G4Polyhedron.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

//////////////////////////////////////////////////////////////////////////
//
// Constructor
//  - Base class constructor

G4CSGSolid::G4CSGSolid(const G4String& name) :
  G4VSolid(name)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4CSGSolid::G4CSGSolid( __void__& a )
  : G4VSolid(a)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor
//

G4CSGSolid::~G4CSGSolid()
{
  delete fpPolyhedron; fpPolyhedron = nullptr;
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//

G4CSGSolid::G4CSGSolid(const G4CSGSolid& rhs)
  : G4VSolid(rhs), fCubicVolume(rhs.fCubicVolume),
    fSurfaceArea(rhs.fSurfaceArea)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4CSGSolid& G4CSGSolid::operator = (const G4CSGSolid& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VSolid::operator=(rhs);

   // Copy data
   //
   fCubicVolume = rhs.fCubicVolume;
   fSurfaceArea = rhs.fSurfaceArea;
   fRebuildPolyhedron = false;
   delete fpPolyhedron; fpPolyhedron = nullptr;

   return *this;
}

G4double G4CSGSolid::GetRadiusInRing(G4double rmin, G4double rmax) const
{
  G4double k = G4QuickRand();
  return (rmin <= 0) ? rmax*std::sqrt(k)
                     : std::sqrt(k*rmax*rmax + (1. - k)*rmin*rmin);
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
  if (fpPolyhedron == nullptr ||
      fRebuildPolyhedron ||
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
    {
      G4AutoLock l(&polyhedronMutex);
      delete fpPolyhedron;
      fpPolyhedron = CreatePolyhedron();
      fRebuildPolyhedron = false;
      l.unlock();
    }
  return fpPolyhedron;
}
