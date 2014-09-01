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
// $Id:$
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4UOrb
//
// Class description:
//
//   Wrapper class for UOrb to make use of UOrb from USolids module.

// History:
// 30.10.13 G.Cosmo, CERN/PH
// --------------------------------------------------------------------
#ifndef G4ORB_HH
#define G4ORB_HH

#include "G4USolid.hh"
#include "UOrb.hh"
#include "G4Polyhedron.hh"

class G4UOrb : public G4USolid
{
  public:  // with description

    G4UOrb(const G4String& pName, G4double pRmax);
       
   ~G4UOrb() ;
    
    void ComputeDimensions(      G4VPVParameterisation* p,
                           const G4int n,
                           const G4VPhysicalVolume* pRep);

    G4VSolid* Clone() const;

    inline UOrb* GetShape() const;

    inline G4double GetRadius() const;
    inline void SetRadius(G4double newRmax);

  public:  // without description

    G4UOrb(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4UOrb(const G4UOrb& rhs);
    G4UOrb& operator=(const G4UOrb& rhs); 
      // Copy constructor and assignment operator.
};

// --------------------------------------------------------------------
// Inline methods
// --------------------------------------------------------------------

inline UOrb* G4UOrb::GetShape() const
{
  return (UOrb*) fShape;
}

inline G4double G4UOrb::GetRadius() const
{
  return GetShape()->GetRadius();
}

inline void G4UOrb::SetRadius(G4double newRmax)
{
  GetShape()->SetRadius(newRmax);
  fRebuildPolyhedron = true;
}

#endif
