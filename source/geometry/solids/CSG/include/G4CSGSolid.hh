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
// $Id: G4CSGSolid.hh 105315 2017-07-20 14:35:13Z gcosmo $
//
//  
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4CSGSolid
//
// Class description:
//
//   An abstract class inherited from G4VSolid for Constructed Solids.
//   Used primarily to structure inheritance tree.

// History:
// 27.03.98 J.Apostolakis   Created first version.
// --------------------------------------------------------------------

#ifndef G4CSGSOLID_HH
#define G4CSGSOLID_HH

#include "G4VSolid.hh"

class G4CSGSolid : public G4VSolid
{
  public:  // with description

    G4CSGSolid(const G4String& pName);
    virtual ~G4CSGSolid();

    virtual std::ostream& StreamInfo(std::ostream& os) const;

    virtual G4Polyhedron* GetPolyhedron () const;

  public:  // without description

    G4CSGSolid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.

    G4CSGSolid(const G4CSGSolid& rhs);
    G4CSGSolid& operator=(const G4CSGSolid& rhs); 
      // Copy constructor and assignment operator.

  protected:

    G4double GetRadiusInRing(G4double rmin, G4double rmax) const;
      // Utility function for subclasses to generate proper circular areas

  protected:

    G4double fCubicVolume;
    G4double fSurfaceArea;
    mutable G4bool fRebuildPolyhedron;
    mutable G4Polyhedron* fpPolyhedron;
};

#endif
