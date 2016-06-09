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
// $Id: G4BREPSolidOpenPCone.hh,v 1.12 2010-10-20 09:14:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4BREPSolidOpenPCone
//
// Class description:
//
//  Definition of a generic BREP open polyconical solid:
//
//  G4BREPSolidOpenPCone ( const G4String& name,
//                               G4double  start_angle,
//                               G4double  opening_angle,                 
//                               G4int     num_z_planes,
//                               G4double  z_start,               
//                               G4double  z_values[],
//                               G4double  RMIN[],
//                               G4double  RMAX[])

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
// ----------------------------------------------------------------------
#ifndef G4BREPSolidOpenPCone_hh
#define G4BREPSolidOpenPCone_hh

#include "G4IntersectionSolid.hh"

class G4BREPSolidOpenPCone : public G4IntersectionSolid
{

 public: // with description
  
  G4BREPSolidOpenPCone( const G4String& name,
                        G4double  start_angle,
                        G4double  opening_angle,                 
                        G4int     num_z_planes, // sections
                        G4double  z_start,               
                        G4double  z_values[],
                        G4double  RMIN[],
                        G4double  RMAX[] );
    // Constructor defining the polycone through its
    // conical/cylindrical surfaces.

  ~G4BREPSolidOpenPCone();
    // Empty destructor.

  void DescribeYourselfTo (G4VGraphicsScene& scene) const;
    // Dispatch function which identifies the solid to the graphics scene.
  
  G4VSolid* Clone() const;
    // Returns a pointer of a dynamically allocated copy of the solid.

  std::ostream& StreamInfo(std::ostream& os) const;
    // Streams solid contents to output stream.

 public:  // without description

  G4BREPSolidOpenPCone(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

  G4BREPSolidOpenPCone(const G4BREPSolidOpenPCone& rhs);
  G4BREPSolidOpenPCone& operator=(const G4BREPSolidOpenPCone& rhs);

 private:

  void InitializeOPCone();

 private:

  struct G4BREPOpenPConeParams
  {
    G4double start_angle, opening_angle;
    G4int num_z_planes;       // sections
    G4double z_start;
    G4double *z_values, *RMIN, *RMAX;
  } constructorParams;
};

#endif
