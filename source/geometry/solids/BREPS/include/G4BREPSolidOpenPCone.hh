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
// $Id: G4BREPSolidOpenPCone.hh,v 1.6 2002-11-06 23:28:48 radoone Exp $
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
  
  virtual G4std::ostream& StreamInfo(G4std::ostream& os) const;
    // Streams solid contents to output stream.

 private:

  G4BREPSolidOpenPCone(const G4BREPSolidOpenPCone&);
  G4BREPSolidOpenPCone& operator=(const G4BREPSolidOpenPCone&);

};

#endif
