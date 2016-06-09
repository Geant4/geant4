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
// $Id: G4BREPSolidOpenPCone.cc,v 1.11 2006/06/29 18:41:21 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4BREPSolidOpenPCone.cc
//
// ----------------------------------------------------------------------

#include "G4BREPSolidOpenPCone.hh"
#include "G4BREPSolidPCone.hh"
#include "G4Tubs.hh"
#include "G4VGraphicsScene.hh"

G4BREPSolidOpenPCone::G4BREPSolidOpenPCone
                                  (const G4String& name,
                                   G4double start_angle,
                                   G4double opening_angle,
                                   G4int    num_z_planes, // sections,
                                   G4double z_start,                 
                                   G4double z_values[],
                                   G4double RMIN[],
                                   G4double RMAX[]
                                   )
 : G4IntersectionSolid ( name,
                         new G4BREPSolidPCone( name,
                                               start_angle, opening_angle, 
                                               num_z_planes, z_start, z_values,
                                               RMIN, RMAX
                                             ),
                         new G4Tubs( "IntersectionTubs",
                                     0., 1., 1., start_angle, opening_angle
                                   )
                       )
   //, constructorParams.z_values( 0 ), constructorParams.RMIN( 0 ), constructorParams.RMAX( 0 )
{

// compute max radius

  G4double MaxRMAX = 0;
  for ( int i = 0; i < num_z_planes; i++ ) 
    if ( RMAX[i] > MaxRMAX ) MaxRMAX = RMAX[i];
  		
  G4double length = z_values[num_z_planes-1] - z_values[0];
  
  ((G4Tubs*)fPtrSolidB)->SetOuterRadius ( MaxRMAX );
  ((G4Tubs*)fPtrSolidB)->SetZHalfLength ( length );

}

G4BREPSolidOpenPCone::G4BREPSolidOpenPCone( __void__& a )
  : G4IntersectionSolid(a)
{
}

G4BREPSolidOpenPCone::~G4BREPSolidOpenPCone()
{
}

void G4BREPSolidOpenPCone::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddSolid ( *this );
}

// Streams solid contents to output stream.
std::ostream& G4BREPSolidOpenPCone::StreamInfo(std::ostream& os) const
{  
  G4IntersectionSolid::StreamInfo( os );

  return os;
}

