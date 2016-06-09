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
// $Id: G4BREPSolidOpenPCone.cc,v 1.15 2010-10-20 09:14:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

G4BREPSolidOpenPCone::
G4BREPSolidOpenPCone ( const G4String& name,
                             G4double sangle,
                             G4double oangle,
                             G4int    nplanes, // sections,
                             G4double zstart,                 
                             G4double zvalues[],
                             G4double radmin[],
                             G4double radmax[]  )
 : G4IntersectionSolid ( name,
                         new G4BREPSolidPCone( name,
                                               sangle, oangle, 
                                               nplanes, zstart, zvalues,
                                               radmin, radmax ),
                         new G4Tubs( "IntersectionTubs",
                                     0., 1., 1., sangle, oangle ) )
{
  // Save local data
  //
  constructorParams.start_angle = sangle;
  constructorParams.opening_angle = oangle;
  constructorParams.num_z_planes = nplanes;
  constructorParams.z_start = zstart;
  constructorParams.z_values = 0;
  constructorParams.RMIN     = 0;
  constructorParams.RMAX     = 0;
  if ( nplanes>0 )
  {
    constructorParams.z_values = new G4double[nplanes];
    constructorParams.RMIN     = new G4double[nplanes];
    constructorParams.RMAX     = new G4double[nplanes];
    for ( G4int i = 0; i < nplanes; ++i )
    { 
      constructorParams.z_values[i] = zvalues[i];
      constructorParams.RMIN[i] = radmin[i];
      constructorParams.RMAX[i] = radmax[i];
    }
  }

  // compute max radius
  //
  InitializeOPCone();
}

G4BREPSolidOpenPCone::G4BREPSolidOpenPCone( __void__& a )
  : G4IntersectionSolid(a)
{
  constructorParams.start_angle    = 0.;
  constructorParams.opening_angle  = 0.;
  constructorParams.num_z_planes   = 0;
  constructorParams.z_start        = 0.;
  constructorParams.z_values = 0;
  constructorParams.RMIN = 0;
  constructorParams.RMAX = 0;
}

G4BREPSolidOpenPCone::~G4BREPSolidOpenPCone()
{
  if( constructorParams.num_z_planes > 0 )
  {
    delete [] constructorParams.z_values;
    delete [] constructorParams.RMIN;
    delete [] constructorParams.RMAX;
  }
}

G4BREPSolidOpenPCone::G4BREPSolidOpenPCone(const G4BREPSolidOpenPCone& rhs)
  : G4IntersectionSolid ( rhs.GetName(),
    new G4BREPSolidPCone( rhs.GetName(), rhs.constructorParams.start_angle,
                          rhs.constructorParams.opening_angle, 
                          rhs.constructorParams.num_z_planes,
                          rhs.constructorParams.z_start,
                          rhs.constructorParams.z_values,
                          rhs.constructorParams.RMIN,
                          rhs.constructorParams.RMAX ),
    new G4Tubs( "IntersectionTubs", 0., 1., 1.,
                rhs.constructorParams.start_angle,
                rhs.constructorParams.opening_angle ) )

{
  SetName(rhs.GetName());
  constructorParams.start_angle = rhs.constructorParams.start_angle;
  constructorParams.opening_angle = rhs.constructorParams.opening_angle;
  constructorParams.num_z_planes = rhs.constructorParams.num_z_planes;
  constructorParams.z_start = rhs.constructorParams.z_start;
  constructorParams.z_values = 0;
  constructorParams.RMIN     = 0;
  constructorParams.RMAX     = 0;
  G4int nplanes = constructorParams.num_z_planes;
  if( nplanes > 0 )
  {
    constructorParams.z_values = new G4double[nplanes];
    constructorParams.RMIN     = new G4double[nplanes];
    constructorParams.RMAX     = new G4double[nplanes];
    for ( G4int i = 0; i < nplanes; ++i )
    { 
      constructorParams.z_values[i] = rhs.constructorParams.z_values[i];
      constructorParams.RMIN[i] = rhs.constructorParams.RMIN[i];
      constructorParams.RMAX[i] = rhs.constructorParams.RMAX[i];
    }
  }
  InitializeOPCone();
}

G4BREPSolidOpenPCone&
G4BREPSolidOpenPCone::operator = (const G4BREPSolidOpenPCone& rhs) 
{
  // Check assignment to self
  //
  if (this == &rhs)  { return *this; }

  // Copy base class data
  //
  G4IntersectionSolid::operator=(rhs);

  // Copy data
  //
  SetName(rhs.GetName());
  constructorParams.start_angle = rhs.constructorParams.start_angle;
  constructorParams.opening_angle = rhs.constructorParams.opening_angle;
  constructorParams.num_z_planes = rhs.constructorParams.num_z_planes;
  constructorParams.z_start = rhs.constructorParams.z_start;
  G4int nplanes = constructorParams.num_z_planes;
  if( nplanes > 0 )
  {
    delete [] constructorParams.z_values;
    delete [] constructorParams.RMIN;
    delete [] constructorParams.RMAX;
    constructorParams.z_values = new G4double[nplanes];
    constructorParams.RMIN     = new G4double[nplanes];
    constructorParams.RMAX     = new G4double[nplanes];
    for( G4int idx = 0; idx < nplanes; ++idx )
    {
      constructorParams.z_values[idx] = rhs.constructorParams.z_values[idx];
      constructorParams.RMIN[idx]     = rhs.constructorParams.RMIN[idx];
      constructorParams.RMAX[idx]     = rhs.constructorParams.RMAX[idx];      
    }
  }
  InitializeOPCone();

  return *this;
}  

void G4BREPSolidOpenPCone::InitializeOPCone()
{
  G4double MaxRMAX = 0;
  for ( G4int i = 0; i < constructorParams.num_z_planes; ++i ) 
    if ( constructorParams.RMAX[i] > MaxRMAX )
      MaxRMAX = constructorParams.RMAX[i];
  		
  G4double length =
    constructorParams.z_values[constructorParams.num_z_planes-1]
  - constructorParams.z_values[0];
  
  ((G4Tubs*)fPtrSolidB)->SetOuterRadius ( MaxRMAX );
  ((G4Tubs*)fPtrSolidB)->SetZHalfLength ( length );
}

void G4BREPSolidOpenPCone::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddSolid ( *this );
}

G4VSolid* G4BREPSolidOpenPCone::Clone() const
{
  return new G4BREPSolidOpenPCone(*this);
}

std::ostream& G4BREPSolidOpenPCone::StreamInfo(std::ostream& os) const
{  
  // Streams solid contents to output stream.

  G4IntersectionSolid::StreamInfo( os );

  return os;
}

