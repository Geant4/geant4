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
// GEANT 4 class source file
//
//
// G4PolyhedraHistorical.cc
//
// Implementation of polyhedra data
//
// --------------------------------------------------------------------

#include "G4PolyhedraHistorical.hh"

G4PolyhedraHistorical::G4PolyhedraHistorical()
  : Start_angle(0.), Opening_angle(0.), numSide(0), Num_z_planes(0),
    Z_values(0), Rmin(0), Rmax(0)
{
}

G4PolyhedraHistorical::G4PolyhedraHistorical( G4int z_planes )
  : Start_angle(0.), Opening_angle(0.), numSide(0),
    Num_z_planes(z_planes)
{
  Z_values = new G4double[z_planes];
  Rmin     = new G4double[z_planes];
  Rmax     = new G4double[z_planes];
  
  for( G4int i = 0; i < z_planes; i++)
  {
    Z_values[i] = 0.0;
    Rmin[i]     = 0.0;
    Rmax[i]     = 0.0;
  }
}

G4PolyhedraHistorical::~G4PolyhedraHistorical()
{
  delete [] Z_values;
  delete [] Rmin;
  delete [] Rmax;
}

G4PolyhedraHistorical::
G4PolyhedraHistorical( const G4PolyhedraHistorical& source )
{
  Start_angle   = source.Start_angle;
  Opening_angle = source.Opening_angle;
  numSide       = source.numSide;
  Num_z_planes  = source.Num_z_planes;
  
  Z_values = new G4double[Num_z_planes];
  Rmin     = new G4double[Num_z_planes];
  Rmax     = new G4double[Num_z_planes];
  
  for( G4int i = 0; i < Num_z_planes; i++)
  {
    Z_values[i] = source.Z_values[i];
    Rmin[i]     = source.Rmin[i];
    Rmax[i]     = source.Rmax[i];
  }
}

G4PolyhedraHistorical&
G4PolyhedraHistorical::operator=( const G4PolyhedraHistorical& right )
{
  if ( &right == this ) return *this;

  Start_angle   = right.Start_angle;
  Opening_angle = right.Opening_angle;
  numSide       = right.numSide;
  Num_z_planes  = right.Num_z_planes;
  
  delete [] Z_values;
  delete [] Rmin;
  delete [] Rmax;
  Z_values = new G4double[Num_z_planes];
  Rmin     = new G4double[Num_z_planes];
  Rmax     = new G4double[Num_z_planes];
  
  for( G4int i = 0; i < Num_z_planes; i++)
  {
    Z_values[i] = right.Z_values[i];
    Rmin[i]     = right.Rmin[i];
    Rmax[i]     = right.Rmax[i];
  }

  return *this;
}
