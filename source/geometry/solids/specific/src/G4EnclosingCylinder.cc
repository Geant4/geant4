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
// $Id: G4EnclosingCylinder.cc,v 1.5 2002-10-30 13:52:23 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4EnclosingCylinder.cc
//
// Implementation of a utility class for a quick check of geometry.
//
// --------------------------------------------------------------------

#include "G4EnclosingCylinder.hh"
#include "G4ReduciblePolygon.hh"

//
// Constructor
//
G4EnclosingCylinder::G4EnclosingCylinder( const G4ReduciblePolygon *rz,
                                                G4bool thePhiIsOpen, 
                                                G4double theStartPhi,
                                                G4double theTotalPhi )
  : startPhi(theStartPhi), totalPhi(theTotalPhi),
    rx1(0.), ry1(0.), dx1(0.), dy1(0.),
    rx2(0.), ry2(0.), dx2(0.), dy2(0.),     
    concave(theTotalPhi > M_PI)
{
  //
  // Obtain largest r and smallest and largest z
  //
  radius = rz->Amax();
  zHi = rz->Bmax();
  zLo = rz->Bmin();
  
  //
  // Save phi info
  //
  phiIsOpen = thePhiIsOpen;
  if ( phiIsOpen )
  {    
    rx1 = cos(startPhi);
    ry1 = sin(startPhi);
    dx1 = +ry1*10*kCarTolerance;
    dy1 = -rx1*10*kCarTolerance;
    
    rx2 = cos(startPhi+totalPhi);
    ry2 = sin(startPhi+totalPhi);
    dx2 = -ry2*10*kCarTolerance;
    dy2 = +rx2*10*kCarTolerance;
  }
  
  //
  // Add safety
  //
  radius += 10*kCarTolerance;
  zLo    -= 10*kCarTolerance;
  zHi    += 10*kCarTolerance;
}

//
// Destructor
//
G4EnclosingCylinder::~G4EnclosingCylinder()
{
}


//
// Outside
//
// Decide very rapidly if the point is outside the cylinder
//
// If one is not certain, return false
//
G4bool G4EnclosingCylinder::MustBeOutside( const G4ThreeVector &p ) const
{
  if (p.perp() > radius) return true;
  if (p.z() < zLo) return true;
  if (p.z() > zHi) return true;

  if (phiIsOpen)
  {
    if (concave)
    {
      if ( ((p.x()-dx1)*ry1 - (p.y()-dy1)*rx1) < 0) return false;
      if ( ((p.x()-dx2)*ry2 - (p.y()-dy2)*rx2) > 0) return false;
    }
    else
    {
      if ( ((p.x()-dx1)*ry1 - (p.y()-dy1)*rx1) > 0) return true;
      if ( ((p.x()-dx2)*ry2 - (p.y()-dy2)*rx2) < 0) return true;
    }
  }
  
  return false;
}
    
    
//
// Misses
//
// Decide very rapidly if the trajectory is going to miss the cylinder
//
// If one is not sure, return false
//
G4bool G4EnclosingCylinder::ShouldMiss( const G4ThreeVector &p,
                                        const G4ThreeVector &v ) const
{
  if (!MustBeOutside(p)) return false;
  
  G4double cross = p.x()*v.y() - p.y()*v.x();
  if (cross > radius) return true;
  
  if (p.perp() > radius)
  {
    G4double dot = p.x()*v.x() + p.y()*v.y();
    if (dot > 0) return true;
  }

  return false;
}
