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
//
// G4PSCellFluxForCylinder3D
#include "G4PSCellFluxForCylinder3D.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for 3D scoring cell flux.
//   The Cell Flux is defined by  a sum of track length divided
//   by the geometry volume, where all of the tracks in the geometry
//   are taken into account. e.g. the unit of Cell Flux is mm/mm3.
//
//
//   If you want to score only tracks passing through the geometry volume,
//  please use G4PSPassageCellFlux3D.
//
//
// Created: 2007-08-14  Tsukasa ASO
// 2010-07-22   Introduce Unit specification.
// 2011-03-24   Give Size and Segmentation for relicated volume in cylinder.
//
///////////////////////////////////////////////////////////////////////////////

G4PSCellFluxForCylinder3D::G4PSCellFluxForCylinder3D(G4String name, G4int ni,
                                                     G4int nj, G4int nk,
                                                     G4int depi, G4int depj,
                                                     G4int depk)
  : G4PSCellFlux3D(name, ni, nj, nk, depi, depj, depk)
{
  nSegment[0] = nSegment[1] = nSegment[2] = 0;
}

G4PSCellFluxForCylinder3D::G4PSCellFluxForCylinder3D(G4String name,
                                                     const G4String& unit,
                                                     G4int ni, G4int nj,
                                                     G4int nk, G4int depi,
                                                     G4int depj, G4int depk)
  : G4PSCellFlux3D(name, unit, ni, nj, nk, depi, depj, depk)
{
  nSegment[0] = nSegment[1] = nSegment[2] = 0;
}

void G4PSCellFluxForCylinder3D::SetCylinderSize(G4ThreeVector cylSize, G4double StartAng, G4double AngSpan)
{
  cylinderSize = cylSize;   // rMin, rMax, halfZ
  fAngle[0] = StartAng;
  fAngle[1] = AngSpan;
}
void G4PSCellFluxForCylinder3D::SetNumberOfSegments(G4int nSeg[3])
{
  nSegment[0] = nSeg[0];  // Z
  nSegment[1] = nSeg[1];  // Phi
  nSegment[2] = nSeg[2];  // R
}

G4double G4PSCellFluxForCylinder3D::ComputeVolume(G4Step*, G4int idx)
{
  G4double dr = (cylinderSize[1] - cylinderSize[0]) / nSegment[2];
  G4double r0 = cylinderSize[0] + dr * (idx);
  G4double r1 = cylinderSize[0] + dr * (idx + 1);
  G4double dRArea = (r1 * r1 - r0 * r0) * pi;

  // cylinderSize is given in Half Size
  G4double fullz    = cylinderSize[2] / nSegment[0] * 2.;
  G4double phiRatio = (fAngle[1] / (CLHEP::twopi*rad)) / nSegment[1];
  G4double v        = dRArea * fullz * phiRatio;

  if(verboseLevel > 9)
  {
    G4cout << " r0= " << r0 / cm << "  r1= " << r1 / cm
           << " fullz=" << fullz / cm << G4endl;
    G4cout << " idx= " << idx << "  v(cm3)= " << v / cm3 << G4endl;
  }

  return v;
}
