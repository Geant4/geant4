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
// G4PSFlatSurfaceFlux3D
#include "G4PSFlatSurfaceFlux3D.hh"

///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring Surface Flux.
//  Current version assumes only for G4Box shape, and the surface
//  is defined at the -Z plane of the box.
//   The surface flux is given in the unit of area.
//    e.g.  sum of 1/cos(T)/mm2,  where T is a incident angle of the
//                                track on the surface.
//
//
// Surface is defined at the -Z surface.
// Direction                  -Z   +Z
//   0  IN || OUT            ->|<-  |        fFlux_InOut
//   1  IN                   ->|    |        fFlux_In
//   2  OUT                    |<-  |        fFlux_Out
//
// Created: 2007-08-14  Tsukasa ASO
// 2010-07-22   Introduce Unit specification.
//
///////////////////////////////////////////////////////////////////////////////

G4PSFlatSurfaceFlux3D::G4PSFlatSurfaceFlux3D(G4String name, G4int direction,
                                             G4int ni, G4int nj, G4int nk,
                                             G4int di, G4int dj, G4int dk)
  : G4PSFlatSurfaceFlux(name, direction)
  , fDepthi(di)
  , fDepthj(dj)
  , fDepthk(dk)
{
  SetNijk(ni, nj, nk);
}

G4PSFlatSurfaceFlux3D::G4PSFlatSurfaceFlux3D(G4String name, G4int direction,
                                             const G4String& unit, G4int ni,
                                             G4int nj, G4int nk, G4int di,
                                             G4int dj, G4int dk)
  : G4PSFlatSurfaceFlux3D(name, direction, ni, nj, nk, di, dj, dk)
{
  SetUnit(unit);
}

G4int G4PSFlatSurfaceFlux3D::GetIndex(G4Step* aStep)
{
  const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
  G4int i                       = touchable->GetReplicaNumber(fDepthi);
  G4int j                       = touchable->GetReplicaNumber(fDepthj);
  G4int k                       = touchable->GetReplicaNumber(fDepthk);

  return i * fNj * fNk + j * fNk + k;
}
