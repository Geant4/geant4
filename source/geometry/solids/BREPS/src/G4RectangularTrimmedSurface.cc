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
// $Id: G4RectangularTrimmedSurface.cc,v 1.8 2010-07-07 14:45:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4RectangularTrimmedSurface.cc
//
// ----------------------------------------------------------------------

#include "G4RectangularTrimmedSurface.hh"
#include "G4FPlane.hh"
#include "G4BSplineSurface.hh"
#include "G4ToroidalSurface.hh"
#include "G4SphericalSurface.hh"

G4RectangularTrimmedSurface::G4RectangularTrimmedSurface()
  : BasisSurface(0), TrimU1(0.), TrimU2(0.), TrimV1(0.), TrimV2(0.)
{
}

G4RectangularTrimmedSurface::~G4RectangularTrimmedSurface()
{
  if (BasisSurface) delete BasisSurface;
}  


const char* G4RectangularTrimmedSurface::Name() const
{
  return "G4RectangularTrimmedSurface";
}

void G4RectangularTrimmedSurface::CalcBBox()
{
  BasisSurface->CalcBBox();
  bbox = BasisSurface->GetBBox();
}


G4int G4RectangularTrimmedSurface::Intersect(const G4Ray& Rayref)
{
  if(BasisSurface->Intersect(Rayref))
  {
    G4double UHit = BasisSurface->GetUHit();
    G4double VHit = BasisSurface->GetVHit();
    
    if((TrimU1<=UHit)&&(TrimU2>=UHit)&&(TrimV1<=VHit)&&(TrimV2>=VHit))
    {
      closest_hit = BasisSurface->GetClosestHit();
      return 1;
    }
  }
  
  return 0;
}
