// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RectangularTrimmedSurface.cc,v 1.5 2000-11-08 14:22:11 gcosmo Exp $
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
  : BasisSurface(0)
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
