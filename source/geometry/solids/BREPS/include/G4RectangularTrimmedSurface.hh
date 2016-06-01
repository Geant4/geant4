// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4RectangularTrimmedSurface.hh,v 2.2 1998/10/20 16:31:32 broglia Exp $
// GEANT4 tag $Name: geant4-00 $
//
#include "G4FCylindricalSurface.hh"

class G4RectangularTrimmedSurface: public G4Surface
{
public:
  G4RectangularTrimmedSurface();
  ~G4RectangularTrimmedSurface();  

  int Intersect(const G4Ray&);  
  void CalcBBox();

  virtual char *Name() const { return "G4RectangularTrimmedSurface"; }    


private:

  G4Surface* BasisSurface;

  G4double TrimU1,TrimU2;
  G4double TrimV1,TrimV2; 
 
  G4Point3D TrimPointU1, TrimPointU2;
  G4Point3D TrimPointV1, TrimPointV2;  
  
};





