#ifndef __CONIC_H
#define __CONIC_H 

#include "G4Curve.hh"
#include "G4Axis2Placement3D.hh"

class G4Conic: public G4Curve
{
public:
 
  G4Conic();
  ~G4Conic();

  // Get/Set to geometric data
  const G4Axis2Placement3D* GetPosition() const;

  // pShift must be added/subtracted from the parameter 
  // no STEP I/O if not 0!!!
  // set by Project members
  G4double GetPShift() const;
  void SetPShift(G4double pShift0); 

  //inline G4Placement GetPosition() {return Position;}
  //virtual const char *Name(){return "G4ConicalCurve";}

protected:
  //void ProjectCurve(const G4Plane&, const G4Plane&);
  //int HitPartOfCurve(G4double, G4double, const G4Point2d&);
  //G4Placement Position;

  // geometric data
  G4Axis2Placement3D position;

private:

  G4double pShift;

};

#include "G4Conic.icc"

#endif
