// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ParametrisedBox.hh,v 1.2 1999-12-15 14:55:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Original class written by Hans-Peter Wellisch.

#ifndef PARAMETRISEDBOX_HH
#define PARAMETRISEDBOX_HH

#include "G4VPVParameterisation.hh"

#include "globals.hh"

class G4VPhysicalVolume;
class G4Box;

class ParametrisedBox: public G4VPVParameterisation
{
public:
  void ComputeTransformation(const G4int n,
			     G4VPhysicalVolume* pRep) const;
  void ComputeDimensions(G4Box& box,
			 const G4int n,
			 const G4VPhysicalVolume* pRep) const;
  void ComputeDimensions(G4Tubs& shape,
			 const G4int n,
			 const G4VPhysicalVolume* pV) const {
    G4VPVParameterisation::ComputeDimensions(shape,n,pV);
  }
  void ComputeDimensions(G4Trd& shape,
			 const G4int n,
			 const G4VPhysicalVolume* pV) const {
    G4VPVParameterisation::ComputeDimensions(shape,n,pV);
  }
  void ComputeDimensions(G4Trap& shape,
			 const G4int n,
			 const G4VPhysicalVolume* pV) const {
    G4VPVParameterisation::ComputeDimensions(shape,n,pV);
  }
  void ComputeDimensions(G4Cons& shape,
			 const G4int n,
			 const G4VPhysicalVolume* pV) const {
    G4VPVParameterisation::ComputeDimensions(shape,n,pV);
  }
  void ComputeDimensions(G4Sphere& shape,
				 const G4int n,
				 const G4VPhysicalVolume* pV) const {
    G4VPVParameterisation::ComputeDimensions(shape,n,pV);
  }
  void ComputeDimensions(G4Torus& shape,
			 const G4int n,
			 const G4VPhysicalVolume* pV) const {
    G4VPVParameterisation::ComputeDimensions(shape,n,pV);
  }
  void ComputeDimensions(G4Para& shape,
			 const G4int n,
			 const G4VPhysicalVolume* pV) const {
    G4VPVParameterisation::ComputeDimensions(shape,n,pV);
  }
  void ComputeDimensions(G4Hype& shape,
			 const G4int n,
			 const G4VPhysicalVolume* pV) const {
    G4VPVParameterisation::ComputeDimensions(shape,n,pV);
  }	
};

#endif
