//
//    ************************************
//    *                                  *
//    *    ThyroidROGeometry.hh          *
//    *                                  *
//    ************************************


#ifndef ThyroidROGeometry_h
#define ThyroidROGeometry_h 1

#include "G4VReadOutGeometry.hh"
#include "G4EllipticalTube.hh"
#include "G4RotationMatrix.hh"

class ThyroidROGeometry : public G4VReadOutGeometry
{
  public:
	ThyroidROGeometry(G4String aString,G4double DetrDx, G4double DetrDy, G4double DetrDz);
  	~ThyroidROGeometry();

  public:
	const G4double m_DetrDx;
	const G4double m_DetrDy;
	const G4double m_DetrDz;
     
  private:
  	G4VPhysicalVolume* Build();
};

#endif
