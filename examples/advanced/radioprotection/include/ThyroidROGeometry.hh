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
