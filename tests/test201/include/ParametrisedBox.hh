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
// $Id: ParametrisedBox.hh,v 1.3 2001-07-11 10:10:22 gunter Exp $
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
