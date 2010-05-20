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
// $Id: ParametrisedBox.hh,v 1.5 2010-05-20 18:10:09 allison Exp $
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
  void ComputeDimensions(G4Orb& shape,
			 const G4int n,
			 const G4VPhysicalVolume* pV) const {
    G4VPVParameterisation::ComputeDimensions(shape,n,pV);
  }	
  void ComputeDimensions(G4Polycone& shape,
			 const G4int n,
			 const G4VPhysicalVolume* pV) const {
    G4VPVParameterisation::ComputeDimensions(shape,n,pV);
  }	
  void ComputeDimensions(G4Polyhedra& shape,
			 const G4int n,
			 const G4VPhysicalVolume* pV) const {
    G4VPVParameterisation::ComputeDimensions(shape,n,pV);
  }	
};

#endif
