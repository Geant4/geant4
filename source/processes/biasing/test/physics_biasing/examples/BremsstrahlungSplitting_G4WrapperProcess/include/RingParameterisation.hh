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
// $Id: RingParameterisation.hh,v 1.1 2007-05-24 21:57:03 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation, based on 
//                         BeamTestScoreParameterisation by T. Aso
//
#ifndef RINGPARAMETERISATION_HH
#define RINGPARAMETERISATION_HH

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4VPVParameterisation.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Sphere;
class G4Material;
class G4VisAttributes;

// Dummy declarations to get rid of warnings ...
class G4Trd;
class G4Trap;
class G4Cons;
class G4Orb;
class G4Tubs;
class G4Torus;
class G4Para;
class G4Hype;
class G4Box;
class G4Polycone;
class G4Polyhedra;

class RingParameterisation : public G4VPVParameterisation {

public:

  RingParameterisation(const G4double rmin, const G4double rmax,
		       const G4double dtheta, const G4double initialTheta );
  
  virtual ~RingParameterisation();
  
  virtual void ComputeTransformation(const G4int, G4VPhysicalVolume*) const {};

  virtual void ComputeDimensions(G4Sphere&, 
				 const G4int, 
				 const G4VPhysicalVolume*) const;
  
private:
  void ComputeDimensions(G4Tubs&, const G4int, const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Trd&, const G4int, const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Trap&, const G4int, const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Cons&, const G4int, const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Orb&, const G4int, const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Torus&, const G4int, const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Para&, const G4int, const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Hype&, const G4int, const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Box&, const G4int, const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Polycone&, const G4int, const G4VPhysicalVolume*) const {}
  void ComputeDimensions(G4Polyhedra&, const G4int, const G4VPhysicalVolume*) const {}

private:
   G4double fRmin;
   G4double fRmax;
   G4double fDtheta;
   G4double fMinTheta;
};

#endif



