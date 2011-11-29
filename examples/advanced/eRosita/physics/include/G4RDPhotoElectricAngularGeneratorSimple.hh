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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:  G4RDPhotoElectricAngularGeneratorSimple
//
// Author:        Andreia Trindade (andreia@lip.pt)
// 
// Creation date: 10 May 2004
//
// Modifications: 
// 10 May 2003       A. Trindade    First implementation acording with new design
//
// Class Description: 
//
// Abstract class for PhotoElectricsstrahlung Angular GeneratorSimple Generation
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4RDPhotoElectricAngularGeneratorSimple_h
#define G4RDPhotoElectricAngularGeneratorSimple_h 1

#include "G4RDVPhotoElectricAngularDistribution.hh"
#include "G4ios.hh"
#include "globals.hh"

class G4RDPhotoElectricAngularGeneratorSimple : public G4RDVPhotoElectricAngularDistribution
{

public:

  G4RDPhotoElectricAngularGeneratorSimple(const G4String& name);

  ~G4RDPhotoElectricAngularGeneratorSimple();
  G4ThreeVector GetPhotoElectronDirection(const G4ThreeVector& direction, 
					  const G4double kineticEnergy, 
					  const G4ThreeVector& polarization, 
					  const G4int shellId) const;

  void PrintGeneratorInformation() const;

protected:

private:

  // hide assignment operator 
     G4RDPhotoElectricAngularGeneratorSimple & operator=(const  G4RDPhotoElectricAngularGeneratorSimple &right);
     G4RDPhotoElectricAngularGeneratorSimple(const  G4RDPhotoElectricAngularGeneratorSimple&);

};

#endif

