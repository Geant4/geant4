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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:  G4VPhotoElectricAngularDistribution
//
// Author:     Pedro Rodrigues  (psilva@lip.pt)
//             Andreia Trindade (andreia@lip.pt)
// 
// Creation date: 10 May 2004
//
// Modifications: 
// 10 May 2004       A. Trindade    First implementation acording with new design
//
// Class Description: 
//
// Abstract class for PhotoElectricsstrahlung Angular Distribution Generation
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4VPhotoElectricAngularDistribution_h
#define G4VPhotoElectricAngularDistribution_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4VPhotoElectricAngularDistribution 
{

public:

  G4VPhotoElectricAngularDistribution(const G4String& name);

  virtual ~G4VPhotoElectricAngularDistribution();

  virtual G4ThreeVector GetPhotoElectronDirection(G4ThreeVector direction, G4double kineticEnergy) = 0;

  virtual void PrintGeneratorInformation() const = 0;

protected:

private:

  // hide assignment operator 
     G4VPhotoElectricAngularDistribution & operator=(const  G4VPhotoElectricAngularDistribution &right);
     G4VPhotoElectricAngularDistribution(const  G4VPhotoElectricAngularDistribution&);

};

#endif

