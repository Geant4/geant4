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
// File name:  G4RDPhotoElectricAngularGeneratorPolarized
//
// Author:        A.C. Farinha, L. Peralta, P. Rodrigues and A. Trindade
// 
// Creation date: 10 January 2006
//
// Class Description: 
//
// Abstract class for PhotoElectricsstrahlung Angular Generator462 Generation
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------
//

#ifndef G4RDPhotoElectricAngularGeneratorPolarized_h
#define G4RDPhotoElectricAngularGeneratorPolarized_h 1

#include "G4RDVPhotoElectricAngularDistribution.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"

class G4RDPhotoElectricAngularGeneratorPolarized : public G4RDVPhotoElectricAngularDistribution
{

public:

  G4RDPhotoElectricAngularGeneratorPolarized(const G4String& name);

  ~G4RDPhotoElectricAngularGeneratorPolarized();

  G4ThreeVector GetPhotoElectronDirection(const G4ThreeVector& direction, 
					  const G4double kineticEnergy, 
					  const G4ThreeVector& polarization, 
					  const G4int shellId) const;

  void PrintGeneratorInformation() const;

protected:
  G4ThreeVector SetPerpendicularVector(const G4ThreeVector& a) const;

private:

  // hide assignment operator 
  G4RDPhotoElectricAngularGeneratorPolarized & operator=(const  G4RDPhotoElectricAngularGeneratorPolarized &right);
  G4RDPhotoElectricAngularGeneratorPolarized(const  G4RDPhotoElectricAngularGeneratorPolarized&);

  void PhotoElectronGetMajorantSurfaceAandCParameters(const G4int shellLevel, 
						      const G4double beta, 
						      G4double *majorantSurfaceParameterA, 
						      G4double *majorantSurfaceParameterC) const;

  void PhotoElectronGeneratePhiAndTheta(const G4int shellLevel, 
					const G4double beta, 
					const G4double aBeta, 
					const G4double cBeta, 
					G4double *pphi, 
					G4double *ptheta) const;

  G4ThreeVector PhotoElectronComputeFinalDirection(const G4RotationMatrix& rotation, 
						   const G4double theta, 
						   const G4double phi) const;

  G4RotationMatrix PhotoElectronRotationMatrix(const G4ThreeVector& direction, 
					       const G4ThreeVector& polarization) const;

  G4double GetMax(const G4double arg1, const G4double arg2) const;

  G4double CrossSectionMajorantFunction(const G4double theta, const G4double cBeta) const;
  G4double DSigmaKshellGavrila1959(const G4double beta, const G4double theta, const G4double phi) const;
  G4double DSigmaL1shellGavrila(const G4double beta, const G4double theta, const G4double phi) const;

  G4double betaArray[3];
  G4double aMajorantSurfaceParameterTable[980][2], cMajorantSurfaceParameterTable[980][2];
};

#endif

