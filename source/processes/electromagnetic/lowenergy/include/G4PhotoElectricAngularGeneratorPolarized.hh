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
// File name:  G4PhotoElectricAngularGeneratorPolarized
//
// Author:        A.C. Farinha, L. Peralta, P. Rodrigues and A. Trindade
// 
// Creation date: 10 January 2006
//
// Class Description: 
//

// -------------------------------------------------------------------
//

#ifndef G4PhotoElectricAngularGeneratorPolarized_h
#define G4PhotoElectricAngularGeneratorPolarized_h 1

#include "G4VEmAngularDistribution.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"

class G4PhotoElectricAngularGeneratorPolarized : public G4VEmAngularDistribution
{

public:

  G4PhotoElectricAngularGeneratorPolarized();

  ~G4PhotoElectricAngularGeneratorPolarized();

  virtual G4ThreeVector& SampleDirection(const G4DynamicParticle* dp,
                                         G4double eKinEnergy,
                                         G4int shellId,
                                         const G4Material* mat = 0);

  void PrintGeneratorInformation() const;

protected:

  G4ThreeVector PerpendicularVector(const G4ThreeVector& a) const;

private:

  // hide assignment operator 
  G4PhotoElectricAngularGeneratorPolarized & operator=(const  G4PhotoElectricAngularGeneratorPolarized &right);
  G4PhotoElectricAngularGeneratorPolarized(const  G4PhotoElectricAngularGeneratorPolarized&);

  void PhotoElectronGetMajorantSurfaceAandCParameters(G4int shellId, 
						      G4double beta, 
						      G4double *majorantSurfaceParameterA, 
						      G4double *majorantSurfaceParameterC) const;

  void PhotoElectronGeneratePhiAndTheta(G4int shellId, 
					G4double beta, 
					G4double aBeta, 
					G4double cBeta, 
					G4double *pphi, 
					G4double *ptheta) const;

  G4ThreeVector PhotoElectronComputeFinalDirection(const G4RotationMatrix& rotation, 
						   G4double theta, 
						   G4double phi) const;

  G4RotationMatrix PhotoElectronRotationMatrix(const G4ThreeVector& direction, 
					       const G4ThreeVector& polarization);

  G4double CrossSectionMajorantFunction(G4double theta, G4double cBeta) const;
  G4double DSigmaKshellGavrila1959(G4double beta, G4double theta, G4double phi) const;
  G4double DSigmaL1shellGavrila(G4double beta, G4double theta, G4double phi) const;

  G4double betaArray[3];
  G4double aMajorantSurfaceParameterTable[980][2];
  G4double cMajorantSurfaceParameterTable[980][2];
};

#endif

