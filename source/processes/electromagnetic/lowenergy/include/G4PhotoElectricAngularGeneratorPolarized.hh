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
// File name:  G4PhotoElectricAngularGeneratorPolarized
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

#ifndef G4PhotoElectricAngularGeneratorPolarized_h
#define G4PhotoElectricAngularGeneratorPolarized_h 1

#include "G4VPhotoElectricAngularDistribution.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"

class G4PhotoElectricAngularGeneratorPolarized : public G4VPhotoElectricAngularDistribution
{

public:

  G4PhotoElectricAngularGeneratorPolarized(const G4String& name);

  ~G4PhotoElectricAngularGeneratorPolarized();

  G4ThreeVector GetPhotoElectronDirection(const G4ThreeVector& direction, const G4double kineticEnergy, const G4ThreeVector& polarization, const G4int shellId) const;

  void PrintGeneratorInformation() const;

protected:
  G4ThreeVector SetPerpendicularVector(const G4ThreeVector& a) const;

private:

  // hide assignment operator 
  G4PhotoElectricAngularGeneratorPolarized & operator=(const  G4PhotoElectricAngularGeneratorPolarized &right);
  G4PhotoElectricAngularGeneratorPolarized(const  G4PhotoElectricAngularGeneratorPolarized&);

private:
  void PhotoElectronGetMajorantSurfaceAandCParameters(const G4int shellLevel, const G4double beta, G4double *majorantSurfaceParameterA, G4double *majorantSurfaceParameterC) const;
  void PhotoElectronGeneratePhiAndTheta(const G4int shellLevel, const G4double beta, const G4double aBeta, 
				         const G4double cBeta, G4double *pphi, G4double *ptheta) const;
  G4ThreeVector PhotoElectronComputeFinalDirection(const G4RotationMatrix& rotation, const G4double theta, const G4double phi) const;
  G4RotationMatrix PhotoElectronRotationMatrix(const G4ThreeVector& direction, const G4ThreeVector& polarization) const;

  G4double GetMax(const G4double arg1, const G4double arg2) const;

  G4double CrossSectionMajorantFunction(const G4double theta, const G4double cBeta) const;
  G4double DSigmaKshellGavrila1959(const G4double beta, const G4double theta, const G4double phi) const;
  G4double DSigmaL1shellGavrila(const G4double beta, const G4double theta, const G4double phi) const;

  G4double betaArray[3];
  G4double aMajorantSurfaceParameterTable[980][2], cMajorantSurfaceParameterTable[980][2];
};

#endif

