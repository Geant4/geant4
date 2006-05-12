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

  G4ThreeVector GetPhotoElectronDirection(G4ThreeVector direction, G4double kineticEnergy, G4ThreeVector polarization, G4int shellId);

  void PrintGeneratorInformation() const;

protected:
  G4ThreeVector SetPerpendicularVector(G4ThreeVector& a);

private:

  // hide assignment operator 
  G4PhotoElectricAngularGeneratorPolarized & operator=(const  G4PhotoElectricAngularGeneratorPolarized &right);
  G4PhotoElectricAngularGeneratorPolarized(const  G4PhotoElectricAngularGeneratorPolarized&);

  void PhotoElectronGetac(G4int level, G4double beta, G4double *pa, G4double *pc);
  void PhotoElectronGenPhiTheta(G4int level,G4double beta, G4double a_beta, 
				G4double c_beta, G4double *pphi, G4double *ptheta);
  G4ThreeVector PhotoElectronGetPlab(G4RotationMatrix rotation, G4double theta, G4double phi);
  G4RotationMatrix PhotoElectronRotationMatrix(G4ThreeVector direction, G4ThreeVector polarization);

  G4double GetMax(G4double arg1, G4double arg2);

  G4double G2Function(G4double theta, G4double c_beta);
  G4double DSigmaKshellGavrila1959(G4double beta, G4double theta, G4double phi);
  G4double DSigmaL1shellGavrila(G4double beta, G4double theta, G4double phi);

  G4double betarray[3];
  G4double a[980][2],c[980][2];
};

#endif

