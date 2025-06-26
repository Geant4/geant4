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
// G4VSIntegration
//
// Class description:
//
// Numerical algorithm for integration of probability density function
// and sampling of the energy. Parameters of the algorithm should
// be defined by consumer class via the InitialiseIntegrator(..) method.
// The method is effective for the case of functions with a peak and long
// tail. Tunning of parameters may increase efficiency of the algorithm.
// The default set of parameters is optimized for the pre-compound model.
//
// Created 03.03.2025 V.Ivanchenko
//
// --------------------------------------------------------------------

#ifndef G4VSIntegration_HH
#define G4VSIntegration_HH 1

#include "globals.hh"

class G4VSIntegration
{
 public:
  G4VSIntegration() = default;

  virtual ~G4VSIntegration() = default;

  // this method should be implemented in a consumer class
  virtual G4double ProbabilityDensityFunction(G4double) = 0;

  virtual const G4String& ModelName() const;

  // initialisation before run called once
  // accuracy - accuracy of integration
  // fact1 - value < 1 to define energy E1: Pmax*fact = P(E1)
  //         and E2: P(E1)*fact = P(E2)
  // fact2 - value > 1 provides tolerance for max cross section
  // deltaE - the default step in energy
  // dmin, dmax - min and max values of step
  void InitialiseIntegrator(G4double accuracy, G4double fact1, G4double fact2,
			    G4double de, G4double dmin, G4double dmax);

  // compute integral of probability density function
  G4double ComputeIntegral(const G4double emin, const G4double  emax);

  // sample value according to probability density function
  // it is assumed that ComputeIntegral(emin, emax) was executed 
  G4double SampleValue();

  G4VSIntegration(const G4VSIntegration&) = delete;
  G4VSIntegration& operator=(const G4VSIntegration&) = delete;
  G4bool operator==(const G4VSIntegration &right) const = delete;
  G4bool operator!=(const G4VSIntegration &right) const = delete;

  void SetVerbose(G4int verb) { fVerbose = verb; }

private:

  G4double fAcc{0.001};    // accuracy of integration
  G4double fMinDelta{0.1}; // minimal step integration
  G4double fMaxDelta{2.0}; // maximal step integration
  G4double fDelta{1.0};    // the default step
  G4double fFactor1{0.25};
  G4double fFactor2{1.05};

  // parameters describing function
  G4double fEmin{0.0};
  G4double fEmax{0.0};
  G4double fE1{0.0};
  G4double fP1{0.0};
  G4double fE2{0.0}; 
  G4double fP2{0.0};
  G4double fPmax{0.0};

  G4int fVerbose{0};
  G4int fWarnLimit{4};
  G4int fnWarn{0};

  G4String dummy{""};
};

#endif
