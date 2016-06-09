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
// $Id: G4MscModel71.hh,v 1.5 2007/05/22 17:34:36 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4MscMOdel71
//
// Author:        Laszlo Urban
//
// Creation date: 01.03.2001
//
// Modifications:
//
// 27-03-03 Move model part from G4MultipleScattering (V.Ivanchenko)
// 27-03-03 Rename (V.Ivanchenko)
//
// 05-08-03 angle distribution has been modified (L.Urban)
// 26-11-03 new data member currentRange (L.Urban)
// 01-03-04 changes in data members + signature changed in SampleCosineTheta
// 11-03-04 changes in data members (L.Urban)
// 23-04-04 changes in data members and in signature of SampleCosineTheta
//          (L.Urban)
// 17-08-04 name of data member facxsi changed to factail (L.Urban)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 15-04-05 optimize internal interface - add SampleSecondaries method (V.Ivanchenko)
// 03-10-05 Model is freezed with the name McsModel71 (V.Ivanchenko)
// 17-02-06 Save table of transport cross sections not mfp (V.Ivanchenko)
// 07-03-06 Create G4UrbanMscModel and move there step limit calculation (V.Ivanchenko)
//

//
// Class Description:
//
// Implementation of the model of multiple scattering based on
// H.W.Lewis Phys Rev 78 (1950) 526 and L.Urban model

// -------------------------------------------------------------------
//

#ifndef G4MscModel71_h
#define G4MscModel71_h 1

#include "G4VEmModel.hh"

class G4ParticleChangeForMSC;
class G4Navigator;

class G4MscModel71 : public G4VEmModel
{

public:

  G4MscModel71(G4double&, G4double&, G4double&, G4double&, G4bool&,
	     const G4String& nam = "MscUni");

  virtual ~G4MscModel71();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double ComputeCrossSectionPerAtom(
                             const G4ParticleDefinition* particle,
                                   G4double KineticEnergy,
                                   G4double AtomicNumber,
                                   G4double AtomicWeight=0., 
				   G4double cut =0.,
				   G4double emax=DBL_MAX);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double length,
				 G4double safety);

  G4double GeomPathLength(G4PhysicsTable* theLambdaTable,
                    const G4MaterialCutsCouple* couple,
                    const G4ParticleDefinition* particle,
                          G4double& kineticEnergy,
                          G4double lambda,
                          G4double range,
                          G4double truePathLength);

  G4double TrueStepLength(G4double geomStepLength);

  G4double SampleCosineTheta(G4double trueStepLength,G4double KineticEnergy);

  G4double SampleDisplacement();

  void SetLateralDisplasmentFlag(G4bool val);

private:

  // hide assignment operator
  G4MscModel71 & operator=(const  G4MscModel71 &right);
  G4MscModel71(const  G4MscModel71&);

  const G4ParticleDefinition* particle;
  G4ParticleChangeForMSC*     fParticleChange;
  G4Navigator*                navigator;

  G4double mass;
  G4double charge;
  G4double massRate;

  G4double taubig;
  G4double tausmall;
  G4double taulim;
  G4double currentTau;
  G4double dtrl;
  G4double NuclCorrPar;
  G4double FactPar;
  G4double factail ;

  G4double sigmafactor;
  G4double b;
  G4double xsi;

  G4double lambda0;
  G4double tPathLength;
  G4double par1,par2,par3 ;

  G4double stepmin ;

  G4double currentKinEnergy;
  G4double currentRange ; 
  G4double currentRadLength;

  G4bool   samplez;
  G4bool   latDisplasment;
  G4bool   isInitialized;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4MscModel71::SetLateralDisplasmentFlag(G4bool val)
{
  latDisplasment = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

