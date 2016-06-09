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
// $Id: G4MscModel.hh,v 1.11 2004/04/29 18:40:52 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4MscModel
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

//
// Class Description:
//
// Implementation of the model of multiple scattering based on
// H.W.Lewis Phys Rev 78 (1950) 526 and L.Urban model

// -------------------------------------------------------------------
//

#ifndef G4MscModel_h
#define G4MscModel_h 1

#include "G4VEmModel.hh"

class G4MscModel : public G4VEmModel
{

public:

  G4MscModel(G4double&, G4double&, G4double&, G4double&, G4bool&,
               const G4String& nam = "MscUni");

  ~G4MscModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double HighEnergyLimit(const G4ParticleDefinition*) {return highKinEnergy;};

  G4double LowEnergyLimit(const G4ParticleDefinition*) {return lowKinEnergy;};

  void SetHighEnergyLimit(G4double e) {highKinEnergy = e;};

  void SetLowEnergyLimit(G4double e) {lowKinEnergy = e;};

  G4double MinEnergyCut(const G4ParticleDefinition*,
                        const G4MaterialCutsCouple*) {return 0.0;};

  G4bool IsInCharge(const G4ParticleDefinition*);

  G4double ComputeDEDX(const G4MaterialCutsCouple*,
                       const G4ParticleDefinition*,
                             G4double,
                             G4double) {return 0.0;};

  G4double CrossSection(const G4MaterialCutsCouple*,
                        const G4ParticleDefinition*,
                              G4double kineticEnergy,
                              G4double cutEnergy,
                              G4double maxEnergy);
  G4DynamicParticle* SampleSecondary(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double,
                                      G4double) {return 0;};

  std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double,
                                      G4double) {return 0;};

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

  G4double MaxSecondaryEnergy(const G4DynamicParticle*) {return 0.0;};

  void SetDynamicParticle(const G4DynamicParticle*);

protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
                                    G4double) {return 0.0;};
private:

  G4double ComputeTransportCrossSection(
                             const G4ParticleDefinition* particle,
                                   G4double KineticEnergy,
                                   G4double AtomicNumber,
                                   G4double AtomicWeight);

  // hide assignment operator
  G4MscModel & operator=(const  G4MscModel &right);
  G4MscModel(const  G4MscModel&);

  const G4ParticleDefinition* particle;
  G4double mass;
  G4double charge;
  G4double massRate;
  G4double highKinEnergy;
  G4double lowKinEnergy;

  G4double taubig;
  G4double tausmall;
  G4double taulim;
  G4double currentTau;
  G4double dtrl;
  G4double NuclCorrPar;
  G4double FactPar;
  G4double facxsi;

  G4double sigmafactor;
  G4double b;
  G4double xsi;

  G4double lambda0;
  G4double tPathLength;
  G4double par1,par2,par3 ;

  G4bool   samplez;

  G4double stepmin ;

  G4double currentKinEnergy;
  G4double currentRange ; 
  G4double currentRadLength;

};

inline void G4MscModel::SetDynamicParticle(const G4DynamicParticle* dp)
{
  particle = dp->GetDefinition();
  mass = dp->GetMass();
  charge = dp->GetCharge()/eplus;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

