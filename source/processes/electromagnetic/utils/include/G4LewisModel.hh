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
// GEANT4 Class header file
//
//
// File name:     G4LewisModel
//
// Author:        Laszlo Urban
//
// Creation date: 01.03.2001
//
// Modifications:
//
// 27-03-03 Move model part from G4MultipleScattering (V.Ivanchenko)
//
//
// Class Description:
//
// Implementation of Lewis model of multiple scattering
// H.W.Lewis Phys Rev 78 (1950) 526

// -------------------------------------------------------------------
//

#ifndef G4LewisModel_h
#define G4LewisModel_h 1

#include "G4VEmModel.hh"

class G4LewisModel : public G4VEmModel
{

public:

  G4LewisModel(G4double&, G4double&, G4double&, G4double&, G4bool&,
               const G4String& nam = "Lewis");

  ~G4LewisModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double HighEnergyLimit(const G4ParticleDefinition*) {return highKinEnergy;};

  G4double LowEnergyLimit(const G4ParticleDefinition*) {return lowKinEnergy;};

  void SetHighEnergyLimit(G4double e) {highKinEnergy = e;};

  void SetLowEnergyLimit(G4double e) {lowKinEnergy = e;};

  G4double MinEnergyCut(const G4ParticleDefinition*,
                        const G4MaterialCutsCouple*) {return 0.0;};

  G4bool IsInCharge(const G4ParticleDefinition*);

  G4double ComputeDEDX(const G4Material*,
                       const G4ParticleDefinition*,
                             G4double kineticEnergy,
                             G4double cutEnergy) {return 0.0;};

  G4double CrossSection(const G4Material*,
                        const G4ParticleDefinition*,
                              G4double kineticEnergy,
                              G4double cutEnergy,
                              G4double maxEnergy);

  G4DynamicParticle* SampleSecondary(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double maxEnergy) {return 0;};

  G4std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double maxEnergy) {return 0;};

  G4double GeomPathLength(G4PhysicsTable* theLambdaTable,
                    const G4MaterialCutsCouple* couple,
		    const G4ParticleDefinition* particle,
		          G4double& kineticEnergy,
			  G4double lambda,
			  G4double range,
			  G4double truePathLength);

  G4double TrueStepLength(G4double geomStepLength);

  G4double SampleCosineTheta(G4double trueStepLength);

  G4double SampleDisplacement();

  G4double MaxSecondaryEnergy(const G4DynamicParticle*) {return 0.0;};

protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
                                    G4double kinEnergy) {return 0.0;};

private:

  G4double ComputeTransportCrossSection(
                             const G4ParticleDefinition* particle,
                                   G4double KineticEnergy,
                                   G4double AtomicNumber,
                                   G4double AtomicWeight);

  // hide assignment operator
  G4LewisModel & operator=(const  G4LewisModel &right);
  G4LewisModel(const  G4LewisModel&);

  const G4ParticleDefinition* particle;
  G4double mass;
  G4double chargeSquare;
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
  G4double alfa1;
  G4double alfa2;
  G4double alfa3;
  G4double b;
  G4double xsi;
  G4double c0;

  G4double lambda0;
  G4double lambda1;
  G4double lambdam;
  G4double alam;
  G4double zm;
  G4double cthm;
  G4double tPathLength;

  G4bool   samplez;

  G4double stepmin ;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
