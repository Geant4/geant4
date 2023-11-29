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
// ===========================================================================
// GEANT4 class source file
//
// Class:                G4IonParametrisedLossModel
//
// Base class:           G4VEmModel (utils)
//
// Author:               Anton Lechner (Anton.Lechner@cern.ch)
//
// First implementation: 10. 11. 2008
//
// Modifications: 03. 02. 2009 - Bug fix iterators (AL)
//                11. 03. 2009 - Introduced new table handler(G4IonDEDXHandler)
//                               and modified method to add/remove tables
//                               (tables are now built in init. phase),
//                               Minor bug fix in ComputeDEDXPerVolume (AL)
//                11. 05. 2009 - Introduced scaling algorithm for heavier ions:
//                               G4IonDEDXScalingICRU73 (AL)
//                12. 11. 2009 - Moved from original ICRU 73 classes to new
//                               class (G4IonStoppingData), which is capable
//                               of reading stopping power data files stored
//                               in G4LEDATA (requires G4EMLOW6.8 or higher).
//                               Simultanesouly, the upper energy limit of
//                               ICRU 73 is increased to 1 GeV/nucleon.
//                             - Removed nuclear stopping from Corrections-
//                               AlongStep since dedicated process was created.
//                             - Added function for switching off scaling
//                               of heavy ions from ICRU 73 data
//                             - Minor fix in ComputeLossForStep function
//                             - Minor fix in ComputeDEDXPerVolume (AL)
//                23. 11. 2009 - Changed energy loss limit from 0.15 to 0.01
//                               to improve accuracy for large steps (AL)
//                24. 11. 2009 - Bug fix: Range calculation corrected if same
//                               materials appears with different cuts in diff.
//                               regions (added UpdateRangeCache function and
//                               modified BuildRangeVector, ComputeLossForStep
//                               functions accordingly, added new cache param.)
//                             - Removed GetRange function (AL)
//                04. 11. 2010 - Moved virtual methods to the source (VI)
//
//
// Class description:
//    Model for computing the energy loss of ions by employing a
//    parameterisation of dE/dx tables (by default ICRU 73 tables). For
//    ion-material combinations and/or projectile energies not covered
//    by this model, the G4BraggIonModel and G4BetheBloch models are
//    employed.
//
// Comments:
//
// ===========================================================================
#include "G4IonParametrisedLossModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4IonStoppingData.hh"
#include "G4VIonDEDXTable.hh"
#include "G4VIonDEDXScalingAlgorithm.hh"
#include "G4IonDEDXScalingICRU73.hh"
#include "G4BraggIonModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4LossTableManager.hh"
#include "G4EmParameters.hh"
#include "G4GenericIon.hh"
#include "G4Electron.hh"
#include "G4DeltaAngle.hh"
#include "Randomize.hh"
#include "G4Exp.hh"

//#define PRINT_TABLE_BUILT
// #########################################################################

G4IonParametrisedLossModel::G4IonParametrisedLossModel(
             const G4ParticleDefinition*,
             const G4String& nam)
  : G4VEmModel(nam),
    braggIonModel(0),
    betheBlochModel(0),
    nmbBins(90),
    nmbSubBins(100),
    particleChangeLoss(0),
    corrFactor(1.0),
    energyLossLimit(0.01),
    cutEnergies(0),
    isInitialised(false)
{
  genericIon = G4GenericIon::Definition();
  genericIonPDGMass = genericIon->GetPDGMass();
  corrections = G4LossTableManager::Instance()->EmCorrections();

  // The Bragg ion and Bethe Bloch models are instantiated
  braggIonModel = new G4BraggIonModel();
  betheBlochModel = new G4BetheBlochModel();

  // The boundaries for the range tables are set
  lowerEnergyEdgeIntegr = 0.025 * MeV;
  upperEnergyEdgeIntegr = betheBlochModel -> HighEnergyLimit();

  // Cache parameters are set
  cacheParticle = 0;
  cacheMass = 0;
  cacheElecMassRatio = 0;
  cacheChargeSquare = 0;

  // Cache parameters are set
  rangeCacheParticle = 0;
  rangeCacheMatCutsCouple = 0;
  rangeCacheEnergyRange = 0;
  rangeCacheRangeEnergy = 0;

  // Cache parameters are set
  dedxCacheParticle = 0;
  dedxCacheMaterial = 0;
  dedxCacheEnergyCut = 0;
  dedxCacheIter = lossTableList.end();
  dedxCacheTransitionEnergy = 0.0;
  dedxCacheTransitionFactor = 0.0;
  dedxCacheGenIonMassRatio = 0.0;

  // default generator
  SetAngularDistribution(new G4DeltaAngle());
}

// #########################################################################

G4IonParametrisedLossModel::~G4IonParametrisedLossModel() {

  // dE/dx table objects are deleted and the container is cleared
  LossTableList::iterator iterTables = lossTableList.begin();
  LossTableList::iterator iterTables_end = lossTableList.end();

  for(;iterTables != iterTables_end; ++iterTables) { delete *iterTables; }
  lossTableList.clear();

  // range table
  RangeEnergyTable::iterator itr = r.begin();
  RangeEnergyTable::iterator itr_end = r.end();
  for(;itr != itr_end; ++itr) { delete itr->second; }
  r.clear();

  // inverse range
  EnergyRangeTable::iterator ite = E.begin();
  EnergyRangeTable::iterator ite_end = E.end();
  for(;ite != ite_end; ++ite) { delete ite->second; }
  E.clear();

}

// #########################################################################

G4double G4IonParametrisedLossModel::MinEnergyCut(
                                       const G4ParticleDefinition*,
                                       const G4MaterialCutsCouple* couple) {

  return couple->GetMaterial()->GetIonisation()->GetMeanExcitationEnergy();
}

// #########################################################################

G4double G4IonParametrisedLossModel::MaxSecondaryEnergy(
                             const G4ParticleDefinition* particle,
                             G4double kineticEnergy) {

  // ############## Maximum energy of secondaries ##########################
  // Function computes maximum energy of secondary electrons which are
  // released by an ion
  //
  // See Geant4 physics reference manual (version 9.1), section 9.1.1
  //
  // Ref.: W.M. Yao et al, Jour. of Phys. G 33 (2006) 1.
  //       C.Caso et al. (Part. Data Group), Europ. Phys. Jour. C 3 1 (1998).
  //       B. Rossi, High energy particles, New York, NY: Prentice-Hall (1952).
  //
  // (Implementation adapted from G4BraggIonModel)

  if(particle != cacheParticle) UpdateCache(particle);

  G4double tau  = kineticEnergy/cacheMass;
  G4double tmax = 2.0 * electron_mass_c2 * tau * (tau + 2.) /
                  (1. + 2.0 * (tau + 1.) * cacheElecMassRatio +
                  cacheElecMassRatio * cacheElecMassRatio);

  return tmax;
}

// #########################################################################

G4double G4IonParametrisedLossModel::GetChargeSquareRatio(
                             const G4ParticleDefinition* particle,
	 	 	     const G4Material* material,
                             G4double kineticEnergy) {    // Kinetic energy

  G4double chargeSquareRatio = corrections ->
                                     EffectiveChargeSquareRatio(particle,
                                                                material,
                                                                kineticEnergy);
  corrFactor = chargeSquareRatio *
                       corrections -> EffectiveChargeCorrection(particle,
                                                                material,
                                                                kineticEnergy);
  return corrFactor;
}

// #########################################################################

G4double G4IonParametrisedLossModel::GetParticleCharge(
                             const G4ParticleDefinition* particle,
			     const G4Material* material,
                             G4double kineticEnergy) {   // Kinetic energy

  return corrections -> GetParticleCharge(particle, material, kineticEnergy);
}

// #########################################################################

void G4IonParametrisedLossModel::Initialise(
                                       const G4ParticleDefinition* particle,
                                       const G4DataVector& cuts) {

  // Cached parameters are reset
  cacheParticle = 0;
  cacheMass = 0;
  cacheElecMassRatio = 0;
  cacheChargeSquare = 0;

  // Cached parameters are reset
  rangeCacheParticle = 0;
  rangeCacheMatCutsCouple = 0;
  rangeCacheEnergyRange = 0;
  rangeCacheRangeEnergy = 0;

  // Cached parameters are reset
  dedxCacheParticle = 0;
  dedxCacheMaterial = 0;
  dedxCacheEnergyCut = 0;
  dedxCacheIter = lossTableList.end();
  dedxCacheTransitionEnergy = 0.0;
  dedxCacheTransitionFactor = 0.0;
  dedxCacheGenIonMassRatio = 0.0;

  // By default ICRU 73 stopping power tables are loaded:
  if(!isInitialised) {
    G4bool icru90 = G4EmParameters::Instance()->UseICRU90Data();
    isInitialised = true;
    AddDEDXTable("ICRU73",
		 new G4IonStoppingData("ion_stopping_data/icru",icru90),
		 new G4IonDEDXScalingICRU73());
  }
  // The cache of loss tables is cleared
  LossTableList::iterator iterTables = lossTableList.begin();
  LossTableList::iterator iterTables_end = lossTableList.end();

  for(;iterTables != iterTables_end; ++iterTables) {
    (*iterTables) -> ClearCache();
  }

  // Range vs energy and energy vs range vectors from previous runs are
  // cleared
  RangeEnergyTable::iterator iterRange = r.begin();
  RangeEnergyTable::iterator iterRange_end = r.end();

  for(;iterRange != iterRange_end; ++iterRange) {
    delete iterRange->second;
  }
  r.clear();

  EnergyRangeTable::iterator iterEnergy = E.begin();
  EnergyRangeTable::iterator iterEnergy_end = E.end();

  for(;iterEnergy != iterEnergy_end; ++iterEnergy) {
    delete iterEnergy->second;
  }
  E.clear();

  // The cut energies 
  cutEnergies = cuts;

  // All dE/dx vectors are built
  const G4ProductionCutsTable* coupleTable=
                     G4ProductionCutsTable::GetProductionCutsTable();
  G4int nmbCouples = (G4int)coupleTable->GetTableSize();

#ifdef PRINT_TABLE_BUILT
    G4cout << "G4IonParametrisedLossModel::Initialise():"
           << " Building dE/dx vectors:"
           << G4endl;
#endif

  for (G4int i = 0; i < nmbCouples; ++i) {

    const G4MaterialCutsCouple* couple = coupleTable->GetMaterialCutsCouple(i);
    const G4Material* material = couple->GetMaterial();

    for(G4int atomicNumberIon = 3; atomicNumberIon < 102; ++atomicNumberIon) {

       LossTableList::iterator iter = lossTableList.begin();
       LossTableList::iterator iter_end = lossTableList.end();

       for(;iter != iter_end; ++iter) {

          if(*iter == 0) {
              G4cout << "G4IonParametrisedLossModel::Initialise():"
                     << " Skipping illegal table."
                     << G4endl;
          }
	  
	  if((*iter)->BuildDEDXTable(atomicNumberIon, material)) {

#ifdef PRINT_TABLE_BUILT
             G4cout << "  Atomic Number Ion = " << atomicNumberIon
                    << ", Material = " << material -> GetName()
                    << ", Table = " << (*iter) -> GetName()
                    << G4endl;
#endif
             break;
	  }
       }
    }
  }

  // The particle change object
  if(! particleChangeLoss) {
    particleChangeLoss = GetParticleChangeForLoss();
    braggIonModel->SetParticleChange(particleChangeLoss, 0);
    betheBlochModel->SetParticleChange(particleChangeLoss, 0);
  }

  // The G4BraggIonModel and G4BetheBlochModel instances are initialised with
  // the same settings as the current model:
  braggIonModel -> Initialise(particle, cuts);
  betheBlochModel -> Initialise(particle, cuts);
}

// #########################################################################

G4double G4IonParametrisedLossModel::ComputeCrossSectionPerAtom(
                                   const G4ParticleDefinition* particle,
                                   G4double kineticEnergy,
				   G4double atomicNumber,
                                   G4double,
                                   G4double cutEnergy,
                                   G4double maxKinEnergy) {

  // ############## Cross section per atom  ################################
  // Function computes ionization cross section per atom
  //
  // See Geant4 physics reference manual (version 9.1), section 9.1.3
  //
  // Ref.: W.M. Yao et al, Jour. of Phys. G 33 (2006) 1.
  //       B. Rossi, High energy particles, New York, NY: Prentice-Hall (1952).
  //
  // (Implementation adapted from G4BraggIonModel)

  G4double crosssection     = 0.0;
  G4double tmax      = MaxSecondaryEnergy(particle, kineticEnergy);
  G4double maxEnergy = std::min(tmax, maxKinEnergy);

  if(cutEnergy < tmax) {

     G4double energy  = kineticEnergy + cacheMass;
     G4double betaSquared  = kineticEnergy *
                                     (energy + cacheMass) / (energy * energy);

     crosssection = 1.0 / cutEnergy - 1.0 / maxEnergy -
                         betaSquared * std::log(maxEnergy / cutEnergy) / tmax;

     crosssection *= twopi_mc2_rcl2 * cacheChargeSquare / betaSquared;
  }

#ifdef PRINT_DEBUG_CS
  G4cout << "########################################################"
         << G4endl
         << "# G4IonParametrisedLossModel::ComputeCrossSectionPerAtom"
         << G4endl
         << "# particle =" << particle -> GetParticleName()
         <<  G4endl
         << "# cut(MeV) = " << cutEnergy/MeV
         << G4endl;

  G4cout << "#"
         << std::setw(13) << std::right << "E(MeV)"
         << std::setw(14) << "CS(um)"
         << std::setw(14) << "E_max_sec(MeV)"
         << G4endl
         << "# ------------------------------------------------------"
         << G4endl;

  G4cout << std::setw(14) << std::right << kineticEnergy / MeV
         << std::setw(14) << crosssection / (um * um)
         << std::setw(14) << tmax / MeV
         << G4endl;
#endif

  crosssection *= atomicNumber;

  return crosssection;
}

// #########################################################################

G4double G4IonParametrisedLossModel::CrossSectionPerVolume(
				   const G4Material* material,
                                   const G4ParticleDefinition* particle,
                                   G4double kineticEnergy,
                                   G4double cutEnergy,
                                   G4double maxEnergy) {

  G4double nbElecPerVolume = material -> GetTotNbOfElectPerVolume();
  G4double cross = ComputeCrossSectionPerAtom(particle,
                                              kineticEnergy,
                                              nbElecPerVolume, 0,
                                              cutEnergy,
                                              maxEnergy);

  return cross;
}

// #########################################################################

G4double G4IonParametrisedLossModel::ComputeDEDXPerVolume(
                                      const G4Material* material,
			              const G4ParticleDefinition* particle,
				      G4double kineticEnergy,
				      G4double cutEnergy) {

  // ############## dE/dx ##################################################
  // Function computes dE/dx values, where following rules are adopted:
  //   A. If the ion-material pair is covered by any native ion data
  //      parameterisation, then:
  //      * This parameterization is used for energies below a given energy
  //        limit,
  //      * whereas above the limit the Bethe-Bloch model is applied, in
  //        combination with an effective charge estimate and high order
  //        correction terms.
  //      A smoothing procedure is applied to dE/dx values computed with
  //      the second approach. The smoothing factor is based on the dE/dx
  //      values of both approaches at the transition energy (high order
  //      correction terms are included in the calculation of the transition
  //      factor).
  //   B. If the particle is a generic ion, the BraggIon and Bethe-Bloch
  //      models are used and a smoothing procedure is applied to values
  //      obtained with the second approach.
  //   C. If the ion-material is not covered by any ion data parameterization
  //      then:
  //      * The BraggIon model is used for energies below a given energy
  //        limit,
  //      * whereas above the limit the Bethe-Bloch model is applied, in
  //        combination with an effective charge estimate and high order
  //        correction terms.
  //      Also in this case, a smoothing procedure is applied to dE/dx values
  //      computed with the second model.

  G4double dEdx = 0.0;
  UpdateDEDXCache(particle, material, cutEnergy);

  LossTableList::iterator iter = dedxCacheIter;

  if(iter != lossTableList.end()) {

     G4double transitionEnergy = dedxCacheTransitionEnergy;

     if(transitionEnergy > kineticEnergy) {

        dEdx = (*iter) -> GetDEDX(particle, material, kineticEnergy);

        G4double dEdxDeltaRays = DeltaRayMeanEnergyTransferRate(material,
                                           particle,
                                           kineticEnergy,
                                           cutEnergy);
        dEdx -= dEdxDeltaRays;
     }
     else {
        G4double massRatio = dedxCacheGenIonMassRatio;

        G4double chargeSquare =
                       GetChargeSquareRatio(particle, material, kineticEnergy);

        G4double scaledKineticEnergy = kineticEnergy * massRatio;
        G4double scaledTransitionEnergy = transitionEnergy * massRatio;

        G4double lowEnergyLimit = betheBlochModel -> LowEnergyLimit();

        if(scaledTransitionEnergy >= lowEnergyLimit) {

           dEdx = betheBlochModel -> ComputeDEDXPerVolume(
                                      material, genericIon,
     				      scaledKineticEnergy, cutEnergy);

           dEdx *= chargeSquare;

           dEdx += corrections -> ComputeIonCorrections(particle,
                                                 material, kineticEnergy);

           G4double factor = 1.0 + dedxCacheTransitionFactor /
                                                           kineticEnergy;

	   dEdx *= factor;
	}
     }
  }
  else {
     G4double massRatio = 1.0;
     G4double chargeSquare = 1.0;

     if(particle != genericIon) {

        chargeSquare = GetChargeSquareRatio(particle, material, kineticEnergy);
        massRatio = genericIonPDGMass / particle -> GetPDGMass();
     }

     G4double scaledKineticEnergy = kineticEnergy * massRatio;

     G4double lowEnergyLimit = betheBlochModel -> LowEnergyLimit();
     if(scaledKineticEnergy < lowEnergyLimit) {
        dEdx = braggIonModel -> ComputeDEDXPerVolume(
                                      material, genericIon,
     				      scaledKineticEnergy, cutEnergy);

        dEdx *= chargeSquare;
     }
     else {
        G4double dEdxLimitParam = braggIonModel -> ComputeDEDXPerVolume(
                                      material, genericIon,
     				      lowEnergyLimit, cutEnergy);

        G4double dEdxLimitBetheBloch = betheBlochModel -> ComputeDEDXPerVolume(
                                      material, genericIon,
     				      lowEnergyLimit, cutEnergy);

        if(particle != genericIon) {
           G4double chargeSquareLowEnergyLimit =
               GetChargeSquareRatio(particle, material,
                                    lowEnergyLimit / massRatio);

           dEdxLimitParam *= chargeSquareLowEnergyLimit;
           dEdxLimitBetheBloch *= chargeSquareLowEnergyLimit;

           dEdxLimitBetheBloch +=
                    corrections -> ComputeIonCorrections(particle,
                                      material, lowEnergyLimit / massRatio);
        }

        G4double factor = (1.0 + (dEdxLimitParam/dEdxLimitBetheBloch - 1.0)
                               * lowEnergyLimit / scaledKineticEnergy);

        dEdx = betheBlochModel -> ComputeDEDXPerVolume(
                                      material, genericIon,
     				      scaledKineticEnergy, cutEnergy);

        dEdx *= chargeSquare;

        if(particle != genericIon) {
           dEdx += corrections -> ComputeIonCorrections(particle,
                                      material, kineticEnergy);
        }

        dEdx *= factor;
     }

  }

  if (dEdx < 0.0) dEdx = 0.0;

  return dEdx;
}

// #########################################################################

void G4IonParametrisedLossModel::PrintDEDXTable(
                   const G4ParticleDefinition* particle,  // Projectile (ion)
                   const G4Material* material,  // Absorber material
                   G4double lowerBoundary,      // Minimum energy per nucleon
                   G4double upperBoundary,      // Maximum energy per nucleon
                   G4int numBins,               // Number of bins
                   G4bool logScaleEnergy) {     // Logarithmic scaling of energy

  G4double atomicMassNumber = particle -> GetAtomicMass();
  G4double materialDensity = material -> GetDensity();

  G4cout << "# dE/dx table for " << particle -> GetParticleName()
         << " in material " << material -> GetName()
         << " of density " << materialDensity / g * cm3
         << " g/cm3"
         << G4endl
         << "# Projectile mass number A1 = " << atomicMassNumber
         << G4endl
         << "# ------------------------------------------------------"
         << G4endl;
  G4cout << "#"
         << std::setw(13) << std::right << "E"
         << std::setw(14) << "E/A1"
         << std::setw(14) << "dE/dx"
         << std::setw(14) << "1/rho*dE/dx"
         << G4endl;
  G4cout << "#"
         << std::setw(13) << std::right << "(MeV)"
         << std::setw(14) << "(MeV)"
         << std::setw(14) << "(MeV/cm)"
         << std::setw(14) << "(MeV*cm2/mg)"
         << G4endl
         << "# ------------------------------------------------------"
         << G4endl;

  G4double energyLowerBoundary = lowerBoundary * atomicMassNumber;
  G4double energyUpperBoundary = upperBoundary * atomicMassNumber;

  if(logScaleEnergy) {

     energyLowerBoundary = std::log(energyLowerBoundary);
     energyUpperBoundary = std::log(energyUpperBoundary);
  }

  G4double deltaEnergy = (energyUpperBoundary - energyLowerBoundary) /
                                                           G4double(nmbBins);

  for(int i = 0; i < numBins + 1; i++) {

      G4double energy = energyLowerBoundary + i * deltaEnergy;
      if(logScaleEnergy) energy = G4Exp(energy);

      G4double dedx = ComputeDEDXPerVolume(material, particle, energy, DBL_MAX);
      G4cout.precision(6);
      G4cout << std::setw(14) << std::right << energy / MeV
             << std::setw(14) << energy / atomicMassNumber / MeV
             << std::setw(14) << dedx / MeV * cm
             << std::setw(14) << dedx / materialDensity / (MeV*cm2/(0.001*g))
             << G4endl;
  }
}

// #########################################################################

void G4IonParametrisedLossModel::PrintDEDXTableHandlers(
                   const G4ParticleDefinition* particle,  // Projectile (ion)
                   const G4Material* material,  // Absorber material
                   G4double lowerBoundary,      // Minimum energy per nucleon
                   G4double upperBoundary,      // Maximum energy per nucleon
                   G4int numBins,               // Number of bins
                   G4bool logScaleEnergy) {     // Logarithmic scaling of energy

  LossTableList::iterator iter = lossTableList.begin();
  LossTableList::iterator iter_end = lossTableList.end();

  for(;iter != iter_end; iter++) {
      G4bool isApplicable =  (*iter) -> IsApplicable(particle, material);
      if(isApplicable) {
	(*iter) -> PrintDEDXTable(particle, material,
                                  lowerBoundary, upperBoundary,
                                  numBins,logScaleEnergy);
        break;
      }
  }
}

// #########################################################################

void G4IonParametrisedLossModel::SampleSecondaries(
                             std::vector<G4DynamicParticle*>* secondaries,
			     const G4MaterialCutsCouple* couple,
			     const G4DynamicParticle* particle,
			     G4double cutKinEnergySec,
			     G4double userMaxKinEnergySec) {
  // ############## Sampling of secondaries #################################
  // The probability density function (pdf) of the kinetic energy T of a
  // secondary electron may be written as:
  //    pdf(T) = f(T) * g(T)
  // where
  //    f(T) = (Tmax - Tcut) / (Tmax * Tcut) * (1 / T^2)
  //    g(T) = 1 - beta^2 * T / Tmax
  // where Tmax is the maximum kinetic energy of the secondary, Tcut
  // is the lower energy cut and beta is the kinetic energy of the
  // projectile.
  //
  // Sampling of the kinetic energy of a secondary electron:
  //  1) T0 is sampled from f(T) using the cumulated distribution function
  //     F(T) = int_Tcut^T f(T')dT'
  //  2) T is accepted or rejected by evaluating the rejection function g(T)
  //     at the sampled energy T0 against a randomly sampled value
  //
  //
  // See Geant4 physics reference manual (version 9.1), section 9.1.4
  //
  //
  // Reference pdf: W.M. Yao et al, Jour. of Phys. G 33 (2006) 1.
  //
  // (Implementation adapted from G4BraggIonModel)

  G4double rossiMaxKinEnergySec = MaxSecondaryKinEnergy(particle);
  G4double maxKinEnergySec =
                        std::min(rossiMaxKinEnergySec, userMaxKinEnergySec);

  if(cutKinEnergySec >= maxKinEnergySec) return;

  G4double kineticEnergy = particle -> GetKineticEnergy();

  G4double energy  = kineticEnergy + cacheMass;
  G4double betaSquared = kineticEnergy * (energy + cacheMass)
    / (energy * energy);

  G4double kinEnergySec;
  G4double grej;

  do {

    // Sampling kinetic energy from f(T) (using F(T)):
    G4double xi = G4UniformRand();
    kinEnergySec = cutKinEnergySec * maxKinEnergySec /
                        (maxKinEnergySec * (1.0 - xi) + cutKinEnergySec * xi);

    // Deriving the value of the rejection function at the obtained kinetic
    // energy:
    grej = 1.0 - betaSquared * kinEnergySec / rossiMaxKinEnergySec;

    if(grej > 1.0) {
        G4cout << "G4IonParametrisedLossModel::SampleSecondary Warning: "
               << "Majorant 1.0 < "
               << grej << " for e= " << kinEnergySec
               << G4endl;
    }

  } while( G4UniformRand() >= grej );

  const G4Material* mat =  couple->GetMaterial();
  G4int Z = SelectRandomAtomNumber(mat);

  const G4ParticleDefinition* electron = G4Electron::Electron();

  G4DynamicParticle* delta = new G4DynamicParticle(electron,
    GetAngularDistribution()->SampleDirection(particle, kinEnergySec,
                                              Z, mat),
                                                   kinEnergySec);

  secondaries->push_back(delta);

  // Change kinematics of primary particle
  G4ThreeVector direction = particle ->GetMomentumDirection();
  G4double totalMomentum = std::sqrt(kineticEnergy*(energy + cacheMass));

  G4ThreeVector finalP = totalMomentum*direction - delta->GetMomentum();
  finalP               = finalP.unit();

  kineticEnergy       -= kinEnergySec;

  particleChangeLoss->SetProposedKineticEnergy(kineticEnergy);
  particleChangeLoss->SetProposedMomentumDirection(finalP);
}

// #########################################################################

void G4IonParametrisedLossModel::UpdateRangeCache(
                     const G4ParticleDefinition* particle,
                     const G4MaterialCutsCouple* matCutsCouple) {

  // ############## Caching ##################################################
  // If the ion-material-cut combination is covered by any native ion data
  // parameterisation (for low energies), range vectors are computed

  if(particle == rangeCacheParticle &&
     matCutsCouple == rangeCacheMatCutsCouple) {
  }
  else{
     rangeCacheParticle = particle;
     rangeCacheMatCutsCouple = matCutsCouple;

     const G4Material* material = matCutsCouple -> GetMaterial();
     LossTableList::iterator iter = IsApplicable(particle, material);

     // If any table is applicable, the transition factor is computed:
     if(iter != lossTableList.end()) {

        // Build range-energy and energy-range vectors if they don't exist
        IonMatCouple ionMatCouple = std::make_pair(particle, matCutsCouple);
        RangeEnergyTable::iterator iterRange = r.find(ionMatCouple);

        if(iterRange == r.end()) BuildRangeVector(particle, matCutsCouple);

        rangeCacheEnergyRange = E[ionMatCouple];
        rangeCacheRangeEnergy = r[ionMatCouple];
     }
     else {
        rangeCacheEnergyRange = 0;
        rangeCacheRangeEnergy = 0;
     }
  }
}

// #########################################################################

void G4IonParametrisedLossModel::UpdateDEDXCache(
                     const G4ParticleDefinition* particle,
                     const G4Material* material,
                     G4double cutEnergy) {

  // ############## Caching ##################################################
  // If the ion-material combination is covered by any native ion data
  // parameterisation (for low energies), a transition factor is computed
  // which is applied to Bethe-Bloch results at higher energies to
  // guarantee a smooth transition.
  // This factor only needs to be calculated for the first step an ion
  // performs inside a certain material.

  if(particle == dedxCacheParticle &&
     material == dedxCacheMaterial &&
     cutEnergy == dedxCacheEnergyCut) {
  }
  else {

     dedxCacheParticle = particle;
     dedxCacheMaterial = material;
     dedxCacheEnergyCut = cutEnergy;

     G4double massRatio = genericIonPDGMass / particle -> GetPDGMass();
     dedxCacheGenIonMassRatio = massRatio;

     LossTableList::iterator iter = IsApplicable(particle, material);
     dedxCacheIter = iter;

     // If any table is applicable, the transition factor is computed:
     if(iter != lossTableList.end()) {

        // Retrieving the transition energy from the parameterisation table
        G4double transitionEnergy =
                 (*iter) -> GetUpperEnergyEdge(particle, material);
        dedxCacheTransitionEnergy = transitionEnergy;

        // Computing dE/dx from low-energy parameterisation at
        // transition energy
        G4double dEdxParam = (*iter) -> GetDEDX(particle, material,
                                           transitionEnergy);

        G4double dEdxDeltaRays = DeltaRayMeanEnergyTransferRate(material,
                                           particle,
                                           transitionEnergy,
                                           cutEnergy);
        dEdxParam -= dEdxDeltaRays;

        // Computing dE/dx from Bethe-Bloch formula at transition
        // energy
        G4double transitionChargeSquare =
              GetChargeSquareRatio(particle, material, transitionEnergy);

        G4double scaledTransitionEnergy = transitionEnergy * massRatio;

        G4double dEdxBetheBloch =
                           betheBlochModel -> ComputeDEDXPerVolume(
                                        material, genericIon,
     		 		        scaledTransitionEnergy, cutEnergy);
        dEdxBetheBloch *= transitionChargeSquare;

        // Additionally, high order corrections are added
        dEdxBetheBloch +=
            corrections -> ComputeIonCorrections(particle,
                                                 material, transitionEnergy);

        // Computing transition factor from both dE/dx values
        dedxCacheTransitionFactor =
                  	 (dEdxParam - dEdxBetheBloch)/dEdxBetheBloch
                             * transitionEnergy;
     }
     else {

        dedxCacheParticle = particle;
        dedxCacheMaterial = material;
        dedxCacheEnergyCut = cutEnergy;

        dedxCacheGenIonMassRatio =
                             genericIonPDGMass / particle -> GetPDGMass();

        dedxCacheTransitionEnergy = 0.0;
        dedxCacheTransitionFactor = 0.0;
     }
  }
}

// #########################################################################

void G4IonParametrisedLossModel::CorrectionsAlongStep(
                             const G4MaterialCutsCouple* couple,
			     const G4DynamicParticle* dynamicParticle,
                             const G4double& length,
			     G4double& eloss) {

  // ############## Corrections for along step energy loss calculation ######
  // The computed energy loss (due to electronic stopping) is overwritten
  // by this function if an ion data parameterization is available for the
  // current ion-material pair.
  // No action on the energy loss (due to electronic stopping) is performed
  // if no parameterization is available. In this case the original
  // generic ion tables (in combination with the effective charge) are used
  // in the along step DoIt function.
  //
  //
  // (Implementation partly adapted from G4BraggIonModel/G4BetheBlochModel)

  const G4ParticleDefinition* particle = dynamicParticle -> GetDefinition();
  const G4Material* material = couple -> GetMaterial();

  G4double kineticEnergy = dynamicParticle -> GetKineticEnergy();

  if(kineticEnergy == eloss) { return; }

  G4double cutEnergy = DBL_MAX;
  std::size_t cutIndex = couple -> GetIndex();
  cutEnergy = cutEnergies[cutIndex];

  UpdateDEDXCache(particle, material, cutEnergy);

  LossTableList::iterator iter = dedxCacheIter;

  // If parameterization for ions is available the electronic energy loss
  // is overwritten
  if(iter != lossTableList.end()) {

     // The energy loss is calculated using the ComputeDEDXPerVolume function
     // and the step length (it is assumed that dE/dx does not change
     // considerably along the step)
     eloss =
        length * ComputeDEDXPerVolume(material, particle,
                                      kineticEnergy, cutEnergy);

#ifdef PRINT_DEBUG
  G4cout.precision(6);
  G4cout << "########################################################"
         << G4endl
         << "# G4IonParametrisedLossModel::CorrectionsAlongStep"
         << G4endl
         << "# cut(MeV) = " << cutEnergy/MeV
         << G4endl;

  G4cout << "#"
         << std::setw(13) << std::right << "E(MeV)"
         << std::setw(14) << "l(um)"
         << std::setw(14) << "l*dE/dx(MeV)"
         << std::setw(14) << "(l*dE/dx)/E"
         << G4endl
         << "# ------------------------------------------------------"
         << G4endl;

  G4cout << std::setw(14) << std::right << kineticEnergy / MeV
         << std::setw(14) << length / um
         << std::setw(14) << eloss / MeV
         << std::setw(14) << eloss / kineticEnergy * 100.0
         << G4endl;
#endif

     // If the energy loss exceeds a certain fraction of the kinetic energy
     // (the fraction is indicated by the parameter "energyLossLimit") then
     // the range tables are used to derive a more accurate value of the
     // energy loss
     if(eloss > energyLossLimit * kineticEnergy) {

        eloss = ComputeLossForStep(couple, particle,
                                   kineticEnergy,length);

#ifdef PRINT_DEBUG
  G4cout << "# Correction applied:"
         << G4endl;

  G4cout << std::setw(14) << std::right << kineticEnergy / MeV
         << std::setw(14) << length / um
         << std::setw(14) << eloss / MeV
         << std::setw(14) << eloss / kineticEnergy * 100.0
         << G4endl;
#endif

     }
  }

  // For all corrections below a kinetic energy between the Pre- and
  // Post-step energy values is used
  G4double energy = kineticEnergy - eloss * 0.5;
  if(energy < 0.0) energy = kineticEnergy * 0.5;

  G4double chargeSquareRatio = corrections ->
                                     EffectiveChargeSquareRatio(particle,
                                                                material,
                                                                energy);
  GetModelOfFluctuations() -> SetParticleAndCharge(particle,
                                                   chargeSquareRatio);

  // A correction is applied considering the change of the effective charge
  // along the step (the parameter "corrFactor" refers to the effective
  // charge at the beginning of the step). Note: the correction is not
  // applied for energy loss values deriving directly from parameterized
  // ion stopping power tables
  G4double transitionEnergy = dedxCacheTransitionEnergy;

  if(iter != lossTableList.end() && transitionEnergy < kineticEnergy) {
     chargeSquareRatio *= corrections -> EffectiveChargeCorrection(particle,
                                                                material,
                                                                energy);

     G4double chargeSquareRatioCorr = chargeSquareRatio/corrFactor;
     eloss *= chargeSquareRatioCorr;
  }
  else if (iter == lossTableList.end()) {

     chargeSquareRatio *= corrections -> EffectiveChargeCorrection(particle,
                                                                material,
                                                                energy);

     G4double chargeSquareRatioCorr = chargeSquareRatio/corrFactor;
     eloss *= chargeSquareRatioCorr;
  }

  // Ion high order corrections are applied if the current model does not
  // overwrite the energy loss (i.e. when the effective charge approach is
  // used)
  if(iter == lossTableList.end()) {

     G4double scaledKineticEnergy = kineticEnergy * dedxCacheGenIonMassRatio;
     G4double lowEnergyLimit = betheBlochModel -> LowEnergyLimit();

     // Corrections are only applied in the Bethe-Bloch energy region
     if(scaledKineticEnergy > lowEnergyLimit)
       eloss += length *
            corrections -> IonHighOrderCorrections(particle, couple, energy);
  }
}

// #########################################################################

void G4IonParametrisedLossModel::BuildRangeVector(
                     const G4ParticleDefinition* particle,
                     const G4MaterialCutsCouple* matCutsCouple) {

  G4double cutEnergy = DBL_MAX;
  std::size_t cutIndex = matCutsCouple -> GetIndex();
  cutEnergy = cutEnergies[cutIndex];

  const G4Material* material = matCutsCouple -> GetMaterial();

  G4double massRatio = genericIonPDGMass / particle -> GetPDGMass();

  G4double lowerEnergy = lowerEnergyEdgeIntegr / massRatio;
  G4double upperEnergy = upperEnergyEdgeIntegr / massRatio;

  G4double logLowerEnergyEdge = std::log(lowerEnergy);
  G4double logUpperEnergyEdge = std::log(upperEnergy);

  G4double logDeltaEnergy = (logUpperEnergyEdge - logLowerEnergyEdge) /
                                                        G4double(nmbBins);

  G4double logDeltaIntegr = logDeltaEnergy / G4double(nmbSubBins);

  G4PhysicsFreeVector* energyRangeVector =
    new G4PhysicsFreeVector(nmbBins+1,
			    lowerEnergy,
			    upperEnergy,
			    /*spline=*/true);

  G4double dedxLow = ComputeDEDXPerVolume(material,
                                          particle,
                                          lowerEnergy,
                                          cutEnergy);

  G4double range = 2.0 * lowerEnergy / dedxLow;

  energyRangeVector -> PutValues(0, lowerEnergy, range);

  G4double logEnergy = std::log(lowerEnergy);
  for(std::size_t i = 1; i < nmbBins+1; ++i) {

      G4double logEnergyIntegr = logEnergy;

      for(std::size_t j = 0; j < nmbSubBins; ++j) {

          G4double binLowerBoundary = G4Exp(logEnergyIntegr);
          logEnergyIntegr += logDeltaIntegr;

          G4double binUpperBoundary = G4Exp(logEnergyIntegr);
          G4double deltaIntegr = binUpperBoundary - binLowerBoundary;

          G4double energyIntegr = binLowerBoundary + 0.5 * deltaIntegr;

          G4double dedxValue = ComputeDEDXPerVolume(material,
                                                    particle,
                                                    energyIntegr,
						    cutEnergy);

          if(dedxValue > 0.0) range += deltaIntegr / dedxValue;

#ifdef PRINT_DEBUG_DETAILS
              G4cout << "   E = "<< energyIntegr/MeV
                     << " MeV -> dE = " << deltaIntegr/MeV
                     << " MeV -> dE/dx = " << dedxValue/MeV*mm
                     << " MeV/mm -> dE/(dE/dx) = " << deltaIntegr /
                                                               dedxValue / mm
                     << " mm -> range = " << range / mm
                     << " mm " << G4endl;
#endif
      }

      logEnergy += logDeltaEnergy;

      G4double energy = G4Exp(logEnergy);

      energyRangeVector -> PutValues(i, energy, range);

#ifdef PRINT_DEBUG_DETAILS
      G4cout << "G4IonParametrisedLossModel::BuildRangeVector() bin = "
             << i <<", E = "
             << energy / MeV << " MeV, R = "
             << range / mm << " mm"
             << G4endl;
#endif

  }
  //vector is filled: activate spline
  energyRangeVector -> FillSecondDerivatives();

  G4double lowerRangeEdge =
                    energyRangeVector -> Value(lowerEnergy);
  G4double upperRangeEdge =
                    energyRangeVector -> Value(upperEnergy);

  G4PhysicsFreeVector* rangeEnergyVector
    = new G4PhysicsFreeVector(nmbBins+1,
			      lowerRangeEdge,
			      upperRangeEdge,
			      /*spline=*/true);

  for(std::size_t i = 0; i < nmbBins+1; ++i) {
      G4double energy = energyRangeVector -> Energy(i);
      rangeEnergyVector ->
             PutValues(i, energyRangeVector -> Value(energy), energy);
  }

  rangeEnergyVector -> FillSecondDerivatives();

#ifdef PRINT_DEBUG_TABLES
  G4cout << *energyLossVector
         << *energyRangeVector
         << *rangeEnergyVector << G4endl;
#endif

  IonMatCouple ionMatCouple = std::make_pair(particle, matCutsCouple);

  E[ionMatCouple] = energyRangeVector;
  r[ionMatCouple] = rangeEnergyVector;
}

// #########################################################################

G4double G4IonParametrisedLossModel::ComputeLossForStep(
                     const G4MaterialCutsCouple* matCutsCouple,
                     const G4ParticleDefinition* particle,
                     G4double kineticEnergy,
                     G4double stepLength) {

  G4double loss = 0.0;

  UpdateRangeCache(particle, matCutsCouple);

  G4PhysicsVector* energyRange = rangeCacheEnergyRange;
  G4PhysicsVector* rangeEnergy = rangeCacheRangeEnergy;

  if(energyRange != 0 && rangeEnergy != 0) {

     G4double lowerEnEdge = energyRange -> Energy( 0 );
     G4double lowerRangeEdge = rangeEnergy -> Energy( 0 );

     // Computing range for pre-step kinetic energy:
     G4double range = energyRange -> Value(kineticEnergy);

     // Energy below vector boundary:
     if(kineticEnergy < lowerEnEdge) {

        range =  energyRange -> Value(lowerEnEdge);
        range *= std::sqrt(kineticEnergy / lowerEnEdge);
     }

#ifdef PRINT_DEBUG
     G4cout << "G4IonParametrisedLossModel::ComputeLossForStep() range = "
            << range / mm << " mm, step = " << stepLength / mm << " mm"
            << G4endl;
#endif

     // Remaining range:
     G4double remRange = range - stepLength;

     // If range is smaller than step length, the loss is set to kinetic
     // energy
     if(remRange < 0.0) loss = kineticEnergy;
     else if(remRange < lowerRangeEdge) {

        G4double ratio = remRange / lowerRangeEdge;
        loss = kineticEnergy - ratio * ratio * lowerEnEdge;
     }
     else {

        G4double energy = rangeEnergy -> Value(range - stepLength);
        loss = kineticEnergy - energy;
     }
  }

  if(loss < 0.0) loss = 0.0;

  return loss;
}

// #########################################################################

G4bool G4IonParametrisedLossModel::AddDEDXTable(
                                const G4String& nam,
                                G4VIonDEDXTable* table,
                                G4VIonDEDXScalingAlgorithm* algorithm) {

  if(table == 0) {
     G4cout << "G4IonParametrisedLossModel::AddDEDXTable() Cannot "
            << " add table: Invalid pointer."
            << G4endl;

     return false;
  }

  // Checking uniqueness of name
  LossTableList::iterator iter = lossTableList.begin();
  LossTableList::iterator iter_end = lossTableList.end();

  for(;iter != iter_end; ++iter) {
     const G4String& tableName = (*iter)->GetName();

     if(tableName == nam) {
        G4cout << "G4IonParametrisedLossModel::AddDEDXTable() Cannot "
               << " add table: Name already exists."
               << G4endl;

        return false;
     }
  }

  G4VIonDEDXScalingAlgorithm* scalingAlgorithm = algorithm;
  if(scalingAlgorithm == 0)
     scalingAlgorithm = new G4VIonDEDXScalingAlgorithm;

  G4IonDEDXHandler* handler =
                      new G4IonDEDXHandler(table, scalingAlgorithm, nam);

  lossTableList.push_front(handler);

  return true;
}

// #########################################################################

G4bool G4IonParametrisedLossModel::RemoveDEDXTable(
				 const G4String& nam) {

  LossTableList::iterator iter = lossTableList.begin();
  LossTableList::iterator iter_end = lossTableList.end();

  for(;iter != iter_end; iter++) {
     G4String tableName = (*iter) -> GetName();

     if(tableName == nam) {
        delete (*iter);

        // Remove from table list
        lossTableList.erase(iter);

        // Range vs energy and energy vs range vectors are cleared
        RangeEnergyTable::iterator iterRange = r.begin();
        RangeEnergyTable::iterator iterRange_end = r.end();

        for(;iterRange != iterRange_end; iterRange++)
                                            delete iterRange -> second;
        r.clear();

        EnergyRangeTable::iterator iterEnergy = E.begin();
        EnergyRangeTable::iterator iterEnergy_end = E.end();

        for(;iterEnergy != iterEnergy_end; iterEnergy++)
                                            delete iterEnergy -> second;
        E.clear();

        return true;
     }
  }

  return false;
}
