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
// $Id: G4DNARuddIonisationModel.cc 104430 2017-05-31 07:43:44Z gcosmo $
// GEANT4 tag $Name:  $
//

#include "G4DNARuddIonisationModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4DNARuddAngle.hh"
#include "G4DeltaAngle.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARuddIonisationModel::G4DNARuddIonisationModel(const G4ParticleDefinition*,
                                                   const G4String& nam) :
    G4VEmModel(nam), isInitialised(false)
{
  fpWaterDensity = 0;

  slaterEffectiveCharge[0] = 0.;
  slaterEffectiveCharge[1] = 0.;
  slaterEffectiveCharge[2] = 0.;
  sCoefficient[0] = 0.;
  sCoefficient[1] = 0.;
  sCoefficient[2] = 0.;

  lowEnergyLimitForZ1 = 0 * eV;
  lowEnergyLimitForZ2 = 0 * eV;
  lowEnergyLimitOfModelForZ1 = 100 * eV;
  lowEnergyLimitOfModelForZ2 = 1 * keV;
  killBelowEnergyForZ1 = lowEnergyLimitOfModelForZ1;
  killBelowEnergyForZ2 = lowEnergyLimitOfModelForZ2;

  verboseLevel = 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if (verboseLevel > 0)
  {
    G4cout << "Rudd ionisation model is constructed " << G4endl;
  }

  // define default angular generator
  SetAngularDistribution(new G4DNARuddAngle());

  //Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);
  fAtomDeexcitation = 0;
  fParticleChangeForGamma = 0;

 // Selection of stationary mode

  statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARuddIonisationModel::~G4DNARuddIonisationModel()
{
  // Cross section

  std::map<G4String, G4DNACrossSectionDataSet*, std::less<G4String> >::iterator pos;
  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
  {
    G4DNACrossSectionDataSet* table = pos->second;
    delete table;
  }

  // The following removal is forbidden since G4VEnergyLossmodel takes care of deletion
  // Coverity however will signal this as an error
  //if (fAtomDeexcitation) {delete  fAtomDeexcitation;}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationModel::Initialise(const G4ParticleDefinition* particle,
                                          const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
  {
    G4cout << "Calling G4DNARuddIonisationModel::Initialise()" << G4endl;
  }

  // Energy limits

  G4String fileProton("dna/sigma_ionisation_p_rudd");
  G4String fileHydrogen("dna/sigma_ionisation_h_rudd");
  G4String fileAlphaPlusPlus("dna/sigma_ionisation_alphaplusplus_rudd");
  G4String fileAlphaPlus("dna/sigma_ionisation_alphaplus_rudd");
  G4String fileHelium("dna/sigma_ionisation_he_rudd");

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();
  G4ParticleDefinition* hydrogenDef = instance->GetIon("hydrogen");
  G4ParticleDefinition* alphaPlusPlusDef = instance->GetIon("alpha++");
  G4ParticleDefinition* alphaPlusDef = instance->GetIon("alpha+");
  G4ParticleDefinition* heliumDef = instance->GetIon("helium");

  G4String proton;
  G4String hydrogen;
  G4String alphaPlusPlus;
  G4String alphaPlus;
  G4String helium;

  G4double scaleFactor = 1 * m*m;

  // LIMITS AND DATA

  // ********************************************************

  proton = protonDef->GetParticleName();
  tableFile[proton] = fileProton;

  lowEnergyLimit[proton] = lowEnergyLimitForZ1;
  highEnergyLimit[proton] = 500. * keV;

  // Cross section

  G4DNACrossSectionDataSet* tableProton = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
      eV,
      scaleFactor );
  tableProton->LoadData(fileProton);
  tableData[proton] = tableProton;

  // ********************************************************

  hydrogen = hydrogenDef->GetParticleName();
  tableFile[hydrogen] = fileHydrogen;

  lowEnergyLimit[hydrogen] = lowEnergyLimitForZ1;
  highEnergyLimit[hydrogen] = 100. * MeV;

  // Cross section

  G4DNACrossSectionDataSet* tableHydrogen = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
      eV,
      scaleFactor );
  tableHydrogen->LoadData(fileHydrogen);

  tableData[hydrogen] = tableHydrogen;

  // ********************************************************

  alphaPlusPlus = alphaPlusPlusDef->GetParticleName();
  tableFile[alphaPlusPlus] = fileAlphaPlusPlus;

  lowEnergyLimit[alphaPlusPlus] = lowEnergyLimitForZ2;
  highEnergyLimit[alphaPlusPlus] = 400. * MeV;

  // Cross section

  G4DNACrossSectionDataSet* tableAlphaPlusPlus = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
      eV,
      scaleFactor );
  tableAlphaPlusPlus->LoadData(fileAlphaPlusPlus);

  tableData[alphaPlusPlus] = tableAlphaPlusPlus;

  // ********************************************************

  alphaPlus = alphaPlusDef->GetParticleName();
  tableFile[alphaPlus] = fileAlphaPlus;

  lowEnergyLimit[alphaPlus] = lowEnergyLimitForZ2;
  highEnergyLimit[alphaPlus] = 400. * MeV;

  // Cross section

  G4DNACrossSectionDataSet* tableAlphaPlus = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
      eV,
      scaleFactor );
  tableAlphaPlus->LoadData(fileAlphaPlus);
  tableData[alphaPlus] = tableAlphaPlus;

  // ********************************************************

  helium = heliumDef->GetParticleName();
  tableFile[helium] = fileHelium;

  lowEnergyLimit[helium] = lowEnergyLimitForZ2;
  highEnergyLimit[helium] = 400. * MeV;

  // Cross section

  G4DNACrossSectionDataSet* tableHelium = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,
      eV,
      scaleFactor );
  tableHelium->LoadData(fileHelium);
  tableData[helium] = tableHelium;

  //

  if (particle==protonDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[proton]);
    SetHighEnergyLimit(highEnergyLimit[proton]);
  }

  if (particle==hydrogenDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[hydrogen]);
    SetHighEnergyLimit(highEnergyLimit[hydrogen]);
  }

  if (particle==heliumDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[helium]);
    SetHighEnergyLimit(highEnergyLimit[helium]);
  }

  if (particle==alphaPlusDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[alphaPlus]);
    SetHighEnergyLimit(highEnergyLimit[alphaPlus]);
  }

  if (particle==alphaPlusPlusDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[alphaPlusPlus]);
    SetHighEnergyLimit(highEnergyLimit[alphaPlusPlus]);
  }

  if( verboseLevel>0 )
  {
    G4cout << "Rudd ionisation model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / eV << " eV - "
    << HighEnergyLimit() / keV << " keV for "
    << particle->GetParticleName()
    << G4endl;
  }

  // Initialize water density pointer
  fpWaterDensity = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

  //

  fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();

  if (isInitialised)
  { return;}
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNARuddIonisationModel::CrossSectionPerVolume(const G4Material* material,
                                                         const G4ParticleDefinition* particleDefinition,
                                                         G4double k,
                                                         G4double,
                                                         G4double)
{
  if (verboseLevel > 3)
  {
    G4cout << "Calling CrossSectionPerVolume() of G4DNARuddIonisationModel"
        << G4endl;
  }

  // Calculate total cross section for model

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  if (
      particleDefinition != G4Proton::ProtonDefinition()
      &&
      particleDefinition != instance->GetIon("hydrogen")
      &&
      particleDefinition != instance->GetIon("alpha++")
      &&
      particleDefinition != instance->GetIon("alpha+")
      &&
      particleDefinition != instance->GetIon("helium")
  )

  return 0;

  G4double lowLim = 0;

  if ( particleDefinition == G4Proton::ProtonDefinition()
      || particleDefinition == instance->GetIon("hydrogen")
  )

  lowLim = lowEnergyLimitOfModelForZ1;

  if ( particleDefinition == instance->GetIon("alpha++")
      || particleDefinition == instance->GetIon("alpha+")
      || particleDefinition == instance->GetIon("helium")
  )

  lowLim = lowEnergyLimitOfModelForZ2;

  G4double highLim = 0;
  G4double sigma=0;

  G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];

  if(waterDensity!= 0.0)
  //  if (material == nistwater || material->GetBaseMaterial() == nistwater)
  {
    const G4String& particleName = particleDefinition->GetParticleName();

    // SI - the following is useless since lowLim is already defined
    /*
     std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
     pos1 = lowEnergyLimit.find(particleName);

     if (pos1 != lowEnergyLimit.end())
     {
     lowLim = pos1->second;
     }
     */

    std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
    pos2 = highEnergyLimit.find(particleName);

    if (pos2 != highEnergyLimit.end())
    {
      highLim = pos2->second;
    }

    if (k <= highLim)
    {
      //SI : XS must not be zero otherwise sampling of secondaries method ignored

      if (k < lowLim) k = lowLim;

      //

      std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
      pos = tableData.find(particleName);

      if (pos != tableData.end())
      {
        G4DNACrossSectionDataSet* table = pos->second;
        if (table != 0)
        {
          sigma = table->FindValue(k);
        }
      }
      else
      {
        G4Exception("G4DNARuddIonisationModel::CrossSectionPerVolume","em0002",
            FatalException,"Model not applicable to particle type.");
      }

    } // if (k >= lowLim && k < highLim)

    if (verboseLevel > 2)
    {
      G4cout << "__________________________________" << G4endl;
      G4cout << "G4DNARuddIonisationModel - XS INFO START" << G4endl;
      G4cout << "Kinetic energy(eV)=" << k/eV << " particle : " << particleDefinition->GetParticleName() << G4endl;
      G4cout << "Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
      G4cout << "Cross section per water molecule (cm^-1)=" << sigma*waterDensity/(1./cm) << G4endl;
      //      G4cout << " - Cross section per water molecule (cm^-1)=" << sigma*material->GetAtomicNumDensityVector()[1]/(1./cm) << G4endl;
      G4cout << "G4DNARuddIonisationModel - XS INFO END" << G4endl;
    }

  } // if (waterMaterial)
  else
  {
    if (verboseLevel > 2)
    {
      G4cout << "Warning : RuddIonisationModel: WATER DENSITY IS NULL" << G4endl;
    }
  }

  return sigma*waterDensity;
  //    return sigma*material->GetAtomicNumDensityVector()[1];

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARuddIonisationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                 const G4MaterialCutsCouple* couple,
                                                 const G4DynamicParticle* particle,
                                                 G4double,
                                                 G4double)
{
  if (verboseLevel > 3)
  {
    G4cout << "Calling SampleSecondaries() of G4DNARuddIonisationModel"
        << G4endl;
  }

  G4double lowLim = 0;
  G4double highLim = 0;

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  if ( particle->GetDefinition() == G4Proton::ProtonDefinition()
      || particle->GetDefinition() == instance->GetIon("hydrogen")
  )

  lowLim = killBelowEnergyForZ1;

  if ( particle->GetDefinition() == instance->GetIon("alpha++")
      || particle->GetDefinition() == instance->GetIon("alpha+")
      || particle->GetDefinition() == instance->GetIon("helium")
  )

  lowLim = killBelowEnergyForZ2;

  G4double k = particle->GetKineticEnergy();

  const G4String& particleName = particle->GetDefinition()->GetParticleName();

  // SI - the following is useless since lowLim is already defined
  /*
   std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
   pos1 = lowEnergyLimit.find(particleName);

   if (pos1 != lowEnergyLimit.end())
   {
   lowLim = pos1->second;
   }
   */

  std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
  pos2 = highEnergyLimit.find(particleName);

  if (pos2 != highEnergyLimit.end())
  {
    highLim = pos2->second;
  }

  if (k >= lowLim && k <= highLim)
  {
    G4ParticleDefinition* definition = particle->GetDefinition();
    G4ParticleMomentum primaryDirection = particle->GetMomentumDirection();
    /*
     G4double particleMass = definition->GetPDGMass();
     G4double totalEnergy = k + particleMass;
     G4double pSquare = k*(totalEnergy+particleMass);
     G4double totalMomentum = std::sqrt(pSquare);
     */

    G4int ionizationShell = RandomSelect(k,particleName);

    // sample deexcitation
    // here we assume that H_{2}O electronic levels are the same of Oxigen.
    // this can be considered true with a rough 10% error in energy on K-shell,

    G4int secNumberInit = 0;// need to know at a certain point the enrgy of secondaries
    G4int secNumberFinal = 0;// So I'll make the diference and then sum the energies

    G4double bindingEnergy = 0;
    bindingEnergy = waterStructure.IonisationEnergy(ionizationShell);

    //SI: additional protection if tcs interpolation method is modified
    if (k<bindingEnergy) return;
    //

    G4int Z = 8;
    if(fAtomDeexcitation)
    {
      G4AtomicShellEnumerator as = fKShell;

      if (ionizationShell <5 && ionizationShell >1)
      {
        as = G4AtomicShellEnumerator(4-ionizationShell);
      }
      else if (ionizationShell <2)
      {
        as = G4AtomicShellEnumerator(3);
      }

      //	DEBUG
      //	if (ionizationShell == 4) {
      //
      //	  G4cout << "Z: " << Z << " as: " << as
      //               << " ionizationShell: " << ionizationShell << " bindingEnergy: "<< bindingEnergy/eV << G4endl;
      //	  G4cout << "Press <Enter> key to continue..." << G4endl;
      //	  G4cin.ignore();
      //	}

      const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
      secNumberInit = fvect->size();
      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, 0, 0);
      secNumberFinal = fvect->size();
    }

    G4double secondaryKinetic = RandomizeEjectedElectronEnergy(definition,k,ionizationShell);

    G4ThreeVector deltaDirection =
    GetAngularDistribution()->SampleDirectionForShell(particle, secondaryKinetic,
        Z, ionizationShell,
        couple->GetMaterial());

    // Ignored for ions on electrons
    /*
     G4double deltaTotalMomentum = std::sqrt(secondaryKinetic*(secondaryKinetic + 2.*electron_mass_c2 ));

     G4double finalPx = totalMomentum*primaryDirection.x() - deltaTotalMomentum*deltaDirection.x();
     G4double finalPy = totalMomentum*primaryDirection.y() - deltaTotalMomentum*deltaDirection.y();
     G4double finalPz = totalMomentum*primaryDirection.z() - deltaTotalMomentum*deltaDirection.z();
     G4double finalMomentum = std::sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz);
     finalPx /= finalMomentum;
     finalPy /= finalMomentum;
     finalPz /= finalMomentum;

     G4ThreeVector direction;
     direction.set(finalPx,finalPy,finalPz);

     fParticleChangeForGamma->ProposeMomentumDirection(direction.unit()) ;
     */
    fParticleChangeForGamma->ProposeMomentumDirection(primaryDirection);

    G4double scatteredEnergy = k-bindingEnergy-secondaryKinetic;
    G4double deexSecEnergy = 0;
    for (G4int j=secNumberInit; j < secNumberFinal; j++)
    {

      deexSecEnergy = deexSecEnergy + (*fvect)[j]->GetKineticEnergy();

    }

    if (!statCode)
    {
      fParticleChangeForGamma->SetProposedKineticEnergy(scatteredEnergy);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(k-scatteredEnergy-secondaryKinetic-deexSecEnergy);
    }
    else
    {
      fParticleChangeForGamma->SetProposedKineticEnergy(k);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(k-scatteredEnergy);
    }

    // debug
    // k-scatteredEnergy-secondaryKinetic-deexSecEnergy = k-(k-bindingEnergy-secondaryKinetic)-secondaryKinetic-deexSecEnergy =
    // = k-k+bindingEnergy+secondaryKinetic-secondaryKinetic-deexSecEnergy=
    // = bindingEnergy-deexSecEnergy
    // SO deexSecEnergy=0 => LocalEnergyDeposit =  bindingEnergy

    G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),deltaDirection,secondaryKinetic);
    fvect->push_back(dp);

    const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
    G4DNAChemistryManager::Instance()->CreateWaterMolecule(eIonizedMolecule,
        ionizationShell,
        theIncomingTrack);
  }

  // SI - not useful since low energy of model is 0 eV

  if (k < lowLim)
  {
    fParticleChangeForGamma->SetProposedKineticEnergy(0.);
    fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(k);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationModel::RandomizeEjectedElectronEnergy(G4ParticleDefinition* particleDefinition,
                                                                  G4double k,
                                                                  G4int shell)
{
  G4double maximumKineticEnergyTransfer = 0.;

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  if (particleDefinition == G4Proton::ProtonDefinition()
      || particleDefinition == instance->GetIon("hydrogen"))
  {
    maximumKineticEnergyTransfer = 4. * (electron_mass_c2 / proton_mass_c2) * k;
  }

  else if (particleDefinition == instance->GetIon("helium")
      || particleDefinition == instance->GetIon("alpha+")
      || particleDefinition == instance->GetIon("alpha++"))
  {
    maximumKineticEnergyTransfer = 4. * (0.511 / 3728) * k;
  }

  G4double crossSectionMaximum = 0.;

  for (G4double value = waterStructure.IonisationEnergy(shell);
      value <= 5. * waterStructure.IonisationEnergy(shell) && k >= value;
      value += 0.1 * eV)
  {
    G4double differentialCrossSection =
        DifferentialCrossSection(particleDefinition, k, value, shell);
    if (differentialCrossSection >= crossSectionMaximum)
      crossSectionMaximum = differentialCrossSection;
  }

  G4double secElecKinetic = 0.;

  do
  {
    secElecKinetic = G4UniformRand()* maximumKineticEnergyTransfer;
  }while(G4UniformRand()*crossSectionMaximum > DifferentialCrossSection(particleDefinition,
          k,
          secElecKinetic+waterStructure.IonisationEnergy(shell),
          shell));

  return (secElecKinetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// The following section is not used anymore but is kept for memory
// GetAngularDistribution()->SampleDirectionForShell is used instead

/*
 void G4DNARuddIonisationModel::RandomizeEjectedElectronDirection(G4ParticleDefinition* particleDefinition,
 G4double k,
 G4double secKinetic,
 G4double & cosTheta,
 G4double & phi )
 {
 G4DNAGenericIonsManager *instance;
 instance = G4DNAGenericIonsManager::Instance();

 G4double maxSecKinetic = 0.;

 if (particleDefinition == G4Proton::ProtonDefinition()
 || particleDefinition == instance->GetIon("hydrogen"))
 {
 maxSecKinetic = 4.* (electron_mass_c2 / proton_mass_c2) * k;
 }

 else if (particleDefinition == instance->GetIon("helium")
 || particleDefinition == instance->GetIon("alpha+")
 || particleDefinition == instance->GetIon("alpha++"))
 {
 maxSecKinetic = 4.* (0.511 / 3728) * k;
 }

 phi = twopi * G4UniformRand();

 // cosTheta = std::sqrt(secKinetic / maxSecKinetic);

 // Restriction below 100 eV from Emfietzoglou (2000)

 if (secKinetic>100*eV) cosTheta = std::sqrt(secKinetic / maxSecKinetic);
 else cosTheta = (2.*G4UniformRand())-1.;

 }
 */
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4DNARuddIonisationModel::DifferentialCrossSection(G4ParticleDefinition* particleDefinition,
                                                            G4double k,
                                                            G4double energyTransfer,
                                                            G4int ionizationLevelIndex)
{
  // Shells ids are 0 1 2 3 4 (4 is k shell)
  // !!Attention, "energyTransfer" here is the energy transfered to the electron which means
  //             that the secondary kinetic energy is w = energyTransfer - bindingEnergy
  //
  //   ds            S                F1(nu) + w * F2(nu)
  //  ---- = G(k) * ----     -------------------------------------------
  //   dw            Bj       (1+w)^3 * [1 + exp{alpha * (w - wc) / nu}]
  //
  // w is the secondary electron kinetic Energy in eV
  //
  // All the other parameters can be found in Rudd's Papers
  //
  // M.Eugene Rudd, 1988, User-Friendly model for the energy distribution of
  // electrons from protons or electron collisions. Nucl. Tracks Rad. Meas.Vol 16 N0 2/3 pp 219-218
  //

  const G4int j = ionizationLevelIndex;

  G4double A1;
  G4double B1;
  G4double C1;
  G4double D1;
  G4double E1;
  G4double A2;
  G4double B2;
  G4double C2;
  G4double D2;
  G4double alphaConst;

  //  const G4double Bj[5] = {12.61*eV, 14.73*eV, 18.55*eV, 32.20*eV, 539.7*eV};
  // The following values are provided by M. dingfelder (priv. comm)
  const G4double Bj[5] = { 12.60 * eV, 14.70 * eV, 18.40 * eV, 32.20 * eV, 540
      * eV };

  if (j == 4)
  {
    //Data For Liquid Water K SHELL from Dingfelder (Protons in Water)
    A1 = 1.25;
    B1 = 0.5;
    C1 = 1.00;
    D1 = 1.00;
    E1 = 3.00;
    A2 = 1.10;
    B2 = 1.30;
    C2 = 1.00;
    D2 = 0.00;
    alphaConst = 0.66;
  } else
  {
    //Data For Liquid Water from Dingfelder (Protons in Water)
    A1 = 1.02;
    B1 = 82.0;
    C1 = 0.45;
    D1 = -0.80;
    E1 = 0.38;
    A2 = 1.07;
    // Value provided by M. Dingfelder (priv. comm)
    B2 = 11.6;
    //
    C2 = 0.60;
    D2 = 0.04;
    alphaConst = 0.64;
  }

  const G4double n = 2.;
  const G4double Gj[5] = { 0.99, 1.11, 1.11, 0.52, 1. };

  G4DNAGenericIonsManager* instance;
  instance = G4DNAGenericIonsManager::Instance();

  G4double wBig = (energyTransfer
      - waterStructure.IonisationEnergy(ionizationLevelIndex));
  if (wBig < 0)
    return 0.;

  G4double w = wBig / Bj[ionizationLevelIndex];
  // Note that the following (j==4) cases are provided by   M. Dingfelder (priv. comm)
  if (j == 4)
    w = wBig / waterStructure.IonisationEnergy(ionizationLevelIndex);

  G4double Ry = 13.6 * eV;

  G4double tau = 0.;

  G4bool isProtonOrHydrogen = false;
  G4bool isHelium = false;

  if (particleDefinition == G4Proton::ProtonDefinition()
      || particleDefinition == instance->GetIon("hydrogen"))
  {
    isProtonOrHydrogen = true;
    tau = (electron_mass_c2 / proton_mass_c2) * k;
  }

  else if (particleDefinition == instance->GetIon("helium")
      || particleDefinition == instance->GetIon("alpha+")
      || particleDefinition == instance->GetIon("alpha++"))
  {
    isHelium = true;
    tau = (0.511 / 3728.) * k;
  }

  G4double S = 4. * pi * Bohr_radius * Bohr_radius * n
      * std::pow((Ry / Bj[ionizationLevelIndex]), 2);
  if (j == 4)
    S = 4. * pi * Bohr_radius * Bohr_radius * n
        * std::pow((Ry / waterStructure.IonisationEnergy(ionizationLevelIndex)),
                   2);

  G4double v2 = tau / Bj[ionizationLevelIndex];
  if (j == 4)
    v2 = tau / waterStructure.IonisationEnergy(ionizationLevelIndex);

  G4double v = std::sqrt(v2);
  G4double wc = 4. * v2 - 2. * v - (Ry / (4. * Bj[ionizationLevelIndex]));
  if (j == 4)
    wc = 4. * v2 - 2. * v
        - (Ry / (4. * waterStructure.IonisationEnergy(ionizationLevelIndex)));

  G4double L1 = (C1 * std::pow(v, (D1))) / (1. + E1 * std::pow(v, (D1 + 4.)));
  G4double L2 = C2 * std::pow(v, (D2));
  G4double H1 = (A1 * std::log(1. + v2)) / (v2 + (B1 / v2));
  G4double H2 = (A2 / v2) + (B2 / (v2 * v2));

  G4double F1 = L1 + H1;
  G4double F2 = (L2 * H2) / (L2 + H2);

  G4double sigma =
      CorrectionFactor(particleDefinition, k) * Gj[j]
          * (S / Bj[ionizationLevelIndex])
          * ((F1 + w * F2)
              / (std::pow((1. + w), 3)
                  * (1. + G4Exp(alphaConst * (w - wc) / v))));

  if (j == 4)
    sigma = CorrectionFactor(particleDefinition, k) * Gj[j]
        * (S / waterStructure.IonisationEnergy(ionizationLevelIndex))
        * ((F1 + w * F2)
            / (std::pow((1. + w), 3)
                * (1. + G4Exp(alphaConst * (w - wc) / v))));

  if ((particleDefinition == instance->GetIon("hydrogen"))
      && (ionizationLevelIndex == 4))
  {
    //    sigma = Gj[j] * (S/Bj[ionizationLevelIndex])
    sigma = Gj[j] * (S / waterStructure.IonisationEnergy(ionizationLevelIndex))
        * ((F1 + w * F2)
            / (std::pow((1. + w), 3)
                * (1. + G4Exp(alphaConst * (w - wc) / v))));
  }
//    if (    particleDefinition == G4Proton::ProtonDefinition()
//            || particleDefinition == instance->GetIon("hydrogen")
//            )
  if (isProtonOrHydrogen)
  {
    return (sigma);
  }

  if (particleDefinition == instance->GetIon("alpha++"))
  {
    slaterEffectiveCharge[0] = 0.;
    slaterEffectiveCharge[1] = 0.;
    slaterEffectiveCharge[2] = 0.;
    sCoefficient[0] = 0.;
    sCoefficient[1] = 0.;
    sCoefficient[2] = 0.;
  }

  else if (particleDefinition == instance->GetIon("alpha+"))
  {
    slaterEffectiveCharge[0] = 2.0;
    // The following values are provided by M. Dingfelder (priv. comm)
    slaterEffectiveCharge[1] = 2.0;
    slaterEffectiveCharge[2] = 2.0;
    //
    sCoefficient[0] = 0.7;
    sCoefficient[1] = 0.15;
    sCoefficient[2] = 0.15;
  }

  else if (particleDefinition == instance->GetIon("helium"))
  {
    slaterEffectiveCharge[0] = 1.7;
    slaterEffectiveCharge[1] = 1.15;
    slaterEffectiveCharge[2] = 1.15;
    sCoefficient[0] = 0.5;
    sCoefficient[1] = 0.25;
    sCoefficient[2] = 0.25;
  }

//    if (    particleDefinition == instance->GetIon("helium")
//            || particleDefinition == instance->GetIon("alpha+")
//            || particleDefinition == instance->GetIon("alpha++")
//            )
  if (isHelium)
  {
    sigma = Gj[j] * (S / Bj[ionizationLevelIndex])
        * ((F1 + w * F2)
            / (std::pow((1. + w), 3)
                * (1. + G4Exp(alphaConst * (w - wc) / v))));

    if (j == 4)
      sigma = Gj[j]
          * (S / waterStructure.IonisationEnergy(ionizationLevelIndex))
          * ((F1 + w * F2)
              / (std::pow((1. + w), 3)
                  * (1. + G4Exp(alphaConst * (w - wc) / v))));

    G4double zEff = particleDefinition->GetPDGCharge() / eplus
        + particleDefinition->GetLeptonNumber();

    zEff -= (sCoefficient[0]
        * S_1s(k, energyTransfer, slaterEffectiveCharge[0], 1.)
        + sCoefficient[1]
            * S_2s(k, energyTransfer, slaterEffectiveCharge[1], 2.)
        + sCoefficient[2]
            * S_2p(k, energyTransfer, slaterEffectiveCharge[2], 2.));

    return zEff * zEff * sigma;
  }

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationModel::S_1s(G4double t,
                                        G4double energyTransferred,
                                        G4double slaterEffectiveChg,
                                        G4double shellNumber)
{
  // 1 - e^(-2r) * ( 1 + 2 r + 2 r^2)
  // Dingfelder, in Chattanooga 2005 proceedings, formula (7)

  G4double r = R(t, energyTransferred, slaterEffectiveChg, shellNumber);
  G4double value = 1. - G4Exp(-2 * r) * ((2. * r + 2.) * r + 1.);

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationModel::S_2s(G4double t,
                                        G4double energyTransferred,
                                        G4double slaterEffectiveChg,
                                        G4double shellNumber)
{
  // 1 - e^(-2 r) * ( 1 + 2 r + 2 r^2 + 2 r^4)
  // Dingfelder, in Chattanooga 2005 proceedings, formula (8)

  G4double r = R(t, energyTransferred, slaterEffectiveChg, shellNumber);
  G4double value = 1.
      - G4Exp(-2 * r) * (((2. * r * r + 2.) * r + 2.) * r + 1.);

  return value;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationModel::S_2p(G4double t,
                                        G4double energyTransferred,
                                        G4double slaterEffectiveChg,
                                        G4double shellNumber)
{
  // 1 - e^(-2 r) * ( 1 + 2 r + 2 r^2 + 4/3 r^3 + 2/3 r^4)
  // Dingfelder, in Chattanooga 2005 proceedings, formula (9)

  G4double r = R(t, energyTransferred, slaterEffectiveChg, shellNumber);
  G4double value = 1.
      - G4Exp(-2 * r)
          * ((((2. / 3. * r + 4. / 3.) * r + 2.) * r + 2.) * r + 1.);

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationModel::R(G4double t,
                                     G4double energyTransferred,
                                     G4double slaterEffectiveChg,
                                     G4double shellNumber)
{
  // tElectron = m_electron / m_alpha * t
  // Dingfelder, in Chattanooga 2005 proceedings, p 4

  G4double tElectron = 0.511 / 3728. * t;
  // The following values are provided by M. Dingfelder (priv. comm)
  G4double H = 2. * 13.60569172 * eV;
  G4double value = std::sqrt(2. * tElectron / H) / (energyTransferred / H)
      * (slaterEffectiveChg / shellNumber);

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationModel::CorrectionFactor(G4ParticleDefinition* particleDefinition,
                                                    G4double k)
{
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  if (particleDefinition == G4Proton::Proton())
  {
    return (1.);
  } else if (particleDefinition == instance->GetIon("hydrogen"))
  {
    G4double value = (std::log10(k / eV) - 4.2) / 0.5;
    // The following values are provided by M. Dingfelder (priv. comm)
    return ((0.6 / (1 + G4Exp(value))) + 0.9);
  } else
  {
    return (1.);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNARuddIonisationModel::RandomSelect(G4double k,
                                             const G4String& particle)
{

  // BEGIN PART 1/2 OF ELECTRON CORRECTION

  // add ONE or TWO electron-water ionisation for alpha+ and helium

  G4int level = 0;

  // Retrieve data table corresponding to the current particle type

  std::map<G4String, G4DNACrossSectionDataSet*, std::less<G4String> >::iterator pos;
  pos = tableData.find(particle);

  if (pos != tableData.end())
  {
    G4DNACrossSectionDataSet* table = pos->second;

    if (table != 0)
    {
      G4double* valuesBuffer = new G4double[table->NumberOfComponents()];

      const size_t n(table->NumberOfComponents());
      size_t i(n);
      G4double value = 0.;

      while (i > 0)
      {
        i--;
        valuesBuffer[i] = table->GetComponent(i)->FindValue(k);
        value += valuesBuffer[i];
      }

      value *= G4UniformRand();

      i = n;

      while (i > 0)
      {
        i--;

        if (valuesBuffer[i] > value)
        {
          delete[] valuesBuffer;
          return i;
        }
        value -= valuesBuffer[i];
      }

      if (valuesBuffer)
        delete[] valuesBuffer;

    }
  } else
  {
    G4Exception("G4DNARuddIonisationModel::RandomSelect",
                "em0002",
                FatalException,
                "Model not applicable to particle type.");
  }

  return level;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationModel::PartialCrossSection(const G4Track& track)
{
  G4double sigma = 0.;

  const G4DynamicParticle* particle = track.GetDynamicParticle();
  G4double k = particle->GetKineticEnergy();

  G4double lowLim = 0;
  G4double highLim = 0;

  const G4String& particleName = particle->GetDefinition()->GetParticleName();

  std::map<G4String, G4double, std::less<G4String> >::iterator pos1;
  pos1 = lowEnergyLimit.find(particleName);

  if (pos1 != lowEnergyLimit.end())
  {
    lowLim = pos1->second;
  }

  std::map<G4String, G4double, std::less<G4String> >::iterator pos2;
  pos2 = highEnergyLimit.find(particleName);

  if (pos2 != highEnergyLimit.end())
  {
    highLim = pos2->second;
  }

  if (k >= lowLim && k <= highLim)
  {
    std::map<G4String, G4DNACrossSectionDataSet*, std::less<G4String> >::iterator pos;
    pos = tableData.find(particleName);

    if (pos != tableData.end())
    {
      G4DNACrossSectionDataSet* table = pos->second;
      if (table != 0)
      {
        sigma = table->FindValue(k);
      }
    } else
    {
      G4Exception("G4DNARuddIonisationModel::PartialCrossSection",
                  "em0002",
                  FatalException,
                  "Model not applicable to particle type.");
    }
  }

  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARuddIonisationModel::Sum(G4double /* energy */,
                                       const G4String& /* particle */)
{
  return 0;
}

