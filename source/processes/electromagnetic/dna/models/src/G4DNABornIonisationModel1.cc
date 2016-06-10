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
// $Id: G4DNABornIonisationModel1.cc 87631 2014-12-14 12:42:05Z matkara $
//

#include "G4DNABornIonisationModel1.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4DNABornAngle.hh"
#include "G4DeltaAngle.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNABornIonisationModel1::G4DNABornIonisationModel1(const G4ParticleDefinition*,
                                                     const G4String& nam) :
    G4VEmModel(nam), isInitialised(false)
{
  verboseLevel = 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if (verboseLevel > 0)
  {
    G4cout << "Born ionisation model is constructed " << G4endl;
  }

  //Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);
  fAtomDeexcitation = 0;
  fParticleChangeForGamma = 0;
  fpMolWaterDensity = 0;

  // define default angular generator
  SetAngularDistribution(new G4DNABornAngle());

  // Selection of computation method

  fasterCode = false;

  // Selection of stationary mode

  statCode = false;

  // Selection of SP scaling

  spScaling = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNABornIonisationModel1::~G4DNABornIonisationModel1()
{
  // Cross section

  std::map<G4String, G4DNACrossSectionDataSet*, std::less<G4String> >::iterator pos;
  for (pos = tableData.begin(); pos != tableData.end(); ++pos)
  {
    G4DNACrossSectionDataSet* table = pos->second;
    delete table;
  }

  // Final state

  eVecm.clear();
  pVecm.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNABornIonisationModel1::Initialise(const G4ParticleDefinition* particle,
                                           const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
  {
    G4cout << "Calling G4DNABornIonisationModel1::Initialise()" << G4endl;
  }

  // Energy limits

  G4String fileElectron("dna/sigma_ionisation_e_born");
  G4String fileProton("dna/sigma_ionisation_p_born");

  G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();
  G4ParticleDefinition* protonDef = G4Proton::ProtonDefinition();

  G4String electron;
  G4String proton;

  G4double scaleFactor = (1.e-22 / 3.343) * m*m;

  char *path = getenv("G4LEDATA");

  // *** ELECTRON

  electron = electronDef->GetParticleName();

  tableFile[electron] = fileElectron;

  lowEnergyLimit[electron] = 11. * eV;
  highEnergyLimit[electron] = 1. * MeV;

  // Cross section

  G4DNACrossSectionDataSet* tableE = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
  tableE->LoadData(fileElectron);

  tableData[electron] = tableE;

  // Final state

  std::ostringstream eFullFileName;

  if (fasterCode) eFullFileName << path << "/dna/sigmadiff_cumulated_ionisation_e_born_hp.dat";
  if (!fasterCode) eFullFileName << path << "/dna/sigmadiff_ionisation_e_born.dat";

  std::ifstream eDiffCrossSection(eFullFileName.str().c_str());

  if (!eDiffCrossSection)
  {
    if (fasterCode) G4Exception("G4DNABornIonisationModel1::Initialise","em0003",
        FatalException,"Missing data file:/dna/sigmadiff_cumulated_ionisation_e_born_hp.dat");

    if (!fasterCode) G4Exception("G4DNABornIonisationModel1::Initialise","em0003",
        FatalException,"Missing data file:/dna/sigmadiff_ionisation_e_born.dat");
  }

  //

  // Clear the arrays for re-initialization case (MT mode)
  // March 25th, 2014 - Vaclav Stepan, Sebastien Incerti

  eTdummyVec.clear();
  pTdummyVec.clear();

  eVecm.clear();
  pVecm.clear();

  for (int j=0; j<5; j++)
  {
    eProbaShellMap[j].clear();
    pProbaShellMap[j].clear();

    eDiffCrossSectionData[j].clear();
    pDiffCrossSectionData[j].clear();

    eNrjTransfData[j].clear();
    pNrjTransfData[j].clear();
  }
  //

  eTdummyVec.push_back(0.);
  while(!eDiffCrossSection.eof())
  {
    double tDummy;
    double eDummy;
    eDiffCrossSection>>tDummy>>eDummy;
    if (tDummy != eTdummyVec.back()) eTdummyVec.push_back(tDummy);

    double tmp;
    for (int j=0; j<5; j++)
    {
      eDiffCrossSection>> tmp;

      eDiffCrossSectionData[j][tDummy][eDummy] = tmp;

      if (fasterCode)
      {
        eNrjTransfData[j][tDummy][eDiffCrossSectionData[j][tDummy][eDummy]]=eDummy;
        eProbaShellMap[j][tDummy].push_back(eDiffCrossSectionData[j][tDummy][eDummy]);
      }

      // SI - only if eof is not reached
      if (!eDiffCrossSection.eof() && !fasterCode) eDiffCrossSectionData[j][tDummy][eDummy]*=scaleFactor;

      if (!fasterCode) eVecm[tDummy].push_back(eDummy);

    }
  }

  // *** PROTON

  proton = protonDef->GetParticleName();

  tableFile[proton] = fileProton;

  lowEnergyLimit[proton] = 500. * keV;
  highEnergyLimit[proton] = 100. * MeV;

  // Cross section

  G4DNACrossSectionDataSet* tableP = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
  tableP->LoadData(fileProton);

  tableData[proton] = tableP;

  // Final state

  std::ostringstream pFullFileName;

  if (fasterCode) pFullFileName << path << "/dna/sigmadiff_cumulated_ionisation_p_born_hp.dat";

  if (!fasterCode) pFullFileName << path << "/dna/sigmadiff_ionisation_p_born.dat";

  std::ifstream pDiffCrossSection(pFullFileName.str().c_str());

  if (!pDiffCrossSection)
  {
    if (fasterCode) G4Exception("G4DNABornIonisationModel1::Initialise","em0003",
        FatalException,"Missing data file:/dna/sigmadiff_cumulated_ionisation_p_born_hp.dat");

    if (!fasterCode) G4Exception("G4DNABornIonisationModel1::Initialise","em0003",
        FatalException,"Missing data file:/dna/sigmadiff_ionisation_p_born.dat");
  }

  pTdummyVec.push_back(0.);
  while(!pDiffCrossSection.eof())
  {
    double tDummy;
    double eDummy;
    pDiffCrossSection>>tDummy>>eDummy;
    if (tDummy != pTdummyVec.back()) pTdummyVec.push_back(tDummy);
    for (int j=0; j<5; j++)
    {
      pDiffCrossSection>>pDiffCrossSectionData[j][tDummy][eDummy];

      if (fasterCode)
      {
        pNrjTransfData[j][tDummy][pDiffCrossSectionData[j][tDummy][eDummy]]=eDummy;
        pProbaShellMap[j][tDummy].push_back(pDiffCrossSectionData[j][tDummy][eDummy]);
      }

      // SI - only if eof is not reached !
      if (!pDiffCrossSection.eof() && !fasterCode) pDiffCrossSectionData[j][tDummy][eDummy]*=scaleFactor;

      if (!fasterCode) pVecm[tDummy].push_back(eDummy);
    }
  }

  //

  if (particle==electronDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[electron]);
    SetHighEnergyLimit(highEnergyLimit[electron]);
  }

  if (particle==protonDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[proton]);
    SetHighEnergyLimit(highEnergyLimit[proton]);
  }

  if( verboseLevel>0 )
  {
    G4cout << "Born ionisation model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / eV << " eV - "
    << HighEnergyLimit() / keV << " keV for "
    << particle->GetParticleName()
    << G4endl;
  }

  // Initialize water density pointer
  fpMolWaterDensity = G4DNAMolecularMaterial::Instance()->
  GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

  //
  fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();

  if (isInitialised)
  { return;}
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNABornIonisationModel1::CrossSectionPerVolume(const G4Material* material,
                                                          const G4ParticleDefinition* particleDefinition,
                                                          G4double ekin,
                                                          G4double,
                                                          G4double)
{
  if (verboseLevel > 3)
  {
    G4cout << "Calling CrossSectionPerVolume() of G4DNABornIonisationModel1"
        << G4endl;

  }

  if (
      particleDefinition != G4Proton::ProtonDefinition()
      &&
      particleDefinition != G4Electron::ElectronDefinition()
  )

  return 0;

  // Calculate total cross section for model

  G4double lowLim = 0;
  G4double highLim = 0;
  G4double sigma=0;

  G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];

  if(waterDensity!= 0.0)
  //  if (material == nistwater || material->GetBaseMaterial() == nistwater)
  {
    const G4String& particleName = particleDefinition->GetParticleName();

    std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
    pos1 = lowEnergyLimit.find(particleName);
    if (pos1 != lowEnergyLimit.end())
    {
      lowLim = pos1->second;
    }

    std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
    pos2 = highEnergyLimit.find(particleName);
    if (pos2 != highEnergyLimit.end())
    {
      highLim = pos2->second;
    }

    if (ekin >= lowLim && ekin < highLim)
    {
      std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
      pos = tableData.find(particleName);

      if (pos != tableData.end())
      {
        G4DNACrossSectionDataSet* table = pos->second;
        if (table != 0)
        {
          sigma = table->FindValue(ekin);

          // ICRU49 electronic SP scaling - ZF, SI

          if (particleDefinition == G4Proton::ProtonDefinition() && ekin < 70*MeV && spScaling)
          {
           G4double A = 1.39241700556072800000E-009 ;
           G4double B = -8.52610412942622630000E-002 ;
           sigma = sigma * std::exp(A*(ekin/eV)+B);
          }
          //

        }
      }
      else
      {
        G4Exception("G4DNABornIonisationModel1::CrossSectionPerVolume","em0002",
            FatalException,"Model not applicable to particle type.");
      }
    }

    if (verboseLevel > 2)
    {
      G4cout << "__________________________________" << G4endl;
      G4cout << "G4DNABornIonisationModel1 - XS INFO START" << G4endl;
      G4cout << "Kinetic energy(eV)=" << ekin/eV << " particle : " << particleName << G4endl;
      G4cout << "Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
      G4cout << "Cross section per water molecule (cm^-1)=" << sigma*waterDensity/(1./cm) << G4endl;
      G4cout << "G4DNABornIonisationModel1 - XS INFO END" << G4endl;
    }
  } // if (waterMaterial)

  return sigma*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNABornIonisationModel1::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                  const G4MaterialCutsCouple* couple,
                                                  const G4DynamicParticle* particle,
                                                  G4double,
                                                  G4double)
{

  if (verboseLevel > 3)
  {
    G4cout << "Calling SampleSecondaries() of G4DNABornIonisationModel1"
        << G4endl;
  }

  G4double lowLim = 0;
  G4double highLim = 0;

  G4double k = particle->GetKineticEnergy();

  const G4String& particleName = particle->GetDefinition()->GetParticleName();

  std::map< G4String,G4double,std::less<G4String> >::iterator pos1;
  pos1 = lowEnergyLimit.find(particleName);

  if (pos1 != lowEnergyLimit.end())
  {
    lowLim = pos1->second;
  }

  std::map< G4String,G4double,std::less<G4String> >::iterator pos2;
  pos2 = highEnergyLimit.find(particleName);

  if (pos2 != highEnergyLimit.end())
  {
    highLim = pos2->second;
  }

  if (k >= lowLim && k < highLim)
  {
    G4ParticleMomentum primaryDirection = particle->GetMomentumDirection();
    G4double particleMass = particle->GetDefinition()->GetPDGMass();
    G4double totalEnergy = k + particleMass;
    G4double pSquare = k * (totalEnergy + particleMass);
    G4double totalMomentum = std::sqrt(pSquare);

    G4int ionizationShell = 0;

    if (!fasterCode) ionizationShell = RandomSelect(k,particleName);

    // SI: The following protection is necessary to avoid infinite loops :
    //  sigmadiff_ionisation_e_born.dat has non zero partial xs at 18 eV for shell 3 (ionizationShell ==2)
    //  sigmadiff_cumulated_ionisation_e_born.dat has zero cumulated partial xs at 18 eV for shell 3 (ionizationShell ==2)
    //  this is due to the fact that the max allowed transfered energy is (18+10.79)/2=17.025 eV and only transfered energies
    //  strictly above this value have non zero partial xs in sigmadiff_ionisation_e_born.dat (starting at trans = 17.12 eV)

    if (fasterCode)
    do
    {
      ionizationShell = RandomSelect(k,particleName);
    }while (k<19*eV && ionizationShell==2 && particle->GetDefinition()==G4Electron::ElectronDefinition());

    // AM: sample deexcitation
    // here we assume that H_{2}O electronic levels are the same as Oxygen.
    // this can be considered true with a rough 10% error in energy on K-shell,

    G4int secNumberInit = 0;// need to know at a certain point the energy of secondaries
    G4int secNumberFinal = 0;// So I'll make the diference and then sum the energies

    G4double bindingEnergy = 0;
    bindingEnergy = waterStructure.IonisationEnergy(ionizationShell);

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

      //	FOR DEBUG ONLY
      //	if (ionizationShell == 4) {
      //
      //	  G4cout << "Z: " << Z << " as: " << as
      //               << " ionizationShell: " << ionizationShell << " bindingEnergy: "<< bindingEnergy/eV << G4endl;
      //        G4cout << "Press <Enter> key to continue..." << G4endl;
      //	  G4cin.ignore();
      //	}

      const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
      secNumberInit = fvect->size();
      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, 0, 0);
      secNumberFinal = fvect->size();
    }

    G4double secondaryKinetic=-1000*eV;

    if (fasterCode == false)
    {
      secondaryKinetic = RandomizeEjectedElectronEnergy(particle->GetDefinition(),k,ionizationShell);
    }
    // SI - 01/04/2014
    else
    {
      secondaryKinetic = RandomizeEjectedElectronEnergyFromCumulatedDcs(particle->GetDefinition(),k,ionizationShell);
    }
    //

    G4ThreeVector deltaDirection =
    GetAngularDistribution()->SampleDirectionForShell(particle, secondaryKinetic,
        Z, ionizationShell,
        couple->GetMaterial());

    if (particle->GetDefinition() == G4Electron::ElectronDefinition())
    {
      G4double deltaTotalMomentum = std::sqrt(secondaryKinetic*(secondaryKinetic + 2.*electron_mass_c2 ));

      G4double finalPx = totalMomentum*primaryDirection.x() - deltaTotalMomentum*deltaDirection.x();
      G4double finalPy = totalMomentum*primaryDirection.y() - deltaTotalMomentum*deltaDirection.y();
      G4double finalPz = totalMomentum*primaryDirection.z() - deltaTotalMomentum*deltaDirection.z();
      G4double finalMomentum = std::sqrt(finalPx*finalPx + finalPy*finalPy + finalPz*finalPz);
      finalPx /= finalMomentum;
      finalPy /= finalMomentum;
      finalPz /= finalMomentum;

      G4ThreeVector direction;
      direction.set(finalPx,finalPy,finalPz);

      fParticleChangeForGamma->ProposeMomentumDirection(direction.unit());
    }

    else fParticleChangeForGamma->ProposeMomentumDirection(primaryDirection);

    // note that secondaryKinetic is the energy of the delta ray, not of all secondaries.
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
    
    // SI - 01/04/2014
    if (secondaryKinetic>0)
    {
      G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),deltaDirection,secondaryKinetic);
      fvect->push_back(dp);
    }
    //

    const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
    G4DNAChemistryManager::Instance()->CreateWaterMolecule(eIonizedMolecule,
        ionizationShell,
        theIncomingTrack);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNABornIonisationModel1::RandomizeEjectedElectronEnergy(G4ParticleDefinition* particleDefinition,
                                                                   G4double k,
                                                                   G4int shell)
{
  // G4cout << "*** SLOW computation for " << " " << particleDefinition->GetParticleName() << G4endl;

  if (particleDefinition == G4Electron::ElectronDefinition())
  {
    G4double maximumEnergyTransfer = 0.;
    if ((k + waterStructure.IonisationEnergy(shell)) / 2. > k)
      maximumEnergyTransfer = k;
    else
      maximumEnergyTransfer = (k + waterStructure.IonisationEnergy(shell)) / 2.;

    // SI : original method
    /*
     G4double crossSectionMaximum = 0.;
     for(G4double value=waterStructure.IonisationEnergy(shell); value<=maximumEnergyTransfer; value+=0.1*eV)
     {
     G4double differentialCrossSection = DifferentialCrossSection(particleDefinition, k/eV, value/eV, shell);
     if(differentialCrossSection >= crossSectionMaximum) crossSectionMaximum = differentialCrossSection;
     }
     */

    // SI : alternative method
    G4double crossSectionMaximum = 0.;

    G4double minEnergy = waterStructure.IonisationEnergy(shell);
    G4double maxEnergy = maximumEnergyTransfer;
    G4int nEnergySteps = 50;

    G4double value(minEnergy);
    G4double stpEnergy(std::pow(maxEnergy / value,
                                1. / static_cast<G4double>(nEnergySteps - 1)));
    G4int step(nEnergySteps);
    while (step > 0)
    {
      step--;
      G4double differentialCrossSection =
          DifferentialCrossSection(particleDefinition,
                                   k / eV,
                                   value / eV,
                                   shell);
      if (differentialCrossSection >= crossSectionMaximum)
        crossSectionMaximum = differentialCrossSection;
      value *= stpEnergy;
    }
    //

    G4double secondaryElectronKineticEnergy = 0.;
    do
    {
      secondaryElectronKineticEnergy = G4UniformRand()* (maximumEnergyTransfer-waterStructure.IonisationEnergy(shell));
    }while(G4UniformRand()*crossSectionMaximum >
        DifferentialCrossSection(particleDefinition, k/eV,
            (secondaryElectronKineticEnergy+waterStructure.IonisationEnergy(shell))/eV,shell));

    return secondaryElectronKineticEnergy;

  }

  else if (particleDefinition == G4Proton::ProtonDefinition())
  {
    G4double maximumKineticEnergyTransfer = 4.
        * (electron_mass_c2 / proton_mass_c2) * k;

    G4double crossSectionMaximum = 0.;
    for (G4double value = waterStructure.IonisationEnergy(shell);
        value <= 4. * waterStructure.IonisationEnergy(shell); value += 0.1 * eV)
    {
      G4double differentialCrossSection =
          DifferentialCrossSection(particleDefinition,
                                   k / eV,
                                   value / eV,
                                   shell);
      if (differentialCrossSection >= crossSectionMaximum)
        crossSectionMaximum = differentialCrossSection;
    }

    G4double secondaryElectronKineticEnergy = 0.;
    do
    {
      secondaryElectronKineticEnergy = G4UniformRand()* maximumKineticEnergyTransfer;
    }while(G4UniformRand()*crossSectionMaximum >=
        DifferentialCrossSection(particleDefinition, k/eV,
            (secondaryElectronKineticEnergy+waterStructure.IonisationEnergy(shell))/eV,shell));

    return secondaryElectronKineticEnergy;
  }

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// The following section is not used anymore but is kept for memory
// GetAngularDistribution()->SampleDirectionForShell is used instead

/*
 void G4DNABornIonisationModel1::RandomizeEjectedElectronDirection(G4ParticleDefinition* particleDefinition,
 G4double k,
 G4double secKinetic,
 G4double & cosTheta,
 G4double & phi )
 {
 if (particleDefinition == G4Electron::ElectronDefinition())
 {
 phi = twopi * G4UniformRand();
 if (secKinetic < 50.*eV) cosTheta = (2.*G4UniformRand())-1.;
 else if (secKinetic <= 200.*eV)
 {
 if (G4UniformRand() <= 0.1) cosTheta = (2.*G4UniformRand())-1.;
 else cosTheta = G4UniformRand()*(std::sqrt(2.)/2);
 }
 else
 {
 G4double sin2O = (1.-secKinetic/k) / (1.+secKinetic/(2.*electron_mass_c2));
 cosTheta = std::sqrt(1.-sin2O);
 }
 }

 else if (particleDefinition == G4Proton::ProtonDefinition())
 {
 G4double maxSecKinetic = 4.* (electron_mass_c2 / proton_mass_c2) * k;
 phi = twopi * G4UniformRand();

 // cosTheta = std::sqrt(secKinetic / maxSecKinetic);

 // Restriction below 100 eV from Emfietzoglou (2000)

 if (secKinetic>100*eV) cosTheta = std::sqrt(secKinetic / maxSecKinetic);
 else cosTheta = (2.*G4UniformRand())-1.;

 }
 }
 */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double G4DNABornIonisationModel1::DifferentialCrossSection(G4ParticleDefinition * particleDefinition,
                                                           G4double k,
                                                           G4double energyTransfer,
                                                           G4int ionizationLevelIndex)
{
  G4double sigma = 0.;

  if (energyTransfer >= waterStructure.IonisationEnergy(ionizationLevelIndex))
  {
    G4double valueT1 = 0;
    G4double valueT2 = 0;
    G4double valueE21 = 0;
    G4double valueE22 = 0;
    G4double valueE12 = 0;
    G4double valueE11 = 0;

    G4double xs11 = 0;
    G4double xs12 = 0;
    G4double xs21 = 0;
    G4double xs22 = 0;

    if (particleDefinition == G4Electron::ElectronDefinition())
    {
      // k should be in eV and energy transfer eV also

      std::vector<double>::iterator t2 = std::upper_bound(eTdummyVec.begin(),
                                                          eTdummyVec.end(),
                                                          k);

      std::vector<double>::iterator t1 = t2 - 1;

      // SI : the following condition avoids situations where energyTransfer >last vector element
      if (energyTransfer <= eVecm[(*t1)].back()
          && energyTransfer <= eVecm[(*t2)].back())
      {
        std::vector<double>::iterator e12 =
            std::upper_bound(eVecm[(*t1)].begin(),
                             eVecm[(*t1)].end(),
                             energyTransfer);
        std::vector<double>::iterator e11 = e12 - 1;

        std::vector<double>::iterator e22 =
            std::upper_bound(eVecm[(*t2)].begin(),
                             eVecm[(*t2)].end(),
                             energyTransfer);
        std::vector<double>::iterator e21 = e22 - 1;

        valueT1 = *t1;
        valueT2 = *t2;
        valueE21 = *e21;
        valueE22 = *e22;
        valueE12 = *e12;
        valueE11 = *e11;

        xs11 = eDiffCrossSectionData[ionizationLevelIndex][valueT1][valueE11];
        xs12 = eDiffCrossSectionData[ionizationLevelIndex][valueT1][valueE12];
        xs21 = eDiffCrossSectionData[ionizationLevelIndex][valueT2][valueE21];
        xs22 = eDiffCrossSectionData[ionizationLevelIndex][valueT2][valueE22];

      }

    }

    if (particleDefinition == G4Proton::ProtonDefinition())
    {
      // k should be in eV and energy transfer eV also
      std::vector<double>::iterator t2 = std::upper_bound(pTdummyVec.begin(),
                                                          pTdummyVec.end(),
                                                          k);
      std::vector<double>::iterator t1 = t2 - 1;

      std::vector<double>::iterator e12 = std::upper_bound(pVecm[(*t1)].begin(),
                                                           pVecm[(*t1)].end(),
                                                           energyTransfer);
      std::vector<double>::iterator e11 = e12 - 1;

      std::vector<double>::iterator e22 = std::upper_bound(pVecm[(*t2)].begin(),
                                                           pVecm[(*t2)].end(),
                                                           energyTransfer);
      std::vector<double>::iterator e21 = e22 - 1;

      valueT1 = *t1;
      valueT2 = *t2;
      valueE21 = *e21;
      valueE22 = *e22;
      valueE12 = *e12;
      valueE11 = *e11;

      xs11 = pDiffCrossSectionData[ionizationLevelIndex][valueT1][valueE11];
      xs12 = pDiffCrossSectionData[ionizationLevelIndex][valueT1][valueE12];
      xs21 = pDiffCrossSectionData[ionizationLevelIndex][valueT2][valueE21];
      xs22 = pDiffCrossSectionData[ionizationLevelIndex][valueT2][valueE22];

    }

    G4double xsProduct = xs11 * xs12 * xs21 * xs22;
    if (xsProduct != 0.)
    {
      sigma = QuadInterpolator(valueE11,
                               valueE12,
                               valueE21,
                               valueE22,
                               xs11,
                               xs12,
                               xs21,
                               xs22,
                               valueT1,
                               valueT2,
                               k,
                               energyTransfer);
    }
  }

  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNABornIonisationModel1::Interpolate(G4double e1,
                                                G4double e2,
                                                G4double e,
                                                G4double xs1,
                                                G4double xs2)
{
  G4double value = 0.;

  // Log-log interpolation by default

  if (e1 != 0 && e2 != 0 && (std::log10(e2) - std::log10(e1)) != 0
      && !fasterCode)
  {
    G4double a = (std::log10(xs2) - std::log10(xs1))
        / (std::log10(e2) - std::log10(e1));
    G4double b = std::log10(xs2) - a * std::log10(e2);
    G4double sigma = a * std::log10(e) + b;
    value = (std::pow(10., sigma));
  }

  // Switch to lin-lin interpolation
  /*
   if ((e2-e1)!=0)
   {
   G4double d1 = xs1;
   G4double d2 = xs2;
   value = (d1 + (d2 - d1)*(e - e1)/ (e2 - e1));
   }
   */

  // Switch to log-lin interpolation for faster code
  if ((e2 - e1) != 0 && xs1 != 0 && xs2 != 0 && fasterCode)
  {
    G4double d1 = std::log10(xs1);
    G4double d2 = std::log10(xs2);
    value = std::pow(10., (d1 + (d2 - d1) * (e - e1) / (e2 - e1)));
  }

  // Switch to lin-lin interpolation for faster code
  // in case one of xs1 or xs2 (=cum proba) value is zero

  if ((e2 - e1) != 0 && (xs1 == 0 || xs2 == 0) && fasterCode)
  {
    G4double d1 = xs1;
    G4double d2 = xs2;
    value = (d1 + (d2 - d1) * (e - e1) / (e2 - e1));
  }

  /*
   G4cout
   << e1 << " "
   << e2 << " "
   << e  << " "
   << xs1 << " "
   << xs2 << " "
   << value
   << G4endl;
   */

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNABornIonisationModel1::QuadInterpolator(G4double e11,
                                                     G4double e12,
                                                     G4double e21,
                                                     G4double e22,
                                                     G4double xs11,
                                                     G4double xs12,
                                                     G4double xs21,
                                                     G4double xs22,
                                                     G4double t1,
                                                     G4double t2,
                                                     G4double t,
                                                     G4double e)
{
  G4double interpolatedvalue1 = Interpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = Interpolate(e21, e22, e, xs21, xs22);
  G4double value = Interpolate(t1,
                               t2,
                               t,
                               interpolatedvalue1,
                               interpolatedvalue2);

  return value;
}

G4double G4DNABornIonisationModel1::GetPartialCrossSection(const G4Material* /*material*/,
                                                           G4int level,
                                                           const G4ParticleDefinition* particle,
                                                           G4double kineticEnergy)
{
  std::map<G4String, G4DNACrossSectionDataSet*, std::less<G4String> >::iterator pos;
  pos = tableData.find(particle->GetParticleName());

  if (pos != tableData.end())
  {
    G4DNACrossSectionDataSet* table = pos->second;
    return table->GetComponent(level)->FindValue(kineticEnergy);
  }

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNABornIonisationModel1::RandomSelect(G4double k,
                                              const G4String& particle)
{
  G4int level = 0;

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
    G4Exception("G4DNABornIonisationModel1::RandomSelect",
                "em0002",
                FatalException,
                "Model not applicable to particle type.");
  }

  return level;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNABornIonisationModel1::RandomizeEjectedElectronEnergyFromCumulatedDcs(G4ParticleDefinition* particleDefinition,
                                                                                   G4double k,
                                                                                   G4int shell)
{
  //G4cout << "*** FAST computation for " << " " << particleDefinition->GetParticleName() << G4endl;

  G4double secondaryElectronKineticEnergy = 0.;

  G4double random = G4UniformRand();

  secondaryElectronKineticEnergy = TransferedEnergy(particleDefinition,
                                                    k / eV,
                                                    shell,
                                                    random) * eV
      - waterStructure.IonisationEnergy(shell);

  //G4cout << RandomTransferedEnergy(particleDefinition, k/eV, shell) << G4endl;
  // SI - 01/04/2014
  if (secondaryElectronKineticEnergy < 0.)
    return 0.;
  //

  return secondaryElectronKineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNABornIonisationModel1::TransferedEnergy(G4ParticleDefinition* particleDefinition,
                                                     G4double k,
                                                     G4int ionizationLevelIndex,
                                                     G4double random)
{
  G4double nrj = 0.;

  G4double valueK1 = 0;
  G4double valueK2 = 0;
  G4double valuePROB21 = 0;
  G4double valuePROB22 = 0;
  G4double valuePROB12 = 0;
  G4double valuePROB11 = 0;

  G4double nrjTransf11 = 0;
  G4double nrjTransf12 = 0;
  G4double nrjTransf21 = 0;
  G4double nrjTransf22 = 0;

  if (particleDefinition == G4Electron::ElectronDefinition())
  {
    // k should be in eV
    std::vector<double>::iterator k2 = std::upper_bound(eTdummyVec.begin(),
                                                        eTdummyVec.end(),
                                                        k);
    std::vector<double>::iterator k1 = k2 - 1;

    /*
     G4cout << "----> k=" << k
     << " " << *k1
     << " " << *k2
     << " " << random
     << " " << ionizationLevelIndex
     << " " << eProbaShellMap[ionizationLevelIndex][(*k1)].back()
     << " " << eProbaShellMap[ionizationLevelIndex][(*k2)].back()
     << G4endl;
     */

    // SI : the following condition avoids situations where random >last vector element
    if (random <= eProbaShellMap[ionizationLevelIndex][(*k1)].back()
        && random <= eProbaShellMap[ionizationLevelIndex][(*k2)].back())
    {
      std::vector<double>::iterator prob12 =
          std::upper_bound(eProbaShellMap[ionizationLevelIndex][(*k1)].begin(),
                           eProbaShellMap[ionizationLevelIndex][(*k1)].end(),
                           random);

      std::vector<double>::iterator prob11 = prob12 - 1;

      std::vector<double>::iterator prob22 =
          std::upper_bound(eProbaShellMap[ionizationLevelIndex][(*k2)].begin(),
                           eProbaShellMap[ionizationLevelIndex][(*k2)].end(),
                           random);

      std::vector<double>::iterator prob21 = prob22 - 1;

      valueK1 = *k1;
      valueK2 = *k2;
      valuePROB21 = *prob21;
      valuePROB22 = *prob22;
      valuePROB12 = *prob12;
      valuePROB11 = *prob11;

      /*
       G4cout << "        " << random << " " << valuePROB11 << " "
       << valuePROB12 << " " << valuePROB21 << " " << valuePROB22 << G4endl;
       */

      nrjTransf11 = eNrjTransfData[ionizationLevelIndex][valueK1][valuePROB11];
      nrjTransf12 = eNrjTransfData[ionizationLevelIndex][valueK1][valuePROB12];
      nrjTransf21 = eNrjTransfData[ionizationLevelIndex][valueK2][valuePROB21];
      nrjTransf22 = eNrjTransfData[ionizationLevelIndex][valueK2][valuePROB22];

      /*
       G4cout << "        " << ionizationLevelIndex << " "
       << random << " " <<valueK1 << " " << valueK2 << G4endl;

       G4cout << "        " << random << " " << nrjTransf11 << " "
       << nrjTransf12 << " " << nrjTransf21 << " " <<nrjTransf22 << G4endl;
       */

    }
    // Avoids cases where cum xs is zero for k1 and is not for k2 (with always k1<k2)
    if (random > eProbaShellMap[ionizationLevelIndex][(*k1)].back())
    {
      std::vector<double>::iterator prob22 =
          std::upper_bound(eProbaShellMap[ionizationLevelIndex][(*k2)].begin(),
                           eProbaShellMap[ionizationLevelIndex][(*k2)].end(),
                           random);

      std::vector<double>::iterator prob21 = prob22 - 1;

      valueK1 = *k1;
      valueK2 = *k2;
      valuePROB21 = *prob21;
      valuePROB22 = *prob22;

      //G4cout << "        " << random << " " << valuePROB21 << " " << valuePROB22 << G4endl;

      nrjTransf21 = eNrjTransfData[ionizationLevelIndex][valueK2][valuePROB21];
      nrjTransf22 = eNrjTransfData[ionizationLevelIndex][valueK2][valuePROB22];

      G4double interpolatedvalue2 = Interpolate(valuePROB21,
                                                valuePROB22,
                                                random,
                                                nrjTransf21,
                                                nrjTransf22);

      // zeros are explicitely set

      G4double value = Interpolate(valueK1, valueK2, k, 0., interpolatedvalue2);

      /*
       G4cout << "        " << ionizationLevelIndex << " "
       << random << " " <<valueK1 << " " << valueK2 << G4endl;

       G4cout << "        " << random << " " << nrjTransf11 << " "
       << nrjTransf12 << " " << nrjTransf21 << " " <<nrjTransf22 << G4endl;

       G4cout << "ici" << " " << value << G4endl;
       */

      return value;
    }
  }
  //
  else if (particleDefinition == G4Proton::ProtonDefinition())
  {
    // k should be in eV

    std::vector<double>::iterator k2 = std::upper_bound(pTdummyVec.begin(),
                                                        pTdummyVec.end(),
                                                        k);

    std::vector<double>::iterator k1 = k2 - 1;

    /*
     G4cout << "----> k=" << k
     << " " << *k1
     << " " << *k2
     << " " << random
     << " " << ionizationLevelIndex
     << " " << pProbaShellMap[ionizationLevelIndex][(*k1)].back()
     << " " << pProbaShellMap[ionizationLevelIndex][(*k2)].back()
     << G4endl;
     */

    // SI : the following condition avoids situations where random > last vector element,
    //      for eg. when the last element is zero
    if (random <= pProbaShellMap[ionizationLevelIndex][(*k1)].back()
        && random <= pProbaShellMap[ionizationLevelIndex][(*k2)].back())
    {
      std::vector<double>::iterator prob12 =
          std::upper_bound(pProbaShellMap[ionizationLevelIndex][(*k1)].begin(),
                           pProbaShellMap[ionizationLevelIndex][(*k1)].end(),
                           random);

      std::vector<double>::iterator prob11 = prob12 - 1;

      std::vector<double>::iterator prob22 =
          std::upper_bound(pProbaShellMap[ionizationLevelIndex][(*k2)].begin(),
                           pProbaShellMap[ionizationLevelIndex][(*k2)].end(),
                           random);

      std::vector<double>::iterator prob21 = prob22 - 1;

      valueK1 = *k1;
      valueK2 = *k2;
      valuePROB21 = *prob21;
      valuePROB22 = *prob22;
      valuePROB12 = *prob12;
      valuePROB11 = *prob11;

      /*
       G4cout << "        " << random << " " << valuePROB11 << " "
       << valuePROB12 << " " << valuePROB21 << " " << valuePROB22 << G4endl;
       */

      nrjTransf11 = pNrjTransfData[ionizationLevelIndex][valueK1][valuePROB11];
      nrjTransf12 = pNrjTransfData[ionizationLevelIndex][valueK1][valuePROB12];
      nrjTransf21 = pNrjTransfData[ionizationLevelIndex][valueK2][valuePROB21];
      nrjTransf22 = pNrjTransfData[ionizationLevelIndex][valueK2][valuePROB22];

      /*
       G4cout << "        " << ionizationLevelIndex << " "
       << random << " " <<valueK1 << " " << valueK2 << G4endl;

       G4cout << "        " << random << " " << nrjTransf11 << " "
       << nrjTransf12 << " " << nrjTransf21 << " " <<nrjTransf22 << G4endl;
       */
    }

    // Avoids cases where cum xs is zero for k1 and is not for k2 (with always k1<k2)

    if (random > pProbaShellMap[ionizationLevelIndex][(*k1)].back())
    {
      std::vector<double>::iterator prob22 =
          std::upper_bound(pProbaShellMap[ionizationLevelIndex][(*k2)].begin(),
                           pProbaShellMap[ionizationLevelIndex][(*k2)].end(),
                           random);

      std::vector<double>::iterator prob21 = prob22 - 1;

      valueK1 = *k1;
      valueK2 = *k2;
      valuePROB21 = *prob21;
      valuePROB22 = *prob22;

      //G4cout << "        " << random << " " << valuePROB21 << " " << valuePROB22 << G4endl;

      nrjTransf21 = pNrjTransfData[ionizationLevelIndex][valueK2][valuePROB21];
      nrjTransf22 = pNrjTransfData[ionizationLevelIndex][valueK2][valuePROB22];

      G4double interpolatedvalue2 = Interpolate(valuePROB21,
                                                valuePROB22,
                                                random,
                                                nrjTransf21,
                                                nrjTransf22);

      // zeros are explicitely set

      G4double value = Interpolate(valueK1, valueK2, k, 0., interpolatedvalue2);

      /*
       G4cout << "        " << ionizationLevelIndex << " "
       << random << " " <<valueK1 << " " << valueK2 << G4endl;

       G4cout << "        " << random << " " << nrjTransf11 << " "
       << nrjTransf12 << " " << nrjTransf21 << " " <<nrjTransf22 << G4endl;

       G4cout << "ici" << " " << value << G4endl;
       */

      return value;
    }
  }
  // End electron and proton cases

  G4double nrjTransfProduct = nrjTransf11 * nrjTransf12 * nrjTransf21
      * nrjTransf22;
  //G4cout << "nrjTransfProduct=" << nrjTransfProduct << G4endl;

  if (nrjTransfProduct != 0.)
  {
    nrj = QuadInterpolator(valuePROB11,
                           valuePROB12,
                           valuePROB21,
                           valuePROB22,
                           nrjTransf11,
                           nrjTransf12,
                           nrjTransf21,
                           nrjTransf22,
                           valueK1,
                           valueK2,
                           k,
                           random);
  }
  //G4cout << nrj << endl;

  return nrj;
}
