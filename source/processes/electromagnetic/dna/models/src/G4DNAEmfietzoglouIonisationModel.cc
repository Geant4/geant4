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
// Based on the work described in
// Rad Res 163, 98-111 (2005)
// D. Emfietzoglou, H. Nikjoo
//
// Authors of the class (2014):
// I. Kyriakou (kyriak@cc.uoi.gr)
// D. Emfietzoglou (demfietz@cc.uoi.gr)
// S. Incerti (incerti@cenbg.in2p3.fr)
//

#include "G4DNAEmfietzoglouIonisationModel.hh"
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

G4DNAEmfietzoglouIonisationModel::G4DNAEmfietzoglouIonisationModel(const G4ParticleDefinition*,
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

  if(verboseLevel > 0)
  {
    G4cout << "Emfietzoglou ionisation model is constructed " << G4endl;
  }

  // Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);
  fAtomDeexcitation = 0;
  fParticleChangeForGamma = 0;
  fpMolWaterDensity = 0;

  // Define default angular generator
  SetAngularDistribution(new G4DNABornAngle());

  SetLowEnergyLimit(10. * eV);
  SetHighEnergyLimit(10. * keV);

  // Selection of computation method

  fasterCode = false;

  // Selection of stationary mode

  statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAEmfietzoglouIonisationModel::~G4DNAEmfietzoglouIonisationModel()
{
  // Cross section

  std::map<G4String, G4DNACrossSectionDataSet*, std::less<G4String> >::iterator pos;
  for(pos = tableData.begin(); pos != tableData.end(); ++pos)
  {
    G4DNACrossSectionDataSet* table = pos->second;
    delete table;
  }

  // Final state

  eVecm.clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAEmfietzoglouIonisationModel::Initialise(const G4ParticleDefinition* particle,
                                                  const G4DataVector& /*cuts*/)
{

  if(verboseLevel > 3)
  {
    G4cout << "Calling G4DNAEmfietzoglouIonisationModel::Initialise()" << G4endl;
  }

  // Energy limits

  G4String fileElectron("dna/sigma_ionisation_e_emfietzoglou");

  G4ParticleDefinition* electronDef = G4Electron::ElectronDefinition();

  G4String electron;

  G4double scaleFactor = (1.e-22 / 3.343) * m*m;

  const char *path = G4FindDataDir("G4LEDATA");

  // *** ELECTRON

  electron = electronDef->GetParticleName();

  tableFile[electron] = fileElectron;

  // Cross section

  G4DNACrossSectionDataSet* tableE = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
  tableE->LoadData(fileElectron);

  tableData[electron] = tableE;

  // Final state

  std::ostringstream eFullFileName;

  if (fasterCode) eFullFileName << path << "/dna/sigmadiff_cumulated_ionisation_e_emfietzoglou.dat";
  if (!fasterCode) eFullFileName << path << "/dna/sigmadiff_ionisation_e_emfietzoglou.dat";

  std::ifstream eDiffCrossSection(eFullFileName.str().c_str());

  if (!eDiffCrossSection)
  {
    if (fasterCode) G4Exception("G4DNAEmfietzoglouIonisationModel::Initialise","em0003",
        FatalException,"Missing data file:/dna/sigmadiff_cumulated_ionisation_e_emfietzoglou.dat");

    if (!fasterCode) G4Exception("G4DNAEmfietzoglouIonisationModel::Initialise","em0003",
        FatalException,"Missing data file:/dna/sigmadiff_ionisation_e_emfietzoglou.dat");
  }

  //

  // Clear the arrays for re-initialization case (MT mode)
  // March 25th, 2014 - Vaclav Stepan, Sebastien Incerti

  eTdummyVec.clear();

  eVecm.clear();

  eProbaShellMap->clear();

  eDiffCrossSectionData->clear();

  eNrjTransfData->clear();

  //

  eTdummyVec.push_back(0.);
  while(!eDiffCrossSection.eof())
  {
    G4double tDummy;
    G4double eDummy;
    eDiffCrossSection>>tDummy>>eDummy;
    if (tDummy != eTdummyVec.back()) eTdummyVec.push_back(tDummy);
    for (G4int j=0; j<5; j++)
    {
      eDiffCrossSection>>eDiffCrossSectionData[j][tDummy][eDummy];

      if (fasterCode)
      {
        eNrjTransfData[j][tDummy][eDiffCrossSectionData[j][tDummy][eDummy]]=eDummy;
        eProbaShellMap[j][tDummy].push_back(eDiffCrossSectionData[j][tDummy][eDummy]);
      }

      // SI - only if eof is not reached
      if (!eDiffCrossSection.eof() && !fasterCode)
      {
        eDiffCrossSectionData[j][tDummy][eDummy]*=scaleFactor;
      }

      if (!fasterCode) eVecm[tDummy].push_back(eDummy);

    }
  }

  //

  if( verboseLevel>0 )
  {
    G4cout << "Emfietzoglou ionisation model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / eV << " eV - "
    << HighEnergyLimit() / keV << " keV for "
    << particle->GetParticleName()
    << G4endl;
  }

  // Initialize water density pointer

  fpMolWaterDensity =
      G4DNAMolecularMaterial::Instance()->
        GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

  // AD

  fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();

  if (isInitialised)
  { return;}
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAEmfietzoglouIonisationModel::
CrossSectionPerVolume(const G4Material* material,
                      const G4ParticleDefinition* particleDefinition,
                      G4double ekin,
                      G4double,
                      G4double)
{
  if(verboseLevel > 3)
  {
    G4cout
        << "Calling CrossSectionPerVolume() of G4DNAEmfietzoglouIonisationModel"
        << G4endl;
  }

  if (particleDefinition != G4Electron::ElectronDefinition()) return 0; // necessary ??

  // Calculate total cross section for model

  G4double sigma=0;

  G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];

  const G4String& particleName = particleDefinition->GetParticleName();

  if (ekin >= LowEnergyLimit() && ekin <= HighEnergyLimit())
  {
    std::map< G4String,G4DNACrossSectionDataSet*,std::less<G4String> >::iterator pos;
    pos = tableData.find(particleName);

    if (pos != tableData.end())
    {
      G4DNACrossSectionDataSet* table = pos->second;
      if (table != 0)
      {
        sigma = table->FindValue(ekin);
      }
    }
    else
    {
      G4Exception("G4DNAEmfietzoglouIonisationModel::CrossSectionPerVolume","em0002",
          FatalException,"Model not applicable to particle type.");
    }
  }

  if (verboseLevel > 2)
  {
    G4cout << "__________________________________" << G4endl;
    G4cout << "G4DNAEmfietzoglouIonisationModel - XS INFO START" << G4endl;
    G4cout << "Kinetic energy(eV)=" << ekin/eV << " particle : " << particleName << G4endl;
    G4cout << "Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
    G4cout << "Cross section per water molecule (cm^-1)=" << sigma*waterDensity/(1./cm) << G4endl;
    G4cout << "G4DNAEmfietzoglouIonisationModel - XS INFO END" << G4endl;
  }

  return sigma*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAEmfietzoglouIonisationModel::
SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                  const G4MaterialCutsCouple* couple,
                  const G4DynamicParticle* particle,
                  G4double,
                  G4double)
{

  if(verboseLevel > 3)
  {
    G4cout << "Calling SampleSecondaries() of G4DNAEmfietzoglouIonisationModel"
           << G4endl;
  }

  G4double k = particle->GetKineticEnergy();

  const G4String& particleName = particle->GetDefinition()->GetParticleName();

  if (k >= LowEnergyLimit() && k <= HighEnergyLimit())
  {
    G4ParticleMomentum primaryDirection = particle->GetMomentumDirection();
    G4double particleMass = particle->GetDefinition()->GetPDGMass();
    G4double totalEnergy = k + particleMass;
    G4double pSquare = k * (totalEnergy + particleMass);
    G4double totalMomentum = std::sqrt(pSquare);

    G4int ionizationShell = 0;

    ionizationShell = RandomSelect(k,particleName);

    G4double bindingEnergy = 0;
    bindingEnergy = waterStructure.IonisationEnergy(ionizationShell);

    // SI : additional protection if tcs interpolation method is modified
    if (k<bindingEnergy) return;
    //

    G4double secondaryKinetic=-1000*eV;

    if (!fasterCode) secondaryKinetic = RandomizeEjectedElectronEnergy(particle->GetDefinition(),k,ionizationShell);

    if (fasterCode)
    secondaryKinetic = RandomizeEjectedElectronEnergyFromCumulatedDcs(particle->GetDefinition(),k,ionizationShell);

    // SI - For atom. deexc. tagging - 23/05/2017

    G4int Z = 8;

    G4ThreeVector deltaDirection =
    GetAngularDistribution()->SampleDirectionForShell(particle, secondaryKinetic,
        Z, ionizationShell,
        couple->GetMaterial());

    if (secondaryKinetic>0)
    {
      G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),deltaDirection,secondaryKinetic);
      fvect->push_back(dp);
    }

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

    // AM: sample deexcitation
    // here we assume that H_{2}O electronic levels are the same as Oxygen.
    // this can be considered true with a rough 10% error in energy on K-shell,

    size_t secNumberInit = 0;// need to know at a certain point the energy of secondaries
    size_t secNumberFinal = 0;// So I'll make the diference and then sum the energies

    G4double scatteredEnergy = k-bindingEnergy-secondaryKinetic;

    // SI: only atomic deexcitation from K shell is considered
    if(fAtomDeexcitation && ionizationShell == 4)
    {
      const G4AtomicShell* shell
        = fAtomDeexcitation->GetAtomicShell(Z, G4AtomicShellEnumerator(0));
      secNumberInit = fvect->size();
      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, 0, 0);
      secNumberFinal = fvect->size();

      if(secNumberFinal > secNumberInit) {
	for (size_t i=secNumberInit; i<secNumberFinal; ++i) {
          //Check if there is enough residual energy
          if (bindingEnergy >= ((*fvect)[i])->GetKineticEnergy())
           {
             //Ok, this is a valid secondary: keep it
	     bindingEnergy -= ((*fvect)[i])->GetKineticEnergy();
           }
          else
           {
 	     //Invalid secondary: not enough energy to create it!
 	     //Keep its energy in the local deposit
             delete (*fvect)[i];
             (*fvect)[i]=0;
           }
	}
      }

    }

    //This should never happen
    if(bindingEnergy < 0.0)
     G4Exception("G4DNAEmfietzoglouIonisatioModel1::SampleSecondaries()",
                 "em2050",FatalException,"Negative local energy deposit");

    //bindingEnergy has been decreased
    //by the amount of energy taken away by deexc. products
    if (!statCode)
    {
      fParticleChangeForGamma->SetProposedKineticEnergy(scatteredEnergy);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(bindingEnergy);
    }
    else
    {
      fParticleChangeForGamma->SetProposedKineticEnergy(k);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(k-scatteredEnergy);
    }
    // TEST //////////////////////////
    // if (secondaryKinetic<0) abort();
    // if (scatteredEnergy<0) abort();
    // if (k-scatteredEnergy-secondaryKinetic-deexSecEnergy<0) abort();
    // if (k-scatteredEnergy<0) abort();
    /////////////////////////////////

    const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
    G4DNAChemistryManager::Instance()->CreateWaterMolecule(eIonizedMolecule,
        ionizationShell,
        theIncomingTrack);
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4DNAEmfietzoglouIonisationModel::
RandomizeEjectedElectronEnergy(G4ParticleDefinition* particleDefinition,
                               G4double k,
                               G4int shell)
{
  // G4cout << "*** SLOW computation for "
  //        << " " << particleDefinition->GetParticleName() << G4endl;

  if(particleDefinition == G4Electron::ElectronDefinition())
  {
    G4double maximumEnergyTransfer = 0.;
    if((k + waterStructure.IonisationEnergy(shell)) / 2. > k)
      maximumEnergyTransfer =  k;
    else
      maximumEnergyTransfer = (k + waterStructure.IonisationEnergy(shell))/ 2.;

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
    while(step > 0)
    {
      step--;
      G4double differentialCrossSection =
          DifferentialCrossSection(particleDefinition,
                                   k / eV,
                                   value / eV,
                                   shell);
      if(differentialCrossSection >= crossSectionMaximum) crossSectionMaximum =
          differentialCrossSection;
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

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// The following section is not used anymore but is kept for memory
// GetAngularDistribution()->SampleDirectionForShell is used instead

/*
 void G4DNAEmfietzoglouIonisationModel::RandomizeEjectedElectronDirection(G4ParticleDefinition* particleDefinition,
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
G4double G4DNAEmfietzoglouIonisationModel::DifferentialCrossSection(G4ParticleDefinition * particleDefinition,
                                                                  G4double k,
                                                                  G4double energyTransfer,
                                                                  G4int ionizationLevelIndex)
{
  G4double sigma = 0.;

  if(energyTransfer >= waterStructure.IonisationEnergy(ionizationLevelIndex)/eV)
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

    if(particleDefinition == G4Electron::ElectronDefinition())
    {
      // Protection against out of boundary access
      if (k==eTdummyVec.back()) k=k*(1.-1e-12);
      //

      // k should be in eV and energy transfer eV also

      std::vector<G4double>::iterator t2 = std::upper_bound(eTdummyVec.begin(),
                                                          eTdummyVec.end(),
                                                          k);

      std::vector<G4double>::iterator t1 = t2 - 1;

      // SI : the following condition avoids situations where energyTransfer >last vector element
      // added strict limitations (09/08/2017)
      if(energyTransfer < eVecm[(*t1)].back() &&
         energyTransfer < eVecm[(*t2)].back())
      {
        std::vector<G4double>::iterator e12 =
            std::upper_bound(eVecm[(*t1)].begin(),
                             eVecm[(*t1)].end(),
                             energyTransfer);
        std::vector<G4double>::iterator e11 = e12 - 1;

        std::vector<G4double>::iterator e22 =
            std::upper_bound(eVecm[(*t2)].begin(),
                             eVecm[(*t2)].end(),
                             energyTransfer);
        std::vector<G4double>::iterator e21 = e22 - 1;

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

        //G4cout << "-------------------" << G4endl;
        //G4cout << "ionizationLevelIndex=" << ionizationLevelIndex << G4endl;
        //G4cout << "valueT1/eV=" << valueT1 << " valueT2/eV=" << valueT2 << G4endl;
        //G4cout << "valueE11/eV=" << valueE11 << " valueE12/eV=" << valueE12
        //       << " valueE21/eV=" << valueE21 << " valueE22/eV=" << valueE22 << G4endl;
        //G4cout << "xs11=" << xs11 / ((1.e-22 / 3.343) * m*m) << G4endl;
        //G4cout << "xs12=" << xs12 / ((1.e-22 / 3.343) * m*m) << G4endl;
        //G4cout << "xs21=" << xs21 / ((1.e-22 / 3.343) * m*m) << G4endl;
        //G4cout << "xs22=" << xs22 / ((1.e-22 / 3.343) * m*m) << G4endl;
        //G4cout << "###################" << G4endl;

      }

    }

    G4double xsProduct = xs11 * xs12 * xs21 * xs22;
    if(xsProduct != 0.)
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

G4double G4DNAEmfietzoglouIonisationModel::Interpolate(G4double e1,
                                                       G4double e2,
                                                       G4double e,
                                                       G4double xs1,
                                                       G4double xs2)
{

  G4double value = 0.;

  // Log-log interpolation by default

  if(e1 != 0 && e2 != 0 && (std::log10(e2) - std::log10(e1)) != 0
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
  if((e2 - e1) != 0 && xs1 != 0 && xs2 != 0 && fasterCode)
  {
    G4double d1 = std::log10(xs1);
    G4double d2 = std::log10(xs2);
    value = std::pow(10., (d1 + (d2 - d1) * (e - e1) / (e2 - e1)));
  }

  // Switch to lin-lin interpolation for faster code
  // in case one of xs1 or xs2 (=cum proba) value is zero

  if((e2 - e1) != 0 && (xs1 == 0 || xs2 == 0) && fasterCode)
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

G4double G4DNAEmfietzoglouIonisationModel::QuadInterpolator(G4double e11,
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNAEmfietzoglouIonisationModel::RandomSelect(G4double k,
                                                     const G4String& particle)
{
  G4int level = 0;

  auto pos = tableData.find(particle);

  if(pos != tableData.cend())
  {
    G4DNACrossSectionDataSet* table = pos->second;

    if(table != 0)
    {
      G4double* valuesBuffer = new G4double[table->NumberOfComponents()];
      const G4int n = (G4int)table->NumberOfComponents();
      G4int i(n);
      G4double value = 0.;

      while(i > 0)
      {
        i--;
        valuesBuffer[i] = table->GetComponent(i)->FindValue(k);
        value += valuesBuffer[i];
      }

      value *= G4UniformRand();

      i = n;

      while(i > 0)
      {
        i--;

        if(valuesBuffer[i] > value)
        {
          delete[] valuesBuffer;
          return i;
        }
        value -= valuesBuffer[i];
      }

      if(valuesBuffer) delete[] valuesBuffer;

    }
  }
  else
  {
    G4Exception("G4DNAEmfietzoglouIonisationModel::RandomSelect",
                "em0002",
                FatalException,
                "Model not applicable to particle type.");
  }

  return level;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNAEmfietzoglouIonisationModel::RandomizeEjectedElectronEnergyFromCumulatedDcs(G4ParticleDefinition* particleDefinition,
                                                                                          G4double k,
                                                                                          G4int shell)
{
  //G4cout << "*** FAST computation for " << " " << particleDefinition->GetParticleName() << G4endl;

  G4double secondaryElectronKineticEnergy = 0.;

  secondaryElectronKineticEnergy = RandomTransferedEnergy(particleDefinition,
                                                          k / eV,
                                                          shell)
                                   * eV
                                   - waterStructure.IonisationEnergy(shell);

  //G4cout << RandomTransferedEnergy(particleDefinition, k/eV, shell) << G4endl;
  if(secondaryElectronKineticEnergy < 0.) return 0.;

  return secondaryElectronKineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNAEmfietzoglouIonisationModel::RandomTransferedEnergy(G4ParticleDefinition* particleDefinition,
                                                                  G4double k,
                                                                  G4int ionizationLevelIndex)
{

  G4double random = G4UniformRand();

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
    // Protection against out of boundary access
    if (k==eTdummyVec.back()) k=k*(1.-1e-12);
    //

    // k should be in eV
    std::vector<G4double>::iterator k2 = std::upper_bound(eTdummyVec.begin(),eTdummyVec.end(), k);

    std::vector<G4double>::iterator k1 = k2-1;

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
    if ( random <= eProbaShellMap[ionizationLevelIndex][(*k1)].back()
        && random <= eProbaShellMap[ionizationLevelIndex][(*k2)].back() )

    {
      std::vector<G4double>::iterator prob12 = std::upper_bound(eProbaShellMap[ionizationLevelIndex][(*k1)].begin(),
          eProbaShellMap[ionizationLevelIndex][(*k1)].end(), random);

      std::vector<G4double>::iterator prob11 = prob12-1;

      std::vector<G4double>::iterator prob22 = std::upper_bound(eProbaShellMap[ionizationLevelIndex][(*k2)].begin(),
          eProbaShellMap[ionizationLevelIndex][(*k2)].end(), random);

      std::vector<G4double>::iterator prob21 = prob22-1;

      valueK1 =*k1;
      valueK2 =*k2;
      valuePROB21 =*prob21;
      valuePROB22 =*prob22;
      valuePROB12 =*prob12;
      valuePROB11 =*prob11;

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

    if ( random > eProbaShellMap[ionizationLevelIndex][(*k1)].back() )

    {
      std::vector<G4double>::iterator prob22 = std::upper_bound(eProbaShellMap[ionizationLevelIndex][(*k2)].begin(),
          eProbaShellMap[ionizationLevelIndex][(*k2)].end(), random);

      std::vector<G4double>::iterator prob21 = prob22-1;

      valueK1 =*k1;
      valueK2 =*k2;
      valuePROB21 =*prob21;
      valuePROB22 =*prob22;

      //G4cout << "        " << random << " " << valuePROB21 << " " << valuePROB22 << G4endl;

      nrjTransf21 = eNrjTransfData[ionizationLevelIndex][valueK2][valuePROB21];
      nrjTransf22 = eNrjTransfData[ionizationLevelIndex][valueK2][valuePROB22];

      G4double interpolatedvalue2 = Interpolate(valuePROB21, valuePROB22, random, nrjTransf21, nrjTransf22);

      // zeros are explicitly set

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

  // End electron

  G4double nrjTransfProduct = nrjTransf11 * nrjTransf12 * nrjTransf21 * nrjTransf22;

  //G4cout << "nrjTransfProduct=" << nrjTransfProduct << G4endl;

  if (nrjTransfProduct != 0.)
  {
    nrj = QuadInterpolator( valuePROB11, valuePROB12,
        valuePROB21, valuePROB22,
        nrjTransf11, nrjTransf12,
        nrjTransf21, nrjTransf22,
        valueK1, valueK2,
        k, random);
  }

  //G4cout << nrj << endl;

  return nrj;
}
