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
// Reference:
//    A.D. Dominguez-Munoz, M.I. Gallardo, M.C. Bordage,
//    Z. Francis, S. Incerti, M.A. Cortes-Giraldo,
//    Radiat. Phys. Chem. 199 (2022) 110363.
//
// Class authors:
//    A.D. Dominguez-Munoz
//    M.A. Cortes-Giraldo (miancortes -at- us.es)
//
// Class creation: 2022-03-03
//
//

#include "G4DNARPWBAIonisationModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4DNABornAngle.hh"
#include "G4Exp.hh"
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4DNARPWBAIonisationModel::G4DNARPWBAIonisationModel(
  const G4ParticleDefinition*, const G4String& nam)
  : G4VEmModel(nam)
{
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if(verboseLevel > 0)
  {
    G4cout << "RPWBA ionisation model is constructed " << G4endl;
  }
  SetDeexcitationFlag(true);
  SetAngularDistribution(new G4DNABornAngle());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNARPWBAIonisationModel::~G4DNARPWBAIonisationModel()
{
  eVecm.clear();
  pVecm.clear();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4DNARPWBAIonisationModel::InEnergyLimit(const G4double& k)
{
  if(lowEnergyLimit == highEnergyLimit)
  {
    G4Exception("G4DNARPWBAIonisationModel::InEnergyLimit", "em0102",
                FatalException, "lowEnergyLimit == highEnergyLimit");
  }
  if(k >= lowEnergyLimit && k <= highEnergyLimit)
  {
    return true;
  }
  else
  {
    return false;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4DNARPWBAIonisationModel::InitialiseForProton(
  const G4ParticleDefinition* part)
{
  if(part != fProtonDef)
  {
    G4Exception("G4DNARPWBAIonisationModel::InitialiseForProton", "em0002",
                FatalException, "Model not applicable to particle type.");
  }
  // Energy limits
  G4String fileProton("dna/sigma_ionisation_p_RPWBA");
  G4double scaleFactor = 1 * cm * cm;
  const char *path    = G4FindDataDir("G4LEDATA");
  lowEnergyLimit      = 100. * MeV;
  highEnergyLimit     = 300. * MeV;

  if(LowEnergyLimit() < lowEnergyLimit || HighEnergyLimit() > highEnergyLimit)
  {
    G4ExceptionDescription ed;
    ed << "Model is applicable from "<<lowEnergyLimit<<" to "<<highEnergyLimit;
    G4Exception("G4DNARPWBAIonisationModel::InitialiseForProton", "em0004",
      FatalException, ed);
  }

  fpTotalCrossSection = make_unique<G4DNACrossSectionDataSet>(
    new G4LogLogInterpolation, eV, scaleFactor);
  fpTotalCrossSection->LoadData(fileProton);

  // Final state

  std::ostringstream pFullFileName;
  fasterCode ? pFullFileName
                 << path << "/dna/sigmadiff_cumulated_ionisation_p_RPWBA.dat"
             : pFullFileName << path << "/dna/sigmadiff_ionisation_p_RPWBA.dat";
  std::ifstream pDiffCrossSection(pFullFileName.str().c_str());
  if(!pDiffCrossSection)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "Missing data file: " + pFullFileName.str();
    G4Exception("G4DNARPWBAIonisationModel::InitialiseForProton", "em0003",
                FatalException, exceptionDescription);
  }

  pTdummyVec.push_back(0.);
  while(!pDiffCrossSection.eof())
  {
    G4double tDummy;
    G4double eDummy;
    pDiffCrossSection >> tDummy >> eDummy;
    if(tDummy != pTdummyVec.back())
    {
      pTdummyVec.push_back(tDummy);
    }

    for(G4int j = 0; j < 5; j++)
    {
      pDiffCrossSection >> pDiffCrossSectionData[j][tDummy][eDummy];

      if(fasterCode)
      {
        pNrjTransfData[j][tDummy][pDiffCrossSectionData[j][tDummy][eDummy]] =
          eDummy;
        pProbaShellMap[j][tDummy].push_back(
          pDiffCrossSectionData[j][tDummy][eDummy]);
      }

      // SI - only if eof is not reached !
      if(!pDiffCrossSection.eof() && !fasterCode)
      {
        pDiffCrossSectionData[j][tDummy][eDummy] *= scaleFactor;
      }

      if(!fasterCode)
      {
        pVecm[tDummy].push_back(eDummy);
      }
    }
  }

  // be careful about this
  // SetLowEnergyLimit(lowEnergyLimit);
  // SetHighEnergyLimit(highEnergyLimit);
}

void G4DNARPWBAIonisationModel::Initialise(const G4ParticleDefinition* particle,
                                           const G4DataVector& /*cuts*/)
{
  if(isInitialised)
  {
    return;
  }
  if(verboseLevel > 3)
  {
    G4cout << "Calling G4DNARPWBAIonisationModel::Initialise()"
           << particle->GetParticleName() << G4endl;
  }

  InitialiseForProton(particle);

  if(verboseLevel > 0)
  {
    G4cout << "RPWBA ionisation model is initialized " << G4endl
           << "Energy range: " << LowEnergyLimit() / MeV << " MeV - "
           << HighEnergyLimit() / MeV << " MeV for "
           << particle->GetParticleName() << G4endl;
  }

  // Initialize water density pointer
  if(G4Material::GetMaterial("G4_WATER") != nullptr)
  {
    fpMolWaterDensity =
      G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(
        G4Material::GetMaterial("G4_WATER"));
  }
  else
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription << "G4_WATER does not exist :";
    G4Exception("G4DNARPWBAIonisationModel::Initialise", "em00020",
                FatalException, exceptionDescription);
  }
  fAtomDeexcitation       = G4LossTableManager::Instance()->AtomDeexcitation();
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised           = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNARPWBAIonisationModel::CrossSectionPerVolume(
  const G4Material* material, const G4ParticleDefinition* particleDefinition,
  G4double ekin, G4double, G4double)
{
  if(particleDefinition != fProtonDef)
  {
    G4Exception("G4DNARPWBAIonisationModel::CrossSectionPerVolume", "em0402",
                FatalException, "Model not applicable to particle type.");
  }
  if(verboseLevel > 3)
  {
    G4cout << "Calling CrossSectionPerVolume() of G4DNARPWBAIonisationModel"
           << G4endl;
  }
  G4double sigma;
  G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];

  if(InEnergyLimit(ekin))
  {
    sigma = fpTotalCrossSection->FindValue(ekin);
  }
  else
  {
    // nput energy is outside this interval the cross section is set to zero
    // should add a warning or exception ?
    return 0;
  }

  if(verboseLevel > 2)
  {
    G4cout << "__________________________________" << G4endl;
    G4cout << "G4DNARPWBAIonisationModel - XS INFO START" << G4endl;
    G4cout << "Kinetic energy(eV)=" << ekin / eV
           << " particle : " << fProtonDef->GetParticleName() << G4endl;
    G4cout << "Cross section per water molecule (cm^2)=" << sigma / cm / cm
           << G4endl;
    G4cout << "Cross section per water molecule (cm^-1)="
           << sigma * waterDensity / (1. / cm) << G4endl;
    G4cout << "G4DNARPWBAIonisationModel - XS INFO END" << G4endl;
  }

  return sigma * waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNARPWBAIonisationModel::SampleSecondaries(
  std::vector<G4DynamicParticle*>* fvect, const G4MaterialCutsCouple* couple,
  const G4DynamicParticle* particle, G4double, G4double)
{
  if(verboseLevel > 3)
  {
    G4cout << "Calling SampleSecondaries() of G4DNARPWBAIonisationModel"
           << G4endl;
  }
  G4double k = particle->GetKineticEnergy();
  if(InEnergyLimit(k))
  {
    G4ParticleMomentum primaryDirection = particle->GetMomentumDirection();
    G4double particleMass  = particle->GetDefinition()->GetPDGMass();
    G4double totalEnergy   = k + particleMass;
    G4double pSquare       = k * (totalEnergy + particleMass);
    G4double totalMomentum = std::sqrt(pSquare);
    G4int ionizationShell;

    if(!fasterCode)
    {
      ionizationShell = RandomSelect(k);
    }
    else
    {
      // fasterCode = true
      do
      {
        ionizationShell = RandomSelect(k);
      } while(k < 19 * eV && ionizationShell == 2 &&
              particle->GetDefinition() == G4Electron::ElectronDefinition());
    }

    G4double bindingEnergy = waterStructure.IonisationEnergy(ionizationShell);

    // SI: additional protection if tcs interpolation method is modified
    if(k < bindingEnergy)
    {
      return;
    }
    //
    G4double secondaryKinetic;
    if(!fasterCode)
    {
      secondaryKinetic = RandomizeEjectedElectronEnergy(k, ionizationShell);
    }
    else
    {
      secondaryKinetic =
        RandomizeEjectedElectronEnergyFromCumulatedDcs(k, ionizationShell);
    }

    G4int Z = 8;  // water Z (6 Oxygen + 2 hydrogen)
    G4ThreeVector deltaDirection =
      GetAngularDistribution()->SampleDirectionForShell(
        particle, secondaryKinetic, Z, ionizationShell, couple->GetMaterial());

    if(secondaryKinetic > 0){
      auto dp = new G4DynamicParticle(G4Electron::Electron(), deltaDirection,
                                      secondaryKinetic);
      fvect->push_back(dp);
    }

    if(particle->GetDefinition() == G4Electron::ElectronDefinition()){
      G4double deltaTotalMomentum = std::sqrt(
        secondaryKinetic * (secondaryKinetic + 2. * electron_mass_c2));

      G4double finalPx = totalMomentum * primaryDirection.x() -
                         deltaTotalMomentum * deltaDirection.x();
      G4double finalPy = totalMomentum * primaryDirection.y() -
                         deltaTotalMomentum * deltaDirection.y();
      G4double finalPz = totalMomentum * primaryDirection.z() -
                         deltaTotalMomentum * deltaDirection.z();
      G4double finalMomentum =
        std::sqrt(finalPx * finalPx + finalPy * finalPy + finalPz * finalPz);
      finalPx /= finalMomentum;
      finalPy /= finalMomentum;
      finalPz /= finalMomentum;
      G4ThreeVector direction;
      direction.set(finalPx, finalPy, finalPz);
      fParticleChangeForGamma->ProposeMomentumDirection(direction.unit());
    }
    else
    {
      fParticleChangeForGamma->ProposeMomentumDirection(primaryDirection);
    }

    // AM: sample deexcitation
    // here we assume that H_{2}O electronic levels are the same as Oxygen.
    // this can be considered true with a rough 10% error in energy on K-shell,

    size_t secNumberInit;  // need to know at a certain point the energy of
                           // secondaries
    size_t
      secNumberFinal;  // So I'll make the diference and then sum the energies

    G4double scatteredEnergy = k - bindingEnergy - secondaryKinetic;

    // SI: only atomic deexcitation from K shell is considered
    if((fAtomDeexcitation != nullptr) && ionizationShell == 4)
    {
      const G4AtomicShell* shell =
        fAtomDeexcitation->GetAtomicShell(Z, G4AtomicShellEnumerator(0));
      secNumberInit = fvect->size();
      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, 0, 0);
      secNumberFinal = fvect->size();

      if(secNumberFinal > secNumberInit){
        for(size_t i = secNumberInit; i < secNumberFinal; ++i){
          if(bindingEnergy >= ((*fvect)[i])->GetKineticEnergy())
          {
            bindingEnergy -= ((*fvect)[i])->GetKineticEnergy();
          }else{
            delete(*fvect)[i];
            (*fvect)[i] = nullptr;
          }
        }
      }
    }

    // This should never happen
    if(bindingEnergy < 0.0)
    {
      G4Exception("G4DNARPWBAIonisatioModel::SampleSecondaries()", "em2050",
                  FatalException, "Negative local energy deposit");
    }

    // bindingEnergy has been decreased
    // by the amount of energy taken away by deexc. products
    if(!statCode){
      fParticleChangeForGamma->SetProposedKineticEnergy(scatteredEnergy);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(bindingEnergy);
    }else{
      fParticleChangeForGamma->SetProposedKineticEnergy(k);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(k - scatteredEnergy);
    }
    const G4Track* theIncomingTrack =
      fParticleChangeForGamma->GetCurrentTrack();
    G4DNAChemistryManager::Instance()->CreateWaterMolecule(
      eIonizedMolecule, ionizationShell, theIncomingTrack);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARPWBAIonisationModel::RandomizeEjectedElectronEnergy(
  const G4double& k, const G4int& shell)
{
  G4double maximumKineticEnergyTransfer =
    4. * (electron_mass_c2 / proton_mass_c2) * k;
  G4double IonisationEnergyInShell = waterStructure.IonisationEnergy(shell);
  G4double kIneV = k / eV;

  G4double crossSectionMaximum = 0.;
  for(G4double value = IonisationEnergyInShell;
      value <= 4. * IonisationEnergyInShell; value += 0.1 * eV)
  {
    G4double differentialCrossSection =
      DifferentialCrossSection(kIneV, value / eV, shell);
    if(differentialCrossSection >= crossSectionMaximum)
    {
      crossSectionMaximum = differentialCrossSection;
    }
  }

  G4double secondaryElectronKineticEnergy;
  do
  {
    secondaryElectronKineticEnergy =
      G4UniformRand() * maximumKineticEnergyTransfer;
  } while(G4UniformRand() * crossSectionMaximum >=
          DifferentialCrossSection(kIneV,
             (secondaryElectronKineticEnergy +
               IonisationEnergyInShell) / eV, shell));

  return secondaryElectronKineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4DNARPWBAIonisationModel::DifferentialCrossSection(
  const G4double& kine, const G4double& energyTransfer,
  const G4int& ionizationLevelIndex)
{
  G4double k     = kine;
  G4double sigma = 0.;
  if(energyTransfer >=
     waterStructure.IonisationEnergy(ionizationLevelIndex) / eV)
  {
    G4double valueT1  = 0;
    G4double valueT2  = 0;
    G4double valueE21 = 0;
    G4double valueE22 = 0;
    G4double valueE12 = 0;
    G4double valueE11 = 0;

    G4double xs11 = 0;
    G4double xs12 = 0;
    G4double xs21 = 0;
    G4double xs22 = 0;

    // Protection against out of boundary access - proton case : 100 MeV

    if(k == pTdummyVec.back())
    {
      k = k * (1. - 1e-12);
    }
    // k should be in eV and energy transfer eV also
    auto t2 = std::upper_bound(pTdummyVec.begin(), pTdummyVec.end(), k);
    auto t1 = t2 - 1;

    auto e12 = std::upper_bound(pVecm[(*t1)].begin(), pVecm[(*t1)].end(),
                                energyTransfer);
    auto e11 = e12 - 1;

    auto e22 = std::upper_bound(pVecm[(*t2)].begin(), pVecm[(*t2)].end(),
                                energyTransfer);
    auto e21 = e22 - 1;

    valueT1  = *t1;
    valueT2  = *t2;
    valueE21 = *e21;
    valueE22 = *e22;
    valueE12 = *e12;
    valueE11 = *e11;

    xs11 = pDiffCrossSectionData[ionizationLevelIndex][valueT1][valueE11];
    xs12 = pDiffCrossSectionData[ionizationLevelIndex][valueT1][valueE12];
    xs21 = pDiffCrossSectionData[ionizationLevelIndex][valueT2][valueE21];
    xs22 = pDiffCrossSectionData[ionizationLevelIndex][valueT2][valueE22];

    G4double xsProduct = xs11 * xs12 * xs21 * xs22;
    if(xsProduct != 0.)
    {
      sigma =
        QuadInterpolator(valueE11, valueE12, valueE21, valueE22, xs11, xs12,
                         xs21, xs22, valueT1, valueT2, k, energyTransfer);
    }
  }
  return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARPWBAIonisationModel::Interpolate(const G4double& e1,
                                                const G4double& e2,
                                                const G4double& e,
                                                const G4double& xs1,
                                                const G4double& xs2)
{
  G4double value = 0.;

  // Log-log interpolation by default

  if(e1 != 0 && e2 != 0 && (std::log10(e2) - std::log10(e1)) != 0 &&
     !fasterCode)
  {
    G4double a =
      (std::log10(xs2) - std::log10(xs1)) / (std::log10(e2) - std::log10(e1));
    G4double b     = std::log10(xs2) - a * std::log10(e2);
    G4double sigma = a * std::log10(e) + b;
    value          = (std::pow(10., sigma));
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
    value       = std::pow(10., (d1 + (d2 - d1) * (e - e1) / (e2 - e1)));
  }

  // Switch to lin-lin interpolation for faster code
  // in case one of xs1 or xs2 (=cum proba) value is zero

  if((e2 - e1) != 0 && (xs1 == 0 || xs2 == 0) && fasterCode)
  {
    G4double d1 = xs1;
    G4double d2 = xs2;
    value       = (d1 + (d2 - d1) * (e - e1) / (e2 - e1));
  }
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARPWBAIonisationModel::QuadInterpolator(
  const G4double& e11, const G4double& e12, const G4double& e21,
  const G4double& e22, const G4double& xs11, const G4double& xs12,
  const G4double& xs21, const G4double& xs22, const G4double& t1,
  const G4double& t2, const G4double& t, const G4double& e)
{
  G4double interpolatedvalue1 = Interpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = Interpolate(e21, e22, e, xs21, xs22);
  G4double value =
    Interpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);

  return value;
}

G4double G4DNARPWBAIonisationModel::GetPartialCrossSection(
  const G4Material* /*material*/, G4int level,
  const G4ParticleDefinition* particle, G4double kineticEnergy)
{
  if(fpTotalCrossSection != nullptr && particle != fProtonDef)
  {
    G4Exception("G4DNARPWBAIonisationModel::GetPartialCrossSection", "em0010",
                FatalException, "Model not applicable to particle type.");
  }
  return fpTotalCrossSection->GetComponent(level)->FindValue(kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNARPWBAIonisationModel::RandomSelect(G4double k)
{
  if(fpTotalCrossSection == nullptr)
  {
    G4Exception("G4DNARPWBAIonisationModel::RandomSelect", "em0010",
                FatalException, "Model not applicable to particle type.");
  }
  else
  {
    auto valuesBuffer = new G4double[fpTotalCrossSection->NumberOfComponents()];
    const G4int n = (G4int)fpTotalCrossSection->NumberOfComponents();
    G4int i(n);
    G4double value = 0.;
    while(i > 0)
    {
      --i;
      valuesBuffer[i] = fpTotalCrossSection->GetComponent(i)->FindValue(k);
      value += valuesBuffer[i];
    }
    value *= G4UniformRand();
    i = n;

    while(i > 0)
    {
      --i;
      if(valuesBuffer[i] > value)
      {
        delete[] valuesBuffer;
        return i;
      }
      value -= valuesBuffer[i];
    }
    delete[] valuesBuffer;
  }
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double
G4DNARPWBAIonisationModel::RandomizeEjectedElectronEnergyFromCumulatedDcs(
  const G4double& k, const G4int& shell)
{
  G4double random                         = G4UniformRand();
  G4double secondaryKineticEnergy =
    TransferedEnergy(k / eV, shell, random) * eV -
    waterStructure.IonisationEnergy(shell);
  if(secondaryKineticEnergy < 0.)
  {
    return 0.;
  }
  return secondaryKineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNARPWBAIonisationModel::TransferedEnergy(G4double k,
                                                     G4int ionizationLevelIndex,
                                                     const G4double& random)
{
  G4double nrj         = 0.;
  G4double valueK1     = 0;
  G4double valueK2     = 0;
  G4double valuePROB21 = 0;
  G4double valuePROB22 = 0;
  G4double valuePROB12 = 0;
  G4double valuePROB11 = 0;
  G4double nrjTransf11 = 0;
  G4double nrjTransf12 = 0;
  G4double nrjTransf21 = 0;
  G4double nrjTransf22 = 0;
  // Protection against out of boundary access - proton case : 100 MeV
  if(k == pTdummyVec.back())
  {
    k = k * (1. - 1e-12);
  }
  // k should be in eV

  auto k2 = std::upper_bound(pTdummyVec.begin(), pTdummyVec.end(), k);
  auto k1 = k2 - 1;

  // SI : the following condition avoids situations where random > last vector
  // element,
  //      for eg. when the last element is zero
  if(random <= pProbaShellMap[ionizationLevelIndex][(*k1)].back() &&
     random <= pProbaShellMap[ionizationLevelIndex][(*k2)].back())
  {
    auto prob12 = std::upper_bound(
      pProbaShellMap[ionizationLevelIndex][(*k1)].begin(),
      pProbaShellMap[ionizationLevelIndex][(*k1)].end(), random);
    auto prob11 = prob12 - 1;
    auto prob22 = std::upper_bound(
      pProbaShellMap[ionizationLevelIndex][(*k2)].begin(),
      pProbaShellMap[ionizationLevelIndex][(*k2)].end(), random);

    auto prob21 = prob22 - 1;

    valueK1     = *k1;
    valueK2     = *k2;
    valuePROB21 = *prob21;
    valuePROB22 = *prob22;
    valuePROB12 = *prob12;
    valuePROB11 = *prob11;

    nrjTransf11 = pNrjTransfData[ionizationLevelIndex][valueK1][valuePROB11];
    nrjTransf12 = pNrjTransfData[ionizationLevelIndex][valueK1][valuePROB12];
    nrjTransf21 = pNrjTransfData[ionizationLevelIndex][valueK2][valuePROB21];
    nrjTransf22 = pNrjTransfData[ionizationLevelIndex][valueK2][valuePROB22];
  }

  // Avoids cases where cum xs is zero for k1 and is not for k2 (with always
  // k1<k2)

  if(random > pProbaShellMap[ionizationLevelIndex][(*k1)].back())
  {
    auto prob22 = std::upper_bound(
      pProbaShellMap[ionizationLevelIndex][(*k2)].begin(),
      pProbaShellMap[ionizationLevelIndex][(*k2)].end(), random);
    auto prob21 = prob22 - 1;

    valueK1     = *k1;
    valueK2     = *k2;
    valuePROB21 = *prob21;
    valuePROB22 = *prob22;
    nrjTransf21 = pNrjTransfData[ionizationLevelIndex][valueK2][valuePROB21];
    nrjTransf22 = pNrjTransfData[ionizationLevelIndex][valueK2][valuePROB22];

    G4double interpolatedvalue2 =
      Interpolate(valuePROB21, valuePROB22, random, nrjTransf21, nrjTransf22);

    G4double value = Interpolate(valueK1, valueK2, k, 0., interpolatedvalue2);
    return value;
  }
  G4double nrjTransfProduct =
    nrjTransf11 * nrjTransf12 * nrjTransf21 * nrjTransf22;

  if(nrjTransfProduct != 0.)
  {
    nrj = QuadInterpolator(valuePROB11, valuePROB12, valuePROB21, valuePROB22,
                           nrjTransf11, nrjTransf12, nrjTransf21, nrjTransf22,
                           valueK1, valueK2, k, random);
  }
  return nrj;
}
