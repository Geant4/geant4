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
// Author: H. N. Tran (Ton Duc Thang University)
// p, H, He, He+ and He++ models are assumed identical
// NIMB 343, 132-137 (2015)
//
// The Geant4-DNA web site is available at http://geant4-dna.org
//

#include "G4DNAIonElasticModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4ParticleTable.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAIonElasticModel::G4DNAIonElasticModel (const G4ParticleDefinition*,
                                            const G4String& nam) :
    G4VEmModel(nam), isInitialised(false)
{
  killBelowEnergy = 100 * eV;
  lowEnergyLimit = 0 * eV;
  highEnergyLimit = 1 * MeV;
  SetLowEnergyLimit(lowEnergyLimit);
  SetHighEnergyLimit(highEnergyLimit);

  verboseLevel = 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if(verboseLevel > 0)
  {
    G4cout << "Ion elastic model is constructed " << G4endl<< "Energy range: "
    << lowEnergyLimit / eV << " eV - "
    << highEnergyLimit / MeV << " MeV"
    << G4endl;
  }
  
  fParticleChangeForGamma = 0;
  fpMolWaterDensity = 0;
  fpTableData = 0;
  fParticle_Mass = -1;

  // Selection of stationary mode

  statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAIonElasticModel::~G4DNAIonElasticModel ()
{
  // For total cross section
  if(fpTableData) delete fpTableData;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4DNAIonElasticModel::Initialise (
    const G4ParticleDefinition* particleDefinition,
    const G4DataVector& /*cuts*/)
{

  if(verboseLevel > 3)
  {
    G4cout << "Calling G4DNAIonElasticModel::Initialise()" << G4endl;
  }

  // Energy limits

  if (LowEnergyLimit() < lowEnergyLimit)
  {
    G4cout << "G4DNAIonElasticModel: low energy limit increased from " <<
    LowEnergyLimit()/eV << " eV to " << lowEnergyLimit/eV << " eV" << G4endl;
    SetLowEnergyLimit(lowEnergyLimit);
  }

  if (HighEnergyLimit() > highEnergyLimit)
  {
    G4cout << "G4DNAIonElasticModel: high energy limit decreased from " <<
    HighEnergyLimit()/MeV << " MeV to " << highEnergyLimit/MeV << " MeV" << G4endl;
    SetHighEnergyLimit(highEnergyLimit);
  }

  // Reading of data files

  G4double scaleFactor = 1e-16*cm*cm;

  const char *path = G4FindDataDir("G4LEDATA");

  if (!path)
  {
    G4Exception("G4IonElasticModel::Initialise","em0006",
        FatalException,"G4LEDATA environment variable not set.");
    return;
  }

  G4String totalXSFile;
  std::ostringstream fullFileName;

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  G4ParticleDefinition* protonDef =
  G4ParticleTable::GetParticleTable()->FindParticle("proton");
  G4ParticleDefinition* hydrogenDef = instance->GetIon("hydrogen");
  G4ParticleDefinition* heliumDef = instance->GetIon("helium");
  G4ParticleDefinition* alphaplusDef = instance->GetIon("alpha+");
  G4ParticleDefinition* alphaplusplusDef = instance->GetIon("alpha++");
  G4String proton, hydrogen, helium, alphaplus, alphaplusplus;

  if (
      (particleDefinition == protonDef && protonDef != 0)
      ||
      (particleDefinition == hydrogenDef && hydrogenDef != 0)
  )
  {
    // For total cross section of p,h
    fParticle_Mass = 1.;   
    totalXSFile = "dna/sigma_elastic_proton_HTran";

    // For final state
    fullFileName << path << "/dna/sigmadiff_cumulated_elastic_proton_HTran.dat";
  }

  if (
      (particleDefinition == instance->GetIon("helium") && heliumDef)
      ||
      (particleDefinition == instance->GetIon("alpha+") && alphaplusDef)
      ||
      (particleDefinition == instance->GetIon("alpha++") && alphaplusplusDef)
  )
  {
    // For total cross section of he,he+,he++
    fParticle_Mass = 4.;    
    totalXSFile = "dna/sigma_elastic_alpha_HTran";

    // For final state
    fullFileName << path << "/dna/sigmadiff_cumulated_elastic_alpha_HTran.dat";
  }

  fpTableData = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, eV,scaleFactor );
  fpTableData->LoadData(totalXSFile);
  std::ifstream diffCrossSection(fullFileName.str().c_str());

  if (!diffCrossSection)
  {
    G4ExceptionDescription description;
    description << "Missing data file:"
    <<fullFileName.str().c_str()<< G4endl;
    G4Exception("G4IonElasticModel::Initialise","em0003",
        FatalException,
        description);
  }

  // Added clear for MT

  eTdummyVec.clear();
  eVecm.clear();
  fDiffCrossSectionData.clear();

  //

  eTdummyVec.push_back(0.);

  while(!diffCrossSection.eof())
  {
    G4double tDummy;
    G4double eDummy;
    diffCrossSection>>tDummy>>eDummy;

    // SI : mandatory eVecm initialization

    if (tDummy != eTdummyVec.back())
    {
      eTdummyVec.push_back(tDummy);
      eVecm[tDummy].push_back(0.);
    }

    diffCrossSection>>fDiffCrossSectionData[tDummy][eDummy];

    if (eDummy != eVecm[tDummy].back()) eVecm[tDummy].push_back(eDummy);

  }

  // End final state
  if( verboseLevel>0 )
  {
    if (verboseLevel > 2)
    {
      G4cout << "Loaded cross section files for ion elastic model" << G4endl;
    }
    G4cout << "Ion elastic model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / eV << " eV - "
    << HighEnergyLimit() / MeV << " MeV"
    << G4endl;
  }

  // Initialize water density pointer
  G4DNAMolecularMaterial::Instance()->Initialize();
  fpMolWaterDensity = G4DNAMolecularMaterial::Instance()->
  GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

  if (isInitialised) return;
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4DNAIonElasticModel::CrossSectionPerVolume (const G4Material* material,
                                             const G4ParticleDefinition* p,
                                             G4double ekin, G4double, G4double)
{
  if(verboseLevel > 3)
  {
    G4cout << "Calling CrossSectionPerVolume() of G4DNAIonElasticModel"
           << G4endl;
  }

  // Calculate total cross section for model

  G4double sigma=0;

  G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];

  const G4String& particleName = p->GetParticleName();

  if (ekin <= highEnergyLimit)
  {
    //SI : XS must not be zero otherwise sampling of secondaries method ignored
    if (ekin < killBelowEnergy) return DBL_MAX;
    //

    if (fpTableData != 0) 
    {
      sigma = fpTableData->FindValue(ekin);
    }
    else
    {
      G4Exception("G4DNAIonElasticModel::ComputeCrossSectionPerVolume","em0002",
          FatalException,"Model not applicable to particle type.");
    }
  }

  if (verboseLevel > 2)
  {
    G4cout << "__________________________________" << G4endl;
    G4cout << "G4DNAIonElasticModel - XS INFO START" << G4endl;
    G4cout << "Kinetic energy(eV)=" << ekin/eV << " particle : " << particleName << G4endl;
    G4cout << "Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
    G4cout << "Cross section per water molecule (cm^-1)=" << sigma*waterDensity/(1./cm) << G4endl;
    G4cout << "G4DNAIonElasticModel - XS INFO END" << G4endl;
  }

  return sigma*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4DNAIonElasticModel::SampleSecondaries (
    std::vector<G4DynamicParticle*>* /*fvect*/,
    const G4MaterialCutsCouple* /*couple*/,
    const G4DynamicParticle* aDynamicParticle, G4double, G4double)
{

  if(verboseLevel > 3)
  {
    G4cout << "Calling SampleSecondaries() of G4DNAIonElasticModel" << G4endl;
  }

  G4double particleEnergy0 = aDynamicParticle->GetKineticEnergy();

  if (particleEnergy0 < killBelowEnergy)
  {
    fParticleChangeForGamma->SetProposedKineticEnergy(0.);
    fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(particleEnergy0);
    return;
  }

  if (particleEnergy0>= killBelowEnergy && particleEnergy0 <= highEnergyLimit)
  {
    G4double water_mass = 18.;

    G4double thetaCM = RandomizeThetaCM(particleEnergy0, aDynamicParticle->GetDefinition());

    //HT:convert to laboratory system

    G4double theta = std::atan(std::sin(thetaCM*CLHEP::pi/180)
        /(fParticle_Mass/water_mass+std::cos(thetaCM*CLHEP::pi/180)));

    G4double cosTheta= std::cos(theta);

    //

    G4double phi = 2. * CLHEP::pi * G4UniformRand();

    G4ThreeVector zVers = aDynamicParticle->GetMomentumDirection();
    G4ThreeVector xVers = zVers.orthogonal();
    G4ThreeVector yVers = zVers.cross(xVers);

    G4double xDir = std::sqrt(1. - cosTheta*cosTheta);
    G4double yDir = xDir;
    xDir *= std::cos(phi);
    yDir *= std::sin(phi);

    G4ThreeVector zPrimeVers((xDir*xVers + yDir*yVers + cosTheta*zVers));

    fParticleChangeForGamma->ProposeMomentumDirection(zPrimeVers.unit());

    G4double depositEnergyCM = 0;

    //HT: deposited energy
    depositEnergyCM = 4. * particleEnergy0 * fParticle_Mass * water_mass *
    (1-std::cos(thetaCM*CLHEP::pi/180))
    / (2 * std::pow((fParticle_Mass+water_mass),2));

    //SI: added protection particleEnergy0 >= depositEnergyCM
    if (!statCode && (particleEnergy0 >= depositEnergyCM) ) 

      fParticleChangeForGamma->SetProposedKineticEnergy(particleEnergy0 - depositEnergyCM);
    
    else fParticleChangeForGamma->SetProposedKineticEnergy(particleEnergy0);
   
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(depositEnergyCM);    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4DNAIonElasticModel::Theta (G4ParticleDefinition * /*particleDefinition*/,
                             G4double k, G4double integrDiff)
{
  G4double theta = 0.;
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

  // Protection against out of boundary access
  if (k==eTdummyVec.back()) k=k*(1.-1e-12);
  //

  std::vector<G4double>::iterator t2 = std::upper_bound(eTdummyVec.begin(),
                                                      eTdummyVec.end(), k);
  std::vector<G4double>::iterator t1 = t2 - 1;

  std::vector<G4double>::iterator e12 = std::upper_bound(eVecm[(*t1)].begin(),
                                                       eVecm[(*t1)].end(),
                                                       integrDiff);
  std::vector<G4double>::iterator e11 = e12 - 1;

  std::vector<G4double>::iterator e22 = std::upper_bound(eVecm[(*t2)].begin(),
                                                       eVecm[(*t2)].end(),
                                                       integrDiff);
  std::vector<G4double>::iterator e21 = e22 - 1;

  valueT1 = *t1;
  valueT2 = *t2;
  valueE21 = *e21;
  valueE22 = *e22;
  valueE12 = *e12;
  valueE11 = *e11;

  xs11 = fDiffCrossSectionData[valueT1][valueE11];
  xs12 = fDiffCrossSectionData[valueT1][valueE12];
  xs21 = fDiffCrossSectionData[valueT2][valueE21];
  xs22 = fDiffCrossSectionData[valueT2][valueE22];

  if(xs11 == 0 && xs12 == 0 && xs21 == 0 && xs22 == 0) return (0.);

  theta = QuadInterpolator(valueE11, valueE12, valueE21, valueE22, xs11, xs12,
                           xs21, xs22, valueT1, valueT2, k, integrDiff);

  return theta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4DNAIonElasticModel::LinLinInterpolate (G4double e1, G4double e2, G4double e,
                                         G4double xs1, G4double xs2)
{
  G4double d1 = xs1;
  G4double d2 = xs2;
  G4double value = (d1 + (d2 - d1) * (e - e1) / (e2 - e1));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4DNAIonElasticModel::LinLogInterpolate (G4double e1, G4double e2, G4double e,
                                         G4double xs1, G4double xs2)
{
  G4double d1 = std::log(xs1);
  G4double d2 = std::log(xs2);
  G4double value = G4Exp(d1 + (d2 - d1) * (e - e1) / (e2 - e1));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4DNAIonElasticModel::LogLogInterpolate (G4double e1, G4double e2, G4double e,
                                         G4double xs1, G4double xs2)
{
  G4double a = (std::log10(xs2) - std::log10(xs1))
      / (std::log10(e2) - std::log10(e1));
  G4double b = std::log10(xs2) - a * std::log10(e2);
  G4double sigma = a * std::log10(e) + b;
  G4double value = (std::pow(10., sigma));
  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4DNAIonElasticModel::QuadInterpolator (G4double e11, G4double e12,
                                        G4double e21, G4double e22,
                                        G4double xs11, G4double xs12,
                                        G4double xs21, G4double xs22,
                                        G4double t1, G4double t2, G4double t,
                                        G4double e)
{
  // Log-Log
  /*
   G4double interpolatedvalue1 = LogLogInterpolate(e11, e12, e, xs11, xs12);
   G4double interpolatedvalue2 = LogLogInterpolate(e21, e22, e, xs21, xs22);
   G4double value = LogLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
   */

  // Lin-Log
  /*
   G4double interpolatedvalue1 = LinLogInterpolate(e11, e12, e, xs11, xs12);
   G4double interpolatedvalue2 = LinLogInterpolate(e21, e22, e, xs21, xs22);
   G4double value = LinLogInterpolate(t1, t2, t, interpolatedvalue1, interpolatedvalue2);
   */

// Lin-Lin
  G4double interpolatedvalue1 = LinLinInterpolate(e11, e12, e, xs11, xs12);
  G4double interpolatedvalue2 = LinLinInterpolate(e21, e22, e, xs21, xs22);
  G4double value = LinLinInterpolate(t1, t2, t, interpolatedvalue1,
                                     interpolatedvalue2);

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4DNAIonElasticModel::RandomizeThetaCM (
    G4double k, G4ParticleDefinition * particleDefinition)
{
  G4double integrdiff = G4UniformRand();
  return Theta(particleDefinition, k / eV, integrdiff);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4DNAIonElasticModel::SetKillBelowThreshold (G4double threshold)
{
  killBelowEnergy = threshold;

  if(killBelowEnergy < 100 * eV)
  {
    G4cout << "*** WARNING : the G4DNAIonElasticModel class is not "
           "activated below 100 eV !"
           << G4endl;
   }
}

