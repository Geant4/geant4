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

#include "G4DNAMillerGreenExcitationModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4Exp.hh"
#include "G4Pow.hh"
#include "G4Alpha.hh"

static G4Pow * gpow = G4Pow::GetInstance();
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAMillerGreenExcitationModel::G4DNAMillerGreenExcitationModel(const G4ParticleDefinition*,
                                                                 const G4String& nam)
:G4VEmModel(nam),isInitialised(false)
{
  fpMolWaterDensity = 0;

  nLevels=0;
  kineticEnergyCorrection[0]=0.;
  kineticEnergyCorrection[1]=0.;
  kineticEnergyCorrection[2]=0.;
  kineticEnergyCorrection[3]=0.;

  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if( verboseLevel>0 )
  {
    G4cout << "Miller & Green excitation model is constructed " << G4endl;
  }
  fParticleChangeForGamma = 0;

  // Selection of stationary mode

  statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4DNAMillerGreenExcitationModel::~G4DNAMillerGreenExcitationModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAMillerGreenExcitationModel::Initialise(const G4ParticleDefinition* particle,
                                                 const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
    G4cout << "Calling G4DNAMillerGreenExcitationModel::Initialise()" << G4endl;

  // Energy limits

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  protonDef = G4Proton::ProtonDefinition();
  hydrogenDef = instance->GetIon("hydrogen");
  alphaPlusPlusDef = G4Alpha::Alpha();
  alphaPlusDef = instance->GetIon("alpha+");
  heliumDef = instance->GetIon("helium");

  G4String proton;
  G4String hydrogen;
  G4String alphaPlusPlus;
  G4String alphaPlus;
  G4String helium;

  // LIMITS AND CONSTANTS

  proton = protonDef->GetParticleName();
  lowEnergyLimit[proton] = 10. * eV;
  highEnergyLimit[proton] = 500. * keV;
  
  kineticEnergyCorrection[0] = 1.;
  slaterEffectiveCharge[0][0] = 0.;
  slaterEffectiveCharge[1][0] = 0.;
  slaterEffectiveCharge[2][0] = 0.;
  sCoefficient[0][0] = 0.;
  sCoefficient[1][0] = 0.;
  sCoefficient[2][0] = 0.;

  hydrogen = hydrogenDef->GetParticleName();
  lowEnergyLimit[hydrogen] = 10. * eV;
  highEnergyLimit[hydrogen] = 500. * keV;
  
  kineticEnergyCorrection[0] = 1.;
  slaterEffectiveCharge[0][0] = 0.;
  slaterEffectiveCharge[1][0] = 0.;
  slaterEffectiveCharge[2][0] = 0.;
  sCoefficient[0][0] = 0.;
  sCoefficient[1][0] = 0.;
  sCoefficient[2][0] = 0.;

  alphaPlusPlus = alphaPlusPlusDef->GetParticleName();
  lowEnergyLimit[alphaPlusPlus] = 1. * keV;
  highEnergyLimit[alphaPlusPlus] = 400. * MeV;

  kineticEnergyCorrection[1] = 0.9382723/3.727417;
  slaterEffectiveCharge[0][1]=0.;
  slaterEffectiveCharge[1][1]=0.;
  slaterEffectiveCharge[2][1]=0.;
  sCoefficient[0][1]=0.;
  sCoefficient[1][1]=0.;
  sCoefficient[2][1]=0.;

  alphaPlus = alphaPlusDef->GetParticleName();
  lowEnergyLimit[alphaPlus] = 1. * keV;
  highEnergyLimit[alphaPlus] = 400. * MeV;

  kineticEnergyCorrection[2] = 0.9382723/3.727417;
  slaterEffectiveCharge[0][2]=2.0;

  // Following values provided by M. Dingfelder
  slaterEffectiveCharge[1][2]=2.00;
  slaterEffectiveCharge[2][2]=2.00;
  //
  sCoefficient[0][2]=0.7;
  sCoefficient[1][2]=0.15;
  sCoefficient[2][2]=0.15;

  helium = heliumDef->GetParticleName();
  lowEnergyLimit[helium] = 1. * keV;
  highEnergyLimit[helium] = 400. * MeV;
  
  kineticEnergyCorrection[3] = 0.9382723/3.727417;
  slaterEffectiveCharge[0][3]=1.7;
  slaterEffectiveCharge[1][3]=1.15;
  slaterEffectiveCharge[2][3]=1.15;
  sCoefficient[0][3]=0.5;
  sCoefficient[1][3]=0.25;
  sCoefficient[2][3]=0.25;

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

  if (particle==alphaPlusPlusDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[alphaPlusPlus]);
    SetHighEnergyLimit(highEnergyLimit[alphaPlusPlus]);
  }

  if (particle==alphaPlusDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[alphaPlus]);
    SetHighEnergyLimit(highEnergyLimit[alphaPlus]);
  }

  if (particle==heliumDef)
  {
    SetLowEnergyLimit(lowEnergyLimit[helium]);
    SetHighEnergyLimit(highEnergyLimit[helium]);
  }

  //

  nLevels = waterExcitation.NumberOfLevels();

  //
  if( verboseLevel>0 )
  {
    G4cout << "Miller & Green excitation model is initialized " << G4endl
           << "Energy range: "
           << LowEnergyLimit() / eV << " eV - "
           << HighEnergyLimit() / keV << " keV for "
           << particle->GetParticleName()
           << G4endl;
  }

  // Initialize water density pointer
  fpMolWaterDensity = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

  if (isInitialised) { return; }
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4DNAMillerGreenExcitationModel::CrossSectionPerVolume(const G4Material* material,
                                                                const G4ParticleDefinition* particleDefinition,
                                                                G4double k,
                                                                G4double,
                                                                G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of G4DNAMillerGreenExcitationModel" << G4endl;

  // Calculate total cross section for model

  if (
       particleDefinition != protonDef
       &&
       particleDefinition != hydrogenDef
       &&
       particleDefinition != alphaPlusPlusDef
       &&
       particleDefinition != alphaPlusDef
       &&
       particleDefinition != heliumDef
     )

     return 0;

  G4double lowLim = 0;
  G4double highLim = 0;
  G4double crossSection = 0.;

  G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];

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

  if (k >= lowLim && k <= highLim)
  {
    crossSection = Sum(k,particleDefinition);

    // add ONE or TWO electron-water excitation for alpha+ and helium
    /*
      if ( particleDefinition == alphaPlusDef
           ||
           particleDefinition == heliumDef
         )
      {

      G4DNAEmfietzoglouExcitationModel * excitationXS = new G4DNAEmfietzoglouExcitationModel();
          excitationXS->Initialise(G4Electron::ElectronDefinition());

      G4double sigmaExcitation=0;
      G4double tmp =0.;

      if (k*0.511/3728 > 8.23*eV && k*0.511/3728 < 10*MeV ) sigmaExcitation =
        excitationXS->CrossSectionPerVolume(material,G4Electron::ElectronDefinition(),k*0.511/3728,tmp,tmp)
        /material->GetAtomicNumDensityVector()[1];

      if ( particleDefinition == alphaPlusDef )
        crossSection = crossSection +  sigmaExcitation ;

      if ( particleDefinition == heliumDef )
        crossSection = crossSection + 2*sigmaExcitation ;

      delete excitationXS;

          // Alternative excitation model

          G4DNABornExcitationModel * excitationXS = new G4DNABornExcitationModel();
          excitationXS->Initialise(G4Electron::ElectronDefinition());

      G4double sigmaExcitation=0;
      G4double tmp=0;

      if (k*0.511/3728 > 9*eV && k*0.511/3728 < 1*MeV ) sigmaExcitation =
        excitationXS->CrossSectionPerVolume(material,G4Electron::ElectronDefinition(),k*0.511/3728,tmp,tmp)
        /material->GetAtomicNumDensityVector()[1];

      if ( particleDefinition == alphaPlusDef )
        crossSection = crossSection +  sigmaExcitation ;

      if ( particleDefinition == heliumDef )
        crossSection = crossSection + 2*sigmaExcitation ;

      delete excitationXS;

      }
    */      

  }

  if (verboseLevel > 2)
  {
    G4cout << "__________________________________" << G4endl;
    G4cout << "G4DNAMillerGreenExcitationModel - XS INFO START" << G4endl;
    G4cout << "Kinetic energy(eV)=" << k/eV << " particle : " << particleDefinition->GetParticleName() << G4endl;
    G4cout << "Cross section per water molecule (cm^2)=" << crossSection/cm/cm << G4endl;
    G4cout << "Cross section per water molecule (cm^-1)=" << crossSection*waterDensity/(1./cm) << G4endl;
    // G4cout << " - Cross section per water molecule (cm^-1)=" << sigma*material->GetAtomicNumDensityVector()[1]/(1./cm) << G4endl;
    G4cout << "G4DNAMillerGreenExcitationModel - XS INFO END" << G4endl;
  }

    return crossSection*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4DNAMillerGreenExcitationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
                                                        const G4MaterialCutsCouple* /*couple*/,
                                                        const G4DynamicParticle* aDynamicParticle,
                                                        G4double,
                                                        G4double)
{

  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4DNAMillerGreenExcitationModel" << G4endl;

  G4double particleEnergy0 = aDynamicParticle->GetKineticEnergy();

  G4int level = RandomSelect(particleEnergy0,aDynamicParticle->GetDefinition());

  // Dingfelder's excitation levels
  const G4double excitation[]={ 8.17*eV, 10.13*eV, 11.31*eV, 12.91*eV, 14.50*eV};
  G4double excitationEnergy = excitation[level];

  G4double newEnergy = 0.;
  
  if (!statCode) newEnergy = particleEnergy0 - excitationEnergy;

  else newEnergy = particleEnergy0;
  
  if (newEnergy>0)
  {
    fParticleChangeForGamma->ProposeMomentumDirection(aDynamicParticle->GetMomentumDirection());
    fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);

    const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
    G4DNAChemistryManager::Instance()->CreateWaterMolecule(eExcitedMolecule,
     level, theIncomingTrack);

  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNAMillerGreenExcitationModel::GetPartialCrossSection(const G4Material*,
                                          G4int level,
                                          const G4ParticleDefinition* particleDefinition,
                                          G4double kineticEnergy)
{
  return PartialCrossSection(kineticEnergy, level, particleDefinition);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNAMillerGreenExcitationModel::PartialCrossSection(G4double k, G4int excitationLevel, 
                                                              const G4ParticleDefinition* particleDefinition)
{
  //                               ( ( z * aj ) ^ omegaj ) * ( t - ej ) ^ nu
  // sigma(t) = zEff^2 * sigma0 * --------------------------------------------
  //                               jj ^ ( omegaj + nu ) + t ^ ( omegaj + nu )
  //
  // where t is the kinetic energy corrected by Helium mass over proton mass for Helium ions
  //
  // zEff is:
  //  1 for protons
  //  2 for alpha++
  //  and  2 - c1 S_1s - c2 S_2s - c3 S_2p for alpha+ and He
  //
  // Dingfelder et al., RPC 59, 255-275, 2000 from Miller and Green (1973)
  // Formula (34) and Table 2

  const G4double sigma0(1.E+8 * barn);
  const G4double nu(1.);
  const G4double aj[]={876.*eV, 2084.* eV, 1373.*eV, 692.*eV, 900.*eV};
  const G4double jj[]={19820.*eV, 23490.*eV, 27770.*eV, 30830.*eV, 33080.*eV};
  const G4double omegaj[]={0.85, 0.88, 0.88, 0.78, 0.78};

  // Dingfelder's excitation levels
  const G4double Eliq[5]={ 8.17*eV, 10.13*eV, 11.31*eV, 12.91*eV, 14.50*eV};

  G4int particleTypeIndex = 0;
 
  if (particleDefinition == protonDef) particleTypeIndex=0;
  if (particleDefinition == hydrogenDef) particleTypeIndex=0;
  if (particleDefinition == alphaPlusPlusDef) particleTypeIndex=1;
  if (particleDefinition == alphaPlusDef) particleTypeIndex=2;
  if (particleDefinition == heliumDef) particleTypeIndex=3;

  G4double tCorrected;
  tCorrected = k * kineticEnergyCorrection[particleTypeIndex];

  // SI - added protection
  if (tCorrected < Eliq[excitationLevel]) return 0;
  //

  G4int z = 10;

  G4double numerator;
  numerator = gpow->powA(z * aj[excitationLevel], omegaj[excitationLevel]) *
              gpow->powA(tCorrected - Eliq[excitationLevel], nu);

  // H case : see S. Uehara et al. IJRB 77, 2, 139-154 (2001) - section 3.3

  if (particleDefinition == hydrogenDef)
      numerator = gpow->powA(z * 0.75*aj[excitationLevel], omegaj[excitationLevel]) *
                  gpow->powA(tCorrected - Eliq[excitationLevel], nu);


  G4double power;
  power = omegaj[excitationLevel] + nu;

  G4double denominator;
  denominator = gpow->powA(jj[excitationLevel], power) + gpow->powA(tCorrected, power);

  G4double zEff = particleDefinition->GetPDGCharge() / eplus + particleDefinition->GetLeptonNumber();

  zEff -= ( sCoefficient[0][particleTypeIndex] * S_1s(k, Eliq[excitationLevel], slaterEffectiveCharge[0][particleTypeIndex], 1.) +
          sCoefficient[1][particleTypeIndex] * S_2s(k, Eliq[excitationLevel], slaterEffectiveCharge[1][particleTypeIndex], 2.) +
          sCoefficient[2][particleTypeIndex] * S_2p(k, Eliq[excitationLevel], slaterEffectiveCharge[2][particleTypeIndex], 2.) );

  if (particleDefinition == hydrogenDef) zEff = 1.;

  G4double cross = sigma0 * zEff * zEff * numerator / denominator;


  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNAMillerGreenExcitationModel::RandomSelect(G4double k,const G4ParticleDefinition* particle)
{
  G4int i = nLevels;
  G4double value = 0.;
  std::deque<G4double> values;

  if ( particle == alphaPlusPlusDef ||
       particle == protonDef||
       particle == hydrogenDef  ||
       particle == alphaPlusDef  ||
       particle == heliumDef
       )
  {
    while (i > 0)
    {
      i--;
      G4double partial = PartialCrossSection(k,i,particle);
      values.push_front(partial);
      value += partial;
    }

    value *= G4UniformRand();

    i = nLevels;

    while (i > 0)
    {
      i--;
      if (values[i] > value) return i;
      value -= values[i];
    }
  }

  /*
  // add ONE or TWO electron-water excitation for alpha+ and helium

  if ( particle == alphaPlusDef
       ||
       particle == heliumDef
     )
  {
    while (i>0)
    {
      i--;

          G4DNAEmfietzoglouExcitationModel * excitationXS = new G4DNAEmfietzoglouExcitationModel();
          excitationXS->Initialise(G4Electron::ElectronDefinition());

      G4double sigmaExcitation=0;

      if (k*0.511/3728 > 8.23*eV && k*0.511/3728 < 10*MeV ) sigmaExcitation = excitationXS->PartialCrossSection(k*0.511/3728,i);

      G4double partial = PartialCrossSection(k,i,particle);

      if (particle == alphaPlusDef) partial = PartialCrossSection(k,i,particle) + sigmaExcitation;
      if (particle == heliumDef) partial = PartialCrossSection(k,i,particle) + 2*sigmaExcitation;

      values.push_front(partial);
      value += partial;
      delete excitationXS;
    }

    value*=G4UniformRand();

    i=5;
    while (i>0)
    {
      i--;

      if (values[i]>value) return i;

      value-=values[i];
    }
  }
  */

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNAMillerGreenExcitationModel::Sum(G4double k, const G4ParticleDefinition* particle)
{
  G4double totalCrossSection = 0.;

  for (G4int i=0; i<nLevels; i++)
  {
    totalCrossSection += PartialCrossSection(k,i,particle);
  }
  return totalCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNAMillerGreenExcitationModel::S_1s(G4double t,
                                               G4double energyTransferred,
                                               G4double _slaterEffectiveCharge,
                                               G4double shellNumber)
{
  // 1 - e^(-2r) * ( 1 + 2 r + 2 r^2)
  // Dingfelder, in Chattanooga 2005 proceedings, formula (7)

  G4double r = R(t, energyTransferred, _slaterEffectiveCharge, shellNumber);
  G4double value = 1. - G4Exp(-2 * r) * ( ( 2. * r + 2. ) * r + 1. );

  return value;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNAMillerGreenExcitationModel::S_2s(G4double t,
                                               G4double energyTransferred,
                                               G4double _slaterEffectiveCharge,
                                               G4double shellNumber)
{
  // 1 - e^(-2 r) * ( 1 + 2 r + 2 r^2 + 2 r^4)
  // Dingfelder, in Chattanooga 2005 proceedings, formula (8)

  G4double r = R(t, energyTransferred, _slaterEffectiveCharge, shellNumber);
  G4double value =  1. - G4Exp(-2 * r) * (((2. * r * r + 2.) * r + 2.) * r + 1.);

  return value;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNAMillerGreenExcitationModel::S_2p(G4double t,
                                               G4double energyTransferred,
                                               G4double _slaterEffectiveCharge,
                                               G4double shellNumber)
{
  // 1 - e^(-2 r) * ( 1 + 2 r + 2 r^2 + 4/3 r^3 + 2/3 r^4)
  // Dingfelder, in Chattanooga 2005 proceedings, formula (9)

  G4double r = R(t, energyTransferred, _slaterEffectiveCharge, shellNumber);
  G4double value =  1. - G4Exp(-2 * r) * (((( 2./3. * r + 4./3.) * r + 2.) * r + 2.) * r  + 1.);

  return value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNAMillerGreenExcitationModel::R(G4double t,
                                            G4double energyTransferred,
                                            G4double _slaterEffectiveCharge,
                                            G4double shellNumber)
{
  // tElectron = m_electron / m_alpha * t
  // Dingfelder, in Chattanooga 2005 proceedings, p 4

  G4double tElectron = 0.511/3728. * t;

  // The following is provided by M. Dingfelder
  G4double H = 2.*13.60569172 * eV;
  G4double value = std::sqrt ( 2. * tElectron / H ) / ( energyTransferred / H ) *  (_slaterEffectiveCharge/shellNumber);

  return value;
}

