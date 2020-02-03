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

// Created by Z. Francis

#include "G4DNASancheExcitationModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAMolecularMaterial.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

//#define SANCHE_VERBOSE

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNASancheExcitationModel::G4DNASancheExcitationModel(const G4ParticleDefinition*,
                                                       const G4String& nam) :
    G4VEmModel(nam), isInitialised(false)
{
  fpWaterDensity = 0;

  SetLowEnergyLimit(2.*eV);
  SetHighEnergyLimit(100*eV);
  nLevels = 9;

  verboseLevel = 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

#ifdef SANCHE_VERBOSE
  if (verboseLevel > 0)
  {
    G4cout << "Sanche Excitation model is constructed "
           << G4endl
           << "Energy range: "
           << LowEnergyLimit() / eV << " eV - "
           << HighEnergyLimit() / eV << " eV"
           << G4endl;
  }
#endif
  
  fParticleChangeForGamma = 0;
  fpWaterDensity = 0;

  // Selection of stationary mode

  statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4DNASancheExcitationModel::~G4DNASancheExcitationModel()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void
G4DNASancheExcitationModel::
Initialise(const G4ParticleDefinition* /*particle*/,
           const G4DataVector& /*cuts*/)
{
  
#ifdef SANCHE_VERBOSE
  if (verboseLevel > 3)
  {
    G4cout << "Calling G4DNASancheExcitationModel::Initialise()"
           << G4endl;
  }
#endif
  
  // Energy limits

  if (LowEnergyLimit() < 2.*eV)
  {
    G4Exception("*** WARNING : the G4DNASancheExcitationModel class is not "
                "validated below 2 eV !",
                "", JustWarning, "");
  }

  if (HighEnergyLimit() > 100.*eV)
  {
    G4cout << "G4DNASancheExcitationModel: high energy limit decreased from " <<
    HighEnergyLimit()/eV << " eV to " << 100. << " eV" << G4endl;
    SetHighEnergyLimit(100.*eV);
  }

  //
#ifdef SANCHE_VERBOSE
  if (verboseLevel > 0)
  {
    G4cout << "Sanche Excitation model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / eV << " eV - "
    << HighEnergyLimit() / eV << " eV"
    << G4endl;
  }
#endif

  // Initialize water density pointer
  fpWaterDensity = G4DNAMolecularMaterial::Instance()->
      GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

  if (isInitialised) {return;}

  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;

  char *path = getenv("G4LEDATA");
  std::ostringstream eFullFileName;
  eFullFileName << path << "/dna/sigma_excitationvib_e_sanche.dat";
  std::ifstream input(eFullFileName.str().c_str());

  if (!input)
  {
    G4Exception("G4DNASancheExcitationModel::Initialise","em0003",
        FatalException,"Missing data file:/dna/sigma_excitationvib_e_sanche.dat");
  }

  // March 25th, 2014 - Vaclav Stepan, Sebastien Incerti
  // Added clear for MT
  tdummyVec.clear();
  //

  G4double t;
  G4double xs;

  while(!input.eof())
  {
    input>>t;
    tdummyVec.push_back(t);

    fEnergyLevelXS.push_back(std::vector<G4double>());
    fEnergyTotalXS.push_back(0);
    std::vector<G4double>& levelXS = fEnergyLevelXS.back();
    levelXS.reserve(9);

    // G4cout<<t;

    for(size_t i = 0 ; i < 9 ;++i)
    {
      input>>xs;
      levelXS.push_back(xs);
      fEnergyTotalXS.back() += xs;
      // G4cout <<"  " << levelXS[i];
    }

    // G4cout << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNASancheExcitationModel::CrossSectionPerVolume(const G4Material* material,
                                                           const G4ParticleDefinition*
#ifdef SANCHE_VERBOSE
                                                           particleDefinition
#endif
                                                           ,
                                                           G4double ekin,
                                                           G4double,
                                                           G4double)
{
#ifdef SANCHE_VERBOSE
  if (verboseLevel > 3)
  {
    G4cout << "Calling CrossSectionPerVolume() of G4DNASancheExcitationModel"
           << G4endl;
  }
#endif

  // Calculate total cross section for model

  G4double sigma = 0.;

  G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];

  if (ekin >= LowEnergyLimit() && ekin <= HighEnergyLimit())
    sigma =  TotalCrossSection(ekin);

#ifdef SANCHE_VERBOSE
  if (verboseLevel > 2)
  {
    G4cout << "__________________________________" << G4endl;
    G4cout << "=== G4DNASancheExcitationModel - XS INFO START" << G4endl;
    G4cout << "=== Kinetic energy(eV)=" << ekin/eV << " particle : " << particleDefinition->GetParticleName() << G4endl;
    G4cout << "=== Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
    G4cout << "=== Cross section per water molecule (cm^-1)=" << sigma*waterDensity/(1./cm) << G4endl;
    G4cout << "=== G4DNASancheExcitationModel - XS INFO END" << G4endl;
  }
#endif

  return sigma*2.*waterDensity;
  // see papers for factor 2 description (liquid phase)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4DNASancheExcitationModel::SampleSecondaries(std::vector<
                                                       G4DynamicParticle*>*,
                                                   const G4MaterialCutsCouple*,
                                                   const G4DynamicParticle* aDynamicElectron,
                                                   G4double,
                                                   G4double)
{
#ifdef SANCHE_VERBOSE
  if (verboseLevel > 3)
  {
    G4cout << "Calling SampleSecondaries() of G4DNASancheExcitationModel"
           << G4endl;
  }
#endif

  G4double electronEnergy0 = aDynamicElectron->GetKineticEnergy();
  G4int level = RandomSelect(electronEnergy0);
  G4double excitationEnergy = VibrationEnergy(level); // levels go from 0 to 8
  G4double newEnergy = electronEnergy0 - excitationEnergy;

  /*
   if (electronEnergy0 < highEnergyLimit)
   {
     if (newEnergy >= lowEnergyLimit)
     {
       fParticleChangeForGamma->ProposeMomentumDirection(aDynamicElectron->GetMomentumDirection());
       fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
       fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
     }

     else
     {
       fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
       fParticleChangeForGamma->ProposeLocalEnergyDeposit(electronEnergy0);
     }
   }
   */

  if (electronEnergy0 <= HighEnergyLimit() && newEnergy>0.)
  {

    if (!statCode)     
    {     
      fParticleChangeForGamma->ProposeMomentumDirection(aDynamicElectron->GetMomentumDirection());
      fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
    }

    else 
    {
      fParticleChangeForGamma->ProposeMomentumDirection(aDynamicElectron->GetMomentumDirection());
      fParticleChangeForGamma->SetProposedKineticEnergy(electronEnergy0);
      fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
    }

  }

  //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNASancheExcitationModel::PartialCrossSection(G4double t,
                                                         G4int level)
{
  // Protection against out of boundary access
  if (t/eV==tdummyVec.back()) t=t*(1.-1e-12);
  //

  std::vector<G4double>::iterator t2 = std::upper_bound(tdummyVec.begin(),
                                                      tdummyVec.end(), t / eV);
  std::vector<G4double>::iterator t1 = t2 - 1;

  size_t i1 = t1 - tdummyVec.begin();
  size_t i2 = t2 - tdummyVec.begin();

  G4double sigma = LinInterpolate((*t1), (*t2),
                                t / eV,
                                fEnergyLevelXS[i1][level],
                                fEnergyLevelXS[i2][level]);

  static const G4double conv_factor =  1e-16 * cm * cm;

  sigma *= conv_factor;
  if (sigma == 0.) sigma = 1e-30;
  return (sigma);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNASancheExcitationModel::TotalCrossSection(G4double t)
{
  // Protection against out of boundary access
  if (t/eV==tdummyVec.back()) t=t*(1.-1e-12);
  //

  std::vector<G4double>::iterator t2 = std::upper_bound(tdummyVec.begin(),
                                                      tdummyVec.end(), t / eV);
  std::vector<G4double>::iterator t1 = t2 - 1;

  size_t i1 = t1 - tdummyVec.begin();
  size_t i2 = t2 - tdummyVec.begin();

  G4double sigma = LinInterpolate((*t1), (*t2),
                                t / eV,
                                fEnergyTotalXS[i1],
                                fEnergyTotalXS[i2]);

  static const G4double conv_factor =  1e-16 * cm * cm;

  sigma *= conv_factor;
  if (sigma == 0.) sigma = 1e-30;
  return (sigma);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNASancheExcitationModel::VibrationEnergy(G4int level)
{
  static G4double energies[9] = { 0.01, 0.024, 0.061, 0.092, 0.204, 0.417, 0.460,
                           0.500, 0.835 };
  return (energies[level] * eV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4DNASancheExcitationModel::RandomSelect(G4double k)
{

  // Level Selection Counting can be done here !

  G4int i = nLevels;
  G4double value = 0.;
  std::deque<G4double> values;

  while (i > 0)
  {
    i--;
    G4double partial = PartialCrossSection(k, i);
    values.push_front(partial);
    value += partial;
  }

  value *= G4UniformRand();

  i = nLevels;

  while (i > 0)
  {
    i--;
    if (values[i] > value)
    {
      //outcount<<i<<"  "<<VibrationEnergy(i)<<G4endl;
      return i;
    }
    value -= values[i];
  }

  //outcount<<0<<"  "<<VibrationEnergy(0)<<G4endl;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNASancheExcitationModel::Sum(G4double k)
{
  G4double totalCrossSection = 0.;

  for (G4int i = 0; i < nLevels; i++)
  {
    totalCrossSection += PartialCrossSection(k, i);
  }

  return totalCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4DNASancheExcitationModel::LinInterpolate(G4double e1,
                                                    G4double e2,
                                                    G4double e,
                                                    G4double xs1,
                                                    G4double xs2)
{
  G4double a = (xs2 - xs1) / (e2 - e1);
  G4double b = xs2 - a * e2;
  G4double value = a * e + b;
  // G4cout<<"interP >>  "<<e1<<"  "<<e2<<"  "<<e<<"  "
  // <<xs1<<"  "<<xs2<<"  "<<a<<"  "<<b<<"  "<<value<<G4endl;

  return value;
}

