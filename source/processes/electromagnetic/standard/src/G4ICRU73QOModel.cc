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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4ICRU73QOModel
//
// Author:        Alexander Bagulya
//
// Creation date: 21.05.2010
//
// Modifications: 
//
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ICRU73QOModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4LossTableManager.hh"
#include "G4AntiProton.hh"
#include "G4DeltaAngle.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ICRU73QOModel::G4ICRU73QOModel(const G4ParticleDefinition* p, const G4String& nam)
  : G4VEmModel(nam),
    particle(nullptr),
    isInitialised(false)
{
  mass = charge = chargeSquare = massRate = ratio = 0.0;
  if(p) { SetParticle(p); }
  SetHighEnergyLimit(10.0*MeV);

  lowestKinEnergy  = 5.0*keV;

  sizeL0 = 67;
  sizeL1 = 22;
  sizeL2 = 14;

  theElectron = G4Electron::Electron();

  for (G4int i = 0; i < 100; ++i)
    {
      indexZ[i] = -1;
    }
  for(G4int i = 0; i < NQOELEM; ++i) 
    {
      if(ZElementAvailable[i] > 0) {
        indexZ[ZElementAvailable[i]] = i;
      }
    }
  fParticleChange = nullptr;
  denEffData = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ICRU73QOModel::Initialise(const G4ParticleDefinition* p,
                                 const G4DataVector&)
{
  if(p != particle) SetParticle(p);

  // always false before the run
  SetDeexcitationFlag(false);

  if(!isInitialised) {
    isInitialised = true;

    if(UseAngularGeneratorFlag() && !GetAngularDistribution()) {
      SetAngularDistribution(new G4DeltaAngle());
    }

    G4String pname = particle->GetParticleName();
    fParticleChange = GetParticleChangeForLoss();
    const G4MaterialTable* mtab = G4Material::GetMaterialTable(); 
    denEffData = (*mtab)[0]->GetIonisation()->GetDensityEffectData();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU73QOModel::ComputeCrossSectionPerElectron(
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double cut,
                                                 G4double maxKinEnergy)
{
  G4double cross = 0.0;
  const G4double tmax      = MaxSecondaryEnergy(p, kineticEnergy);
  const G4double maxEnergy = std::min(tmax, maxKinEnergy);
  const G4double cutEnergy = std::max(cut, lowestKinEnergy*massRate);
  if(cutEnergy < maxEnergy) {

    const G4double energy  = kineticEnergy + mass;
    const G4double energy2 = energy*energy;
    const G4double beta2   = kineticEnergy*(kineticEnergy + 2.0*mass)/energy2;
    cross = (maxEnergy - cutEnergy)/(cutEnergy*maxEnergy) 
      - beta2*G4Log(maxEnergy/cutEnergy)/tmax;

    cross *= CLHEP::twopi_mc2_rcl2*chargeSquare/beta2;
  }
 
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU73QOModel::ComputeCrossSectionPerAtom(
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double Z, G4double,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  G4double cross = Z*ComputeCrossSectionPerElectron
                                         (p,kineticEnergy,cutEnergy,maxEnergy);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU73QOModel::CrossSectionPerVolume(
                                           const G4Material* material,
                                           const G4ParticleDefinition* p,
                                                 G4double kineticEnergy,
                                                 G4double cutEnergy,
                                                 G4double maxEnergy)
{
  G4double eDensity = material->GetElectronDensity();
  G4double cross = eDensity*ComputeCrossSectionPerElectron
                                         (p,kineticEnergy,cutEnergy,maxEnergy);
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU73QOModel::ComputeDEDXPerVolume(const G4Material* material,
                                               const G4ParticleDefinition* p,
                                               G4double kineticEnergy,
                                               G4double cut)
{
  SetParticle(p);
  const G4double tmax  = MaxSecondaryEnergy(p, kineticEnergy);
  const G4double tkin  = kineticEnergy/massRate;
  const G4double cutEnergy = std::max(cut, lowestKinEnergy*massRate);
  G4double dedx  = 0.0;
  if(tkin > lowestKinEnergy) { dedx = DEDX(material, tkin); }
  else { dedx = DEDX(material, lowestKinEnergy)*sqrt(tkin/lowestKinEnergy); }

  if (cutEnergy < tmax) {

    const G4double tau = kineticEnergy/mass;
    const G4double x = cutEnergy/tmax;
    dedx += (G4Log(x)*(tau + 1.)*(tau + 1.)/(tau * (tau + 2.0)) + 1.0 - x) * 
      CLHEP::twopi_mc2_rcl2 *chargeSquare * material->GetElectronDensity();
  }
  dedx = std::max(dedx, 0.0);
  return dedx;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU73QOModel::DEDX(const G4Material* material,
                               G4double kineticEnergy) 
{
  G4double eloss = 0.0;
  const std::size_t numberOfElements = material->GetNumberOfElements();
  const G4double* theAtomicNumDensityVector =
                                 material->GetAtomicNumDensityVector();
  
  // Bragg's rule calculation
  const G4ElementVector* theElementVector =
                           material->GetElementVector() ;
  
  //  loop for the elements in the material
  for (std::size_t i=0; i<numberOfElements; ++i)
    {
      const G4Element* element = (*theElementVector)[i] ;
      eloss += DEDXPerElement(element->GetZasInt(), kineticEnergy)
        * theAtomicNumDensityVector[i] * element->GetZ();
    }      
  return eloss;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU73QOModel::DEDXPerElement(G4int AtomicNumber,
                                         G4double kineticEnergy)
{
  G4int Z = std::min(AtomicNumber, 97);
  G4int nbOfShells = std::max(GetNumberOfShells(Z), 1);

  G4double v = CLHEP::c_light * std::sqrt( 2.0*kineticEnergy/proton_mass_c2 );

  G4double fBetheVelocity = CLHEP::fine_structure_const*CLHEP::c_light/v;

  G4double tau   = kineticEnergy/proton_mass_c2;
  G4double gam   = tau + 1.0;
  G4double bg2   = tau * (tau+2.0);
  G4double beta2 = bg2/(gam*gam);

  G4double l0Term = 0, l1Term = 0, l2Term = 0;
  
  for (G4int nos = 0; nos < nbOfShells; ++nos){
    
    G4double NormalizedEnergy = (2.0*CLHEP::electron_mass_c2*beta2) / 
      GetShellEnergy(Z,nos);

    G4double shStrength = GetShellStrength(Z,nos);

    G4double l0 = GetL0(NormalizedEnergy);
    l0Term += shStrength  * l0; 

    G4double l1 = GetL1(NormalizedEnergy);
    l1Term += shStrength * l1; 

    G4double l2 = GetL2(NormalizedEnergy);
    l2Term += shStrength * l2;

  }
  G4double dedx  = 2*CLHEP::twopi_mc2_rcl2*chargeSquare*factorBethe[Z]*
    (l0Term + charge*fBetheVelocity*l1Term 
     + chargeSquare*fBetheVelocity*fBetheVelocity*l2Term)/beta2;
  return dedx;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU73QOModel::GetOscillatorEnergy(G4int Z,
                                              G4int nbOfTheShell) const
{ 
  G4int idx = denEffData->GetElementIndex(Z, kStateUndefined);
  if(idx == -1) { idx = denEffData->GetElementIndex(Z-1, kStateUndefined); }
  G4double PlasmaEnergy = denEffData->GetPlasmaEnergy(idx);
 
  G4double PlasmaEnergy2 = PlasmaEnergy * PlasmaEnergy;

  G4double plasmonTerm = 0.66667 
    * G4AtomicShells::GetNumberOfElectrons(Z,nbOfTheShell)  
    * PlasmaEnergy2 / (Z*Z) ; 

  static const G4double exphalf = G4Exp(0.5);
  G4double ionTerm = exphalf * 
    (G4AtomicShells::GetBindingEnergy(Z,nbOfTheShell)) ;
  G4double ionTerm2 = ionTerm*ionTerm ;
   
  G4double oscShellEnergy = std::sqrt( ionTerm2 + plasmonTerm );

  return  oscShellEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int G4ICRU73QOModel::GetNumberOfShells(G4int Z) const
{
  G4int nShell = 0;

  if(indexZ[Z] >= 0) { 
    nShell = nbofShellsForElement[indexZ[Z]]; 
  } else { 
    nShell = G4AtomicShells::GetNumberOfShells(Z); 
  }
  return nShell;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU73QOModel::GetShellEnergy(G4int Z, G4int nbOfTheShell) const
{
  G4double shellEnergy = 0.;

  G4int idx = indexZ[Z];

  if(idx >= 0) { 
    shellEnergy = ShellEnergy[startElemIndex[idx] + nbOfTheShell]*CLHEP::eV; 
  } else { 
    shellEnergy = GetOscillatorEnergy(Z, nbOfTheShell); 
  }

  return  shellEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU73QOModel::GetShellStrength(G4int Z, G4int nbOfTheShell) const
{
  G4double shellStrength = 0.;

  G4int idx = indexZ[Z];

  if(idx >= 0) { 
    shellStrength = SubShellOccupation[startElemIndex[idx] + nbOfTheShell] / Z; 
  } else { 
    shellStrength = G4double(G4AtomicShells::GetNumberOfElectrons(Z,nbOfTheShell))/Z; 
  }
  
  return shellStrength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU73QOModel::GetL0(G4double normEnergy) const 
{
  G4int n;
  
  for(n = 0; n < sizeL0; n++) {
    if( normEnergy < L0[n][0] ) break;
  }
  if(0 == n) { n = 1; }
  if(n >= sizeL0) { n = sizeL0 - 1; }

  G4double l0    = L0[n][1];
  G4double l0p   = L0[n-1][1];
  G4double bethe = l0p + (l0 - l0p) * ( normEnergy - L0[n-1][0]) / 
                  (L0[n][0] - L0[n-1][0]);

  return bethe ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU73QOModel::GetL1(G4double normEnergy) const
{
  G4int n;

  for(n = 0; n < sizeL1; n++) {
    if( normEnergy < L1[n][0] ) break;
  }
  if(0 == n) n = 1 ;
  if(n >= sizeL1) n = sizeL1 - 1 ;

  G4double l1    = L1[n][1];
  G4double l1p   = L1[n-1][1];
  G4double barkas= l1p + (l1 - l1p) * ( normEnergy - L1[n-1][0]) / 
                  (L1[n][0] - L1[n-1][0]);

  return barkas;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU73QOModel::GetL2(G4double normEnergy) const
{
  G4int n;
  for(n = 0; n < sizeL2; n++) {
    if( normEnergy < L2[n][0] ) break;
  }
  if(0 == n) n = 1 ;
  if(n >= sizeL2) n = sizeL2 - 1 ;

  G4double l2    = L2[n][1];
  G4double l2p   = L2[n-1][1];
  G4double bloch = l2p + (l2 - l2p) * ( normEnergy - L2[n-1][0]) / 
                  (L2[n][0] - L2[n-1][0]);

  return bloch;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ICRU73QOModel::CorrectionsAlongStep(const G4MaterialCutsCouple*,
                                           const G4DynamicParticle*,
                                           const G4double&,
                                           G4double&)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ICRU73QOModel::SampleSecondaries(std::vector<G4DynamicParticle*>* vdp,
                                        const G4MaterialCutsCouple* couple,
                                        const G4DynamicParticle* dp,
                                        G4double minEnergy,
                                        G4double maxEnergy)
{
  const G4double tmax = MaxSecondaryKinEnergy(dp);
  const G4double xmax = std::min(tmax, maxEnergy);
  const G4double xmin = std::max(minEnergy, lowestKinEnergy*massRate);
  if(xmin >= xmax) { return; }

  G4double kineticEnergy = dp->GetKineticEnergy();
  const G4double energy  = kineticEnergy + mass;
  const G4double energy2 = energy*energy;
  const G4double beta2   = kineticEnergy*(kineticEnergy + 2.0*mass)/energy2;
  G4double grej    = 1.0;
  G4double deltaKinEnergy, f;

  G4ThreeVector direction = dp->GetMomentumDirection();

  // sampling follows ...
  do {
    G4double x = G4UniformRand();
    deltaKinEnergy = xmin*xmax/(xmin*(1.0 - x) + xmax*x);

    f = 1.0 - beta2*deltaKinEnergy/tmax;

    if(f > grej) {
        G4cout << "G4ICRU73QOModel::SampleSecondary Warning! "
               << "Majorant " << grej << " < "
               << f << " for e= " << deltaKinEnergy
               << G4endl;
    }

    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while( grej*G4UniformRand() >= f );

  G4ThreeVector deltaDirection;

  if(UseAngularGeneratorFlag()) {
    const G4Material* mat =  couple->GetMaterial();
    G4int Z = SelectRandomAtomNumber(mat);

    deltaDirection = 
      GetAngularDistribution()->SampleDirection(dp, deltaKinEnergy, Z, mat);

  } else {

    G4double deltaMomentum =
           sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
    G4double totMomentum = energy*sqrt(beta2);
    G4double cost = deltaKinEnergy * (energy + electron_mass_c2) /
                                   (deltaMomentum * totMomentum);
    if(cost > 1.0) { cost = 1.0; }
    G4double sint = sqrt((1.0 - cost)*(1.0 + cost));

    G4double phi = twopi * G4UniformRand() ;

    deltaDirection.set(sint*cos(phi),sint*sin(phi), cost) ;
    deltaDirection.rotateUz(direction);
  }
  // create G4DynamicParticle object for delta ray
  auto delta = new G4DynamicParticle(theElectron,deltaDirection,deltaKinEnergy);

  // Change kinematics of primary particle
  kineticEnergy -= deltaKinEnergy;
  G4ThreeVector finalP = dp->GetMomentum() - delta->GetMomentum();
  finalP               = finalP.unit();
  
  fParticleChange->SetProposedKineticEnergy(kineticEnergy);
  fParticleChange->SetProposedMomentumDirection(finalP);

  vdp->push_back(delta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4ICRU73QOModel::MaxSecondaryEnergy(const G4ParticleDefinition* pd,
                                             G4double kinEnergy)
{
  if(pd != particle) { SetParticle(pd); }
  G4double tau  = kinEnergy/mass;
  G4double tmax = 2.0*electron_mass_c2*tau*(tau + 2.) /
                  (1. + 2.0*(tau + 1.)*ratio + ratio*ratio);
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4int G4ICRU73QOModel::ZElementAvailable[NQOELEM] = 
  {1,2,4,6,7,8,10,13,14,-18,
   22,26,28,29,32,36,42,47,
   50,54,73,74,78,79,82,92};

const G4int G4ICRU73QOModel::nbofShellsForElement[NQOELEM] = 
  {1,1,2,3,3,3,3,4,5,4,
   5,5,5,5,6,4,6,6,
   7,6,6,8,7,7,9,9};

const G4int G4ICRU73QOModel::startElemIndex[NQOELEM] = 
  {0,1,2,4,7,10,13,16,20,25,
   29,34,39,44,49,55,59,65,
   71,78,84,90,98,105,112,121};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// SubShellOccupation = Z * ShellStrength
const G4double G4ICRU73QOModel::SubShellOccupation[NQODATA] = 
  {
    1.000, // H 0
    2.000, // He 1
    1.930, 2.070, // Be 2-3
    1.992, 1.841, 2.167, // C 4-6    
    1.741, 1.680, 3.579, // N 7-9
    1.802, 1.849, 4.349, // O 10-12
    1.788, 2.028, 6.184, // Ne 13-15
    1.623, 2.147, 6.259, 2.971, // Al 16-19
    1.631, 2.094, 6.588, 2.041, 1.646, // Si 20-24
    1.535, 8.655, 1.706, 6.104, // Ar 25-28
    1.581, 8.358, 8.183, 2.000, 1.878, // Ti 29-33
    1.516, 8.325, 8.461, 6.579, 1.119, // Fe 34-38
    1.422, 7.81, 8.385, 8.216, 2.167, // Ni 39-43
    1.458, 8.049, 8.79, 9.695, 1.008, // Cu 44-48
    1.442, 7.791, 7.837, 10.122, 2.463, 2.345, // Ge 49-54
    1.645, 7.765, 19.192, 7.398, // Kr 55-58
    1.313, 6.409, 19.229, 8.633, 5.036, 1.380, // Mo 59-64
    1.295, 6.219, 18.751, 8.748, 10.184, 1.803, // Ag 65-70
    1.277, 6.099, 20.386, 8.011, 10.007, 2.272, 1.948, // Sn 71-77
    1.563, 6.312, 21.868, 5.762, 11.245, 7.250, // Xe 78-83
    0.9198, 6.5408, 18.9727, 24.9149, 15.0161, 6.6284, // Ta 84-89
    1.202, 5.582, 19.527, 18.741, 8.411, 14.387, 4.042, 2.108, // W 90-97
    1.159, 5.467, 18.802, 33.905, 8.300, 9.342, 1.025, // Pt 98-104
    1.124, 5.331, 18.078, 34.604, 8.127, 10.414, 1.322, // Au 105-111
    2.000, 8.000, 18.000, 18.000, 14.000, 8.000, 10.000, 2.000, 2.000, // Pb 112-120
    2.000, 8.000, 18.000, 32.000, 18.000, 8.000, 2.000, 1.000, 3.000  // U 121-129
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// ShellEnergy in eV
const G4double G4ICRU73QOModel::ShellEnergy[NQODATA] = 
  {
    19.2, // H
    41.8, // He
    209.11, 21.68, // Be
    486.2, 60.95, 23.43, // C
    732.61, 100.646, 23.550, // N    
    965.1, 129.85, 31.60, // O
    1525.9, 234.9, 56.18, // Ne
    2701, 476.5, 150.42, 16.89, // Al
    3206.1, 586.4, 186.8, 23.52, 14.91, // Si
    5551.6, 472.43, 124.85, 22.332, // Ar
    8554.6, 850.58, 93.47, 39.19, 19.46, // Ti
    12254.7, 1279.29, 200.35, 49.19, 17.66, // Fe
    14346.9, 1532.28, 262.71, 74.37, 23.03, // Ni
    15438.5, 1667.96, 294.1, 70.69, 16.447, // Cu
    19022.1, 2150.79, 455.79, 179.87, 57.89, 20.95, // Ge
    24643, 2906.4, 366.85, 22.24, // Kr
    34394, 4365.3, 589.36, 129.42, 35.59, 18.42, // Mo
    43664.3, 5824.91, 909.79, 175.47, 54.89, 19.63, // Ag
    49948, 6818.2, 1036.1, 172.65, 70.89, 33.87, 14.54, // Sn
    58987, 8159, 1296.6, 356.75, 101.03, 16.52, // Xe
    88926, 18012, 3210, 575, 108.7, 30.8, // Ta
    115025.9, 17827.44, 3214.36, 750.41, 305.21, 105.50, 38.09, 21.25, // W
    128342, 20254, 3601.8, 608.1, 115.0, 42.75, 17.04, // Pt
    131872, 20903, 3757.4, 682.1, 105.2, 44.89, 17.575, // Au
    154449, 25067, 5105.0, 987.44, 247.59, 188.1, 40.61, 19.2, 15.17, // Pb
    167282, 27868, 6022.7, 1020.4, 244.81, 51.33, 13, 11.06, 14.43  // U
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Data for L0 from: Sigmund P., Haagerup U. Phys. Rev. A34 (1986) 892-910
const G4double G4ICRU73QOModel::L0[67][2] =
{
  {0.00,        0.000001},
  {0.10,        0.000001},
  {0.12,        0.00001},
  {0.14,        0.00005},
  {0.16,        0.00014},
  {0.18,        0.00030},
  {0.20,        0.00057},
  {0.25,        0.00189},
  {0.30,        0.00429},
  {0.35,        0.00784},
  {0.40,        0.01248},
  {0.45,        0.01811},
  {0.50,        0.02462},
  {0.60,        0.03980},
  {0.70,        0.05731},
  {0.80,        0.07662},
  {0.90,        0.09733},
  {1.00,        0.11916},
  {1.20,        0.16532},
  {1.40,        0.21376},
  {1.60,        0.26362},
  {1.80,        0.31428},
  {2.00,        0.36532},
  {2.50,        0.49272},
  {3.00,        0.61765},
  {3.50,        0.73863},
  {4.00,        0.85496},
  {4.50,        0.96634},
  {5.00,        1.07272},
  {6.00,        1.27086},
  {7.00,        1.45075},
  {8.00,        1.61412},
  {9.00,        1.76277},
  {10.00,       1.89836},
  {12.00,       2.13625},
  {14.00,       2.33787},
  {16.00,       2.51093},
  {18.00,       2.66134},
  {20.00,       2.79358},
  {25.00,       3.06539},
  {30.00,       3.27902},
  {35.00,       3.45430},
  {40.00,       3.60281},
  {45.00,       3.73167},
  {50.00,       3.84555},
  {60.00,       4.04011},
  {70.00,       4.20264},
  {80.00,       4.34229},
  {90.00,       4.46474},
  {100.00,      4.57378},
  {120.00,      4.76155},
  {140.00,      4.91953},
  {160.00,      5.05590},
  {180.00,      5.17588},
  {200.00,      5.28299},
  {250.00,      5.50925},
  {300.00,      5.69364},
  {350.00,      5.84926},
  {400.00,      5.98388},
  {450.00,      6.10252},
  {500.00,      6.20856},
  {600.00,      6.39189},
  {700.00,      6.54677},
  {800.00,      6.68084},
  {900.00,      6.79905},
  {1000.00,     6.90474}
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Data for L1 from: Mikkelsen H.H., Sigmund P. Phys. Rev. A40 (1989) 101-116
const G4double G4ICRU73QOModel::L1[22][2] =
{
  {0.00,       -0.000001},
  {0.10,       -0.00001},
  {0.20,       -0.00049},
  {0.30,       -0.00084},
  {0.40,        0.00085},
  {0.50,        0.00519},
  {0.60,        0.01198},
  {0.70,        0.02074},
  {0.80,        0.03133},
  {0.90,        0.04369},
  {1.00,        0.06035},
  {2.00,        0.24023},
  {3.00,        0.44284},
  {4.00,        0.62012},
  {5.00,        0.77031},
  {6.00,        0.90390},
  {7.00,        1.02705},
  {8.00,        1.10867},
  {9.00,        1.17546},
  {10.00,       1.21599},
  {15.00,        1.24349},
  {20.00,        1.16752}
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Data for L2 from: Mikkelsen H.H. Nucl. Instr. Meth. B58 (1991) 136-148
const G4double G4ICRU73QOModel::L2[14][2] =
{
  {0.00,        0.000001},
  {0.10,        0.00001},
  {0.20,        0.00000},
  {0.40,       -0.00120},
  {0.60,       -0.00036},
  {0.80,        0.00372},
  {1.00,        0.01298},
  {2.00,        0.08296},
  {4.00,        0.21953},
  {6.00,        0.23903},
  {8.00,        0.20893},
  {10.00,       0.10879},
  {20.00,      -0.88409},          
  {40.00,      -1.13902}
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Correction obtained by V.Ivanchenko using G4BetheBlochModel
const G4double G4ICRU73QOModel::factorBethe[99] = { 1.0, 
0.9637, 0.9872, 0.9469, 0.9875, 0.91, 0.989, 0.9507, 0.9773, 0.8621, 0.979,   // 1 - 10
0.8357, 0.868, 0.9417, 0.9466, 0.8911, 0.905, 0.944, 0.9607, 0.928, 0.96,   // 11 - 20
0.9098, 0.976, 0.8425, 0.8099, 0.7858, 0.947, 0.7248, 0.9106, 0.9246, 0.6821,   // 21 - 30
0.7223, 0.9784, 0.774, 0.7953, 0.829, 0.9405, 0.8318, 0.8583, 0.8563, 0.8481,   // 31 - 40
0.8207, 0.9033, 0.8063, 0.7837, 0.7818, 0.744, 0.875, 0.7693, 0.7871, 0.8459,   // 41 - 50
0.8231, 0.8462, 0.853, 0.8736, 0.856, 0.8762, 0.8629, 0.8323, 0.8064, 0.7828,   // 51 - 60
0.7533, 0.7273, 0.7093, 0.7157, 0.6823, 0.6612, 0.6418, 0.6395, 0.6323, 0.6221,   // 61 - 70
0.6497, 0.6746, 0.8568, 0.8541, 0.6958, 0.6962, 0.7051, 0.863, 0.8588, 0.7226,   // 71 - 80
0.7454, 0.78, 0.7783, 0.7996, 0.8216, 0.8632, 0.8558, 0.8792, 0.8745, 0.8676,   // 81 - 90
0.8321, 0.8272, 0.7999, 0.7934, 0.7787, 0.7851, 0.7692, 0.7598}; 
