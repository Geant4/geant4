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
// File name:     G4BetheHeitlerModel
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 15.03.2005
//
// Modifications by Vladimir Ivanchenko, Michel Maire, Mihaly Novak
//
// Class Description:
//
// -------------------------------------------------------------------
//

#include "G4BetheHeitlerModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4Pow.hh"
#include "G4Exp.hh"
#include "G4ModifiedTsai.hh"

const G4int G4BetheHeitlerModel::gMaxZet = 120; 
std::vector<G4BetheHeitlerModel::ElementData*> G4BetheHeitlerModel::gElementData;

G4BetheHeitlerModel::G4BetheHeitlerModel(const G4ParticleDefinition*, 
                                         const G4String& nam)
: G4VEmModel(nam), 
  fG4Calc(G4Pow::GetInstance()), fTheGamma(G4Gamma::Gamma()),
  fTheElectron(G4Electron::Electron()), fThePositron(G4Positron::Positron()),
  fParticleChange(nullptr) 
{
  SetAngularDistribution(new G4ModifiedTsai());
}

G4BetheHeitlerModel::~G4BetheHeitlerModel()
{
  if (IsMaster()) {
    // clear ElementData container
    for (std::size_t iz = 0; iz < gElementData.size(); ++iz) {
      if (gElementData[iz]) delete gElementData[iz];
    }
    gElementData.clear(); 
  }
}

void G4BetheHeitlerModel::Initialise(const G4ParticleDefinition* p, 
                                     const G4DataVector& cuts)
{
  if (IsMaster()) {
    InitialiseElementData();
  }
  if (!fParticleChange) { fParticleChange = GetParticleChangeForGamma(); }
  if (IsMaster()) { 
    InitialiseElementSelectors(p, cuts); 
  }
}

void G4BetheHeitlerModel::InitialiseLocal(const G4ParticleDefinition*, 
                                          G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}

// Calculates the microscopic cross section in GEANT4 internal units.
// A parametrized formula from L. Urban is used to estimate
// the total cross section.
// It gives a good description of the data from 1.5 MeV to 100 GeV.
// below 1.5 MeV: sigma=sigma(1.5MeV)*(GammaEnergy-2electronmass)
//                                   *(GammaEnergy-2electronmass) 
G4double 
G4BetheHeitlerModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*, 
                                                G4double gammaEnergy, G4double Z, 
                                                G4double, G4double, G4double)
{
  G4double xSection = 0.0 ;
  // short versions
  static const G4double kMC2  = CLHEP::electron_mass_c2;
  // zero cross section below the kinematical limit: Eg<2mc^2
  if (Z < 0.9 || gammaEnergy <= 2.0*kMC2) { return xSection; }
  //
  static const G4double gammaEnergyLimit = 1.5*CLHEP::MeV;
  // set coefficients a, b c
  static const G4double a0 =  8.7842e+2*CLHEP::microbarn;
  static const G4double a1 = -1.9625e+3*CLHEP::microbarn; 
  static const G4double a2 =  1.2949e+3*CLHEP::microbarn;
  static const G4double a3 = -2.0028e+2*CLHEP::microbarn; 
  static const G4double a4 =  1.2575e+1*CLHEP::microbarn; 
  static const G4double a5 = -2.8333e-1*CLHEP::microbarn;
  
  static const G4double b0 = -1.0342e+1*CLHEP::microbarn;
  static const G4double b1 =  1.7692e+1*CLHEP::microbarn;
  static const G4double b2 = -8.2381   *CLHEP::microbarn;
  static const G4double b3 =  1.3063   *CLHEP::microbarn;
  static const G4double b4 = -9.0815e-2*CLHEP::microbarn;
  static const G4double b5 =  2.3586e-3*CLHEP::microbarn;
  
  static const G4double c0 = -4.5263e+2*CLHEP::microbarn;
  static const G4double c1 =  1.1161e+3*CLHEP::microbarn; 
  static const G4double c2 = -8.6749e+2*CLHEP::microbarn;
  static const G4double c3 =  2.1773e+2*CLHEP::microbarn; 
  static const G4double c4 = -2.0467e+1*CLHEP::microbarn;
  static const G4double c5 =  6.5372e-1*CLHEP::microbarn;
  // check low energy limit of the approximation (1.5 MeV)
  G4double gammaEnergyOrg = gammaEnergy;
  if (gammaEnergy < gammaEnergyLimit) { gammaEnergy = gammaEnergyLimit; }
  // compute gamma energy variables
  const G4double x  = G4Log(gammaEnergy/kMC2);
  const G4double x2 = x *x; 
  const G4double x3 = x2*x;
  const G4double x4 = x3*x;
  const G4double x5 = x4*x;
  //
  const G4double F1 = a0 + a1*x + a2*x2 + a3*x3 + a4*x4 + a5*x5;
  const G4double F2 = b0 + b1*x + b2*x2 + b3*x3 + b4*x4 + b5*x5;
  const G4double F3 = c0 + c1*x + c2*x2 + c3*x3 + c4*x4 + c5*x5;     
  // compute the approximated cross section 
  xSection = (Z + 1.)*(F1*Z + F2*Z*Z + F3);
  // check if we are below the limit of the approximation and apply correction
  if (gammaEnergyOrg < gammaEnergyLimit) {
    const G4double dum = (gammaEnergyOrg-2.*kMC2)/(gammaEnergyLimit-2.*kMC2);
    xSection *= dum*dum;
  }
  // make sure that the cross section is never negative
  xSection = std::max(xSection, 0.); 
  return xSection;
}

// The secondaries e+e- energies are sampled using the Bethe - Heitler
// cross sections with Coulomb correction.
// A modified version of the random number techniques of Butcher & Messel
// is used (Nuc Phys 20(1960),15).
//
// GEANT4 internal units.
//
// Note 1 : Effects due to the breakdown of the Born approximation at
//          low energy are ignored.
// Note 2 : The differential cross section implicitly takes account of 
//          pair creation in both nuclear and atomic electron fields.
//          However triplet prodution is not generated.
void G4BetheHeitlerModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                            const G4MaterialCutsCouple* couple,
                                            const G4DynamicParticle* aDynamicGamma,
                                            G4double, G4double)
{
  // set some constant values
  const G4double    gammaEnergy = aDynamicGamma->GetKineticEnergy();
  const G4double    eps0        = CLHEP::electron_mass_c2/gammaEnergy;
  //
  // check kinematical limit: gamma energy(Eg) must be at least 2 e- rest mass
  if (eps0 > 0.5) { return; }
  //
  // select target element of the material (probs. are based on partial x-secs)
  const G4Element* anElement = SelectTargetAtom(couple, fTheGamma, gammaEnergy,
                                          aDynamicGamma->GetLogKineticEnergy());

  // 
  // get the random engine
  CLHEP::HepRandomEngine* rndmEngine = G4Random::getTheEngine();
  //
  // 'eps' is the total energy transferred to one of the e-/e+ pair in initial
  // gamma energy units Eg. Since the corresponding DCS is symmetric on eps=0.5,
  // the kinematical limits for eps0=mc^2/Eg <= eps <= 0.5 
  // 1. 'eps' is sampled uniformly on the [eps0, 0.5] inteval if Eg<Egsmall 
  // 2. otherwise, on the [eps_min, 0.5] interval according to the DCS (case 2.) 
  G4double eps;
  // case 1.
  static const G4double Egsmall = 2.*CLHEP::MeV;
  if (gammaEnergy < Egsmall) {
    eps = eps0 + (0.5-eps0)*rndmEngine->flat();
  } else {
  // case 2.
    // get the Coulomb factor for the target element (Z) and gamma energy (Eg)
    // F(Z) = 8*ln(Z)/3           if Eg <= 50 [MeV] => no Coulomb correction
    // F(Z) = 8*ln(Z)/3 + 8*fc(Z) if Eg  > 50 [MeV] => fc(Z) is the Coulomb cor.
    //
    // The screening variable 'delta(eps)' = 136*Z^{-1/3}*eps0/[eps(1-eps)]
    // Due to the Coulomb correction, the DCS can go below zero even at 
    // kinematicaly allowed eps > eps0 values. In order to exclude this eps 
    // range with negative DCS, the minimum eps value will be set to eps_min = 
    // max[eps0, epsp] with epsp is the solution of SF(delta(epsp)) - F(Z)/2 = 0 
    // with SF being the screening function (SF1=SF2 at high value of delta). 
    // The solution is epsp = 0.5 - 0.5*sqrt[ 1 - 4*136*Z^{-1/3}eps0/deltap] 
    // with deltap = Exp[(42.038-F(Z))/8.29]-0.958. So the limits are:
    // - when eps=eps_max = 0.5            => delta_min = 136*Z^{-1/3}*eps0/4
    // - epsp = 0.5 - 0.5*sqrt[ 1 - delta_min/deltap]
    // - and eps_min = max[eps0, epsp]  
    static const G4double midEnergy = 50.*CLHEP::MeV;
    const  G4int           iZet = std::min(gMaxZet, anElement->GetZasInt());   
    const  G4double deltaFactor = 136.*eps0/anElement->GetIonisation()->GetZ3();
    G4double           deltaMax = gElementData[iZet]->fDeltaMaxLow;
    G4double                 FZ = 8.*anElement->GetIonisation()->GetlogZ3();
    if (gammaEnergy > midEnergy) { 
      FZ      += 8.*(anElement->GetfCoulomb()); 
      deltaMax = gElementData[iZet]->fDeltaMaxHigh;
    }
    const G4double deltaMin = 4.*deltaFactor; 
    // 
    // compute the limits of eps
    const G4double epsp     = 0.5 - 0.5*std::sqrt(1. - deltaMin/deltaMax) ;
    const G4double epsMin   = std::max(eps0,epsp);
    const G4double epsRange = 0.5 - epsMin;
    //
    // sample the energy rate (eps) of the created electron (or positron)
    G4double F10, F20;
    ScreenFunction12(deltaMin, F10, F20); 
    F10 -= FZ;
    F20 -= FZ; 
    const G4double NormF1   = std::max(F10 * epsRange * epsRange, 0.); 
    const G4double NormF2   = std::max(1.5 * F20                , 0.);
    const G4double NormCond = NormF1/(NormF1 + NormF2); 
    // we will need 3 uniform random number for each trial of sampling 
    G4double rndmv[3];
    G4double greject = 0.;
    do {
      rndmEngine->flatArray(3, rndmv);
      if (NormCond > rndmv[0]) {
        eps = 0.5 - epsRange * fG4Calc->A13(rndmv[1]);
        const G4double delta = deltaFactor/(eps*(1.-eps));
        greject = (ScreenFunction1(delta)-FZ)/F10;
      } else { 
        eps = epsMin + epsRange*rndmv[1];
        const G4double delta = deltaFactor/(eps*(1.-eps));
        greject = (ScreenFunction2(delta)-FZ)/F20;
      }
      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while (greject < rndmv[2]);
  } //  end of eps sampling
  //
  // select charges randomly
  G4double eTotEnergy, pTotEnergy;
  if (rndmEngine->flat() > 0.5) {
    eTotEnergy = (1.-eps)*gammaEnergy;
    pTotEnergy = eps*gammaEnergy; 
  } else {
    pTotEnergy = (1.-eps)*gammaEnergy;
    eTotEnergy = eps*gammaEnergy;
  }
  //
  // sample pair kinematics
  const G4double eKinEnergy = std::max(0.,eTotEnergy - CLHEP::electron_mass_c2);
  const G4double pKinEnergy = std::max(0.,pTotEnergy - CLHEP::electron_mass_c2);
  //
  G4ThreeVector eDirection, pDirection;
  //
  GetAngularDistribution()->SamplePairDirections(aDynamicGamma, 
                                                 eKinEnergy, pKinEnergy,
                                                 eDirection, pDirection);
  // create G4DynamicParticle object for the particle1
  auto aParticle1= new G4DynamicParticle(fTheElectron,eDirection,eKinEnergy);
  // create G4DynamicParticle object for the particle2
  auto aParticle2= new G4DynamicParticle(fThePositron,pDirection,pKinEnergy);
  // Fill output vector
  fvect->push_back(aParticle1);
  fvect->push_back(aParticle2);
  // kill incident photon
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);   
}

// should be called only by the master and at initialisation
void G4BetheHeitlerModel::InitialiseElementData() 
{
  G4int size = (G4int)gElementData.size();
  if (size < gMaxZet+1) {
    gElementData.resize(gMaxZet+1, nullptr);
  }
  // create for all elements that are in the detector
  const G4ElementTable* elemTable = G4Element::GetElementTable();
  std::size_t numElems = (*elemTable).size();
  for (std::size_t ie = 0; ie < numElems; ++ie) {
    const G4Element* elem = (*elemTable)[ie];
    const G4int        iz = std::min(gMaxZet, elem->GetZasInt());
    if (!gElementData[iz]) { // create it if doesn't exist yet
      G4double FZLow     = 8.*elem->GetIonisation()->GetlogZ3();
      G4double FZHigh    = FZLow + 8.*elem->GetfCoulomb();
      auto elD           = new ElementData();
      elD->fDeltaMaxLow  = G4Exp((42.038 - FZLow )/8.29) - 0.958;
      elD->fDeltaMaxHigh = G4Exp((42.038 - FZHigh)/8.29) - 0.958;
      gElementData[iz]   = elD;
    }
  }
}

