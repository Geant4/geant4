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
// File name:     G4KleinNishinaModel
//
// Author:        Vladimir Ivanchenko on base of G4KleinNishinaCompton
//
// Creation date: 13.06.2010
//
// Modifications:
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4KleinNishinaModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4AtomicShells.hh"
#include "G4LossTableManager.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4KleinNishinaModel::G4KleinNishinaModel(const G4String& nam)
  : G4VEmModel(nam), 
    lv1(0.,0.,0.,0.),
    lv2(0.,0.,0.,0.),
    bst(0.,0.,0.)
{
  theGamma = G4Gamma::Gamma();
  theElectron = G4Electron::Electron();
  lowestSecondaryEnergy = 10*eV;
  limitFactor       = 4;
  fProbabilities.resize(9,0.0);
  SetDeexcitationFlag(true);
  fParticleChange = nullptr;
  fAtomDeexcitation = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4KleinNishinaModel::~G4KleinNishinaModel() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4KleinNishinaModel::Initialise(const G4ParticleDefinition* p,
                                     const G4DataVector& cuts)
{
  fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();
  if(IsMaster()) { InitialiseElementSelectors(p, cuts); }
  if(nullptr == fParticleChange) { 
    fParticleChange = GetParticleChangeForGamma(); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4KleinNishinaModel::InitialiseLocal(const G4ParticleDefinition*,
                                          G4VEmModel* masterModel)
{
  SetElementSelectors(masterModel->GetElementSelectors());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4KleinNishinaModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                                G4double gammaEnergy,
                                                G4double Z, G4double,
                                                G4double, G4double)
{
  G4double xSection = 0.0 ;
  if (gammaEnergy <= LowEnergyLimit()) { return xSection; }

  static const G4double a = 20.0 , b = 230.0 , c = 440.0;

static const G4double
  d1= 2.7965e-1*CLHEP::barn, d2=-1.8300e-1*CLHEP::barn, 
  d3= 6.7527   *CLHEP::barn, d4=-1.9798e+1*CLHEP::barn,
  e1= 1.9756e-5*CLHEP::barn, e2=-1.0205e-2*CLHEP::barn, 
  e3=-7.3913e-2*CLHEP::barn, e4= 2.7079e-2*CLHEP::barn,
  f1=-3.9178e-7*CLHEP::barn, f2= 6.8241e-5*CLHEP::barn, 
  f3= 6.0480e-5*CLHEP::barn, f4= 3.0274e-4*CLHEP::barn;
  
  G4double p1Z = Z*(d1 + e1*Z + f1*Z*Z), p2Z = Z*(d2 + e2*Z + f2*Z*Z),
           p3Z = Z*(d3 + e3*Z + f3*Z*Z), p4Z = Z*(d4 + e4*Z + f4*Z*Z);

  G4double T0  = 15.0*keV; 
  if (Z < 1.5) { T0 = 40.0*keV; } 

  G4double X   = max(gammaEnergy, T0) / electron_mass_c2;
  xSection = p1Z*G4Log(1.+2.*X)/X
               + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
                
  //  modification for low energy. (special case for Hydrogen)
  static const G4double dT0 = keV;
  if (gammaEnergy < T0) {
    X = (T0+dT0) / electron_mass_c2 ;
    G4double sigma = p1Z*G4Log(1.+2*X)/X
                    + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
    G4double   c1 = -T0*(sigma-xSection)/(xSection*dT0);             
    G4double   c2 = 0.150; 
    if (Z > 1.5) { c2 = 0.375-0.0556*G4Log(Z); }
    G4double    y = G4Log(gammaEnergy/T0);
    xSection *= G4Exp(-y*(c1+c2*y));          
  }

  if(xSection < 0.0) { xSection = 0.0; }
  //  G4cout << "e= " << GammaEnergy << " Z= " << Z 
  //  << " cross= " << xSection << G4endl;
  return xSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4KleinNishinaModel::SampleSecondaries(
                             std::vector<G4DynamicParticle*>* fvect,
                             const G4MaterialCutsCouple* couple,
                             const G4DynamicParticle* aDynamicGamma,
                             G4double,
                             G4double)
{
  // primary gamma
  G4double energy = aDynamicGamma->GetKineticEnergy();

  // do nothing below the threshold
  if(energy <= LowEnergyLimit()) { return; }

  G4ThreeVector direction = aDynamicGamma->GetMomentumDirection();

  // select atom
  const G4Element* elm = SelectRandomAtom(couple, theGamma, energy);

  // select shell first
  G4int nShells = elm->GetNbOfAtomicShells();
  if(nShells > (G4int)fProbabilities.size()) { fProbabilities.resize(nShells); }
  G4double totprob = 0.0;
  G4int i;
  for(i=0; i<nShells; ++i) {
    //G4double bindingEnergy = elm->GetAtomicShell(i);
    totprob += elm->GetNbOfShellElectrons(i);
    //totprob += elm->GetNbOfShellElectrons(i)/(bindingEnergy*bindingEnergy);
    fProbabilities[i] = totprob; 
  }

  // Loop on sampling
  static const G4int nlooplim = 1000;
  G4int nloop = 0;

  G4double bindingEnergy, ePotEnergy, eKinEnergy;
  G4double gamEnergy0, gamEnergy1;

  CLHEP::HepRandomEngine* rndmEngineMod = G4Random::getTheEngine();
  G4double rndm[4];

  do {
    ++nloop;

    // 4 random numbers to select e-
    rndmEngineMod->flatArray(4, rndm);
    G4double xprob = totprob*rndm[0];

    // select shell
    for(i=0; i<nShells; ++i) { if(xprob <= fProbabilities[i]) { break; } }
   
    bindingEnergy = elm->GetAtomicShell(i);
    lv1.set(0.0,0.0,energy,energy);
    /*
    G4cout << "nShells= " << nShells << " i= " << i 
       << " Egamma= " << energy << " Ebind= " << bindingEnergy
       << G4endl;
    */
    // for rest frame of the electron
    G4double x = -G4Log(rndm[1]);
    eKinEnergy = bindingEnergy*x;
    ePotEnergy = bindingEnergy*(1.0 + x);

    // for rest frame of the electron
    G4double eTotMomentum = sqrt(eKinEnergy*(eKinEnergy + 2*electron_mass_c2));
    G4double phi = rndm[2]*twopi;
    G4double costet = 2*rndm[3] - 1;
    G4double sintet = sqrt((1 - costet)*(1 + costet));
    lv2.set(eTotMomentum*sintet*cos(phi),eTotMomentum*sintet*sin(phi),
            eTotMomentum*costet,eKinEnergy + electron_mass_c2);
    bst = lv2.boostVector();
    lv1.boost(-bst);

    gamEnergy0 = lv1.e();
   
    // In the rest frame of the electron
    // The scattered gamma energy is sampled according to Klein-Nishina formula
    // The random number techniques of Butcher & Messel are used 
    // (Nuc Phys 20(1960),15). 
    G4double E0_m = gamEnergy0/electron_mass_c2;

    //G4cout << "Nloop= "<< nloop << " Ecm(keV)= " << gamEnergy0/keV << G4endl;
    //
    // sample the energy rate of the scattered gamma 
    //

    G4double epsilon, epsilonsq, onecost, sint2, greject ;

    G4double eps0       = 1./(1 + 2*E0_m);
    G4double epsilon0sq = eps0*eps0;
    G4double alpha1     = - G4Log(eps0);
    G4double alpha2     = alpha1 + 0.5*(1 - epsilon0sq);

    do {
      ++nloop;
      // false interaction if too many iterations
      if(nloop > nlooplim) { return; }

      // 3 random numbers to sample scattering
      rndmEngineMod->flatArray(3, rndm);

      if ( alpha1 > alpha2*rndm[0] ) {
        epsilon   = G4Exp(-alpha1*rndm[1]);   // epsilon0**r
        epsilonsq = epsilon*epsilon; 

      } else {
        epsilonsq = epsilon0sq + (1.- epsilon0sq)*rndm[1];
        epsilon   = sqrt(epsilonsq);
      }

      onecost = (1.- epsilon)/(epsilon*E0_m);
      sint2   = onecost*(2.-onecost);
      greject = 1. - epsilon*sint2/(1.+ epsilonsq);

      // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    } while (greject < rndm[2]);
    gamEnergy1 = epsilon*gamEnergy0;
 
    // before scattering total 4-momentum in e- system
    lv2.set(0.0,0.0,0.0,electron_mass_c2);
    lv2 += lv1;
 
    //
    // scattered gamma angles. ( Z - axis along the parent gamma)
    //
    if(sint2 < 0.0) { sint2 = 0.0; }
    costet = 1. - onecost; 
    sintet = sqrt(sint2);
    phi  = twopi * rndmEngineMod->flat();

    // e- recoil
    //
    // in  rest frame of the electron
    G4ThreeVector gamDir = lv1.vect().unit();
    G4ThreeVector v = G4ThreeVector(sintet*cos(phi),sintet*sin(phi),costet);
    v.rotateUz(gamDir);
    lv1.set(gamEnergy1*v.x(),gamEnergy1*v.y(),gamEnergy1*v.z(),gamEnergy1);
    lv2 -= lv1;
    //G4cout<<"Egam(keV)= " << lv1.e()/keV
    //          <<" Ee(keV)= " << (lv2.e()-electron_mass_c2)/keV << G4endl;
    lv2.boost(bst);
    eKinEnergy = lv2.e() - electron_mass_c2 - ePotEnergy;   
    //G4cout << "Nloop= " << nloop << " eKinEnergy= " << eKinEnergy << G4endl;

    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  } while ( eKinEnergy < 0.0 );

  //
  // update G4VParticleChange for the scattered gamma
  //
   
  lv1.boost(bst);
  gamEnergy1 = lv1.e();
  if(gamEnergy1 > lowestSecondaryEnergy) {
    G4ThreeVector gamDirection1 = lv1.vect().unit();
    gamDirection1.rotateUz(direction);
    fParticleChange->ProposeMomentumDirection(gamDirection1);
  } else { 
    fParticleChange->ProposeTrackStatus(fStopAndKill);
    gamEnergy1 = 0.0;
  }
  fParticleChange->SetProposedKineticEnergy(gamEnergy1);

  //
  // kinematic of the scattered electron
  //

  if(eKinEnergy > lowestSecondaryEnergy) {
    G4ThreeVector eDirection = lv2.vect().unit();
    eDirection.rotateUz(direction);
    auto dp = new G4DynamicParticle(theElectron,eDirection,eKinEnergy);
    fvect->push_back(dp);
  } else { eKinEnergy = 0.0; }

  G4double edep = energy - gamEnergy1 - eKinEnergy;
  G4double esec = 0.0;
  
  // sample deexcitation
  //
  if(nullptr != fAtomDeexcitation) {
    G4int index = couple->GetIndex();
    if(fAtomDeexcitation->CheckDeexcitationActiveRegion(index)) {
      G4int Z = elm->GetZasInt();
      auto as = (G4AtomicShellEnumerator)(i);
      const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
      G4int nbefore = (G4int)fvect->size();
      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, index);
      G4int nafter = (G4int)fvect->size();
      //G4cout << "N1= " << nbefore << "  N2= " << nafter << G4endl;
      for (G4int j=nbefore; j<nafter; ++j) {
        G4double e = ((*fvect)[j])->GetKineticEnergy();
        if(esec + e > edep) {
          // correct energy in order to have energy balance
          e = edep - esec;
          ((*fvect)[j])->SetKineticEnergy(e);
          esec += e;
          /*            
            G4cout << "### G4KleinNishinaModel Edep(eV)= " << edep/eV 
                   << " Esec(eV)= " << esec/eV 
                   << " E["<< j << "](eV)= " << e/eV
                   << " N= " << nafter
                   << " Z= " << Z << " shell= " << i 
                   << "  Ebind(keV)= " << bindingEnergy/keV 
                   << "  Eshell(keV)= " << shell->BindingEnergy()/keV 
                   << G4endl;
          */
          // delete the rest of secondaries (should not happens)
          for (G4int jj=nafter-1; jj>j; --jj) { 
            delete (*fvect)[jj]; 
            fvect->pop_back(); 
          }
          break;              
        }
        esec += e; 
      }
      edep -= esec;
    }
  }
  if(std::abs(energy - gamEnergy1 - eKinEnergy - esec - edep) > eV) {
    G4cout << "### G4KleinNishinaModel dE(eV)= " 
           << (energy - gamEnergy1 - eKinEnergy - esec - edep)/eV 
           << " shell= " << i 
           << "  E(keV)= " << energy/keV 
           << "  Ebind(keV)= " << bindingEnergy/keV 
           << "  Eg(keV)= " << gamEnergy1/keV 
           << "  Ee(keV)= " << eKinEnergy/keV 
           << "  Esec(keV)= " << esec/keV 
           << "  Edep(keV)= " << edep/keV 
           << G4endl;
  }
  // energy balance
  if(edep > 0.0) { 
    fParticleChange->ProposeLocalEnergyDeposit(edep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

