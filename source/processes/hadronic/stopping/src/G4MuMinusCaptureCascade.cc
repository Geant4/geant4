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
// $Id: G4MuMinusCaptureCascade.cc,v 1.16 2008-05-05 09:09:06 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//   G4MuonMinusCaptureAtRest physics process
//                   
//   E-mail: Vladimir.Ivantchenko@cern.ch
//
//   Created:   02.04.00 V.Ivanchenko
//
// Modified:  
// 06.04.01 V.Ivanchenko Bug in theta distribution fixed
// 13.02.07 V.Ivanchenko Fixes in decay - add random distribution of e- 
//                       direction; factor 2 in potential energy 
//
//----------------------------------------------------------------------

#include "G4MuMinusCaptureCascade.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4MuonMinus.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4GHEKinematicsVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuMinusCaptureCascade::G4MuMinusCaptureCascade()
{ 
  theElectron = G4Electron::Electron();
  theGamma = G4Gamma::Gamma();
  Emass = theElectron->GetPDGMass();
  MuMass = G4MuonMinus::MuonMinus()->GetPDGMass();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuMinusCaptureCascade::~G4MuMinusCaptureCascade()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuMinusCaptureCascade::GetKShellEnergy(G4double Z)
{ 
  // Calculate the Energy of K Mesoatom Level for this Element using
  // the Energy of Hydrogen Atom taken into account finite size of the
  // nucleus (V.Ivanchenko)
  const G4int ListK = 28;
  static G4double ListZK[ListK] = {
      1., 2.,  4.,  6.,  8., 11., 14., 17., 18., 21., 24.,
     26., 29., 32., 38., 40., 41., 44., 49., 53., 55.,
     60., 65., 70., 75., 81., 85., 92.};
  static G4double ListKEnergy[ListK] = {
     0.00275, 0.011, 0.043, 0.098, 0.173, 0.326,
     0.524, 0.765, 0.853, 1.146, 1.472,
     1.708, 2.081, 2.475, 3.323, 3.627, 
     3.779, 4.237, 5.016, 5.647, 5.966,
     6.793, 7.602, 8.421, 9.249, 10.222,
    10.923,11.984};

  // Energy with finit size corrections
  G4double KEnergy = GetLinApprox(ListK,ListZK,ListKEnergy,Z);

  return KEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuMinusCaptureCascade::AddNewParticle(G4ParticleDefinition* aParticle,
                                             G4ThreeVector& Momentum,
                                             G4double mass,
                                             G4int* nParticle,
                                             G4GHEKinematicsVector* Cascade)
{
  // Store particle in the HEK vector and increment counter
  Cascade[*nParticle].SetZero();
  Cascade[*nParticle].SetMass( mass );
  Cascade[*nParticle].SetMomentumAndUpdate(Momentum.x(), Momentum.y(), Momentum.z());
  Cascade[*nParticle].SetParticleDef( aParticle );
  (*nParticle)++;

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4MuMinusCaptureCascade::DoCascade(const G4double Z, const G4double massA, 
					 G4GHEKinematicsVector* Cascade)
{
  // Inicialization - cascade start from 14th level
  // N.C.Mukhopadhyay Phy. Rep. 30 (1977) 1.
  G4int nPart = 0;
  G4double EnergyLevel[14];

  G4double mass = MuMass * massA / (MuMass + massA) ;

  const G4double KEnergy = 13.6 * eV * Z * Z * mass/ electron_mass_c2;

  EnergyLevel[0] = GetKShellEnergy(Z);
  for( G4int i = 2; i < 15; i++ ) {
    EnergyLevel[i-1] = KEnergy / (i*i) ;
  }

  G4int nElec  = G4int(Z);
  G4int nAuger = 1;
  G4int nLevel = 13;
  G4double DeltaE;
  G4double pGamma = Z*Z*Z*Z;

  // Capture on 14-th level
  G4double ptot = std::sqrt(EnergyLevel[13]*(EnergyLevel[13] + 2.0*Emass));
  G4ThreeVector moment = ptot * GetRandomVec();

  AddNewParticle(theElectron,moment,Emass,&nPart,Cascade);

  // Emit new photon or electron
  // Simplified model for probabilities
  // N.C.Mukhopadhyay Phy. Rep. 30 (1977) 1.
  do {

    // case of Auger electrons
    if((nAuger < nElec) && ((pGamma + 10000.0) * G4UniformRand() < 10000.0) ) {
        nAuger++;
        DeltaE = EnergyLevel[nLevel-1] - EnergyLevel[nLevel];
        nLevel--;

        ptot = std::sqrt(DeltaE * (DeltaE + 2.0*Emass));
        moment = ptot * GetRandomVec();

        AddNewParticle(theElectron, moment, Emass, &nPart, Cascade);

    } else {

      // Case of photon cascade, probabilities from
      // C.S.Wu and L.Wilets, Ann. Rev. Nuclear Sci. 19 (1969) 527.

      G4double var = (10.0 + G4double(nLevel - 1) ) * G4UniformRand();
      G4int iLevel = nLevel - 1 ;
      if(var > 10.0) iLevel -= G4int(var-10.0) + 1;
      if( iLevel < 0 ) iLevel = 0;
      DeltaE = EnergyLevel[iLevel] - EnergyLevel[nLevel];
      nLevel = iLevel;
      moment = DeltaE * GetRandomVec();
      AddNewParticle(theGamma, moment, 0.0, &nPart, Cascade);
    }

  } while( nLevel > 0 );

  return nPart;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuMinusCaptureCascade::DoBoundMuonMinusDecay(G4double Z, 
                                                    G4int* nCascade, 
                                                    G4GHEKinematicsVector* Cascade)
{
  // Simulation on Decay of mu- on a K-shell of the muonic atom
  G4double xmax = ( 1.0 + Emass*Emass/ (MuMass*MuMass) );
  G4double xmin = 2.0*Emass/MuMass;
  G4double KEnergy = GetKShellEnergy(Z);
  /*
  G4cout << "G4MuMinusCaptureCascade::DoBoundMuonMinusDecay" 
	 << " XMAX= " << xmax
	 << " Ebound= " << KEnergy
	 << G4endl;
  */
  G4double pmu = std::sqrt(KEnergy*(KEnergy + 2.0*MuMass));
  G4double emu = KEnergy + MuMass;
  G4ThreeVector moment = GetRandomVec();
  G4LorentzVector MU(pmu*moment,emu);
  G4ThreeVector bst = MU.boostVector();

  G4double Eelect, Pelect, x, ecm;
  G4LorentzVector EL, NN;
  // Calculate electron energy
  do {
    do {
      x = xmin + (xmax-xmin)*G4UniformRand();
    } while (G4UniformRand() > (3.0 - 2.0*x)*x*x );
    Eelect = x*MuMass*0.5;
    Pelect = 0.0;
    if(Eelect > Emass) { 
      Pelect = std::sqrt( Eelect*Eelect - Emass*Emass );
    } else {
      Pelect = 0.0;
      Eelect = Emass;
    }
    G4ThreeVector e_mom = GetRandomVec();
    EL = G4LorentzVector(Pelect*e_mom,Eelect);
    EL.boost(bst);
    Eelect = EL.e() - Emass - 2.0*KEnergy;
    //
    // Calculate rest frame parameters of 2 neutrinos
    //
    NN = MU - EL;
    ecm = NN.mag2();
  } while (Eelect < 0.0 || ecm < 0.0);

  //
  // Create electron
  //
  moment = std::sqrt(Eelect * (Eelect + 2.0*Emass))*(EL.vect().unit());
  AddNewParticle(theElectron, moment, Emass, nCascade, Cascade);
  //
  // Create Neutrinos
  //
  ecm = 0.5*std::sqrt(ecm);
  bst = NN.boostVector();
  G4ThreeVector p1 = ecm * GetRandomVec();
  G4LorentzVector N1 = G4LorentzVector(p1,ecm);
  N1.boost(bst);
  G4ThreeVector p1lab = N1.vect();
  AddNewParticle(G4AntiNeutrinoE::AntiNeutrinoE(),p1lab,0.0,nCascade,Cascade);
  NN -= N1;
  G4ThreeVector p2lab = NN.vect();
  AddNewParticle(G4NeutrinoMu::NeutrinoMu(),p2lab,0.0,nCascade,Cascade);

  return;
}
















