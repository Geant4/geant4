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
// $Id: G4KL3DecayChannel.cc 95906 2016-03-02 10:56:50Z gcosmo $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      30 May 1997 H.Kurashige
// ------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4KL3DecayChannel.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"

G4KL3DecayChannel::G4KL3DecayChannel()
  :G4VDecayChannel(),
   pLambda(0.0), pXi0(0.0)
{
}


G4KL3DecayChannel::G4KL3DecayChannel(
			const G4String& theParentName, 
			G4double        theBR,
			const G4String& thePionName,
			const G4String& theLeptonName,
			const G4String& theNutrinoName)
                   :G4VDecayChannel("KL3 Decay",theParentName,
				   theBR,  3,
				   thePionName,theLeptonName,theNutrinoName)
{
  static const G4String K_plus("kaon+");
  static const G4String K_minus("kaon-");
  static const G4String K_L("kaon0L");
  static const G4String Mu_plus("mu+");
  static const G4String Mu_minus("mu-");
  static const G4String E_plus("e+");
  static const G4String E_minus("e-");
  
  // check modes
  if ( ((theParentName == K_plus)&&(theLeptonName == E_plus)) ||
       ((theParentName == K_minus)&&(theLeptonName == E_minus))   ) {
    // K+- (Ke3)
    pLambda = 0.0286;
    pXi0    = -0.35;
   } else if ( ((theParentName == K_plus)&&(theLeptonName == Mu_plus)) ||
       ((theParentName == K_minus)&&(theLeptonName == Mu_minus))   ) {
    // K+- (Kmu3)
    pLambda = 0.033;
    pXi0    = -0.35;
  } else if ( (theParentName == K_L) && 
              ((theLeptonName == E_plus) ||(theLeptonName == E_minus))  ){
    // K0L (Ke3)
    pLambda = 0.0300;
    pXi0    = -0.11;
  } else if ( (theParentName == K_L) && 
              ((theLeptonName == Mu_plus) ||(theLeptonName == Mu_minus))  ){
    // K0L (Kmu3)
    pLambda = 0.034;
    pXi0    = -0.11;
  } else {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>2) {
      G4cout << "G4KL3DecayChannel:: constructor :";
      G4cout << "illegal arguments " << G4endl;;
      DumpInfo();
    }
#endif
    // set values for K0L (Ke3) temporarily
    pLambda = 0.0300;
    pXi0    = -0.11;
  }
}

G4KL3DecayChannel::~G4KL3DecayChannel()
{
}

G4KL3DecayChannel::G4KL3DecayChannel(const G4KL3DecayChannel &right):
  G4VDecayChannel(right),
  //massK(right.massK),
  pLambda(right.pLambda), 
  pXi0(right.pXi0)
{
}

G4KL3DecayChannel & G4KL3DecayChannel::operator=(const G4KL3DecayChannel & right)
{
  if (this != &right) { 
    kinematics_name = right.kinematics_name;
    verboseLevel = right.verboseLevel;
    rbranch = right.rbranch;

    // copy parent name
    parent_name = new G4String(*right.parent_name);

    // clear daughters_name array
    ClearDaughtersName();

    // recreate array
    numberOfDaughters = right.numberOfDaughters;
    if ( numberOfDaughters >0 ) {
      if (daughters_name !=0) ClearDaughtersName();
      daughters_name = new G4String*[numberOfDaughters];
      //copy daughters name
      for (G4int index=0; index < numberOfDaughters; index++) {
          daughters_name[index] = new G4String(*right.daughters_name[index]);
      }
    }
    //massK = right.massK;
    pLambda = right.pLambda; 
    pXi0 = right.pXi0;
  }
  return *this;
}


G4DecayProducts* G4KL3DecayChannel::DecayIt(G4double) 
{
  // this version neglects muon polarization 
  //              assumes the pure V-A coupling
  //              gives incorrect energy spectrum for Nutrinos
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4KL3DecayChannel::DecayIt " << G4endl;
#endif

  // fill parent particle and its mass
  CheckAndFillParent();
  G4double massK = G4MT_parent->GetPDGMass();

  // fill daughter particles and their mass
  CheckAndFillDaughters();
  G4double daughterM[3];
  daughterM[idPi] = G4MT_daughters[idPi]->GetPDGMass();
  daughterM[idLepton] = G4MT_daughters[idLepton]->GetPDGMass();
  daughterM[idNutrino] = G4MT_daughters[idNutrino]->GetPDGMass();

  // determine momentum/energy of daughters 
  //  according to DalitzDensity 
  G4double daughterP[3], daughterE[3];
  G4double w;
  G4double r;
  const size_t MAX_LOOP = 10000;
  for (size_t loop_counter=0; loop_counter <MAX_LOOP; ++loop_counter){
    r = G4UniformRand();
    PhaseSpace(massK, &daughterM[0], &daughterE[0], &daughterP[0]);
    w = DalitzDensity(massK,daughterE[idPi],daughterE[idLepton],daughterE[idNutrino],
                      daughterM[idPi],daughterM[idLepton],daughterM[idNutrino]);
    if ( r <= w) break;
  }

  // output message
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << *daughters_name[0] << ":" << daughterP[0]/GeV << "[GeV/c]" <<G4endl;
    G4cout << *daughters_name[1] << ":" << daughterP[1]/GeV << "[GeV/c]" <<G4endl;
    G4cout << *daughters_name[2] << ":" << daughterP[2]/GeV << "[GeV/c]" <<G4endl;
  }
#endif
   //create parent G4DynamicParticle at rest
  G4ThreeVector* direction = new G4ThreeVector(1.0,0.0,0.0);
  G4DynamicParticle * parentparticle = new G4DynamicParticle( G4MT_parent, *direction, 0.0);
  delete direction;

  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  //create daughter G4DynamicParticle 
  G4double costheta, sintheta, phi, sinphi, cosphi; 
  G4double costhetan, sinthetan, phin, sinphin, cosphin;
 
  // pion
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
  phi  = twopi*G4UniformRand()*rad;
  sinphi = std::sin(phi);
  cosphi = std::cos(phi);
  direction = new G4ThreeVector(sintheta*cosphi,sintheta*sinphi,costheta);
  G4ThreeVector momentum0 =  (*direction)*daughterP[0]; 
  G4DynamicParticle * daughterparticle 
       = new G4DynamicParticle( G4MT_daughters[0], momentum0);
  products->PushProducts(daughterparticle);

  // neutrino
  costhetan = (daughterP[1]*daughterP[1]-daughterP[2]*daughterP[2]-daughterP[0]*daughterP[0])/(2.0*daughterP[2]*daughterP[0]);
  sinthetan = std::sqrt((1.0-costhetan)*(1.0+costhetan));
  phin  = twopi*G4UniformRand()*rad;
  sinphin = std::sin(phin);
  cosphin = std::cos(phin);
  direction->setX( sinthetan*cosphin*costheta*cosphi - sinthetan*sinphin*sinphi + costhetan*sintheta*cosphi); 
  direction->setY( sinthetan*cosphin*costheta*sinphi + sinthetan*sinphin*cosphi + costhetan*sintheta*sinphi); 
  direction->setZ( -sinthetan*cosphin*sintheta + costhetan*costheta);

  G4ThreeVector momentum2 =  (*direction)*daughterP[2]; 
  daughterparticle = new G4DynamicParticle( G4MT_daughters[2], momentum2);
  products->PushProducts(daughterparticle);

  //lepton
  G4ThreeVector momentum1 = (momentum0 + momentum2) * (-1.0);
  daughterparticle = 
       new G4DynamicParticle( G4MT_daughters[1], momentum1);
  products->PushProducts(daughterparticle);

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
     G4cout << "G4KL3DecayChannel::DecayIt ";
     G4cout << "  create decay products in rest frame " <<G4endl;
     G4cout << "  decay products address=" << products << G4endl;
     products->DumpInfo();
  }
#endif
  delete direction;
  return products;
}

void G4KL3DecayChannel::PhaseSpace(G4double parentM,
				   const G4double* M,
				   G4double*       E,
				   G4double*       P )
// algorism of this code is originally written in GDECA3 of GEANT3
{
  
  //sum of daughters'mass
  G4double sumofdaughtermass = 0.0;
  G4int index;
  const G4int N_DAUGHTER=3;
  
  for (index=0; index<N_DAUGHTER; index++){
    sumofdaughtermass += M[index];
  }

  //calculate daughter momentum
  //  Generate two 
  G4double rd1, rd2, rd;
  G4double momentummax=0.0, momentumsum = 0.0;
  G4double energy;
  const size_t MAX_LOOP=10000;
  for (size_t loop_counter=0; loop_counter <MAX_LOOP; ++loop_counter){
    rd1 = G4UniformRand();
    rd2 = G4UniformRand();
    if (rd2 > rd1) {
      rd  = rd1;
      rd1 = rd2;
      rd2 = rd;
    } 
    momentummax = 0.0;
    momentumsum = 0.0;
    // daughter 0
    energy = rd2*(parentM - sumofdaughtermass);
    P[0] = std::sqrt(energy*energy + 2.0*energy*M[0]);
    E[0] = energy;
    if ( P[0] >momentummax )momentummax =  P[0];
    momentumsum  +=  P[0];
    // daughter 1
    energy = (1.-rd1)*(parentM - sumofdaughtermass);
    P[1] = std::sqrt(energy*energy + 2.0*energy*M[1]);
    E[1] = energy;
    if ( P[1] >momentummax )momentummax =  P[1];
    momentumsum  +=  P[1];
    // daughter 2
    energy = (rd1-rd2)*(parentM - sumofdaughtermass);
    P[2] = std::sqrt(energy*energy + 2.0*energy*M[2]);
    E[2] = energy;
    if ( P[2] >momentummax )momentummax =  P[2];
    momentumsum  +=  P[2];
    if (momentummax <=  momentumsum - momentummax ) break;
  }
#ifdef G4VERBOSE
  if (GetVerboseLevel()>2) {
     G4cout << "G4KL3DecayChannel::PhaseSpace    ";
     G4cout << "Kon mass:" << parentM/GeV << "GeV/c/c" << G4endl;
     for (index=0; index<3; index++){
       G4cout << index << " : " << M[index]/GeV << "GeV/c/c  ";
       G4cout << " : " << E[index]/GeV << "GeV  ";
       G4cout << " : " << P[index]/GeV << "GeV/c " << G4endl;
     }
  }
#endif
}


G4double G4KL3DecayChannel::DalitzDensity(G4double massK, G4double Epi, G4double El, G4double Enu,
                                          G4double massPi, G4double massL , G4double massNu )
{
  // KL3 decay   Dalitz Plot Density
  //               see Chounet et al Phys. Rep. 4, 201
  //  arguments
  //    Epi: kinetic enregy of pion
  //    El:  kinetic enregy of lepton (e or mu)
  //    Enu: kinetic energy of nutrino
  //  constants
  //    pLambda : linear energy dependence of f+
  //    pXi0    : = f+(0)/f-
  //    pNorm   : normalization factor
  //  variables
  //    Epi: total energy of pion
  //    El:  total energy of lepton (e or mu)
  //    Enu: total energy of nutrino

  // calcurate total energy
  Epi = Epi + massPi;
  El  = El  + massL;
  Enu = Enu + massNu;
  
  G4double Epi_max = (massK*massK+massPi*massPi-massL*massL)/2.0/massK;
  G4double E  = Epi_max - Epi;
  G4double q2 = massK*massK + massPi*massPi - 2.0*massK*Epi;

  G4double F    = 1.0 + pLambda*q2/massPi/massPi;
  G4double Fmax = 1.0;
  if (pLambda >0.0) Fmax = (1.0 + pLambda*(massK*massK/massPi/massPi+1.0));

  G4double Xi = pXi0*(1.0 + pLambda*q2/massPi/massPi);

  G4double coeffA = massK*(2.0*El*Enu-massK*E)+massL*massL*(E/4.0-Enu);
  G4double coeffB = massL*massL*(Enu-E/2.0);
  G4double coeffC = massL*massL*E/4.0;

  G4double RhoMax = (Fmax*Fmax)*(massK*massK*massK/8.0);

  G4double Rho = (F*F)*(coeffA + coeffB*Xi + coeffC*Xi*Xi);
 
#ifdef G4VERBOSE
  if (GetVerboseLevel()>2) {
    G4cout << "G4KL3DecayChannel::DalitzDensity  " <<G4endl;
    G4cout << " Pi[" << massPi/GeV <<"GeV/c/c] :" << Epi/GeV << "GeV" <<G4endl;
    G4cout << " L[" << massL/GeV <<"GeV/c/c] :" << El/GeV << "GeV" <<G4endl;
    G4cout << " Nu[" << massNu/GeV <<"GeV/c/c] :" << Enu/GeV << "GeV" <<G4endl;
    G4cout << " F :" << F  << " Fmax :" << Fmax << "  Xi :" << Xi << G4endl;
    G4cout << " A :" << coeffA << "  B :" << coeffB << "  C :"<< coeffC <<G4endl; 
    G4cout << " Rho :" << Rho  << "   RhoMax :" << RhoMax << G4endl;
  }
#endif
  return (Rho/RhoMax);
}


