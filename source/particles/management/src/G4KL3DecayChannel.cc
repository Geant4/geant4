// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4KL3DecayChannel.cc,v 1.1 1999-01-07 16:10:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      30 May 1997 H.Kurashige
// ------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4KL3DecayChannel.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"


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
  //#ifdef G4VERBOSE
  //if (GetVerboseLevel()>1) {
  //  G4cout << "G4KL3DecayChannel:: constructor ";
  //  G4cout << "addr[" << this << "]" << endl;
  //}
  //#endif
  // check modes
  if ( ((theParentName == "kaon+")&&(theLeptonName == "e+")) ||
       ((theParentName == "kaon-")&&(theLeptonName == "e-"))   ) {
    // K+- (Ke3)
    pLambda = 0.0286;
    pXi0    = -0.35;
   } else if ( ((theParentName == "kaon+")&&(theLeptonName == "mu+")) ||
       ((theParentName == "kaon-")&&(theLeptonName == "mu-"))   ) {
    // K+- (Kmu3)
    pLambda = 0.033;
    pXi0    = -0.35;
  } else if ( (theParentName == "kaon0L") && 
              ((theLeptonName == "e+") ||(theLeptonName == "e-"))  ){
    // K0L (Ke3)
    pLambda = 0.0300;
    pXi0    = -0.11;
  } else if ( (theParentName == "kaon0L") && 
              ((theLeptonName == "mu+") ||(theLeptonName == "mu-"))  ){
    // K0L (Kmu3)
    pLambda = 0.034;
    pXi0    = -0.11;
  } else {
    //#ifdef G4VERBOSE
    //if (GetVerboseLevel()>0) {
    //   G4cout << "G4KL3DecayChannel:: constructor :";
    //   G4cout << "illegal arguments " << endl;;
    //   DumpInfo();
    // }
    //#endif
    // set values for K0L (Ke3) temporarily
    pLambda = 0.0300;
    pXi0    = -0.11;
  }
}

G4KL3DecayChannel::~G4KL3DecayChannel()
{
}

G4DecayProducts* G4KL3DecayChannel::DecayIt(G4double) 
{
  // this version neglects muon polarization 
  //              assumes the pure V-A coupling
  //              gives incorrect energy spectrum for Nutrinos
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4KL3DecayChannel::DecayIt " << endl;
#endif
  // fill parent particle and its mass
  if (parent == NULL) {
    FillParent();
  }
  massK = parent->GetPDGMass();

  // fill daughter particles and their mass
  if (daughters == NULL) {
    FillDaughters();
  }
  daughterM[idPi] = daughters[idPi]->GetPDGMass();
  daughterM[idLepton] = daughters[idLepton]->GetPDGMass();
  daughterM[idNutrino] = daughters[idNutrino]->GetPDGMass();

  // determine momentum/energy of daughters 
  //  according to DalitzDensity 
  G4double daughterP[3], daughterE[3];
  G4double w;
  G4double r;
  do {
    r = G4UniformRand();
    PhaseSpace(massK, &daughterM[0], &daughterE[0], &daughterP[0]);
    w = DalitzDensity(daughterE[idPi],daughterE[idLepton],daughterE[idNutrino]);
  } while ( r > w);

  // output message
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << *daughters_name[0] << ":" << daughterP[0]/GeV << "[GeV/c]" <<endl;
    G4cout << *daughters_name[1] << ":" << daughterP[1]/GeV << "[GeV/c]" <<endl;
    G4cout << *daughters_name[2] << ":" << daughterP[2]/GeV << "[GeV/c]" <<endl;
  }
#endif
   //create parent G4DynamicParticle at rest
  G4ParticleMomentum* direction = new G4ParticleMomentum(1.0,0.0,0.0);
  G4DynamicParticle * parentparticle = new G4DynamicParticle( parent, *direction, 0.0);
  delete direction;

  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  //create daughter G4DynamicParticle 
  G4double costheta, sintheta, phi, sinphi, cosphi; 
  G4double costhetan, sinthetan, phin, sinphin, cosphin;
 
  // pion
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = sqrt((1.0-costheta)*(1.0+costheta));
  phi  = 2.0*M_PI*G4UniformRand()*rad;
  sinphi = sin(phi);
  cosphi = cos(phi);
  direction = new G4ParticleMomentum(sintheta*cosphi,sintheta*sinphi,costheta);
  G4ThreeVector momentum0 =  (*direction)*daughterP[0]; 
  G4DynamicParticle * daughterparticle 
       = new G4DynamicParticle( daughters[0], momentum0);
  products->PushProducts(daughterparticle);

  // neutrino
  costhetan = (daughterP[1]*daughterP[1]-daughterP[2]*daughterP[2]-daughterP[0]*daughterP[0])/(2.0*daughterP[2]*daughterP[0]);
  sinthetan = sqrt((1.0-costhetan)*(1.0+costhetan));
  phin  = 2.0*M_PI*G4UniformRand()*rad;
  sinphin = sin(phin);
  cosphin = cos(phin);
  direction->setX( sinthetan*cosphin*costheta*cosphi - sinthetan*sinphin*sinphi + costhetan*sintheta*cosphi); 
  direction->setY( sinthetan*cosphin*costheta*sinphi + sinthetan*sinphin*cosphi + costhetan*sintheta*sinphi); 
  direction->setZ( -sinthetan*cosphin*sintheta + costhetan*costheta);

  G4ThreeVector momentum2 =  (*direction)*daughterP[2]; 
  daughterparticle = new G4DynamicParticle( daughters[2], momentum2);
  products->PushProducts(daughterparticle);

  //lepton
  G4ThreeVector momentum1 = (momentum0 + momentum2) * (-1.0);
  daughterparticle = 
       new G4DynamicParticle( daughters[1], momentum1);
  products->PushProducts(daughterparticle);

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
     G4cout << "G4KL3DecayChannel::DecayIt ";
     G4cout << "  create decay products in rest frame " <<endl;
     G4cout << "  decay products address=" << products << endl;
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
  for (index=0; index<3; index++){
    sumofdaughtermass += M[index];
  }

  //calculate daughter momentum
  //  Generate two 
  G4double rd1, rd2, rd;
  G4double momentummax=0.0, momentumsum = 0.0;
  G4double energy;

  do {
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
    P[0] = sqrt(energy*energy + 2.0*energy*M[0]);
    E[0] = energy;
    if ( P[0] >momentummax )momentummax =  P[0];
    momentumsum  +=  P[0];
    // daughter 1
    energy = (1.-rd1)*(parentM - sumofdaughtermass);
    P[1] = sqrt(energy*energy + 2.0*energy*M[1]);
    E[1] = energy;
    if ( P[1] >momentummax )momentummax =  P[1];
    momentumsum  +=  P[1];
    // daughter 2
    energy = (rd1-rd2)*(parentM - sumofdaughtermass);
    P[2] = sqrt(energy*energy + 2.0*energy*M[2]);
    E[2] = energy;
    if ( P[2] >momentummax )momentummax =  P[2];
    momentumsum  +=  P[2];
  } while (momentummax >  momentumsum - momentummax );

#ifdef G4VERBOSE
  if (GetVerboseLevel()>2) {
     G4cout << "G4KL3DecayChannel::PhaseSpace    ";
     G4cout << "Kon mass:" << parentM/GeV << "GeV/c/c" << endl;
     for (index=0; index<3; index++){
       G4cout << index << " : " << M[index]/GeV << "GeV/c/c  ";
       G4cout << " : " << E[index]/GeV << "GeV  ";
       G4cout << " : " << P[index]/GeV << "GeV/c " << endl;
     }
  }
#endif
}


G4double G4KL3DecayChannel::DalitzDensity(G4double Epi, G4double El, G4double Enu)
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

  // mass of daughters
  G4double massPi = daughterM[idPi];
  G4double massL  = daughterM[idLepton]; 
  G4double massNu = daughterM[idNutrino];

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
    G4cout << "G4KL3DecayChannel::DalitzDensity  " <<endl;
    G4cout << " Pi[" << massPi/GeV <<"GeV/c/c] :" << Epi/GeV << "GeV" <<endl;
    G4cout << " L[" << massL/GeV <<"GeV/c/c] :" << El/GeV << "GeV" <<endl;
    G4cout << " Nu[" << massNu/GeV <<"GeV/c/c] :" << Enu/GeV << "GeV" <<endl;
    G4cout << " F :" << F  << " Fmax :" << Fmax << "  Xi :" << Xi << endl;
    G4cout << " A :" << coeffA << "  B :" << coeffB << "  C :"<< coeffC <<endl; 
    G4cout << " Rho :" << Rho  << "   RhoMax :" << RhoMax << endl;
  }
#endif
  return (Rho/RhoMax);
}


