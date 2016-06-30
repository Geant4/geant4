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
// $Id: G4GeneralPhaseSpaceDecay.cc 95909 2016-03-02 11:01:47Z gcosmo $
// ----------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, A. Feliciello, 21st May 1998
//
//      Note: this class is a generalization of the 
//            G4PhaseSpaceDecayChannel one
// ----------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4GeneralPhaseSpaceDecay.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4ios.hh"


G4GeneralPhaseSpaceDecay::G4GeneralPhaseSpaceDecay(G4int Verbose) : 
                          G4VDecayChannel("Phase Space", Verbose),
                          parentmass(0.), theDaughterMasses(0)
{
  if (GetVerboseLevel()>1) G4cout << "G4GeneralPhaseSpaceDecay:: constructor " << G4endl;
}

G4GeneralPhaseSpaceDecay::G4GeneralPhaseSpaceDecay(const G4String& theParentName,
			                           G4double        theBR,
			                           G4int           theNumberOfDaughters,
			                           const G4String& theDaughterName1,
			                           const G4String& theDaughterName2,
			                           const G4String& theDaughterName3) :
                                   G4VDecayChannel("Phase Space",
					           theParentName,theBR,
					           theNumberOfDaughters,
					           theDaughterName1,
					           theDaughterName2,
					           theDaughterName3),
				   theDaughterMasses(0)
{
  if (GetVerboseLevel()>1) G4cout << "G4GeneralPhaseSpaceDecay:: constructor " << G4endl;
  
  //   Set the parent particle (resonance) mass to the (default) PDG vale
  if (G4MT_parent != NULL)
  {
      parentmass = G4MT_parent->GetPDGMass();
  } else {
	  parentmass=0.;
  }

}

G4GeneralPhaseSpaceDecay::G4GeneralPhaseSpaceDecay(const G4String& theParentName,
                                                   G4double        theParentMass,
			                           G4double        theBR,
			                           G4int           theNumberOfDaughters,
			                           const G4String& theDaughterName1,
			                           const G4String& theDaughterName2,
			                           const G4String& theDaughterName3) :
                                   G4VDecayChannel("Phase Space",
					           theParentName,theBR,
					           theNumberOfDaughters,
					           theDaughterName1,
					           theDaughterName2,
					           theDaughterName3),
			           parentmass(theParentMass),
				   theDaughterMasses(0)
{
  if (GetVerboseLevel()>1) G4cout << "G4GeneralPhaseSpaceDecay:: constructor " << G4endl;
}

G4GeneralPhaseSpaceDecay::G4GeneralPhaseSpaceDecay(const G4String& theParentName,
                                                   G4double        theParentMass,
			                           G4double        theBR,
			                           G4int           theNumberOfDaughters,
			                           const G4String& theDaughterName1,
			                           const G4String& theDaughterName2,
			                           const G4String& theDaughterName3,
						   const G4double *masses) :
                                   G4VDecayChannel("Phase Space",
					           theParentName,theBR,
					           theNumberOfDaughters,
					           theDaughterName1,
					           theDaughterName2,
					           theDaughterName3),
			           parentmass(theParentMass),
				   theDaughterMasses(masses)
{
  if (GetVerboseLevel()>1) G4cout << "G4GeneralPhaseSpaceDecay:: constructor " << G4endl;
}

G4GeneralPhaseSpaceDecay::G4GeneralPhaseSpaceDecay(const G4String& theParentName,
                                                   G4double        theParentMass,
			                           G4double        theBR,
			                           G4int           theNumberOfDaughters,
			                           const G4String& theDaughterName1,
			                           const G4String& theDaughterName2,
			                           const G4String& theDaughterName3,
			                           const G4String& theDaughterName4,
						   const G4double *masses) :
                                   G4VDecayChannel("Phase Space",
					           theParentName,theBR,
					           theNumberOfDaughters,
					           theDaughterName1,
					           theDaughterName2,
					           theDaughterName3,
					           theDaughterName4),
			           parentmass(theParentMass),
				   theDaughterMasses(masses)
{
  if (GetVerboseLevel()>1) G4cout << "G4GeneralPhaseSpaceDecay:: constructor " << G4endl;
}

G4GeneralPhaseSpaceDecay::~G4GeneralPhaseSpaceDecay()
{
}

G4DecayProducts *G4GeneralPhaseSpaceDecay::DecayIt(G4double) 
{
  if (GetVerboseLevel()>1) G4cout << "G4GeneralPhaseSpaceDecay::DecayIt ";
  G4DecayProducts * products = NULL;
 
  CheckAndFillParent();
  CheckAndFillDaughters();
  
  switch (numberOfDaughters){
  case 0:
    if (GetVerboseLevel()>0) {
      G4cout << "G4GeneralPhaseSpaceDecay::DecayIt ";
      G4cout << " daughters not defined " <<G4endl;
    }
    break;
  case 1:
    products =  OneBodyDecayIt();
    break;
  case 2:
    products =  TwoBodyDecayIt();
    break;
  case 3:
    products =  ThreeBodyDecayIt();
    break;
  default:
    products =  ManyBodyDecayIt();
    break;
  }
  if ((products == NULL) && (GetVerboseLevel()>0)) {
    G4cout << "G4GeneralPhaseSpaceDecay::DecayIt ";
    G4cout << *parent_name << " can not decay " << G4endl;
    DumpInfo();
  }
  return products;
}

G4DecayProducts *G4GeneralPhaseSpaceDecay::OneBodyDecayIt()
{
  if (GetVerboseLevel()>1) G4cout << "G4GeneralPhaseSpaceDecay::OneBodyDecayIt()"<<G4endl;

//  G4double daughtermass = daughters[0]->GetPDGMass();

  //create parent G4DynamicParticle at rest
  G4ParticleMomentum dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle(G4MT_parent, dummy, 0.0);

  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  //create daughter G4DynamicParticle at rest
  G4DynamicParticle * daughterparticle = new G4DynamicParticle(G4MT_daughters[0], dummy, 0.0);
  products->PushProducts(daughterparticle);

  if (GetVerboseLevel()>1) 
    {
     G4cout << "G4GeneralPhaseSpaceDecay::OneBodyDecayIt ";
     G4cout << "  create decay products in rest frame " <<G4endl;
     products->DumpInfo();
    }
  return products;
}

G4DecayProducts *G4GeneralPhaseSpaceDecay::TwoBodyDecayIt()
{
  if (GetVerboseLevel()>1) G4cout << "G4GeneralPhaseSpaceDecay::TwoBodyDecayIt()"<<G4endl;
  
  //daughters'mass
  G4double daughtermass[2]; 
  G4double daughtermomentum;
  if ( theDaughterMasses )
  { 
     daughtermass[0]= *(theDaughterMasses);
     daughtermass[1] = *(theDaughterMasses+1);
  } else {   
     daughtermass[0] = G4MT_daughters[0]->GetPDGMass();
     daughtermass[1] = G4MT_daughters[1]->GetPDGMass();
  }
  
//  G4double sumofdaughtermass =  daughtermass[0] + daughtermass[1];

  //create parent G4DynamicParticle at rest
  G4ParticleMomentum dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( G4MT_parent, dummy, 0.0);

  //create G4Decayproducts  @@GF why dummy parentparticle?
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  //calculate daughter momentum
  daughtermomentum = Pmx(parentmass,daughtermass[0],daughtermass[1]);
  G4double costheta = 2.*G4UniformRand()-1.0;
  G4double sintheta = std::sqrt((1.0 - costheta)*(1.0 + costheta));
  G4double phi  = twopi*G4UniformRand()*rad;
  G4ParticleMomentum direction(sintheta*std::cos(phi),sintheta*std::sin(phi),costheta);

  //create daughter G4DynamicParticle
  G4double Etotal= std::sqrt(daughtermass[0]*daughtermass[0] + daughtermomentum*daughtermomentum); 
  G4DynamicParticle * daughterparticle = new G4DynamicParticle( G4MT_daughters[0],Etotal, direction*daughtermomentum);
  products->PushProducts(daughterparticle);
  Etotal= std::sqrt(daughtermass[1]*daughtermass[1] + daughtermomentum*daughtermomentum);
  daughterparticle = new G4DynamicParticle( G4MT_daughters[1],Etotal, direction*(-1.0*daughtermomentum));
  products->PushProducts(daughterparticle);

  if (GetVerboseLevel()>1) 
    {
     G4cout << "G4GeneralPhaseSpaceDecay::TwoBodyDecayIt ";
     G4cout << "  create decay products in rest frame " <<G4endl;
     products->DumpInfo();
    }
  return products;
}

G4DecayProducts *G4GeneralPhaseSpaceDecay::ThreeBodyDecayIt()
// algorism of this code is originally written in GDECA3 of GEANT3
{
  if (GetVerboseLevel()>1) G4cout << "G4GeneralPhaseSpaceDecay::ThreeBodyDecayIt()"<<G4endl;

  //daughters'mass
  G4double daughtermass[3]; 
  G4double sumofdaughtermass = 0.0;
  for (G4int index=0; index<3; index++)
    {
     if ( theDaughterMasses )
     { 
         daughtermass[index]= *(theDaughterMasses+index);
     } else {   
         daughtermass[index] = G4MT_daughters[index]->GetPDGMass();
     }   
     sumofdaughtermass += daughtermass[index];
    }
  
  //create parent G4DynamicParticle at rest
  G4ParticleMomentum dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( G4MT_parent, dummy, 0.0);

  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  //calculate daughter momentum
  //  Generate two 
  G4double rd1, rd2, rd;
  G4double daughtermomentum[3];
  G4double momentummax=0.0, momentumsum = 0.0;
  G4double energy;
  const G4int maxNumberOfLoops = 10000;
  G4int loopCounter = 0;

  do 
    {
     rd1 = G4UniformRand();
     rd2 = G4UniformRand();
     if (rd2 > rd1) 
       {
        rd  = rd1;
        rd1 = rd2;
        rd2 = rd;
       } 
     momentummax = 0.0;
     momentumsum = 0.0;
     // daughter 0

     energy = rd2*(parentmass - sumofdaughtermass);
     daughtermomentum[0] = std::sqrt(energy*energy + 2.0*energy* daughtermass[0]);
     if ( daughtermomentum[0] >momentummax )momentummax =  daughtermomentum[0];
     momentumsum  +=  daughtermomentum[0];

     // daughter 1
     energy = (1.-rd1)*(parentmass - sumofdaughtermass);
     daughtermomentum[1] = std::sqrt(energy*energy + 2.0*energy* daughtermass[1]);
     if ( daughtermomentum[1] >momentummax )momentummax =  daughtermomentum[1];
     momentumsum  +=  daughtermomentum[1];

     // daughter 2
     energy = (rd1-rd2)*(parentmass - sumofdaughtermass);
     daughtermomentum[2] = std::sqrt(energy*energy + 2.0*energy* daughtermass[2]);
     if ( daughtermomentum[2] >momentummax )momentummax =  daughtermomentum[2];
     momentumsum  +=  daughtermomentum[2];
    } while ( ( momentummax >  momentumsum - momentummax ) &&  /* Loop checking, 02.11.2015, A.Ribon */ 
              ++loopCounter < maxNumberOfLoops );
    if ( loopCounter >= maxNumberOfLoops ) {
      G4ExceptionDescription ed;
      ed << " Failed sampling after maxNumberOfLoops attempts : forced exit" << G4endl;
      G4Exception( " G4GeneralPhaseSpaceDecay::ThreeBodyDecayIt ", "HAD_PHASESPACE_001", FatalException, ed );
    }

  // output message
  if (GetVerboseLevel()>1) {
    G4cout << "     daughter 0:" << daughtermomentum[0]/GeV << "[GeV/c]" <<G4endl;
    G4cout << "     daughter 1:" << daughtermomentum[1]/GeV << "[GeV/c]" <<G4endl;
    G4cout << "     daughter 2:" << daughtermomentum[2]/GeV << "[GeV/c]" <<G4endl;
    G4cout << "   momentum sum:" << momentumsum/GeV << "[GeV/c]" <<G4endl;
  }

  //create daughter G4DynamicParticle 
  G4double costheta, sintheta, phi, sinphi, cosphi; 
  G4double costhetan, sinthetan, phin, sinphin, cosphin; 
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
  phi  = twopi*G4UniformRand()*rad;
  sinphi = std::sin(phi);
  cosphi = std::cos(phi);
  G4ParticleMomentum direction0(sintheta*cosphi,sintheta*sinphi,costheta);
  G4double Etotal=std::sqrt( daughtermass[0]*daughtermass[0] + daughtermomentum[0]*daughtermomentum[0]);
  G4DynamicParticle * daughterparticle 
         = new G4DynamicParticle( G4MT_daughters[0], Etotal, direction0*daughtermomentum[0]);
  products->PushProducts(daughterparticle);

  costhetan = (daughtermomentum[1]*daughtermomentum[1]-daughtermomentum[2]*daughtermomentum[2]-daughtermomentum[0]*daughtermomentum[0])/(2.0*daughtermomentum[2]*daughtermomentum[0]);
  sinthetan = std::sqrt((1.0-costhetan)*(1.0+costhetan));
  phin  = twopi*G4UniformRand()*rad;
  sinphin = std::sin(phin);
  cosphin = std::cos(phin);
  G4ParticleMomentum direction2;
  direction2.setX( sinthetan*cosphin*costheta*cosphi - sinthetan*sinphin*sinphi + costhetan*sintheta*cosphi); 
  direction2.setY( sinthetan*cosphin*costheta*sinphi + sinthetan*sinphin*cosphi + costhetan*sintheta*sinphi); 
  direction2.setZ( -sinthetan*cosphin*sintheta + costhetan*costheta);
  Etotal=std::sqrt( daughtermass[2]*daughtermass[2] + daughtermomentum[2]*daughtermomentum[2]/direction2.mag2());
  daughterparticle = new G4DynamicParticle( G4MT_daughters[2],Etotal, direction2*(daughtermomentum[2]/direction2.mag()));
  products->PushProducts(daughterparticle);
  G4ThreeVector mom=(direction0*daughtermomentum[0] + direction2*(daughtermomentum[2]/direction2.mag()))*(-1.0);
  Etotal= std::sqrt( daughtermass[1]*daughtermass[1] + mom.mag2() );
  daughterparticle = 
       new G4DynamicParticle(G4MT_daughters[1], Etotal, mom);
  products->PushProducts(daughterparticle);

  if (GetVerboseLevel()>1) {
     G4cout << "G4GeneralPhaseSpaceDecay::ThreeBodyDecayIt ";
     G4cout << "  create decay products in rest frame " <<G4endl;
     products->DumpInfo();
  }
  return products;
}

G4DecayProducts *G4GeneralPhaseSpaceDecay::ManyBodyDecayIt()
// algorism of this code is originally written in FORTRAN by M.Asai
//*****************************************************************
//  NBODY
//   N-body phase space Monte-Carlo generator
//  Makoto Asai 
//   Hiroshima Institute of Technology
//    (asai@kekvax.kek.jp)
//  Revised release : 19/Apr/1995
//
{
  //return value
  G4DecayProducts *products;

  if (GetVerboseLevel()>1) G4cout << "G4GeneralPhaseSpaceDecay::ManyBodyDecayIt()"<<G4endl;

  //daughters'mass
  G4double *daughtermass = new G4double[numberOfDaughters]; 
  G4double sumofdaughtermass = 0.0;
  for (G4int index=0; index<numberOfDaughters; index++){
    daughtermass[index] = G4MT_daughters[index]->GetPDGMass();
    sumofdaughtermass += daughtermass[index];
  }
  
  //Calculate daughter momentum
  G4double *daughtermomentum = new G4double[numberOfDaughters];
  G4ParticleMomentum direction;  
  G4DynamicParticle **daughterparticle;
  G4double *sm = new G4double[numberOfDaughters];
  G4double tmas;
  G4double weight = 1.0;
  G4int    numberOfTry = 0; 
  G4int index1;

  do {
    //Generate rundom number in descending order 
    G4double temp;
    G4double *rd = new G4double[numberOfDaughters];
    rd[0] = 1.0;
    for(index1 =1; index1 < numberOfDaughters -1; index1++)
      rd[index1] = G4UniformRand(); 
    rd[ numberOfDaughters -1] = 0.0;
    for(index1 =1; index1 < numberOfDaughters -1; index1++) {
      for(G4int index2 = index1+1; index2 < numberOfDaughters; index2++) {
        if (rd[index1] < rd[index2]){
          temp         = rd[index1];
          rd[index1]   = rd[index2];
          rd[index2]   = temp;
        }
      }
    }
    
    //calcurate virtual mass 
    tmas = parentmass -  sumofdaughtermass;
    temp =  sumofdaughtermass; 
    for(index1 =0; index1 < numberOfDaughters; index1++) {
      sm[index1] = rd[index1]*tmas + temp;
      temp -= daughtermass[index1];
      if (GetVerboseLevel()>1) {
        G4cout << index1 << "  rundom number:" << rd[index1]; 
        G4cout << "   virtual mass:" << sm[index1]/GeV << "[GeV/c/c]" <<G4endl; 
      }
    }
    delete [] rd;

    //Calculate daughter momentum
    weight = 1.0;
    index1 =numberOfDaughters-1;
    daughtermomentum[index1]= Pmx( sm[index1-1],daughtermass[index1-1],sm[index1]);
    if (GetVerboseLevel()>1) {
      G4cout << "     daughter " << index1 << ":" << *daughters_name[index1];
      G4cout << " momentum:" << daughtermomentum[index1]/GeV << "[GeV/c]" <<G4endl;
    }
    for(index1 =numberOfDaughters-2; index1>=0; index1--) {
      // calculate 
      daughtermomentum[index1]= Pmx( sm[index1],daughtermass[index1], sm[index1 +1]);
      if(daughtermomentum[index1] < 0.0) {
        // !!! illegal momentum !!!
        if (GetVerboseLevel()>0) {
          G4cout << "G4GeneralPhaseSpaceDecay::ManyBodyDecayIt ";
          G4cout << "     can not calculate daughter momentum " <<G4endl;
          G4cout << "     parent:" << *parent_name;
          G4cout << " mass:" << parentmass/GeV << "[GeV/c/c]" <<G4endl;
          G4cout << "     daughter " << index1 << ":" << *daughters_name[index1];
          G4cout << " mass:" << daughtermass[index1]/GeV << "[GeV/c/c]" ;
          G4cout << " mass:" << daughtermomentum[index1]/GeV << "[GeV/c]" <<G4endl;
        }
	delete [] sm;
	delete [] daughtermass;
	delete [] daughtermomentum;
	return NULL;   // Error detection

      } else {
	// calculate weight of this events
        weight *=  daughtermomentum[index1]/sm[index1];
        if (GetVerboseLevel()>1) {
          G4cout << "     daughter " << index1 << ":" << *daughters_name[index1];
          G4cout << " momentum:" << daughtermomentum[index1]/GeV << "[GeV/c]" <<G4endl;
        }
      }
    }
    if (GetVerboseLevel()>1) {
      G4cout << "    weight: " << weight <<G4endl;
    }
    
    // exit if number of Try exceeds 100
    if (numberOfTry++ >100) {
      if (GetVerboseLevel()>0) {
        G4cout << "G4GeneralPhaseSpaceDecay::ManyBodyDecayIt: ";
	G4cout << " can not determine Decay Kinematics " << G4endl;
      }
      delete [] sm;
      delete [] daughtermass;
      delete [] daughtermomentum;
      return NULL;  // Error detection
    }
  } while ( weight > G4UniformRand());  /* Loop checking, 02.11.2015, A.Ribon */
  if (GetVerboseLevel()>1) {
      G4cout << "Start calulation of daughters momentum vector "<<G4endl;
  }
  
  G4double costheta, sintheta, phi;
  G4double beta;
  daughterparticle = new G4DynamicParticle*[numberOfDaughters];

  index1 = numberOfDaughters -2;
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
  phi  = twopi*G4UniformRand()*rad;
  direction.setZ(costheta);
  direction.setY(sintheta*std::sin(phi));
  direction.setX(sintheta*std::cos(phi));
  daughterparticle[index1] = new G4DynamicParticle( G4MT_daughters[index1], direction*daughtermomentum[index1] );
  daughterparticle[index1+1] = new G4DynamicParticle( G4MT_daughters[index1+1], direction*(-1.0*daughtermomentum[index1]) );

  for (index1 = numberOfDaughters -3;  index1 >= 0; index1--) {
    //calculate momentum direction
    costheta = 2.*G4UniformRand()-1.0;
    sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
    phi  = twopi*G4UniformRand()*rad;
    direction.setZ(costheta);
    direction.setY(sintheta*std::sin(phi));
    direction.setX(sintheta*std::cos(phi));

    // boost already created particles 
    beta = daughtermomentum[index1];
    beta /= std::sqrt( daughtermomentum[index1]*daughtermomentum[index1] + sm[index1+1]*sm[index1+1] );
    for (G4int index2 = index1+1; index2<numberOfDaughters; index2++) {
      G4LorentzVector p4;
      // make G4LorentzVector for secondaries
      p4 = daughterparticle[index2]->Get4Momentum();

      // boost secondaries to  new frame 
      p4.boost( direction.x()*beta, direction.y()*beta, direction.z()*beta);

      // change energy/momentum
      daughterparticle[index2]->Set4Momentum(p4);
    }
    //create daughter G4DynamicParticle 
    daughterparticle[index1]= new G4DynamicParticle( G4MT_daughters[index1], direction*(-1.0*daughtermomentum[index1]));
  }

  //create G4Decayproducts
  G4DynamicParticle *parentparticle; 
  direction.setX(1.0);  direction.setY(0.0); direction.setZ(0.0);
  parentparticle = new G4DynamicParticle( G4MT_parent, direction, 0.0);
  products = new G4DecayProducts(*parentparticle);
  delete parentparticle;
  for (index1 = 0; index1<numberOfDaughters; index1++) {
    products->PushProducts(daughterparticle[index1]);
  }
  if (GetVerboseLevel()>1) { 
    G4cout << "G4GeneralPhaseSpaceDecay::ManyBodyDecayIt ";
    G4cout << "  create decay products in rest frame " << G4endl;
    products->DumpInfo();
  }

  delete [] daughterparticle;
  delete [] daughtermomentum;
  delete [] daughtermass;
  delete [] sm;
  
  return products;
}





