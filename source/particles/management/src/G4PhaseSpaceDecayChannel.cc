// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PhaseSpaceDecayChannel.cc,v 1.1 1999-01-07 16:10:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      27 July 1996 H.Kurashige
//      20 Oct 1996 H.Kurashige
//      30 May 1997 H.Kurashige
// ------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4DecayProducts.hh"
#include "G4VDecayChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"


G4PhaseSpaceDecayChannel::G4PhaseSpaceDecayChannel(G4int Verbose)
 :G4VDecayChannel("Phase Space", Verbose)
{

}

G4PhaseSpaceDecayChannel::G4PhaseSpaceDecayChannel(
			       const G4String& theParentName,
			       G4double        theBR,
			       G4int           theNumberOfDaughters,
			       const G4String& theDaughterName1,
			       const G4String& theDaughterName2,
			       const G4String& theDaughterName3,
			       const G4String& theDaughterName4    )
                         :G4VDecayChannel("Phase Space",
					  theParentName,theBR,
					  theNumberOfDaughters,
					  theDaughterName1,
					  theDaughterName2,
					  theDaughterName3,
					  theDaughterName4)
{

}

G4PhaseSpaceDecayChannel::~G4PhaseSpaceDecayChannel()
{
}

G4DecayProducts *G4PhaseSpaceDecayChannel::DecayIt(G4double parentMass)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4PhaseSpaceDecayChannel::DecayIt ";
#endif
  
  G4DecayProducts * products = NULL;
 
  if (parent == NULL) FillParent();  
  if (daughters == NULL) FillDaughters();

  if (parentMass >0.0) parent_mass = parentMass;  
  switch (numberOfDaughters){
  case 0:
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4PhaseSpaceDecayChannel::DecayIt ";
      G4cout << " daughters not defined " <<endl;
    }
#endif
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
#ifdef G4VERBOSE
  if ((products == NULL) && (GetVerboseLevel()>0)) {
    G4cout << "G4PhaseSpaceDecayChannel::DecayIt ";
    G4cout << *parent_name << " can not decay " << endl;
    DumpInfo();
  }
#endif
  return products;
}

G4DecayProducts *G4PhaseSpaceDecayChannel::OneBodyDecayIt()
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4PhaseSpaceDecayChannel::OneBodyDecayIt()"<<endl;
#endif
  G4double parentmass = parent_mass;
  G4double daughtermass = daughters_mass[0];

  //create parent G4DynamicParticle at rest
  G4ParticleMomentum dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( parent, dummy, 0.0);
  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  //create daughter G4DynamicParticle at rest
  G4DynamicParticle * daughterparticle = new G4DynamicParticle( daughters[0], dummy, 0.0);
  products->PushProducts(daughterparticle);

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
     G4cout << "G4PhaseSpaceDecayChannel::OneBodyDecayIt ";
     G4cout << "  create decay products in rest frame " <<endl;
     products->DumpInfo();
  }
#endif
  return products;
}

G4DecayProducts *G4PhaseSpaceDecayChannel::TwoBodyDecayIt()
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4PhaseSpaceDecayChannel::TwoBodyDecayIt()"<<endl;
#endif
  // parent mass
  G4double parentmass = parent_mass;
  
  //daughters'mass
  G4double daughtermass[2]; 
  G4double daughtermomentum;
  daughtermass[0] = daughters_mass[0];
  daughtermass[1] = daughters_mass[1];
  G4double sumofdaughtermass =  daughtermass[0] + daughtermass[1];

  //create parent G4DynamicParticle at rest
  G4ParticleMomentum dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( parent, dummy, 0.0);
  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  //calculate daughter momentum
  daughtermomentum = Pmx(parentmass,daughtermass[0],daughtermass[1]);
  G4double costheta = 2.*G4UniformRand()-1.0;
  G4double sintheta = sqrt((1.0 - costheta)*(1.0 + costheta));
  G4double phi  = 2.0*M_PI*G4UniformRand()*rad;
  G4ParticleMomentum direction(sintheta*cos(phi),sintheta*sin(phi),costheta);

  //create daughter G4DynamicParticle 
  G4DynamicParticle * daughterparticle = new G4DynamicParticle( daughters[0], direction*daughtermomentum);
  products->PushProducts(daughterparticle);
  daughterparticle = new G4DynamicParticle( daughters[1], direction*(-1.0*daughtermomentum));
  products->PushProducts(daughterparticle);

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
     G4cout << "G4PhaseSpaceDecayChannel::TwoBodyDecayIt ";
     G4cout << "  create decay products in rest frame " <<endl;
     products->DumpInfo();
  }
#endif
  return products;
}

G4DecayProducts *G4PhaseSpaceDecayChannel::ThreeBodyDecayIt()
// algorism of this code is originally written in GDECA3 of GEANT3
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4PhaseSpaceDecayChannel::ThreeBodyDecayIt()"<<endl;
#endif
  // parent mass
  G4double parentmass = parent_mass;
  //daughters'mass
  G4double daughtermass[3]; 
  G4double sumofdaughtermass = 0.0;
  for (G4int index=0; index<3; index++){
    daughtermass[index] = daughters_mass[index];
    sumofdaughtermass += daughtermass[index];
  }
  
   //create parent G4DynamicParticle at rest
  G4ParticleMomentum dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( parent, dummy, 0.0);
  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  //calculate daughter momentum
  //  Generate two 
  G4double rd1, rd2, rd;
  G4double daughtermomentum[3];
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
    energy = rd2*(parentmass - sumofdaughtermass);
    daughtermomentum[0] = sqrt(energy*energy + 2.0*energy* daughtermass[0]);
    if ( daughtermomentum[0] >momentummax )momentummax =  daughtermomentum[0];
    momentumsum  +=  daughtermomentum[0];
    // daughter 1
    energy = (1.-rd1)*(parentmass - sumofdaughtermass);
    daughtermomentum[1] = sqrt(energy*energy + 2.0*energy* daughtermass[1]);
    if ( daughtermomentum[1] >momentummax )momentummax =  daughtermomentum[1];
    momentumsum  +=  daughtermomentum[1];
    // daughter 2
    energy = (rd1-rd2)*(parentmass - sumofdaughtermass);
    daughtermomentum[2] = sqrt(energy*energy + 2.0*energy* daughtermass[2]);
    if ( daughtermomentum[2] >momentummax )momentummax =  daughtermomentum[2];
    momentumsum  +=  daughtermomentum[2];
  } while (momentummax >  momentumsum - momentummax );
  // output message
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "     daughter 0:" << daughtermomentum[0]/GeV << "[GeV/c]" <<endl;
    G4cout << "     daughter 1:" << daughtermomentum[1]/GeV << "[GeV/c]" <<endl;
    G4cout << "     daughter 2:" << daughtermomentum[2]/GeV << "[GeV/c]" <<endl;
    G4cout << "   momentum sum:" << momentumsum/GeV << "[GeV/c]" <<endl;
  }
#endif

  //create daughter G4DynamicParticle 
  G4double costheta, sintheta, phi, sinphi, cosphi; 
  G4double costhetan, sinthetan, phin, sinphin, cosphin; 
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = sqrt((1.0-costheta)*(1.0+costheta));
  phi  = 2.0*M_PI*G4UniformRand()*rad;
  sinphi = sin(phi);
  cosphi = cos(phi);
  G4ParticleMomentum direction0(sintheta*cosphi,sintheta*sinphi,costheta);
  G4DynamicParticle * daughterparticle 
         = new G4DynamicParticle( daughters[0], direction0*daughtermomentum[0]);
  products->PushProducts(daughterparticle);

  costhetan = (daughtermomentum[1]*daughtermomentum[1]-daughtermomentum[2]*daughtermomentum[2]-daughtermomentum[0]*daughtermomentum[0])/(2.0*daughtermomentum[2]*daughtermomentum[0]);
  sinthetan = sqrt((1.0-costhetan)*(1.0+costhetan));
  phin  = 2.0*M_PI*G4UniformRand()*rad;
  sinphin = sin(phin);
  cosphin = cos(phin);
  G4ParticleMomentum direction2;
  direction2.setX( sinthetan*cosphin*costheta*cosphi - sinthetan*sinphin*sinphi + costhetan*sintheta*cosphi); 
  direction2.setY( sinthetan*cosphin*costheta*sinphi + sinthetan*sinphin*cosphi + costhetan*sintheta*sinphi); 
  direction2.setZ( -sinthetan*cosphin*sintheta + costhetan*costheta);
  daughterparticle = new G4DynamicParticle( daughters[2], direction2*(daughtermomentum[2]/direction2.mag()));
  products->PushProducts(daughterparticle);

  daughterparticle = 
       new G4DynamicParticle( 
	        daughters[1],
	       (direction0*daughtermomentum[0] + direction2*(daughtermomentum[2]/direction2.mag()))*(-1.0)   
		);
  products->PushProducts(daughterparticle);

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
     G4cout << "G4PhaseSpaceDecayChannel::ThreeBodyDecayIt ";
     G4cout << "  create decay products in rest frame " <<endl;
     products->DumpInfo();
  }
#endif
  return products;
}

G4DecayProducts *G4PhaseSpaceDecayChannel::ManyBodyDecayIt()
// algorithm of this code is originally written in FORTRAN by M.Asai
//******************************************************************
//  NBODY
//   N-body phase space Monte-Carlo generator
//  Makoto Asai 
//   Hiroshima Institute of Technology
//    (asai@kekvax.kek.jp)
//  Revised release : 19/Apr/1995
//
{
  G4int index, index2;

  //return value
  G4DecayProducts *products;

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4PhaseSpaceDecayChannel::ManyBodyDecayIt()"<<endl;
#endif
  // parent mass
  G4double parentmass = parent_mass;
  //daughters'mass
  G4double *daughtermass = new G4double[numberOfDaughters]; 
  G4double sumofdaughtermass = 0.0;
  for (index=0; index<numberOfDaughters; index++){
    daughtermass[index] = daughters_mass[index];
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
 
  do {
    //Generate rundom number in descending order 
    G4double temp;
    G4double *rd = new G4double[numberOfDaughters];
    rd[0] = 1.0;
    for(index =1; index < numberOfDaughters -1; index++) rd[index] = G4UniformRand(); 
    rd[ numberOfDaughters -1] = 0.0;
    for(index =1; index < numberOfDaughters -1; index++) {
      for(index2 = index+1; index2 < numberOfDaughters; index2++) {
        if (rd[index] < rd[index2]){
          temp        = rd[index];
          rd[index]  = rd[index2];
          rd[index2] = temp;
        }
      }
    }
    
    //calcurate virtual mass 
    tmas = parentmass -  sumofdaughtermass;
    temp =  sumofdaughtermass; 
    for(index =0; index < numberOfDaughters; index++) {
      sm[index] = rd[index]*tmas + temp;
      temp -= daughtermass[index];
      if (GetVerboseLevel()>1) {
        G4cout << index << "  rundom number:" << rd[index]; 
        G4cout << "   virtual mass:" << sm[index]/GeV << "[GeV/c/c]" <<endl; 
      }
    }
    delete [] rd;

    //Calculate daughter momentum
    weight = 1.0;
    index =numberOfDaughters-1;
    daughtermomentum[index]= Pmx( sm[index-1],daughtermass[index-1], sm[index]);
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cout << "     daughter " << index << ":" << *daughters_name[index];
      G4cout << " momentum:" << daughtermomentum[index]/GeV << "[GeV/c]" <<endl;
    }
#endif
    for(index =numberOfDaughters-2; index>=0; index--) {
      // calculate 
      daughtermomentum[index]= Pmx( sm[index],daughtermass[index], sm[index +1]);
      if(daughtermomentum[index] < 0.0) {
        // !!! illegal momentum !!!
#ifdef G4VERBOSE
        if (GetVerboseLevel()>0) {
          G4cout << "G4PhaseSpaceDecayChannel::ManyBodyDecayIt ";
          G4cout << "     can not calculate daughter momentum " <<endl;
          G4cout << "     parent:" << *parent_name;
          G4cout << " mass:" << parentmass/GeV << "[GeV/c/c]" <<endl;
          G4cout << "     daughter " << index << ":" << *daughters_name[index];
          G4cout << " mass:" << daughtermass[index]/GeV << "[GeV/c/c]" ;
          G4cout << " mass:" << daughtermomentum[index]/GeV << "[GeV/c]" <<endl;
        }
#endif
	delete [] sm;
	delete [] daughtermass;
	delete [] daughtermomentum;
	return NULL;   // Error detection

      } else {
	// calculate weight of this events
        weight *=  daughtermomentum[index]/sm[index];
        if (GetVerboseLevel()>1) {
          G4cout << "     daughter " << index << ":" << *daughters_name[index];
          G4cout << " momentum:" << daughtermomentum[index]/GeV << "[GeV/c]" <<endl;
        }
      }
    }
    if (GetVerboseLevel()>1) {
      G4cout << "    weight: " << weight <<endl;
    }
    
    // exit if number of Try exceeds 100
    if (numberOfTry++ >100) {
      if (GetVerboseLevel()>0) {
        G4cout << "G4PhaseSpaceDecayChannel::ManyBodyDecayIt: ";
	G4cout << " can not determine Decay Kinematics " << endl;
      }
      delete [] sm;
      delete [] daughtermass;
      delete [] daughtermomentum;
      return NULL;  // Error detection
    }
  } while ( weight > G4UniformRand());

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
      G4cout << "Start calulation of daughters momentum vector "<<endl;
  }
#endif  
  G4double costheta, sintheta, phi;
  G4double beta;
  daughterparticle = new G4DynamicParticle*[numberOfDaughters];

  index = numberOfDaughters -2;
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = sqrt((1.0-costheta)*(1.0+costheta));
  phi  = 2.0*M_PI*G4UniformRand()*rad;
  direction.setZ(costheta);
  direction.setY(sintheta*sin(phi));
  direction.setX(sintheta*cos(phi));
  daughterparticle[index] = new G4DynamicParticle( daughters[index], direction*daughtermomentum[index] );
  daughterparticle[index+1] = new G4DynamicParticle( daughters[index+1], direction*(-1.0*daughtermomentum[index]) );

  for (index = numberOfDaughters -3;  index >= 0; index--) {
    //calculate momentum direction
    costheta = 2.*G4UniformRand()-1.0;
    sintheta = sqrt((1.0-costheta)*(1.0+costheta));
    phi  = 2.0*M_PI*G4UniformRand()*rad;
    direction.setZ(costheta);
    direction.setY(sintheta*sin(phi));
    direction.setX(sintheta*cos(phi));

    // boost already created particles 
    beta = daughtermomentum[index];
    beta /= sqrt( daughtermomentum[index]*daughtermomentum[index] + sm[index+1]*sm[index+1] );
    for (G4int index2 = index+1; index2<numberOfDaughters; index2++) {
      G4LorentzVector p4;
      // make G4LorentzVector for secondaries
      p4 = daughterparticle[index2]->Get4Momentum();

      // boost secondaries to  new frame 
      p4.boost( direction.x()*beta, direction.y()*beta, direction.z()*beta);

      // change energy/momentum
      daughterparticle[index2]->Set4Momentum(p4);
    }
    //create daughter G4DynamicParticle 
    daughterparticle[index]= new G4DynamicParticle( daughters[index], direction*(-1.0*daughtermomentum[index]));
  }

  //create G4Decayproducts
  G4DynamicParticle *parentparticle; 
  direction.setX(1.0);  direction.setY(0.0); direction.setZ(0.0);
  parentparticle = new G4DynamicParticle( parent, direction, 0.0);
  products = new G4DecayProducts(*parentparticle);
  delete parentparticle;
  for (index = 0; index<numberOfDaughters; index++) {
    products->PushProducts(daughterparticle[index]);
  }
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) { 
    G4cout << "G4PhaseSpaceDecayChannel::ManyBodyDecayIt ";
    G4cout << "  create decay products in rest frame " << endl;
    products->DumpInfo();
  }
#endif
  delete [] daughterparticle;
  delete [] daughtermomentum;
  delete [] daughtermass;
  delete [] sm;
  
  return products;
}





