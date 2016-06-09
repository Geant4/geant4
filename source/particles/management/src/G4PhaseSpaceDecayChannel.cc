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
// $Id$
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      History: first implementation, based on object model of
//      27 July 1996 H.Kurashige
//      20 Oct 1996 H.Kurashige
//      30 May 1997 H.Kurashige
// ------------------------------------------------------------

#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
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
					  theDaughterName4),
			  current_parent_mass(0.0)
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
  
  G4DecayProducts * products = 0;
 
  if (parent == 0) FillParent();  
  if (daughters == 0) FillDaughters();

  if (parentMass >0.0) current_parent_mass = parentMass;
  else                 current_parent_mass = parent_mass;

  switch (numberOfDaughters){
  case 0:
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4PhaseSpaceDecayChannel::DecayIt ";
      G4cout << " daughters not defined " <<G4endl;
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
  if ((products == 0) && (GetVerboseLevel()>0)) {
    G4cout << "G4PhaseSpaceDecayChannel::DecayIt ";
    G4cout << *parent_name << " can not decay " << G4endl;
    DumpInfo();
  }
#endif
  return products;
}

G4DecayProducts *G4PhaseSpaceDecayChannel::OneBodyDecayIt()
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4PhaseSpaceDecayChannel::OneBodyDecayIt()"<<G4endl;
#endif

  //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
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
     G4cout << "  create decay products in rest frame " <<G4endl;
     products->DumpInfo();
  }
#endif
  return products;
}

G4DecayProducts *G4PhaseSpaceDecayChannel::TwoBodyDecayIt()
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4PhaseSpaceDecayChannel::TwoBodyDecayIt()"<<G4endl;
#endif
  // parent mass
  G4double parentmass = current_parent_mass;
  
  //daughters'mass
  G4double daughtermass[2]; 
  G4double daughtermomentum;
  daughtermass[0] = daughters_mass[0];
  daughtermass[1] = daughters_mass[1];

  //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( parent, dummy, 0.0);
  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  //calculate daughter momentum
  daughtermomentum = Pmx(parentmass,daughtermass[0],daughtermass[1]);
  if (daughtermomentum <0.0) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << "G4PhaseSpaceDecayChannel::TwoBodyDecayIt " 
             << "sum of daughter mass is larger than parent mass" << G4endl;
      G4cerr << "parent :" << parent->GetParticleName() << "  " << current_parent_mass/GeV << G4endl;
      G4cerr << "daughter 1 :" << daughters[0]->GetParticleName() << "  " << daughtermass[0]/GeV << G4endl;
      G4cerr << "daughter 2:" << daughters[1]->GetParticleName() << "  " << daughtermass[1]/GeV << G4endl;
    }
#endif
    G4Exception("G4PhaseSpaceDecayChannel::TwoBodyDecayIt",
                "PART112", JustWarning,
                "Can not create decay products: sum of daughter mass is larger than parent mass");
    return products;
  }

  G4double costheta = 2.*G4UniformRand()-1.0;
  G4double sintheta = std::sqrt((1.0 - costheta)*(1.0 + costheta));
  G4double phi  = twopi*G4UniformRand()*rad;
  G4ThreeVector direction(sintheta*std::cos(phi),sintheta*std::sin(phi),costheta);

  //create daughter G4DynamicParticle 
  G4DynamicParticle * daughterparticle = new G4DynamicParticle( daughters[0], direction*daughtermomentum);
  products->PushProducts(daughterparticle);
  daughterparticle = new G4DynamicParticle( daughters[1], direction*(-1.0*daughtermomentum));
  products->PushProducts(daughterparticle);

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
     G4cout << "G4PhaseSpaceDecayChannel::TwoBodyDecayIt ";
     G4cout << "  create decay products in rest frame " <<G4endl;
     products->DumpInfo();
  }
#endif
  return products;
}

G4DecayProducts *G4PhaseSpaceDecayChannel::ThreeBodyDecayIt()
// algorism of this code is originally written in GDECA3 of GEANT3
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4PhaseSpaceDecayChannel::ThreeBodyDecayIt()"<<G4endl;
#endif
  // parent mass
  G4double parentmass = current_parent_mass;
  //daughters'mass
  G4double daughtermass[3]; 
  G4double sumofdaughtermass = 0.0;
  for (G4int index=0; index<3; index++){
    daughtermass[index] = daughters_mass[index];
    sumofdaughtermass += daughtermass[index];
  }
  
   //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( parent, dummy, 0.0);


  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  if (sumofdaughtermass >parentmass) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cerr << "G4PhaseSpaceDecayChannel::ThreeBodyDecayIt " 
             << "sum of daughter mass is larger than parent mass" << G4endl;
      G4cerr << "parent :" << parent->GetParticleName() << "  " << current_parent_mass/GeV << G4endl;
      G4cerr << "daughter 1 :" << daughters[0]->GetParticleName() << "  " << daughtermass[0]/GeV << G4endl;
      G4cerr << "daughter 2:" << daughters[1]->GetParticleName() << "  " << daughtermass[1]/GeV << G4endl;
      G4cerr << "daughter 3:" << daughters[2]->GetParticleName() << "  " << daughtermass[2]/GeV << G4endl;
    }
#endif
    G4Exception("G4PhaseSpaceDecayChannel::ThreeBodyDecayIt",
                "PART112", JustWarning,
                "Can not create decay products: sum of daughter mass is larger than parent mass");
    return products;
  }



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
  } while (momentummax >  momentumsum - momentummax );
  // output message
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "     daughter 0:" << daughtermomentum[0]/GeV << "[GeV/c]" <<G4endl;
    G4cout << "     daughter 1:" << daughtermomentum[1]/GeV << "[GeV/c]" <<G4endl;
    G4cout << "     daughter 2:" << daughtermomentum[2]/GeV << "[GeV/c]" <<G4endl;
    G4cout << "   momentum sum:" << momentumsum/GeV << "[GeV/c]" <<G4endl;
  }
#endif

  //create daughter G4DynamicParticle 
  G4double costheta, sintheta, phi, sinphi, cosphi; 
  G4double costhetan, sinthetan, phin, sinphin, cosphin; 
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
  phi  = twopi*G4UniformRand()*rad;
  sinphi = std::sin(phi);
  cosphi = std::cos(phi);
  G4ThreeVector direction0(sintheta*cosphi,sintheta*sinphi,costheta);
  G4DynamicParticle * daughterparticle 
         = new G4DynamicParticle( daughters[0], direction0*daughtermomentum[0]);
  products->PushProducts(daughterparticle);

  costhetan = (daughtermomentum[1]*daughtermomentum[1]-daughtermomentum[2]*daughtermomentum[2]-daughtermomentum[0]*daughtermomentum[0])/(2.0*daughtermomentum[2]*daughtermomentum[0]);
  sinthetan = std::sqrt((1.0-costhetan)*(1.0+costhetan));
  phin  = twopi*G4UniformRand()*rad;
  sinphin = std::sin(phin);
  cosphin = std::cos(phin);
  G4ThreeVector direction2;
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
     G4cout << "  create decay products in rest frame " <<G4endl;
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
  if (GetVerboseLevel()>1) G4cout << "G4PhaseSpaceDecayChannel::ManyBodyDecayIt()"<<G4endl;
#endif
  // parent mass
  G4double parentmass = current_parent_mass;
  //daughters'mass
  G4double *daughtermass = new G4double[numberOfDaughters]; 
  G4double sumofdaughtermass = 0.0;
  for (index=0; index<numberOfDaughters; index++){
    daughtermass[index] = daughters_mass[index];
    sumofdaughtermass += daughtermass[index];
  }
  
  //Calculate daughter momentum
  G4double *daughtermomentum = new G4double[numberOfDaughters];
  G4ThreeVector direction;  
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
        G4cout << "   virtual mass:" << sm[index]/GeV << "[GeV/c/c]" <<G4endl; 
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
      G4cout << " momentum:" << daughtermomentum[index]/GeV << "[GeV/c]" <<G4endl;
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
          G4cout << "     can not calculate daughter momentum " <<G4endl;
          G4cout << "     parent:" << *parent_name;
          G4cout << " mass:" << parentmass/GeV << "[GeV/c/c]" <<G4endl;
          G4cout << "     daughter " << index << ":" << *daughters_name[index];
          G4cout << " mass:" << daughtermass[index]/GeV << "[GeV/c/c]" ;
          G4cout << " mass:" << daughtermomentum[index]/GeV << "[GeV/c]" <<G4endl;
        }
#endif
	delete [] sm;
	delete [] daughtermass;
	delete [] daughtermomentum;

	G4Exception("G4PhaseSpaceDecayChannel::ManyBodyDecayIt",
                "PART112", JustWarning,
                "Can not create decay products: sum of daughter mass is larger than parent mass");

	return 0;   // Error detection

      } else {
	// calculate weight of this events
        weight *=  daughtermomentum[index]/sm[index];
#ifdef G4VERBOSE
        if (GetVerboseLevel()>1) {
          G4cout << "     daughter " << index << ":" << *daughters_name[index];
          G4cout << " momentum:" << daughtermomentum[index]/GeV << "[GeV/c]" <<G4endl;
        }
#endif
      }
    }

#ifdef G4VERBOSE
    if (GetVerboseLevel()>1) {
      G4cout << "    weight: " << weight <<G4endl;
    }
#endif
    
    // exit if number of Try exceeds 100
    if (numberOfTry++ >100) {
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0) {
        G4cout << "G4PhaseSpaceDecayChannel::ManyBodyDecayIt: ";
	G4cout << " can not determine Decay Kinematics " << G4endl;
	G4cout << "parent : " << *parent_name << G4endl;
	G4cout << "daughters : "; 
	for(index=0; index<numberOfDaughters; index++) {
	  G4cout <<  *daughters_name[index] << " , ";
	}
	G4cout << G4endl;
      }
#endif
      G4Exception("G4PhaseSpaceDecayChannel::ManyBodyDecayIt: ",
                  "PART113", JustWarning,
		  " Cannot decay :  Decay Kinematics cannot be calculated ");
      
      delete [] sm;
      delete [] daughtermass;
      delete [] daughtermomentum;
      return 0;  // Error detection
    }
  } while ( weight > G4UniformRand());

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
      G4cout << "Start calulation of daughters momentum vector "<<G4endl;
  }
#endif 
 
  G4double costheta, sintheta, phi;
  G4double beta;
  daughterparticle = new G4DynamicParticle*[numberOfDaughters];

  index = numberOfDaughters -2;
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
  phi  = twopi*G4UniformRand()*rad;
  direction.setZ(costheta);
  direction.setY(sintheta*std::sin(phi));
  direction.setX(sintheta*std::cos(phi));
  daughterparticle[index] = new G4DynamicParticle( daughters[index], direction*daughtermomentum[index] );
  daughterparticle[index+1] = new G4DynamicParticle( daughters[index+1], direction*(-1.0*daughtermomentum[index]) );

  for (index = numberOfDaughters -3;  index >= 0; index--) {
    //calculate momentum direction
    costheta = 2.*G4UniformRand()-1.0;
    sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
    phi  = twopi*G4UniformRand()*rad;
    direction.setZ(costheta);
    direction.setY(sintheta*std::sin(phi));
    direction.setX(sintheta*std::cos(phi));

    // boost already created particles 
    beta = daughtermomentum[index];
    beta /= std::sqrt( daughtermomentum[index]*daughtermomentum[index] + sm[index+1]*sm[index+1] );
    for (index2 = index+1; index2<numberOfDaughters; index2++) {
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
    G4cout << "  create decay products in rest frame " << G4endl;
    products->DumpInfo();
  }
#endif

  delete [] daughterparticle;
  delete [] daughtermomentum;
  delete [] daughtermass;
  delete [] sm;
  
  return products;
}


G4double G4PhaseSpaceDecayChannel::Pmx(G4double e, G4double p1, G4double p2)
{
   // calcurate momentum of daughter particles in two-body decay
   G4double ppp = (e+p1+p2)*(e+p1-p2)*(e-p1+p2)*(e-p1-p2)/(4.0*e*e);
   if (ppp>0) return std::sqrt(ppp);
   else       return -1.;
}



