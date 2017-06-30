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
// $Id: G4PhaseSpaceDecayChannel.cc 102916 2017-03-02 12:58:58Z gcosmo $
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

//G4ThreadLocal G4double G4PhaseSpaceDecayChannel::current_parent_mass = 0.0;

G4PhaseSpaceDecayChannel::G4PhaseSpaceDecayChannel(G4int Verbose)
  :G4VDecayChannel("Phase Space", Verbose),
   useGivenDaughterMass(false)
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
			  useGivenDaughterMass(false)
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
 
  CheckAndFillParent();
  CheckAndFillDaughters();

  if (parentMass >0.0) current_parent_mass.Put( parentMass );
  else                 current_parent_mass.Put(G4MT_parent_mass);

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
  // parent mass
  G4double parentmass = current_parent_mass.Get();

  //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( G4MT_parent, dummy, 0.0, parentmass);
  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  //create daughter G4DynamicParticle at rest
  G4DynamicParticle * daughterparticle = new G4DynamicParticle( G4MT_daughters[0], dummy, 0.0);
  if (useGivenDaughterMass) daughterparticle->SetMass(givenDaughterMasses[0]);
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
  G4double parentmass = current_parent_mass.Get();
  
  //daughters'mass, width
  G4double daughtermass[2], daughterwidth[2]; 
  daughtermass[0] = G4MT_daughters_mass[0];
  daughtermass[1] = G4MT_daughters_mass[1];
  daughterwidth[0] = G4MT_daughters_width[0];
  daughterwidth[1] = G4MT_daughters_width[1];

  //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( G4MT_parent, dummy, 0.0, parentmass);
  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  if (!useGivenDaughterMass) {
    G4bool withWidth =   (daughterwidth[0]>1.0e-3*daughtermass[0]) 
                      || (daughterwidth[1]>1.0e-3*daughtermass[1]);
    if (withWidth) {
      G4double sumofdaughterwidthsq = daughterwidth[0]*daughterwidth[0]+daughterwidth[1]*daughterwidth[1];
      // check parent mass and daughter mass
      G4double maxDev = (parentmass - daughtermass[0] - daughtermass[1] )/std::sqrt(sumofdaughterwidthsq);
      if (maxDev <= -1.0*rangeMass ){
#ifdef G4VERBOSE
	if (GetVerboseLevel()>0) {
	  G4cout << "G4PhaseSpaceDecayChannel::TwoBodyDecayIt " 
		 << "sum of daughter mass is larger than parent mass" << G4endl;
	  G4cout << "parent :" << G4MT_parent->GetParticleName() << "  " << current_parent_mass.Get()/GeV << G4endl;
	  G4cout << "daughter 1 :" << G4MT_daughters[0]->GetParticleName() << "  " << daughtermass[0]/GeV << G4endl;
	  G4cout << "daughter 2:" << G4MT_daughters[1]->GetParticleName() << "  " << daughtermass[1]/GeV << G4endl;
	}
#endif
	G4Exception("G4PhaseSpaceDecayChannel::TwoBodyDecayIt",
		    "PART112", JustWarning,
		    "Can not create decay products: sum of daughter mass is larger than parent mass");
	return products;
      }
      G4double dm1=daughtermass[0];
      if (daughterwidth[0] > 0.) dm1= DynamicalMass(daughtermass[0],daughterwidth[0], maxDev);
      G4double dm2= daughtermass[1];
      if (daughterwidth[1] > 0.) dm2= DynamicalMass(daughtermass[1],daughterwidth[1], maxDev);
      while (dm1+dm2>parentmass){ // Loop checking, 09.08.2015, K.Kurashige
	dm1= DynamicalMass(daughtermass[0],daughterwidth[0], maxDev);
	dm2= DynamicalMass(daughtermass[1],daughterwidth[1], maxDev);
      }
      daughtermass[0] = dm1;
      daughtermass[1] = dm2;
    }
  } else {
    // use given daughter mass;
    daughtermass[0] = givenDaughterMasses[0];
    daughtermass[1] = givenDaughterMasses[1];
  }
  if (parentmass < daughtermass[0] + daughtermass[1] ){
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4PhaseSpaceDecayChannel::TwoBodyDecayIt " 
	     << "sum of daughter mass is larger than parent mass" << G4endl;
      G4cout << "parent :" << G4MT_parent->GetParticleName() << "  " << current_parent_mass.Get()/GeV << G4endl;
      G4cout << "daughter 1 :" << G4MT_daughters[0]->GetParticleName() << "  " << daughtermass[0]/GeV << G4endl;
      G4cout << "daughter 2:" << G4MT_daughters[1]->GetParticleName() << "  " << daughtermass[1]/GeV << G4endl;
      if (useGivenDaughterMass) {
	G4cout << "Daughter Mass is given " << G4endl;
      }
    }
#endif
    G4Exception("G4PhaseSpaceDecayChannel::TwoBodyDecayIt",
		"PART112", JustWarning,
		"Can not create decay products: sum of daughter mass is larger than parent mass");
    return products;
  }  
  
  //calculate daughter momentum
  G4double daughtermomentum = Pmx(parentmass,daughtermass[0],daughtermass[1]);

  G4double costheta = 2.*G4UniformRand()-1.0;
  G4double sintheta = std::sqrt((1.0 - costheta)*(1.0 + costheta));
  G4double phi  = twopi*G4UniformRand()*rad;
  G4ThreeVector direction(sintheta*std::cos(phi),sintheta*std::sin(phi),costheta);

  //create daughter G4DynamicParticle 
  G4double Ekin = std::sqrt(daughtermomentum*daughtermomentum + daughtermass[0]*daughtermass[0]) - daughtermass[0];
  G4DynamicParticle * daughterparticle = new G4DynamicParticle( G4MT_daughters[0], direction, Ekin, daughtermass[0]);
  products->PushProducts(daughterparticle);
  Ekin = std::sqrt(daughtermomentum*daughtermomentum + daughtermass[1]*daughtermass[1]) - daughtermass[1];
  daughterparticle = new G4DynamicParticle( G4MT_daughters[1], -1.0*direction, Ekin, daughtermass[1]);
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
  G4double parentmass = current_parent_mass.Get();
  //daughters'mass
  G4double daughtermass[3], daughterwidth[3]; 
  G4double sumofdaughtermass = 0.0;
  G4double sumofdaughterwidthsq = 0.0;
  G4bool withWidth = false; 
  for (G4int index=0; index<3; index++){
    daughtermass[index] = G4MT_daughters_mass[index];
    sumofdaughtermass += daughtermass[index];
    daughterwidth[index] = G4MT_daughters_width[index];
    sumofdaughterwidthsq += daughterwidth[index]*daughterwidth[index];
    withWidth = withWidth ||(daughterwidth[index]>1.0e-3*daughtermass[index]);
  }
  
   //create parent G4DynamicParticle at rest
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( G4MT_parent, dummy, 0.0, parentmass);


  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  if (!useGivenDaughterMass) {
    if (withWidth){
      G4double maxDev = (parentmass - sumofdaughtermass )/std::sqrt(sumofdaughterwidthsq) ;
      if (maxDev <= -1.0*rangeMass ){
#ifdef G4VERBOSE
	if (GetVerboseLevel()>0) {
	  G4cout << "G4PhaseSpaceDecayChannel::ThreeBodyDecayIt " 
		 << "sum of daughter mass is larger than parent mass" << G4endl;
	  G4cout << "parent :" << G4MT_parent->GetParticleName() << "  " << current_parent_mass.Get()/GeV << G4endl;
	  G4cout << "daughter 1 :" << G4MT_daughters[0]->GetParticleName() << "  " << daughtermass[0]/GeV << G4endl;
	  G4cout << "daughter 2:" << G4MT_daughters[1]->GetParticleName() << "  " << daughtermass[1]/GeV << G4endl;
	  G4cout << "daughter 3:" << G4MT_daughters[2]->GetParticleName() << "  " << daughtermass[2]/GeV << G4endl;
	}
#endif
	G4Exception("G4PhaseSpaceDecayChannel::ThreeBodyDecayIt",
		    "PART112", JustWarning,
		    "Can not create decay products: sum of daughter mass is larger than parent mass");
	return products;
      }
      G4double dm1=daughtermass[0];
      if (daughterwidth[0] > 0.) dm1= DynamicalMass(daughtermass[0],daughterwidth[0], maxDev);
      G4double dm2= daughtermass[1];
      if (daughterwidth[1] > 0.) dm2= DynamicalMass(daughtermass[1],daughterwidth[1], maxDev);
      G4double dm3= daughtermass[2];
      if (daughterwidth[2] > 0.) dm3= DynamicalMass(daughtermass[2],daughterwidth[2], maxDev);
      while (dm1+dm2+dm3>parentmass){ // Loop checking, 09.08.2015, K.Kurashige
	dm1= DynamicalMass(daughtermass[0],daughterwidth[0], maxDev);
	dm2= DynamicalMass(daughtermass[1],daughterwidth[1], maxDev);
	dm3= DynamicalMass(daughtermass[2],daughterwidth[2], maxDev);
      }
      daughtermass[0] = dm1;
      daughtermass[1] = dm2;
      daughtermass[2] = dm3;
      sumofdaughtermass = dm1+dm2+dm3;
    }
  } else {
    // use given daughter mass;
    daughtermass[0] = givenDaughterMasses[0];
    daughtermass[1] = givenDaughterMasses[1];
    daughtermass[2] = givenDaughterMasses[2];
    sumofdaughtermass = daughtermass[0] + daughtermass[1] + daughtermass[2];
  } 
 
  if (sumofdaughtermass >parentmass) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4PhaseSpaceDecayChannel::ThreeBodyDecayIt " 
	     << "sum of daughter mass is larger than parent mass" << G4endl;
      G4cout << "parent :" << G4MT_parent->GetParticleName() << "  " << current_parent_mass.Get()/GeV << G4endl;
      G4cout << "daughter 1 :" << G4MT_daughters[0]->GetParticleName() << "  " << daughtermass[0]/GeV << G4endl;
      G4cout << "daughter 2:" << G4MT_daughters[1]->GetParticleName() << "  " << daughtermass[1]/GeV << G4endl;
      G4cout << "daughter 3:" << G4MT_daughters[2]->GetParticleName() << "  " << daughtermass[2]/GeV << G4endl;
    }
    if (useGivenDaughterMass) {
      G4cout << "Daughter Mass is given " << G4endl;
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
    if (momentummax <=  momentumsum - momentummax ) break;
  }

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
  G4double Ekin = std::sqrt(daughtermomentum[0]*daughtermomentum[0] + daughtermass[0]*daughtermass[0]) - daughtermass[0];
  G4DynamicParticle * daughterparticle = new G4DynamicParticle( G4MT_daughters[0], 
								direction0, 
								Ekin, daughtermass[0]);
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
  G4ThreeVector pmom = daughtermomentum[2]*direction2/direction2.mag();
  Ekin = std::sqrt(pmom.mag2() + daughtermass[2]*daughtermass[2]) - daughtermass[2];
  daughterparticle = new G4DynamicParticle( G4MT_daughters[2], 
					    pmom/pmom.mag(),
					    Ekin, daughtermass[2]);
  products->PushProducts(daughterparticle);

  pmom = (direction0*daughtermomentum[0] + direction2*(daughtermomentum[2]/direction2.mag()))*(-1.0);
  Ekin = std::sqrt(pmom.mag2() + daughtermass[1]*daughtermass[1]) - daughtermass[1];
  daughterparticle = new G4DynamicParticle( 
					   G4MT_daughters[1],
					   pmom/pmom.mag(),
					   Ekin, daughtermass[1]);
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

#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) G4cout << "G4PhaseSpaceDecayChannel::ManyBodyDecayIt()"<<G4endl;
#endif


  // parent mass
  G4double parentmass = current_parent_mass.Get();

  //parent particle
  G4ThreeVector dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( G4MT_parent, dummy, 0.0, parentmass);

  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  //daughters'mass
  G4double *daughtermass = new G4double[numberOfDaughters]; 
    
  G4double sumofdaughtermass = 0.0;
  for (index=0; index<numberOfDaughters; index++){
    if (!useGivenDaughterMass) {
      daughtermass[index] = G4MT_daughters_mass[index];
    } else {
      daughtermass[index] = givenDaughterMasses[index];
    }
    sumofdaughtermass += daughtermass[index];
  }
  if (sumofdaughtermass >parentmass) {
#ifdef G4VERBOSE
    if (GetVerboseLevel()>0) {
      G4cout << "G4PhaseSpaceDecayChannel::ManyBodyDecayIt " 
	     << "sum of daughter mass is larger than parent mass" << G4endl;
      G4cout << "parent :" << G4MT_parent->GetParticleName() << "  " << current_parent_mass.Get()/GeV << G4endl;
      G4cout << "daughter 1 :" << G4MT_daughters[0]->GetParticleName() << "  " << daughtermass[0]/GeV << G4endl;
      G4cout << "daughter 2:" << G4MT_daughters[1]->GetParticleName() << "  " << daughtermass[1]/GeV << G4endl;
      G4cout << "daughter 3:" << G4MT_daughters[2]->GetParticleName() << "  " << daughtermass[2]/GeV << G4endl;
      G4cout << "daughter 4:" << G4MT_daughters[3]->GetParticleName() << "  " << daughtermass[3]/GeV << G4endl;
    }
#endif
    G4Exception("G4PhaseSpaceDecayChannel::ManyBodyDecayIt",
		"PART112", JustWarning,
		"Can not create decay products: sum of daughter mass is larger than parent mass");
    delete [] daughtermass;
    return products;
  }  

  //Calculate daughter momentum
  G4double *daughtermomentum = new G4double[numberOfDaughters];
  G4ThreeVector direction;  
  G4DynamicParticle **daughterparticle;
  G4double *sm = new G4double[numberOfDaughters];
  G4double tmas;
  G4double weight = 1.0;
  G4int    numberOfTry = 0;
  const size_t MAX_LOOP=10000;

  for (size_t loop_counter=0; loop_counter <MAX_LOOP; ++loop_counter){
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
    G4bool smOK=true;
    for(index =0; index< numberOfDaughters-1 && smOK; index--) {
      smOK = (sm[index]-daughtermass[index]- sm[index+1] >=0.); 
    }
    if (!smOK) continue;

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
	  if (useGivenDaughterMass) {
	    G4cout << "Daughter Mass is given " << G4endl;
	  }
        }
#endif
	delete [] sm;
	delete [] daughtermass;
	delete [] daughtermomentum;
        delete products;
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
      delete products;
      return 0;  // Error detection
    }
    if ( weight < G4UniformRand()) break; 
  }

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
  daughterparticle[index] = new G4DynamicParticle( G4MT_daughters[index], direction*daughtermomentum[index] );
  daughterparticle[index+1] = new G4DynamicParticle( G4MT_daughters[index+1], direction*(-1.0*daughtermomentum[index]) );

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
    daughterparticle[index]= new G4DynamicParticle( G4MT_daughters[index], direction*(-1.0*daughtermomentum[index]));
  }

  //add daughters to G4Decayproducts
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

G4bool  G4PhaseSpaceDecayChannel::SetDaughterMasses( G4double masses[])
{
  for (G4int idx=0; idx<numberOfDaughters; idx++){
    givenDaughterMasses[idx] = masses[idx];
  }
  useGivenDaughterMass = true;
  return useGivenDaughterMass; 
}

G4bool  G4PhaseSpaceDecayChannel::SampleDaughterMasses()
{
  useGivenDaughterMass = false;
  return useGivenDaughterMass; 
}

G4bool  G4PhaseSpaceDecayChannel::IsOKWithParentMass(G4double parentMass) {
  if (!useGivenDaughterMass) return G4VDecayChannel::IsOKWithParentMass(parentMass);
 
  CheckAndFillParent();
  CheckAndFillDaughters();

  G4double sumOfDaughterMassMin=0.0;
  for (G4int index=0; index < numberOfDaughters;  index++) { 
    sumOfDaughterMassMin += givenDaughterMasses[index];
  }
  return (parentMass >= sumOfDaughterMassMin); 
}

G4double G4PhaseSpaceDecayChannel::Pmx(G4double e, G4double p1, G4double p2)
{
   // calcurate momentum of daughter particles in two-body decay
   G4double ppp = (e+p1+p2)*(e+p1-p2)*(e-p1+p2)*(e-p1-p2)/(4.0*e*e);
   if (ppp>0) return std::sqrt(ppp);
   else       return -1.;
}



