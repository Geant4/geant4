// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DecayProducts.cc,v 1.2 1999-02-06 10:10:13 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      History: first implementation, based on object model of
//      10 July 1996 H.Kurashige
//      21 Oct  1996 H.Kurashige
//      12 Dec 1997 H.Kurashige
// ------------------------------------------------------------

#include "G4ios.hh"
#include "globals.hh"
#include "G4DecayProducts.hh"

#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"

G4Allocator<G4DecayProducts> aDecayProductsAllocator;


G4DecayProducts::G4DecayProducts()
                :numberOfProducts(0),theParentParticle(NULL)
{ 

}

G4DecayProducts::G4DecayProducts(const G4DynamicParticle &aParticle)
                :numberOfProducts(0),theParentParticle(NULL)
{
  theParentParticle = new G4DynamicParticle(aParticle);
}

G4DecayProducts::G4DecayProducts(const G4DecayProducts &right) 
                :numberOfProducts(0)
{
  // copy parent (Deep Copy)
  theParentParticle = new G4DynamicParticle(*right.theParentParticle);

  //copy daughters (Deep Copy)
  for (G4int index=0; index < right.numberOfProducts; index++)
  {
    PushProducts( new G4DynamicParticle(*right.theProductVector[index]) );
  }
}

G4DecayProducts & G4DecayProducts::operator=(const G4DecayProducts &right)
{
  G4int index;

  if (this != &right)
  { 
    // recreate parent
    if (theParentParticle != NULL) delete theParentParticle;
    theParentParticle = new G4DynamicParticle(*right.theParentParticle);

    // delete G4DynamicParticle objects
    for (index=0; index < numberOfProducts; index++)
    {
      delete theProductVector[index];
    }

    //copy daughters (Deep Copy)
    for (index=0; index < right.numberOfProducts; index++) {
      PushProducts( new G4DynamicParticle(*right.theProductVector[index]) );
    } 
    numberOfProducts = right.numberOfProducts;
    
  }
  return *this;
}

G4DecayProducts::~G4DecayProducts()
{
  //delete parent
  if (theParentParticle != NULL) delete theParentParticle;
  
  // delete G4DynamicParticle object
  for (G4int index=0; index < numberOfProducts; index++)
  {
      delete theProductVector[index];
  }
  numberOfProducts = 0;    
}

G4DynamicParticle* G4DecayProducts::PopProducts()
{
   if ( numberOfProducts >0 ) {
     numberOfProducts -= 1;   
     return  theProductVector[numberOfProducts];
   } else {
     return NULL;
   }
}

G4int G4DecayProducts::PushProducts(G4DynamicParticle *aParticle)
{
   if ( numberOfProducts < MaxNumberOfProducts) {
     theProductVector[numberOfProducts]= aParticle;
     numberOfProducts += 1; 
   } else {
#ifdef G4VERBOSE
     G4cerr << "G4DecayProducts::PushProducts ";
     G4cerr << " exceeds MaxNumberOfProducts(=" <<MaxNumberOfProducts << ")";
     G4cerr << endl;
#endif
   }
   return numberOfProducts;
}

G4DynamicParticle* G4DecayProducts::operator[](G4int anIndex) const
{
   if ((numberOfProducts > anIndex) && (anIndex >=0) ) {
     return  theProductVector[anIndex];
   } else {
     return NULL;
   }
}

void  G4DecayProducts::SetParentParticle(const G4DynamicParticle &aParticle)
{
  if (theParentParticle != NULL) delete theParentParticle;
  theParentParticle = new G4DynamicParticle(aParticle);
}


void G4DecayProducts::Boost(G4double totalEnergy, const G4ThreeVector &momentumDirection)
{
  // calcurate new beta
  G4double   mass = theParentParticle->GetMass();
  G4double   totalMomentum  = sqrt( (totalEnergy - mass)*(totalEnergy + mass) );
  G4double   betax = momentumDirection.x()*totalMomentum/totalEnergy;  
  G4double   betay = momentumDirection.y()*totalMomentum/totalEnergy;  
  G4double   betaz = momentumDirection.z()*totalMomentum/totalEnergy;  
  this->Boost(betax, betay, betaz);
}

void G4DecayProducts::Boost(G4double newbetax, G4double newbetay, G4double newbetaz)
{ 
  G4double   mass = theParentParticle->GetMass();
  G4double   energy  = theParentParticle->GetTotalEnergy();
  G4double   momentum  = 0.0;

  G4ThreeVector direction(0.0,0.0,1.0);    
  G4LorentzVector p4;

  if (energy - mass > DBL_MIN) {
    // calcurate  beta of initial state
    momentum  = theParentParticle->GetTotalMomentum();
    direction = theParentParticle->GetMomentumDirection();
    G4double betax = -1.0*direction.x()*momentum/energy;  
    G4double betay = -1.0*direction.y()*momentum/energy;  
    G4double betaz = -1.0*direction.z()*momentum/energy;  
    
    for (G4int index=0; index < numberOfProducts; index++) {
       // make G4LorentzVector for secondaries
       p4 = theProductVector[index]->Get4Momentum();

       // boost secondaries to theParentParticle's rest frame 
       p4.boost(betax, betay, betaz);

       // boost secondaries to  new frame 
       p4.boost(newbetax, newbetay, newbetaz);

       // change energy/momentum
       theProductVector[index]->Set4Momentum(p4);
     }
   } else {
     for (G4int index=0; index < numberOfProducts; index++) {
       // make G4LorentzVector for secondaries
       p4 = theProductVector[index]->Get4Momentum();

       // boost secondaries to  new frame 
       p4.boost(newbetax, newbetay, newbetaz);

       // change energy/momentum
       theProductVector[index]->Set4Momentum(p4);
      }
   }
   // make G4LorentzVector for parent in its rest frame
   mass = theParentParticle->GetMass();
   G4LorentzVector parent4( 0.0, 0.0, 0.0, mass);

   // boost parent to new frame 
   parent4.boost(newbetax, newbetay, newbetaz);

   // change energy/momentum
   theParentParticle->Set4Momentum(parent4);
}

G4bool G4DecayProducts::IsChecked()
{
  G4bool returnValue = true;
  // check parent 
  //   energy/momentum
  G4double   parent_mass = theParentParticle->GetMass();
  G4double   parent_energy  = theParentParticle->GetTotalEnergy();
  G4ThreeVector direction = theParentParticle->GetMomentumDirection();
  G4ThreeVector parent_momentum = direction*(theParentParticle->GetTotalMomentum());
  // check momentum dirction is a unit vector
  if ( (parent_momentum.mag() >0.0) && (abs(direction.mag()-1.0) >1.0e-6 ) ) {
#ifdef G4VERBOSE
    G4cout << " Momentum Direction Vector of Parent is not normalized ";
    G4cout << "  (=" << direction.mag() << ")" << endl;
#endif
    returnValue = false;
    parent_momentum = parent_momentum * (1./direction.mag());
  }

  //daughters
  G4double   mass, energy;
  G4ThreeVector momentum;
  G4double   total_energy = parent_energy;
  G4ThreeVector total_momentum =  parent_momentum;
  for (G4int index=0; index < numberOfProducts; index++) 
  {
    mass = theProductVector[index]->GetMass();
    energy  = theProductVector[index]->GetTotalEnergy();
    direction = theProductVector[index]->GetMomentumDirection();
    momentum = direction*(theProductVector[index]->GetTotalMomentum());
    // check momentum dirction is a unit vector
    if ( (momentum.mag()>0.0) && (abs(direction.mag()-1.0) > 1.0e-6)) {
#ifdef G4VERBOSE
      G4cout << " Momentum Direction Vector of Daughter [" << index;
      G4cout << "]  is not normalized (=" << direction.mag() << ")" << endl;
#endif
      returnValue = false;
      momentum = momentum * (1./direction.mag());
    }
    // whether daughter stops or not
    if (energy - mass < DBL_MIN ) {
#ifdef G4VERBOSE
      G4cout << "Daughter [" << index << "] has no kinetic energy "<< endl;
#endif
      returnValue = false;
    }
    total_energy -= energy; 
    total_momentum -= momentum;
  }
  // check energy/momentum conservation
  if ( (abs(total_energy) >1.0e-6) || (total_momentum.mag() >1.0e-6 ) ){ 
#ifdef G4VERBOSE
    G4cout << " Energy/Momentum is not conserved   ";
    G4cout << total_energy << "  " << total_momentum << endl;
#endif
    returnValue = false;
  }
  return returnValue;
}

void G4DecayProducts::DumpInfo()
{
   G4cout << " ----- List of DecayProducts  -----" << endl;
   G4cout << " ------ Parent Particle ----------" << endl;
   if (theParentParticle != NULL) theParentParticle->DumpInfo();
   G4cout << " ------ Daughter Particles  ------" << endl;  
   for (G4int index=0; index < numberOfProducts; index++) 
   {
      G4cout << " ----------" << index+1 << " -------------" << endl;  
      theProductVector[index]-> DumpInfo();
   }
   G4cout << " ----- End List of DecayProducts  -----" << endl;
   G4cout << endl;
} 
