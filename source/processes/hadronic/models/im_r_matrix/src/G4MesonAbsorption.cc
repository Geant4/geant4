#include "G4MesonAbsorption.hh"
#include "G4LorentzRotation.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"

// first prototype

G4CollisionInitialState * G4MesonAbsorption::
GetCollision(G4KineticTrack * projectile, vector<G4KineticTrack *> targets)
{
}


vector<G4KineticTrack *> * G4MesonAbsorption::
Scatter(G4KineticTrack * projectile, vector<G4KineticTrack *> & targets)
{
  // Only 2-body absorption for the time being.
  // If insufficient, use 2-body absorption and clusterization to emulate 3-,4-,etc body absorption.
  G4LorentzVector thePro = projectile->Get4Momentum();
  G4LorentzVector theT1 = targets[0]->Get4Momentum();
  G4LorentzVector theT2 = targets[1]->Get4Momentum();
  G4LorentzVector incoming = thePro + theT1 + theT2;
  
  // boost all to MMS
  G4LorentzRotation toSPS((-1)*(thePro + theT1 + theT2).boostVector());
  theT1 = toSPS * theT1;
  theT2 = toSPS * theT2;
  thePro = toSPS * thePro;
  G4LorentzRotation fromSPS(toSPS.inverse());
 
  // rotate to z
  G4LorentzRotation toZ;
  G4LorentzVector Ptmp=projectile->Get4Momentum();
  toZ.rotateZ(-1*Ptmp.phi());
  toZ.rotateY(-1*Ptmp.theta());
  theT1 = toZ * theT1;
  theT2 = toZ * theT2;
  thePro = toZ * thePro;
  G4LorentzRotation toLab(toZ.inverse());
  
  // calculate the momenta.
  G4double M = (thePro+theT1+theT2).mag();
  G4double m1 = targets[0]->GetDefinition()->GetPDGMass();
  G4double m2 = targets[1]->GetDefinition()->GetPDGMass();
  G4double m = sqrt(M*M-m1*m1-m2*m2);
  G4double p = sqrt((m*m*m*m - 4.*m1*m1 * m2*m2)/(4.*(M*M)));
  G4double costh = 2.*G4UniformRand()-1.;
  G4double phi = 2.*pi*G4UniformRand();
  G4ThreeVector pFinal(p*sin(acos(costh))*cos(phi), p*sin(acos(costh))*sin(phi), p*costh);
  
  // G4cout << "testing p "<<p-pFinal.mag()<<G4endl;
  // construct the final state particles lorentz momentum.
  G4double eFinal1 = sqrt(m1*m1+pFinal.mag2());
  G4LorentzVector final1(pFinal, eFinal1);
  G4double eFinal2 = sqrt(m2*m2+pFinal.mag2());
  G4LorentzVector final2(-1.*pFinal, eFinal2);
  
  // rotate back.
  final1 = toLab * final1;
  final2 = toLab * final2;
  
  // boost back to LAB
  final1 = fromSPS * final1;
  final2 = fromSPS * final2;
  
  // make the final state
  G4KineticTrack * f1 = 
      new G4KineticTrack(targets[0]->GetDefinition(), 0., 
                        targets[0]->GetPosition(), final1);
  G4KineticTrack * f2 = 
      new G4KineticTrack(targets[1]->GetDefinition(), 0., 
                        targets[1]->GetPosition(), final2);
  vector<G4KineticTrack *> * result = new vector<G4KineticTrack *>;
  result->push_back(f1);
  result->push_back(f2);

  return result;
}

