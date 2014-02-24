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
#include "G4Absorber.hh"
#include "G4KineticTrack.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4LorentzRotation.hh"

G4Absorber::G4Absorber(G4double cutOnP)
{
  theCutOnP = cutOnP;
  theAbsorbers = new G4KineticTrackVector;
  theProducts = new G4KineticTrackVector;
}


G4Absorber::~G4Absorber()
{
  delete theAbsorbers;
  delete theProducts;
}


bool G4Absorber::WillBeAbsorbed(const G4KineticTrack & kt)
{
 // FixMe: actually only for pions
//  if(kt.Get4Momentum().vect().mag() < theCutOnP) 
// Cut on kinetic Energy...
  if (kt.Get4Momentum().e() - kt.GetActualMass() < theCutOnP) 
  {
      if(kt.GetDefinition() == G4PionPlus::PionPlus()  ||
	 kt.GetDefinition() == G4PionZero::PionZero()  ||
	 kt.GetDefinition() == G4PionMinus::PionMinus())
      {
	 return true;
      }   
  }
  return false;
}



G4bool G4Absorber::Absorb(G4KineticTrack & kt, G4KineticTrackVector & tgt)
{
  if(!FindAbsorbers(kt, tgt))
    return false;
  return FindProducts(kt);
}


G4bool G4Absorber::FindAbsorbers(G4KineticTrack & kt,
				 G4KineticTrackVector & tgt)
{
//  Find a closest ( in space) pair of Nucleons capable to absorb pi+/pi-
//    pi+ can be absorbed on np or nn resulting in pp or np
//    pi- can be absorbed on np or pp resulting in nn or np

// @GF: FindAbsorbers is unused, logic is seriously wrong

  G4KineticTrack * kt1 = NULL;
  G4KineticTrack * kt2 = NULL;
  G4double dist1 = DBL_MAX;		// dist to closest nucleon
  G4double dist2 = DBL_MAX;		// dist to next close 
  G4double charge1 = 0;
//  G4double charge2 = 0;		// charge2 is only assigned to, never used
  G4double charge0 = kt.GetDefinition()->GetPDGCharge();
  G4ThreeVector pos = kt.GetPosition();

  std::vector<G4KineticTrack *>::iterator iter;
  for(iter = tgt.begin(); iter != tgt.end(); ++iter)
  {
    G4KineticTrack * curr = *iter;
    G4double dist = (pos-curr->GetPosition()).mag();
    if(dist >= dist2)
      continue;
    if(dist < dist1)
    {
      if(dist1 == DBL_MAX) // accept 1st as a candidate, 
      {
	kt1 = curr;
	charge1 = kt1->GetDefinition()->GetPDGCharge();
	dist1 = dist;
	continue;
      }	
      if(dist2 == DBL_MAX) // accept the candidate and shift kt1 to kt2
      {				// @GF: should'nt we check if compatible?
	kt2 = kt1;
//	charge2 = charge1;
	dist2 = dist1;
	kt1 = curr;
	charge1 = kt1->GetDefinition()->GetPDGCharge();
	dist1 = dist;
	continue;
      }
// test the compatibility with charge conservation for new config
      G4double charge = curr->GetDefinition()->GetPDGCharge();
      if((charge0+charge1+charge < 0.) ||	//test config (curr,kt1) 
	 (charge0+charge1+charge) > 2*eplus)
      {  // incompatible: change kt1 with curr.
	kt1 = curr;
	charge1 = charge;
	dist1 = dist;
      }
      else
      { // compatible: change kt1 with curr and kt2 with kt1
	kt2 = kt1;
//	charge2 = charge1;
	dist2 = dist1;
	kt1 = curr;
	charge1 = charge;
	dist1 = dist;
      }
      continue;
    }
// here if dist1 < dist < dist2
    if(dist2 == DBL_MAX) // accept the candidate
    {
      kt2 = curr;
//      charge2 = kt2->GetDefinition()->GetPDGCharge();
      dist2 = dist;
      continue;
    }	
// test the compatibility with charge conservation
    G4double charge = curr->GetDefinition()->GetPDGCharge();
    if((charge0+charge1+charge < 0.) ||
       (charge0+charge1+charge) > 2*eplus)
      continue;   // incomatible: do nothing
// compatible: change kt2 with curr
    kt2 = curr;
//    charge2 = charge;
    dist2 = dist;
  }

  theAbsorbers->clear(); // do not delete tracks in theAbsorbers vector!
  if((kt1 == NULL) || (kt2 == NULL))
    return false;

  theAbsorbers->push_back(kt1);
  theAbsorbers->push_back(kt2);
  return true;
}



G4bool G4Absorber::FindProducts(G4KineticTrack & kt)
{
// Choose the products type
  const G4ParticleDefinition * prod1;
  const G4ParticleDefinition * prod2;
  G4KineticTrack * abs1 = (*theAbsorbers)[0];
  G4KineticTrack * abs2 = (*theAbsorbers)[1];

  G4double charge = kt.GetDefinition()->GetPDGCharge();
  if(charge == eplus)
  { // a neutron become proton
    prod1 = G4Proton::Proton();
    if(abs1->GetDefinition() == G4Neutron::Neutron())
      prod2 = abs2->GetDefinition();
    else
      prod2 = G4Proton::Proton();
  }
  else if(charge == -eplus)
  { // a proton become neutron
    prod1 = G4Neutron::Neutron();
    if(abs1->GetDefinition() == G4Proton::Proton())
      prod2 = abs2->GetDefinition();
    else
      prod2 = G4Neutron::Neutron();
  }
  else  // charge = 0: leave particle types unchenged
  {
    prod1 = abs1->GetDefinition();
    prod2 = abs2->GetDefinition();
  }

// Translate to the CMS frame
  G4LorentzVector momLab = kt.Get4Momentum()+abs1->Get4Momentum()+
    abs2->Get4Momentum();
  G4LorentzRotation toCMSFrame((-1)*momLab.boostVector());
  G4LorentzRotation toLabFrame(momLab.boostVector());
  G4LorentzVector momCMS = toCMSFrame*momLab;

// Evaluate the final momentum of products
  G4double ms1 = prod1->GetPDGMass();
  G4double ms2 = prod2->GetPDGMass();
  G4double e0 = momCMS.e();
  G4double squareP = (e0*e0*e0*e0-2*e0*e0*(ms1*ms1+ms2*ms2)+
    (ms2*ms2-ms1*ms1)*(ms2*ms2-ms1*ms1))/(4*e0*e0);
//  if(squareP < 0)  // should never happen
//    squareP = 0;
  G4ThreeVector mom1CMS = GetRandomDirection();
  mom1CMS = std::sqrt(squareP)*mom1CMS;
  G4LorentzVector final4Mom1CMS(mom1CMS, std::sqrt(squareP+ms1*ms1));
  G4LorentzVector final4Mom2CMS((-1)*mom1CMS, std::sqrt(squareP+ms2*ms2));

// Go back to the lab frame
  G4LorentzVector mom1 = toLabFrame*final4Mom1CMS;
  G4LorentzVector mom2 = toLabFrame*final4Mom2CMS;

// ------ debug
/*
  G4LorentzVector temp = mom1+mom2;

  cout << (1/MeV)*momLab.x() << " " << (1/MeV)*momLab.y() << " "
       << (1/MeV)*momLab.z() << " " << (1/MeV)*momLab.t() << " "
       << (1/MeV)*momLab.vect().mag() << " " << (1/MeV)*momLab.mag() << " "
       << (1/MeV)*temp.x() << " " << (1/MeV)*temp.y() << " "
       << (1/MeV)*temp.z() << " " << (1/MeV)*temp.t() << " "
       << (1/MeV)*temp.vect().mag() << " " << (1/MeV)*temp.mag() << " "
       << (1/MeV)*std::sqrt(squareP) << endl;

*/
// ------ end debug

// Build two new kinetic tracks and add to products
  G4KineticTrack * kt1 = new G4KineticTrack(prod1, 0., abs1->GetPosition(),
					    mom1);
  G4KineticTrack * kt2 = new G4KineticTrack(prod2, 0., abs2->GetPosition(),
					    mom2);
// ------ debug
/*
  G4LorentzVector initialMom1 = abs1->Get4Momentum();
  G4LorentzVector initialMom2 = abs2->Get4Momentum();
  G4LorentzVector pion4MomCMS = toCMSFrame*kt.Get4Momentum();
  cout << (1/MeV)*initialMom1.x() << " " << (1/MeV)*initialMom1.y() << " "
       << (1/MeV)*initialMom1.z() << " " << (1/MeV)*initialMom1.e() << " "
       << (1/MeV)*initialMom1.vect().mag() << " "
       << (1/MeV)*initialMom2.x() << " " << (1/MeV)*initialMom2.y() << " "
       << (1/MeV)*initialMom2.z() << " " << (1/MeV)*initialMom2.e() << " "
       << (1/MeV)*initialMom2.vect().mag() << " "
       << (1/MeV)*mom1.x() << " " << (1/MeV)*mom1.y() << " "
       << (1/MeV)*mom1.z() << " " << (1/MeV)*mom1.e() << " "
       << (1/MeV)*mom1.vect().mag() << " "
       << (1/MeV)*mom2.x() << " " << (1/MeV)*mom2.y() << " "
       << (1/MeV)*mom2.z() << " " << (1/MeV)*mom2.e() << " "
       << (1/MeV)*mom2.vect().mag() << " "
       << (1/MeV)*pion4MomCMS.x() << " " << (1/MeV)*pion4MomCMS.y() << " "
       << (1/MeV)*pion4MomCMS.z() << " " << (1/MeV)*pion4MomCMS.e() << " "
       << (1/MeV)*pion4MomCMS.vect().mag() << " "
       << (1/MeV)*final4Mom1CMS.x() << " " << (1/MeV)*final4Mom1CMS.y() << " "
       << (1/MeV)*final4Mom1CMS.z() << " " << (1/MeV)*final4Mom1CMS.e() << " "
       << (1/MeV)*final4Mom1CMS.vect().mag() << " "
       << (1/MeV)*final4Mom2CMS.x() << " " << (1/MeV)*final4Mom2CMS.y() << " "
       << (1/MeV)*final4Mom2CMS.z() << " " << (1/MeV)*final4Mom2CMS.e() << " "
       << (1/MeV)*final4Mom2CMS.vect().mag() << endl;
*/
// ------ end debug

  theProducts->clear();
  theProducts->push_back(kt1);
  theProducts->push_back(kt2);
  return true;
}



G4ThreeVector G4Absorber::GetRandomDirection()
{
  G4double theta = 2.0*G4UniformRand()-1.0;
  theta = std::acos(theta);
  G4double phi = G4UniformRand()*2*pi;
  G4ThreeVector direction(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta));
  return direction;
}







