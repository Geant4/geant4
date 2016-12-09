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
// $Id: G4FermiFragmentsPool.cc,v 1.5 2006-06-29 20:13:13 gunter Exp $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara
//
// Modifications:
// J.M.Quesada,  July 2009, bug fixed in excitation energies: 
// ALL of them are in MeV instead of keV (as they were expressed previously)
// source:  http://www.nndc.bnl.gov/chart
// Unknown excitation energies in He5  and Li5 have been suppressed
// Long lived levels (half-lives of the order ps-fs have been included)   
//
// J. M. Quesada,  April 2010: excitation energies according to tabulated values 
// in PhotonEvaporatoion2.0. Fake photons eliminated. 
//
// 01.04.2011 General cleanup by V.Ivanchenko - more clean usage of static
//
// 04.05.2011 J. M. Quesada: added detailed printout for testing

#include "G4FermiFragmentsPool.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4StableFermiFragment.hh"
#include "G4B9FermiFragment.hh"
#include "G4Be8FermiFragment.hh"
#include "G4He5FermiFragment.hh"
#include "G4Li5FermiFragment.hh"

G4FermiFragmentsPool* G4FermiFragmentsPool::theInstance = nullptr;
 
G4FermiFragmentsPool* G4FermiFragmentsPool::Instance()
{
  if(0 == theInstance) {
    static G4FermiFragmentsPool instance;
    theInstance = &instance;
  }
  return theInstance;
}

G4FermiFragmentsPool::G4FermiFragmentsPool()
{
  maxZ = 9;
  maxA = 17;
  verbose = 0;
  Initialise();
}

G4FermiFragmentsPool::~G4FermiFragmentsPool()
{
  size_t nn;
  for(size_t i=0; i<17; ++i) {
    nn = list1[i].size();
    if(0 < nn) { for(size_t j=0; j<nn; ++j) { delete (list1[i])[j]; }}
    nn = list2[i].size();
    if(0 < nn) { for(size_t j=0; j<nn; ++j) { delete (list2[i])[j]; }}
    nn = list3[i].size();
    if(0 < nn) { for(size_t j=0; j<nn; ++j) { delete (list3[i])[j]; }}
    nn = list4[i].size();
    if(0 < nn) { for(size_t j=0; j<nn; ++j) { delete (list4[i])[j]; }}
  }
  nn = fragment_pool.size();
  if(0 < nn) { for(size_t j=0; j<nn; ++j) { delete fragment_pool[j]; }}
}

G4bool G4FermiFragmentsPool::IsAvailable(G4int Z, G4int A) const
{
  G4bool res = true;
  if     (2 == Z && 5 == A) { res = false; }
  else if(3 == Z && 5 == A) { res = false; }
  else if(4 == Z && 8 == A) { res = false; }
  else if(5 == Z && 9 == A) { res = false; }
  return res;
}

G4int G4FermiFragmentsPool::GetMaxZ() const
{
  return maxZ;
}

G4int G4FermiFragmentsPool::GetMaxA() const
{
  return maxA;
}

const G4FermiPhaseSpaceDecay* 
G4FermiFragmentsPool::GetFermiPhaseSpaceDecay() const
{
  return &thePhaseSpace;
}

void G4FermiFragmentsPool::Initialise()
{
  // JMQ 30/06/09 unknown levels have been supressed
  // JMQ 01/07/09 corrected excitation energies for 64-66, according to 
  // http://www.nndc.bnl.gov/chart
  // JMQ 19/04/10 new level, fragment numbering shifted accordingly from here onwards
  //                                                 A  Z  Pol  ExcitE
  fragment_pool.push_back(new G4StableFermiFragment(  1, 0,  2,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment(  1, 1,  2,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment(  2, 1,  3,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment(  3, 1,  2,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment(  3, 2,  2,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment(  4, 2,  1,  0.00*MeV )); 
  fragment_pool.push_back(new G4He5FermiFragment   (  5, 2,  4,  0.00*MeV )); 
  fragment_pool.push_back(new G4Li5FermiFragment   (  5, 3,  4,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment(  6, 2,  1,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment(  6, 3,  3,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment(  6, 3,  1,  3.562880*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment(  7, 3,  4,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment(  7, 3,  2,  0.4776120*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment(  7, 4,  4,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment(  7, 4,  2,  0.4290800*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment(  8, 3,  5,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment(  8, 3,  3,  0.9808000*MeV )); 
  fragment_pool.push_back(new G4Be8FermiFragment   (  8, 4,  1,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment(  9, 4,  4,  0.00*MeV )); 
  fragment_pool.push_back(new G4B9FermiFragment    (  9, 5,  4,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 10, 4,  1,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 10, 4,  5,  3.368030*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 10, 4,  8,  5.958390*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 10, 4,  1,  6.179300*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 10, 4,  5,  6.263300*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 10, 5,  7,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 10, 5,  3,  0.7183800*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 10, 5,  1,  1.740050*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 10, 5,  3,  2.154270*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 10, 5,  5,  3.587130*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 10, 6,  3,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 10, 6,  5,  3.353700*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 5,  4,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 5,  2,  2.124693*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 5,  6,  4.444980*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 5,  4,  5.020300*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 5,  8,  6.741850*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 5,  2,  6.791800*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 5,  6,  7.285510*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 5,  4,  7.977840*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 5,  6,  8.560100*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 6,  4,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 6,  2,  2.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 6,  6,  4.318800*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 6,  4,  4.804200*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 6,  2,  6.339200*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 6,  8,  6.478200*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 6,  6,  6.904800*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 6,  4,  7.499700*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 6,  4,  8.104500*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 11, 6,  6,  8.420000*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 12, 5,  3,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 12, 5,  5,  0.9531400*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 12, 5,  5,  1.673650*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 12, 5,  3,  2.620800*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 12, 6,  1,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 12, 6,  5,  4.438910*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 13, 6,  2,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 13, 6,  2,  3.089443*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 13, 6,  4,  3.684507*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 13, 6,  6,  3.853807*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 13, 7,  2,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 6,  1,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 6,  3,  6.093800*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 6,  1,  6.589400*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 6,  7,  6.728200*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 6,  1,  6.902600*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 6,  5,  7.012000*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 6,  5,  7.341000*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 7,  3,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 7,  1,  2.312798*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 7,  3,  3.948100*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 7,  1,  4.915100*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 7,  5,  5.105890*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 7,  3,  5.691440*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 7,  7,  5.834250*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 7,  3,  6.203500*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 7,  7,  6.446170*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 14, 7,  5,  7.029120*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  2,  0.00*MeV )); 
  // JMQ 010709 two very close levels instead of only one, with their own spins
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  6,  5.270155*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  2,  5.298822*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  4,  6.323780*MeV )); 
  //JMQ 010709 new level and corrected energy and spins
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  6,  7.155050*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  4,  7.300830*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  8,  7.567100*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  2,  8.312620*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  4,  8.571400*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  2,  9.049710*MeV )); 
  //JMQ 010709 new levels for N15
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  4,  9.151900*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  6,  9.154900*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  2,  9.222100*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  6,  9.760000*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  8,  9.829000*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  4,  9.925000*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 7,  4, 10.06600*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 8,  2,  0.00*MeV )); 
  //JMQ 010709 new level and spins
  fragment_pool.push_back(new G4StableFermiFragment( 15, 8,  2,  5.183000*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 8,  6,  5.240900*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 8,  4,  6.176300*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 8,  4,  6.793100*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 8,  6,  6.859400*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 15, 8,  8,  7.275900*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 16, 7,  5,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 16, 7,  1,  0.1204200*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 16, 7,  7,  0.2982200*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 16, 7,  3,  0.3972700*MeV )); 
  //JMQ 010709   some energies and spins have been changed 
  fragment_pool.push_back(new G4StableFermiFragment( 16, 8,  1,  0.00*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 16, 8,  1,  6.049400*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 16, 8,  7,  6.129890*MeV )); 
  fragment_pool.push_back(new G4StableFermiFragment( 16, 8,  5,  6.917100*MeV )); 
  //JMQ 180510 fixed fragment 111
  fragment_pool.push_back(new G4StableFermiFragment( 16, 8,  3,  7.116850*MeV )); 

  G4int nfrag = fragment_pool.size();
  std::vector<const G4VFermiFragment*> newvec;
  newvec.reserve(4);

  // list of fragments ordered by A
  newvec.resize(1);
  for(G4int i=0; i<nfrag; ++i) {
    newvec[0] = fragment_pool[i];
    const G4FermiConfiguration* conf = new G4FermiConfiguration(newvec);
    G4int A = fragment_pool[i]->GetA();
    list1[A].push_back(conf);
  }
  if(verbose > 0) { 
    G4cout << "### G4FermiFragmentPool: " << nfrag 
	   << " fragments" << G4endl;
    for(G4int A=1; A<maxA; ++A) {
      G4cout << "  A= " << A << " : Z= ";
      for(size_t j=0; j<list1[A].size(); ++j) { 
	G4cout << (list1[A])[j]->GetZ() << "  "; 
      }
      G4cout << G4endl;
    }
  }

  // list of fragment pairs ordered by A
  G4int counter = 0;
  G4int tot = 0;
  newvec.resize(2);
  for(G4int i=0; i<nfrag; ++i) {
    G4int Z1 = fragment_pool[i]->GetZ();
    G4int A1 = fragment_pool[i]->GetA();
    for(G4int j=0; j<nfrag; ++j) {
      G4int Z2 = fragment_pool[j]->GetZ();
      G4int A2 = fragment_pool[j]->GetA();
      G4int Z = Z1 + Z2;
      G4int A = A1 + A2;
      if(Z < maxZ && A < maxA && IsAvailable(Z, A)) {
	newvec[0] = fragment_pool[i];
	newvec[1] = fragment_pool[j];
	if(!IsExist(Z, A, newvec)) { 
	  const G4FermiConfiguration* conf = new G4FermiConfiguration(newvec);
	  list2[A].push_back(conf); 
	  ++counter;
	}
      }
    }
  }
  if(verbose > 0) { 
    G4cout << G4endl;
    G4cout << "### Pairs of fragments: " << counter << G4endl;
    for(G4int A=2; A<maxA; ++A) {
      G4cout << "  A= " << A<<G4endl; 
      for(size_t j=0; j<list2[A].size(); ++j) {
        const std::vector<const G4VFermiFragment*>* vec 
	  = (list2[A])[j]->GetFragmentList(); 
	G4int a1=(*vec)[0]->GetA();
	G4int z1=(*vec)[0]->GetZ();
	G4int a2=(*vec)[1]->GetA();
	G4int z2=(*vec)[1]->GetZ();
 	G4cout << "("<<a1<<","<<z1<<")("<<a2<<","<<z2<<") % "; 
      }
      G4cout<<G4endl;
      G4cout<<"----------------------------------------------------------------"
	    << G4endl;
    }
  }

  // list of fragment triples ordered by A
  tot += counter;
  counter = 0;
  newvec.resize(3);
  for(G4int A1=2; A1<maxA; ++A1) {
    size_t nz = list2[A1].size();
    for(size_t idx=0; idx<nz; ++idx) {
      const G4FermiConfiguration* conf2 = (list2[A1])[idx];
      G4int Z1 = conf2->GetZ();
      const std::vector<const G4VFermiFragment*>* vec2 = 
	conf2->GetFragmentList(); 
      for(G4int j=0; j<nfrag; ++j) {
	G4int Z2 = fragment_pool[j]->GetZ();
	G4int A2 = fragment_pool[j]->GetA();
	G4int Z = Z1 + Z2;
	G4int A = A1 + A2;
	if(Z < maxZ && A < maxA && IsAvailable(Z, A)) {
	  newvec[0] = (*vec2)[0];
	  newvec[1] = (*vec2)[1];
	  newvec[2] = fragment_pool[j];
	  if(!IsExist(Z, A, newvec)) { 
	    const G4FermiConfiguration* conf3 = new G4FermiConfiguration(newvec);
	    list3[A].push_back(conf3);
	    ++counter;
	  }
	}
      }
    }
  }
  if(verbose > 0) { 
    G4cout << G4endl;
    G4cout << "### Triples of fragments: " << counter << G4endl;
    for(G4int A=3; A<maxA; ++A) {
      G4cout << "  A= " << A<<G4endl;
      for(size_t j=0; j<list3[A].size(); ++j) { 
	const std::vector<const G4VFermiFragment*>* vec 
	  = (list3[A])[j]->GetFragmentList(); 
	G4int a1=(*vec)[0]->GetA();
	G4int z1=(*vec)[0]->GetZ();
	G4int a2=(*vec)[1]->GetA();
	G4int z2=(*vec)[1]->GetZ();
	G4int a3=(*vec)[2]->GetA();
	G4int z3=(*vec)[2]->GetZ();
 	G4cout << "("<<a1<<","<<z1<<")("<<a2<<","<<z2<<")("<<a3<<","<<z3<<") % ";
      }
      G4cout<<G4endl;
      G4cout<<"----------------------------------------------------------------"
	    << G4endl;
    }
  }

  // list of fragment quartets (3 + 1) ordered by A
  tot += counter;
  counter = 0;
  newvec.resize(4);
  for(G4int A1=3; A1<maxA; ++A1) {
    size_t nz = list3[A1].size();
    for(size_t idx=0; idx<nz; ++idx) {
      const G4FermiConfiguration* conf3 = (list3[A1])[idx];
      G4int Z1 = conf3->GetZ();
      const std::vector<const G4VFermiFragment*>* vec3 = 
	conf3->GetFragmentList(); 
      for(G4int j=0; j<nfrag; ++j) {
	G4int Z2 = fragment_pool[j]->GetZ();
	G4int A2 = fragment_pool[j]->GetA();
	G4int Z = Z1 + Z2;
	G4int A = A1 + A2;
	if(Z < maxZ && A < maxA && IsAvailable(Z, A)) {
	  newvec[0] = (*vec3)[0];
	  newvec[1] = (*vec3)[1];
	  newvec[2] = (*vec3)[2];
	  newvec[3] = fragment_pool[j];
	  if(!IsExist(Z, A, newvec)) { 
	    const G4FermiConfiguration* conf4 = new G4FermiConfiguration(newvec);
	    list4[A].push_back(conf4);
	    ++counter;
	  }
	}
      }
    }
  }
  // list of fragment quartets (2 + 2) ordered by A
  for(G4int A1=2; A1<maxA; ++A1) {
    size_t nz1 = list2[A1].size();
    for(size_t id1=0; id1<nz1; ++id1) {
      const G4FermiConfiguration* conf1 = (list2[A1])[id1];
      G4int Z1 = conf1->GetZ();
      const std::vector<const G4VFermiFragment*>* vec1 = 
	conf1->GetFragmentList(); 
      for(G4int A2=2; A2<maxA; ++A2) {
	size_t nz2 = list2[A2].size();
	for(size_t id2=0; id2<nz2; ++id2) {
	  const G4FermiConfiguration* conf2 = (list2[A2])[id2];
	  G4int Z2 = conf2->GetZ();
	  const std::vector<const G4VFermiFragment*>* vec2 = 
	    conf2->GetFragmentList(); 
	  G4int Z = Z1 + Z2;
	  G4int A = A1 + A2;
	  if(Z < maxZ && A < maxA && IsAvailable(Z, A)) {
	    newvec[0] = (*vec1)[0];
	    newvec[1] = (*vec1)[1];
	    newvec[2] = (*vec2)[0];
	    newvec[3] = (*vec2)[1];
	    if(!IsExist(Z, A, newvec)) { 
	      const G4FermiConfiguration* conf4 = 
		new G4FermiConfiguration(newvec);
	      list4[A].push_back(conf4);
	      ++counter;
	    }
	  }
	}
      }
    }
  }
  if(verbose > 0) { 
    tot += counter;
    G4cout << G4endl;
    G4cout << "### Quartets of fragments: " << counter << G4endl;
    for(G4int A=4; A<maxA; ++A) {
      G4cout << "  A= " << A<<G4endl;
      for(size_t j=0; j<list4[A].size(); ++j) { 
	const std::vector<const G4VFermiFragment*>* vec 
	  = (list4[A])[j]->GetFragmentList(); 
	G4int a1=(*vec)[0]->GetA();
	G4int z1=(*vec)[0]->GetZ();
	G4int a2=(*vec)[1]->GetA();
	G4int z2=(*vec)[1]->GetZ();
	G4int a3=(*vec)[2]->GetA();
	G4int z3=(*vec)[2]->GetZ();
	G4int a4=(*vec)[3]->GetA();
	G4int z4=(*vec)[3]->GetZ();

 	G4cout << "("<<a1<<","<<z1<<")("<<a2<<","<<z2<<")("<<a3<<","<<z3<<")("
	       <<a4<<","<<z4<<") % "; 
      }
      G4cout<<G4endl;
      G4cout<<"----------------------------------------------------------------"
	    << G4endl;
    }
    G4cout << "Total number: " << tot << G4endl;
  }
}

G4bool G4FermiFragmentsPool::IsApplicable(G4int Z, G4int A, G4double mass) const
{
  if(Z >= maxZ || A >= maxA || A <= 0) { return false; }
  // look into pair list
  size_t nz = list2[A].size();
  if(0 < nz) {
    for(size_t j=0; j<nz; ++j) {
      const G4FermiConfiguration* conf = (list2[A])[j];
      if(Z == conf->GetZ() && mass >= conf->GetMass()) { return true; }
    }
  }
  // look into triple list
  nz = list3[A].size();
  if(0 < nz) {
    for(size_t j=0; j<nz; ++j) {
      const G4FermiConfiguration* conf = (list3[A])[j];
      if(Z == conf->GetZ() && mass >= conf->GetMass()) { return true; }
    }
  }
  // look into quartet list
  nz = list4[A].size();
  if(0 < nz) {
    for(size_t j=0; j<nz; ++j) {
      const G4FermiConfiguration* conf = (list4[A])[j];
      if(Z == conf->GetZ() && mass >= conf->GetMass()) { return true; }
    }
  }

  // search in the pool and if found then return vector with one element
  nz = list1[A].size();
  if(0 < nz) {
    for(size_t j=0; j<nz; ++j) {
      const G4FermiConfiguration* conf = (list1[A])[j];
      if(Z == conf->GetZ() && mass >= conf->GetMass()) {
	if(!(*(conf->GetFragmentList()))[0]->IsStable()) { return true; }
      }
    }
  }
  return false;
}

const std::vector<const G4FermiConfiguration*>* 
G4FermiFragmentsPool::GetConfigurationList(G4int Z, G4int A, G4double mass) const
{
  std::vector<const G4FermiConfiguration*>* v = 
    new std::vector<const G4FermiConfiguration*>;
  if(Z >= maxZ || A >= maxA) { return v; }

  //G4cout << "G4FermiFragmentsPool::GetConfigurationList:"
  // << " Z= " << Z << " A= " << A << " Mass(GeV)= " << mass/GeV<< G4endl;

  // look into pair list
  size_t nz = list2[A].size();
  if(0 < nz) {
    for(size_t j=0; j<nz; ++j) {
      const G4FermiConfiguration* conf = (list2[A])[j];
      if(Z == conf->GetZ() && mass >= conf->GetMass()) { 
	v->push_back(conf); 
      }
      //G4cout << "Pair dM(MeV)= " << mass - conf->GetMass() << G4endl; }
    }
  }
  // look into triple list
  nz = list3[A].size();
  if(0 < nz) {
    for(size_t j=0; j<nz; ++j) {
      const G4FermiConfiguration* conf = (list3[A])[j];
      if(Z == conf->GetZ() && mass >= conf->GetMass()) { 
	v->push_back(conf); 
      }
      //G4cout << "Triple dM(MeV)= " << mass - conf->GetMass() << G4endl; }
    }
  }
  // look into quartet list
  nz = list4[A].size();
  if(0 < nz) {
    for(size_t j=0; j<nz; ++j) {
      const G4FermiConfiguration* conf = (list4[A])[j];
      if(Z == conf->GetZ() && mass >= conf->GetMass()) { 
	v->push_back(conf);
      }
      //  G4cout << "Quartet dM(MeV)= " << mass - conf->GetMass() << G4endl; }
    }
  }
  // return if vector not empty
  if(0 < v->size()) { 
    if(verbose > 0) { 
      G4double ExEn= mass - G4NucleiProperties::GetNuclearMass(A,Z);
      G4cout<<"Total number of configurations = "<<v->size()<<" for A= "
	    <<A<<"   Z= "<<Z<<"   E*= "<< ExEn<<" MeV"<<G4endl;
      size_t size_vector_conf = v->size();
      for(size_t jc=0; jc<size_vector_conf; ++jc) {     
	const std::vector<const G4VFermiFragment*>* v_frag = 
	  (*v)[jc]->GetFragmentList();
	size_t size_vector_fragments = v_frag->size();
	G4cout<<size_vector_fragments<<"-body configuration "<<jc+1<<": ";
	for(size_t jf=0;jf<size_vector_fragments;++jf){
	  G4int af= (*v_frag)[jf]->GetA();
	  G4int zf= (*v_frag)[jf]->GetZ();
	  G4double ex=(*v_frag)[jf]->GetExcitationEnergy();
	  G4cout<<"(a="<<af<<", z="<<zf<<", ex="<<ex<<")  ";
	}
	G4cout<<G4endl;
	G4cout<<"-----------------------------------------------------"<<G4endl;
      }
    }
    return v; 
  }

  // search in the pool and if found then return vector with one element
  nz = list1[A].size();
  if(0 < nz) {
    for(size_t j=0; j<nz; ++j) {
      const G4FermiConfiguration* conf = (list1[A])[j];

      //  G4cout << "Single dM(MeV)= " << mass - conf->GetMass() << G4endl; }
      if(Z == conf->GetZ() && mass >= conf->GetMass()) {
	if(!(*(conf->GetFragmentList()))[0]->IsStable()) {
	  v->push_back(conf);

	  if(verbose > 0) { 
	    G4double ExEn= mass -G4NucleiProperties::GetNuclearMass(A,Z);
	    G4cout<<"Only 1 configurations for A= "
		  <<A<<"   Z= "<<Z<<"   E*= "<< ExEn<<" MeV"<<G4endl;
	    const std::vector<const G4VFermiFragment*>* v_frag 
	      = (*v)[0]->GetFragmentList();
	    size_t size_vector_fragments=v_frag->size();
	    G4cout<<"1 Fragment configuration: ";
	    for(size_t jf=0;jf<size_vector_fragments;++jf){
	      G4int af= (*v_frag)[jf]->GetA();
	      G4int zf= (*v_frag)[jf]->GetZ();
	      G4double ex=(*v_frag)[jf]->GetExcitationEnergy();
	      G4cout<<"(a="<<af<<", z="<<zf<<", ex="<<ex<<")  ";
	    }
	    G4cout<<G4endl;
	    G4cout<<"-----------------------------------------------------"<<G4endl;    
	  }
	  return v;
	}
      }
    }
  }
  //failer
  if(verbose > 0) { 
    G4cout << "G4FermiFragmentsPool::GetConfigurationList: WARNING: not "
	   << "able decay fragment Z= " << Z << " A= " << A
	   << " Mass(GeV)= " << mass/GeV<< G4endl;
  }
  return v;
}

G4bool G4FermiFragmentsPool::IsExist(G4int Z, G4int A, 
	      std::vector<const G4VFermiFragment*>& newconf) const
{
  size_t nn = newconf.size();
  G4double mass = 0.0;
  for(size_t i=0; i<nn; ++i) { mass +=  newconf[i]->GetTotalEnergy(); }
  // look into pair list
  if(2 == nn) {
    size_t nz = list2[A].size();
    if(0 < nz) {
      for(size_t j=0; j<nz; ++j) {
	const G4FermiConfiguration* conf = (list2[A])[j];
	if(Z == conf->GetZ() && A == conf->GetA() && 
	   std::fabs(mass - conf->GetMass()) < keV) {return true; }
      }
    }
    return false;
  }
  // look into triple list
  if(3 == nn) {
    size_t nz = list3[A].size();
    if(0 < nz) {
      for(size_t j=0; j<nz; ++j) {
	const G4FermiConfiguration* conf = (list3[A])[j];
	if(Z == conf->GetZ() && A == conf->GetA() && 
	   std::fabs(mass - conf->GetMass()) < keV) { return true; }
      }
    }
    return false;
  }
  // look into quartet list
  if(4 == nn) {
    size_t nz = list4[A].size();
    if(0 < nz) {
      for(size_t j=0; j<nz; ++j) {
	const G4FermiConfiguration* conf = (list4[A])[j];
	if(Z == conf->GetZ() && A == conf->GetA() && 
	   std::fabs(mass - conf->GetMass()) < keV) { return true; }
      }
    }
    return false;
  }
  return false;
}

const G4VFermiFragment* 
G4FermiFragmentsPool::GetFragment(G4int Z, G4int A) const
{
  const G4VFermiFragment* f = 0;
  if(Z < maxZ && A < maxA) { 
    size_t nz = list1[A].size();
    for(size_t j=0; j<nz; ++j) {
      const G4FermiConfiguration* conf = (list1[A])[j];
      if(Z == conf->GetZ()) { 
	f = (*(conf->GetFragmentList()))[0]; 
	break; 
      }
    }
  }
  return f;
}

void G4FermiFragmentsPool::DumpFragment(const G4VFermiFragment* f) const
{
  if(f) {
    G4cout << "Z= " << f->GetZ() << " A= " << f->GetA() 
	   << " Mass(GeV)= " << f->GetFragmentMass()/GeV
	   << " Eexc(MeV)= " << f->GetExcitationEnergy() << G4endl;
  }
}

void G4FermiFragmentsPool::Dump() const
{
  G4cout << "##### List of Fragments in the Fermi Fragment Pool #####" 
	 << G4endl;
  G4int nfrag = fragment_pool.size();
  for(G4int i=0; i<nfrag; ++i) {
    DumpFragment(fragment_pool[i]);
  }
  G4cout << G4endl;
}
