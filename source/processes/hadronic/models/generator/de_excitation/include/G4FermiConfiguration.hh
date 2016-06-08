//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4FermiConfiguration.hh,v 1.9 2002/06/06 17:26:31 larazb Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#ifndef G4FermiConfiguration_h
#define G4FermiConfiguration_h 1

#include "g4std/deque"

#include "globals.hh"
#include "Randomize.hh"
#include "G4VFermiFragment.hh"
#include "G4StableFermiFragment.hh"
#include "G4B9FermiFragment.hh"
#include "G4Be8FermiFragment.hh"
#include "G4He5FermiFragment.hh"
#include "G4Li5FermiFragment.hh"
#include "G4ParticleMomentum.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Fragment.hh"


static const G4int NumberOfFragments = 100;

class G4FermiConfiguration 
{
public:
  G4FermiConfiguration();
  ~G4FermiConfiguration();
  
  G4FermiConfiguration(const G4FermiConfiguration &right);
  
  const G4FermiConfiguration & operator=(const G4FermiConfiguration &right);
  G4bool operator==(const G4FermiConfiguration &right) const;
  G4bool operator!=(const G4FermiConfiguration &right) const;
  
public:

  void Initialize(const G4int max);

  G4bool SplitNucleus(const G4int A, const G4int Z);

  G4double DecayProbability(const G4int A, const G4double TotalE);

  G4FragmentVector * GetFragments(const G4Fragment & theNucleus);


private:

  G4double CoulombBarrier(void);


  G4std::deque<G4LorentzVector*>* FragmentsMomentum(G4double KineticEnergy);
  
  G4double RNKSI(const G4int K);

  G4ParticleMomentum IsotropicVector(const G4double Magnitude = 1.0);


  // Kappa = V/V_0 it is used in calculation of Coulomb energy
  static const G4double Kappa;

  // Nuclear radius r0 (is a model parameter)
  static const G4double r0;



  static G4StableFermiFragment Fragment00;
  static G4StableFermiFragment Fragment01;
  static G4StableFermiFragment Fragment02;
  static G4StableFermiFragment Fragment03;
  static G4StableFermiFragment Fragment04;
  static G4StableFermiFragment Fragment05;
  static G4He5FermiFragment Fragment06;    // He5
  static G4Li5FermiFragment Fragment07;    // Li5
  static G4StableFermiFragment Fragment08;
  static G4StableFermiFragment Fragment09;

  static G4StableFermiFragment Fragment10;
  static G4StableFermiFragment Fragment11;
  static G4StableFermiFragment Fragment12;
  static G4StableFermiFragment Fragment13;
  static G4StableFermiFragment Fragment14;
  static G4StableFermiFragment Fragment15;
  static G4StableFermiFragment Fragment16;
  static G4Be8FermiFragment Fragment17;     // Be8
  static G4StableFermiFragment Fragment18;
  static G4B9FermiFragment Fragment19;  // B9

  static G4StableFermiFragment Fragment20;
  static G4StableFermiFragment Fragment21;
  static G4StableFermiFragment Fragment22;
  static G4StableFermiFragment Fragment23;
  static G4StableFermiFragment Fragment24;
  static G4StableFermiFragment Fragment25;
  static G4StableFermiFragment Fragment26;
  static G4StableFermiFragment Fragment27;
  static G4StableFermiFragment Fragment28;
  static G4StableFermiFragment Fragment29;

  static G4StableFermiFragment Fragment30;
  static G4StableFermiFragment Fragment31;
  static G4StableFermiFragment Fragment32;
  static G4StableFermiFragment Fragment33;
  static G4StableFermiFragment Fragment34;
  static G4StableFermiFragment Fragment35;
  static G4StableFermiFragment Fragment36;
  static G4StableFermiFragment Fragment37;
  static G4StableFermiFragment Fragment38;
  static G4StableFermiFragment Fragment39;

  static G4StableFermiFragment Fragment40;
  static G4StableFermiFragment Fragment41;
  static G4StableFermiFragment Fragment42;
  static G4StableFermiFragment Fragment43;
  static G4StableFermiFragment Fragment44;
  static G4StableFermiFragment Fragment45;
  static G4StableFermiFragment Fragment46;
  static G4StableFermiFragment Fragment47;
  static G4StableFermiFragment Fragment48;
  static G4StableFermiFragment Fragment49;

  static G4StableFermiFragment Fragment50;
  static G4StableFermiFragment Fragment51;
  static G4StableFermiFragment Fragment52;
  static G4StableFermiFragment Fragment53;
  static G4StableFermiFragment Fragment54;
  static G4StableFermiFragment Fragment55;
  static G4StableFermiFragment Fragment56;
  static G4StableFermiFragment Fragment57;
  static G4StableFermiFragment Fragment58;
  static G4StableFermiFragment Fragment59;

  static G4StableFermiFragment Fragment60;
  static G4StableFermiFragment Fragment61;
  static G4StableFermiFragment Fragment62;
  static G4StableFermiFragment Fragment63;
  static G4StableFermiFragment Fragment64;
  static G4StableFermiFragment Fragment65;
  static G4StableFermiFragment Fragment66;
  static G4StableFermiFragment Fragment67;
  static G4StableFermiFragment Fragment68;
  static G4StableFermiFragment Fragment69;

  static G4StableFermiFragment Fragment70;
  static G4StableFermiFragment Fragment71;
  static G4StableFermiFragment Fragment72;
  static G4StableFermiFragment Fragment73;
  static G4StableFermiFragment Fragment74;
  static G4StableFermiFragment Fragment75;
  static G4StableFermiFragment Fragment76;
  static G4StableFermiFragment Fragment77;
  static G4StableFermiFragment Fragment78;
  static G4StableFermiFragment Fragment79;

  static G4StableFermiFragment Fragment80;
  static G4StableFermiFragment Fragment81;
  static G4StableFermiFragment Fragment82;
  static G4StableFermiFragment Fragment83;
  static G4StableFermiFragment Fragment84;
  static G4StableFermiFragment Fragment85;
  static G4StableFermiFragment Fragment86;
  static G4StableFermiFragment Fragment87;
  static G4StableFermiFragment Fragment88;
  static G4StableFermiFragment Fragment89;

  static G4StableFermiFragment Fragment90;
  static G4StableFermiFragment Fragment91;
  static G4StableFermiFragment Fragment92;
  static G4StableFermiFragment Fragment93;
  static G4StableFermiFragment Fragment94;
  static G4StableFermiFragment Fragment95;
  static G4StableFermiFragment Fragment96;
  static G4StableFermiFragment Fragment97;
  static G4StableFermiFragment Fragment98;
  static G4StableFermiFragment Fragment99;


  static G4VFermiFragment * theListOfFragments[NumberOfFragments];


  //  G4VFermiFragment * theConfiguration[MaxConfigSize];

  //  G4int Index[MaxConfigSize];
  G4std::vector<G4int> Index;


    struct DeleteFragment 
    {
	template<typename T>
	void operator()(const T* ptr) const
	    {
		delete ptr;
	    }
    };


};


#endif


