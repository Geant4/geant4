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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4Mars5GeV.hh,v 1.1 2001-12-13 14:59:46 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
// 
// ------------------------------------------------------------
//   First Implemention    17 Nov. 1998  M.Asai, H.Kurahige
// 
//   History:
//   modified as hadronic model 28 Oct 2001 N.Kanaya
// ------------------------------------------------------------
//  Class Description
//  This is a Event Biasing mechanism based on MARS code
//   This model is applicable to 
//   proton/neutron/pi+-/K+-/gamma/anti_proton
//   with energy < 5.0GeV
//*
//  Original code is MARS13 written by Nikolai Mokhov (FNAL)
//**************************************************************
//*   MARS13: 9. hA EVENT GENERATOR:
//*     Copyright Nikolai Mokhov (Fermilab)
//*
//*     LAST CHANGE: 14-NOV-1998
//**************************************************************
//*     Copyright Nikolai Mokhov (Fermilab)
//*
//*     MARS13(98)
//*
//*     INCLUSIVE HADRON(photon)-NUCLEUS VERTEX AT E < 5 GEV !!!
//*     THREE WEIGHTED HADRONS IN FINAL STATE:     !!!
//*     IP+A -> N/P(CASC)+ PI+/PI-(K+/K-) + PI0
//

#ifndef G4Mars5GeV_h
#define G4Mars5GeV_h 1
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4ElementTable.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4LightMedia.hh"
#include "G4Step.hh"
#include "G4TrackStatus.hh"
#include "G4InelasticInteraction.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

class G4Mars5GeV : public G4InelasticInteraction
{
public: // with description
  // constructors
  G4Mars5GeV(); 
  
  //destructor
  ~G4Mars5GeV() { ; };
  
  // virtual methods derived from G4VEvtBiasMechanism
  G4VParticleChange* ApplyYourself( const G4Track &aTrack,
				    G4Nucleus &targetNucleus );

 private:
  void Treem5(); 
  // This is the mothod which invoke MARS 
  // secondary information will be filled up 
  
  G4bool CoulombBarrier(G4int pType, G4double  pE);
  // Check if coulomb barrier exists

  void CreateNucleon(G4int ib, G4int pType, G4double  pE);
  void CreatePion(G4int ib, G4int pType, G4double  pE);
  void CreatePionZero(G4int ib, G4int pType, G4double  pE);
  // Create secondary particles and add them into the list 
   
  void     AddSecondary(); 
  // Add a secondary particle  into the list
   
  G4double SelBS(G4int pType, G4double aNucl, G4double zNucl);
  // Calculate weight of secondary
  
  G4double D2N2(G4int    pType,       G4double incidentE, 
	        G4double prodE,       G4double tin,
	        G4int    reacType,    G4int    proType,
                G4double ai,          G4double z);
  // Calculate Hadron Inclusive Yield  

  G4double Rkaon(G4int ib, G4int jp, G4double eRaw);
  // Calculate energy dependent K/pi ratio

  void Trans(G4ThreeVector* d1, G4ThreeVector* d2);
  // Direction cosine transformation using selec2(cs,ss,ch,sh)

  G4double GetMaxWeight() const;
  G4double GetMinWeight() const;
  void     SetMaxWeight(G4double);
  void     SetMinWeight(G4double);

 public:
  enum {FastVectorSize = 16};
  typedef  G4FastVector<G4DynamicParticle ,FastVectorSize> G4MarsSecondaryVector;

 private:
  // weight max/min
  G4double             maxWeight;
  G4double             minWeight; 

  // information of secondary 
  G4int                 numberOfSecondaries;
  G4double              weightOfSecondaries[FastVectorSize];
  G4MarsSecondaryVector secondaries;
    
 private:
  G4bool IsApplicable(G4int marsEncoding) const;

 private:
  // Particle Change

  // Particle Table
  G4ParticleTable* theParticleTable;

  // particle encoding for MARS
  enum { MarsUndefined =0,
         MarsP, MarsN, MarsPIplus, MarsPIminus, MarsKplus, MarsKminus,
         MarsMUplus, MarsMUminus, MarsGAM, MarsEminus, MarsEplus, MarsAP,
         MarsPI0, MarsD, MarsT, MarsHe3, MarsHe4 };

  G4int GetMarsEncoding(const G4ParticleDefinition* )const;
  const G4String& GetParticleName(G4int marsEncoding) const;
  G4ParticleDefinition* GetParticleDefinition(G4int marsEncoding) const;   
  G4double ProtonMass;
 
 private:
  // incident information
  G4double                   incidentWeight;
  const G4DynamicParticle*   incidentParticle; 
  G4int                      incidentMarsEncoding;

  void GetTargetNuclei(const G4Material*);  //fill up fANucle and fZnucl
  G4double            fANucl, fZNucl;       // target nucleus

 private:  
  // these class is to define common blocks in original code

  class Selec1 
  {
   public:
    G4double   Einc,  EN,   V,  V10;
    G4int      Treac, Tprod;
  } selec1;
  class Selec2
  {
   public:
    G4double   Cs, Ss, Ch, Sh;
  } selec2;
  class Selec3
  {
   public:
    G4double   Eth, Emax, Sqs, X, Pt2, Pt, P;
  } selec3;
};

#include "G4ParticleDefinition.hh"
inline
 G4double G4Mars5GeV::GetMaxWeight() const
{
 return maxWeight;
}
inline
 G4double G4Mars5GeV::GetMinWeight() const
{
 return minWeight;
}

inline
 void  G4Mars5GeV::SetMaxWeight(G4double value)
{
 maxWeight = value;
}
inline
 void  G4Mars5GeV::SetMinWeight(G4double value)
{
 minWeight = value;
}

inline
 G4int G4Mars5GeV::GetMarsEncoding(const G4ParticleDefinition* particle) const{
  const G4String& name = particle->GetParticleName();
  G4int encoding = MarsUndefined;
  if (name == "proton") {
	encoding = MarsP;
  } else if  (name == "neutron") {
	encoding = MarsN; 
  } else if  (name == "pi+") {
	encoding = MarsPIplus; 
  } else if  (name == "pi-") {
	encoding = MarsPIminus; 
  } else if  (name == "kaon+") {
	encoding = MarsKplus; 
  } else if  (name == "kaon-") {
	encoding = MarsKminus; 
  } else if  (name == "mu+") {
	encoding = MarsMUplus; 
  } else if  (name == "mu-") {
	encoding = MarsMUminus; 
  } else if  (name == "gamma") {
	encoding = MarsGAM; 
  } else if  (name == "e+") {
	encoding = MarsEplus; 
  } else if  (name == "e-") {
	encoding = MarsEminus; 
  } else if  (name == "anti_proton") {
	encoding = MarsAP; 
  } else if  (name == "pi0") {
	encoding = MarsPI0; 
  } else if  (name == "deuteron") {
	encoding = MarsD; 
  } else if  (name == "triton") {
	encoding = MarsT; 
  } else if  (name == "He3") {
	encoding = MarsHe3; 
  } else if  (name == "alpha") {
	encoding = MarsHe4; 
  }
  return encoding;
}

inline 
 const G4String& G4Mars5GeV::GetParticleName(G4int encoding) const
{
  static G4String name;
  name = "None";
  switch (encoding) 
  {
    case MarsP:
      name = "proton";
      break;
    case MarsN:
      name = "neutron";
      break;
    case MarsPIplus:
      name = "pi+";
      break;
    case MarsPIminus:
      name = "pi-";
      break;
    case MarsKplus:
      name = "kaon+";
      break;
    case MarsKminus:
      name = "kaon-";
      break;
    case MarsMUplus:
      name = "mu+";
      break;
    case MarsMUminus:
      name = "mu-";
      break;
    case MarsGAM:
      name = "gamma";
      break;
    case MarsEplus:
      name = "e+";
      break;
    case MarsEminus:
      name = "e-";
      break;
    case MarsAP:
      name = "anti_proton";
      break;
    case MarsPI0:
      name = "pi0";
      break;
    case MarsD:
      name = "deuteron";
      break;
    case MarsT:
      name = "triton";
      break;
    case MarsHe3:
      name = "He3";
      break;
    case MarsHe4:
      name = "alpha";
      break;
    default:
      break;
  }
  return name;
}   

inline 
G4ParticleDefinition* G4Mars5GeV::GetParticleDefinition(G4int encoding) const
{
  G4String name = GetParticleName(encoding);
  G4ParticleDefinition* particle = 0;
  if (name != "None") {
	particle = theParticleTable->FindParticle(name);
  }
  return particle;
}   

inline
 G4bool G4Mars5GeV::IsApplicable(G4int marsEncoding) const
{
  return ( ((marsEncoding!=MarsUndefined) && (marsEncoding<=MarsKminus) )|| 
	   (marsEncoding==MarsGAM) ||
	   (marsEncoding==MarsAP)      );
}

#endif
