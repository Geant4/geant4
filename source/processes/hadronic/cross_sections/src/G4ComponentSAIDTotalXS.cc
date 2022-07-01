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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:    G4ComponentSAIDTotalXS
//
// Authors:  G.Folger, V.Ivanchenko, D.Wright
//
// Modifications:
//

#include "G4ComponentSAIDTotalXS.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsFreeVector.hh"

const G4String G4ComponentSAIDTotalXS::fnames[13] = {
  "","pp","np","pip","pim",
  "pin","pie",
  "gp_pi0p","gp_pi+n","gn_pi-p","gn_pi0n","gp_etap","gp_etapp"
};

#ifdef G4MULTITHREADED
  G4Mutex G4ComponentSAIDTotalXS::saidXSMutex = G4MUTEX_INITIALIZER;
#endif

G4ComponentSAIDTotalXS::G4ComponentSAIDTotalXS() 
  : G4VComponentCrossSection("xsSAID")
{
  for(G4int i=0; i<numberOfSaidXS; ++i) {
    elastdata[i] = nullptr;
    inelastdata[i] = nullptr;
  }
}

G4ComponentSAIDTotalXS::~G4ComponentSAIDTotalXS()
{
  for(G4int i=0; i<numberOfSaidXS; ++i) {
    if(elastdata[i]) {
      delete elastdata[i];
      elastdata[i] = nullptr;
    }
    if(inelastdata[i]) {
      delete inelastdata[i];
      inelastdata[i] = nullptr;
    }
  }
}

G4double 
G4ComponentSAIDTotalXS::GetTotalElementCrossSection(
      const G4ParticleDefinition* part,
      G4double, G4int Z, G4double N)
{
  PrintWarning(part,0,Z,G4lrint(N),
	       "G4ComponentSAIDTotalXS::GetTotalElementCrossSection",
	       "Method is not implemented");
  return 0.0; 
}

G4double 
G4ComponentSAIDTotalXS::GetTotalIsotopeCrossSection(
      const G4ParticleDefinition* part,
      G4double kinEnergy, G4int Z, G4int N)
{
  return GetInelasticIsotopeCrossSection(part,kinEnergy,Z,N)
    + GetElasticIsotopeCrossSection(part,kinEnergy,Z,N);
}

G4double 
G4ComponentSAIDTotalXS::GetInelasticElementCrossSection(
      const G4ParticleDefinition* part,
      G4double, G4int Z, G4double N)
{
  PrintWarning(part,0,Z,G4lrint(N),
	       "G4ComponentSAIDTotalXS::GetTotalElementCrossSection",
	       "Method is not implemented");
  return 0.0; 
}

G4double 
G4ComponentSAIDTotalXS::GetInelasticIsotopeCrossSection(
      const G4ParticleDefinition* part,
      G4double kinEnergy, G4int Z, G4int N)
{
  G4double cross = 0.0;
  G4SAIDCrossSectionType tp = GetType(part,0,Z,N);
  if(saidUnknown != tp) {
    G4int idx = G4int(tp);
    if(!inelastdata[idx]) { Initialise(tp); }
    if(inelastdata[idx]) { 
      cross = (inelastdata[idx])->Value(kinEnergy);
    }
  }
  return cross;
}

G4double 
G4ComponentSAIDTotalXS::GetElasticElementCrossSection(
      const G4ParticleDefinition* part,
      G4double, G4int Z, G4double N)
{
  PrintWarning(part,0,Z,G4lrint(N),
	       "G4ComponentSAIDTotalXS::GetTotalElementCrossSection",
	       "Method is not implemented");
  return 0.0; 
}

G4double 
G4ComponentSAIDTotalXS::GetElasticIsotopeCrossSection(
      const G4ParticleDefinition* part,
      G4double kinEnergy, G4int Z, G4int N)
{
  G4double cross = 0.0;
  G4SAIDCrossSectionType tp = GetType(part,0,Z,N);
  if(saidUnknown != tp) {
    G4int idx = G4int(tp);
    if(!elastdata[idx]) { Initialise(tp); }
    if(elastdata[idx]) { 
      cross = (elastdata[idx])->Value(kinEnergy); 
    }
  }
  return cross;
}

G4double 
G4ComponentSAIDTotalXS::GetChargeExchangeCrossSection(
     const G4ParticleDefinition* prim,
     const G4ParticleDefinition* sec,
     G4double kinEnergy, G4int Z, G4int N)
{
  G4double cross = 0.0;
  G4SAIDCrossSectionType tp = GetType(prim,sec,Z,N);
  if(saidUnknown != tp) {
    G4int idx = G4int(tp);
    if(!inelastdata[idx]) { Initialise(tp); }
    if(inelastdata[idx]) { 
      cross = (inelastdata[idx])->Value(kinEnergy); 
    }
  }
  return cross;
}

void G4ComponentSAIDTotalXS::Description(std::ostream&) const
{
}

G4SAIDCrossSectionType 
G4ComponentSAIDTotalXS::GetType(const G4ParticleDefinition* prim,
				const G4ParticleDefinition* sec,
				G4int Z, G4int N)
{
  G4SAIDCrossSectionType type = saidUnknown;
  if(1 == N && prim) {
    G4int code = prim->GetPDGEncoding();

    // only gamma + N x-sections available
    if(0 == Z && sec && 22 == code) {
      G4int code1 = sec->GetPDGEncoding();
      if(-211 == code1)     { type = saidGN_PINP; }
      else if(111 == code1) { type = saidGN_PI0N; }

      // x-sections off proton
    } else if(1 == Z) {
      if(sec) {
	G4int code1 = sec->GetPDGEncoding();
        if(-211 == code) {
          if(111 == code1)      { type = saidPINP_PI0N; }
          else if(221 == code1) { type = saidPINP_ETAN; }

	} else if(22 == code) {
          if(111 == code1)      { type = saidGP_PI0P; }
          else if(211 == code1) { type = saidGP_PIPN; }
          else if(221 == code1) { type = saidGP_ETAP; }
          else if(331 == code1) { type = saidGP_ETAPP; }
	}
      } else {
        if(2212 == code)        { type = saidPP; } 
        else if(2112 == code)   { type = saidNP; } 
        else if(211 == code)    { type = saidPIPP; } 
        else if(-211 == code)   { type = saidPINP; } 
      }
    }
  }
  //G4cout << "G4ComponentSAIDTotalXS::Type= " << type << G4endl;
  return type;
}

void G4ComponentSAIDTotalXS::Initialise(G4SAIDCrossSectionType tp)
{
  G4int idx = G4int(tp);
#ifdef G4MULTITHREADED
  G4MUTEXLOCK(&saidXSMutex);
  if(!inelastdata[idx]) { 
#endif
    // check environment variable 
    // Build the complete string identifying the file with the data set
    const char* path = G4FindDataDir("G4SAIDXSDATA");
    if (!path){
      G4Exception("G4ComponentSAIDTotalXS::Initialise(..)","had013",
		  FatalException,
		  "Environment variable G4SAIDXSDATA is not defined");
      return;
    }
    if(idx <= 4) {
      elastdata[idx] = new G4PhysicsFreeVector(true);
      inelastdata[idx] = new G4PhysicsFreeVector(true);
      ReadData(idx,elastdata[idx],path,"_el.dat");
      ReadData(idx,inelastdata[idx],path,"_in.dat");
    } else {
      inelastdata[idx] = new G4PhysicsFreeVector();
      ReadData(idx,inelastdata[idx],path,".dat");
    }
#ifdef G4MULTITHREADED
  }
  G4MUTEXUNLOCK(&saidXSMutex);
#endif
}

void G4ComponentSAIDTotalXS::ReadData(G4int index, 
				      G4PhysicsVector* v,
				      const G4String& ss1, 
				      const G4String& ss2)
{
  std::ostringstream ost;
  ost << ss1 << "/" << fnames[index] << ss2;
  std::ifstream filein(ost.str().c_str());
  if (!(filein)) {
    G4ExceptionDescription ed;
    ed << "Data file <" << ost.str().c_str()
       << "> is not opened!";
    G4Exception("G4ComponentSAIDTotalXS::ReadData(..)","had014",
                FatalException, ed, "Check G4SAIDXSDATA");
  } else {
    if(GetVerboseLevel() > 1) {
      G4cout << "File " << ost.str() 
             << " is opened by G4ComponentSAIDTotalXS" << G4endl;
    }
    // retrieve data from DB
    v->Retrieve(filein, true);
    v->ScaleVector(CLHEP::MeV,CLHEP::millibarn);
    v->FillSecondDerivatives();
  } 
}

void 
G4ComponentSAIDTotalXS::PrintWarning(const G4ParticleDefinition* prim,
				     const G4ParticleDefinition* sec,
				     G4int Z, G4int N,
				     const G4String& ss1,
				     const G4String& ss2)
{
  G4cout << ss1 << ": " << ss2 << G4endl;
  G4cout << "For Z= " << Z << " N= " << N << " of ";
  if(prim) { G4cout << prim->GetParticleName() << " "; }
  if(sec) { G4cout << " x-section to " << sec->GetParticleName(); }
  G4cout << G4endl;
}
