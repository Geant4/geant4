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
// $Id: G4ExcitedLambdaConstructor.cc 83749 2014-09-12 12:14:59Z gcosmo $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//      History: first implementation, based on object model of
//      10 oct 1998  H.Kurashige
// ---------------------------------------------------------------


#include "G4ExcitedLambdaConstructor.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"


G4ExcitedLambdaConstructor::G4ExcitedLambdaConstructor():
    G4ExcitedBaryonConstructor(NStates, LambdaIsoSpin)
{

}

G4ExcitedLambdaConstructor::~G4ExcitedLambdaConstructor()
{
}

G4DecayTable* G4ExcitedLambdaConstructor::CreateDecayTable(
						 const G4String&  parentName,  
						 G4int iIso3, 
						 G4int iState,
						 G4bool fAnti)
{
  // create decay table
  G4DecayTable* decayTable =  new G4DecayTable();

  G4double br;
  if ( (br=bRatio[iState][NK]) >0.0) {
    AddNKMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][NKStar]) >0.0) {
    AddNKStarMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][SigmaPi]) >0.0) {
    AddSigmaPiMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][SigmaStarPi]) >0.0) {
    AddSigmaStarPiMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][LambdaGamma]) >0.0) {
    AddLambdaGammaMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][LambdaEta]) >0.0) {
    AddLambdaEtaMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][LambdaOmega]) >0.0) {
    AddLambdaOmegaMode( decayTable, parentName, br, iIso3, fAnti);
  }

  return  decayTable;
}

G4DecayTable*  G4ExcitedLambdaConstructor::AddLambdaGammaMode( 
				   G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int , G4bool fAnti)
{
  G4VDecayChannel* mode;

  // 
  G4String lambda = "lambda";  
  if (fAnti) lambda = "anti_" + lambda;

  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           lambda,"gamma");
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}
G4DecayTable*  G4ExcitedLambdaConstructor::AddLambdaEtaMode( 
				   G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int , G4bool fAnti)
{
  G4VDecayChannel* mode;

  // 
  G4String lambda = "lambda";  
  if (fAnti) lambda = "anti_" + lambda;

  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           lambda,"eta");
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedLambdaConstructor::AddLambdaOmegaMode( 
				   G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int , G4bool fAnti)
{
  G4VDecayChannel* mode;

  // 
  G4String lambda = "lambda";  
  if (fAnti) lambda = "anti_" + lambda;

  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           lambda,"omega");
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedLambdaConstructor::AddNKMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int , G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterN;
  G4String daughterK;

  // ------------ N K- ------------ 
  // determine daughters
  daughterN  = "proton";
  if (!fAnti) {
    daughterK = "kaon-";
  } else {
    daughterK = "kaon+";
  }
  if (fAnti) daughterN = "anti_" + daughterN;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/2.0, 2,
                                           daughterN,daughterK);
  // add decay table
  decayTable->Insert(mode);

  // ------------ N K0 ------------ 
  // determine daughters
  daughterN  = "neutron";
  if (!fAnti) {
    daughterK = "anti_kaon0";
  } else {
    daughterK = "kaon0";
  }
  if (fAnti) daughterN = "anti_" + daughterN;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/2.0, 2,
                                           daughterN,daughterK);
  // add decay table
  decayTable->Insert(mode);


  return decayTable;
}

G4DecayTable*  G4ExcitedLambdaConstructor::AddNKStarMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int , G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterN;
  G4String daughterK;

  // ------------ N K- ------------ 
  // determine daughters
  daughterN  = "proton";
  if (!fAnti) {
    daughterK = "k_star-";
  } else {
    daughterK = "k_star+";
  }
  if (fAnti) daughterN = "anti_" + daughterN;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/2.0, 2,
                                           daughterN,daughterK);
  // add decay table
  decayTable->Insert(mode);

  // ------------ N K0 ------------ 
  // determine daughters
  daughterN  = "neutron";
  if (!fAnti) {
    daughterK = "anti_k_star0";
  } else {
    daughterK = "k_star0";
  }
  if (fAnti) daughterN = "anti_" + daughterN;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/2.0, 2,
                                           daughterN,daughterK);
  // add decay table
  decayTable->Insert(mode);


  return decayTable;
}

G4DecayTable*  G4ExcitedLambdaConstructor::AddSigmaPiMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int , G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterSigma;
  G4String daughterPi;

  // ------------ Sigma+ pi - ------------ 
  // determine daughters
  daughterSigma = "sigma+";
  if (!fAnti) {
    daughterPi = "pi-";
  } else {
    daughterPi = "pi+";
  }
  if (fAnti) daughterSigma = "anti_" + daughterSigma;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/3.0, 2,
                                           daughterSigma,daughterPi);
  // add decay table
  decayTable->Insert(mode);

   // ------------ Sigma0 Pi0 ------------ 
  // determine daughters
  daughterSigma  = "sigma0";
  daughterPi = "pi0";

  if (fAnti) daughterSigma = "anti_" + daughterSigma;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/3.0, 2,
                                           daughterSigma,daughterPi);

  // add decay table
  decayTable->Insert(mode);

  // ------------ Sigma- pi + ------------ 
  // determine daughters
  daughterSigma = "sigma-";
  if (!fAnti) {
    daughterPi = "pi+";
  } else {
    daughterPi = "pi-";
  }
  if (fAnti) daughterSigma = "anti_" + daughterSigma;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/3.0, 2,
                                           daughterSigma,daughterPi);
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}


G4DecayTable*  G4ExcitedLambdaConstructor::AddSigmaStarPiMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int , G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterSigma;
  G4String daughterPi;

  // ------------ Sigma+ pi - ------------ 
  // determine daughters
  daughterSigma = "sigma(1385)+";
  if (!fAnti) {
    daughterPi = "pi-";
  } else {
    daughterPi = "pi+";
  }
  if (fAnti) daughterSigma = "anti_" + daughterSigma;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/3.0, 2,
                                           daughterSigma,daughterPi);
  // add decay table
  decayTable->Insert(mode);

   // ------------ Sigma0 Pi0 ------------ 
  // determine daughters
  daughterSigma  = "sigma(1385)0";
  daughterPi = "pi0";

  if (fAnti) daughterSigma = "anti_" + daughterSigma;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/3.0, 2,
                                           daughterSigma,daughterPi);

  // add decay table
  decayTable->Insert(mode);

  // ------------ Sigma- pi + ------------ 
  // determine daughters
  daughterSigma = "sigma(1385)-";
  if (!fAnti) {
    daughterPi = "pi+";
  } else {
    daughterPi = "pi-";
  }
  if (fAnti) daughterSigma = "anti_" + daughterSigma;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/3.0, 2,
                                           daughterSigma,daughterPi);
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

const char* G4ExcitedLambdaConstructor::name[] = {
  "lambda(1405)","lambda(1520)","lambda(1600)","lambda(1670)","lambda(1690)", 
  "lambda(1800)","lambda(1810)","lambda(1820)","lambda(1830)","lambda(1890)",
  "lambda(2100)","lambda(2110)"
};

const G4double G4ExcitedLambdaConstructor::mass[] = {
  1.4051*GeV,1.5195*GeV, 1.600*GeV, 1.670*GeV,  1.690*GeV, 
   1.800*GeV, 1.810*GeV, 1.820*GeV, 1.830*GeV,  1.890*GeV, 
   2.100*GeV, 2.110*GeV
};

const G4double G4ExcitedLambdaConstructor::width[] = {
   50.5*MeV,  15.6*MeV, 150.0*MeV,  35.0*MeV,  60.0*MeV,
  300.0*MeV, 150.0*MeV,  80.0*MeV,  95.0*MeV, 100.0*MeV,
  200.0*MeV, 200.0*MeV
};

const G4int G4ExcitedLambdaConstructor::iSpin[] = {
    1,   3,   1,   1,   3,
    1,   1,   5,   5,   3,
    7,   5
};

const G4int G4ExcitedLambdaConstructor::iParity[] = {
  -1,  -1,   +1,  -1,  -1,
  -1,  +1,   +1,  -1,  +1,
  -1,  +1 
};

const G4int G4ExcitedLambdaConstructor::encodingOffset[] = {
  10000,     0, 20000, 30000, 10000,
  40000, 50000,     0, 10000, 20000, 
      0, 20000
};

const G4double G4ExcitedLambdaConstructor::bRatio[ G4ExcitedLambdaConstructor::NStates ][ G4ExcitedLambdaConstructor::NumberOfDecayModes] = 
{
   {   0.0,  0.0,  1.0,  0.0,  0.0,   0.0,   0.0}, 
   {  0.45,  0.0, 0.43, 0.11, 0.01,   0.0,   0.0}, 
   {  0.35,  0.0, 0.65,  0.0,  0.0,   0.0,   0.0}, 
   {  0.20,  0.0, 0.50,  0.0,  0.0,  0.30,   0.0}, 
   {  0.25,  0.0, 0.45, 0.30,  0.0,   0.0,   0.0}, 
   {  0.40, 0.20, 0.20, 0.20,  0.0,   0.0,   0.0}, 
   {  0.35, 0.45, 0.15, 0.05,  0.0,   0.0,   0.0}, 
   {  0.73,  0.0, 0.16, 0.11,  0.0,   0.0,   0.0}, 
   {  0.10,  0.0, 0.70, 0.20,  0.0,   0.0,   0.0}, 
   {  0.37, 0.21, 0.11, 0.31,  0.0,   0.0,   0.0}, 
   {  0.35, 0.20, 0.05, 0.30,  0.0,  0.02,  0.08}, 
   {  0.25, 0.45, 0.30,  0.0,  0.0,   0.0,   0.0}
};
















