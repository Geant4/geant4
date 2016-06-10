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
// $Id: G4ExcitedXiConstructor.cc 72955 2013-08-14 14:23:14Z gcosmo $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//      History: first implementation, based on object model of
//      10 oct 1998  H.Kurashige
// ---------------------------------------------------------------


#include "G4ExcitedXiConstructor.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"

G4ExcitedXiConstructor::G4ExcitedXiConstructor():
    G4ExcitedBaryonConstructor(NStates, XiIsoSpin)
{

}

G4ExcitedXiConstructor::~G4ExcitedXiConstructor()
{
}

G4DecayTable* G4ExcitedXiConstructor::CreateDecayTable(
                                                 const G4String&  parentName,  
                                                 G4int iIso3, 
                                                 G4int iState,
                                                 G4bool fAnti)
{

  // create decay table
  G4DecayTable* decayTable =  new G4DecayTable();

  G4double br;
  if ( (br=bRatio[iState][XiPi]) >0.0) {
    AddXiPiMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][XiGamma]) >0.0) {
    AddXiGammaMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][LambdaK]) >0.0) {
    AddLambdaKMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][SigmaK]) >0.0) {
    AddSigmaKMode( decayTable, parentName, br, iIso3, fAnti);
  }

  return  decayTable;
}

G4DecayTable*  G4ExcitedXiConstructor::AddXiGammaMode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;
  // 
  G4String daughterH;  
  if (iIso3== +1) {
     daughterH = "xi0";
   } else if (iIso3==-1) {
     daughterH = "xi-";
   }
  if (fAnti) daughterH = "anti_" + daughterH;

  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           daughterH,"gamma");
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedXiConstructor::AddLambdaKMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4bool fAnti)


{
  G4VDecayChannel* mode;

  G4String lambda = "lambda";
  G4String daughterK;
  G4double r = 0.;

  // ------------ Lambda K- ------------ 
  // determine daughters
  if (iIso3 == +1) {
    if (!fAnti) {
      daughterK = "kaon0";
    } else {
      daughterK = "anti_kaon0";
    }  
    r = br;
  } else if (iIso3 == -1) {
    if (!fAnti) {
      daughterK = "kaon-";
    } else {
      daughterK = "kaon+";
    }  
    r = br;
  }
  if (fAnti) lambda = "anti_" + lambda;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        lambda,daughterK);
    // add decay table
    decayTable->Insert(mode);
  }

  return decayTable;
}

G4DecayTable*  G4ExcitedXiConstructor::AddSigmaKMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterH;
  G4String daughterK;
  G4double r = 0.;

  // ------------ Sigma K- ------------ 
  // determine daughters
  if (iIso3== +1) {
    daughterH  = "sigma+";
    r= br/2.;
  } else if (iIso3== -1) {
    daughterH  = "sigma0";
    r = br/2.;
  }
  if (!fAnti) {
    daughterK = "kaon-";
  } else {
    daughterK = "kaon+";
  }
  if (fAnti) daughterH = "anti_" + daughterH;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterH,daughterK);
    // add decay table
    decayTable->Insert(mode);
  }

  // ------------ Sigma K0 ------------ 
  // determine daughters
  if (iIso3 == +1) {
    daughterH  = "sigma0";
    r= br/2.;
  } else if (iIso3 == -1) {
    daughterH  = "sigma-";
     r = br/2.;
  }
  if (!fAnti) {
    daughterK = "anti_kaon0";
  } else {
    daughterK = "kaon0";
  }
  if (fAnti) daughterH = "anti_" + daughterH;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterH,daughterK);
    // add decay table
    decayTable->Insert(mode);
  }

  return decayTable;
}

G4DecayTable*  G4ExcitedXiConstructor::AddXiPiMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterXi;
  G4String daughterPi;
  G4double r = 0.;

  // ------------ Xi pi-  ------------ 
  // determine daughters
  if (iIso3== +1) {
    r = 0.;
  } else if (iIso3 == -1) {
    daughterXi = "xi0";
    r = br/2.;
  }
  if (!fAnti) {
    daughterPi = "pi-";
  } else {
    daughterPi = "pi+";
  }
  if (fAnti) daughterXi = "anti_" + daughterXi;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterXi,daughterPi);
    // add decay table
    decayTable->Insert(mode);
  }
  // ------------ Xi Pi0 ------------ 
  // determine daughters
  if (iIso3== +1) {
    daughterXi = "xi0";
    r = br/2.;
  } else if (iIso3 == -1) {
    daughterXi = "xi-";
    r = br/2.;
  }
  daughterPi = "pi0";
  if (fAnti) daughterXi = "anti_" + daughterXi;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterXi,daughterPi);
    // add decay table
    decayTable->Insert(mode);
  }

  // ------------ XI pi + ------------ 
  // determine daughters
  if (iIso3== +1) {
    daughterXi = "xi-";
    r = br/2.;
  } else if (iIso3==-1) {
    r = 0.;
  }
  if (!fAnti) {
    daughterPi = "pi+";
  } else {
    daughterPi = "pi-";
  }
  if (fAnti) daughterXi = "anti_" + daughterXi;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterXi,daughterPi);
    // add decay table
    decayTable->Insert(mode);
  }

  return decayTable;
}

G4double G4ExcitedXiConstructor::GetMass(G4int iState, G4int iso3)
{ 
  G4double fm = mass[iState];
  if ( iState==0 ) {
    if (iso3== -1) fm = 1.5350*GeV; // xi-
  }
  return fm;
}

G4double G4ExcitedXiConstructor::GetWidth(G4int iState, G4int iso3)
{
  G4double fw = width[iState];
  if ( iState==0 ) {
    if (iso3== -1) fw = 9.9*MeV; // xi-
  }
  return fw;
}

const char* G4ExcitedXiConstructor::name[] = {
   "xi(1530)", "xi(1690)", "xi(1820)", "xi(1950)", "xi(2030)"
};

const G4double G4ExcitedXiConstructor::mass[] = {
 1.5318*GeV, 1.690*GeV, 1.823*GeV, 1.950*GeV,  2.025*GeV 
};

const G4double G4ExcitedXiConstructor::width[] = {
    9.1*MeV,  50.0*MeV,  24.0*MeV,  60.0*MeV,  20.0*MeV
};

const G4int G4ExcitedXiConstructor::iSpin[] = {
    3,   3,   3,   3,   5
};

const G4int G4ExcitedXiConstructor::iParity[] = {
   +1,  +1,   -1,  -1,  +1
};


const G4int G4ExcitedXiConstructor::encodingOffset[] = {
      0, 20000, 10000, 30000, 10000
};

const G4double G4ExcitedXiConstructor::bRatio[ G4ExcitedXiConstructor::NStates ][ G4ExcitedXiConstructor::NumberOfDecayModes] = 
{
   {  0.98, 0.02,  0.0,  0.0}, 
   {  0.10,  0.0, 0.70, 0.20}, 
   {  0.15,  0.0, 0.70, 0.15}, 
   {  0.25,  0.0, 0.50, 0.25}, 
   {  0.10,  0.0, 0.20, 0.70}
};











