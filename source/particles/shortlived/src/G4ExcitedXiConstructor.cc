// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedXiConstructor.cc,v 1.3 1999-12-15 14:51:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
//      History: first implementation, based on object model of
//      10 oct 1998  H.Kurashige
// ---------------------------------------------------------------


#include "G4ExcitedXiConstructor.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedTable.hh"
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
  G4double r;

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
  G4double r;

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
  G4double r;

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


const char* G4ExcitedXiConstructor::name[] = {
   "xi(1530)", "xi(1690)", "xi(1820)", "xi(1950)", "xi(2030)"
};

const G4double G4ExcitedXiConstructor::mass[] = {
  1.532*GeV, 1.700*GeV, 1.823*GeV, 1.950*GeV,  2.025*GeV 
};

const G4double G4ExcitedXiConstructor::width[] = {
    9.0*MeV,  50.0*MeV,  24.0*MeV,  60.0*MeV,  20.0*MeV
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











