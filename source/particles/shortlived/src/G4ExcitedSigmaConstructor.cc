// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedSigmaConstructor.cc,v 1.2 1999-06-08 07:33:30 kurasige Exp $
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


#include "G4ExcitedSigmaConstructor.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"

G4ExcitedSigmaConstructor::G4ExcitedSigmaConstructor():
    G4ExcitedBaryonConstructor(NStates, SigmaIsoSpin)
{

}

G4ExcitedSigmaConstructor::~G4ExcitedSigmaConstructor()
{
}

G4DecayTable* G4ExcitedSigmaConstructor::CreateDecayTable(
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

  if ( (br=bRatio[iState][LambdaPi]) >0.0) {
    AddLambdaPiMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][SigmaEta]) >0.0) {
    AddSigmaEtaMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][LambdaStarPi]) >0.0) {
    AddLambdaStarPiMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][DeltaK]) >0.0) {
    AddDeltaKMode( decayTable, parentName, br, iIso3, fAnti);
  }

  return  decayTable;
}

G4DecayTable*  G4ExcitedSigmaConstructor::AddSigmaEtaMode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;
  // 
  G4String daughterH;  
  if (iIso3== +2) {
     daughterH = "sigma+";
   } else if (iIso3== 0) {
     daughterH = "sigma0";
   } else if (iIso3== -2) {
     daughterH = "sigma-";
   }
  if (fAnti) daughterH = "anti_" + daughterH;

  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           daughterH,"eta");
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedSigmaConstructor::AddNKMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4bool fAnti)


{
  G4VDecayChannel* mode;

  G4String daughterN;
  G4String daughterK;
  G4double r;

  // ------------ N K- ------------ 
  // determine daughters
  if (iIso3== +2) {
    r=0.;
  } else if (iIso3== 0) {
    daughterN  = "proton";
    r = br/2.;
  } else if (iIso3== -2) {
    daughterN  = "neutron";
    r = br;
  }
  if (!fAnti) {
    daughterK = "kaon-";
  } else {
    daughterK = "kaon+";
  }  
  if (fAnti) daughterN = "anti_" + daughterN;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterN,daughterK);
    // add decay table
    decayTable->Insert(mode);
  }

  // ------------ N K0 ------------ 
  // determine daughters
  if (iIso3== +2) {
    daughterN  = "proton";
    r=br;
  } else if (iIso3== 0) {
    daughterN  = "neutron";
     r = br/2.;
  } else if (iIso3== -2) {
    r = 0.;
  }
  if (!fAnti) {
    daughterK = "anti_kaon0";
  } else {

    daughterK = "kaon0";
  }
  if (fAnti) daughterN = "anti_" + daughterN;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterN,daughterK);
    // add decay table
    decayTable->Insert(mode);
  }

  return decayTable;
}

G4DecayTable*  G4ExcitedSigmaConstructor::AddDeltaKMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterN;
  G4String daughterK;
  G4double r;

  // ------------ N K- ------------ 
  // determine daughters
  if (iIso3== +2) {
    daughterN  = "delta++";
    r=0.75*br;
  } else if (iIso3== 0) {
    daughterN  = "delta+";
    r = br/2.;
  } else if (iIso3== -2) {
    daughterN  = "delta0";
    r = 0.25*br;
  }
  if (!fAnti) {
    daughterK = "kaon-";
  } else {
    daughterK = "kaon+";
  }
  if (fAnti) daughterN = "anti_" + daughterN;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterN,daughterK);
    // add decay table
    decayTable->Insert(mode);
  }

  // ------------ N K0 ------------ 
  // determine daughters
  if (iIso3== +2) {
    daughterN  = "delta+";
    r=0.25*br;
  } else if (iIso3== 0) {
    daughterN  = "delta0";
     r = br/2.;
  } else if (iIso3== -2) {
    daughterN  = "delta-";
    r=0.75*br;
  }
  if (!fAnti) {
    daughterK = "anti_kaon0";
  } else {
    daughterK = "kaon0";
  }
  if (fAnti) daughterN = "anti_" + daughterN;
  if (r>0.) {


    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterN,daughterK);
    // add decay table
    decayTable->Insert(mode);
  }

  return decayTable;
}


G4DecayTable*  G4ExcitedSigmaConstructor::AddNKStarMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterN;
  G4String daughterK;
  G4double r;

  // ------------ N K- ------------ 
  // determine daughters
  if (iIso3== +2) {

    r=0.;
  } else if (iIso3== 0) {
    daughterN  = "proton";
    r = br/2.;
  } else if (iIso3== -2) {
    daughterN  = "neutron";
    r = br;
  }
  if (!fAnti) {
    daughterK = "k_star-";
  } else {
    daughterK = "k_star+";
  }
  if (fAnti) daughterN = "anti_" + daughterN;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterN,daughterK);
    // add decay table
    decayTable->Insert(mode);
  }

  // ------------ N K0 ------------ 

  // determine daughters
  if (iIso3== +2) {
    daughterN  = "proton";
    r=br;
  } else if (iIso3== 0) {
    daughterN  = "neutron";
    r = br/2.;
  } else if (iIso3== -2) {
    r = 0.;
  }
  if (!fAnti) {
    daughterK = "anti_k_star0";
  } else {
    daughterK = "k_star0";
  }
  if (fAnti) daughterN = "anti_" + daughterN;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                           daughterN,daughterK);
  // add decay table
  decayTable->Insert(mode);


  return decayTable;
}

G4DecayTable*  G4ExcitedSigmaConstructor::AddSigmaPiMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterSigma;
  G4String daughterPi;
  G4double r;

  // ------------ Sigma+ pi - ------------ 
  // determine daughters
  if (iIso3== +2) {
    r = 0.;
  } else if (iIso3== 0) {
    daughterSigma = "sigma+";
    r = br/2.;
  } else if (iIso3== -2) {
    daughterSigma = "sigma0";
    r = br/2.;
  }
  if (!fAnti) {
    daughterPi = "pi-";
  } else {
    daughterPi = "pi+";
  }
  if (fAnti) daughterSigma = "anti_" + daughterSigma;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterSigma,daughterPi);
    // add decay table
    decayTable->Insert(mode);
  }
  // ------------ Sigma0 Pi0 ------------ 
  // determine daughters
  if (iIso3== +2) {
    daughterSigma = "sigma+";
    r = br/2.;
  } else if (iIso3== 0) {
    r = 0.;
  } else if (iIso3== -2) {
    daughterSigma = "sigma-";


    r = br/2.;
  }
  daughterPi = "pi0";
  if (fAnti) daughterSigma = "anti_" + daughterSigma;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterSigma,daughterPi);
    // add decay table
    decayTable->Insert(mode);
  }

  // ------------ Sigma- pi + ------------ 
  // determine daughters
  if (iIso3== +2) {
    daughterSigma = "sigma0";
    r = br/2.;
  } else if (iIso3== 0) {
    daughterSigma = "sigma-";
    r = br/2.;
  } else if (iIso3== -2) {
    r = 0.;
  }
  if (!fAnti) {
    daughterPi = "pi+";
  } else {
    daughterPi = "pi-";
  }
  if (fAnti) daughterSigma = "anti_" + daughterSigma;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterSigma,daughterPi);
    // add decay table
    decayTable->Insert(mode);
  }

  return decayTable;
}


G4DecayTable*  G4ExcitedSigmaConstructor::AddSigmaStarPiMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;


  G4String daughterSigma;
  G4String daughterPi;
  G4double r;

  // ------------ Sigma+ pi - ------------ 
  // determine daughters
  if (iIso3== +2) {
    r = 0.;
  } else if (iIso3== 0) {
    daughterSigma = "sigma(1385)+";
    r = br/2.;
  } else if (iIso3== -2) {
    daughterSigma = "sigma(1385)0";
    r = br/2.;
  }
  if (!fAnti) {
    daughterPi = "pi-";
  } else {
    daughterPi = "pi+";
  }
  if (fAnti) daughterSigma = "anti_" + daughterSigma;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterSigma,daughterPi);
    // add decay table
    decayTable->Insert(mode);
  }
  // ------------ Sigma0 Pi0 ------------ 
  // determine daughters
  if (iIso3== +2) {
    daughterSigma = "sigma(1385)+";
    r = br/2.;
  } else if (iIso3== 0) {
    r = 0.;
  } else if (iIso3== -2) {
    daughterSigma = "sigma(1385)-";
    r = br/2.;
  }
  daughterPi = "pi0";
  if (fAnti) daughterSigma = "anti_" + daughterSigma;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterSigma,daughterPi);

    // add decay table
    decayTable->Insert(mode);
  }

  // ------------ Sigma- pi + ------------ 
  // determine daughters
  if (iIso3== +2) {
    daughterSigma = "sigma(1385)0";
    r = br/2.;
  } else if (iIso3== 0) {
    daughterSigma = "sigma(1385)-";
    r = br/2.;
  } else if (iIso3== -2) {
    r = 0.;
  }
  if (!fAnti) {
    daughterPi = "pi+";
  } else {
    daughterPi = "pi-";
  }
  if (fAnti) daughterSigma = "anti_" + daughterSigma;
  if (r>0.) {
    // create decay channel  [parent    BR     #daughters]

    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                        daughterSigma,daughterPi);
    // add decay table
    decayTable->Insert(mode);
  }

  return decayTable;
}

G4DecayTable*  G4ExcitedSigmaConstructor::AddLambdaPiMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterLambda = "lambda";
  G4String daughterPi;

  // determine daughters
  if (iIso3== +2) {
    if (!fAnti) {
      daughterPi = "pi+";
    } else {
      daughterPi = "pi-";
    }
  } else if (iIso3== 0) {
    daughterPi = "pi0";
  } else if (iIso3== -2) {
    if (!fAnti) {
      daughterPi = "pi-";
    } else {
      daughterPi = "pi+";
    }
  }
  if (fAnti) daughterLambda = "anti_" + daughterLambda;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           daughterLambda, daughterPi);
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedSigmaConstructor::AddLambdaStarPiMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterLambda = "lambda(1405)";
  G4String daughterPi;

  // determine daughters
  if (iIso3== +2) {
    if (!fAnti) {
      daughterPi = "pi+";
    } else {
      daughterPi = "pi-";
    }
  } else if (iIso3== 0) {
    daughterPi = "pi0";
  } else if (iIso3== -2) {
    if (!fAnti) {
      daughterPi = "pi-";
    } else {
      daughterPi = "pi+";
    }
  }

  if (fAnti) daughterLambda = "anti_" + daughterLambda;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           daughterLambda,daughterPi);
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

const G4String G4ExcitedSigmaConstructor::name[] = {
  "sigma(1385)","sigma(1660)","sigma(1670)","sigma(1750)","sigma(1775)",
  "sigma(1915)","sigma(1940)","sigma(2030)"
};

const G4double G4ExcitedSigmaConstructor::mass[] = {
  1.384*GeV, 1.660*GeV, 1.670*GeV, 1.750*GeV,  1.775*GeV, 
  1.915*GeV, 1.940*GeV, 2.030*GeV
};

const G4double G4ExcitedSigmaConstructor::width[] = {
   36.0*MeV, 100.0*MeV,  60.0*MeV,  90.0*MeV, 120.0*MeV,
  120.0*MeV, 330.0*MeV, 180.0*MeV
};

const G4int G4ExcitedSigmaConstructor::iSpin[] = {
    3,   1,   3,   1,   5,
    5,   3,   7
};

const G4int G4ExcitedSigmaConstructor::iParity[] = {
  +1,  +1,   -1,  -1,  -1,
  +1,  -1,   +1
};


const G4int G4ExcitedSigmaConstructor::encodingOffset[] = {
      0, 10000, 10000, 20000,     0,
  10000, 20000,     0
};

const G4double G4ExcitedSigmaConstructor::bRatio[ G4ExcitedSigmaConstructor::NStates ][ G4ExcitedSigmaConstructor::NumberOfDecayModes] = 
{
   {   0.0,  0.0, 0.12,  0.0, 0.88,   0.0,   0.0,   0.0}, 
   {  0.30,  0.0, 0.35,  0.0, 0.35,   0.0,   0.0,   0.0}, 
   {  0.15,  0.0, 0.70,  0.0, 0.15,   0.0,   0.0,   0.0}, 
   {  0.40,  0.0, 0.05,  0.0,  0.0,  0.55,   0.0,   0.0}, 
   {  0.40,  0.0, 0.04, 0.10, 0.23,   0.0,  0.23,   0.0}, 
   {  0.15,  0.0, 0.40, 0.05, 0.40,   0.0,   0.0,   0.0}, 
   {  0.10, 0.15, 0.15, 0.15, 0.15,   0.0,  0.15,  0.15}, 
   {  0.20, 0.04, 0.10, 0.10, 0.20,   0.0,  0.18,  0.18} 

};











