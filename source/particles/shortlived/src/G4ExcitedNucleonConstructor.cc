// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedNucleonConstructor.cc,v 1.4 1999-12-15 14:51:17 gunter Exp $
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


#include "G4ExcitedNucleonConstructor.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"


G4ExcitedNucleonConstructor::G4ExcitedNucleonConstructor():
    G4ExcitedBaryonConstructor(NStates, NucleonIsoSpin)
{

}

G4ExcitedNucleonConstructor::~G4ExcitedNucleonConstructor()
{
}

G4int G4ExcitedNucleonConstructor::GetEncoding(G4int iIsoSpin3, G4int idxState)
{
  G4int encoding;
  // Delta has exceptinal encoding
  if ((idxState==1)||(idxState==6)||(idxState==8)||(idxState==9)||(idxState==12) )  {
    encoding = GetEncodingOffset(idxState);
    if ((iIsoSpin3==3)||(iIsoSpin3==-3)) {
      // normal encoding 
      encoding += 1000*GetQuarkContents(0, iIsoSpin3);
      encoding +=  100*GetQuarkContents(1, iIsoSpin3);
      encoding +=   10*GetQuarkContents(2, iIsoSpin3);
    } else if (iIsoSpin3== +1){
      // 1st <--> 2nd quark 
      encoding += 1000*GetQuarkContents(0, iIsoSpin3);
      encoding +=   10*GetQuarkContents(1, iIsoSpin3);
      encoding +=  100*GetQuarkContents(2, iIsoSpin3);
    } else if (iIsoSpin3== -1){
      // 1st <--> 0th quark 
      encoding +=  100*GetQuarkContents(0, iIsoSpin3);
      encoding += 1000*GetQuarkContents(1, iIsoSpin3);
      encoding +=   10*GetQuarkContents(2, iIsoSpin3);
    }
    encoding += GetiSpin(idxState) +1;
  } else {
    encoding = G4ExcitedBaryonConstructor::GetEncoding(iIsoSpin3, idxState);
  }
  return encoding;
}

G4DecayTable* G4ExcitedNucleonConstructor::CreateDecayTable(
						 const G4String&  parentName,  
						 G4int iIso3, 
						 G4int iState,
						 G4bool fAnti)
{
  // create decay table
  G4DecayTable* decayTable =  new G4DecayTable();

  G4double br;
  if ( (br=bRatio[iState][NGamma]) >0.0) {
    AddNGammaMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][NPi]) >0.0) {
    AddNPiMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][NEta]) >0.0) {
    AddNEtaMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][NOmega]) >0.0) {
    AddNOmegaMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][NRho]) >0.0) {
    AddNRhoMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][N2Pi]) >0.0) {
    AddN2PiMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][DeltaPi]) >0.0) {
    AddDeltaPiMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][NStarPi]) >0.0) {
    AddNStarPiMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][LambdaK]) >0.0) {
    AddLambdaKMode( decayTable, parentName, br, iIso3, fAnti);
  }

  return  decayTable;
}

G4DecayTable*  G4ExcitedNucleonConstructor::AddNGammaMode( 
				   G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  // 
  G4String daughterN;
  if (iIso3 == +1) { 
    daughterN = "proton";  
  } else {
    daughterN = "neutron"; 
  } 
  if (fAnti) daughterN = "anti_" + daughterN;

  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           daughterN,"gamma");
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedNucleonConstructor::AddNPiMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterN;
  G4String daughterPi;

  // ------------ N pi0 ------------ 
  // determine daughters
  if (iIso3 == +1) {
    daughterN  = "proton";
    daughterPi = "pi0";
  } else {
    daughterN  = "neutron";
    daughterPi = "pi0";
  }
  if (fAnti) daughterN = "anti_" + daughterN;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/2.0, 2,
                                           daughterN,daughterPi);
  // add decay table
  decayTable->Insert(mode);

  // -------------N pi +/- -------------- 
  // determine daughters
  if (iIso3 == +1) {
    daughterN  = "neutron";
    if (!fAnti) {
      daughterPi = "pi+";
    } else {
      daughterPi = "pi-";
    }
  } else {
    daughterN  = "proton";
    if (!fAnti) {
      daughterPi = "pi-";
    } else {
      daughterPi = "pi+";
    }
  }
  if (fAnti) daughterN = "anti_" + daughterN;

  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/2.0, 2,
                                           daughterN,daughterPi);
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedNucleonConstructor::AddNEtaMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterN;

  // ------------ N eta------------ 
  // determine daughters
  if (iIso3 == +1) {
    daughterN  = "proton";
  } else {
    daughterN  = "neutron";
  }
  if (fAnti) daughterN = "anti_" + daughterN;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           daughterN, "eta");
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedNucleonConstructor::AddNOmegaMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterN;

  // ------------ N omega------------ 
  // determine daughters
  if (iIso3 == +1) {
    daughterN  = "proton";
  } else {
    daughterN  = "neutron";
  }
  if (fAnti) daughterN = "anti_" + daughterN;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           daughterN, "omega");
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedNucleonConstructor::AddNRhoMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterN;
  G4String daughterRho;

  // ------------ N rho0 ------------ 
  // determine daughters
  if (iIso3 == +1) {
    daughterN  = "proton";
    daughterRho = "rho0";
  } else {
    daughterN  = "neutron";
    daughterRho = "rho0";
  }
  if (fAnti) daughterN = "anti_" + daughterN;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/2.0, 2,
                                           daughterN,daughterRho);
  // add decay table
  decayTable->Insert(mode);

  // -------------N rho+/- -------------- 
  // determine daughters
  if (iIso3 == +1) {
    daughterN  = "neutron";
    if (!fAnti) {
      daughterRho = "rho+";
    } else {
      daughterRho = "rho-";
    }
  } else {
    daughterN  = "proton";
    if (!fAnti) {
      daughterRho = "rho-";
    } else {
      daughterRho = "rho+";
    }
  }
  if (fAnti) daughterN = "anti_" + daughterN;

  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/2.0, 2,
                                           daughterN,daughterRho);
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedNucleonConstructor::AddN2PiMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int iIso3, G4bool fAnti)
{
  // Decay Modes
  //   N* --> N + pi + pi
  //     Only I=0 states are included for 2-pi system 
  
  G4VDecayChannel* mode;

  G4String daughterN;
  G4String daughterPi1;
  G4String daughterPi2;

  // -------------N pi+ pi- -------------- 
  // determine daughters
  if (iIso3 == +1) {
    daughterN  = "proton";
    daughterPi1 = "pi+";
    daughterPi2 = "pi-";
  } else {
    daughterN  = "neutron";
    daughterPi1 = "pi+";
    daughterPi2 = "pi-";
  }
  if (fAnti) daughterN = "anti_" + daughterN;

  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/2.0, 3,
                                           daughterN,daughterPi1,daughterPi2);
  // add decay table
  decayTable->Insert(mode);

  // -------------N pi0 pi0 --------------  
  // determine daughters
  if (iIso3 == +1) {
    daughterN  = "proton";
    daughterPi1 = "pi0";
    daughterPi2 = "pi0";
  } else {
    daughterN  = "neutron";
    daughterPi1 = "pi0";
    daughterPi2 = "pi0";
  }
  if (fAnti) daughterN = "anti_" + daughterN;

  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/2.0, 3,
                                           daughterN,daughterPi1,daughterPi2);
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedNucleonConstructor::AddNStarPiMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterN;
  G4String daughterPi;

  // ------------ N pi0 ------------ 
  // determine daughters
  if (iIso3 == +1) {
    daughterN  = "N(1440)+";
    daughterPi = "pi0";
  } else {
    daughterN  = "N(1440)0";
    daughterPi = "pi0";
  }
  if (fAnti) daughterN = "anti_" + daughterN;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/2.0, 2,
                                           daughterN,daughterPi);
  // add decay table
  decayTable->Insert(mode);

  // -------------N pi +/- -------------- 
  // determine daughters
  if (iIso3 == +1) {
    daughterN  = "N(1440)0";
    if (!fAnti) {
      daughterPi = "pi+";
    } else {
      daughterPi = "pi-";
    }
  } else {
    daughterN  = "N(1440)+";
    if (!fAnti) {
      daughterPi = "pi-";
    } else {
      daughterPi = "pi+";
    }
  }
  if (fAnti) daughterN = "anti_" + daughterN;

  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/2.0, 2,
                                           daughterN,daughterPi);
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedNucleonConstructor::AddDeltaPiMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterDelta;
  G4String daughterPi;
  G4double r;

  // ------------ Delta pi+/- ------------ 
  // determine daughters
  if (iIso3 == +1) {
    daughterDelta  = "delta0";
    if (!fAnti) {
      daughterPi = "pi+";
    } else {
      daughterPi = "pi-";
    }
    r = br/6.0;
  } else {
    daughterDelta  = "delta+";
    if (!fAnti) {
      daughterPi = "pi-";
    } else {
      daughterPi = "pi+";
    }
    r = br/6.0;
  }
  if (fAnti) daughterDelta = "anti_" + daughterDelta;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                           daughterDelta,daughterPi);
  // add decay table
  decayTable->Insert(mode);

  // ------------ Delta pi+/- ------------ 
  // determine daughters
  if (iIso3 == +1) {
    daughterDelta  = "delta++";
    if (!fAnti) {
      daughterPi = "pi-";
    } else {
      daughterPi = "pi+";
    }
    r = br/2.0;
  } else {
    daughterDelta  = "delta-";
    if (!fAnti) {
      daughterPi = "pi+";
    } else {
      daughterPi = "pi-";
    }
    r = br/2.0;
  }
  if (fAnti) daughterDelta = "anti_" + daughterDelta;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                           daughterDelta,daughterPi);
  // add decay table
  decayTable->Insert(mode);

  // ------------ Delta pi0 ------------ 
  // determine daughters
  if (iIso3 == +1) {
    daughterDelta  = "delta+";
    daughterPi = "pi0";
    r = br/3.0;
  } else {
    daughterDelta  = "delta0";
    daughterPi = "pi0";
    r = br/3.0;
  }
  if (fAnti) daughterDelta = "anti_" + daughterDelta;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                           daughterDelta,daughterPi);
  // add decay table
  decayTable->Insert(mode);


  return decayTable;
}

G4DecayTable*  G4ExcitedNucleonConstructor::AddLambdaKMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String lambda = "lambda";
  G4String daughterK;

  // ------------ N pi0 ------------ 
  // determine daughters
  if (iIso3 == +1) {
    if (!fAnti) {
      daughterK = "kaon+";
    } else {
      daughterK = "kaon-";
    }
  } else {
    if (!fAnti) {
      daughterK = "kaon0";
    } else {
      daughterK = "anti_kaon0";
    }
  }
  if (fAnti) lambda = "anti_" + lambda;
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           lambda, daughterK);
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}



const char* G4ExcitedNucleonConstructor::name[] = {
   "N(1440)", "N(1520)", "N(1535)", "N(1650)", "N(1675)",
   "N(1680)", "N(1700)", "N(1710)", "N(1720)", "N(1900)", 
   "N(1990)", "N(2090)", "N(2190)", "N(2220)", "N(2250)"
};

const G4double G4ExcitedNucleonConstructor::mass[] = {
  1.440*GeV, 1.520*GeV, 1.535*GeV, 1.650*GeV,  1.675*GeV, 
  1.680*GeV, 1.700*GeV, 1.710*GeV, 1.720*GeV,  1.870*GeV, 
  1.990*GeV, 2.040*GeV, 2.190*GeV, 2.220*GeV,  2.250*GeV
};

const G4double G4ExcitedNucleonConstructor::width[] = {
  350.0*MeV, 120.0*MeV, 150.0*MeV, 150.0*MeV, 150.0*MeV,
  130.0*MeV, 100.0*MeV, 100.0*MeV, 150.0*MeV, 500.0*MeV,
  550.0*MeV, 250.0*MeV, 450.0*MeV
};

const G4int G4ExcitedNucleonConstructor::iSpin[] = {
    1,   3,   1,   1,   5,
    5,   3,   1,   3,   3,
    7,   3,   7,   9,   9
};

const G4int G4ExcitedNucleonConstructor::iParity[] = {
  +1,  -1,   -1,  -1,  -1,
  +1,  -1,   +1,  +1,  +1,
  +1,  -1,   +1,  -1,  -1 
};

const G4int G4ExcitedNucleonConstructor::encodingOffset[] = {
   10000,       0,  20000,  30000,       0,
   10000,   20000,  40000,  30000,   40000, 
   10000,   50000,       0,     0,   10000
};

const G4double G4ExcitedNucleonConstructor::bRatio[ G4ExcitedNucleonConstructor::NStates ][ G4ExcitedNucleonConstructor::NumberOfDecayModes] = 
{
   {  0.0, 0.70,  0.0,  0.0,  0.0,  0.05,  0.25,  0.0,  0.0}, 
   {  0.0, 0.60,  0.0,  0.0,  0.0,  0.15,  0.25,  0.0,  0.0}, 
   {0.001, 0.55, 0.35,  0.0,  0.0,  0.05,  0.00, 0.05,  0.0},
   {  0.0, 0.65, 0.05,  0.0,  0.0,  0.05,  0.10, 0.05, 0.10},
   {  0.0, 0.45,  0.0,  0.0,  0.0,  0.00,  0.55,  0.0,  0.0},
   {  0.0, 0.65,  0.0,  0.0,  0.0,  0.20,  0.15,  0.0,  0.0},
   {  0.0, 0.10, 0.05,  0.0, 0.05,  0.45,  0.35,  0.0,  0.0},
   {  0.0, 0.15, 0.20,  0.0, 0.05,  0.20,  0.20, 0.10, 0.10},
   {  0.0, 0.15, 0.00,  0.0, 0.25,  0.45,  0.10, 0.00, 0.05},
   {  0.0, 0.35,  0.0, 0.55, 0.05,  0.00,  0.05,  0.0,  0.0},
   {  0.0, 0.05,  0.0,  0.0, 0.15,  0.25,  0.30, 0.15, 0.10},
   {  0.0, 0.60, 0.05,  0.0, 0.25,  0.05,  0.05,  0.0,  0.0},
   {  0.0, 0.35,  0.0, 0.00, 0.30,  0.15,  0.15, 0.05,  0.0},
   {  0.0, 0.35,  0.0,  0.0, 0.25,  0.20,  0.20,  0.0,  0.0},
   {  0.0, 0.30,  0.0, 0.00, 0.25,  0.20,  0.20, 0.05,  0.0}
};
















