// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedDeltaConstructor.cc,v 1.1 1999-01-07 16:10:41 gunter Exp $
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


#include "G4ExcitedDeltaConstructor.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ShortLivedTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"


G4ExcitedDeltaConstructor::G4ExcitedDeltaConstructor():
    G4ExcitedBaryonConstructor(NStates, DeltaIsoSpin)
{

}

G4ExcitedDeltaConstructor::~G4ExcitedDeltaConstructor()
{
}

G4int G4ExcitedDeltaConstructor::GetEncoding(G4int iIsoSpin3, G4int idxState)
{
  G4int encoding;
  // Delta has exceptinal encoding
  if ((idxState==1)||(idxState==3)||(idxState==4)||(idxState==5)||(idxState==7))  {
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
G4DecayTable* G4ExcitedDeltaConstructor::CreateDecayTable(
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

  if ( (br=bRatio[iState][NRho]) >0.0) {
    AddNRhoMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][DeltaPi]) >0.0) {
    AddDeltaPiMode( decayTable, parentName, br, iIso3, fAnti);
  }

  if ( (br=bRatio[iState][NStarPi]) >0.0) {
    AddNStarPiMode( decayTable, parentName, br, iIso3, fAnti);
  }

  return  decayTable;
}

G4DecayTable*  G4ExcitedDeltaConstructor::AddNGammaMode( 
				   G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  // 
  G4String daughterN;
  if (iIso3 == +1) { 
    daughterN = "proton";  
  } else if (iIso3 == -1) {
    daughterN = "neutron"; 
  } else {
    // can not decay into N+gamma
    return decayTable;
  }

  if (fAnti) daughterN = "anti_" + daughterN;

  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           daughterN,"gamma");
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedDeltaConstructor::AddNPiMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterN;
  G4String daughterPi;
  G4double r;

  // ------------ N pi0 ------------ 
  // determine daughters
  if ((iIso3 == +1)||(iIso3 == -1)) {
    if (iIso3 == +1) {
      daughterN  = "proton";
      daughterPi = "pi0";
      r = br*2./3.;
    } else  if (iIso3 == -1) {
      daughterN  = "neutron";
      daughterPi = "pi0";
      r = br/3.;
    } 
    if (fAnti) daughterN = "anti_" + daughterN;
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					daughterN,daughterPi);
    // add decay table
    decayTable->Insert(mode);
  }

  // -------------N pi +/- -------------- 
  // determine daughters
  if (iIso3 == +3) {
    daughterN  = "proton";
    if (!fAnti) {
      daughterPi = "pi+";
    } else {
      daughterPi = "pi-";
    }
    r = br;
  } else  if (iIso3 == +1) {
    daughterN  = "neutron";
    if (!fAnti) {
      daughterPi = "pi+";
    } else {
      daughterPi = "pi-";
    }
    r = br/3.;
  } else if (iIso3 == -1) {
    daughterN  = "proton";
    if (!fAnti) {
      daughterPi = "pi-";
    } else {
      daughterPi = "pi+";
    }
    r = br*2./3.;
  } else if (iIso3 == -3) {
    daughterN  = "neutron";
    if (!fAnti) {
      daughterPi = "pi-";
    } else {
      daughterPi = "pi+";
    }
    r = br;
  }
  if (fAnti) daughterN = "anti_" + daughterN;

  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                           daughterN,daughterPi);
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}


G4DecayTable*  G4ExcitedDeltaConstructor::AddNRhoMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterN;
  G4String daughterRho;
  G4double r;

  // ------------ N Rho0 ------------ 
  // determine daughters
  if ((iIso3 == +1)||(iIso3 == -1)) {
    if (iIso3 == +1) {
      daughterN  = "proton";
      daughterRho = "rho0";
      r = br*2./3.;
    } else  if (iIso3 == -1) {
      daughterN  = "neutron";
      daughterRho = "rho0";
      r = br/3.;
    } 
    if (fAnti) daughterN = "anti_" + daughterN;
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					daughterN,daughterRho);
    // add decay table
    decayTable->Insert(mode);
  }

  // -------------N Rho +/- -------------- 
  // determine daughters
  if (iIso3 == +3) {
    daughterN  = "proton";
    if (!fAnti) {
      daughterRho = "rho+";
    } else {
      daughterRho = "rho-";
    }
    r = br;
  } else  if (iIso3 == +1) {
    daughterN  = "neutron";
    if (!fAnti) {
      daughterRho = "rho+";
    } else {
      daughterRho = "rho-";
    }
    r = br/3.;
  } else if (iIso3 == -1) {
    daughterN  = "proton";
    if (!fAnti) {
      daughterRho = "rho-";
    } else {
      daughterRho = "rho+";
    }
    r = br*2./3.;
  } else if (iIso3 == -3) {
    daughterN  = "neutron";
    if (!fAnti) {
      daughterRho = "rho-";
    } else {
      daughterRho = "rho+";
    }
    r = br;
  }
  if (fAnti) daughterN = "anti_" + daughterN;

  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                           daughterN,daughterRho);
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedDeltaConstructor::AddNStarPiMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterN;
  G4String daughterPi;
  G4double r;

  // ------------ N pi0 ------------ 
  // determine daughters
  if ((iIso3 == +1)||(iIso3 == -1)) {
    if (iIso3 == +1) {
      daughterN  = "N(1440)+";
      daughterPi = "pi0";
      r = br*2./3.;
    } else  if (iIso3 == -1) {
      daughterN  = "N(1440)0";
      daughterPi = "pi0";
      r = br/3.;
    } 
    if (fAnti) daughterN = "anti_" + daughterN;
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					daughterN,daughterPi);
    // add decay table
    decayTable->Insert(mode);
  }

  // -------------N pi +/- -------------- 
  // determine daughters
  if (iIso3 == +3) {
    daughterN  = "N(1440)+";
    if (!fAnti) {
      daughterPi = "pi+";
    } else {
      daughterPi = "pi-";
    }
    r = br;
  } else  if (iIso3 == +1) {
    daughterN  = "N(1440)0";
    if (!fAnti) {
      daughterPi = "pi+";
    } else {
      daughterPi = "pi-";
    }
    r = br/3.;
  } else if (iIso3 == -1) {
    daughterN  = "N(1440)+";
    if (!fAnti) {
      daughterPi = "pi-";
    } else {
      daughterPi = "pi+";
    }
    r = br*2./3.;
  } else if (iIso3 == -3) {
    daughterN  = "N(1440)0";
    if (!fAnti) {
      daughterPi = "pi-";
    } else {
      daughterPi = "pi+";
    }
    r = br;
  }
  if (fAnti) daughterN = "anti_" + daughterN;

  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
                                           daughterN,daughterPi);
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedDeltaConstructor::AddDeltaPiMode( 
				    G4DecayTable* decayTable, const G4String& nameParent,
				    G4double br, G4int iIso3, G4bool fAnti)
{
  G4VDecayChannel* mode;

  G4String daughterDelta;
  G4String daughterPi;
  G4double r;

  // ------------ Delta pi +------------ 
  // determine daughters
  if (iIso3 == +3) {
    daughterDelta  = "delta+";
    r = br*0.4;
  } else  if (iIso3 == +1) {
    daughterDelta  = "delta0";
    r = br*8./15.0;
  } else  if (iIso3 == -1) {
    daughterDelta  = "delta-";
    r = br*6./15.;
  } else {
     r  = 0.;
  }
  if (!fAnti) {
    daughterPi = "pi+";
  } else {
    daughterPi = "pi-";
  }
   if (fAnti) daughterDelta = "anti_" + daughterDelta;
  if (r>0.0) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					daughterDelta,daughterPi);
    // add decay table
    decayTable->Insert(mode);
  }

  // ------------ Delta pi0 ------------ 
  // determine daughters
  if (iIso3 == +3) {
    daughterDelta  = "delta++";
    r = br*0.6;
  } else  if (iIso3 == +1) {
    daughterDelta  = "delta+";
    r = br*1./15.0;
  } else  if (iIso3 == -1) {
    daughterDelta  = "delta0";
    r = br*1./15.;
  } else {
    daughterDelta  = "delta-";
    r = br*0.6;
  }
  daughterPi = "pi0";
  if (fAnti) daughterDelta = "anti_" + daughterDelta;
  
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
				      daughterDelta,daughterPi);
  // add decay table
  decayTable->Insert(mode);

  // ------------ Delta pi - ------------- 
  // determine daughters
  if (iIso3 == +3) {
    r= 0.;
  } else  if (iIso3 == +1) {
    daughterDelta  = "delta++";
    r = br*6./15.0;
  } else  if (iIso3 == -1) {
    daughterDelta  = "delta+";
    r = br*8./15.;
  } else {
    daughterDelta  = "delta0";
    r = br*0.4;
  }
  if (!fAnti) {
    daughterPi = "pi-";
  } else {
    daughterPi = "pi+";
  }
   if (fAnti) daughterDelta = "anti_" + daughterDelta;
  if (r>0.0) {
    // create decay channel  [parent    BR     #daughters]
    mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					daughterDelta,daughterPi);
    // add decay table
    decayTable->Insert(mode);
  }

  return decayTable;
}

const G4String G4ExcitedDeltaConstructor::name[] = 
{
  "delta(1600)", "delta(1620)", "delta(1700)", "delta(1900)", "delta(1905)",
  "delta(1910)", "delta(1920)", "delta(1930)", "delta(1950)"
};

const G4double G4ExcitedDeltaConstructor::mass[] = 
{
  1.700*GeV, 1.675*GeV, 1.750*GeV, 1.850*GeV,  1.880*GeV, 
  1.900*GeV, 1.920*GeV, 1.930*GeV, 1.950*GeV
};

const G4double G4ExcitedDeltaConstructor::width[] = {
  200.0*MeV, 180.0*MeV, 300.0*MeV, 240.0*MeV, 280.0*MeV,
  250.0*MeV, 105.0*MeV, 250.0*MeV, 250.0*MeV
};

const G4int G4ExcitedDeltaConstructor::iSpin[] = 
{
    3,   1,   3,   1,   5,
    1,   3,   5,   7
};

const G4int G4ExcitedDeltaConstructor::iParity[] = {
  +1,  -1,   -1,  -1,  +1,
  +1,  +1,   -1,  +1
};

const G4int G4ExcitedDeltaConstructor::encodingOffset[] = {
  30000,     0, 10000, 10000,     0,
  20000, 20000, 10000,     0
};

const G4double G4ExcitedDeltaConstructor::bRatio[ G4ExcitedDeltaConstructor::NStates ][ G4ExcitedDeltaConstructor::NumberOfDecayModes] = 
{
//   NGamma    Npi     NRho DeltaPi    N*Pi
   {    0.0,   0.15,    0.0,   0.55,   0.30 },
   {    0.0,   0.25,    0.0,   0.60,   0.15 },
   {    0.0,   0.20,   0.10,   0.55,   0.15 },
   {    0.0,   0.30,   0.15,   0.30,   0.25 },
   {    0.0,   0.20,   0.60,   0.10,   0.10 },
   {    0.0,   0.35,   0.40,   0.15,   0.10 },
   {    0.0,   0.15,   0.30,   0.30,   0.25 },
   {    0.0,   0.20,   0.25,   0.25,   0.30 },
   {   0.01,   0.44,   0.15,   0.20,   0.20 }
};
















