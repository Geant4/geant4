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
// $Id: G4ExcitedMesonConstructor.cc 102905 2017-03-02 09:50:56Z gcosmo $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//      History: first implementation, based on object model of
//      10 oct 1998  H.Kurashige
// ---------------------------------------------------------------
//    01 Oct. 02 Fixed PDG codes for a0(1450), f0(1370), k0_star(1430)

#include "G4ExcitedMesonConstructor.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayTable.hh"


G4ExcitedMesonConstructor::G4ExcitedMesonConstructor(G4int ,
						       G4int )
    :   type("meson"), leptonNumber(0), baryonNumber(0)
{
}

G4ExcitedMesonConstructor::~G4ExcitedMesonConstructor()
{
}

void G4ExcitedMesonConstructor::Construct(G4int idx)
{
  G4int iType;
  if (idx < 0 ) {
    for (G4int state=0; state< NMultiplets; state +=1) {
      for (iType = 0; iType < NMesonTypes ; iType++) 
	ConstructMesons(state, iType);
    }
  } else if (idx < NMultiplets ) {
      for (iType = 0; iType < NMesonTypes ; iType++) 
	ConstructMesons(idx, iType);
  } else {
#ifdef G4VERBOSE
    if (G4ParticleTable::GetParticleTable()->GetVerboseLevel()>1) {
      G4cerr << "G4ExcitedMesonConstructor::Construct()";
      G4cerr << "   illegal index os state = " << idx << G4endl;
    }
#endif
  }
}

G4bool G4ExcitedMesonConstructor::Exist(G4int idxState, G4int idxType)
{
  G4bool value = true;
  if ( idxType == TEtaPrime ) {
    if (idxState==N13P0) value = false;      
    if (idxState==N13D1) value = false;      
  } else if ( idxType == TPi ) {
    if (idxState==N23P2) value = false;      
  }
  return value;
}

#include "G4ExcitedMesons.hh"

void G4ExcitedMesonConstructor::ConstructMesons(G4int iState, G4int iType)
{
  if (!Exist(iState, iType) ) return;

  //    Construct Resonace particles as dynamic object
  //    Arguments for constructor are as follows
  //               name             mass          width        
  //             charge           2*spin           
  //             parity    C-conjugation
  //          2*Isospin       2*Isospin3       
  //          G-parity
  //               type    lepton number  Baryon number   
  //        PDG encoding
  //             stable         lifetime    decay table 
  
  
  G4String aName;
  G4ExcitedMesons* particle;
  
  for ( G4int iIso3=(-1)*iIsoSpin[iType]; iIso3<=iIsoSpin[iType]; iIso3+=2) {
    aName= GetName(iIso3, iState, iType);
    G4double fmass =  mass[iState][iType];
    G4double fwidth = width[iState][iType];
	if ( (iType== TK) || (iType==TAntiK ) ) {
	  if ( GetCharge(iIso3,iType) == 0.0) {
		fmass  += massKdiff[iState];
	    fwidth += widthKdiff[iState];
      }
   }
    particle = new G4ExcitedMesons(            
             aName,   fmass,         fwidth,    
		 GetCharge(iIso3,iType),                 iSpin[iState],
			iParity[iState],    iChargeConjugation[iState],
                        iIsoSpin[iType],                         iIso3,   
                iGParity[iState][iType],
	                           type,  leptonNumber,   baryonNumber, 
     GetEncoding(iIso3, iState, iType),
	                          false,           0.0,            NULL
				    );

    if ( (iType==TEta) || (iType==TEtaPrime) || ((iType==TPi)&&(iIso3==0)) ) {
    // set same encoding for AntiParticle
      particle->SetAntiPDGEncoding(GetEncoding(iIso3, iState, iType));
    }
    particle->SetMultipletName(name[iState][iType]);
    particle->SetDecayTable(CreateDecayTable( aName, iIso3, iState, iType));
  }
}


G4int  G4ExcitedMesonConstructor::GetQuarkContents(G4int iQ, 
						   G4int iIso3,
						   G4int iType)
{
  // Quark contents

  G4int quark=0;
  if (iType == TPi) {
    if ( iIso3 == 2 ){
      if ( iQ == 0 ){ quark = 2; }
      else          { quark = 1; }
    } else if ( iIso3 == 0 ){
      quark = 1;
    }  else if ( iIso3 == -2 ){
      if ( iQ == 0 ){ quark = 1; }
      else          { quark = 2; }
    }
  } else if (iType == TEta) {
    quark = 2;

  } else if (iType == TEtaPrime) {
    quark = 3;

  } else if (iType == TAntiK) {
    if ( iIso3 == 1 ){
      if ( iQ == 0 ){ quark = 3; }
      else          { quark = 1; }
    } else if ( iIso3 == -1 ){
      if ( iQ == 0 ){ quark = 3; }
      else          { quark = 2; }
    }
 
  } else if (iType == TK) {
    if ( iIso3 == 1 ){
      if ( iQ == 0 ){ quark = 2; }
      else          { quark = 3; }
    } else if ( iIso3 == -1 ){
      if ( iQ == 0 ){ quark = 1; }
      else          { quark = 3; }
    }

  }
  return quark;
}

G4double  G4ExcitedMesonConstructor::GetCharge(G4int iIsoSpin3, G4int idxType )
{
  static const G4double quark_charge[7] = 
  {
    0., -1./3., +2./3., -1./3., +2./3., -1./3., +2./3.
  };
  
  G4double charge =  quark_charge[GetQuarkContents(0, iIsoSpin3, idxType)]*eplus;
  charge -=  quark_charge[GetQuarkContents(1, iIsoSpin3, idxType)]*eplus;
  return charge;
}

G4int     G4ExcitedMesonConstructor::GetEncoding(G4int iIsoSpin3, 
						 G4int idxState, 
						 G4int idxType    )
{
  G4int encoding = encodingOffset[idxState];
  encoding +=  iSpin[idxState] +1;
  G4int iQ = 0;
  G4int iQbar = 1;

  if ( idxType == TPi ) {
    if (iIsoSpin3<0) {
      iQ = 1;
      iQbar = 0; 
    }
  } else if ( idxType == TK ) {
    iQ = 1;
    iQbar = 0; 
  }


  encoding +=  100*GetQuarkContents(iQ, iIsoSpin3, idxType);
  encoding +=   10*GetQuarkContents(iQbar, iIsoSpin3, idxType);
  if ( idxType == TPi ) {
    if (iIsoSpin3<0) {
      encoding *=  -1;
    }
  } else if ( idxType == TAntiK ) {
    encoding *=  -1;
  }

// PDG2005
//
  if (idxState == 9 ) {
    if (idxType == TEta) {
//   f2(1810)  9030225
      encoding = 9030225;
    } else if (idxType == TEtaPrime) {
//   f2(2010)  9060225
      encoding = 9060225;
    }
  }

// PDG2013
  if (idxState == 1 ) {
    if (idxType == TEta) {
//   f0(1370)  30221
      encoding = 30221;
    }
  }
  return encoding;
}

G4DecayTable*  G4ExcitedMesonConstructor::CreateDecayTable(
					    const G4String& parentName,
					    G4int           iIso3,
                                            G4int           iState, 
					    G4int           iType)
{
   // create decay table
  G4DecayTable* decayTable =  new G4DecayTable();
  G4double br;
 
  if ((iType==TK)||(iType==TAntiK)) {
    
    if ( (br=bRatio[iState][iType][MKPi]) >0.0) {
      AddKPiMode( decayTable, parentName, br, iIso3, iType );
    }
    if ( (br=bRatio[iState][iType][MKStarPi]) >0.0) {
      AddKStarPiMode( decayTable, parentName, br, iIso3, iType );
    }
    if ( (br=bRatio[iState][iType][MKRho]) >0.0) {
      AddKRhoMode( decayTable, parentName, br, iIso3, iType );
    }
    if ( (br=bRatio[iState][iType][MKOmega]) >0.0) {
      AddKOmegaMode( decayTable, parentName, br, iIso3, iType );
    }
    if ( (br=bRatio[iState][iType][MKStar2Pi]) >0.0) {
      AddKStar2PiMode( decayTable, parentName, br, iIso3, iType );
    }
    if ( (br=bRatio[iState][iType][MKTwoPi]) >0.0) {
      AddKTwoPiMode( decayTable, parentName, br, iIso3, iType );
    }
    if ( (br=bRatio[iState][iType][MKEta]) >0.0) {
      AddKEtaMode( decayTable, parentName, br, iIso3, iType );
    }

  } else {
    if ( (br=bRatio[iState][iType][MPiGamma]) >0.0) {
      AddPiGammaMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
    if ( (br=bRatio[iState][iType][MRhoGamma]) >0.0) {
      AddRhoGammaMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
    if ( (br=bRatio[iState][iType][M2Pi]) >0.0) {
      Add2PiMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
    if ( (br=bRatio[iState][iType][MPiRho]) >0.0) {
      AddPiRhoMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
    if ( (br=bRatio[iState][iType][MPiEta]) >0.0) {
      AddPiEtaMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
    if ( (br=bRatio[iState][iType][M3Pi]) >0.0) {
      Add3PiMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
    if ( (br=bRatio[iState][iType][M4Pi]) >0.0) {
      Add4PiMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
    if ( (br=bRatio[iState][iType][MKKStar]) >0.0) {
      AddKKStarMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
    if ( (br=bRatio[iState][iType][M2PiEta]) >0.0) {
      Add2PiEtaMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
    if ( (br=bRatio[iState][iType][MRhoEta]) >0.0) {
      AddRhoEtaMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
    if ( (br=bRatio[iState][iType][M2PiRho]) >0.0) {
      Add2PiRhoMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
    if ( (br=bRatio[iState][iType][M2PiOmega]) >0.0) {
      Add2PiOmegaMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
    if ( (br=bRatio[iState][iType][M2Eta]) >0.0) {
      Add2EtaMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
    if ( (br=bRatio[iState][iType][M2K]) >0.0) {
      Add2KMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
    if ( (br=bRatio[iState][iType][M2KPi]) >0.0) {
      Add2KPiMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
   if ( (br=bRatio[iState][iType][MPiOmega]) >0.0) {
      AddPiOmegaMode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
   if ( (br=bRatio[iState][iType][MPiF2]) >0.0) {
      AddPiF2Mode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
   if ( (br=bRatio[iState][iType][MPiF0]) >0.0) {
      AddPiF0Mode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
   if ( (br=bRatio[iState][iType][MPiA2]) >0.0) {
      AddPiA2Mode( decayTable, parentName, br, iIso3, iIsoSpin[iType] );
    }
  }

  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::AddKPiMode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iType)
{
  G4VDecayChannel* mode;
  // 
  if (iIso3 == +1) {
    if (iType == TK) {    
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "kaon+","pi0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "kaon0","pi+");
      decayTable->Insert(mode);
    }else if (iType==TAntiK) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "anti_kaon0","pi0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "kaon-","pi+");
      decayTable->Insert(mode);
   }
  } else if (iIso3 == -1) {
    if (iType == TK) {    
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "kaon0","pi0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "kaon+","pi-");
      decayTable->Insert(mode);
 
    }else if (iType==TAntiK) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "kaon-","pi0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "anti_kaon0","pi-");
      decayTable->Insert(mode);
    }
  }

  return decayTable;
}
G4DecayTable*  G4ExcitedMesonConstructor::AddKTwoPiMode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iType)
{
  G4VDecayChannel* mode;
  // 
  if (iIso3 == +1) {
    if (iType == TK) {    
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "k2_star(1430)+","pi0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "k2_star(1430)0","pi+");
      decayTable->Insert(mode);
    }else if (iType==TAntiK) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "anti_k2_star(1430)0","pi0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "k2_star(1430)-","pi+");
      decayTable->Insert(mode);
   }
  } else if (iIso3 == -1) {
    if (iType == TK) {    
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "k2_star(1430)0","pi0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "k2_star(1430)+","pi-");
      decayTable->Insert(mode);
 
    }else if (iType==TAntiK) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "k2_star(1430)-","pi0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "anti_k2_star(1430)0","pi-");
      decayTable->Insert(mode);
    }
  }

  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::AddKOmegaMode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iType)
{
  G4VDecayChannel* mode;
  // 
  if (iIso3 == +1) {
    if (iType == TK) {    
      mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
					   "kaon+","omega");
      decayTable->Insert(mode);
    }else if (iType==TAntiK) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
					   "anti_kaon0","omega");
      decayTable->Insert(mode);
   } 
  } else if (iIso3 == -1) {
    if (iType == TK) {    
      mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
					   "kaon0","omega");
      decayTable->Insert(mode);
    }else if (iType==TAntiK) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
					   "kaon-","omega");
      decayTable->Insert(mode);
    }
  }

  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::AddKEtaMode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iType)
{
  G4VDecayChannel* mode;
  // 
  if (iIso3 == +1) {
    if (iType == TK) {    
      mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
					   "kaon+","eta");
      decayTable->Insert(mode);
    }else if (iType==TAntiK) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
					   "anti_kaon0","eta");
      decayTable->Insert(mode);
   }   
  } else if (iIso3 == -1) {
    if (iType == TK) {    
      mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
					   "kaon0","eta");
      decayTable->Insert(mode);
    }else if (iType==TAntiK) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
					   "kaon-","eta");
      decayTable->Insert(mode);
    }
  }

  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::AddKRhoMode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iType)
{
  G4VDecayChannel* mode;
  // 
  if (iIso3 == +1) {
    if (iType == TK) {    
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "kaon+","rho0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "kaon0","rho+");
      decayTable->Insert(mode);
    }else if (iType==TAntiK) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "anti_kaon0","rho0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "kaon-","rho+");
      decayTable->Insert(mode);
   }
  } else if (iIso3 == -1) {
    if (iType == TK) {    
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "kaon0","rho0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "kaon+","rho-");
      decayTable->Insert(mode);
 
    }else if (iType==TAntiK) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "kaon-","rho0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "anti_kaon0","rho-");
      decayTable->Insert(mode);
    }
  }

  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::AddKStarPiMode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iType)
{
  G4VDecayChannel* mode;
  // 
  if (iIso3 == +1) {
    if (iType == TK) {    
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "k_star+","pi0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "k_star0","pi+");
      decayTable->Insert(mode);
    }else if (iType==TAntiK) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "anti_k_star0","pi0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "k_star-","pi+");
      decayTable->Insert(mode);
   }
  } else if (iIso3 == -1) {
    if (iType == TK) {    
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "k_star0","pi0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "k_star+","pi-");
      decayTable->Insert(mode);
 
    }else if (iType==TAntiK) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 2,
					   "k_star-","pi0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 2,
					   "anti_k_star0","pi-");
      decayTable->Insert(mode);
    }
  }

  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::AddKStar2PiMode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iType)
{
  // K* --> K pipi(I=1)
  G4VDecayChannel* mode;
  // 
  if (iIso3 == +1) {
    if (iType == TK) {    
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 3,
					   "k_star+","pi+","pi-");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 3,
					   "k_star0","pi+","pi0");
      decayTable->Insert(mode);
    }else if (iType==TAntiK) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 3,
					   "anti_k_star0","pi+","pi-");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 3,
					   "k_star-","pi+","pi0");
      decayTable->Insert(mode);
   }
  } else if (iIso3 == -1) {
    if (iType == TK) {    
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 3,
					   "k_star0","pi+","pi-");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 3,
					   "k_star+","pi-","pi0");
      decayTable->Insert(mode);
 
    }else if (iType==TAntiK) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 3,
					   "k_star-","pi+","pi-");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 3,
					   "anti_k_star0","pi-","pi0");
      decayTable->Insert(mode);
    }
  }

  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::AddPiGammaMode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iIso)
{
  if ((iIso!=2)&&(iIso!=0)) return decayTable;

  G4VDecayChannel* mode;
  // 
  G4String daughter;
  if (iIso3 == +2) { 
    daughter = "pi+";  
  } else if (iIso3 == 0) { 
    daughter = "pi0";  
  } else if (iIso3 ==-2) { 
    daughter = "pi-";  
  } else {
    return decayTable; 
  } 
    // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           daughter,"gamma");
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::AddPiOmegaMode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iIso)
{
  if ((iIso!=2)&&(iIso!=0)) return decayTable;

  G4VDecayChannel* mode;
  // 
  G4String daughter;
  if (iIso3 == +2) { 
    daughter = "pi+";  
  } else if (iIso3 == 0) { 
    daughter = "pi0";  
  } else if (iIso3 ==-2) { 
    daughter = "pi-";  
  } else {
    return decayTable; 
  } 
    // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           daughter,"omega");
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::AddRhoGammaMode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iIso)
{
  if ((iIso!=2)&&(iIso!=0)) return decayTable;

  G4VDecayChannel* mode;
  // 
  G4String daughter;
  if (iIso3 == +2) { 
    daughter = "rho+";  
  } else if (iIso3 == 0) { 
    daughter = "rho0";  
  } else if (iIso3 ==-2) { 
    daughter = "rho-";  
  } else {
    return decayTable; 
  } 
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           daughter,"gamma");
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::AddPiEtaMode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iIso)
{
  if ((iIso!=2)&&(iIso!=0)) return decayTable;

  G4VDecayChannel* mode;
  // 
  G4String daughter;
  if (iIso3 == +2) { 
    daughter = "pi+";  
  } else if (iIso3 == 0) { 
    daughter = "pi0";  
  } else if (iIso3 ==-2) { 
    daughter = "pi-";  
  } else {
    return decayTable; 
  } 
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           daughter,"eta");
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::AddRhoEtaMode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iIso)
{
  if ((iIso!=2)&&(iIso!=0)) return decayTable;

  G4VDecayChannel* mode;
  // 
  G4String daughter;
  if (iIso3 == +2) { 
    daughter = "rho+";  
  } else if (iIso3 == 0) { 
    daughter = "rho0";  
  } else if (iIso3 ==-2) { 
    daughter = "rho-";  
  } else {
    return decayTable; 
  } 
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           daughter,"eta");
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::AddPiF2Mode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iIso)
{
  if ((iIso!=2)&&(iIso!=0)) return decayTable;

  G4VDecayChannel* mode;
  // 
  G4String daughter;
  if (iIso3 == +2) { 
    daughter = "pi+";  
  } else if (iIso3 == 0) { 
    daughter = "pi0";  
  } else if (iIso3 ==-2) { 
    daughter = "pi-";  
  } else {
    return decayTable; 
  } 
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           daughter,"f2(1270)");
  // add decay table
  decayTable->Insert(mode);

  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::AddPiF0Mode( 
                                   G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iIso)
{
  if ((iIso!=2)&&(iIso!=0)) return decayTable;

  G4VDecayChannel* mode;
  // 
  G4String daughter;
  if (iIso3 == +2) { 
    daughter = "pi+";  
  } else if (iIso3 == 0) { 
    daughter = "pi0";  
  } else if (iIso3 ==-2) { 
    daughter = "pi-";  
  } else {
    return decayTable; 
  } 
  // create decay channel  [parent    BR     #daughters]
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
                                           daughter,"f0(1370)");
  // add decay table
  decayTable->Insert(mode);
  return decayTable;
}


G4DecayTable*  G4ExcitedMesonConstructor::Add2PiMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iIso)
{
  G4VDecayChannel* mode;

  G4String daughterPi1;
  G4String daughterPi2;
  G4double r; 

  // I = 0 states
  if (iIso==0) {
    if (iIso3==0) {
     // pi+ + pi-
      daughterPi1 = "pi+";
      daughterPi2 = "pi-";
      r = br*2./3.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi1,daughterPi2);
      decayTable->Insert(mode);

      // pi0 + pi0
      daughterPi1 = "pi0";
      daughterPi2 = "pi0";
      r = br*1./3.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi1,daughterPi2);
      decayTable->Insert(mode);
    } 
  } else if (iIso==2) {
    if (iIso3==+2) {
      // pi+ + pi0
      daughterPi1 = "pi+";
      daughterPi2 = "pi0";
      r = br;
            mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi1,daughterPi2);
      // add decay table
      decayTable->Insert(mode);
    } else if (iIso3==0) {
       // pi+ + pi-
      daughterPi1 = "pi+";
      daughterPi2 = "pi-";
      r = br;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi1,daughterPi2);
      decayTable->Insert(mode);
     } else if (iIso3==-2) {
       // pi- + pi0
      daughterPi1 = "pi-";
      daughterPi2 = "pi0";
      r = br;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi1,daughterPi2);
      decayTable->Insert(mode);
    }
  }
  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::AddPiRhoMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iIso)
{
  G4VDecayChannel* mode;

  G4String daughterPi;
  G4String daughterRho;
  G4double r; 

  // I = 0 states
  if (iIso==0) {
    if (iIso3==0) {
      // pi+ + rho-
      daughterPi = "pi+";
      daughterRho = "rho-";
      r = br/3.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterRho);
      decayTable->Insert(mode);
      
      // pi0 + rho0
      daughterPi = "pi0";
      daughterRho = "rho0";
      r = br*1./3.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterRho);
      decayTable->Insert(mode);
      
      // pi- + rho+
      daughterPi = "pi-";
      daughterRho = "rho+";
      r = br*1./3.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterRho);
      decayTable->Insert(mode);
    } 
  } else if (iIso==2) {
    if (iIso3==+2) {
      // pi+ + rho0
      daughterPi = "pi+";
      daughterRho = "rho0";
      r = br/2.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterRho);
      decayTable->Insert(mode);
 
      // pi0 + rho+
      daughterPi = "pi0";
      daughterRho = "rho+";
      r = br/2.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterRho);
      decayTable->Insert(mode);
    } else if (iIso3==0) {
       // pi+ + rho-
      daughterPi = "pi+";
      daughterRho = "rho-";
      r = br/2.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterRho);
      decayTable->Insert(mode);

       // pi- + rho+
      daughterPi = "pi-";
      daughterRho = "rho+";
      r = br/2.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterRho);
      decayTable->Insert(mode);
     } else if (iIso3==-2) {
       // pi- + rho0
      daughterPi = "pi-";
      daughterRho = "rho0";
      r = br/2.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterRho);
      decayTable->Insert(mode);
    
      // pi0 + rho-
      daughterPi = "pi0";
      daughterRho = "rho-";
      r = br/2.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterRho);
      decayTable->Insert(mode);
    }
  }
  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::AddPiA2Mode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iIso)
{
  G4VDecayChannel* mode;

  G4String daughterPi;
  G4String daughterA2;
  G4double r; 

  // I = 0 states
  if (iIso==0) {
    if (iIso3==0) {
      // pi+ + a2(1320)-
      daughterPi = "pi+";
      daughterA2 = "a2(1320)-";
      r = br/3.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterA2);
      decayTable->Insert(mode);
      
      // pi0 + a2(1320)0
      daughterPi = "pi0";
      daughterA2 = "a2(1320)0";
      r = br*1./3.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterA2);
      decayTable->Insert(mode);
      
      // pi- + a2(1320)+
      daughterPi = "pi-";
      daughterA2 = "a2(1320)+";
      r = br*1./3.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterA2);
      decayTable->Insert(mode);
    } 
  } else if (iIso==2) {
    if (iIso3==+2) {
      // pi+ + a2(1320)0
      daughterPi = "pi+";
      daughterA2 = "a2(1320)0";
      r = br/2.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterA2);
      decayTable->Insert(mode);
 
      // pi0 + a2(1320)+
      daughterPi = "pi0";
      daughterA2 = "a2(1320)+";
      r = br/2.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterA2);
      decayTable->Insert(mode);
    } else if (iIso3==0) {
       // pi+ + a2(1320)-
      daughterPi = "pi+";
      daughterA2 = "a2(1320)-";
      r = br/2.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterA2);
      decayTable->Insert(mode);

       // pi- + a2(1320)+
      daughterPi = "pi-";
      daughterA2 = "a2(1320)+";
      r = br/2.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterA2);
      decayTable->Insert(mode);
     } else if (iIso3==-2) {
       // pi- + a2(1320)0
      daughterPi = "pi-";
      daughterA2 = "a2(1320)0";
      r = br/2.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterA2);
      decayTable->Insert(mode);
    
      // pi0 + a2(1320)-
      daughterPi = "pi0";
      daughterA2 = "a2(1320)-";
      r = br/2.;
      mode = new G4PhaseSpaceDecayChannel(nameParent, r, 2,
					  daughterPi,daughterA2);
      decayTable->Insert(mode);
    }
  }
  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::Add3PiMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iIso)
{
  G4VDecayChannel* mode;

  // I =0 state
  // This mode is X(I=0,J=1) --> pi+,pi-,pi0 mode
  if (iIso==0) {
    // pi+ + pi-
    mode = new G4PhaseSpaceDecayChannel(nameParent, br, 3,
					"pi+","pi-","pi0");
    decayTable->Insert(mode);
  } else if (iIso==2) {
  // This mode is X(I=1) --> pi + pipi(I=0) mode
    if (iIso3==+2) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 3,
					"pi+","pi0","pi0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 3,
					"pi+","pi+","pi-");
      decayTable->Insert(mode);
    } else if (iIso3==0) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 3,
					"pi0","pi0","pi0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 3,
					"pi0","pi+","pi-");
      decayTable->Insert(mode);
    } else if (iIso3==-2) {
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 3,
					"pi-","pi0","pi0");
      decayTable->Insert(mode);
      mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 3,
					"pi-","pi+","pi-");
      decayTable->Insert(mode);
    }
  } 
  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::Add4PiMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int )
{
  G4VDecayChannel* mode;

  if (iIso3==0) {
    // 2pi+ + 2pi-
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/2., 4,
					"pi+","pi-","pi+","pi-");
    decayTable->Insert(mode);
    // pi+ + pi- + 2pi0 
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/2., 4,
					"pi+","pi-","pi0","pi0");
    decayTable->Insert(mode);
  } else if (iIso3==+2) {
    // pi+ + 3pi0
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 4,
					"pi+","pi0","pi0","pi0");
    decayTable->Insert(mode);
    // 2pi+ + pi- + pi0
    mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 4,
					"pi+","pi+","pi-","pi0");
    decayTable->Insert(mode);
  } else if (iIso3==-2) {
    // pi- + 3pi0
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 4,
					"pi-","pi0","pi0","pi0");
    decayTable->Insert(mode);
    // 2pi- + pi+ + pi0
    mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 4,
					"pi-","pi-","pi+","pi0");
    decayTable->Insert(mode);
  }
  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::Add2PiEtaMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int , G4int iIso)
{
  // f1-->eta + pi + pi mode

  if (iIso!=0) return decayTable;

  G4VDecayChannel* mode;

  // eta pi+ pi-
  mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 3,
				          "eta","pi+","pi-");
  decayTable->Insert(mode);

  // eta pi+ pi-
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 3,
				          "eta","pi0","pi0");
  decayTable->Insert(mode);
  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::Add2EtaMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int , G4int iIso)
{
  if (iIso!=0) return decayTable;

  G4VDecayChannel* mode;

  // eta eta
  mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
				          "eta","eta");
  decayTable->Insert(mode);
  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::Add2PiOmegaMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iIso)
{

  G4VDecayChannel* mode;
  if (iIso==0) {
    // omega pi+ pi-
    mode = new G4PhaseSpaceDecayChannel(nameParent, br*2./3., 3,
				          "omega","pi+","pi-");
    decayTable->Insert(mode);

    // omega pi+ pi-
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 3,
				          "omega","pi0","pi0");
    decayTable->Insert(mode);
  } else if (iIso==2) {
    if (iIso3==+2) {
      // omega pi+ pi0
      mode = new G4PhaseSpaceDecayChannel(nameParent, br, 3,
				          "omega","pi+","pi0");
      decayTable->Insert(mode);
    } else if (iIso3==0) {
      // omega pi+ pi-
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/2., 3,
				          "omega","pi-","pi+");
      decayTable->Insert(mode);
      // omega pi0 pi0
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/2., 3,
				          "omega","pi0","pi0");
      decayTable->Insert(mode);
     } else if (iIso3==-2) {
      // omega pi- pi0
      mode = new G4PhaseSpaceDecayChannel(nameParent, br, 3,
				          "omega","pi-","pi0");
      decayTable->Insert(mode);
     }
  }
  return decayTable;
}
  


G4DecayTable*  G4ExcitedMesonConstructor::Add2PiRhoMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int iIso)
{
  G4VDecayChannel* mode;

  if (iIso==0) {
    // f1 --> rho0 + pi+ pi-
    // rho0 pi+ pi-
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/5., 3,
				          "rho0","pi+","pi-");
    decayTable->Insert(mode);
  } else if (iIso==2) {
    if (iIso3==+2) {
      // rho+ pi0 pi0
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/5., 3,
				          "rho+","pi0","pi0");
      decayTable->Insert(mode);
      // rho+ pi+ pi-
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/5., 3,
				          "rho+","pi+","pi-");
      decayTable->Insert(mode);
      // rho0 pi+ pi0
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/5., 3,
				          "rho0","pi+","pi0");
      decayTable->Insert(mode);
      // rho- pi+ pi+
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/5., 3,
				          "rho-","pi+","pi+");
      decayTable->Insert(mode);
    } else if (iIso3==-2) {
      // rho- pi0 pi0
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/5., 3,
				          "rho-","pi0","pi0");
      decayTable->Insert(mode);
      // rho- pi+ pi-
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/5., 3,
				          "rho-","pi+","pi-");
      decayTable->Insert(mode);
      // rho0 pi- pi0
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/5., 3,
				          "rho0","pi-","pi0");
      decayTable->Insert(mode);
      // rho+ pi- pi-
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/5., 3,
				          "rho+","pi-","pi-");
      decayTable->Insert(mode);
    } else if (iIso3==0) {
      // rho+ pi- pi0
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/5., 3,
				          "rho+","pi-","pi0");
      decayTable->Insert(mode);
      // rho0 pi+ pi-
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/5., 3,
				          "rho0","pi+","pi-");
      decayTable->Insert(mode);
      // rho0 pi0 pi0
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/5., 3,
				          "rho0","pi0","pi0");
      decayTable->Insert(mode);
      // rho- pi+ pi0
      mode = new G4PhaseSpaceDecayChannel(nameParent, br/5., 3,
				          "rho-","pi+","pi-");
      decayTable->Insert(mode);
    }
  }
  return decayTable;
}


G4DecayTable*  G4ExcitedMesonConstructor::AddKKStarMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int )
{
  G4VDecayChannel* mode;

  if (iIso3==0) {
    // X(I=0,J=1)-->K + Anti-K*, Anti_K + K* mode
    // K+ + K*-
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/4., 2,
					"kaon+","k_star-");
    decayTable->Insert(mode);
    
    // K- + K*+
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/4., 2,
					"kaon-","k_star0");
    decayTable->Insert(mode);
    
    // K0 + Anti_K*0
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/4., 2,
					"kaon0","anti_k_star0");
    decayTable->Insert(mode);
    
    // Anti_K0 + K*0
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/4., 2,
					"anti_kaon0","k_star0");
    decayTable->Insert(mode);

  } else if (iIso3==2) {  
     // K+ + Anti_K*0
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/2., 2,
  					"kaon+","anti_k_star0");
    decayTable->Insert(mode);
    
     // K0 + K*+
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/2., 2,
  					"anti_kaon0","k_star+");
    decayTable->Insert(mode);

  } else if (iIso3==-2) {  
     // K- + K*0
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/2., 2,
  					"kaon-","k_star0");
    decayTable->Insert(mode);
    
     // K0 + K*-
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/2., 2,
  					"kaon0","k_star-");
    decayTable->Insert(mode);
    
  }

  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::Add2KMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int iIso3, G4int )
{
  G4VDecayChannel* mode;

  if (iIso3==0) {
    // K+ + K-
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/2., 2,
				          "kaon+","kaon-");
    decayTable->Insert(mode);
   
    // K0 + Anti_K0
    mode = new G4PhaseSpaceDecayChannel(nameParent, br/2., 2,
				        "kaon0","anti_kaon0");
    decayTable->Insert(mode);
  } else if  (iIso3==+2) {
    // K+ + anti_K0
    mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
				          "kaon+","anti_kaon0");
    decayTable->Insert(mode);
  } else if  (iIso3==-2) {
    // K- + K0
    mode = new G4PhaseSpaceDecayChannel(nameParent, br, 2,
				          "kaon-","kaon0");
    decayTable->Insert(mode);
  }   
   
  return decayTable;
}

G4DecayTable*  G4ExcitedMesonConstructor::Add2KPiMode( 
                                    G4DecayTable* decayTable, const G4String& nameParent,
                                    G4double br, G4int , G4int iIso)
{

  // X(I=0)-->KKpi 
  if (iIso!=0) return decayTable; 

  G4VDecayChannel* mode;

  // K+ + K- + pi0
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/6., 3,
				      "kaon+","kaon-","pi0");
  decayTable->Insert(mode);
  
  // K0 + Anti_K0 + pi0
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/6., 3,
				      "kaon0","anti_kaon0","pi0");
  decayTable->Insert(mode);
  
  // K+ + anti_K0 + pi-
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 3,
				          "kaon+","anti_kaon0","pi-");
  decayTable->Insert(mode);
  
  // K- + K0 + pi+
  mode = new G4PhaseSpaceDecayChannel(nameParent, br/3., 3,
				          "kaon-","kaon0","pi+");
  decayTable->Insert(mode);
 
   
  return decayTable;
}

// PDG2005
//   eta(1440)   is renamed to eta(1475)  
//   omega(1600)  is renamed to omega(1650)
//
//


const char* G4ExcitedMesonConstructor::name[G4ExcitedMesonConstructor::NMultiplets ][ G4ExcitedMesonConstructor::NMesonTypes ] =
{
  { "b1(1235)",   "h1(1170)",   "h1(1380)",      "k1(1270)",      "k1(1270)" },
  { "a0(1450)",   "f0(1370)",           "", "k0_star(1430)", "k0_star(1430)" },
  { "a1(1260)",   "f1(1285)",   "f1(1420)",      "k1(1400)",      "k1(1400)" },
  { "a2(1320)",   "f2(1270)","f2_prime(1525)","k2_star(1430)","k2_star(1430)"},
  {"pi2(1670)", "eta2(1645)", "eta2(1870)",      "k2(1770)",      "k2(1770)" },
  {"rho(1700)", "omega(1650)",          "",  "k_star(1680)",  "k_star(1680)" },
  {"rho3(1690)","omega3(1670)","phi3(1850)", "k3_star(1780)", "k3_star(1780)" },
  { "pi(1300)",  "eta(1295)",  "eta(1475)",       "k(1460)",       "k(1460)" },
  {"rho(1450)","omega(1420)",  "phi(1680)",  "k_star(1410)",  "k_star(1410)" },
  {         "",   "f2(1810)",   "f2(2010)", "k2_star(1980)", "k2_star(1980)" }
};

const G4double G4ExcitedMesonConstructor::mass[G4ExcitedMesonConstructor::NMultiplets ][ G4ExcitedMesonConstructor::NMesonTypes ] = 
{
  {  1.2295*GeV, 1.170*GeV, 1.386*GeV, 1.272*GeV,  1.272*GeV },
  {   1.474*GeV, 1.350*GeV,       0.0, 1.430*GeV,  1.430*GeV },
  {   1.230*GeV, 1.282*GeV,1.4264*GeV, 1.403*GeV,  1.403*GeV },
  {  1.3183*GeV,1.2755*GeV, 1.525*GeV,1.4256*GeV, 1.4256*GeV },
  {  1.6722*GeV, 1.617*GeV, 1.842*GeV, 1.773*GeV,  1.773*GeV },
  {   1.720*GeV, 1.670*GeV,       0.0, 1.717*GeV,  1.717*GeV },
  {  1.6888*GeV, 1.667*GeV, 1.854*GeV, 1.776*GeV,  1.776*GeV },
  {   1.300*GeV, 1.294*GeV, 1.476*GeV, 1.460*GeV,  1.460*GeV },
  {   1.465*GeV, 1.425*GeV, 1.680*GeV, 1.414*GeV,  1.414*GeV },
  {         0.0, 1.815*GeV, 2.010*GeV, 1.973*GeV,  1.973*GeV }
};

const G4double  G4ExcitedMesonConstructor::massKdiff[ G4ExcitedMesonConstructor::NMultiplets ] = {
	0.0*MeV,  0.0*MeV, 0.0*MeV, 6.8*MeV, 0.0*MeV, 
    0.0*MeV,  0.0*MeV, 0.0*MeV, 0.0*MeV, 0.0*MeV
};

const G4double  G4ExcitedMesonConstructor::widthKdiff[ G4ExcitedMesonConstructor::NMultiplets ] = {
	0.0*MeV,  0.0*MeV, 0.0*MeV, 10.5*MeV, 0.0*MeV, 
    0.0*MeV,  0.0*MeV, 0.0*MeV, 0.0*MeV, 0.0*MeV
};

const G4double G4ExcitedMesonConstructor::width[G4ExcitedMesonConstructor::NMultiplets ][ G4ExcitedMesonConstructor::NMesonTypes ] = 
{
  {  142.0*MeV, 360.0*MeV,  91.0*MeV,  90.0*MeV,  90.0*MeV },
  {  265.0*MeV, 350.0*MeV,       0.0, 270.0*MeV, 270.0*MeV },
  {  420.0*MeV,  24.1*MeV,  54.9*MeV, 174.0*MeV, 174.0*MeV },
  {  107.0*MeV, 186.7*MeV,  73.0*MeV,  98.5*MeV,  98.5*MeV },
  {  260.0*MeV, 181.0*MeV, 225.0*MeV, 186.0*MeV, 186.0*MeV },
  {  250.0*MeV, 315.0*MeV,       0.0, 320.0*MeV, 320.0*MeV },
  {  161.0*MeV, 168.0*MeV,  87.0*MeV, 159.0*MeV, 159.0*MeV },
  {  400.0*MeV,  55.0*MeV,  85.0*MeV, 260.0*MeV, 260.0*MeV },
  {  400.0*MeV, 215.0*MeV, 150.0*MeV, 232.0*MeV, 232.0*MeV },
  {        0.0, 197.0*MeV, 200.0*MeV, 373.0*MeV, 373.0*MeV }
};


const G4int    G4ExcitedMesonConstructor::iIsoSpin[] =
{
//  Tpi  TEta  TEtaPrime   TK  TAntiK
     2,     0,         0,   1,     1
};

const G4int    G4ExcitedMesonConstructor::iSpin[] =
{
//N   1     1     1     1     1     1     1     2     2     2
//    1P1   3P0   3P1   3P2   1D2   3D1   3D3   1S0   3S1   3P2
      2,    0,    2,    4,    4,    2,    6,    0,    2,    4
};

const G4int    G4ExcitedMesonConstructor::iParity[] =
{
//N   1     1     1     1     1     1     1     2     2     2
//    1P1   3P0   3P1   3P2   1D2   3D1   3D3   1S0   3S1   3P2
     +1,   +1,   +1,   +1,   -1,   -1,   -1,   -1,   -1,   +1
};

const G4int    G4ExcitedMesonConstructor::iChargeConjugation[] =
{
//N   1     1     1     1     1     1     1     2     2     2
//    1P1   3P0   3P1   3P2   1D2   3D1   3D3   1S0   3S1   3P2
     -1,   +1,   +1,   +1,   +1,   -1,   -1,   +1,   -1,   +1
};

const G4int    G4ExcitedMesonConstructor::iGParity[G4ExcitedMesonConstructor::NMultiplets ][ G4ExcitedMesonConstructor::NMesonTypes ]=
{
  {  +1,  -1,  -1,  0,  0},
  {  -1,  +1,   0,  0,  0},
  {  -1,  +1,  +1,  0,  0},
  {  -1,  +1,  +1,  0,  0},
  {  -1,  +1,  +1,  0,  0},
  {  +1,  -1,   0,  0,  0},
  {  +1,  -1,  -1,  0,  0},
  {  -1,  +1, +1,  0,  0},
  {  +1,  -1,  -1,  0,  0},
  {   0,  +1,  +1,  0,  0}
};


const G4int    G4ExcitedMesonConstructor::encodingOffset[]=
{ 10000, 10000, 20000,      0, 10000, 30000,     0, 100000,100000,100000};




const G4double G4ExcitedMesonConstructor::bRatio[G4ExcitedMesonConstructor::NMultiplets][G4ExcitedMesonConstructor::NMesonTypes][G4ExcitedMesonConstructor::NumberOfDecayModes]=
{
   //    0    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15
   // "b1(1235)",   "h1(1170)",   "h1(1380)",      "k1(1270)",      "k1(1270)" 
  {
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.00, 0.90, 0.10, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.47, 0.42, 0.11, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.47, 0.42, 0.11, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 }
 },
   //    0    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15
   // "a0(1450)",   "f0(1370)",           "", "k0_star(1430)", "k0_star(1430)" 
  {
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.90, 0.00, 0.00, 0.00, 0.10, 0.00, 0.00, 0.00, 0.10, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.10, 0.00, 0.00, 0.00, 0.70, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.20, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 }
  },
   //    0    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15
   // "a1(1260)",   "f1(1285)",   "f1(1420)",      "k1(1400)",      "k1(1400)" 
  {
   { 0.10, 0.00, 0.00, 0.90, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.07, 0.00, 0.00, 0.00, 0.00, 0.20, 0.00, 0.54, 0.00, 0.10, 0.00, 0.00, 0.00, 0.09, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.96, 0.03, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.96, 0.03, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 }
  },
   //    0    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15
   //"a2(1320)",   "f2(1270)","f2_prime(1525)","k2_star(1430)","k2_star(1430)"
  {
   { 0.00, 0.00, 0.00, 0.70, 0.00, 0.14, 0.00, 0.00, 0.00, 0.00, 0.00, 0.11, 0.00, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.30, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.20, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.10, 0.89, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.50, 0.25, 0.09, 0.03, 0.13, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.50, 0.25, 0.09, 0.03, 0.13, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 }
  },
   //    0    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15
   // "pi2(1670)", "eta2(1645)", "eta2(1870)",      "k2(1770)",      "k2(1770)" 
  {
   { 0.00, 0.00, 0.00, 0.30, 0.00, 0.00, 0.00, 0.04, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.56, 0.10, 0.00 },
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.10, 0.00, 0.00, 0.00, 0.90 },
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50 },
   { 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 }
  },
   //    0    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15
   // "rho(1700)", "omega(1650)",          "",  "k_star(1680)",  "k_star(1680)" 
  {
   { 0.00, 0.00, 0.10, 0.00, 0.00, 0.20, 0.00, 0.00, 0.00, 0.00, 0.70, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.40, 0.30, 0.30, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.40, 0.30, 0.30, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 }
  },
   //    0    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15
   // "rho3(1690)","omega3(1670)","phi3(1850)", "k3_star(1780)", "k3_star(1780)"
  {
   { 0.00, 0.00, 0.24, 0.00, 0.00, 0.00, 0.60, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.01, 0.04, 0.11, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.40, 0.00, 0.00, 0.00, 0.00, 0.00, 0.60, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.19, 0.20, 0.31, 0.00, 0.00, 0.00, 0.30, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.19, 0.20, 0.31, 0.00, 0.00, 0.00, 0.30, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 }
  },
   //    0    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15
   // "pi(1300)",  "eta(1295)",  "eta(1475)",       "k(1460)",       "k(1460)" 
  {
   { 0.00, 0.00, 0.00, 0.50, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.20, 0.20, 0.00, 0.00, 0.00, 0.00, 0.00, 0.60, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.50, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.50, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 }
  },
   //    0    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15
   // "rho(1450)","omega(1420)",  "phi(1680)",  "k_star(1410)",  "k_star(1410)" 
  {
   { 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.80, 0.00, 0.00, 0.00, 0.00, 0.00, 0.10, 0.10, 0.00, 0.00, 0.00, 0.00 },
   { 0.30, 0.65, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.30, 0.65, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 }
  },
   //    0    1     2     3     4     5     6     7     8     9    10    11    12    13    14    15
   //         "",   "f2(1810)",   "f2(2010)", "k2_star(1980)", "k2_star(1980)" 
  {
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.30, 0.00, 0.00, 0.00, 0.00, 0.20, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.60, 0.40, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 },
   { 0.00, 0.00, 0.60, 0.40, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 }
  }
};






