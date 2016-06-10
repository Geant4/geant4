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
// This code plementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedMesonConstructor.hh 67971 2013-03-13 10:13:24Z gcosmo $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//      History: first implementation, based on object model of
//      10 oct 1998  H.Kurashige
// ---------------------------------------------------------------
#ifndef G4ExcitedMesonConstructor_h
#define G4ExcitedMesonConstructor_h 1

#include "globals.hh"
#include "G4ios.hh"
class     G4DecayTable;

class G4ExcitedMesonConstructor
{ 
  //This class is a utility class for construction 
  //short lived particles

  public:
    G4ExcitedMesonConstructor(G4int nStates = 0, G4int isoSpin=0);
    virtual  ~G4ExcitedMesonConstructor();
  
  public:
    virtual  void Construct(G4int indexOfState = -1);
 
  protected:
    void ConstructMesons(G4int indexOfState, G4int indexOfType);
     
    G4String GetName(G4int iIso3, G4int iState, G4int idxType);
    G4double GetCharge(G4int iIsoSpin3);
    G4int    GetEncoding(G4int iIsoSpin3, G4int idxState, G4int idxType);
    G4int    GetQuarkContents(G4int iQ, G4int iIso3,  G4int iType);

  public:
    enum { NMultiplets = 10 };
  protected:    
    enum { 
      N11P1 = 0, N13P0 = 1, N13P1 = 2, N13P2 = 3,
      N11D2 = 4, N13D1 = 5, N13D3 = 6,
      N21S0 = 7, N23S1 = 8, N23P2 = 9
    };
    
  public:
    enum { NMesonTypes = 5 };
  protected:    
    enum { TPi=0, TEta=1, TEtaPrime=2, TK=3, TAntiK=4 }; 

  protected:    
    const G4String type;
    const G4int    leptonNumber;
    const G4int    baryonNumber;

    G4bool Exist(G4int idxState, G4int idxType);
    G4double GetCharge(G4int iIsoSpin3, G4int idxType);
    static const char* name[ NMultiplets ][ NMesonTypes ];
    static const G4double mass[ NMultiplets ][ NMesonTypes ];
    static const G4double massKdiff[ NMultiplets ];
    static const G4double width[ NMultiplets ][ NMesonTypes ];
    static const G4double widthKdiff[ NMultiplets ];
    static const G4int    iIsoSpin[ NMesonTypes ];    
    static const G4int    iSpin[ NMultiplets ];
    static const G4int    iParity[ NMultiplets ];
    static const G4int    iGParity[ NMultiplets ][ NMesonTypes ];
    static const G4int    iChargeConjugation[ NMultiplets ];
    static const G4int    encodingOffset[ NMultiplets ];
   
  public:
    enum     { NumberOfDecayModes = 19 };
  protected:    
    enum     { MPiGamma = 0, MRhoGamma=1, M2Pi=2,      MPiRho=3, 
               M3Pi= 4,      MPiEta=5,    M4Pi=6,      MKKStar=7,
               M2PiEta=8,    MRhoEta=9,   M2PiRho=10,  M2PiOmega=11,
               M2Eta=12,     M2K=13,      M2KPi=14,    MPiOmega=15,
               MPiF2=16,     MPiF0=17,    MPiA2=18 };
    enum     { MKPi = 0,     MKStarPi=1,  MKRho=2,     MKOmega=3,
               MKStar2Pi=4,  MKTwoPi=5,   MKEta=6}; 
               
               
    static const G4double bRatio[ NMultiplets ][ NMesonTypes ][ NumberOfDecayModes];

    G4DecayTable* CreateDecayTable(const G4String&,
					G4int , G4int, G4int);
    
    G4DecayTable* AddKPiMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iType);
    G4DecayTable* AddKStarPiMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iType);
    G4DecayTable* AddKStar2PiMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iType);
    G4DecayTable* AddKRhoMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iType);
    G4DecayTable* AddKTwoPiMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iType);
    G4DecayTable* AddKOmegaMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iType);
    G4DecayTable* AddKEtaMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iType);
    G4DecayTable* AddPiGammaMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3,G4int iIso);
    G4DecayTable* AddRhoGammaMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* Add2PiMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* AddPiRhoMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* AddPiEtaMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* AddPiF2Mode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* AddPiF0Mode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* AddPiA2Mode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* Add3PiMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* Add4PiMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* AddKKStarMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* Add2PiEtaMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* AddRhoEtaMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* Add2PiRhoMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* Add2PiOmegaMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* AddPiOmegaMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* Add2EtaMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* Add2KMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
    G4DecayTable* Add2KPiMode( G4DecayTable* table, const G4String& name,
				        G4double br, G4int iIso3, G4int iIso);
	       
    	       

};


inline 
  G4String G4ExcitedMesonConstructor::GetName(G4int iIso3, 
					      G4int iState, 
					      G4int iType)
{
  G4String particle = name[iState][iType];
  if (iType == TPi) {
    if ( iIso3 == +2 ){
      particle += "+";
    } else if ( iIso3 == -2 ){
      particle += "-";
    } else {
      particle += "0";
    }
  } else if (iType == TK) {
    if ( iIso3 == +1 ){
      particle += "+";
    } else if ( iIso3 == -1 ){
      particle += "0";
    }
  }  else if (iType == TAntiK) {
    if ( iIso3 == +1 ){
      particle += "0";
      particle = "anti_" + particle;
    } else if ( iIso3 == -1 ){
      particle += "-";
    }
  } 
  return particle;
}

#endif
