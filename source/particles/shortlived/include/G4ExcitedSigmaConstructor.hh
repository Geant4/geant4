// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedSigmaConstructor.hh,v 1.3 1999-12-15 14:51:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
//      History: first implementation, based on object model of
//      10 oct 1998  H.Kurashige
// ---------------------------------------------------------------
#ifndef G4ExcitedSigmaConstructor_h
#define G4ExcitedSigmaConstructor_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ExcitedBaryonConstructor.hh"

class G4ExcitedSigmaConstructor: public G4ExcitedBaryonConstructor
{
  //This class is a utility class for construction 
  //short lived particles

  public:
    G4ExcitedSigmaConstructor();
    virtual  ~G4ExcitedSigmaConstructor();

  protected:  
    virtual  G4bool   Exist( G4int ){return true;}

    virtual  G4int    GetQuarkContents(G4int, G4int);
    virtual  G4String GetName(G4int iIso3, G4int iState);

    virtual  G4double GetMass(G4int iState);
    virtual  G4double GetWidth(G4int iState);
    virtual  G4int    GetiSpin(G4int iState);
    virtual  G4int    GetiParity(G4int iState);
    virtual  G4int    GetEncodingOffset(G4int iState);

    virtual  G4DecayTable* CreateDecayTable(const G4String& name,
					    G4int iIso3, G4int iState,
					    G4bool fAnti = false);
  private:
    G4DecayTable* AddNKMode( G4DecayTable* table, const G4String& name,
				    G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddNKStarMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddSigmaPiMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddSigmaStarPiMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddLambdaPiMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddSigmaEtaMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddLambdaStarPiMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddDeltaKMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);

  public:
    enum     { NStates = 8  };
  private:
    enum     { SigmaIsoSpin = 2 };

  private:
    static const char* name[ NStates ];
    static const G4double mass[ NStates ];
    static const G4double width[ NStates ];
    static const G4int    iSpin[ NStates ];
    static const G4int    iParity[ NStates ];
    static const G4int    encodingOffset[ NStates ];
   
  public:
    enum     { NumberOfDecayModes = 8 };
  private:
    enum     { NK=0,  NKStar=1,   SigmaPi=2, SigmaStarPi=3,  LambdaPi=4,
	       SigmaEta=5,  LambdaStarPi=6, DeltaK=7};
   private:
   static const G4double bRatio[ NStates ][ NumberOfDecayModes];
};

inline
 G4double G4ExcitedSigmaConstructor::GetMass(G4int iState)
{ 
  return mass[iState]; 
}

inline
 G4double G4ExcitedSigmaConstructor::GetWidth(G4int iState)
{
  return width[iState];
}

inline
 G4int    G4ExcitedSigmaConstructor::GetiSpin(G4int iState)
{
  return iSpin[iState];
}

inline
 G4int    G4ExcitedSigmaConstructor::GetiParity(G4int iState)
{
  return iParity[iState];
}

inline
 G4int    G4ExcitedSigmaConstructor::GetEncodingOffset(G4int iState)
{
  return encodingOffset[iState];
}

inline
 G4int  G4ExcitedSigmaConstructor::GetQuarkContents(G4int iQ, G4int iIso3)
{
  G4int quark=0;
  if ( iQ == 0 ){
    // s-quark
    quark = 3;
  } else if ( iQ == 1 ){
    if (iIso3 == -2) {
      // d-quark
      quark = 1;
    } else {
      // u-quark
      quark = 2;
    }
  }  else if ( iQ == 2 ){
    if (iIso3 == +2) {
      // u-quark
      quark = 2;
    } else {
      // d-quark
      quark = 1;
    }
  } 
  return quark;
}

inline 
 G4String  G4ExcitedSigmaConstructor::GetName(G4int iIso3, G4int iState)
 {
   G4String particle = name[iState];
   if (iIso3 == +2) {
     particle += "+";
   } else if (iIso3 == 0) {
     particle += "0";
   } else if (iIso3 == -2) {
     particle += "-";
  }
  return particle;
}
#endif










