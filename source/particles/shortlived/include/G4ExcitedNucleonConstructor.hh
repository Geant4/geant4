// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedNucleonConstructor.hh,v 1.1 1999-01-07 16:10:38 gunter Exp $
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
#ifndef G4ExcitedNucleonConstructor_h
#define G4ExcitedNucleonConstructor_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ExcitedBaryonConstructor.hh"

class G4ExcitedNucleonConstructor: public G4ExcitedBaryonConstructor
{
  //This class is a utility class for construction 
  //short lived particles

  public:
    G4ExcitedNucleonConstructor();
    virtual  ~G4ExcitedNucleonConstructor();

  protected:  
    virtual  G4int    GetEncoding(G4int iIsoSpin3, G4int idxState);

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
    G4DecayTable* AddNGammaMode( G4DecayTable* table, const G4String& name,
				    G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddNPiMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddNEtaMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddNOmegaMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddNRhoMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddN2PiMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddDeltaPiMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddNStarPiMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddLambdaKMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);

  public:
    enum     { NStates = 15  };
  private:
   enum     { NucleonIsoSpin = 1  };

  private:
    static const G4String name[ NStates ];
    static const G4double mass[ NStates ];
    static const G4double width[ NStates ];
    static const G4int    iSpin[ NStates ];
    static const G4int    iParity[ NStates ];
    static const G4int    encodingOffset[ NStates ];
   
  public:
    enum     { NumberOfDecayModes = 9 };
  private:
    enum     { NGamma=0,  NPi=1,   NEta=2, NOmega=3,  NRho=4,
	       N2Pi=5,  DeltaPi=6, NStarPi=7, LambdaK = 8 };
  private:
    static const G4double bRatio[ NStates ][ NumberOfDecayModes];
};

inline
 G4double G4ExcitedNucleonConstructor::GetMass(G4int iState)
{ 
  return mass[iState]; 
}

inline
 G4double G4ExcitedNucleonConstructor::GetWidth(G4int iState)
{
  return width[iState];
}

inline
 G4int    G4ExcitedNucleonConstructor::GetiSpin(G4int iState)
{
  return iSpin[iState];
}

inline
 G4int    G4ExcitedNucleonConstructor::GetiParity(G4int iState)
{
  return iParity[iState];
}

inline
 G4int    G4ExcitedNucleonConstructor::GetEncodingOffset(G4int iState)
{
  return encodingOffset[iState];
}

inline
 G4int  G4ExcitedNucleonConstructor::GetQuarkContents(G4int iQ, G4int iIso3)
{
  // Quark contents
  //    iIso3 = -1 : udd
  //    iIso3 = +1 : uud
  G4int quark=0;
  if ( iQ == 0 ){
    // u-quark
    quark = 2;
  } else if ( iQ == 2 ){
    // d-quark
    quark = 1;
  } else {
    if ( iIso3 == -1 ){
    // d-quark
      quark = 1;
    } else {
    // u-quark
      quark = 2;
    }
  }
  return quark;
}

inline 
 G4String  G4ExcitedNucleonConstructor::GetName(G4int iIso3, G4int iState)
{
  G4String particle = name[iState];
  if ( iIso3 == -1 ){
    particle += "0";
  } else {
    particle += "+";
  }
  return particle;
}
#endif










