// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitedXiConstructor.hh,v 1.1 1999-01-07 16:10:39 gunter Exp $
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
#ifndef G4ExcitedXiConstructor_h
#define G4ExcitedXiConstructor_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ExcitedBaryonConstructor.hh"

class G4ExcitedXiConstructor: public G4ExcitedBaryonConstructor
{
  //This class is a utility class for construction 
  //short lived particles

  public:
    G4ExcitedXiConstructor();
    virtual  ~G4ExcitedXiConstructor();

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
    G4DecayTable* AddXiPiMode( G4DecayTable* table, const G4String& name,
				    G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddXiGammaMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddLambdaKMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);
    G4DecayTable* AddSigmaKMode( G4DecayTable* table, const G4String& name,
				     G4double br, G4int iIso3, G4bool fAnti);

  public:   
    enum     { NStates = 5  };
  private:
    enum     { XiIsoSpin = 1 };

  private:
    static const G4String name[ NStates ];
    static const G4double mass[ NStates ];
    static const G4double width[ NStates ];
    static const G4int    iSpin[ NStates ];
    static const G4int    iParity[ NStates ];
    static const G4int    encodingOffset[ NStates ];

  public:   
    enum     { NumberOfDecayModes = 4 };
  private:
    enum     { XiPi=0,  XiGamma=1,  LambdaK=2, SigmaK=3 };
  private:
   static const G4double bRatio[ NStates ][ NumberOfDecayModes];
};

inline
 G4double G4ExcitedXiConstructor::GetMass(G4int iState)
{ 
  return mass[iState]; 
}

inline
 G4double G4ExcitedXiConstructor::GetWidth(G4int iState)
{
  return width[iState];
}

inline
 G4int    G4ExcitedXiConstructor::GetiSpin(G4int iState)
{
  return iSpin[iState];
}

inline
 G4int    G4ExcitedXiConstructor::GetiParity(G4int iState)
{
  return iParity[iState];
}

inline
 G4int    G4ExcitedXiConstructor::GetEncodingOffset(G4int iState)
{
  return encodingOffset[iState];
}

inline
 G4int  G4ExcitedXiConstructor::GetQuarkContents(G4int iQ, G4int iIso3)
{
  G4int quark=0;
  if ( iQ == 0 ){
    // s-quark
    quark = 3;
  } else if ( iQ == 1 ){
    // s-quark
    quark = 3;
  }  else if ( iQ == 2 ){
    if (iIso3 == +1) {
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
 G4String  G4ExcitedXiConstructor::GetName(G4int iIso3, G4int iState)
{
  G4String particle = name[iState];
  if (iIso3 == +1) {
    particle += "0";
  } else if (iIso3 ==-1) {
    particle += "-";
  }
  return particle;
}
#endif










