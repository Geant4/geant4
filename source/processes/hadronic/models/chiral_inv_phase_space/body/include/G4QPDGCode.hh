// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QPDGCode.hh,v 1.3 2000-09-10 13:58:55 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4QPDGCode_h
#define G4QPDGCode_h 1

// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QPDGCode ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for Hadron definition in CHIPS Model
// ------------------------------------------------------------

#include <iostream.h>
#include "globals.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4QContent.hh"

class G4QPDGCode
{
public:
  // Constructors
  G4QPDGCode(G4int PDGCode = 0);                     // Construction by PDGCode
  G4QPDGCode(G4QContent QCont);                      // Construction by Quark Content
  G4QPDGCode(const G4QPDGCode& rhs);                 // Copy Constructor by value
  G4QPDGCode(G4QPDGCode* rhs);                       // Copy Constructor by pointer

  ~G4QPDGCode();                                     // Destructor

  // Operators
  const G4QPDGCode& operator=(const G4QPDGCode& rhs);
  G4int             operator==(const G4QPDGCode& rhs) const;
  G4int             operator==(const G4int& rhs) const;
  G4int             operator!=(const G4QPDGCode& rhs) const;
  G4int             operator!=(const G4int& rhs) const;
  G4QPDGCode        operator+=(const G4int& rhs);
  G4QPDGCode        operator+=(const G4QPDGCode& rhs);
  G4QPDGCode        operator-=(const G4int& rhs);
  G4QPDGCode        operator-=(const G4QPDGCode& rhs);
  G4QPDGCode        operator*=(const G4int& rhs);
  G4QPDGCode        operator/=(const G4int& rhs);

  // Selectors
  G4double   GetMass();                                // GS Mass for the QHadron
  G4double   GetMass2();                               // Squared GS Mass for the QHadron
  G4double   GetWidth();                               // Width for the QHadron
  G4double   GetNuclMass(G4int Z, G4int N, G4int S);   // Wrapper for nuclear mass calculation
  G4double   GetNuclMass(G4int PDGCode);               // Wrapper for nuclear mass calculation
  G4QContent GetQuarkContent()                  const; // Get QC for the particle
  G4int      GetBaryNum()                       const; // Get Baryon Number of the Hadron
  G4int      GetSpin()                          const; // Returns 2s+1 for hadrons or 1 for nuclei
  G4int      GetCharge()                        const; // Get Charge of the Hadron
  G4int      GetPDGCode()                       const; // Get PDG code of the Hadron
  G4int      GetQCode()                         const; // Get Q code of the Hadron
  G4QContent GetExQContent(G4int i, G4int o)    const; // Get Q Content for the Quark Exchange  
  G4int      GetRelCrossIndex(G4int i, G4int o) const; // Relative Cross Index for q_i->q_o
  G4int      GetNumOfComb(G4int i, G4int o)     const; // Get number of Combinations for q_i->q_o
  G4int      GetTotNumOfComb(G4int i)           const; // Get a total # of Combinations for q_i 

  // Modifiers
  void  SetPDGCode(G4int newPDGCode);               // Set PDG code of the Hadron
  void  InitByQCont(G4QContent QCont);              // Init existing QPDG by Quark Content
  void  InitByQCode(G4int QCode);                   // Init existing QPDG by Q Code

  // General
  G4bool     TestRealNeutral();
  void       NegPDGCode();

private:
  // Encapsulated functions
  G4bool   TestRealNeutral(const G4int& PDGCode);
  G4int    MakeQCode(const G4int& PDGCode);         // Make Q Code, using PDG Code
  G4int    MakePDGCode(const G4int& QCode);         // Make PDG Code, using Q Code

private:
  // the Body
  G4int              thePDGCode;
  G4int              theQCode;
};

// Not member operators
ostream&   operator<<(ostream& lhs, G4QPDGCode& rhs);
ostream&   operator<<(ostream& lhs, const G4QPDGCode& rhs);
G4int      operator+(const G4QPDGCode& lhs, const G4QPDGCode& rhs);
G4int      operator+(const G4QPDGCode& lhs, const G4int&      rhs);
G4int      operator+(const G4int&      lhs, const G4QPDGCode& rhs);
G4int      operator-(const G4QPDGCode& lhs, const G4QPDGCode& rhs);
G4int      operator-(const G4QPDGCode& lhs, const G4int&      rhs);
G4int      operator-(const G4int&      lhs, const G4QPDGCode& rhs);
G4int      operator*(const G4QPDGCode& lhs, const G4QPDGCode& rhs);
G4int      operator*(const G4QPDGCode& lhs, const G4int&      rhs);
G4int      operator*(const G4int&      lhs, const G4QPDGCode& rhs);
G4int      operator/(const G4QPDGCode& lhs, const G4QPDGCode& rhs);
G4int      operator/(const G4QPDGCode& lhs, const G4int&      rhs);
G4int      operator/(const G4int&      lhs, const G4QPDGCode& rhs);
G4int      operator%(const G4QPDGCode& lhs, const G4int&      rhs);
// Not member functions
//----------------------------------------------------------------------------------------

inline G4int G4QPDGCode::operator==(const G4QPDGCode& rhs) const {return this==&rhs;}
inline G4int G4QPDGCode::operator==(const G4int&      rhs) const {return thePDGCode==rhs;}
inline G4int G4QPDGCode::operator!=(const G4QPDGCode& rhs) const {return this!=&rhs;}
inline G4int G4QPDGCode::operator!=(const G4int&      rhs) const {return thePDGCode!=rhs;}

inline G4QPDGCode G4QPDGCode::operator+=(const G4QPDGCode& rhs)
{
  thePDGCode+=rhs.GetPDGCode();
  if(!thePDGCode) theQCode = -2;
  else theQCode = MakeQCode(thePDGCode);
  return *this;
}
inline G4QPDGCode G4QPDGCode::operator+=(const G4int& rhs)
{
  thePDGCode+=rhs;
  if(!thePDGCode) theQCode = -2;
  else theQCode = MakeQCode(thePDGCode);
  return *this;
}
inline G4QPDGCode G4QPDGCode::operator-=(const G4QPDGCode& rhs)
{
  thePDGCode-=rhs.GetPDGCode();
  if(!thePDGCode) theQCode = -2;
  else theQCode = MakeQCode(thePDGCode);
  return *this;
}
inline G4QPDGCode G4QPDGCode::operator-=(const G4int& rhs)
{
  thePDGCode-=rhs;
  if(!thePDGCode) theQCode = -2;
  else theQCode = MakeQCode(thePDGCode);
  return *this;
}
inline G4QPDGCode G4QPDGCode::operator*=(const G4int& rhs)
{
  thePDGCode*=rhs;
  if(!thePDGCode) theQCode = -2;
  else theQCode = MakeQCode(thePDGCode);
  return *this;
}
inline G4QPDGCode G4QPDGCode::operator/=(const G4int& rhs)
{
  thePDGCode/=rhs;
  if(!thePDGCode) theQCode = -2;
  else theQCode = MakeQCode(thePDGCode);
  return *this;
}
 
inline G4double G4QPDGCode::GetMass2() {G4double m=GetMass(); return m*m;}
inline G4double G4QPDGCode::GetNuclMass(G4int PDG)
{
  if(PDG>90000000)
  {
    G4int zns=PDG-90000000;
    G4int zs=zns/1000;
    G4int n =zns%1000;
    G4int s = zs/1000;
    G4int z = zs%1000;
    return GetNuclMass(z,n,s);
  }
  return 0.;
}
inline G4int    G4QPDGCode::GetPDGCode() const {return thePDGCode;}
inline G4int    G4QPDGCode::GetQCode()   const {return theQCode;}
inline G4int    G4QPDGCode::GetCharge()  const {return GetQuarkContent().GetCharge();}
inline G4int    G4QPDGCode::GetBaryNum() const {return GetQuarkContent().GetBaryonNumber();}
inline G4int    G4QPDGCode::GetSpin()    const 
{
  if(thePDGCode<90000000)               return thePDGCode%10;
  else if(GetQuarkContent().GetTot()%2) return 3;
  else                                  return 1;
}
inline void     G4QPDGCode::NegPDGCode()     {thePDGCode=-thePDGCode;}
inline G4bool   G4QPDGCode::TestRealNeutral(){return TestRealNeutral(thePDGCode);}

// Redefinition of the PDG instance
inline void  G4QPDGCode::SetPDGCode(G4int newPDGCode)
//           ========================================
{
  thePDGCode=newPDGCode;
  if(!thePDGCode) theQCode = -2;
  else theQCode = MakeQCode(newPDGCode);
}

// Init existing QPDG by Quark Content
inline void  G4QPDGCode::InitByQCont(G4QContent QCont)
{//          =========================================
  thePDGCode = QCont.GetSPDGCode();
  if(!thePDGCode) theQCode = -2;
  else theQCode   = MakeQCode(thePDGCode);
}

// Init existing QPDG by Quark Content
inline void  G4QPDGCode::InitByQCode(G4int QCode)
{//          ====================================
  theQCode   = QCode;
  thePDGCode = MakePDGCode(QCode);
}

#endif



