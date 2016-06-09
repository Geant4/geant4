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
// $Id: G4QPDGCode.hh,v 1.26 2006/06/29 20:06:23 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//      ---------------- G4QPDGCode ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for Hadron definition in CHIPS Model
// ------------------------------------------------------------

#ifndef G4QPDGCode_h
#define G4QPDGCode_h 1

#include <iostream>
#include "globals.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4ParticleTypes.hh"

#include "G4QContent.hh"

class G4QPDGCode
{
public:
  // Constructors
  G4QPDGCode(G4int PDGCode = 0);                     // Construction by PDGCode
  G4QPDGCode(G4bool f, G4int QCode);                 // Construction by QCode, f-verb
  G4QPDGCode(G4QContent QCont);                      // Construction by Quark Content
  G4QPDGCode(const G4QPDGCode& rhs);                 // Copy Constructor by value
  G4QPDGCode(G4QPDGCode* rhs);                       // Copy Constructor by pointer

  ~G4QPDGCode();                                     // Public Destructor

  // Operators
  const G4QPDGCode& operator=(const G4QPDGCode& rhs);
  G4bool            operator==(const G4QPDGCode& rhs) const;
  G4bool            operator==(const G4int& rhs) const;
  G4bool            operator!=(const G4QPDGCode& rhs) const;
  G4bool            operator!=(const G4int& rhs) const;
  G4QPDGCode        operator+=(const G4int& rhs);
  G4QPDGCode        operator+=(const G4QPDGCode& rhs);
  G4QPDGCode        operator-=(const G4int& rhs);
  G4QPDGCode        operator-=(const G4QPDGCode& rhs);
  G4QPDGCode        operator*=(const G4int& rhs);
  G4QPDGCode        operator/=(const G4int& rhs);

  // Selectors
  G4int      GetNQHadr();                              // Return # of predefined hadrons
  G4double   GetMass();                                // GS Mass for the QHadron
  G4double   GetMass2();                               // Squared GS Mass for the QHadron
  G4double   GetWidth();                               // Width for the QHadron
  G4double   GetNuclMass(G4int Z, G4int N, G4int S);   // Wrapper forNuclearMassCalculation
  G4double   GetNuclMass(G4int PDGCode);               // Wrapper forNuclearMassCalculation
  G4QContent GetQuarkContent()                  const; // Get QC for the particle
  G4int      GetBaryNum()                       const; // Get Baryon Number of the Hadron
  G4int      GetSpin()                          const; // Returns 2s+1 for hadrons, 1 for A
  G4int      GetCharge()                        const; // Get Charge of the Hadron
  G4int      GetPDGCode()                       const; // Get PDG code of the Hadron
  G4int      GetQCode()                         const; // Get Q code of the Hadron
  G4QContent GetExQContent(G4int i, G4int o)    const; // Get Q Content for Quark Exchange
  G4int      GetRelCrossIndex(G4int i, G4int o) const; // Relative Cross Index for q_i->q_o
  G4int      GetNumOfComb(G4int i, G4int o)     const; // Get #ofCombinations for q_i->q_o
  G4int      GetTotNumOfComb(G4int i)           const; // Get total#ofCombinations for q_i 

  // Modifiers
  void       SetPDGCode(G4int newPDGCode);             // Set PDG code of the Hadron
  void       InitByQCont(G4QContent QCont);            // Init ExistingQPDG by QuarkContent
  void       InitByQCode(G4int QCode);                 // Init ExistingQPDG by Q Code

  // General
  G4bool     TestRealNeutral();
  void       NegPDGCode();
  void       ConvertPDGToZNS(G4int PDG, G4int& z, G4int& n, G4int& s);

private:
  // Encapsulated functions
  G4bool   TestRealNeutral(const G4int& PDGCode);
  G4int    MakeQCode(const G4int& PDGCode);              // Make Q Code, using PDG Code
  G4int    MakePDGCode(const G4int& QCode);              // Make PDG Code, using Q Code
  G4double CalculateNuclMass(G4int Z, G4int N, G4int S); // Nuclear Mass Calculation
  G4double QHaM(G4int nQ);                      // Definition of hadronic masses in Q-order
private:
  // Static parameter
  static const G4int nQHM=90;
  // the Body
  G4int              thePDGCode;
  G4int              theQCode;
};

// Not member operators
std::ostream&   operator<<(std::ostream& lhs, G4QPDGCode& rhs);
std::ostream&   operator<<(std::ostream& lhs, const G4QPDGCode& rhs);
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

inline G4bool G4QPDGCode::operator==(const G4QPDGCode& rhs) const {return this==&rhs;}
inline G4bool G4QPDGCode::operator==(const G4int&      rhs) const {return thePDGCode==rhs;}
inline G4bool G4QPDGCode::operator!=(const G4QPDGCode& rhs) const {return this!=&rhs;}
inline G4bool G4QPDGCode::operator!=(const G4int&      rhs) const {return thePDGCode!=rhs;}

inline G4int  G4QPDGCode::GetNQHadr() {return nQHM;} // Return # of predefined hadrons

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
  if(PDG>80000000)
  {
    G4int szn=PDG-90000000;
    G4int ds=0;
    G4int dz=0;
    G4int dn=0;
    if(szn<-100000)
    {
      G4int ns=(-szn)/1000000+1;
      szn+=ns*1000000;
      ds+=ns;
    }
    else if(szn<-100)
    {
      G4int nz=(-szn)/1000+1;
      szn+=nz*1000;
      dz+=nz;
    }
    else if(szn<0)
    {
      G4int nn=-szn;
      szn=0;
      dn+=nn;
    }
    G4int sz =szn/1000;
    G4int n  =szn%1000;
    if(n>700)
    {
      n-=1000;
      dz--;
    }
    G4int z  =sz%1000-dz;
    if(z>700)
    {
      z-=1000;
      ds--;
    }
    G4int s  =sz/1000-ds;
    return GetNuclMass(z,n,s);
  }
  return 0.;
}
inline G4int   G4QPDGCode::GetPDGCode() const {return thePDGCode;}
inline G4int   G4QPDGCode::GetQCode()   const {return theQCode;}
inline G4int   G4QPDGCode::GetCharge()  const {return GetQuarkContent().GetCharge();}
inline G4int   G4QPDGCode::GetBaryNum() const {return GetQuarkContent().GetBaryonNumber();}
inline G4int   G4QPDGCode::GetSpin()    const 
{
  if(thePDGCode<80000000)               return thePDGCode%10;
  else if(GetQuarkContent().GetTot()%2) return 3; // @@ Take into account higher resonances
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



