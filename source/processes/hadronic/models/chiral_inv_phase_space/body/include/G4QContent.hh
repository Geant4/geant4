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
// $Id$
//
//      ---------------- G4QContent ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for Quasmon initiated Contents used by the CHIPS Model
// ----------------------------------------------------------------------
// Short description: This is the basic class of the CHIPS model. It
// describes the quark content of the Quasmon, which is a generalized
// hadronic state. All Quasmons are bags, characterized by the quark
// Content (QContent), but the spin is not fixed and only light (u,d,s)
// quarks are considered (SU(3)). The hadrons are the ground states for
// the corresponding quasmons. The Chipolino (G4QChipolino) or nuclear
// cluster are examples for another Quark Content.
// --------------------------------------------------------------------
// @@ In future total spin & c,b,t of the Hadron can be added @@ M.K.@@
// --------------------------------------------------------------------

#ifndef G4QContent_h
#define G4QContent_h 1

#include <iostream>
#include "globals.hh"
#include "Randomize.hh"
#include "G4QException.hh"

class G4QContent
{
public:
  G4QContent(G4int d=0, G4int u=0, G4int s=0, G4int ad=0, G4int au=0, G4int as=0);
  G4QContent(std::pair<G4int,G4int> PP); // Create HadronQC from two partons
  G4QContent(const G4QContent& rhs);     // Copy constructor by value
  G4QContent(G4QContent* rhs);           // Copy constructor by pointer

  ~G4QContent();                         // Public Destructor

  // Overloaded operators
  const G4QContent& operator=(const G4QContent& rhs);
  G4bool            operator==(const G4QContent& rhs) const;
  G4bool            operator!=(const G4QContent& rhs) const;
  G4QContent        operator+=(G4QContent& rhs);
  G4QContent        operator-=(G4QContent& rhs);
  G4QContent        operator*=(G4int& rhs);
  G4QContent        operator+=(const G4QContent& rhs);
  G4QContent        operator-=(const G4QContent& rhs);
  G4QContent        operator*=(const G4int& rhs);

  // Selectors  
  G4int             GetCharge() const;
  G4int             GetBaryonNumber() const;
  G4int             GetStrangeness() const;
  G4int             GetSPDGCode() const;
  G4int             GetZNSPDGCode() const;
  G4int             NOfCombinations(const G4QContent& rhs) const;//@@ can be an "operator/"
  G4int             GetQ() const;
  G4int             GetAQ() const;
  G4int             GetTot() const; 
  G4bool            CheckNegative() const;

  G4int GetP() const;       // A#of protons
  G4int GetN() const;       // A#of neutrons
  G4int GetL() const;       // A#of lambdas
  G4int GetAP() const;      // A#of anti-protons
  G4int GetAN() const;      // A#of anti-neutrons
  G4int GetAL() const;      // A#of anti-lambdas

  G4int GetD() const;       // A#of d-quarks
  G4int GetU() const;       // A#of u-quarks
  G4int GetS() const;       // A#of s-quarks
  G4int GetAD() const;      // A#of anti-d-quarks
  G4int GetAU() const;      // A#of anti-u-quarks
  G4int GetAS() const;      // A#of anti-s-quarks

  G4int GetNetD() const;    // A net#of d-quarks
  G4int GetNetU() const;    // A net#of u-quarks
  G4int GetNetS() const;    // A net#of s-quarks
  G4int GetNetAD() const;   // A net#of anti-d-quarks
  G4int GetNetAU() const;   // A net#of anti-u-quarks
  G4int GetNetAS() const;   // A net#of anti-s-quarks

  G4int GetDD() const;      // A#of dd-di-quarks
  G4int GetUU() const;      // A#of uu-di-quarks
  G4int GetSS() const;      // A#of ss-di-quarks
  G4int GetUD() const;      // A#of ud-di-quarks
  G4int GetDS() const;      // A#of ds-di-quarks
  G4int GetUS() const;      // A#of us-di-quarks
  G4int GetADAD() const;    // A#of anti-dd-di-quarks
  G4int GetAUAU() const;    // A#of anti-uu-di-quarks
  G4int GetASAS() const;    // A#of anti-ss-di-quarks
  G4int GetAUAD() const;    // A#of anti-ud-di-quarks
  G4int GetADAS() const;    // A#of anti-ds-di-quarks
  G4int GetAUAS() const;    // A#of anti-us-di-quarks

  std::pair<G4int,G4int> MakePartonPair() const;// Make PDG's of PartonPairs for Mes&Baryon
  G4int AddParton(G4int pPDG) const; // Add parton with pPDG to the hadron with this QC

  // Modifiers
  void              Anti();
  G4QContent        IndQ (G4int ind=0);
  G4QContent        IndAQ(G4int ind=0);
  G4QContent        SplitChipo(G4double mQ);
  G4bool            SubtractHadron(G4QContent h);
  G4bool            SubtractPi0();
  G4bool            SubtractPion();
  G4bool            SubtractKaon(G4double mQ);

  void SetD (G4int n=0);
  void SetU (G4int n=0);
  void SetS (G4int n=0);
  void SetAD(G4int n=0);
  void SetAU(G4int n=0);
  void SetAS(G4int n=0);

  void IncD (G4int n=1);
  void IncU (G4int n=1);
  void IncS (G4int n=1);
  void IncAD(G4int n=1);
  void IncAU(G4int n=1);
  void IncAS(G4int n=1);
  void IncQAQ(const G4int& nQAQ=1, const G4double& sProb = 1.);

  void DecD (G4int n=1);
  void DecU (G4int n=1);
  void DecS (G4int n=1);
  void DecAD(G4int n=1);
  void DecAU(G4int n=1);
  void DecAS(G4int n=1);
  G4int DecQAQ(const G4int& nQAQ=1);

private:
  G4QContent        GetThis() const;

  // Body
private:            //                       C    S
  G4int nD;         // a#of      d-quarks (-1/3)( 0)
  G4int nU;         // a#of      u-quarks (+2/3)( 0)
  G4int nS;         // a#of      s-quarks (-1/3)(+1)
  G4int nAD;        // a#of anti-d-quarks (+1/3)( 0)
  G4int nAU;        // a#of anti-u-quarks (-2/3)( 0)
  G4int nAS;        // a#of anti-s-quarks (+1/3)(-1)
};

// Not member operators
std::ostream&   operator<<(std::ostream& lhs, G4QContent& rhs);
std::ostream&   operator<<(std::ostream& lhs, const G4QContent& rhs);
G4QContent operator+(const G4QContent& lhs, const G4QContent& rhs);
G4QContent operator-(const G4QContent& lhs, const G4QContent& rhs);
G4QContent operator*(const G4QContent& lhs, const G4int&      rhs);
G4QContent operator*(const G4int&      lhs, const G4QContent& rhs);
//G4int      operator/(const G4QContent& lhs, const G4QContent& rhs); // Combinations
// Not member functions
//----------------------------------------------------------------------------------------

inline G4bool     G4QContent::operator==(const G4QContent& rhs) const {return this==&rhs;}
inline G4bool     G4QContent::operator!=(const G4QContent& rhs) const {return this!=&rhs;}
inline G4int      G4QContent::GetQ() const {return nU+nD+nS;}
inline G4int      G4QContent::GetAQ() const {return nAU+nAD+nAS;}
inline G4int      G4QContent::GetTot() const {return nU+nD+nS+nAU+nAD+nAS;}
inline G4int      G4QContent::GetStrangeness() const {return nS-nAS;}
// @@ Temporary for tests
inline G4bool G4QContent::CheckNegative() const
                                            {return nU<0||nD<0||nS<0||nAU<0||nAD<0||nAS<0;}

inline G4int G4QContent::GetU() const{return nU;}
inline G4int G4QContent::GetD() const{return nD;}
inline G4int G4QContent::GetS() const{return nS;}
inline G4int G4QContent::GetAU() const{return nAU;}
inline G4int G4QContent::GetAD() const{return nAD;}
inline G4int G4QContent::GetAS() const{return nAS;}

inline G4int G4QContent::GetNetU() const{return nU-nAU;}
inline G4int G4QContent::GetNetD() const{return nD-nAD;}
inline G4int G4QContent::GetNetS() const{return nS-nAS;}
inline G4int G4QContent::GetNetAU() const{return nAU-nU;}
inline G4int G4QContent::GetNetAD() const{return nAD-nD;}
inline G4int G4QContent::GetNetAS() const{return nAS-nS;}

inline G4int G4QContent::GetUU() const{return nU*(nU-1)/2;}
inline G4int G4QContent::GetDD() const{return nD*(nD-1)/2;}
inline G4int G4QContent::GetSS() const{return nS*(nS-1)/2;}
inline G4int G4QContent::GetUD() const{return nU*nD;}
inline G4int G4QContent::GetUS() const{return nU*nS;}
inline G4int G4QContent::GetDS() const{return nD*nS;}
inline G4int G4QContent::GetAUAU() const{return nAU*(nAU-1)/2;}
inline G4int G4QContent::GetADAD() const{return nAD*(nAD-1)/2;}
inline G4int G4QContent::GetASAS() const{return nAS*(nAS-1)/2;}
inline G4int G4QContent::GetAUAD() const{return nAU*nAD;}
inline G4int G4QContent::GetAUAS() const{return nAU*nAS;}
inline G4int G4QContent::GetADAS() const{return nAD*nAS;}

// Convert particle to anti-particle
inline G4int G4QContent::GetZNSPDGCode() const
{
  G4int kD=nD-nAD;                           // A net # of d quarks
  G4int kU=nU-nAU;                           // A net # of u quarks
  G4int kS=nS-nAS;                           // A net # of s quarks
  // if(kD>=0&&kU>=0&&kS>=0&&kD+kU+kS>0)        // => "Normal nucleus" case
  //{
  //  G4int b=(kU+kD-kS-kS)/3;
  //  G4int d=kU-kD;
  //  G4int n=(b-d)/2;
  //  return 90000000+1000*(1000*kS+n+d)+n;
  //}
  //else if(kD<=0&&kU<=0&&kS<=0&&kD+kU+kS<0)   // => "Normal anti-nucleus" case
  //{
  //  G4int b=(kS+kS-kD-kU)/3;
  //  G4int d=kD-kU;
  //  G4int n=(b-d)/2;
  //  return -90000000-1000*(1000*kS+n+d)-n;   // @@ double notation for anti-nuclei
  //}
  //else
  //{
    G4int b=(kU+kD-kS-kS)/3;                 // Baryon number-n*{LAMBDA=kS)
    if(!b && !kS) return GetSPDGCode();      // Not a nucleus
    G4int d=kU-kD;                           // Isotopic shift
    G4int n=(b-d)/2;                         // A#of neutrons
    return 90000000+1000*(1000*kS+n+d)+n;
  //}
  //return 0;
}

// Convert particle to anti-particle
inline void G4QContent::Anti()
{
  G4int r=nD;
  nD = nAD;
  nAD= r;
  r  = nU;
  nU = nAU;
  nAU= r;
  r  = nS;
  nS = nAS;
  nAS= r;
}

// Add Quark Content
inline G4QContent G4QContent::operator+=(const G4QContent& rhs)
{
  nD += rhs.nD;
  nU += rhs.nU;
  nS += rhs.nS;
  nAD+= rhs.nAD;
  nAU+= rhs.nAU;
  nAS+= rhs.nAS;
  return *this;
} 

// Add Quark Content
inline G4QContent G4QContent::operator+=(G4QContent& rhs)
{
  nD += rhs.nD;
  nU += rhs.nU;
  nS += rhs.nS;
  nAD+= rhs.nAD;
  nAU+= rhs.nAU;
  nAS+= rhs.nAS;
  return *this;
} 

// Multiply Quark Content by integer number
inline G4QContent G4QContent::operator*=(const G4int& rhs)
{
  nU *= rhs;
  nD *= rhs;
  nS *= rhs;
  nAU*= rhs;
  nAD*= rhs;
  nAS*= rhs;
  return *this;
} 

// Multiply Quark Content by integer number
inline G4QContent G4QContent::operator*=(G4int& rhs)
{
  nU *= rhs;
  nD *= rhs;
  nS *= rhs;
  nAU*= rhs;
  nAD*= rhs;
  nAS*= rhs;
  return *this;
} 

inline void  G4QContent::SetU(G4int n) {nU=n;}
inline void  G4QContent::SetD(G4int n) {nD=n;}
inline void  G4QContent::SetS(G4int n) {nS=n;}
inline void  G4QContent::SetAU(G4int n){nAU=n;}
inline void  G4QContent::SetAD(G4int n){nAD=n;}
inline void  G4QContent::SetAS(G4int n){nAS=n;}

inline void  G4QContent::IncU(G4int n) {nU+=n;}
inline void  G4QContent::IncD(G4int n) {nD+=n;}
inline void  G4QContent::IncS(G4int n) {nS+=n;}
inline void  G4QContent::IncAU(G4int n){nAU+=n;}
inline void  G4QContent::IncAD(G4int n){nAD+=n;}
inline void  G4QContent::IncAS(G4int n){nAS+=n;}

inline void  G4QContent::DecU(G4int n) {nU-=n;}
inline void  G4QContent::DecD(G4int n) {nD-=n;}
inline void  G4QContent::DecS(G4int n) {nS-=n;}
inline void  G4QContent::DecAU(G4int n){nAU-=n;}
inline void  G4QContent::DecAD(G4int n){nAD-=n;}
inline void  G4QContent::DecAS(G4int n){nAS-=n;}

// Private member functions
inline G4QContent G4QContent::GetThis()const{return G4QContent(nD,nU,nS,nAD,nAU,nAS);}
#endif
