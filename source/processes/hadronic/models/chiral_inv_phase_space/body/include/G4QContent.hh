// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QContent.hh,v 1.1 1999-11-17 11:04:14 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4QContent_h
#define G4QContent_h 1

// ----------------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QContent ----------------
//             by Mikhail Kossov, Sept 1999.
//  class header for Quasmon initiated Contents used by the CHIPS Model
// ----------------------------------------------------------------------

#include <iostream.h>
#include "globals.hh"
#include "Randomize.hh"

class G4QContent
{
public:
  G4QContent(G4int u=0, G4int au=0, G4int d=0, G4int ad=0, G4int s=0, G4int as=0);
  G4QContent(const G4QContent& rhs);

  ~G4QContent();

  const G4QContent& operator=(const G4QContent& rhs);
  int        operator==(const G4QContent& rhs) const;
  int        operator!=(const G4QContent& rhs) const;
  G4QContent operator+=(const G4QContent& rhs);
  G4QContent operator-=(const G4QContent& rhs);
  G4QContent operator*=(const G4int& rhs);
  void           IncQAQ(const G4int& nQAQ=1);
  G4int          DecQAQ(const G4int& nQAQ=1);
  
  const G4int  GetCharge() const;
  const G4int  GetStrangeness() const;
  const G4int  NOfCombinations(const G4QContent& rhs); //@@ can be "operator/"
  const G4int  GetTot();
  const G4bool CheckNegative();

  const G4int& GetU();      
  const G4int& GetD(); 
  const G4int& GetS(); 
  const G4int& GetAU();  
  const G4int& GetAD();
  const G4int& GetAS();      

  G4int GetUU();      
  G4int GetDD(); 
  G4int GetSS(); 
  G4int GetUD();      
  G4int GetUS(); 
  G4int GetDS(); 
  G4int GetAUAU();  
  G4int GetADAD();
  G4int GetASAS();      
  G4int GetAUAD();  
  G4int GetAUAS();
  G4int GetADAS();      

  void SetU(G4int n=0);
  void SetD(G4int n=0);
  void SetS(G4int n=0);
  void SetAU(G4int n=0);
  void SetAD(G4int n=0);
  void SetAS(G4int n=0);

  void IncU(G4int n=1);
  void IncD(G4int n=1);
  void IncS(G4int n=1);
  void IncAU(G4int n=1);
  void IncAD(G4int n=1);
  void IncAS(G4int n=1);

  void DecU(G4int n=1);
  void DecD(G4int n=1);
  void DecS(G4int n=1);
  void DecAU(G4int n=1);
  void DecAD(G4int n=1);
  void DecAS(G4int n=1);

private:            //                       C    S
  G4int nU;         // a#of      u-quarks (+2/3)( 0)
  G4int nD;         // a#of      d-quarks (-1/3)( 0)
  G4int nS;         // a#of      s-quarks (-1/3)(+1)
  G4int nAU;        // a#of anti-u-quarks (-2/3)( 0)
  G4int nAD;        // a#of anti-d-quarks (+1/3)( 0)
  G4int nAS;        // a#of anti-s-quarks (+1/3)(-1)
};

ostream& operator<<(ostream& lhs, G4QContent& rhs);

inline int G4QContent::operator==(const G4QContent& rhs) const {return this==&rhs;}	
inline int G4QContent::operator!=(const G4QContent& rhs) const {return this!=&rhs;}
inline const G4int  G4QContent::GetTot() {return nU+nD+nS+nAU+nAD+nAS;}
inline const G4bool G4QContent::CheckNegative(){return nU<0||nD<0||nS<0||nAU<0||nAD<0||nAS<0;}

inline const G4int& G4QContent::GetU(){return nU;}
inline const G4int& G4QContent::GetD(){return nD;}
inline const G4int& G4QContent::GetS(){return nS;}
inline const G4int& G4QContent::GetAU(){return nAU;}
inline const G4int& G4QContent::GetAD(){return nAD;}
inline const G4int& G4QContent::GetAS(){return nAS;}

inline G4int G4QContent::GetUU(){return nU*(nU-1)/2;}
inline G4int G4QContent::GetDD(){return nD*(nD-1)/2;}
inline G4int G4QContent::GetSS(){return nS*(nS-1)/2;}
inline G4int G4QContent::GetUD(){return nU*nD;}
inline G4int G4QContent::GetUS(){return nU*nS;}
inline G4int G4QContent::GetDS(){return nD*nS;}
inline G4int G4QContent::GetAUAU(){return nAU*(nAU-1)/2;}
inline G4int G4QContent::GetADAD(){return nAD*(nAD-1)/2;}
inline G4int G4QContent::GetASAS(){return nAS*(nAS-1)/2;}
inline G4int G4QContent::GetAUAD(){return nAU*nAD;}
inline G4int G4QContent::GetAUAS(){return nAU*nAS;}
inline G4int G4QContent::GetADAS(){return nAD*nAS;}

inline G4QContent G4QContent::operator+=(const G4QContent& rhs)
{
  nU += rhs.nU;
  nD += rhs.nD;
  nS += rhs.nS;
  nAU+= rhs.nAU;
  nAD+= rhs.nAD;
  nAS+= rhs.nAS;
  return *this;
} 

inline G4QContent G4QContent::operator-=(const G4QContent& rhs)
{
  nU -= rhs.nU;
  if (nU<0)
  {
    nAU -= nU;
    nU   = 0;
  }
  nD -= rhs.nD;
  if (nD<0)
  {
    nAD -= nD;
    nD   = 0;
  }
  nS -= rhs.nS;
  if (nS<0)
  {
    nAS -= nS;
    nS   = 0;
  }
  nAU -= rhs.nAU;
  if (nAU<0)
  {
    nU -= nAU;
    nAU = 0;
  }
  nAD -= rhs.nAD;
  if (nAD<0)
  {
    nD -= nAD;
    nAD = 0;
  }
  nAS -= rhs.nAS;
  if (nAS<0)
  {
    nS -= nAS;
    nAS = 0;
  }
  return *this;
} 
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

#endif







