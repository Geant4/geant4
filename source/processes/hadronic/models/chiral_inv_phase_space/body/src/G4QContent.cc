// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QContent.cc,v 1.2 1999-12-15 14:52:10 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QContent ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Quasmon initiated Contents used by CHIPS Model
// ------------------------------------------------------------------

//#define debug

#include "G4QContent.hh"

G4std::ostream& operator<<(G4std::ostream& lhs, G4QContent& rhs)
{
  lhs << "{" << rhs.GetU() << "," << rhs.GetD() << "," << rhs.GetS() << ","
      << rhs.GetAU() << "," << rhs.GetAD() << "," << rhs.GetAS() << "}";
  return lhs;
}

G4QContent::G4QContent(G4int u, G4int au, G4int d, G4int ad, G4int s, G4int as)
{
  nU      = u;
  nD      = d;
  nS      = s;
  nAU     = au;
  nAD     = ad;
  nAS     = as;
}

G4QContent::G4QContent(const G4QContent &right)
{
  nU                  = right.nU;
  nD                  = right.nD;
  nS                  = right.nS;
  nAU                 = right.nAU;
  nAD                 = right.nAD;
  nAS                 = right.nAS;
}

const G4QContent & G4QContent::operator=(const G4QContent &right)
{
  nU                  = right.nU;
  nD                  = right.nD;
  nS                  = right.nS;
  nAU                 = right.nAU;
  nAD                 = right.nAD;
  nAS                 = right.nAS;
		
  return *this;
}

G4QContent::~G4QContent() {}

G4int G4QContent::DecQAQ(const G4int& nQAQ)
{
#ifdef debug
  cout<<"DecQC:in U="<<nU<<",D="<<nD<<",S="<<nS<<",AU="<<nAU<<",AD="<<nAD<<",AS="<<nAS<<G4endl;
#endif
  G4int nUP=0;             // U/AU min factor
  if (nU>=nAU) nUP+=nAU;
  else         nUP+= nU;

  G4int nDP=0;             // D/AD min factor
  if (nD>=nAD) nDP+=nAD;
  else         nDP+= nD;

  G4int nSP=0;             // S/AS min factor
  if (nS>=nAS) nSP+=nAS;
  else         nSP+= nS;

  G4int nLP  =nUP+nDP;     // a#of light quark pairs
  G4int nTotP=nLP+nSP;     // total#of existing pairs
  G4int nReal=nQAQ;        // demanded #of pairs for reduction
  G4int nRet =nTotP-nQAQ;  // full job is not possible
  if (nRet<0) nReal=nTotP; // make the task easier
  if (!nReal) return nRet; // Now nothing to do
  // ---------- Decrimenting by nReal pairs
#ifdef debug
  cout << "DecQC: demanded "<<nQAQ<<" pairs, executed "<<nReal<<" pairs"<<G4endl;
#endif
  for (int i=0; i<nReal; i++)
  {
    int j = nTotP*G4UniformRand();
    if (nUP && j<nUP)      // U-antiU pair
	{
#ifdef debug
      cout << "DecQC: decrementing UAU pair UP="<<nUP<<",nU="<<nU<<",nAU="<<nAU<<G4endl;
#endif
      nU--;
      nAU--;
      nUP--;
      nLP--;
      nTotP--;
	} 
    else if (nDP && j<nLP) // D-antiD pair
	{
#ifdef debug
      cout << "DecQC: decrementing DAD pair DP="<<nDP<<",nD="<<nD<<",nAD="<<nAD<<G4endl;
#endif
      nD--;
      nAD--;
      nDP--;
      nLP--;
      nTotP--;
	} 
    else if (nSP)          // S-antiS pair
	{
#ifdef debug
      cout << "DecQC: decrementing SAS pair SP="<<nSP<<",nS="<<nS<<",nAS="<<nAS<<G4endl;
#endif
      nS--;
      nAS--;
      nSP--;
      nTotP--;
	}
    else G4cerr<<"***DecQC:i="<<i<<",j="<<j<<",UP="<<nUP<<",DP="<<nDP<<",T="<<nTotP<<G4endl;
  }
#ifdef debug
  cout<<"DecQC:out U="<<nU<<",D="<<nD<<",S="<<nS<<",AU="<<nAU<<",AD="<<nAD<<",AS="<<nAS<<G4endl;
#endif
  return nRet;
}

void G4QContent::IncQAQ(const G4int& nQAQ)
{
  for (int i=0; i<nQAQ; i++)
  {
    int j = 3*G4UniformRand();
    if (j = 0)
    {
      nU++;
      nAU++;
	}
    else if (j = 1)
	{
      nD++;
      nAD++;
	}
    else
	{
      nS++;
      nAS++;
    }      
  }
}

const G4int G4QContent::GetCharge() const
{
  static const G4int cU  = 2;
  static const G4int cD  =-1;
  static const G4int cS  =-1;
  static const G4int cAU =-2;
  static const G4int cAD = 1;
  static const G4int cAS = 1;

  G4int c=0;
  if(nU) c+=nU*cU;
  if(nD) c+=nD*cD;
  if(nS) c+=nS*cS;
  if(nAU)c+=nAU*cAU;
  if(nAD)c+=nAD*cAD;
  if(nAS)c+=nAS*cAS;

  if(c%3) G4cerr << "***G4Content: Charge="<<c<<"/3 is not an integer value" << G4endl;

  return c/3;
}

const G4int G4QContent::GetStrangeness() const
{
  G4int  s =0;
  if(nS) s+=nS;
  if(nAS)s-=nAS;

  return s;
}

// === Calculate a number of combinations of rhc out of lhc ===
const G4int G4QContent::NOfCombinations(const G4QContent& rhs)
{
  G4int c=1;

  if(rhs.nU)
  {
    int j=nU;
    if (j<=0) return 0;
    for (int i=1; i<=rhs.nU; i++)
	{
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  };

  if(rhs.nD)
  {
    int j=nD;
    if (j<=0) return 0;
    for (int i=1; i<=rhs.nD; i++)
	{
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  };

  if(rhs.nS)
  {
    int j=nS;
    if (j<=0) return 0;
    for (int i=1; i<=rhs.nS; i++)
	{
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  };

 if(rhs.nAU)
  {
    int j=nAU;
    if (j<=0) return 0;
    for (int i=1; i<=rhs.nAU; i++)
	{
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  };

 if(rhs.nAD)
  {
    int j=nAD;
    if (j<=0) return 0;
    for (int i=1; i<=rhs.nAD; i++)
	{
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  };

 if(rhs.nAS)
  {
    int j=nAS;
    if (j<=0) return 0;
    for (int i=1; i<=rhs.nAS; i++)
	{
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  };
  
  return c;
} 
