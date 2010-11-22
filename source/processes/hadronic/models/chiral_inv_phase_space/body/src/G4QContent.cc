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
// $Id: G4QContent.cc,v 1.50 2010-11-22 07:08:01 dennis Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QContent ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Quasmon initiated Contents used by CHIPS Model
// --------------------------------------------------------------------
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

//#define debug
//#define pdebug
//#define ppdebug
//#define erdebug

#include "G4QContent.hh"
#include <cmath>
using namespace std;

// Constructors

// Initialize by  the full quark content
G4QContent::G4QContent(G4int d, G4int u, G4int s, G4int ad, G4int au, G4int as):
  nD(d),nU(u),nS(s),nAD(ad),nAU(au),nAS(as)
{
  // Bug report 1131 identifies conditional to have no effect.
  // Remove it. 
  // D.H. Wright 11 November 2010 
  /*
  if(d<0||u<0||s<0||ad<0||au<0||as<0)
  {
#ifdef erdebug
    G4cerr<<"***G4QContent:"<<d<<","<<u<<","<<s<<","<<ad<<","<<au<<","<<as<<G4endl;
#endif
    if(d<0) ad-=d;
    if(u<0) au-=u;
    if(s<0) as-=s;
    if(ad<0) d-=ad;
    if(au<0) u-=au;
    if(as<0) s-=as;
  }
  */
}

// Initialize by a pair of partons
G4QContent::G4QContent(std::pair<G4int,G4int> PP): nD(0),nU(0),nS(0),nAD(0),nAU(0),nAS(0)
{
  G4int P1=PP.first;
  G4int P2=PP.second;
  if(!P1 || !P2)
  {
    G4cerr<<"***G4QContent::Constr(pair): Zero parton P1="<<P1<<", P2="<<P2<<G4endl;
    G4Exception("G4QContent::Constructor(pair):","72",FatalException,"WrongPartonPair");
  }
#ifdef pdebug
  G4cout<<"G4QContent::PairConstr: P1="<<P1<<", P2="<<P2<<G4endl;
#endif
  G4bool suc=true;
  G4int A1=P1;
  if     (P1 > 7) A1=  P1/100;
  else if(P1 <-7) A1=(-P1)/100;
  else if(P1 < 0) A1= -P1;
  G4int P11=0;
  G4int P12=0;
  if(A1>7)
  {
    P11=A1/10;
    P12=A1%10;
  }
  if(P1>0)
  {
    if(!P11)
    {
      if     (A1==1) ++nD; 
      else if(A1==2) ++nU; 
      else if(A1==3) ++nS;
      else           suc=false;
    }
    else
    {
      if     (P11==1) ++nD; 
      else if(P11==2) ++nU; 
      else if(P11==3) ++nS;
      else            suc=false;
      if     (P12==1) ++nD; 
      else if(P12==2) ++nU; 
      else if(P12==3) ++nS;
      else            suc=false;
    }
  }
  else // negative parton
  {
    if(!P11)
    {
      if     (A1==1) ++nAD; 
      else if(A1==2) ++nAU; 
      else if(A1==3) ++nAS;
      else           suc=false;
    }
    else
    {
      if     (P11==1) ++nAD; 
      else if(P11==2) ++nAU; 
      else if(P11==3) ++nAS;
      else            suc=false;
      if     (P12==1) ++nAD; 
      else if(P12==2) ++nAU; 
      else if(P12==3) ++nAS;
      else            suc=false;
    }
  }
#ifdef pdebug
  G4cout<<"G4QContent::PCo:1:"<<nD<<","<<nU<<","<<nS<<","<<nAD<<","<<nAU<<","<<nAS<<G4endl;
#endif
  G4int A2=P2;
  if     (P2 > 7) A2=  P2/100;
  else if(P2 <-7) A2=(-P2)/100;
  else if(P2 < 0) A2= -P2;
  G4int P21=0;
  G4int P22=0;
  if(A2>7)
  {
    P21=A2/10;
    P22=A2%10;
  }
  if(P2>0)
  {
    if(!P21)
    {
      if     (A2==1) ++nD; 
      else if(A2==2) ++nU; 
      else if(A2==3) ++nS;
      else           suc=false;
    }
    else
    {
      if     (P21==1) ++nD; 
      else if(P21==2) ++nU; 
      else if(P21==3) ++nS;
      else            suc=false;
      if     (P22==1) ++nD; 
      else if(P22==2) ++nU; 
      else if(P22==3) ++nS;
      else            suc=false;
    }
  }
  else // negative parton
  {
    if(!P21)
    {
      if     (A2==1) ++nAD; 
      else if(A2==2) ++nAU; 
      else if(A2==3) ++nAS;
      else           suc=false;
    }
    else
    {
      if     (P21==1) ++nAD; 
      else if(P21==2) ++nAU; 
      else if(P21==3) ++nAS;
      else            suc=false;
      if     (P22==1) ++nAD; 
      else if(P22==2) ++nAU; 
      else if(P22==3) ++nAS;
      else            suc=false;
    }
  }
  if(!suc)
  {
    G4cerr<<"***G4QContent::Constr(pair): Impossible partons, P1="<<P1<<",P2="<<P2<<G4endl;
    G4Exception("G4QContent::Constructor(pair):","72",FatalException,"ImpossibPartonPair");
  }
#ifdef pdebug
  G4cout<<"G4QContent::PCo:2:"<<nD<<","<<nU<<","<<nS<<","<<nAD<<","<<nAU<<","<<nAS<<G4endl;
#endif
}

G4QContent::G4QContent(const G4QContent& right)
{
  nU  = right.nU;
  nD  = right.nD;
  nS  = right.nS;
  nAU = right.nAU;
  nAD = right.nAD;
  nAS = right.nAS;
}

G4QContent::G4QContent(G4QContent* right)
{
  nU  = right->nU;
  nD  = right->nD;
  nS  = right->nS;
  nAU = right->nAU;
  nAD = right->nAD;
  nAS = right->nAS;
}

// Assignment operator (copy stile for possible Vector extention)
const G4QContent& G4QContent::operator=(const G4QContent &right)
{//               ==============================================
  if(this != &right)                          // Beware of self assignment
  {
    nU  = right.nU;
    nD  = right.nD;
    nS  = right.nS;
    nAU = right.nAU;
    nAD = right.nAD;
    nAS = right.nAS;
  }
  return *this;
}

// Standard output for QC {d,u,s,ad,au,as}
ostream& operator<<(ostream& lhs, G4QContent& rhs)
{//      =========================================
  lhs << "{" << rhs.GetD() << "," << rhs.GetU() << "," << rhs.GetS() << ","
      << rhs.GetAD() << "," << rhs.GetAU() << "," << rhs.GetAS() << "}";
  return lhs;
}

// Standard output for const QC {d,u,s,ad,au,as}
ostream& operator<<(ostream& lhs, const G4QContent& rhs)
{//      ===============================================
  lhs << "{" << rhs.GetD() << "," << rhs.GetU() << "," << rhs.GetS() << ","
      << rhs.GetAD() << "," << rhs.GetAU() << "," << rhs.GetAS() << "}";
  return lhs;
}

// Subtract Quark Content
G4QContent G4QContent::operator-=(const G4QContent& rhs)
//         =============================================
{
#ifdef debug
  G4cout<<"G4QC::-=(const): is called:"<<G4endl;
#endif
  G4int rD=rhs.nD;
  G4int rU=rhs.nU;
  G4int rS=rhs.nS;
  G4int rAD=rhs.nAD;
  G4int rAU=rhs.nAU;
  G4int rAS=rhs.nAS;
  ///////////G4int rQ =rD+rU+rS;
  ///////////G4int rAQ=rAD+rAU+rAS;
  ///////////G4int nQ =nD+nU+nS;
  ///////////G4int nAQ=nAD+nAU+nAS;
  if(nU<rU||nAU<rAU||nD<rD||nAD<rAD)
  {
    G4int dU=rU-nU;
    G4int dAU=rAU-nAU;
    if(dU>0||dAU>0)
    {
      G4int kU=dU;
      if(kU<dAU) kU=dAU;                  // Get biggest difference
      G4int mU=rU;
      if(rAU<mU) mU=rAU;                  // Get a#of possible SS pairs
      if(kU<=mU)                          // Total compensation
      {
        rU-=kU;
        rAU-=kU;
        rD+=kU;
        rAD+=kU;
      }
      else                                // Partial compensation
      {
        rU-=mU;
        rAU-=mU;
        rD+=mU;
        rAD+=mU;
      }
    }
    G4int dD=rD-nD;
    G4int dAD=rAD-nAD;
    if(dD>0||dAD>0)
    {
      G4int kD=dD;
      if(kD<dAD) kD=dAD;                  // Get biggest difference
      G4int mD=rD;
      if(rAD<mD) mD=rAD;                  // Get a#of possible SS pairs
      if(kD<=mD)                          // Total compensation
      {
        rD-=kD;
        rAD-=kD;
        rU+=kD;
        rAU+=kD;
      }
      else                                // Partial compensation
      {
        rD-=mD;
        rAD-=mD;
        rU+=mD;
        rAU+=mD;
      }
    }
  }
#ifdef debug
  G4cout<<"G4QC::-=:comp: "<<rD<<","<<rU<<","<<rS<<","<<rAD<<","<<rAU<<","<<rAS<<G4endl;
#endif
  if(rS==1 && rAS==1 && (nS<1 || nAS<1))  // Eta case, switch quark pairs (?)
  {
    rS =0;
    rAS=0;
    if(nU>rU&&nAU>rAU)
    {
      rU +=1;
      rAU+=1;
    }
    else
    {
      rD +=1;
      rAD+=1;
    }
  }
  nD -= rD;
  if (nD<0)
  {
    nAD -= nD;
    nD   = 0;
  }
  nU -= rU;
  if (nU<0)
  {
    nAU -= nU;
    nU   = 0;
  }
  nS -= rS;
  if (nS<0)
  {
    nAS -= nS;
    nS   = 0;
  }
  nAD -= rAD;
  if (nAD<0)
  {
    nD -= nAD;
    nAD = 0;
  }
  nAU -= rAU;
  if (nAU<0)
  {
    nU -= nAU;
    nAU = 0;
  }
  nAS -= rAS;
  if (nAS<0)
  {
    nS -= nAS;
    nAS = 0;
  }
  return *this;
} 

// Subtract Quark Content
G4QContent G4QContent::operator-=(G4QContent& rhs)
//         =======================================
{
#ifdef debug
  G4cout<<"G4QC::-=: is called:"<<G4endl;
#endif
  G4int rD=rhs.nD;
  G4int rU=rhs.nU;
  G4int rS=rhs.nS;
  G4int rAD=rhs.nAD;
  G4int rAU=rhs.nAU;
  G4int rAS=rhs.nAS;
  G4int rQ =rD+rU+rS;
  G4int rAQ=rAD+rAU+rAS;
  G4int nQ =nD+nU+nS;
  G4int nAQ=nAD+nAU+nAS;
  if(nQ<rQ||nAQ<rAQ)
  {
    G4int dU=rU-nU;
    G4int dAU=rAU-nAU;
    if(dU>0||dAU>0)
    {
      G4int kU=dU;
      if(kU<dAU) kU=dAU;                  // Get biggest difference
      G4int mU=rU;
      if(rAU<mU) mU=rAU;                  // Get a#of possible SS pairs
      if(kU<=mU)                          // Total compensation
      {
        rU-=kU;
        rAU-=kU;
        rD+=kU;
        rAD+=kU;
      }
      else                                // Partial compensation
      {
        rU-=mU;
        rAU-=mU;
        rD+=mU;
        rAD+=mU;
      }
    }
    G4int dD=rD-nD;
    G4int dAD=rAD-nAD;
    if(dD>0||dAD>0)
    {
      G4int kD=dD;
      if(kD<dAD) kD=dAD;                  // Get biggest difference
      G4int mD=rD;
      if(rAD<mD) mD=rAD;                  // Get a#of possible SS pairs
      if(kD<=mD)                          // Total compensation
      {
        rD-=kD;
        rAD-=kD;
        rU+=kD;
        rAU+=kD;
      }
      else                                // Partial compensation
      {
        rD-=mD;
        rAD-=mD;
        rU+=mD;
        rAU+=mD;
      }
    }
  }
  if(rS==1 && rAS==1 && (nS<1 || nAS<1))  // Eta case, switch quark pairs (?)
  {
    rS =0;
    rAS=0;
    if(nU>rU&&nAU>rAU)
    {
      rU +=1;
      rAU+=1;
    }
    else
    {
      rD +=1;
      rAD+=1;
    }
  }
  nD -= rD;
  if (nD<0)
  {
    nAD -= nD;
    nD   = 0;
  }
  nU -= rU;
  if (nU<0)
  {
    nAU -= nU;
    nU   = 0;
  }
  nS -= rS;
  if (nS<0)
  {
    nAS -= nS;
    nS   = 0;
  }
  nAD -= rAD;
  if (nAD<0)
  {
    nD -= nAD;
    nAD = 0;
  }
  nAU -= rAU;
  if (nAU<0)
  {
    nU -= nAU;
    nAU = 0;
  }
  nAS -= rAS;
  if (nAS<0)
  {
    nS -= nAS;
    nAS = 0;
  }
  return *this;
} 

// Overloading of QC addition
G4QContent operator+(const G4QContent& lhs, const G4QContent& rhs)
{//        =======================================================
  G4QContent s  = lhs;
  return     s += rhs;
}

// Overloading of QC subtraction
G4QContent operator-(const G4QContent& lhs, const G4QContent& rhs)
{//        =======================================================
  G4QContent s  = lhs;
  return     s -= rhs;
}

// Overloading of QC multiplication by Int
G4QContent operator*(const G4QContent& lhs, const G4int&      rhs)
{//        =======================================================
  G4QContent s  = lhs;
  return     s *= rhs;
}

// Overloading of Int times QC multiplication
G4QContent operator*(const G4int&      lhs, const G4QContent& rhs)
{//        =======================================================
  G4QContent s  = rhs;
  return     s *= lhs;
}

// Destructor - - - - - - - 
G4QContent::~G4QContent() {}

// Subtract neutral pion from Quark Content (with possible hidden strangeness)
G4bool G4QContent::SubtractPi0()
{//    =========================
#ifdef debug
  G4cout<<"G4QC::SubtractPi0: U="<<nU<<", AU="<<nAU<<", D="<<nD<<", AD="<<nAD<<G4endl;
#endif
  G4int tot=GetTot();
  G4int ab =GetBaryonNumber();
  if(ab){if(tot<3*ab+2) return false;}
  else if(tot<4) return false;

  if(nU>0 && nAU>0)
  {
    nU--;
    nAU--;
    return true;
  }
  else if(nD>0 && nAD>0)
  {
    nD--;
    nAD--;
    return true;
  }
  return false;
}

// Subtract charged pion from Quark Content
G4bool G4QContent::SubtractPion()
{//    ==========================
#ifdef debug
  G4cout<<"G4QC::SubtractPion: U="<<nU<<", AU="<<nAU<<", D="<<nD<<", AD="<<nAD<<G4endl;
#endif
  G4int tot=GetTot();
  G4int ab =GetBaryonNumber();
  if(ab){if(tot<3*ab+2) return false;}
  else if(tot<4) return false;

  if((nU>nAU || (ab && nU>0))&& nAD>0)
  {
    nU--;
    nAD--;
    return true;
  }
  else if((nAU>nU || (ab && nAU>0)) && nD>0)
  {
    nAU--;
    nD--;
    return true;
  }
  return false;
}

// Subtract Hadron from Quark Content
G4bool G4QContent::SubtractHadron(G4QContent h)
{//    ========================================
#ifdef debug
  G4cout<<"G4QC::SubtractHadron "<<h<<" is called for QC="<<GetThis()<<G4endl;
#endif
  if(h.GetU()<=nU  && h.GetD()<=nD  && h.GetS()<=nS&&
     h.GetAU()<=nAU&& h.GetAD()<=nAD&& h.GetAS()<=nAS) return true;
  return false;
}

// Subtract Kaon from Quark Content
G4bool G4QContent::SubtractKaon(G4double mQ)
{//    =====================================
#ifdef debug
  G4cout<<"G4QC::SubtractKaon is called: QC="<<GetThis()<<G4endl;
#endif
  if(mQ<640.) return false;
  G4int tot=GetTot();
  G4int ab =GetBaryonNumber();
  if(ab){if(tot<3*ab+2) return false;}
  else if(tot<4) return false;

  if((nS>nAS || (ab && nS>0)) && (nAD>0 || nAU>0))
  {
    nS--;
    if (nAU>0) nAU--;
    else       nAD--;
    return true;
  }
  else if((nAS>nS || (ab && nAS>0)) && (nD>0 || nU>0))
  {
    nAS--;
    if (nU>0) nU--;
    else      nD--;
    return true;
  }
#ifdef debug
  G4cout<<"G4QCont::SubtractKaon Can't SubtractKaon: QC="<<GetThis()<<G4endl;
#endif
  return false;
}

// Split any hadronic system in two hadrons
G4QContent G4QContent::SplitChipo (G4double mQ)
{//        ====================================
  G4QContent Pi(0,1,0,1,0,0);
  if      (nU>0&&nAU>0) Pi=G4QContent(0,1,0,0,1,0);
  else if (nD>0&&nAD>0) Pi=G4QContent(1,0,0,1,0,0);
  else if (nD>=nU&&nAU>=nAD) Pi=G4QContent(1,0,0,0,1,0);
  G4int bn=GetBaryonNumber();
  G4int b =abs(bn);
  if(!b && mQ<545.&&nS>0&&nAS>0) // Cancel strange sea
  {
    G4int      ss= nS;
    if(nAS<nS) ss=nAS;
    nS -=ss;
    nAS-=ss;
  }
  if       (!b)DecQAQ(-4);
  else if(b==1)DecQAQ(-5);
  else         DecQAQ(0);
  G4int tot=GetTot();
  G4int q=GetQ();
  G4int aq=GetAQ();
  G4QContent r=Pi;     // Pion prototype of returned value
  if((tot!=4||q!=2) && (tot!=5||(q!=1&&aq!=1)) && (tot!=6||abs(b)!=2))
  {
    //#ifdef erdebug
    G4cerr<<"***G4QCont::SplitChipo: QC="<<GetThis()<<G4endl;
    //#endif
  }
  else if(tot==4)      // Mesonic (eight possibilities)
  {
    r=GetThis();
    if     (r.SubtractPi0())    ; // Try any trivial algorithm of splitting
    else if(r.SubtractPion())   ;
    else if(r.SubtractKaon(mQ)) ;
    else
    {
      //#ifdef debug
      G4cerr<<"***G4QCont::SplitChipo:MesTot="<<tot<<",b="<<b<<",q="<<q<<",a="<<aq<<G4endl;
      //#endif
    }
  }
  else if(b==1&&tot==5)      // Baryonic (four possibilities)
  {
    if(nU==3)
    {
      r.SetU(1);
      r+=IndAQ();
    }
    else if(nD==3)
    {
      r.SetD(1);
      r+=IndAQ();
    }
    else if(nS==3)
    {
      r.SetS(1);
      r+=IndAQ();
    }
    else if(nAU==3)
    {
      r.SetAU(1);
      r+=IndQ();
    }
    else if(nAD==3)
    {
      r.SetAD(1);
      r+=IndQ();
    }
    else if(nAS==3)
    {
      r.SetAS(1);
      r+=IndQ();
    }
    else if(q==1&&nU)
    {
      r.SetU(1);
      if(nAU) r.SetAU(1);
      else    r.SetAD(1);
    }
    else if(q==1&&nD)
    {
      r.SetD(1);
      if(nAD) r.SetAD(1);
      else    r.SetAU(1);
    }
    else if(q==1&&nS)
    {
      r.SetS(1);
      if(nAS) r.SetAS(1);
      else    r.SetAU(1);
    }
    else if(aq==1&&nAU)
    {
      r.SetAU(1);
      if(nU) r.SetU(1);
      else   r.SetD(1);
    }
    else if(aq==1&&nAD)
    {
      r.SetAD(1);
      if(nD) r.SetD(1);
      else   r.SetU(1);
    }
    else if(aq==1&&nAS)
    {
      r.SetAS(1);
      if(nS) r.SetS(1);
      else   r.SetU(1);
    }
    else
    {
      //#ifdef erdebug
      G4cerr<<"***G4QCont::SplitChipo: Baryonic tot=5,b=1,qCont="<<GetThis()<<G4endl;
      //#endif
    }
  }
  else if(tot==b*3)          // MultyBaryon cace
  {
    r=GetThis();
    if (bn>0)                // baryonium
    {
      G4QContent la(1,1,1,0,0,0);
      G4QContent nt(2,1,0,0,0,0);
      G4QContent pr(1,2,0,0,0,0);
      G4QContent ks(0,1,2,0,0,0);
      if (nD>nU) ks=G4QContent(1,0,2,0,0,0);
      G4QContent dm(3,0,0,0,0,0);
      G4QContent dp(0,3,0,0,0,0);
      G4QContent om(0,0,3,0,0,0);
      if     (nU>=nD&&nU>=nS)
      {
        if     (r.SubtractHadron(pr)) r-=pr; // These functions only check
        else if(r.SubtractHadron(dp)) r-=dp;
        else if(r.SubtractHadron(nt)) r-=nt;
        else if(r.SubtractHadron(la)) r-=la;
        else if(r.SubtractHadron(dm)) r-=dm;
        else
        {
          //#ifdef erdebug
          G4cerr<<"***G4QCont::SplitChipo:Dibar (1) tot=6, b=2, qCont="<<GetThis()<<G4endl;
          //#endif
        }
      }
      else if(nD>=nU&&nD>=nS)
      {
        if     (r.SubtractHadron(nt)) r-=nt; // These functions only check
        else if(r.SubtractHadron(dm)) r-=dm;
        else if(r.SubtractHadron(pr)) r-=pr;
        else if(r.SubtractHadron(dp)) r-=dp;
        else if(r.SubtractHadron(la)) r-=la;
        else
        {
          //#ifdef erdebug
          G4cerr<<"***G4QContent::SplitChipo:Dib(2) tot=6, b=2, qCont="<<GetThis()<<G4endl;
          //#endif
        }
      }
      else
      {
        if     (r.SubtractHadron(la)) r-=la; // These functions only check
        else if(r.SubtractHadron(ks)) r-=ks;
        else if(r.SubtractHadron(om)) r-=om;
        else if(r.SubtractHadron(pr)) r-=pr;
        else if(r.SubtractHadron(nt)) r-=nt;
        else
        {
          //#ifdef erdebug
          G4cerr<<"***G4QContent::SplitChipo:Dib(3) tot=6, b=2, qCont="<<GetThis()<<G4endl;
          //#endif
        }
      }
    }
    else                     // Anti-baryonium
    {
      G4QContent la(0,0,0,1,1,1);
      G4QContent pr(0,0,0,1,2,0);
      G4QContent nt(0,0,0,2,1,0);
      G4QContent ks(0,1,2,0,0,0);
      if(nAD>nAU)ks=G4QContent(0,0,0,1,0,2);
      G4QContent dm(0,0,0,3,0,0);
      G4QContent dp(0,0,0,0,3,0);
      G4QContent om(0,0,0,0,0,3);
      if     (nAU>=nAD&&nAU>=nAS)
      {
        if     (r.SubtractHadron(pr)) r-=pr; // These functions only check
        else if(r.SubtractHadron(dp)) r-=dp;
        else if(r.SubtractHadron(nt)) r-=nt;
        else if(r.SubtractHadron(la)) r-=la;
        else if(r.SubtractHadron(dm)) r-=dm;
        else
        {
          //#ifdef erdebug
          G4cerr<<"***G4QContent::SplitChipo:ADib(1) tot=6,b=2, qCont="<<GetThis()<<G4endl;
          //#endif
        }
      }
      else if(nAD>=nAU&&nAD>=nAS)
      {
        if     (r.SubtractHadron(nt)) r-=nt; // These functions only check
        else if(r.SubtractHadron(dm)) r-=dm;
        else if(r.SubtractHadron(pr)) r-=pr;
        else if(r.SubtractHadron(dp)) r-=dp;
        else if(r.SubtractHadron(la)) r-=la;
        else
        {
          //#ifdef erdebug
          G4cerr<<"***G4QContent::SplitChipo:ADib(2) tot=6,b=2, qCont="<<GetThis()<<G4endl;
          //#endif
        }
      }
      else
      {
        if     (r.SubtractHadron(la)) r-=la; // These functions only check
        else if(r.SubtractHadron(ks)) r-=ks;
        else if(r.SubtractHadron(om)) r-=om;
        else if(r.SubtractHadron(pr)) r-=pr;
        else if(r.SubtractHadron(nt)) r-=nt;
        else
        {
          //#ifdef erdebug
          G4cerr<<"***G4QContent::SplitChipo:ADib(3) tot=6,b=2, qCont="<<GetThis()<<G4endl;
          //#endif
        }
      }
    }
  }
  else // More than Dibaryon (@@ can use the same algorithm as for dibaryon)
  {
    //#ifdef erdebug
    G4cerr<<"*G4QContent::SplitChipolino:UnknownHadron with QuarkCont="<<GetThis()<<G4endl;
    //#endif
  }
  return r;
}// End of G4QContent::SplitChipolino

// Return one-quark QC using index (a kind of iterator)
G4QContent G4QContent::IndQ (G4int index)
{//        ==============================
#ifdef debug
  G4cout << "G4QC::IndQ is called"<<G4endl;
#endif
  if(index<nD) return G4QContent(1,0,0,0,0,0);
  else if(index<nD+nU) return G4QContent(0,1,0,0,0,0);
  else if(index<nD+nU+nS) return G4QContent(0,0,1,0,0,0);
  //#ifdef erdebug
  else G4cerr<<"***G4QC::IndQ:index="<<index<<" for the QuarkContent="<<GetThis()<<G4endl;
  //throw G4QException("***G4QC::IndQ: Index exceeds the total number of quarks");
  //#endif
  return G4QContent(0,0,0,0,0,0);
}

// Return one-antiquark QC using index (a kind of iterator)
G4QContent G4QContent::IndAQ (G4int index)
{//        ==============================
#ifdef debug
  G4cout << "G4QC::IndAQ is called"<<G4endl;
#endif
  if(index<nAD) return G4QContent(0,0,0,1,0,0);
  else if(index<nAD+nAU) return G4QContent(0,0,0,0,1,0);
  else if(index<nAD+nAU+nAS) return G4QContent(0,0,0,0,0,1);
  //#ifdef erdebug
  else G4cerr<<"***G4QC::IndAQ:index="<<index<<" for the QuarkContent="<<GetThis()<<G4endl;
  //throw G4QException("***G4QC::IndAQ: Index exceeds the total number of antiquarks");
  //#endif
  return G4QContent(0,0,0,0,0,0);
}

// Reduces number (if nQAQ<0:all) of valence Q-Qbar pairs, returns #of pairs to reduce more
G4int G4QContent::DecQAQ(const G4int& nQAQ)
{//   =====================================
#ifdef debug
  G4cout<<"G4QCont::DecQC: n="<<nQAQ<<","<<GetThis()<<G4endl;
#endif
  G4int ban = GetBaryonNumber();
  G4int tot = GetTot();    // Total number of quarks in QC
  if (tot==ban*3) return 0;// Nothing to reduce & nothing to reduce more
  G4int nUP=0;             // U/AU min factor (a#of u which can be subtracted)
  if (nU>=nAU) nUP=nAU;
  else         nUP= nU;

  G4int nDP=0;             // D/AD min factor (a#of d which can be subtracted)
  if (nD>=nAD) nDP=nAD;
  else         nDP= nD;

  G4int nSP=0;             // S/AS min factor (a#of s which can be subtracted)
  if (nS>=nAS) nSP=nAS;
  else         nSP= nS;

  G4int nLP  =nUP+nDP;     // a#of light quark pairs which can be subtracted
  G4int nTotP=nLP+nSP;     // total#of existing pairs which can be subtracted
  G4int nReal=nQAQ;        // demanded #of pairs for reduction (by demand)
  G4int nRet =nTotP-nQAQ;  // a#of additional pairs for further reduction
  if (nQAQ<0)              // === Limited reduction case @@ not tuned for baryons !!
  {
    G4int res=tot+nQAQ;
#ifdef debug
 G4cout<<"G4QC::DecQC: tot="<<tot<<", nTP="<<nTotP<<", res="<<res<<G4endl;
#endif
    if(res<0)
    {
      IncQAQ(1,0.);        // Increment by one not strange pair to get the minimum
      return 1;
    }
    res -=nTotP+nTotP;
    if(res<0) nReal=nTotP+res/2;
    else nReal=nTotP;
    nRet = tot/2-nReal;
  }
  else if(!nQAQ)
  {
    nReal=nTotP;
    nRet =0;
  }
  else if(nRet<0) nReal=nTotP;

  if (!nReal) return nRet; // Now nothing to be done
  // ---------- Decrimenting by nReal pairs
#ifdef debug
  G4cout<<"G4QC::DecQC: demanded "<<nQAQ<<" pairs, executed "<<nReal<<" pairs"<<G4endl;
#endif
  ///////////G4int nt = tot - nTotP - nTotP;
  for (G4int i=0; i<nReal; i++)
  {
    G4double base = nTotP;
    //if (nRet && nSP==1 && !nQAQ) base = nLP;              // Keep S-Sbar pair if possible
    G4int j = static_cast<int>(base*G4UniformRand());       // Random integer "SortOfQuark"
    if (nUP && j<nUP && (nRet>2 || nUP>1 || (nD<2 && nS<2)))// --- U-Ubar pair
    {
#ifdef debug
      G4cout<<"G4QC::DecQC: decrementing UAU pair UP="<<nUP<<", QC="<<GetThis()<<G4endl;
#endif
      nU--;
      nAU--;
      nUP--;
      nLP--;
      nTotP--;
    } 
    else if (nDP && j<nLP && (nRet>2 || nDP>1 || (nU<2 && nS<2)))// --- D-Ubar pair
    {
#ifdef debug
      G4cout<<"G4QC::DecQC: decrementing DAD pair DP="<<nDP<<", QC="<<GetThis()<<G4endl;
#endif
      nD--;
      nAD--;
      nDP--;
      nLP--;
      nTotP--;
    } 
    else if (nSP&& (nRet>2 || nSP>1 || (nU<2 && nD<2)))          // --- S-Sbar pair
    {
#ifdef debug
      G4cout<<"G4QC::DecQC: decrementing SAS pair SP="<<nSP<<", QC="<<GetThis()<<G4endl;
#endif
      nS--;
      nAS--;
      nSP--;
      nTotP--;
    }
    else if (nUP)                                  // --- U-Ubar pair cancelation (final)
    {
#ifdef debug
      G4cout<<"G4QC::DecQC:Decrement UAU pair (final) UP="<<nUP<<",QC="<<GetThis()<<G4endl;
#endif
      nU--;
      nAU--;
      nUP--;
      nLP--;
      nTotP--;
    } 
    else if (nDP)                                 // --- D-Ubar pair cancelation (final)
    {
#ifdef debug
      G4cout<<"G4QC::DecQC:Decrement DAD pair (final) DP="<<nDP<<",QC="<<GetThis()<<G4endl;
#endif
      nD--;
      nAD--;
      nDP--;
      nLP--;
      nTotP--;
    } 
    else if (nSP)                                 // --- S-Sbar pair cancelation (final)
    {
#ifdef debug
      G4cout<<"G4QC::DecQC: decrementing SAS pair SP="<<nSP<<", QC="<<GetThis()<<G4endl;
#endif
      nS--;
      nAS--;
      nSP--;
      nTotP--;
    }
    else G4cout<<"***G4QC::DecQC:i="<<i<<",j="<<j<<",D="<<nDP<<",U="<<nUP<<",S="<<nSP
               <<",T="<<nTotP<<",nRet="<<nRet<<", QC="<<GetThis()<<G4endl;
  }
#ifdef debug
  G4cout<<"G4QC::DecQC: >>>OUT<<< nRet="<<nRet<<", QC="<<GetThis()<<G4endl;
#endif
  return nRet;
}

// Increment quark pairs
void G4QContent::IncQAQ(const G4int& nQAQ, const G4double& sProb)
{//  ============================================================
  G4int tot = GetTot();
  G4QContent mQC = GetThis();
  for (int i=0; i<nQAQ; i++)
  {
    G4int j = static_cast<int>((2.+sProb)*G4UniformRand()); // 0-U, 1-D, 2-S
#ifdef debug
    G4cout<<"IncQC:out QC="<<GetThis()<<",j="<<j<<" for i="<<i<<G4endl;
#endif
    //if      (!j)
    if      ( !j && (nU<=nD || nU<=nS))
    {
      nU++;
      nAU++;
      tot+=2;
    }
    //else if (j==1)
    else if (j==1 && (nD<=nU || nD<=nS))
    {
      nD++;
      nAD++;
      tot+=2;
    }
    //else
    else if (j>1&& (nS<=nU || nS<=nD))
    {
      nS++;
      nAS++;
      tot+=2;
    }
    else if (!j)
    {
      nD++;
      nAD++;
      tot+=2;
    }
    else if (j==1)
    {
      nU++;
      nAU++;
      tot+=2;
    }      
    else
    {
      nS++;
      nAS++;
      tot+=2;
    }
    //else if (nD<=nU)
    //{
    //  nD++;
    //  nAD++;
    //  tot+=2;
    //}
    //else
    //{
    //  nU++;
    //  nAU++;
    //  tot+=2;
    //}      
  }
}

// Calculate a#of protons in the QC (for nuclei)
G4int G4QContent::GetP() const
{//   ========================
  G4int rD=nD-nAD;                                   // Constituent d-quarks
  G4int rU=nU-nAU;                                   // Constituent u-quarks
  G4int rS=nS-nAS;                                   // Constituent s-quarks
  G4int dQ=rD-rU;                                    // Isotopic assimetry
  G4int b3=rD+rU+rS;                                 // (Baryon number) * 3
  return (b3-3*(rS+dQ))/6;
}

// Calculate a#of neutrons in the QC (for nuclei)
G4int G4QContent::GetN() const
{//   ========================
  G4int rD=nD-nAD;                                   // Constituent d-quarks
  G4int rU=nU-nAU;                                   // Constituent u-quarks
  G4int rS=nS-nAS;                                   // Constituent s-quarks
  G4int dQ=rD-rU;                                    // Isotopic assimetry
  G4int b3=rD+rU+rS;                                 // (Baryon number) * 3
  return (b3-3*(rS-dQ))/6;
}

// Calculate a#of lambdas in the QC (for nuclei)
G4int G4QContent::GetL() const
{//   ========================
  G4int rS=nS-nAS;                                   // Constituent s-quarks
  return rS;
}

// Calculate a#of anti-protons in the QC (for anti-nuclei)
G4int G4QContent::GetAP() const
{//   =========================
  G4int rD=nAD-nD;                                   // Constituent anti-d-quarks
  G4int rU=nAU-nU;                                   // Constituent anti-u-quarks
  G4int rS=nAS-nS;                                   // Constituent anti-s-quarks
  G4int dQ=rD-rU;                                    // Isotopic assimetry
  G4int b3=rD+rU+rS;                                 // - (Baryon number) * 3
  return (b3-3*(rS+dQ))/6;
}

// Calculate a#of anti-neutrons in the QC (for anti-nuclei)
G4int G4QContent::GetAN() const
{//   =========================
  G4int rD=nAD-nD;                                   // Constituent anti-d-quarks
  G4int rU=nAU-nU;                                   // Constituent anti-u-quarks
  G4int rS=nAS-nS;                                   // Constituent anti-s-quarks
  G4int dQ=rD-rU;                                    // Isotopic assimetry
  G4int b3=rD+rU+rS;                                 // - (Baryon number) * 3
  return (b3-3*(rS-dQ))/6;
}

// Calculate a#of anti-lambdas in the QC (for anti-nuclei)
G4int G4QContent::GetAL() const
{//   =========================
  G4int rS=nAS-nS;                                   // Constituent anti-s-quarks
  return rS;
}

// Calculate charge for the QC
G4int G4QContent::GetCharge() const
{//   =============================
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
  //#ifdef erdebug
  if(c%3) G4cerr<<"***G4QCont:GetCharge:c="<<c<<"/3 isn't integer, QC="<<GetThis()<<G4endl;
  //#endif
  return c/3;
}

// Calculate a Baryon Number for the QC
G4int G4QContent::GetBaryonNumber() const
{//   ===================================
  G4int b=nU+nD+nS-nAU-nAD-nAS;
  //#ifdef erdebug
  if(b%3)
  {
     G4cerr<<"-Warning-G4QContent::GetBaryonNumber="<<b<<"/3 isn't anIntegerValue"<<G4endl;
     G4Exception("G4QContent::GetBaryonNumber:","72",FatalException,"Wrong Baryon Number");
  }
  //#endif
  return b/3;
}

// Gives the PDG of the lowest (in mass) S-wave hadron or Chipolino (=10) for double hadron
G4int G4QContent::GetSPDGCode() const
{//   ===============================
  G4int p = 0;           // Prototype of output SPDG
  G4int n = GetTot();    // Total number of quarks
  if(!n) return 22;      // Photon does not have any Quark Content
  G4int mD=nD;           // A # of D quarks or anti U quarks
  if (nD<=0) mD=nAD;
  G4int mU=nU;           // A # of U quarks or anti U quarks
  if (nU<=0) mU=nAU;
  G4int mS=nS;           // A # of S quarks or anti U quarks
  if (nS<=0) mS= nAS;
  // ---------------------- Cancelation of q-qbar pairs in case of an excess
  if ( nU>nAU && nAU>0)
  {
    mU=nU-nAU;
    n-=nAU+nAU;
  }
  if (nAU>nU  &&  nU>0)
  {
    mU=nAU-nU;
    n-=nU+nU;
  }
  if ( nD>nAD && nAD>0)
  {
    mD=nD-nAD;
    n-=nAD+nAD;
  }
  if (nAD>nD  &&  nD>0)
  {
    mD=nAD-nD;
    n-=nD+nD;
  }
  if ( nS>nAS && nAS>0)
  {
    mS=nS-nAS;
    n-=nAS+nAS;
  }
  if (nAS>nS  &&  nS>0)
  {
    mS= nAS-nS;
    n-=nS+nS;
  }
  // ---------------------- Cancelation of q-qbar pairs in case of an equality
  if (nAD==nD  &&  nD>0)
  {
    G4int dD=nD+nD;
    if(n>dD)
    {
      mD=0;
      n-=dD;
    }
    else if (n==dD)
    {
      mD=2;
      n=2;
    }
    else
    {
#ifdef debug
      G4cout<<"***G4QC::SPDG:CanD U="<<mU<<",D="<<mD<<",S="<<mS<<",QC="<<GetThis()<<G4endl;
#endif
      return 0;
    }
  }
  if (nAU==nU  &&  nU>0)
  {
    G4int dU=nU+nU;
    if(n>dU)
    {
      mU=0;
      n-=dU;
    }
    else if (n==dU)
    {
      mU=2;
      n=2;
    }
    else
    {
#ifdef debug
      G4cout<<"***G4QC::SPDG:CanU U="<<mU<<",D="<<mD<<",S="<<mS<<",QC="<<GetThis()<<G4endl;
#endif
      return 0;
    }
  }
  if (nAS==nS  &&  nS>0) //@@ Starts with S-quarks - should be randomized and mass limited
  {
    G4int dS=nS+nS;
    if(n>dS)
    {
      mS=0;
      n-=dS;
    }
    else if (n==dS)
    {
      mS=2;
      n=2;
    }
    else
    {
#ifdef debug
      G4cout<<"***G4QC::SPDG:CanS U="<<mU<<",D="<<mD<<",S="<<mS<<",QC="<<GetThis()<<G4endl;
#endif
      return 0;
    }
  }

  G4int b=GetBaryonNumber();
  G4int c=GetCharge();
  G4int s=GetStrangeness();
#ifdef debug
  G4cout<<"G4QC::SPDGC:bef. b="<<b<<",n="<<n<<",c="<<c<<",s="<<s<<",Q="<<GetThis()<<G4endl;
#endif
  if (b)                                         // ==================== Baryon case
  {
    
    G4int ab=abs(b);
    if(ab>=2 && n>=6)                            // Multi-Baryonium (NuclearFragment)
    {
      G4int mI=nU-nAU-nD+nAD;
      //if     (abs(mI)>3||mS>3||(b>0&&s<-1)||(b<0&&s>1)) return  0;
      //else if(abs(mI)>2||mS>2||(b>0&&s< 0)||(b<0&&s>0)) return 10;
      if ( (b > 0 && s < -1) || (b < 0 && s > 1) ) return 10;
      else if (abs(mI) > 2 || mS > 2 
                           || (b > 0 && s < 0) 
                           || (b < 0 && s > 0)) return GetZNSPDGCode();
      else if(mU>=mS&&mD>=mS&&mU+mD+mS==3*b)     // Possible Unary Nuclear Cluster
      {
        G4int mZ=(mU+mD-mS-mS+3*mI)/6;
        p = 90000000+1000*(1000*mS+mZ)+mZ-mI;
        if(b>0) return  p;
        else    return -p;
      }
      else return 10;
    }
    // Normal One Baryon States: Heavy quark should come first
    if(n>5) return GetZNSPDGCode();            //B+M+M Tripolino etc
    if(n==5) return 10;                        //B+M Chipolino
    if(mS>0)                                   // Strange Baryons
    {
      p=3002;
      if      (mS==3)            p+=332;       // Decuplet
      else if (mS==2)
      {
        if      (mU==1 && mD==0) p+=320;
        else if (mU==0 && mD==1) p+=310;
        else
        {
#ifdef debug
          G4cout<<"**G4QC::SPDG:ExoticBSS,U="<<mU<<",D="<<mD<<",S="<<mS<<GetThis()<<G4endl;
#endif
          return GetZNSPDGCode();
        }
      }
      else if (mS==1)
      {
        if      (mU==2 && mD==0) p+=220;
        else if (mU==1 && mD==1) p+=120;       // Lambda (M_Lambda<M_Sigma0) PDG_Sigma=3212
        else if (mU==0 && mD==2) p+=110;
        else
        {
#ifdef debug
          G4cout<<"***G4QC::SPDG:ExoticBS,U="<<mU<<",D="<<mD<<",S="<<mS<<GetThis()<<G4endl;
#endif
          return GetZNSPDGCode();
        }
      }
      else                                     // Superstrange case
      {
#ifdef debug
        G4cout<<"***G4QC::GetSPDG:ExoBarS,U="<<mU<<",D="<<mD<<",S="<<mS<<GetThis()<<G4endl;
#endif
        return GetZNSPDGCode();
      }
    }
    else if (mU>0)                               // Not Strange Baryons
    {
      p=2002;
      if      (mU==3 && mD==0) p+=222;           // Decuplet
      else if (mU==2 && mD==1) p+=210;
      else if (mU==1 && mD==2) p+=110;           // There is a higher Delta S31 (PDG=1212)
      else
      {
#ifdef debug
        G4cout<<"***G4QC::SPDG:ExoBaryonU,U="<<mU<<",D="<<mD<<",S="<<mS<<GetThis()<<G4endl;
#endif
        return GetZNSPDGCode();
      }
    }
    else if (mD==3) p=1114;                      // Decuplet
    else
    {
#ifdef debug
      G4cout<<"**G4QC::SPDG:ExoticBaryonD,U="<<mU<<",D="<<mD<<",S="<<mS<<GetThis()<<G4endl;
#endif
      return GetZNSPDGCode();
    }
    if (b<0) p=-p;
  }
  else             // ====================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Meson case
  {
#ifdef debug
    G4cout<<"G4QC::SPDG:mDUS="<<mD<<","<<mU<<","<<mS<<",b,c,s="<<b<<","<<c<<","<<s<<G4endl;
#endif
    if(n>4)                                      // Super Exotics
    {
#ifdef debug
      G4cout<<"G4QC::SPDG:n>4 SEx:U="<<mU<<",D="<<mD<<",S="<<mS<<",QC="<<GetThis()<<G4endl;
#endif
      return 0;
    }
    if(n==4) return 10;                          // M+M Chipolino
    if(abs(s)>1)
    {
#ifdef debug
      G4cout<<"**G4QC::SPDG:Stran="<<s<<",QC="<<GetThis()<<" - Superstrange Meson"<<G4endl;
#endif
      return 0;
    }
    // Heavy quark should come first
    if(mS>0)                                     // Strange Mesons
    {
      p=301;
      if      (mS==2)
      {
        //if (G4UniformRand()<0.333) p=221;        // eta
        //else                       p+=30;        // eta'
        p=221;
      }
      else if (mU==1 && mD==0) p+=20;
      else if (mU==0 && mD==1) p+=10;
      else
      {
#ifdef debug
        G4cout<<"*G4QC::SPDG:ExMS U="<<mU<<",D="<<mD<<",S="<<mS<<",QC="<<GetThis()<<G4endl;
#endif
        return 0;
      }
    }
    else if (mU>0)                               // Isotopic Mesons
    {
      p=201;
      //if      (mU==2 && mD==0) p=221; // Performance problems
      if      (mU==2 && mD==0) p=111;
      else if (mU==1 && mD==1) p+=10;
      else
      {
#ifdef debug
        G4cout<<"*G4QC::SPDG:ExMU U="<<mU<<",D="<<mD<<",S="<<mS<<",QC="<<GetThis()<<G4endl;
#endif
        return 0;
      }
    }
    else if (mD==2) p=111;
    else
    {
#ifdef debug
      G4cout<<"***G4QC::SPDG:ExMD U="<<mU<<",D="<<mD<<",S="<<mS<<",QC="<<GetThis()<<G4endl;
#endif
      return 0;
    }
    if (c<0 || (c==0 && mS>0 && s>0)) p=-p;
  }
#ifdef debug
  G4cout<<"G4QC::GetSPDG:output SPDGcode="<<p<<" for the QuarkContent="<<GetThis()<<G4endl;
#endif
  return p;
}

// === Calculate a number of combinations of rhc out of lhc ==
G4int G4QContent::NOfCombinations(const G4QContent& rhs) const
{//   ========================================================
  G4int c=1; // Default number of combinations?
#ifdef ppdebug
  G4cout<<"G4QContent::NOfComb: This="<<GetThis()<<", selectQC="<<rhs<<G4endl;
#endif
  G4int mD=rhs.GetD();
  G4int mU=rhs.GetU();
  G4int mS=rhs.GetS();
  G4int mAD=rhs.GetAD();
  G4int mAU=rhs.GetAU();
  G4int mAS=rhs.GetAS();
  G4int mN=mD+mU+mS-mAD-mAU-mAS;
  ////////////G4int PDG=abs(GetSPDGCode());
  if (( ((nD < mD || nAD < mAD) && !(mD-mAD)) || 
        ((nU < mU || nAU < mAU) && !(mU-mAU)) ||
        ((nS < mS || nAS < mAS) && !(mS-mAS)) ) && !mN) return 1;
  if(mD>0)
  {
    int j=nD;
    if (j<=0) return 0;
    if(mD>1||j>1) for (int i=1; i<=mD; i++)
    {
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  }
  if(mU>0)
  {
    int j=nU;
    if (j<=0) return 0;
    if(mU>1||j>1) for (int i=1; i<=mU; i++)
    {
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  }
  if(mS>0)
  {
    int j=nS;
    if (j<=0) return 0;
    if(mS>1||j>1) for (int i=1; i<=mS; i++)
    {
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  }
  if(mAD>0)
  {
    int j=nAD;
    if (j<=0) return 0;
    if(mAD>1||j>1) for (int i=1; i<=mAD; i++)
    {
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  }
  if(mAU>0)
  {
    int j=nAU;
    if (j<=0) return 0;
    if(mAU>1||j>1) for (int i=1; i<=mAU; i++)
    {
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  }
  if(mAS>0)
  {
    int j=nAS;
    if (j<=0) return 0;
    if(mAS>1||j>1) for (int i=1; i<=mAS; i++)
    {
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  } 
  return c;
} 

// Make PDG's of PartonPairs for Mesons & Baryons (only)
std::pair<G4int,G4int> G4QContent::MakePartonPair() const
{
  G4double S=0.;
  S+=nD;
  G4double dP=S;
  S+=nU;
  G4double uP=S;
  S+=nS;
  G4double sP=S;
  S+=nAD;
  G4double dA=S;
  S+=nAU;
  G4double uA=S;
  S+=nAS;
  if(!S)
  {
    G4int f= static_cast<int>(1.+2.3*G4UniformRand()); // Random flavor @@ a Parameter
    return std::make_pair(f,-f);
  }
  G4int f=0;
  G4double R=S*G4UniformRand();
  if     (R<dP) f=1;
  else if(R<uP) f=2;
  else if(R<sP) f=3;
  else if(R<dA) f=-1;
  else if(R<uA) f=-2;
  else          f=-3;

  if (f < 0) {      // anti-quark
    if(nD || nU || nS) // a Meson
    {
      if     (nD) return std::make_pair(1,f);
      else if(nU) return std::make_pair(2,f);
      else        return std::make_pair(3,f);
    }
    else               // Anti-Baryon
    {
      // @@ Can be improved, taking into acount weights (i,i): w=3, (i,j#i): w=4(s=0 + s=1)
      G4int AD=nAD;
      if(f==-1) AD--;
      G4int AU=nAU;
      if(f==-1) AU--;
      G4int AS=nAS;
      if(f==-1) AS--;
      if       (AS)
      {
        if     (AS==2) return std::make_pair(-3303,f); // 3301 does not exist
        else if(AU)    return std::make_pair(-3201,f); // @@ only lightest
        else           return std::make_pair(-3101,f); // @@ only lightest
      }
      else if(AU)
      {
        if     (AU==2) return std::make_pair(-2203,f); // 2201 does not exist
        else           return std::make_pair(-2101,f); // @@ only lightest
      }
      else return std::make_pair(-1103,f); // 1101 does not exist
    }

  } else {        // quark (f is a PDG code of the quark)

    if(nAD || nAU || nAS) // a Meson
    {
      if     (nAD) return std::make_pair(f,-1);
      else if(nAU) return std::make_pair(f,-2);
      else         return std::make_pair(f,-3);
    }
    else               // Anti-Baryon
    {
      // @@ Can be improved, taking into acount weights (i,i): w=3, (i,j#i): w=4(s=0 + s=1)
      //  G4int AD=nD;      AD is unused
      //  if(f==-1) AD--;   dead code: f >= 0 in this branch  (DHW 11 Nov 2010)
      G4int AU=nU;
      //  if(f==-1) AU--;   dead code: f >= 0 in this branch  (DHW 11 Nov 2010)
      G4int AS=nS;
      //  if(f==-1) AS--;   dead code: f >= 0 in this branch  (DHW 11 Nov 2010)

      if (AS) {
        if     (AS==2) return std::make_pair(f,3303); // 3301 does not exist
        else if(AU)    return std::make_pair(f,3201); // @@ only lightest
        else           return std::make_pair(f,3101); // @@ only lightest

      } else if(AU) {
        if     (AU==2) return std::make_pair(f,2203); // 2201 does not exist
        else           return std::make_pair(f,2101); // @@ only lightest
      }
      else return std::make_pair(f,1103); // 1101 does not exist
    }
  }  // test on f
}


// Add parton (pPDG) to the hadron (this QC) & get parton PDG
// (Baryons,Mesons,Anti-Baryons)

G4int G4QContent::AddParton(G4int pPDG) const
{
#ifdef debug
  G4cout<<"G4QContent::AddParton: This="<<GetThis()<<", pPDG="<<pPDG<<G4endl;
#endif
  if(!pPDG || pPDG==9 || pPDG==21)
  {
#ifdef debug
    G4cout<<"-Warning-G4QContent::AddParton: ImpossibleToAdd PartonWithPDG="<<pPDG<<G4endl;
#endif
    return 0;
  }
  G4int aPDG = std::abs(pPDG);
  if( (aPDG>3 && aPDG<1101) || pPDG>3303) // @@ 1101 does not exist
  {
#ifdef debug
    G4cout<<"-Warning-G4QContent::AddParton: Impossible Parton with PDG="<<pPDG<<G4endl;
#endif
    return 0;
  }
  G4int HBN = GetBaryonNumber();
  if( HBN > 1 || HBN <-1)
  {
#ifdef debug
    G4cout<<"-Warning-G4QContent::AddParton: Impossible Hadron with BaryonN="<<HBN<<G4endl;
#endif
    return 0;
  }
  G4int AD=nAD;
  G4int AU=nAU;
  G4int AS=nAS;
  G4int QD=nD;
  G4int QU=nU;
  G4int QS=nS;
  if(aPDG>99)                             // Parton is DiQuark/antiDiQuark
  {
    G4int rPDG=aPDG/100;
    G4int P1=rPDG/10;                     // First quark
    G4int P2=rPDG%10;                     // Second quark
#ifdef debug
    G4cout<<"G4QContent::AddParton: DiQuark/AntiDiQuark, P1="<<P1<<", P2="<<P2<<G4endl;
#endif
    if(pPDG>0)                            // -- DiQuark
    {
#ifdef debug
      G4cout<<"G4QContent::AddParton: DiQuark, P1="<<P1<<", P2="<<P2<<",HBN="<<HBN<<G4endl;
#endif
      if     (P1==3 && P2==3)             // ----> ss DiQuark
      {
        if(HBN<0 && AS>1) AS-=2;          // >> Annihilation of ss-DiQuark with anti-Baryon
        else if(!HBN && AS==1)
        {
          AS=0;
          ++QS;
        }
        else if(HBN || (!HBN && !AS)) return 0;
      }
      else if(P1==3 && P2==2)             // ----> su DiQuark
      {
        if(HBN<0 && AS && AU)             // >> Annihilation of su-DiQuark with anti-Baryon
        {
          --AS;
          --AU;
        }
        else if(!HBN && (AS || AU))
        {
          if(AS)
          {
            --AS;
            ++QU;
          }
          else
          {
            --AU;
            ++QS;
          }
        }
        else if(HBN || (!HBN && !AS && !AU)) return 0;
      }
      else if(P1==3 && P2==1)             // ----> sd DiQuark
      {
        if(HBN<0 && AS && AD)             // >> Annihilation of sd-DiQuark with anti-Baryon
        {
          --AS;
          --AD;
        }
        else if(!HBN && (AS || AD))
        {
          if(AS)
          {
            --AS;
            ++QD;
          }
          else
          {
            --AD;
            ++QS;
          }
        }
        else if(HBN || (!HBN && !AS && !AD)) return 0;
      }
      else if(P1==2 && P2==2)             // ----> uu DiQuark
      {
        if(HBN<0 && AU>1) AU-=2;          // >> Annihilation of uu-DiQuark with anti-Baryon
        else if(!HBN && AU==1)
        {
          AU=0;
          ++QU;
        }
        else if(HBN || (!HBN && !AU)) return 0;
      }
      else if(P1==2 && P2==1)             // ----> ud DiQuark
      {
        if(HBN<0 && AD && AU)             // >> Annihilation of ud-DiQuark with anti-Baryon
        {
          --AD;
          --AU;
        }
        else if(!HBN && (AD || AU))
        {
          if(AD)
          {
            --AD;
            ++QU;
          }
          else
          {
            --AU;
            ++QD;
          }
        }
        else if(HBN || (!HBN && !AU && !AD)) return 0;
      }
      else                                // ----> dd DiQuark
      {
        if(HBN<0 && AD>1) AD-=2;          // >> Annihilation of dd-DiQuark with anti-Baryon
        else if(!HBN && AD==1)
        {
          AD=0;
          ++QD;
        }
        else if(HBN || (!HBN && !AD)) return 0;
      }
#ifdef debug
      G4cout<<"G4QContent::AddParton: DQ, QC="<<QD<<","<<QU<<","<<QS<<","<<AD<<","<<AU<<","
            <<AS<<G4endl;
#endif
      if     (HBN<0)                      // ....... Hadron is an Anti-Baryon
      {
        if     (AD) return -1;            // >>>>>>> Answer is anti-d
        else if(AU) return -2;            // >>>>>>> Answer is anti-u
        else        return -3;            // >>>>>>> Answer is anti-s
      }
      else                                // ... Hadron is aMeson with annihilatedAntiQuark
      {					      
        if    (QS)                       // --------- There is an s-quark
        {					      
          if  (QS==2) return 3303;       // >>>>>>> Answer is ss (3301 does not exist)
          else if(QU) return 3201;       // >>>>>>> Answer is su (@@ only lightest)
          else        return 3101;       // >>>>>>> Answer is sd (@@ only lightest)
        }
        else if(QU)                      // --------- There is an u quark
        {					      
          if  (QU==2) return 2203;       // >>>>>>> Answer is uu (2201 does not exist)
          else        return 2101;       // >>>>>>> Answer is ud (@@ only lightest)
        }
        else          return 1103;       // >>>>>>> Answer is dd (1101 does not exist)
      }
    }
    else                                  // -- antiDiQuark
    {
#ifdef debug
      G4cout<<"G4QContent::AddParton: AntiDiQuark,P1="<<P1<<",P2="<<P2<<",B="<<HBN<<G4endl;
#endif
      if     (P1==3 && P2==3)             // ----> anti-s-anti-s DiQuark
      {
        if(HBN>0 && QS>1) QS-=2;          // >> Annihilation of anti-ss-DiQuark with Baryon
        else if(!HBN && QS==1)
        {
          QS=0;
          ++AS;
        }
        else if(HBN || (!HBN && !QS)) return 0;
      }
      else if(P1==3 && P2==2)             // ----> anti-s-anti-u DiQuark
      {
        if(HBN>0 && QS && QU)             // >> Annihilation of anti-su-DiQuark with Baryon
        {
          --QS;
          --QU;
        }
        else if(!HBN && (QS || QU))
        {
          if(QS)
          {
            --QS;
            ++AU;
          }
          else
          {
            --QU;
            ++AS;
          }
        }
        else if(HBN || (!HBN && !QS && !QU)) return 0;
      }
      else if(P1==3 && P2==1)             // ----> anti-s-anti-d DiQuark
      {
        if(HBN>0 && QS && QD)             // >> Annihilation of anti-sd-DiQuark with Baryon
        {
          --QS;
          --QD;
        }
        else if(!HBN && (QS || QD))
        {
          if(QS)
          {
            --QS;
            ++AD;
          }
          else
          {
            --QD;
            ++AS;
          }
        }
        else if(HBN || (!HBN && !QS && !QD)) return 0;
      }
      else if(P1==2 && P2==2)             // ----> anti-u-anti-u DiQuark
      {
        if(HBN>0 && QU>1) QU-=2;          // >> Annihilation of anti-uu-DiQuark with Baryon
        else if(!HBN && QU==1)
        {
          QU=0;
          ++AU;
        }
        else if(HBN || (!HBN && !QU)) return 0;
      }
      else if(P1==2 && P2==1)             // ----> anti-u-anti-d DiQuark
      {
        if(HBN>0 && QU && QD)             // >> Annihilation of anti-ud-DiQuark with Baryon
        {
          --QU;
          --QD;
        }
        else if(!HBN && (QU || QD))
        {
          if(QU)
          {
            --QU;
            ++AD;
          }
          else
          {
            --QD;
            ++AU;
          }
        }
        else if(HBN || (!HBN && !QU && !QD)) return 0;
      }
      else                                // ----> anti-d=anti-d DiQuark
      {
        if(HBN>0 && QD>1) QD-=2;          // >> Annihilation of anti-dd-DiQuark with Baryon
        else if(!HBN && QD==1)
        {
          QD=0;
          ++AD;
        }
        else if(HBN || (!HBN && !QD)) return 0;
      }
#ifdef debug
      G4cout<<"G4QContent::AddParton:ADQ, QC="<<QD<<","<<QU<<","<<QS<<","<<AD<<","<<AU<<","
            <<AS<<G4endl;
#endif
      if     (HBN>0)                      // ....... Hadron is an Baryon
      {
        if     (QD) return 1;             // >>>>>>> Answer is d
        else if(QU) return 2;             // >>>>>>> Answer is u
        else        return 3;             // >>>>>>> Answer is s
      }
      else                                // ....... Meson with annihilated Anti-Quark
      {					      
        if    (AS)                       // --------- There is an anti-s quark
        {					      
          if  (AS==2) return -3303;      // >>>>>>> Answer is anti-ss (3301 does not exist)
          else if(AU) return -3201;      // >>>>>>> Answer is anti-su (@@ only lightest)
          else        return -3101;      // >>>>>>> Answer is anti-sd (@@ only lightest)
        }
        else if(AU)                      // --------- There is an anti-u quark
        {					      
          if  (AU==2) return -2203;      // >>>>>>> Answer is anti-uu (2201 does not exist)
          else        return -2101;      // >>>>>>> Answer is anti-ud (@@ only lightest)
        }
        else          return -1103;      // >>>>>>> Answer is anti-dd (1101 does not exist)
      }
    }
  }
  else                                    // Parton is Quark/antiQuark
  {
    if(pPDG>0)                            // -- Quark
    {
#ifdef debug
      G4cout<<"G4QContent::AddParton: Quark, A="<<AD<<","<<AU<<","<<AS<<",B="<<HBN<<G4endl;
#endif
      if     (aPDG==1)                    // ----> d quark
      {
        if(HBN<0 && AD) AD--;             // ====> Annihilation of d-quark with anti-Baryon
        else if(HBN || (!HBN && !AD)) return 0;
      }
      else if(aPDG==2)                    // ----> u quark
      {
        if(HBN<0 && AU) AU--;             // ====> Annihilation of u-quark with anti-Baryon
        else if(HBN || (!HBN && !AU)) return 0;
      }
      else                                // ----> s quark
      {
        if(HBN<0 && AS) AS--;             // ====> Annihilation of s-quark with anti-Baryon
        else if(HBN || (!HBN && !AS)) return 0;
      }
#ifdef debug
      G4cout<<"G4QContent::AddParton: Q, QC="<<QD<<","<<QU<<","<<QS<<","<<AD<<","<<AU<<","
            <<AS<<G4endl;
#endif
      if     (!HBN)                       // ....... Hadron is a Meson (passingThrougAbove)
      {					      
        if     (QD) return 1;             // >>>>>>> Answer is d
        else if(QU) return 2;             // >>>>>>> Answer is u
        else        return 3;             // >>>>>>> Answer is s
      }					      
      else                                // ....... AntiBaryon with annihilated AntiQuark
      {					      
        if    (AS)                        // --------- There is an anti-s quark
        {					      
          if  (AS==2) return -3303;       // >>>>>>> Answer is ss (3301 does not exist)
          else if(AU) return -3201;       // >>>>>>> Answer is su (@@ only lightest)
          else        return -3101;       // >>>>>>> Answer is sd (@@ only lightest)
        }					      
        else if(AU)			      
        {					      
          if  (AU==2) return -2203;       // >>>>>>> Answer is uu (2201 does not exist)
          else        return -2101;       // >>>>>>> Answer is ud (@@ only lightest)
        }					      
        else          return -1103;       // >>>>>>> Answer is dd (1101 does not exist)
      }					      
    }
    else                                  // -- antiQuark
    {
#ifdef debug
      G4cout<<"G4QContent::AddParton: antiQ, Q="<<QD<<","<<QU<<","<<QS<<",B="<<HBN<<G4endl;
#endif
      if     (aPDG==1)                    // ----> anti-d quark
      {
        if(HBN>0 && QD) QD--;             // ====> Annihilation of anti-d-quark with Baryon
        else if(HBN || (!HBN && !QD)) return 0;
      }
      else if(aPDG==2)                    // ----> anti-u quark
      {
        if(HBN>0 && QU) QU--;             // ====> Annihilation of anti-u-quark with Baryon
        else if(HBN || (!HBN && !QU)) return 0;
      }
      else                                // ----> anti-s quark
      {
        if(HBN>0 && QS) QS--;             // ====> Annihilation of anti-s-quark with Baryon
        else if(HBN || (!HBN && !QS)) return 0;
      }
#ifdef debug
      G4cout<<"G4QContent::AddParton: AQ, QC="<<QD<<","<<QU<<","<<QS<<","<<AD<<","<<AU<<","
            <<AS<<G4endl;
#endif
      if     (!HBN)                       // ....... Hadron is a Meson (passingThrougAbove)
      {					      
        if     (AD) return -1;            // >>>>>>> Answer is anti-d
        else if(AU) return -2;            // >>>>>>> Answer is anti-u
        else        return -3;            // >>>>>>> Answer is anti-s
      }					      
      else                                // ....... Baryon with annihilated Quark
      {					      
        if    (QS)                        // --------- There is an anti-s quark
        {					      
          if  (QS==2) return 3303;        // >>>>>>> Answer is ss (3301 does not exist)
          else if(QU) return 3201;        // >>>>>>> Answer is su (@@ only lightest)
          else        return 3101;        // >>>>>>> Answer is sd (@@ only lightest)
        }					      
        else if(QU)			      
        {					      
          if  (QU==2) return 2203;        // >>>>>>> Answer is uu (2201 does not exist)
          else        return 2101;        // >>>>>>> Answer is ud (@@ only lightest)
        }					      
        else          return 1103;        // >>>>>>> Answer is dd (1101 does not exist)
      }					      
    }
  }
}
