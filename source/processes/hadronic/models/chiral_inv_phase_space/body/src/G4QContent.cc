// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QContent.cc,v 1.3 2000-08-17 13:53:19 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QContent ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Quasmon initiated Contents used by CHIPS Model
// --------------------------------------------------------------------
// @@ In future total spin & c,b,t of the Hadron can be added @@ M.K.@@

//#define debug
//#define pdebug

#include "G4QContent.hh"

// Constructors
G4QContent::G4QContent(G4int d, G4int u, G4int s, G4int ad, G4int au, G4int as):
  nD(d),nU(u),nS(s),nAD(ad),nAU(au),nAS(as){}

G4QContent::G4QContent(const G4QContent &right)
{
  nU  = right.nU;
  nD  = right.nD;
  nS  = right.nS;
  nAU = right.nAU;
  nAD = right.nAD;
  nAS = right.nAS;
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

// Assignment operator (copy stile for possible Vector extention)
const G4QContent& G4QContent::operator=(const G4QContent &right)
{//               ==============================================
  nU  = right.nU;
  nD  = right.nD;
  nS  = right.nS;
  nAU = right.nAU;
  nAD = right.nAD;
  nAS = right.nAS;
		
  return *this;
}

// Destructor - - - - - - - 
G4QContent::~G4QContent() {}

// Subtract neutral pion from Quark Content (with possible hidden strangeness)
G4bool G4QContent::SubtractPi0()
{//    =========================
#ifdef debug
  cout << "G4QC::SubtractPi0 is called: U="<<nU<<", AU="<<nAU<<", D="<<nD<<", AD="<<nAD<<endl;
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
  cout << "G4QC::SubtractPion is called: U="<<nU<<", AU="<<nAU<<", D="<<nD<<", AD="<<nAD<<endl;
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
  cout << "G4QC::SubtractHadron "<<h<<" is called for QC="<<GetThis()<<endl;
#endif
  if(h.GetU()<=nU&&h.GetD()<=nD&&h.GetS()<=nS&&h.GetAU()<=nAU&&h.GetAD()<=nAD&&h.GetAS()<=nAS)
    return true;
  return false;
}

// Subtract Kaon from Quark Content
G4bool G4QContent::SubtractKaon(G4double mQ)
{//    =====================================
#ifdef debug
  cout<<"G4QC::SubtractKaon is called:S="<<nS<<",AS="<<nAS<<",UDAUAD="<<nU<<nD<<nAU<<nAD<<endl;
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
  cout<<"QC::Can't SubtractKaon:S="<<nS<<",AS="<<nAS<<",UD,AUAD="<<nU<<nD<<","<<nAU<<nAD<<endl;
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
    cerr<<"***G4QCont::SplitChipo: QC="<<GetThis()<<endl;
	//G4Exception("G4Quasmon::FillHadronVector: Cipolino is too fat");
    return Pi;
  }
  else if(tot==4)      // Mesonic (eight possibilities)
  {
    r=GetThis();
    if     (r.SubtractPi0())  return r;
    else if(r.SubtractPion()) return r;
    else if(r.SubtractKaon(mQ)) return r;
    else
    {
      cerr<<"***G4QCont::SplitChipo: Mesonic tot="<<tot<<", b="<<b<<",q="<<q<<",aq="<<aq<<endl;
      return Pi;
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
      cerr<<"***G4QCont::SplitChipo: Baryonic tot=5,b=1,qCont="<<GetThis()<<endl;
      return Pi;
    }
    return r;
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
        if     (r.SubtractHadron(pr)) return r-pr;
        else if(r.SubtractHadron(dp)) return r-dp;
        else if(r.SubtractHadron(nt)) return r-nt;
        else if(r.SubtractHadron(la)) return r-la;
        else if(r.SubtractHadron(dm)) return r-dm;
        else
        {
          cerr<<"***G4QContent::SplitChipo: Dibaryon (1) tot=6, b=2, qCont="<<GetThis()<<endl;
          return Pi;
        }
      }
      else if(nD>=nU&&nD>=nS)
      {
        if     (r.SubtractHadron(nt)) return r-nt;
        else if(r.SubtractHadron(dm)) return r-dm;
        else if(r.SubtractHadron(pr)) return r-pr;
        else if(r.SubtractHadron(dp)) return r-dp;
        else if(r.SubtractHadron(la)) return r-la;
        else
        {
          cerr<<"***G4QContent::SplitChipo: Dibaryon (2) tot=6, b=2, qCont="<<GetThis()<<endl;
          return Pi;
        }
      }
      else
      {
        if     (r.SubtractHadron(la)) return r-la;
        else if(r.SubtractHadron(ks)) return r-ks;
        else if(r.SubtractHadron(om)) return r-om;
        else if(r.SubtractHadron(pr)) return r-pr;
        else if(r.SubtractHadron(nt)) return r-nt;
        else
        {
          cerr<<"***G4QContent::SplitChipo: Dibaryon (3) tot=6, b=2, qCont="<<GetThis()<<endl;
          return Pi;
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
        if     (r.SubtractHadron(pr)) return r-pr;
        else if(r.SubtractHadron(dp)) return r-dp;
        else if(r.SubtractHadron(nt)) return r-nt;
        else if(r.SubtractHadron(la)) return r-la;
        else if(r.SubtractHadron(dm)) return r-dm;
        else
        {
          cerr<<"***G4QContent::SplitChipo: ADibaryon (1) tot=6, b=2, qCont="<<GetThis()<<endl;
          return Pi;
        }
      }
      else if(nAD>=nAU&&nAD>=nAS)
      {
        if     (r.SubtractHadron(nt)) return r-nt;
        else if(r.SubtractHadron(dm)) return r-dm;
        else if(r.SubtractHadron(pr)) return r-pr;
        else if(r.SubtractHadron(dp)) return r-dp;
        else if(r.SubtractHadron(la)) return r-la;
        else
        {
          cerr<<"***G4QContent::SplitChipo: ADibaryon (2) tot=6, b=2, qCont="<<GetThis()<<endl;
          return Pi;
        }
      }
      else
      {
        if     (r.SubtractHadron(la)) return r-la;
        else if(r.SubtractHadron(ks)) return r-ks;
        else if(r.SubtractHadron(om)) return r-om;
        else if(r.SubtractHadron(pr)) return r-pr;
        else if(r.SubtractHadron(nt)) return r-nt;
        else
        {
          cerr<<"***G4QContent::SplitChipo: ADibaryon (3) tot=6, b=2, qCont="<<GetThis()<<endl;
          return Pi;
        }
      }
	}
  }
  else                 // More than Dibaryon (@@ can use the same algorithm as for dibaryon)
  {
    cerr<<"***G4QContent::SplitChipolino: Strange Hadron with QuarkContent="<<GetThis()<<endl;
    return r;
  }
}// End of G4QContent::SplitChipolino

// Return one-quark QC using index (a kind of iterator)
G4QContent G4QContent::IndQ (G4int index)
{//        ==============================
#ifdef debug
  cout << "G4QC::IndQ is called"<<endl;
#endif
  if(index<nD) return G4QContent(1,0,0,0,0,0);
  else if(index<nD+nU) return G4QContent(0,1,0,0,0,0);
  else if(index<nD+nU+nS) return G4QContent(0,0,1,0,0,0);
  else cerr<<"***G4QC::IndQ:index="<<index<<" for the QuarkContent="<<GetThis()<<endl;
  G4Exception("***G4QC::IndQ: Index exceeds the total number of quarks");
  return G4QContent(0,0,0,0,0,0);
}

// Return one-antiquark QC using index (a kind of iterator)
G4QContent G4QContent::IndAQ (G4int index)
{//        ==============================
#ifdef debug
  cout << "G4QC::IndAQ is called"<<endl;
#endif
  if(index<nAD) return G4QContent(0,0,0,1,0,0);
  else if(index<nAD+nAU) return G4QContent(0,0,0,0,1,0);
  else if(index<nAD+nAU+nAS) return G4QContent(0,0,0,0,0,1);
  else cerr<<"***G4QC::IndAQ:index="<<index<<" for the QuarkContent="<<GetThis()<<endl;
  G4Exception("***G4QC::IndAQ: Index exceeds the total number of antiquarks");
  return G4QContent(0,0,0,0,0,0);
}

// Reduce a fixet number (if<0:all?) of valence Q-Qbar pairs, returns a#of pairs to reduce more
G4int G4QContent::DecQAQ(const G4int& nQAQ)
{//   =====================================
#ifdef debug
  cout<<"DecQC: n="<<nQAQ<<", U="<<nU<<", D="<<nD<<", S="<<nS
      <<", AU="<<nAU<<", AD="<<nAD<<", AS="<<nAS<<endl;
#endif
  G4int ban = GetBaryonNumber();
  G4int tot = GetTot();    // Total number of quarks in QC
  if (tot==ban*3) return 0;// Nothing to reduce
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
  G4int nRet =nTotP-nQAQ;  // a#of additional pairs for reduction
  if (nQAQ<0)              // === Limited reduction case @@ not tuned for baryons !!
  {
    G4int res=tot+nQAQ;
#ifdef debug
	cout<<"DecQC: tot="<<tot<<", nTP="<<nTotP<<", res="<<res<<endl;
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
  cout << "DecQC: demanded "<<nQAQ<<" pairs, executed "<<nReal<<" pairs"<<endl;
#endif
  G4int nt = tot - nTotP - nTotP;
  for (int i=0; i<nReal; i++)
  {
    G4double base = nTotP;
    //if (nRet && nSP==1 && !nQAQ) base = nLP;                   // Keep S-Sbar pair if possible
    G4int j = static_cast<int>(base*G4UniformRand());            // Random integer "SortOfQuark"
    if (nUP && j<nUP && (nRet>2 || nUP>1 || (nD<2 && nS<2)))     // --- U-Ubar pair
	{
#ifdef debug
      cout << "DecQC: decrementing UAU pair UP="<<nUP<<",nU="<<nU<<",nAU="<<nAU<<endl;
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
      cout << "DecQC: decrementing DAD pair DP="<<nDP<<",nD="<<nD<<",nAD="<<nAD<<endl;
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
      cout << "DecQC: decrementing SAS pair SP="<<nSP<<",nS="<<nS<<",nAS="<<nAS<<endl;
#endif
      nS--;
      nAS--;
      nSP--;
      nTotP--;
	}
    else if (nUP)                                  // --- U-Ubar pair cancelation (final)
	{
#ifdef debug
      cout << "DecQC: decrementing UAU pair (final) UP="<<nUP<<",nU="<<nU<<",nAU="<<nAU<<endl;
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
      cout << "DecQC: decrementing DAD pair (final) DP="<<nDP<<",nD="<<nD<<",nAD="<<nAD<<endl;
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
      cout << "DecQC: decrementing SAS pair SP="<<nSP<<",nS="<<nS<<",nAS="<<nAS<<endl;
#endif
      nS--;
      nAS--;
      nSP--;
      nTotP--;
	}
    else cerr<<"***DecQC:i="<<i<<",j="<<j<<",D="<<nDP<<",U="<<nUP<<",S="<<nSP<<",T="<<nTotP
             <<",nRet="<<nRet<<", QC="<<GetThis()<<endl;
  }
#ifdef debug
  cout<<"DecQC:out U="<<nU<<",D="<<nD<<",S="<<nS<<",AU="<<nAU<<",AD="<<nAD<<",AS="<<nAS<<endl;
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
    cout<<"IncQC:out QC="<<GetThis()<<",j="<<j<<" for i="<<i<<endl;
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
  if(b3>0) return (b3-3*(rS+dQ))/6;                  // Matter case
  else     return 0;                                 // Meson  case or anti-matter case
}

// Calculate a#of neutrons in the QC (for nuclei)
G4int G4QContent::GetN() const
{//   ========================
  G4int rD=nD-nAD;                                   // Constituent d-quarks
  G4int rU=nU-nAU;                                   // Constituent u-quarks
  G4int rS=nS-nAS;                                   // Constituent s-quarks
  G4int dQ=rD-rU;                                    // Isotopic assimetry
  G4int b3=rD+rU+rS;                                 // (Baryon number) * 3
  if(b3>0) return (b3-3*(rS-dQ))/6;                  // Matter case
  else     return 0;                                 // Meson  case or anti-matter case
}

// Calculate a#of lambdas in the QC (for nuclei)
G4int G4QContent::GetL() const
{//   ========================
  G4int rS=nS-nAS;                                   // Constituent s-quarks
  G4int b3=nD-nAD+nU-nAU+rS;                         // (Baryon number) * 3
  if(b3>0) return rS;                                // Matter case
  else     return 0;                                 // Meson  case or anti-matter case
}

// Calculate a#of anti-protons in the QC (for anti-nuclei)
G4int G4QContent::GetAP() const
{//   =========================
  G4int rD=nAD-nD;                                   // Constituent anti-d-quarks
  G4int rU=nAU-nU;                                   // Constituent anti-u-quarks
  G4int rS=nAS-nS;                                   // Constituent anti-s-quarks
  G4int dQ=rD-rU;                                    // Isotopic assimetry
  G4int b3=rD+rU+rS;                                 // - (Baryon number) * 3
  if(b3>0) return (b3-3*(rS+dQ))/6;                  // Anti-matter case
  else     return 0;                                 // Meson  case or matter case
}

// Calculate a#of anti-neutrons in the QC (for anti-nuclei)
G4int G4QContent::GetAN() const
{//   =========================
  G4int rD=nAD-nD;                                   // Constituent anti-d-quarks
  G4int rU=nAU-nU;                                   // Constituent anti-u-quarks
  G4int rS=nAS-nS;                                   // Constituent anti-s-quarks
  G4int dQ=rD-rU;                                    // Isotopic assimetry
  G4int b3=rD+rU+rS;                                 // - (Baryon number) * 3
  if(b3>0) return (b3-3*(rS-dQ))/6;                  // Anti-matter case
  else     return 0;                                 // Meson  case or matter case
}

// Calculate a#of anti-lambdas in the QC (for anti-nuclei)
G4int G4QContent::GetAL() const
{//   =========================
  G4int rS=nAS-nS;                                   // Constituent anti-s-quarks
  G4int b3=nAD-nD+nAU-nU+rS;                         // - (Baryon number) * 3
  if(b3>0) return rS;                                // Anti-matter case
  else     return 0;                                 // Meson  case or matter case
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

  if(c%3)cerr<<"***G4QContent:GetCharge: c="<<c<<"/3 isn't integer, QC="<<GetThis()<<endl;

  return c/3;
}

// Calculate a Baryon Number for the QC
G4int G4QContent::GetBaryonNumber() const
{//   ===================================
  G4int b=nU+nD+nS-nAU-nAD-nAS;

  if(b%3) cerr<<"***G4Content: BaryonNumber="<<b<<"/3 is not an integer value"<<endl;

  return b/3;
}

// Gives the PDG of the lowest (in mass) S-wave hadron or Chipolino (=10) for double hadron
G4int G4QContent::GetSPDGCode() const
{//   ===============================
  static const G4int NUCPDG  = 90000000;
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
      cerr<<"***G4QC::SPDG: CancelationU U="<<mU<<",D="<<mD<<",S="<<mS<<",QC="<<GetThis()<<endl;
      return 0;
    }
  }
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
      cerr<<"***G4QC::SPDG: CancelationD U="<<mU<<",D="<<mD<<",S="<<mS<<",QC="<<GetThis()<<endl;
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
      cerr<<"***G4QC::SPDG: CancelationS U="<<mU<<",D="<<mD<<",S="<<mS<<",QC="<<GetThis()<<endl;
      return 0;
    }
  }

  G4int b=GetBaryonNumber();
  G4int c=GetCharge();
  G4int s=GetStrangeness();
#ifdef pdebug
  cout<<"G4QContent::SPDGC: before b="<<b<<",n="<<n<<",c="<<c<<",s="<<s<<",Q="<<GetThis()<<endl;
#endif
  if (b)                                         // ==================== Baryon case
  {
    G4int ab=abs(b);
    if(ab>=2 && n>=6)                            // Multi-Baryonium (NuclearFragment)
	{
      G4int mI=nU-nAU-nD+nAD;
      //if     (abs(mI)>3||mS>3||(b>0&&s<-1)||(b<0&&s>1)) return  0;
      //else if(abs(mI)>2||mS>2||(b>0&&s< 0)||(b<0&&s>0)) return 10;
      if     (b>0&&s==-1||b<0&&s==1) return 10;
      else if(abs(mI)>2||mS>2||b>0&&s< 0||b<0&&s>0) return GetZNSPDGCode();
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
    if(n>5)  return  0;                          //B+M+M Tripolino etc
    if(n==5) return 10;                          //B+M Chipolino
    if(mS>0)                                     // Strange Baryons
	{
      p=3002;
      if      (mS==3)            p+=332;         // Decuplet
      else if (mS==2)
      {
        if      (mU==1 && mD==0) p+=320;
        else if (mU==0 && mD==1) p+=310;
        else
		{
          cerr<<"***G4QContent::SPDG: Exotic BaryonSS with U="<<mU<<",D="<<mD<<",S="<<mS<<endl;
          return 0;
        }
	  }
      else if (mS==1)
      {
        if      (mU==2 && mD==0) p+=220;
        else if (mU==1 && mD==1) p+=120;         // Lambda (M_Lambda<M_Sigma0) PDG_Sigma=3212
        else if (mU==0 && mD==2) p+=110;
        else
        {
          cerr<<"***G4QContent::GetSPDG:Exotic BaryonS with U="<<mU<<",D="<<mD<<",S="<<mS<<endl;
          return 0;
        }
	  }
      else                                       // Superstrange case
      {
        //#ifdef debug
        cerr<<"***G4QContent::GetSPDG:ExoBaryonS?U="<<mU<<",D="<<mD<<",S="<<mS<<GetThis()<<endl;
        //#endif
        return 0;
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
        cerr<<"***G4QContent::GetSPDG: Exotic BaryonU with U="<<mU<<",D="<<mD<<",S="<<mS<<endl;
        return 0;
      }
	}
    else if (mD==3) p=1114;                      // Decuplet
    else
    {
      cerr<<"***G4QContent::GetSPDG: Exotic BaryonD with U="<<mU<<",D="<<mD<<",S="<<mS<<endl; 
      return 0;
	}
    if (b<0) p=-p;
  }
  else             // ====================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Meson case
  {
#ifdef debug
    cout<<"G4QCont::SPDG: mD="<<mD<<",mU="<<mU<<",mS="<<mS<<", b="<<b<<",c="<<c<<",s="<<s<<endl;
#endif
    if(n>4)                                      // Super Exotics
    {
#ifdef debug
      cerr<<"***G4QC::SPDG:n>4 SuperExo: U="<<mU<<",D="<<mD<<",S="<<mS<<",QC="<<GetThis()<<endl;
#endif
      return 0;
	}
    if(n==4) return 10;                          // M+M Chipolino
    if(abs(s)>1)
	{
      cerr<<"***G4QC::SPDG: Strangeness="<<s<<", QC="<<GetThis()<<" - Superstrange Meson"<<endl;
      return 0;
	}
    // Heavy quark should come first
    if(mS>0)                                     // Strange Mesons
	{
      p=301;
      if      (mS==2)
      {
        if (G4UniformRand()<0.333) p=221;        // eta
        else                       p+=30;        // eta'
      }
      else if (mU==1 && mD==0) p+=20;
      else if (mU==0 && mD==1) p+=10;
      else
      {
        cerr<<"*G4QC::SPDG:Exotic MesonS U="<<mU<<",D="<<mD<<",S="<<mS<<",QC="<<GetThis()<<endl;
        return 0;
      }
	}
    else if (mU>0)                               // Isotopic Mesons
	{
      p=201;
      if      (mU==2 && mD==0)
      {
        G4double rn=G4UniformRand();
        if     (rn<0.500) p=111;
        //else if(rn<0.667) p=331;
        else              p=221;
      }
      else if (mU==1 && mD==1) p+=10;
      else
	  {
        cerr<<"*G4QC::SPDG:Exotic MesonU U="<<mU<<",D="<<mD<<",S="<<mS<<",QC="<<GetThis()<<endl;
        return 0;
      }
	}
    else if (mD==2)
    {
      G4double rn=G4UniformRand();
      if     (rn<0.500) p=111;
      //else if(rn<0.667) p=331;
      else              p=221;
    }
    else
    {
#ifdef debug
      cerr<<"***G4QC::SPDG:Exotic MesonD U="<<mU<<",D="<<mD<<",S="<<mS<<",QC="<<GetThis()<<endl;
#endif
      return 0;
    }
    if (c<0 || (c==0 && mS>0 && s>0)) p=-p;
  }
#ifdef debug
  cout<<"G4QContent::GetSPDG: output SPDGcode="<<p<<" for the QuarkContent="<<GetThis()<<endl;
#endif
  return p;
}

// === Calculate a number of combinations of rhc out of lhc ==
G4int G4QContent::NOfCombinations(const G4QContent& rhs) const
{//   ========================================================
  G4int c=1; // Default number of combinations?

  G4int mD=rhs.GetD();
  if(mD>0)
  {
    int j=nD;
    if (j<=0) return 0;
    for (int i=1; i<=mD; i++)
	{
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  };

  G4int mU=rhs.GetU();
  if(mU>0)
  {
    int j=nU;
    if (j<=0) return 0;
    for (int i=1; i<=mU; i++)
	{
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  };

  G4int mS=rhs.GetS();
  if(mS>0)
  {
    int j=nS;
    if (j<=0) return 0;
    for (int i=1; i<=mS; i++)
	{
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  };

  G4int mAD=rhs.GetAD();
  if(mAD>0)
  {
    int j=nAD;
    if (j<=0) return 0;
    for (int i=1; i<=mAD; i++)
	{
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  };

  G4int mAU=rhs.GetAU();
  if(mAU>0)
  {
    int j=nAU;
    if (j<=0) return 0;
    for (int i=1; i<=mAU; i++)
	{
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  };

  G4int mAS=rhs.GetAS();
  if(mAS>0)
  {
    int j=nAS;
    if (j<=0) return 0;
    for (int i=1; i<=mAS; i++)
	{
      if(!j) return 0;
      c*=j/i;
      j--;
    }
  };
  
  return c;
} 











