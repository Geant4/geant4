// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QChipolino.cc,v 1.1 2000-08-17 13:55:49 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QChipolino ----------------
//             by Mikhail Kossov, Sept 1999.
//      class for Quasmon initiated Chipolinos generated by CHIPS Model
// ------------------------------------------------------------------

//#define debug
//#define pdebug

#include "G4QChipolino.hh"

G4QChipolino::G4QChipolino(G4QContent& QCont)
{
  // @@ Does not work as static const ??
  G4QPDGCode Pi0(111);
  G4double   mPi0  = Pi0.GetMass();
  G4QContent Pi0QC = Pi0.GetQuarkContent();
  G4int ban =QCont.GetBaryonNumber();
  G4int tban=abs(3*ban);
  G4int tot=QCont.GetTot();   // Initial total number of quarks in QC
  G4int tod=tot%2;            // tot is even - meson or dibaryon-nucleus
  if(!tod&&(tot<4||ban&&tot<tban)||tod&&tot<tban+2) QCont.IncQAQ(1,0.); // Add quark-pair
  G4QContent rQC=QCont;       // Copy for possible reduction ("annihilation" of q-qbar pairs)
  tot=rQC.GetTot();           // New total number of quarks in QC  (temporary)
  if   (tot%2)rQC.DecQAQ(-tban-2); // Reduce pairs, keep only 5 quarks  (baryon case)
  else if(ban)rQC.DecQAQ(-tban); // Reduce pairs, keep only 6 quarks  (dibaryon case)
  else        rQC.DecQAQ(-4); // Reduce pairs, keep only 4 quarks  (meson case)
  tot=rQC.GetTot();           // Final total number of quarks      (updated)
#ifdef debug
  cout<<"G4QChipolino is called with QC="<<QCont<<",rQC="<<rQC<<",tot="<<tot<<endl;
#endif
  minM=1000000.;              // Prototype of minimal mass         (@@ just a big number)
  G4double m1=mPi0;
  theQPDG1   = Pi0;
  G4double m2=mPi0;
  theQPDG2   = Pi0;
  theQCont1  = Pi0QC;
  if      (!tot)              // Should not be here, just in case     (strange input)
  {
    cerr<<"***G4QChipolino: shouldn't be here 1 QC="<<rQC<<endl;
    G4Exception("G4QChipolino cannot be constructed as n_tot=0");
  }
  else if (tot==2 || tot==3)  // Should not be here (basic octet/singlet states)
  {
    cerr<<"***G4QChipolino: shouldn't be here 2 QC="<<rQC<<endl;
    theQCont1= rQC;
    theQPDG1.InitByQCont(rQC);
    theQCont = rQC+Pi0QC;
  }
  else if (tot==4)                 // Two possible combinations for the meson
  {
    G4QContent bQC=rQC.IndQ();
#ifdef debug
    cout<<"G4QChipolino: tot=4,rQC="<<rQC<<",bQC="<<bQC<<endl;
#endif
    for(int j=0; j<2; j++)
    {
      G4QContent aQC=rQC.IndAQ(j);
      G4QContent cQC=bQC+aQC;
      G4QPDGCode cQPDG(cQC);
      G4double                    M1=cQPDG.GetMass(); 
      if(cQPDG.GetPDGCode()==221) M1=mPi0;
      G4QContent oQC=rQC-cQC;
#ifdef debug
	  cout<<"G4QChipolino: aQC="<<aQC<<", cQC="<<cQC<<", oQC="<<oQC<<endl;
#endif
      G4QPDGCode oQPDG(oQC);
      G4double                    M2=oQPDG.GetMass();
      if(oQPDG.GetPDGCode()==221) M2=mPi0;
      G4double m=M1+M2;
#ifdef debug
	  cout<<"G4QChipolino: c="<<cQPDG<<",cM="<<M1<<",o="<<oQPDG<<",oM="<<M2
          <<",cM+0M="<<m<<", curMinM="<<minM<<endl;
#endif
      if(m<minM)
      {
        minM=m;
        theQPDG1  = cQPDG;
        theQCont1 = cQC;
        theQPDG2  = oQPDG;
      }
    }
  }
  else if (tot==5)                 // Four possible combinations for the baryon
  {
    G4int nQ=rQC.GetQ();
    G4int nA=rQC.GetAQ();
    G4bool fl=nA>nQ;               // Flag of antibaryon case
#ifdef pdebug
	cout<<"G4QChipolino: Baryon case nQ="<<nQ<<",nA="<<nA<<",QC="<<rQC
        <<",fl="<<fl<<endl;
#endif
    G4QContent bQC;
    if (fl) bQC=rQC.IndQ();       // Antibaryon case
    else    bQC=rQC.IndAQ();      // Baryon case - QC of antiquark
    for (int i=0; i<4; i++)
    {
      G4QContent cQC;
      if (fl) cQC=bQC+rQC.IndAQ(i);
      else    cQC=bQC+rQC.IndQ(i);// Make mesonout of anti-quark
      G4QPDGCode cQPDG(cQC);      // Make QPDG particle
      G4double                    M1=cQPDG.GetMass(); // Get meson mass
      if(cQPDG.GetPDGCode()==221) M1=mPi0; // Make pi0 out of eta
      G4QContent oQC=rQC-cQC;     // Make residual baryon
      G4QPDGCode oQPDG(oQC);      // Make QPDG of residual baryon
      G4double                    M2=oQPDG.GetMass(); // Get baryon mass
      if(oQPDG.GetPDGCode()==221) M2=mPi0; // @@ Never !!
      G4double m=M1+M2;
      if(m<minM)
      {
        minM=m;
        theQPDG1  = cQPDG;
        theQCont1 = cQC;
        theQPDG2  = oQPDG;
      }
    }
#ifdef pdebug
	cout<<"G4QChipolino: Baryon case minM="<<minM<<", M="<<theQCont1<<theQPDG1
        <<", B="<<theQPDG2<<endl;
#endif
  }
  else if (tot==6)                 // Four possible combinations for the di-baryon
  {
    if(ban)
	{
      G4int nQ=rQC.GetQ();
      G4int nA=rQC.GetAQ();
      G4bool fl=nA>nQ;             // Flag of anti-dibaryon case
#ifdef debug
	  cout<<"G4QChipolino: Di-Baryon case nQ="<<nQ<<",nA="<<nA<<",QC="<<rQC<<",fl="<<fl<<endl;
#endif
      for (int i=0; i<4; i++)
      {
        G4QContent aQC;
        if (fl) aQC=rQC.IndAQ(i);
        else    aQC=rQC.IndQ(i);
        for (int j=i+1; j<5; j++)
        {
          G4QContent bQC;
          if (fl) bQC=aQC+rQC.IndAQ(j);
          else    bQC=aQC+rQC.IndQ(j);
          for (int k=j+1; k<6; k++)
          {
            G4QContent cQC;
            if (fl) cQC=bQC+rQC.IndAQ(k);
            else    cQC=bQC+rQC.IndQ(k);
            G4QPDGCode cQPDG(cQC);
            G4double                    M1=cQPDG.GetMass();
            if(cQPDG.GetPDGCode()==221) M1=mPi0;
            G4QContent oQC=rQC-cQC;
            G4QPDGCode oQPDG=(oQC);
            G4double                    M2=oQPDG.GetMass();
            if(oQPDG.GetPDGCode()==221) M2=mPi0;
            G4double m=M1+M2;
            if(m<minM)
            {
              minM=m;
              theQPDG1  = cQPDG;
              theQCont1 = cQC;
              theQPDG2  = oQPDG;
            }
		  }
	    }
      }
    }
    else                       // Baryon-AntiBaryon
	{
      theQCont1 = rQC.IndQ(0)+rQC.IndQ(1)+rQC.IndQ(2);
      theQPDG1.InitByQCont(theQCont1);
      theQPDG2.InitByQCont(rQC.IndAQ(0)+rQC.IndAQ(1)+rQC.IndAQ(2));
	}
  }
  else if(((rQC.GetU() )>(rQC.GetS() -4) && (rQC.GetD() )>(rQC.GetS() -4)) ||
          ((rQC.GetAU())>(rQC.GetAS()-4) && (rQC.GetAD())>(rQC.GetAS()-4)) )
  {
    G4int kD=rQC.GetD();
    G4int kU=rQC.GetU();
    G4int kS=rQC.GetS();
    G4int mD=rQC.GetAD();
    G4int mU=rQC.GetAU();
    G4int mS=rQC.GetAS();
    G4int nQ=rQC.GetQ();
    G4int nA=rQC.GetAQ();
    G4bool fl=nA>nQ;           // Flag of anti-fragment case
#ifdef debug
	cout<<"G4QChipolino: NucFragment case nQ="<<nQ<<",nAQ="<<nA<<", QC="<<rQC<<",fl="<<fl<<endl;
#endif
	if((fl&&kS>1)||(!fl&&mS>1))
    {
      cerr<<"***G4QChipolino: ***Overfowed by strange quarks*** rQC="<<rQC<<endl;
	  G4Exception("G4QChipolino: Nuclear Fragment-Chipolino is overflowed by strange quarks");
    }
    else if(fl)                // ===> Anti-fragment
	{
      //cerr<<"***G4QChipolino: ***Anti-nuclear fragments*** rQC="<<rQC<<endl;
	  //G4Exception("G4QChipolino: Antinuclear fragments are not yet supported");
      if(!mS)                                                           // No strange quarks
	  {
        G4int nI=mU-mD;                                                 // Isotopic shift
        G4int nN=(mU+mD-nI*3)/6;
        if(!kS)                                                         // No kaons
		{
          if((nI>=0&&nN>=0)||(nI<0&&nN>=-nI))                           // Delta isn't necessary
          {
            if(nI>0)                                                    // Excess of antiprotons
            {
              theQPDG1 = G4QPDGCode(-(90000000+1000*(nN+nI-1)+nN));     // A Fragment-AProton
              theQPDG2 = G4QPDGCode(-2212);                             // An Anti-Proton
            }
            else                                                        // Excess of a-neutrons
            {
              theQPDG1 = G4QPDGCode(-(90000000+1000*(nN+nI)+nN-1));     // A Fragment-ANeutron
              theQPDG2 = G4QPDGCode(-2112);                             // An Anti-Neutron
            }
	      }
          else if((nI>=0&&nN>-2)||(nI<0&&nN>-nI-2))                     // Delta can be a part
		  {
            if(nI>0)                                                    // Excess of au-quarks
            {
              theQPDG1=G4QPDGCode(-(90000000+1000*(nN+nI-2)+nN+1));     // A Fragment-AProton
              theQPDG2=G4QPDGCode(-2224);                               // An Anti-Delta++
            }
            else                                                        // Excess of ad-quarks
            {
              theQPDG1=G4QPDGCode(-(90000000+1000*(nN+nI+1)+nN-2));     // A Fragment-ANeutron
              theQPDG2=G4QPDGCode(-1114);                               // An Anti-Delta-
            }
		  }
          else
		  {
            cerr<<"***G4QChipolino: **A**Isotopic asymmetry (without S)*** rQC="<<rQC<<endl;
	        G4Exception("G4QChipolino: Exotic Isotopic asymmety of AntiMultyBaryon Quasmon");
		  }
        }
        else if(kS<2)                                                   // NucFrag+K is possible
		{
          G4int    nN =(mU+mD-4-nI*3)/6;
          if(nI>0)                                                      // Excess of au-quarks
          {
            nN+=1;
            theQPDG1 = G4QPDGCode(-(90000000+1000*(nN+nI-1)+nN));       // An Anti-Fragment
            theQPDG2 = G4QPDGCode(-321);                                // A K- meson
		  }
          else
		  {
            theQPDG1 = G4QPDGCode(-(90000000+1000*(nN+nI+1)+nN));       // An AntiFragment
            theQPDG2 = G4QPDGCode(-311);                                // An Anti-K0 meson
          }
		}
        else
	    {
          cerr<<"***G4QChipolino: ***Too many kaons are needed*** rQC="<<rQC<<endl;
	      G4Exception("G4QChipolino: Too much Kaons are needed together with AntiNuclearFragm");
		}
	  }
      else                     // Fragment with strangeness
	  {
        if(mS<=mU&&mS<=mD)     // Fragment consisting of Neutrons, Protons & Lambrdas only
		{
          G4int nI=mU-mD;                                               // Isotopic shift
          G4int nN=(mU+mD-mS-mS-nI*3)/6;
          if((nI>=0&&nN>=0)||(nI<0&&nN>=-nI))                           // Delta isn't necessary
          {
            if(nI>0)                                                    // Excess of protons
            {
              theQPDG1 = G4QPDGCode(-(90000000+1000*(kS*1000+nN+nI-1)+nN));// A Fragment-AProton
              theQPDG2 = G4QPDGCode(-2212);                             // An Anti-Proton
            }
            else                                                        // Excess of neutrons
            {
              theQPDG1 = G4QPDGCode(-(90000000+1000*(kS*1000+nN+nI)+nN-1));//A Fragment-ANeutron
              theQPDG2 = G4QPDGCode(-2112);                             // An Anti-Neutron
            }
		  }
          else if((nI>=0&&nN>-2)||(nI<0&&nN>-nI-2))                     // Delta can be a part
		  {
            if(nI>0)                                                    // Excess of au-quarks
            {
              theQPDG1=G4QPDGCode(-(90000000+1000*(kS*1000+nN+nI-2)+nN+1));// A Fragment-AProton
              theQPDG2=G4QPDGCode(-2224);                               // An Anti-Delta++
            }
            else                                                        // Excess of ad-quarks
            {
              theQPDG1=G4QPDGCode(-(90000000+1000*(kS*1000+nN+nI+1)+nN-2));//A Fragment-ANeutron
              theQPDG2=G4QPDGCode(-1114);                               // An Anti-Delta-
            }
		  }
          else
		  {
            cerr<<"***G4QChipolino: **A**Isotopic assimetry (with S)*** rQC="<<rQC<<endl;
	        G4Exception("G4QChipolino:Exotic Isotopics of Strange AntiMultyBaryon Quasmon");
		  }
		}
        else                                                            // Excess of s-quarks
		{
          G4int       lam=mU;                                           // A#of Anti-Lambdas
          if (lam>mD) lam=mD;
          G4int lD=mD-lam;                                              // Residual ad-quarks
          G4int lU=mU-lam;                                              // Residual au-quarks
          G4int lS=mS-lam;                                              // Residual as-quarks
          if(lD+lU+lS!=3||lD<0||lU<0||lS<0)
		  {
            cerr<<"***G4QChipolino:*AntiFragm* rQC="<<rQC<<",s="<<lS<<",u="<<lU<<",d"<<lD<<endl;
	        G4Exception("G4QChipolino: Exotic superstrange AntiMultyBaryon");
		  }
          if     ( !lD && lU==2) theQPDG2=G4QPDGCode(-3222);            // Anti-Sigma+
          else if( !lU && lD==2) theQPDG2=G4QPDGCode(-3112);            // Anti-Sigma-
          else if( !lD && lU==1) theQPDG2=G4QPDGCode(-3322);            // Anti-Ksi0
          else if( !lU && lD==1) theQPDG2=G4QPDGCode(-3312);            // Anti-Ksi-
          else                   theQPDG2=G4QPDGCode(-3334);            // Anti-Omega-
          theQPDG1=G4QPDGCode(-(90+lam)*1000000);                       // Anti Strange Matter
		}
        theQCont1  = rQC-theQPDG2.GetQuarkContent();                    // QCont of Fragment-H
        theQCont   = rQC;                                               // QCont of Chipolino
	  }
	}
	else                       // ===> Nuclear Fragment
	{
      if(!kS)                                                           // No strange quarks
	  {
        G4int nI=kU-kD;                                                 // Isotopic shift
        G4int nN=(kU+kD-nI*3)/6;
        if(!mS)                                                         // No kaons
		{
          if((nI>=0&&nN>=0)||(nI<0&&nN>=-nI))                           // Delta isn't necessary
          {
            if(nI>0)                                                    // Excess of protons
            {
              theQPDG1 = G4QPDGCode(90000000+1000*(nN+nI-1)+nN);        // A Fragment-Proton
              theQPDG2 = G4QPDGCode(2212);                              // A Proton
            }
            else                                                        // Excess of neutrons
            {
              theQPDG1 = G4QPDGCode(90000000+1000*(nN+nI)+nN-1);        // A Fragment-Neutron
              theQPDG2 = G4QPDGCode(2112);                              // A Neutron
            }
	      }
          else if((nI>=0&&nN>-2)||(nI<0&&nN>-nI-2))                     // Delta can be a part
		  {
            if(nI>0)                                                    // Excess of u-quarks
            {
              theQPDG1=G4QPDGCode(90000000+1000*(nN+nI-2)+nN+1);        // A Fragment-Proton
              theQPDG2=G4QPDGCode(2224);                                // A Delta++
            }
            else                                                        // Excess of d-quarks
            {
              theQPDG1=G4QPDGCode(90000000+1000*(nN+nI+1)+nN-2);        // A Fragment-Neutron
              theQPDG2=G4QPDGCode(1114);                                // A Delta-
            }
		  }
          else
		  {
            cerr<<"***G4QChipolino: ***Isotopic assimetry (without S)*** rQC="<<rQC<<endl;
	        G4Exception("G4QChipolino: Exotic Isotopic assimety of MultyBaryon Quasmon");
		  }
        }
        else if(mS<2)                                                   // NucFrag+K is possible
		{
          G4int    nN =(kU+kD-4-nI*3)/6;
          if(nI>0)                                                      // Excess of u-quarks
          {
            nN+=1;
            theQPDG1 = G4QPDGCode(90000000+1000*(nN+nI-1)+nN);          // A Fragment
            theQPDG2 = G4QPDGCode(321);                                 // A K+ meson
		  }
          else
		  {
            theQPDG1 = G4QPDGCode(90000000+1000*(nN+nI+1)+nN);            // A Fragment
            theQPDG2 = G4QPDGCode(311);                                 // A K0 meson
          }
		}
        else
	    {
          cerr<<"***G4QChipolino: ***Too many kaons are needed*** rQC="<<rQC<<endl;
	      G4Exception("G4QChipolino: More than one Kaon is needed together with Nuclear Fragm");
		}
	  }
      else                     // Fragment with strangeness
	  {
        if(kS<=kU&&kS<=kD)     // Fragment consisting of Neutrons, Protons & Lambrdas only
		{
          G4int nI=kU-kD;                                               // Isotopic shift
          G4int nN=(kU+kD-kS-kS-nI*3)/6;
          if((nI>=0&&nN>=0)||(nI<0&&nN>=-nI))                           // Delta isn't necessary
          {
            if(nI>0)                                                    // Excess of protons
            {
              theQPDG1 = G4QPDGCode(90000000+1000*(kS*1000+nN+nI-1)+nN);// A Fragment-Proton
              theQPDG2 = G4QPDGCode(2212);                              // A Proton
            }
            else                                                        // Excess of neutrons
            {
              theQPDG1 = G4QPDGCode(90000000+1000*(kS*1000+nN+nI)+nN-1);// A Fragment-Neutron
              theQPDG2 = G4QPDGCode(2112);                              // A Neutron
            }
		  }
          else if((nI>=0&&nN>-2)||(nI<0&&nN>-nI-2))                     // Delta can be a part
		  {
            if(nI>0)                                                    // Excess of u-quarks
            {
              theQPDG1=G4QPDGCode(90000000+1000*(kS*1000+nN+nI-2)+nN+1);// A Fragment-Proton
              theQPDG2=G4QPDGCode(2224);                                // A Delta++
            }
            else                                                        // Excess of d-quarks
            {
              theQPDG1=G4QPDGCode(90000000+1000*(kS*1000+nN+nI+1)+nN-2);// A Fragment-Neutron
              theQPDG2=G4QPDGCode(1114);                                // A Delta-
            }
		  }
          else
		  {
            cerr<<"***G4QChipolino: ***Isotopic assimetry (with S)*** rQC="<<rQC<<endl;
	        G4Exception("G4QChipolino:Exotic Isotopic Assimety of Strange MultyBaryon Quasmon");
		  }
		}
        else                                                            // Excess of s-quarks
		{
          G4int       lam=kU;                                           // A#of Lambda
          if (lam>kD) lam=kD;
          G4int lD=kD-lam;                                              // Residual d-quarks
          G4int lU=kU-lam;                                              // Residual u-quarks
          G4int lS=kS-lam;                                              // Residual s-quarks
          if(lD+lU+lS!=3||lD<0||lU<0||lS<0)
		  {
            cerr<<"***G4QChipolino: *Fragment* rQC="<<rQC<<",s="<<lS<<",u="<<lU<<",d"<<lD<<endl;
	        G4Exception("G4QChipolino: Exotic superstrange Multy Baryon");
		  }
          if     ( !lD && lU==2) theQPDG2=G4QPDGCode(3222);             // Sigma+
          else if( !lU && lD==2) theQPDG2=G4QPDGCode(3112);             // Sigma-
          else if( !lD && lU==1) theQPDG2=G4QPDGCode(3322);             // Ksi0
          else if( !lU && lD==1) theQPDG2=G4QPDGCode(3312);             // Ksi-
          else                   theQPDG2=G4QPDGCode(3334);             // Omega-
          theQPDG1=G4QPDGCode((90+lam)*1000000);                        // A Strange Matter
		}
        theQCont1  = rQC-theQPDG2.GetQuarkContent();                    // QCont of Fragment-H
        theQCont   = rQC;                                               // QCont of Chipolino
	  }
	}
  }
  else
  {
    cerr<<"***G4QChipolino: ***Exotics*** rQC="<<rQC<<endl;
	G4Exception("G4QChipolino: can not be constructed for the exotic baryon or meson");
  }
}


G4QChipolino::G4QChipolino(const G4QChipolino &right)
{
  theQPDG1  = right.theQPDG1;
  theQPDG2  = right.theQPDG2;
  theQCont  = right.theQCont;
  theQCont1 = right.theQCont1;
}

const G4QChipolino& G4QChipolino::operator=(const G4QChipolino &right)
{
  theQPDG1  = right.theQPDG1;
  theQPDG2  = right.theQPDG2;
  theQCont  = right.theQCont;
  theQCont1 = right.theQCont1;

  return *this;
}

G4QChipolino::~G4QChipolino() {}

// Standard output for G4QChipolino
ostream& operator<<(ostream& lhs, G4QChipolino& rhs)
{//      ===========================================
  lhs<<"{1="<<rhs.GetQPDG1()<<",2="<<rhs.GetQPDG2()<< "}";
  return lhs;
}






