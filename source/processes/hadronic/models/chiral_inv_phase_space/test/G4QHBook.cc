// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QHBook.cc,v 1.1 2000-08-17 14:17:14 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
//      ---------------- G4QHBook ----------------
//      by Mikhail Kossov and P. Degtiarenko, Nov 1999.
//      class for booking and filling histograms and ntuples
//      for the main CHIPStest routine
// ------------------------------------------------------------------

//#define pdebug

#include "G4QHBook.hh"

G4QHBook::G4QHBook() :
  nEvnt(0),
  histNevt(1,"N events generated",1,0.0,2.0),
  tupleEvtA(25,"EvtA",33),
  tupleEvtQ(27,"EvtQ",33),
  tupleIncl(20,"Incl",43),
  tuple3pi(22,"3pi",40)
{
#ifdef pdebug
  cout<<"G4QHBook::G4QHBook() - Start"<<endl;
#endif
  tuple3pi
    .setTag("Nevt").setTag("MtotD").setTag("MtotR")
    .setTag("Mprot").setTag("Mneut")
    .setTag("Mdeut").setTag("Mtrit").setTag("MHe3").setTag("MHe4")
    .setTag("Mgam").setTag("Mpim")
    .setTag("Mpip").setTag("Mpi0").setTag("MKp").setTag("MK0")
    .setTag("MKm").setTag("MaK0").setTag("Meta").setTag("Metap")
    .setTag("Mrhom").setTag("Mrhop").setTag("Mrho0").setTag("Momega")
    .setTag("Mphi").setTag("MKS0").setTag("MKSC")
    .setTag("MaKS0").setTag("MaKSC")
    .setTag("Mf2").setTag("Ma2m").setTag("Ma2p").setTag("Ma20").setTag("Mf2p")
    .setTag("m3pi").setTag("m12").setTag("m13").setTag("m23")
    .setTag("pdg1").setTag("pdg2").setTag("pdg3")
    .book();
  tupleIncl
    .setTag("Nevt").setTag("MtotD").setTag("MtotR")
    .setTag("Mprot").setTag("Mneut")
    .setTag("Mdeut").setTag("Mtrit").setTag("MHe3").setTag("MHe4")
    .setTag("Mgam").setTag("Mpim")
    .setTag("Mpip").setTag("Mpi0").setTag("MKp").setTag("MK0")
    .setTag("MKm").setTag("MaK0").setTag("Meta").setTag("Metap")
    .setTag("Mrhom").setTag("Mrhop").setTag("Mrho0").setTag("Momega")
    .setTag("Mphi").setTag("MKS0").setTag("MKSC")
    .setTag("MaKS0").setTag("MaKSC")
    .setTag("Mf2").setTag("Ma2m").setTag("Ma2p").setTag("Ma20").setTag("Mf2p")
    .setTag("ND").setTag("PDG").setTag("NS").setTag("NZ").setTag("NN")
    .setTag("m").setTag("Px").setTag("Py").setTag("Pz").setTag("E")
    .book();
  tupleEvtA
    .setTag("Nevt").setTag("MtotD").setTag("MtotR")
    .setTag("Mprot").setTag("Mneut")
    .setTag("Mdeut").setTag("Mtrit").setTag("MHe3").setTag("MHe4")
    .setTag("Mgam").setTag("Mpim")
    .setTag("Mpip").setTag("Mpi0").setTag("MKp").setTag("MK0")
    .setTag("MKm").setTag("MaK0").setTag("Meta").setTag("Metap")
    .setTag("Mrhom").setTag("Mrhop").setTag("Mrho0").setTag("Momega")
    .setTag("Mphi").setTag("MKS0").setTag("MKSC")
    .setTag("MaKS0").setTag("MaKSC")
    .setTag("Mf2").setTag("Ma2m").setTag("Ma2p").setTag("Ma20").setTag("Mf2p")
    .book();
  tupleEvtQ
    .setTag("Nevt").setTag("MtotD").setTag("MtotR")
    .setTag("Mprot").setTag("Mneut")
    .setTag("Mdeut").setTag("Mtrit").setTag("MHe3").setTag("MHe4")
    .setTag("Mgam").setTag("Mpim")
    .setTag("Mpip").setTag("Mpi0").setTag("MKp").setTag("MK0")
    .setTag("MKm").setTag("MaK0").setTag("Meta").setTag("Metap")
    .setTag("Mrhom").setTag("Mrhop").setTag("Mrho0").setTag("Momega")
    .setTag("Mphi").setTag("MKS0").setTag("MKSC")
    .setTag("MaKS0").setTag("MaKSC")
    .setTag("Mf2").setTag("Ma2m").setTag("Ma2p").setTag("Ma20").setTag("Mf2p")
    .book();
#ifdef pdebug
  cout<<"G4QHBook::G4QHBook() - End"<<endl;
#endif
}

void G4QHBook::FillEvt(const G4QHadronVector& hadrons)
{
#ifdef pdebug
  cout<<"G4QHBook::FillEvt - Start"<<endl;
#endif
  // Fill Histo No.1 - N events
  histNevt.fill(1.);

  nEvnt++;
  G4int tNH = hadrons.entries();
  G4int MtotD=tNH;
  G4int MtotR=tNH;
  // Correction on the number of d=-1 particles
  for (G4int ind=0; ind<tNH; ind++)
  {
	G4int    d = hadrons[ind]->GetNFragments();
	if(d==-1) {MtotD--; MtotR--;}
  }
  G4int MtotDm=MtotD;
  G4int MtotRm=MtotR;
#ifdef pdebug
  cout<<"G4QHBook::FillEvt - # of generated hadrons = "<<tNH<<",tD="<<MtotD<<endl;
#endif
  // Calculate multiplicities & fill ntuples for All tracks
  G4int Mprot=0;
  G4int Mneut=0;
  G4int Mdeut=0;
  G4int Mtrit=0;
  G4int MHe3=0;
  G4int MHe4=0;
  G4int Mgam=0;
  G4int Mpim=0;
  G4int Mpip=0;
  G4int Mpi0=0;
  G4int MKp=0;
  G4int MK0=0;
  G4int MKm=0;
  G4int MaK0=0;
  G4int Meta=0;
  G4int Metap=0;
  G4int Mrhom=0;
  G4int Mrhop=0;
  G4int Mrho0=0;
  G4int Momega=0;
  G4int Mphi=0;
  G4int MKS0=0;
  G4int MKSC=0; 
  G4int MaKS0=0;
  G4int MaKSC=0;
  G4int Mf2=0;
  G4int Ma2m=0;
  G4int Ma2p=0;
  G4int Ma20=0;
  G4int Mf2p=0;

  G4int picount=0;
  G4int pdgm[3];
  G4LorentzVector lorVm[3];

  for (ind=0; ind<tNH; ind++)
  {
    G4double m = hadrons[ind]->GetMass();
    G4int    c = hadrons[ind]->GetPDGCode();
	G4int    d = hadrons[ind]->GetNFragments();
	
	if(d>-1)
	{
	  if( picount<3 && 
        (c==111 || c==211 || c==-211 || c==311 || c==-311 || c==221 || c==331) )
	  {
 	    pdgm[picount]  = c;
	    lorVm[picount] = hadrons[ind]->Get4Momentum();
	    picount++;
	  }

	  MtotR         = MtotR - d;
	  if(d>0) MtotD = MtotD - 1;

	  if     (c==2212 || c==90001000)      Mprot++;
	  else if(c==2112 || c==90000001)      Mneut++;
	  else if(c==90001001)                 Mdeut++;
	  else if(c==90001002)                 Mtrit++;
	  else if(c==90002001)                 MHe3++;
	  else if(c==90002002)                 MHe4++;
	  else if(c== 111)      Mpi0++;
	  else if(c== 211)      Mpip++;
	  else if(c==-211)      Mpim++;
	  else if(c== 311)      MK0++;
	  else if(c== 321)      MKp++;
	  else if(c==-311)      MaK0++;
	  else if(c==-321)      MKm++;
	  else if(c== 22)       Mgam++;
	  else if(c== 221)      Meta++;
	  else if(c== 331)      Metap++;
	  else if(c== 113)      Mrho0++;
	  else if(c== 213)      Mrhop++;
	  else if(c==-213)      Mrhom++;
	  else if(c== 223)      Momega++;
	  else if(c== 333)      Mphi++;
	  else if(c== 313)      MKS0++;
	  else if(c== 323)      MKSC++;
	  else if(c==-313)      MaKS0++;
	  else if(c==-323)      MaKSC++;
	  else if(c== 215)      Ma2p++;
	  else if(c==-215)      Ma2m++;
	  else if(c== 115)      Ma20++;
	  else if(c== 225)      Mf2++;
	  else if(c== 335)      Mf2p++;
	}
  }
  if(MtotD==3 && ( (Mpi0+Mpip+Mpim)==3 ||
                   (Mpi0==1 && MK0==1 && MaK0==1) ||
                   (Mpi0==1 && Meta==1 && Metap==1)
                 )
    )
  {
	G4LorentzVector lorV12  = lorVm[0] + lorVm[1];
	G4LorentzVector lorV13  = lorVm[0] + lorVm[2];
	G4LorentzVector lorV23  = lorVm[1] + lorVm[2];
	G4LorentzVector lorV123 = lorVm[0] + lorV23;
	G4double m3pi = lorV123.m();
	G4double m12  = lorV12.m();
	G4double m13  = lorV13.m();
    G4double m23  = lorV23.m();

    float xtup3pi[40];
    G4int pidnt=0;
    xtup3pi[pidnt++] = nEvnt;
    xtup3pi[pidnt++] = MtotD;
    xtup3pi[pidnt++] = MtotR;
    xtup3pi[pidnt++] = Mprot;
    xtup3pi[pidnt++] = Mneut;
    xtup3pi[pidnt++] = Mdeut;
    xtup3pi[pidnt++] = Mtrit;
    xtup3pi[pidnt++] = MHe3;
    xtup3pi[pidnt++] = MHe4;
    xtup3pi[pidnt++] = Mgam;
    xtup3pi[pidnt++] = Mpip;
    xtup3pi[pidnt++] = Mpim;
    xtup3pi[pidnt++] = Mpi0;
    xtup3pi[pidnt++] = MKp;
    xtup3pi[pidnt++] = MK0;
    xtup3pi[pidnt++] = MKm;
    xtup3pi[pidnt++] = MaK0;
    xtup3pi[pidnt++] = Meta;
    xtup3pi[pidnt++] = Metap;
    xtup3pi[pidnt++] = Mrhom;
    xtup3pi[pidnt++] = Mrhop;
    xtup3pi[pidnt++] = Mrho0;
    xtup3pi[pidnt++] = Momega;
    xtup3pi[pidnt++] = Mphi;
    xtup3pi[pidnt++] = MKS0;
    xtup3pi[pidnt++] = MKSC;
    xtup3pi[pidnt++] = MaKS0;
    xtup3pi[pidnt++] = MaKSC;
    xtup3pi[pidnt++] = Mf2;
    xtup3pi[pidnt++] = Ma2m;
    xtup3pi[pidnt++] = Ma2p;
    xtup3pi[pidnt++] = Ma20;
    xtup3pi[pidnt++] = Mf2p;
    xtup3pi[pidnt++] = m3pi;
    xtup3pi[pidnt++] = m12;
    xtup3pi[pidnt++] = m13;
    xtup3pi[pidnt++] = m23;
    xtup3pi[pidnt++] = pdgm[0];
    xtup3pi[pidnt++] = pdgm[1];
    xtup3pi[pidnt++] = pdgm[2];
	tuple3pi.fill(xtup3pi);
  }   

    float xtupev[33];
    G4int evdnt=0;
    xtupev[evdnt++] = nEvnt;
    xtupev[evdnt++] = MtotD;
    xtupev[evdnt++] = MtotR;
    xtupev[evdnt++] = Mprot;
    xtupev[evdnt++] = Mneut;
    xtupev[evdnt++] = Mdeut;
    xtupev[evdnt++] = Mtrit;
    xtupev[evdnt++] = MHe3;
    xtupev[evdnt++] = MHe4;
    xtupev[evdnt++] = Mgam;
    xtupev[evdnt++] = Mpip;
    xtupev[evdnt++] = Mpim;
    xtupev[evdnt++] = Mpi0;
    xtupev[evdnt++] = MKp;
    xtupev[evdnt++] = MK0;
    xtupev[evdnt++] = MKm;
    xtupev[evdnt++] = MaK0;
    xtupev[evdnt++] = Meta;
    xtupev[evdnt++] = Metap;
    xtupev[evdnt++] = Mrhom;
    xtupev[evdnt++] = Mrhop;
    xtupev[evdnt++] = Mrho0;
    xtupev[evdnt++] = Momega;
    xtupev[evdnt++] = Mphi;
    xtupev[evdnt++] = MKS0;
    xtupev[evdnt++] = MKSC;
    xtupev[evdnt++] = MaKS0;
    xtupev[evdnt++] = MaKSC;
    xtupev[evdnt++] = Mf2;
    xtupev[evdnt++] = Ma2m;
    xtupev[evdnt++] = Ma2p;
    xtupev[evdnt++] = Ma20;
    xtupev[evdnt++] = Mf2p;
    tupleEvtA.fill(xtupev);
    
    for (ind=0; ind<tNH; ind++)
    {
    	float xtupin[43];
    	G4double m = hadrons[ind]->GetMass();
    	G4LorentzVector lorV = hadrons[ind]->Get4Momentum();
    	G4int c=hadrons[ind]->GetPDGCode();
    	G4int d=hadrons[ind]->GetNFragments();
    	G4int ns=0;
    	G4int nz=0;
    	G4int nn=0;
		if(c>90000000)
		{
		  ns=(c/1000000)%10;
		  nz=(c/1000)%1000;
		  nn=c%1000;
		}
    	G4double px = lorV.x();
    	G4double py = lorV.y();
    	G4double pz = lorV.z();
		G4double pt2 = px*px+py*py;
		G4double p2  = pt2+pz*pz;
    	G4double e  = lorV.t();
    	G4int indnt=0;
    	xtupin[indnt++] = nEvnt;
    	xtupin[indnt++] = MtotD;
    	xtupin[indnt++] = MtotR;
        xtupin[indnt++] = Mprot;
        xtupin[indnt++] = Mneut;
        xtupin[indnt++] = Mdeut;
        xtupin[indnt++] = Mtrit;
        xtupin[indnt++] = MHe3;
        xtupin[indnt++] = MHe4;
    	xtupin[indnt++] = Mgam;
    	xtupin[indnt++] = Mpip;
    	xtupin[indnt++] = Mpim;
    	xtupin[indnt++] = Mpi0;
    	xtupin[indnt++] = MKp;
    	xtupin[indnt++] = MK0;
    	xtupin[indnt++] = MKm;
    	xtupin[indnt++] = MaK0;
    	xtupin[indnt++] = Meta;
    	xtupin[indnt++] = Metap;
    	xtupin[indnt++] = Mrhom;
    	xtupin[indnt++] = Mrhop;
    	xtupin[indnt++] = Mrho0;
    	xtupin[indnt++] = Momega;
    	xtupin[indnt++] = Mphi;
    	xtupin[indnt++] = MKS0;
    	xtupin[indnt++] = MKSC;
    	xtupin[indnt++] = MaKS0;
    	xtupin[indnt++] = MaKSC;
    	xtupin[indnt++] = Mf2;
    	xtupin[indnt++] = Ma2m;
    	xtupin[indnt++] = Ma2p;
    	xtupin[indnt++] = Ma20;
    	xtupin[indnt++] = Mf2p;
    	xtupin[indnt++] = d;
    	xtupin[indnt++] = c;
    	xtupin[indnt++] = ns;
    	xtupin[indnt++] = nz;
    	xtupin[indnt++] = nn;
    	xtupin[indnt++] = m;
    	xtupin[indnt++] = px;
    	xtupin[indnt++] = py;
    	xtupin[indnt++] = pz;
    	xtupin[indnt++] = e;
		//if(pz/sqrt(p2)>0.984808) tupleIncl.fill(xtupin);
		tupleIncl.fill(xtupin);
    }
    
    // Calculate multiplicities & fill ntuple for Q (Level 1) tracks
    MtotD=tNH;
    MtotR=tNH;
    //MtotD=MtotDm;
    //MtotR=MtotRm;
    Mprot=0;
    Mneut=0;
    Mdeut=0;
    Mtrit=0;
    MHe3=0;
    MHe4=0;
    Mgam=0;
    Mpim=0;
    Mpip=0;
    Mpi0=0;
    MKp=0;
    MK0=0;
    MKm=0;
    MaK0=0;
    Meta=0;
    Metap=0;
    Mrhom=0;
    Mrhop=0;
    Mrho0=0;
    Momega=0;
    Mphi=0;
    MKS0=0;
    MKSC=0; 
    MaKS0=0;
    MaKSC=0;
    Mf2=0;
    Ma2m=0;
    Ma2p=0;
    Ma20=0;
    Mf2p=0;
    
    G4int k=0;
    G4int i=0;
    while (i<tNH)
    {
    	G4int       d = hadrons[i]->GetNFragments();
    	G4int       m = d;
		if(d==-1)   m = 0;
    	G4int       c = hadrons[i]->GetPDGCode();
    	G4double mass = hadrons[i]->GetMass();
    
    	// Ignore all resonances with PDG ending with 5
    	//if(c%10 != 5)
    	//{
    	  MtotR -= m;
    	  if(m) MtotD--;
    
    	  if(m) while(m)
    	  {
    		i++;
    		G4int n = hadrons[i]->GetNFragments();
		    if(n==-1)   n = 0;
    		if(n) 
    	    {
    		  MtotR -= n;
    		  MtotD--;	
    		  m+=n;
    		}
    		m--;
    	  }
    
		  if(d>-1)
    	  {
			if     (c==2212 || c==90001000)      Mprot++;
	        else if(c==2112 || c==90000001)      Mneut++;
	        else if(c==90001001)                 Mdeut++;
	        else if(c==90001002)                 Mtrit++;
	        else if(c==90002001)                 MHe3++;
	        else if(c==90002002)                 MHe4++;
	        else if(c== 111)      Mpi0++;
            else if(c== 211)      Mpip++;
            else if(c==-211)      Mpim++;
            else if(c== 311)      MK0++;
            else if(c== 321)      MKp++;
            else if(c==-311)      MaK0++;
            else if(c==-321)      MKm++;
    	    else if(c== 22)       Mgam++;
    	    else if(c== 221)      Meta++;
    	    else if(c== 331)      Metap++;
    	    else if(c== 113)      Mrho0++;
    	    else if(c== 213)      Mrhop++;
    	    else if(c==-213)      Mrhom++;
    	    else if(c== 223)      Momega++;
    	    else if(c== 333)      Mphi++;
    	    else if(c== 313)      MKS0++;
    	    else if(c== 323)      MKSC++;
    	    else if(c==-313)      MaKS0++;
    	    else if(c==-323)      MaKSC++;	
    	    else if(c== 215)      Ma2p++;
    	    else if(c==-215)      Ma2m++;
    	    else if(c== 115)      Ma20++;
    	    else if(c== 225)      Mf2++;
    	    else if(c== 335)      Mf2p++;
		  }
    	//}
    	//else
    	//{
    	//  MtotR--;
    	//  MtotD--;
    	//}
    	i++;
    }
    evdnt=0;
    xtupev[evdnt++] = nEvnt;
    xtupev[evdnt++] = MtotD;
    xtupev[evdnt++] = MtotR;
    xtupev[evdnt++] = Mprot;
    xtupev[evdnt++] = Mneut;
    xtupev[evdnt++] = Mdeut;
    xtupev[evdnt++] = Mtrit;
    xtupev[evdnt++] = MHe3;
    xtupev[evdnt++] = MHe4;
    xtupev[evdnt++] = Mgam;
    xtupev[evdnt++] = Mpip;
    xtupev[evdnt++] = Mpim;
    xtupev[evdnt++] = Mpi0;
    xtupev[evdnt++] = MKp;
    xtupev[evdnt++] = MK0;
    xtupev[evdnt++] = MKm;
    xtupev[evdnt++] = MaK0;
    xtupev[evdnt++] = Meta;
    xtupev[evdnt++] = Metap;
    xtupev[evdnt++] = Mrhom;
    xtupev[evdnt++] = Mrhop;
    xtupev[evdnt++] = Mrho0;
    xtupev[evdnt++] = Momega;
    xtupev[evdnt++] = Mphi;
    xtupev[evdnt++] = MKS0;
    xtupev[evdnt++] = MKSC;
    xtupev[evdnt++] = MaKS0;
    xtupev[evdnt++] = MaKSC;
    xtupev[evdnt++] = Mf2;
    xtupev[evdnt++] = Ma2m;
    xtupev[evdnt++] = Ma2p;
    xtupev[evdnt++] = Ma20;
    xtupev[evdnt++] = Mf2p;
    tupleEvtQ.fill(xtupev);
	//}
#ifdef pdebug
  cout<<"G4QHBook::FillEvt - End"<<endl;
#endif
}

G4QHBook::~G4QHBook() {}




