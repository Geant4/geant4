// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4QHBook.cc,v 1.2 2000-09-11 09:03:53 mkossov Exp $
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

//#define sdebug

#include "G4QHBook.hh"

G4QHBook::G4QHBook() :
  nEvnt(0),
  histNevt( "histnevt.out", ios::out ),
  tuplEvtA( "tuplevta.out", ios::out ),
  tuplEvtQ( "tuplevtq.out", ios::out ),
  tuplIncl( "tuplincl.out", ios::out ),
  tuple3pi( "tuple3pi.out", ios::out )
{
#ifdef sdebug
  cout<<"G4QHBook::G4QHBook() - Start"<<endl;
#endif
  histNevt.setf( ios::scientific, ios::floatfield );
  tuplEvtA.setf( ios::scientific, ios::floatfield );
  tuplEvtQ.setf( ios::scientific, ios::floatfield );
  tuplIncl.setf( ios::scientific, ios::floatfield );
  tuple3pi.setf( ios::scientific, ios::floatfield );
#ifdef sdebug
  cout<<"G4QHBook::G4QHBook() - End"<<endl;
#endif
}

void G4QHBook::FillEvt(const G4QHadronVector* hadrons)
{
#ifdef sdebug
  cout<<"G4QHBook::FillEvt - Start"<<endl;
#endif
  nEvnt++;
  G4int tNH = hadrons->entries();
  G4int MtotD=tNH;
  G4int MtotR=tNH;
  // Correction on the number of d=-1 particles
  for (G4int ind=0; ind<tNH; ind++)
  {
	G4int    d = hadrons->at(ind)->GetNFragments();
	if(d==-1) {MtotD--; MtotR--;}
  }
  G4int MtotDm=MtotD;
  G4int MtotRm=MtotR;
#ifdef sdebug
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
    G4QHadron* curH= hadrons->at(ind);
    G4double m = curH->GetMass();
    G4int    c = curH->GetPDGCode();
	G4int    d = curH->GetNFragments();
	
	if(d>-1)
	{
	  if( picount<3 && 
        (c==111 || c==211 || c==-211 || c==311 || c==-311 || c==221 || c==331) )
	  {
 	    pdgm[picount]  = c;
	    lorVm[picount] = curH->Get4Momentum();
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
#ifdef sdebug
  cout<<"G4QHBook::FillEvt Multiplicities are collected Pi0= "<<Mpi0<<endl;
#endif
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

	tuple3pi<<nEvnt<<" "<<MtotD<<" "<<MtotR <<" "<<Mprot<<" "<<Mneut<<" "
	        <<Mdeut<<" "<<Mtrit<<" "<<MHe3  <<" "<<MHe4 <<" "<<Mgam <<" "
	        <<Mpip <<" "<<Mpim <<" "<<Mpi0  <<" "<<MKp  <<" "<<MK0  <<" "
	        <<MKm  <<" "<<MaK0 <<" "<<Meta  <<" "<<Metap<<" "<<Mrhom<<" "
	        <<Mrhop<<" "<<Mrho0<<" "<<Momega<<" "<<Mphi <<" "<<MKS0 <<" "
	        <<MKSC <<" "<<MaKS0<<" "<<MaKSC <<" "<<Mf2  <<" "<<Ma2m <<" "
	        <<Ma2p <<" "<<Ma20 <<" "<<Mf2p  <<" "
			<<m3pi <<" "<<m12  <<" "<<m13   <<" "<<m23  <<" "
		    <<pdgm[0]<<" "<<pdgm[1]<<" "<<pdgm[2]<<" "<< endl;
  }   

	tuplEvtA<<nEvnt<<" "<<MtotD<<" "<<MtotR <<" "<<Mprot<<" "<<Mneut<<" "
	        <<Mdeut<<" "<<Mtrit<<" "<<MHe3  <<" "<<MHe4 <<" "<<Mgam <<" "
	        <<Mpip <<" "<<Mpim <<" "<<Mpi0  <<" "<<MKp  <<" "<<MK0  <<" "
	        <<MKm  <<" "<<MaK0 <<" "<<Meta  <<" "<<Metap<<" "<<Mrhom<<" "
	        <<Mrhop<<" "<<Mrho0<<" "<<Momega<<" "<<Mphi <<" "<<MKS0 <<" "
	        <<MKSC <<" "<<MaKS0<<" "<<MaKSC <<" "<<Mf2  <<" "<<Ma2m <<" "
	        <<Ma2p <<" "<<Ma20 <<" "<<Mf2p  <<" "<< endl;
    
    for (ind=0; ind<tNH; ind++)
    {
      float xtupin[43];
      G4QHadron* cH= hadrons->at(ind);
      G4double m = cH->GetMass();
      G4LorentzVector lorV = cH->Get4Momentum();
      G4int c=cH->GetPDGCode();
      G4int d=cH->GetNFragments();
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
	  //if(pz/sqrt(p2)>0.984808) 
      tuplIncl<<nEvnt<<" "<<MtotD<<" "<<MtotR <<" "<<Mprot<<" "<<Mneut<<" "
              <<Mdeut<<" "<<Mtrit<<" "<<MHe3  <<" "<<MHe4 <<" "<<Mgam <<" "
              <<Mpip <<" "<<Mpim <<" "<<Mpi0  <<" "<<MKp  <<" "<<MK0  <<" "
              <<MKm  <<" "<<MaK0 <<" "<<Meta  <<" "<<Metap<<" "<<Mrhom<<" "
              <<Mrhop<<" "<<Mrho0<<" "<<Momega<<" "<<Mphi <<" "<<MKS0 <<" "
              <<MKSC <<" "<<MaKS0<<" "<<MaKSC <<" "<<Mf2  <<" "<<Ma2m <<" "
              <<Ma2p <<" "<<Ma20 <<" "<<Mf2p  <<" "<<  d  <<" "<<  c  <<" "
              << ns  <<" "<< nz  <<" "<< nn   <<" "<<  m  <<" "<< px  <<" "
              << py  <<" "<< pz  <<" "<<  e   <<" "<< endl;
    }
#ifdef sdebug
    cout<<"G4QHBook::FillEvt 4M are collected for all hadrons"<<endl;
#endif
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
      G4QHadron* cHd= hadrons->at(i);
      G4int       d = cHd->GetNFragments();
      G4int       m = d;
	  if(d==-1)   m = 0;
      G4int       c = cHd->GetPDGCode();
      G4double mass = cHd->GetMass();
#ifdef sdebug
	  cout<<"G4QHBook::FillEvt i="<<i<<",d="<<d<<",c="<<c<<",m="<<m<<",mass="<<mass<<endl;
#endif
      // Ignore all resonances with PDG ending with 5
      //if(c%10 != 5)
      //{
    	MtotR -= m;
    	if(m) MtotD--;
    	if(m) while(m)
    	{
    	  i++;
#ifdef sdebug
    	  cout<<"G4QHBook::FillEvt Before i="<<i<<",tNH="<<tNH<<endl;
#endif
          G4QHadron* cHm= hadrons->at(i);
    	  G4int n = cHm->GetNFragments();
#ifdef sdebug
    	  cout<<"G4QHBook::FillEvt m="<<m<<",n="<<n<<endl;
#endif
		  if(n==-1) n = 0;
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
#ifdef sdebug
	  cout<<"G4QHBook::FillEvt<L> tNH="<<tNH<<", i="<<i<<", m="<<m<<",d="<<d<<endl;
#endif
    }
	tuplEvtQ<<nEvnt<<" "<<MtotD<<" "<<MtotR <<" "<<Mprot<<" "<<Mneut<<" "
	        <<Mdeut<<" "<<Mtrit<<" "<<MHe3  <<" "<<MHe4 <<" "<<Mgam <<" "
	        <<Mpip <<" "<<Mpim <<" "<<Mpi0  <<" "<<MKp  <<" "<<MK0  <<" "
	        <<MKm  <<" "<<MaK0 <<" "<<Meta  <<" "<<Metap<<" "<<Mrhom<<" "
	        <<Mrhop<<" "<<Mrho0<<" "<<Momega<<" "<<Mphi <<" "<<MKS0 <<" "
	        <<MKSC <<" "<<MaKS0<<" "<<MaKSC <<" "<<Mf2  <<" "<<Ma2m <<" "
	        <<Ma2p <<" "<<Ma20 <<" "<<Mf2p  <<" "<< endl;
#ifdef sdebug
  cout<<"G4QHBook::FillEvt - End"<<endl;
#endif
}

G4QHBook::~G4QHBook() 
{
#ifdef sdebug
  cout<<"~G4QHBook(): Writing a number of events generated to histNevt"<<endl;
#endif
  histNevt<<nEvnt<<" "<< endl;
}




