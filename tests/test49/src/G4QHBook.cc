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
// $Id: G4QHBook.cc,v 1.2 2009-08-27 14:13:54 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QHBook ----------------
//      made for process level tests by Mikhail Kossov - Feb 2005
//      class for booking and filling histograms and ntuples
// ------------------------------------------------------------------
//#define pdebug

#include "G4QHBook.hh"

G4QHBook::G4QHBook() :
  nEvnt(0),
  histNevt( "histnevt.out", std::ios::out ),
  tuplEvtA( "tuplevta.out", std::ios::out ),
  tuplEvtQ( "tuplevtq.out", std::ios::out ),
  tuplIncl( "tuplincl.out", std::ios::out ),
  tuple3pi( "tuple3pi.out", std::ios::out )
{
#ifdef pdebug
  G4cout<<"G4QHBook::G4QHBook() - Start"<<G4endl;
#endif
  histNevt.setf( std::ios::scientific, std::ios::floatfield );
  tuplEvtA.setf( std::ios::scientific, std::ios::floatfield );
  tuplEvtQ.setf( std::ios::scientific, std::ios::floatfield );
  tuplIncl.setf( std::ios::scientific, std::ios::floatfield );
  tuple3pi.setf( std::ios::scientific, std::ios::floatfield );
#ifdef pdebug
  G4cout<<"G4QHBook::G4QHBook() - End"<<G4endl;
#endif
}

//G4QHBook::G4QHBook() :
//  nEvnt(0),
//  histNevt(1,"N events generated",1,0.0,2.0),
//  tupleEvtA(25,"EvtA",33),
//  tupleEvtQ(27,"EvtQ",33),
//  tupleIncl(20,"Incl",43),
//  tuple3pi(22,"3pi",40)
//{
//#ifdef pdebug
//  G4cout<<"G4QHBook::G4QHBook() - Start"<<G4endl;
//#endif
//  tuple3pi
//    .setTag("Nevt").setTag("MtotD").setTag("MtotR")
//    .setTag("Mprot").setTag("Mneut")
//    .setTag("Mdeut").setTag("Mtrit").setTag("MHe3").setTag("MHe4")
//    .setTag("Mgam").setTag("Mpim")
//    .setTag("Mpip").setTag("Mpi0").setTag("MKp").setTag("MK0")
//    .setTag("MKm").setTag("MaK0").setTag("Meta").setTag("Metap")
//    .setTag("Mrhom").setTag("Mrhop").setTag("Mrho0").setTag("Momega")
//    .setTag("Mphi").setTag("MKS0").setTag("MKSC")
//    .setTag("MaKS0").setTag("MaKSC")
//    .setTag("Mf2").setTag("Ma2m").setTag("Ma2p").setTag("Ma20").setTag("Mf2p")
//    .setTag("m3pi").setTag("m12").setTag("m13").setTag("m23")
//    .setTag("pdg1").setTag("pdg2").setTag("pdg3")
//    .book();
//  tupleIncl
//    .setTag("Nevt").setTag("MtotD").setTag("MtotR")
//    .setTag("Mprot").setTag("Mneut")
//    .setTag("Mdeut").setTag("Mtrit").setTag("MHe3").setTag("MHe4")
//    .setTag("Mgam").setTag("Mpim")
//    .setTag("Mpip").setTag("Mpi0").setTag("MKp").setTag("MK0")
//    .setTag("MKm").setTag("MaK0").setTag("Meta").setTag("Metap")
//    .setTag("Mrhom").setTag("Mrhop").setTag("Mrho0").setTag("Momega")
//    .setTag("Mphi").setTag("MKS0").setTag("MKSC")
//    .setTag("MaKS0").setTag("MaKSC")
//    .setTag("Mf2").setTag("Ma2m").setTag("Ma2p").setTag("Ma20").setTag("Mf2p")
//    .setTag("ND").setTag("PDG").setTag("NS").setTag("NZ").setTag("NN")
//    .setTag("m").setTag("Px").setTag("Py").setTag("Pz").setTag("E")
//    .book();
//  tupleEvtA
//    .setTag("Nevt").setTag("MtotD").setTag("MtotR")
//    .setTag("Mprot").setTag("Mneut")
//    .setTag("Mdeut").setTag("Mtrit").setTag("MHe3").setTag("MHe4")
//    .setTag("Mgam").setTag("Mpim")
//    .setTag("Mpip").setTag("Mpi0").setTag("MKp").setTag("MK0")
//    .setTag("MKm").setTag("MaK0").setTag("Meta").setTag("Metap")
//    .setTag("Mrhom").setTag("Mrhop").setTag("Mrho0").setTag("Momega")
//    .setTag("Mphi").setTag("MKS0").setTag("MKSC")
//    .setTag("MaKS0").setTag("MaKSC")
//    .setTag("Mf2").setTag("Ma2m").setTag("Ma2p").setTag("Ma20").setTag("Mf2p")
//    .book();
//  tupleEvtQ
//    .setTag("Nevt").setTag("MtotD").setTag("MtotR")
//    .setTag("Mprot").setTag("Mneut")
//    .setTag("Mdeut").setTag("Mtrit").setTag("MHe3").setTag("MHe4")
//    .setTag("Mgam").setTag("Mpim")
//    .setTag("Mpip").setTag("Mpi0").setTag("MKp").setTag("MK0")
//    .setTag("MKm").setTag("MaK0").setTag("Meta").setTag("Metap")
//    .setTag("Mrhom").setTag("Mrhop").setTag("Mrho0").setTag("Momega")
//    .setTag("Mphi").setTag("MKS0").setTag("MKSC")
//    .setTag("MaKS0").setTag("MaKSC")
//    .setTag("Mf2").setTag("Ma2m").setTag("Ma2p").setTag("Ma20").setTag("Mf2p")
//    .book();
//#ifdef pdebug
//  G4cout<<"G4QHBook::G4QHBook() - End"<<G4endl;
//#endif
//}

void G4QHBook::FillEvt(const G4VParticleChange* hadrons, const G4DynamicParticle* pdProj)
{
  // *** For ntuple case *** Fill Histo No.1 - N events
  //histNevt.fill(1.);
  nEvnt++;
#ifdef pdebug
  G4cout<<"G4QHBook::FillEvt - Start, nEvnt="<<nEvnt<<G4endl;
#endif
  G4int tNH = hadrons->GetNumberOfSecondaries();
  G4TrackStatus lead = hadrons->GetTrackStatus();
  G4int begi=0;
  if(lead==fAlive) begi=-1;   // Means that leading particle is alive
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

  const G4DynamicParticle* cHd = 0;
  G4ParticleDefinition* pd = 0;

  for (G4int jnd=begi; jnd<tNH; jnd++)
  {
    if(jnd<0) pd  = pdProj->GetDefinition();
    else
    {
      cHd = hadrons->GetSecondary(jnd)->GetDynamicParticle();
      pd  = cHd->GetDefinition();
    }
    G4int c = pd->GetPDGEncoding();
    if(!c)
    {
      G4int chrg=static_cast<G4int>(pd->GetPDGCharge());
      G4int bary=static_cast<G4int>(pd->GetBaryonNumber());
      //c=90000000+chrg*999+bary; // OLD CHIPS notation
      c=1000000000+chrg*10000+bary*10; // New PDG 2006 definition
    }
	   if(picount<3 && (c==111 || c==211 || c==-211 || c==311 || c==-311 || c==221 || c==331))
	   {
 	    pdgm[picount]  = c;
      G4double mass  = pd->GetPDGMass();                      // Mass of the particle
      G4ThreeVector mom;        // Direction of the momentum
      G4double ener;            // Kinetic energy
      if(jnd<0)
      {
        mom = pdProj->GetMomentumDirection();        // Direction of the momentum
        ener = pdProj->GetKineticEnergy();                // Kinetic energy
      }
      else
      {
        mom = cHd->GetMomentumDirection();        // Direction of the momentum
        ener = cHd->GetKineticEnergy();                // Kinetic energy
      }
      G4double ten = ener + mass;                             // Total energy
	     G4double p = std::sqrt(ener*(ten + mass));                   // Abs value of the momentum
	     mom *= p;                                               // 3-momentum
	     lorVm[picount] =  G4LorentzVector(mom, ten);
	     picount++;
	   }
    // count number of particles
	   if     (c==2212 || c==1000010010)      Mprot++;
	   else if(c==2112 || c==1000000010)      Mneut++;
	   else if(c==1000010020)                 Mdeut++;
	   else if(c==1000010030)                 Mtrit++;
	   else if(c==1000020030)                 MHe3++;
	   else if(c==1000020040)                 MHe4++;
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
  G4int fNH=tNH-begi;
  // Fill the output file for the meson triplets (necessary for P-antiP) simulation only
  if(fNH==3 && ( (Mpi0+Mpip+Mpim)==3 ||
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

    //float xtup3pi[40];
    //G4int pidnt=0;
    //xtup3pi[pidnt++] = nEvnt;
    //xtup3pi[pidnt++] = fNH;
    //xtup3pi[pidnt++] = fNH;
    //xtup3pi[pidnt++] = Mprot;
    //xtup3pi[pidnt++] = Mneut;
    //xtup3pi[pidnt++] = Mdeut;
    //xtup3pi[pidnt++] = Mtrit;
    //xtup3pi[pidnt++] = MHe3;
    //xtup3pi[pidnt++] = MHe4;
    //xtup3pi[pidnt++] = Mgam;
    //xtup3pi[pidnt++] = Mpip;
    //xtup3pi[pidnt++] = Mpim;
    //xtup3pi[pidnt++] = Mpi0;
    //xtup3pi[pidnt++] = MKp;
    //xtup3pi[pidnt++] = MK0;
    //xtup3pi[pidnt++] = MKm;
    //xtup3pi[pidnt++] = MaK0;
    //xtup3pi[pidnt++] = Meta;
    //xtup3pi[pidnt++] = Metap;
    //xtup3pi[pidnt++] = Mrhom;
    //xtup3pi[pidnt++] = Mrhop;
    //xtup3pi[pidnt++] = Mrho0;
    //xtup3pi[pidnt++] = Momega;
    //xtup3pi[pidnt++] = Mphi;
    //xtup3pi[pidnt++] = MKS0;
    //xtup3pi[pidnt++] = MKSC;
    //xtup3pi[pidnt++] = MaKS0;
    //xtup3pi[pidnt++] = MaKSC;
    //xtup3pi[pidnt++] = Mf2;
    //xtup3pi[pidnt++] = Ma2m;
    //xtup3pi[pidnt++] = Ma2p;
    //xtup3pi[pidnt++] = Ma20;
    //xtup3pi[pidnt++] = Mf2p;
    //xtup3pi[pidnt++] = pdgm[0];
    //xtup3pi[pidnt++] = pdgm[1];
    //xtup3pi[pidnt++] = pdgm[2];
	   //tuple3pi.fill(xtup3pi);
	   tuple3pi<<nEvnt<<" "<<fNH  <<" "<<fNH   <<" "<<Mprot<<" "<<Mneut<<" "
	           <<Mdeut<<" "<<Mtrit<<" "<<MHe3  <<" "<<MHe4 <<" "<<Mgam <<" "
	           <<Mpip <<" "<<Mpim <<" "<<Mpi0  <<" "<<MKp  <<" "<<MK0  <<" "
	           <<MKm  <<" "<<MaK0 <<" "<<Meta  <<" "<<Metap<<" "<<Mrhom<<" "
	           <<Mrhop<<" "<<Mrho0<<" "<<Momega<<" "<<Mphi <<" "<<MKS0 <<" "
	           <<MKSC <<" "<<MaKS0<<" "<<MaKSC <<" "<<Mf2  <<" "<<Ma2m <<" "
	           <<Ma2p <<" "<<Ma20 <<" "<<Mf2p  <<" "
			         <<m3pi <<" "<<m12  <<" "<<m13   <<" "<<m23  <<" "
		          <<pdgm[0]<<" "<<pdgm[1]<<" "<<pdgm[2]<<" "<< G4endl;
  }   
  // *** For Ntuple case ***
  //float xtupev[33];
  //G4int evdnt=0;
  //xtupev[evdnt++] = nEvnt;
  //xtupev[evdnt++] = fNH;
  //xtupev[evdnt++] = fNH;
  //xtupev[evdnt++] = Mprot;
  //xtupev[evdnt++] = Mneut;
  //xtupev[evdnt++] = Mdeut;
  //xtupev[evdnt++] = Mtrit;
  //xtupev[evdnt++] = MHe3;
  //xtupev[evdnt++] = MHe4;
  //xtupev[evdnt++] = Mgam;
  //xtupev[evdnt++] = Mpip;
  //xtupev[evdnt++] = Mpim;
  //xtupev[evdnt++] = Mpi0;
  //xtupev[evdnt++] = MKp;
  //xtupev[evdnt++] = MK0;
  //xtupev[evdnt++] = MKm;
  //xtupev[evdnt++] = MaK0;
  //xtupev[evdnt++] = Meta;
  //xtupev[evdnt++] = Metap;
  //xtupev[evdnt++] = Mrhom;
  //xtupev[evdnt++] = Mrhop;
  //xtupev[evdnt++] = Mrho0;
  //xtupev[evdnt++] = Momega;
  //xtupev[evdnt++] = Mphi;
  //xtupev[evdnt++] = MKS0;
  //xtupev[evdnt++] = MKSC;
  //xtupev[evdnt++] = MaKS0;
  //xtupev[evdnt++] = MaKSC;
  //xtupev[evdnt++] = Mf2;
  //xtupev[evdnt++] = Ma2m;
  //xtupev[evdnt++] = Ma2p;
  //xtupev[evdnt++] = Ma20;
  //xtupev[evdnt++] = Mf2p;
  //tupleEvtA.fill(xtupev);
  tuplEvtA<<nEvnt<<" "<<fNH  <<" "<<fNH   <<" "<<Mprot<<" "<<Mneut<<" "
          <<Mdeut<<" "<<Mtrit<<" "<<MHe3  <<" "<<MHe4 <<" "<<Mgam <<" "
          <<Mpip <<" "<<Mpim <<" "<<Mpi0  <<" "<<MKp  <<" "<<MK0  <<" "
          <<MKm  <<" "<<MaK0 <<" "<<Meta  <<" "<<Metap<<" "<<Mrhom<<" "
          <<Mrhop<<" "<<Mrho0<<" "<<Momega<<" "<<Mphi <<" "<<MKS0 <<" "
          <<MKSC <<" "<<MaKS0<<" "<<MaKSC <<" "<<Mf2  <<" "<<Ma2m <<" "
          <<Ma2p <<" "<<Ma20 <<" "<<Mf2p  <<" "<< G4endl;
    
  for (G4int knd=0; knd<tNH; knd++)
  {
    //float xtupin[43];
    G4int d=0; // a fake number to support the universal way of writing
    if(knd<0) pd  = pdProj->GetDefinition();
    else
    {
      cHd = hadrons->GetSecondary(knd)->GetDynamicParticle();
      pd  = cHd->GetDefinition();
    }
    G4double m  = pd->GetPDGMass();                      // Mass of the particle
    G4ThreeVector mom;        // Direction of the momentum
    G4double ener;            // Kinetic energy
    if(knd<0)
    {
      mom = pdProj->GetMomentumDirection();        // Direction of the momentum
      ener = pdProj->GetKineticEnergy();                // Kinetic energy
    }
    else
    {
      mom = cHd->GetMomentumDirection();        // Direction of the momentum
      ener = cHd->GetKineticEnergy();                // Kinetic energy
    }
    G4double e = ener + m;                             // Total energy
	   G4double p = std::sqrt(ener*(e + m));                   // Abs value of the momentum
	   mom *= p;                                            // 3-momentum
    G4LorentzVector lorV =  G4LorentzVector(mom, e);
    G4int   c = pd->GetPDGEncoding();
    G4int ns=pd->GetQuarkContent(3)-pd->GetAntiQuarkContent(3);
    G4int nz=static_cast<G4int>(pd->GetPDGCharge());
    G4int nn=pd->GetBaryonNumber()-nz-ns;
    //if(!c)
    //{
    //  nz=static_cast<G4int>(pd->GetPDGCharge());
    //  nn=static_cast<G4int>(pd->GetBaryonNumber())-nz;
				//  //c=90000000+nz*1000+nn;     // Old CHIPS notation
				//  c=1000000000+nz*10010+nn*10; // New PDG 2006 definition
    //}
    G4double px = lorV.x();
    G4double py = lorV.y();
    G4double pz = lorV.z();
	   //G4double pt2 = px*px+py*py;
	   //G4double p2  = pt2+pz*pz;
	   //if(pz/sqrt(p2)>0.984808) 
	   // *** For Ntuple case ***
    //float xtupin[43];
    //G4int indnt=0;
    //xtupin[indnt++] = nEvnt;
    //xtupin[indnt++] = fNH;
    //xtupin[indnt++] = fNH;
    //xtupin[indnt++] = Mprot;
    //xtupin[indnt++] = Mneut;
    //xtupin[indnt++] = Mdeut;
    //xtupin[indnt++] = Mtrit;
    //xtupin[indnt++] = MHe3;
    //xtupin[indnt++] = MHe4;
    //xtupin[indnt++] = Mgam;
    //xtupin[indnt++] = Mpip;
    //xtupin[indnt++] = Mpim;
    //xtupin[indnt++] = Mpi0;
    //xtupin[indnt++] = MKp;
    //xtupin[indnt++] = MK0;
    //xtupin[indnt++] = MKm;
    //xtupin[indnt++] = MaK0;
    //xtupin[indnt++] = Meta;
    //xtupin[indnt++] = Metap;
    //xtupin[indnt++] = Mrhom;
    //xtupin[indnt++] = Mrhop;
    //xtupin[indnt++] = Mrho0;
    //xtupin[indnt++] = Momega;
    //xtupin[indnt++] = Mphi;
    //xtupin[indnt++] = MKS0;
    //xtupin[indnt++] = MKSC;
    //xtupin[indnt++] = MaKS0;
    //xtupin[indnt++] = MaKSC;
    //xtupin[indnt++] = Mf2;
    //xtupin[indnt++] = Ma2m;
    //xtupin[indnt++] = Ma2p;
    //xtupin[indnt++] = Ma20;
    //xtupin[indnt++] = Mf2p;
    //xtupin[indnt++] = m;
	   //tupleIncl.fill(xtupin);
    tuplIncl<<nEvnt<<" "<<fNH  <<" "<<fNH   <<" "<<Mprot<<" "<<Mneut<<" "
            <<Mdeut<<" "<<Mtrit<<" "<<MHe3  <<" "<<MHe4 <<" "<<Mgam <<" "
            <<Mpip <<" "<<Mpim <<" "<<Mpi0  <<" "<<MKp  <<" "<<MK0  <<" "
            <<MKm  <<" "<<MaK0 <<" "<<Meta  <<" "<<Metap<<" "<<Mrhom<<" "
            <<Mrhop<<" "<<Mrho0<<" "<<Momega<<" "<<Mphi <<" "<<MKS0 <<" "
            <<MKSC <<" "<<MaKS0<<" "<<MaKSC <<" "<<Mf2  <<" "<<Ma2m <<" "
            <<Ma2p <<" "<<Ma20 <<" "<<Mf2p  <<" "<<  d  <<" "<<  c  <<" "
            << ns  <<" "<< nz  <<" "<< nn   <<" "<<  m  <<" "<< px  <<" "
            << py  <<" "<< pz  <<" "<<  e   <<" "<< G4endl;
  }
  
  // Calculate multiplicities & fill ntuple for Q (Level 1) tracks
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
  G4int i=begi;
  while (i<tNH)
  {
    if(i<0) pd  = pdProj->GetDefinition();
    else
    {
      cHd = hadrons->GetSecondary(i)->GetDynamicParticle();
      pd  = cHd->GetDefinition();
    }
    G4int c = pd->GetPDGEncoding();
    if(!c)
    {
      G4int chrg=static_cast<G4int>(pd->GetPDGCharge());
      G4int bary=static_cast<G4int>(pd->GetBaryonNumber());
      //c=90000000+chrg*999+bary;     // Old CHIPS notation
      c=1000000000+chrg*10000+bary*10; // New PDG 2006 definition
    }
	   if     (c==2212 || c==1000010010)      Mprot++;
	   else if(c==2112 || c==1000000010)      Mneut++;
	   else if(c==1000010020)                 Mdeut++;
	   else if(c==1000010030)                 Mtrit++;
	   else if(c==1000020030)                 MHe3++;
	   else if(c==1000020040)                 MHe4++;
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
  	 i++;
  }
  // *** For Ntuple case ***    
  //    evdnt=0;
  //    xtupev[evdnt++] = nEvnt;
  //    xtupev[evdnt++] = fNH;
  //    xtupev[evdnt++] = fNH;
  //    xtupev[evdnt++] = Mprot;
  //    xtupev[evdnt++] = Mneut;
  //    xtupev[evdnt++] = Mdeut;
  //    xtupev[evdnt++] = Mtrit;
  //    xtupev[evdnt++] = MHe3;
  //    xtupev[evdnt++] = MHe4;
  //    xtupev[evdnt++] = Mgam;
  //    xtupev[evdnt++] = Mpip;
  //    xtupev[evdnt++] = Mpim;
  //    xtupev[evdnt++] = Mpi0;
  //    xtupev[evdnt++] = MKp;
  //    xtupev[evdnt++] = MK0;
  //    xtupev[evdnt++] = MKm;
  //    xtupev[evdnt++] = MaK0;
  //    xtupev[evdnt++] = Meta;
  //    xtupev[evdnt++] = Metap;
  //    xtupev[evdnt++] = Mrhom;
  //    xtupev[evdnt++] = Mrhop;
  //    xtupev[evdnt++] = Mrho0;
  //    xtupev[evdnt++] = Momega;
  //    xtupev[evdnt++] = Mphi;
  //    xtupev[evdnt++] = MKS0;
  //    xtupev[evdnt++] = MKSC;
  //    xtupev[evdnt++] = MaKS0;
  //    xtupev[evdnt++] = MaKSC;
  //    xtupev[evdnt++] = Mf2;
  //    xtupev[evdnt++] = Ma2m;
  //    xtupev[evdnt++] = Ma2p;
  //    xtupev[evdnt++] = Ma20;
  //    xtupev[evdnt++] = Mf2p;
  //    tupleEvtQ.fill(xtupev);
  tuplEvtQ<<nEvnt<<" "<<fNH  <<" "<<fNH <<" "<<Mprot<<" "<<Mneut<<" "
          <<Mdeut<<" "<<Mtrit<<" "<<MHe3  <<" "<<MHe4 <<" "<<Mgam <<" "
          <<Mpip <<" "<<Mpim <<" "<<Mpi0  <<" "<<MKp  <<" "<<MK0  <<" "
          <<MKm  <<" "<<MaK0 <<" "<<Meta  <<" "<<Metap<<" "<<Mrhom<<" "
          <<Mrhop<<" "<<Mrho0<<" "<<Momega<<" "<<Mphi <<" "<<MKS0 <<" "
          <<MKSC <<" "<<MaKS0<<" "<<MaKSC <<" "<<Mf2  <<" "<<Ma2m <<" "
	         <<Ma2p <<" "<<Ma20 <<" "<<Mf2p  <<" "<< G4endl;
#ifdef pdebug
  G4cout<<"G4QHBook::FillEvt - End"<<G4endl;
#endif
}

G4QHBook::~G4QHBook() 
{
#ifdef pdebug
  G4cout<<"~G4QHBook(): Writing number of events histNevt = "<<nEvnt<<G4endl;
#endif
  histNevt<<nEvnt<<" "<< G4endl;
}
