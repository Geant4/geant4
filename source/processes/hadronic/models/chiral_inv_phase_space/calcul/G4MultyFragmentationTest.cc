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
//#define debug
//#define pdebug

//#include "G4UIterminal.hh"
#include "G4ios.hh"
#include "G4QEnvironment.hh"

int main()
{
  //G4double theEnergyLossPerFermi = 1.*GeV;
  G4int nop = 152;                               // clusters (A<6)
  G4double fractionOfSingleQuasiFreeNucleons = 0.5; // It is A-dependent (C=.85, U=.40)
  G4double fractionOfPairedQuasiFreeNucleons = 0.05;
  G4double clusteringCoefficient = 5.;
  G4double temperature = 180.;
  G4double halfTheStrangenessOfSee = 0.3; // = s/d = s/u
  G4double etaToEtaPrime = 0.3;
  G4double fusionToExchange = 100.;
  //G4double theInnerCoreDensityCut = 50.;
#ifdef debug
	 G4cout<<"G4MulFragTst: nop="<<nop<<",fN="<<fractionOfSingleQuasiFreeNucleons
        <<",fD="<<fractionOfPairedQuasiFreeNucleons<<",wC="<<clusteringCoefficient
        <<",vn="<<fusionToExchange<<",T="<<temperature
        <<",ss="<<halfTheStrangenessOfSee<<",se="<<etaToEtaPrime<<G4endl;
#endif
  //-------- Initialize CHIPS
  G4QCHIPSWorld* theW=G4QCHIPSWorld::Get();
  theW->GetParticles(nop);           // Create CHIPS World of nop particles
  G4QNucleus::SetParameters(fractionOfSingleQuasiFreeNucleons,
                            fractionOfPairedQuasiFreeNucleons,
			                         clusteringCoefficient, 
					                       fusionToExchange);
  G4Quasmon::SetParameters(temperature, halfTheStrangenessOfSee, etaToEtaPrime);

  G4int tPDG = 90006005;    // PDG of the Target Nucleus
  G4int nEvents= 100000;
  //#ifdef debug
	 G4cout<<"G4MulFragTst: targetPDG="<<tPDG<<", nEv="<<nEvents<<G4endl;
  //#endif
  G4QHadronVector projHV;

  G4int PDG1=211;
  G4LorentzVector LV1(-485.221,-216.252,267.887,611.103);
  //#ifdef debug
	 G4cout<<"G4MulFragTst: #1 projPDG="<<PDG1<<", 4M="<<LV1<<G4endl;
  //#endif

  G4int PDG2=-211;
  G4LorentzVector LV2(-201.944,475.785,952.224,1092.41);
  //#ifdef debug
	 G4cout<<"G4MulFragTst: #2 projPDG="<<PDG2<<", 4M="<<LV2<<G4endl;
  //#endif

  G4int PDG3=111;
  G4LorentzVector LV3(213.293,183.683,996.566,1044.32);
  //#ifdef debug
	 G4cout<<"G4MulFragTst: #3 projPDG="<<PDG3<<", 4M="<<LV3<<G4endl;
  //#endif

  //std::ifstream inFile("chipstest.in", std::ios::in);
  //G4double temperature;
  //G4double ssin2g;
  //G4double eteps;
  //G4double momb;
  //G4double enb;
  //G4double cP;
  //G4double fN;
  //G4double fD;
  //G4double rM;
  //G4double sA;
  //G4int    nop;
  //G4int    pPDG;
  //G4int    tPDG;
  //G4int    nEvt;
  ////G4int    nofdecays;
  ////G4int    decmask=0;
  //inFile>>temperature>>ssin2g>>eteps>>nop>>momb>>enb>>pPDG>>tPDG>>nEvt>>fN>>fD>>cP>>rM>>sA;
  //G4cout<<"CHIPStest Par's: Temp="<<temperature<<",SSs="<<ssin2g<<",P/S="<<eteps<<",nop="
  //      <<nop<<",p="<<momb<<",e="<<enb<<",projCode="<<pPDG<<",targCode="<<tPDG<<",nEv="
  //      <<nEvt<<",fN="<<fN<<",fD="<<fD<<",cP="<<cP<<",rM="<<rM<<",sA="<<sA<<G4endl;
  //-------- Initialize CHIPS
  //G4QCHIPSWorld* theW=G4QCHIPSWorld::Get();
  //theW->GetParticles(nop);           // Create CHIPS World of nop particles
  //G4Exception("***CHIPStest: TMP");
  //G4QNucleus::SetParameters(fN,fD,cP,rM);
  //G4Quasmon::SetParameters(temperature,ssin2g,eteps);
  //G4QEnvironment::SetParameters(sA);                 // SolAngle (pbar-A secondary capture)
  //-------- write results onto a file --------
  //std::ofstream outFile( "chipstest.out", std::ios::out);
  //outFile.setf( std::ios::scientific, std::ios::floatfield );

  //--------- Example how to write in the outFile -----------
  //outFile << " " << endl;

  // *********** Now momb is a momentum of the incident particle *************
  //G4double mp=G4QPDGCode(pPDG).GetMass();             // @@ just for the check
  //G4QContent pQC=G4QPDGCode(pPDG).GetQuarkContent();
  //G4int    cp=pQC.GetCharge();
  //G4int    bnp=pQC.GetBaryonNumber();
  //G4double mp1=G4QPDGCode(PDG1).GetMass();             // @@ just for the check
  G4QContent pQC1=G4QPDGCode(PDG1).GetQuarkContent();
  G4int    cp1=pQC1.GetCharge();
  G4int    bnp1=pQC1.GetBaryonNumber();
  //G4double mp2=G4QPDGCode(PDG2).GetMass();             // @@ just for the check
  G4QContent pQC2=G4QPDGCode(PDG2).GetQuarkContent();
  G4int    cp2=pQC2.GetCharge();
  G4int    bnp2=pQC2.GetBaryonNumber();
  //G4double mp3=G4QPDGCode(PDG3).GetMass();             // @@ just for the check
  G4QContent pQC3=G4QPDGCode(PDG3).GetQuarkContent();
  G4int    cp3=pQC3.GetCharge();
  G4int    bnp3=pQC3.GetBaryonNumber();
  //momb=momb;
  //G4double ep=sqrt(mp*mp+momb*momb);                  // @@ just for the check
  //if(enb>0.) ep=enb;
  G4double mt=G4QPDGCode(tPDG).GetMass();             // @@ just for the check
  G4QContent tQC=G4QPDGCode(tPDG).GetQuarkContent();
  G4int    ct=tQC.GetCharge();
  G4int    bnt=tQC.GetBaryonNumber();
  G4int    totC=cp1+cp2+cp3+ct;
  G4int    totBN=bnp1+bnp2+bnp3+bnt;
  G4LorentzVector preSumLV(0.,0.,0.,mt);
  preSumLV+=LV1+LV2+LV3;
  G4double fEvt=nEvents;
  G4double sumE=0.;
  G4double sumK=0.;
  G4double sumG=0.;
  G4double sumT=0.;
  G4double sumN=0.;
  G4double sum0=0.;
  G4double sumP=0.;
  G4double sum1N=0.;
  G4double sumNN=0.;
  G4double sumPP=0.;
  G4double sumAL=0.;
  G4double time=clock()/CLOCKS_PER_SEC;
  // Main LOOP over events ======================================
  for (G4int ir=0; ir<nEvents; ir++)
  {
    // Randomization loop: cycle random generator, using 2 lower digits in nEvents
	   G4int    iRandCount = nEvents%100;
	   G4double vRandCount = 0.;
	   while (iRandCount>0)
	   {
	     vRandCount = G4UniformRand();
 	    iRandCount--;
	   }
    if(!(ir%1000) && ir) G4cout<<"G4MultFragTst: "<<ir<<" events are simulated"<<G4endl;
    //G4cout<<"G4MultFragTst: "<<ir<<" events are simulated"<<G4endl;
    G4LorentzVector totSum = preSumLV;
    G4QHadron* H1 = new G4QHadron(PDG1,LV1);
    G4QHadron* H2 = new G4QHadron(PDG2,LV2);
    G4QHadron* H3 = new G4QHadron(PDG3,LV3);
    G4int           totCharge = totC;
    G4int           totBaryN = totBN;
    G4QHadronVector projHV;
    projHV.push_back(H1);                                 // DESTROYED over 3 line
    projHV.push_back(H2);                                 // DESTROYED over 2 line
    projHV.push_back(H3);                                 // DESTROYED over 1 line
    G4QEnvironment* pan= new G4QEnvironment(projHV,tPDG); // DELETED over 8 lines
    std::for_each(projHV.begin(), projHV.end(), DeleteQHadron());
    projHV.clear();
#ifdef debug
    G4cout<<"CHIPStest:===>>> Now call Fragment (HadronizeQuasmon) function" << G4endl;
#endif
    G4QHadronVector* output; // Prototype of the output
    try
    {
      output = pan->Fragment();// DESTROYED in the end of the LOOP work space
    }
    catch (G4QException& error)
    {
#ifdef pdebug
      G4cout<<"***CHIPStest: Exception is catched"<<G4endl;
#endif
      G4cerr<<"***CHIPStest Abort: "<<error.GetMessage()<<G4endl;
      abort();
    }
#ifdef debug
    G4cout<<"CHIPStest:--->>>Now come out of Fragment (HadronizeQuasmon) function"<<G4endl;
#endif
    delete pan; // Destruct theHadronVector (& theCandidateVector) of the Quasmon
#ifdef debug
    G4cout << "CHIPStest: >>> Here the histograms are filled" << G4endl;
#endif
    G4int tNH = output->size();
    G4int npt=0;
    G4int nGamma=0;
    G4double EGamma=0;
    G4int nP0=0;
    G4int nPP=0;
    G4int nPN=0;
    G4int nKaons=0;
    G4int nEta=0;
    G4int nAlphas=0;
    G4int nPhotons=0;
    G4int nProtons=0;
    G4int nNeutrons=0;
    G4int nSpNeut=0;
    G4int nSpAlph=0;
    G4int nOmega=0;
    G4int nDec=0;
    G4int dirN=0;
#ifdef pdebug
    G4cout<<"----------DONE^^^^^^^************^^^^^^^^^^^:ir="<<ir<<": #ofH="<<tNH<<G4endl;
    if(!(ir%100)) G4cerr<<"#"<<ir<<G4endl;
#endif
    G4bool alarm=false;
    G4bool rad=false;
    G4bool hyp=false;
    G4bool badPDG=false;
    for (G4int ind=0; ind<tNH; ind++)
    {
      G4QHadron* curH=output->operator[](ind);
      G4double      m = curH->GetMass();            // Mass of the particle
      G4LorentzVector lorV = curH->Get4Momentum();  // 4-momentum of the particle
      if(std::fabs(m-lorV.m())>.005)
	     {
		      G4cerr<<"***CHIPStest: m="<<lorV.m()<<" # "<<m<<", d="<<lorV.m()-m<<G4endl;
        alarm=true;
	     }
      if(!(lorV.e()>=0||lorV.e()<0)   || !(lorV.px()>=0||lorV.px()<0) ||
         !(lorV.py()>=0||lorV.py()<0) || !(lorV.pz()>=0||lorV.pz()<0))
	     {
		      G4cerr<<"***CHIPStest: NAN in LorentzVector="<<lorV<<G4endl;
        alarm=true;
	     }
      G4int         d=curH->GetNFragments(); // In how many particles this particle decayed
      G4ThreeVector p = lorV.vect();         // 3-momentum of the particle
      G4double      e = lorV.e();            // Energy of the particle
      G4int         c=curH->GetPDGCode();    // PDG Code of the particle
      //if(!d&&(c==90000002||c==90002000||c==92000000||c==221||c==331))
      if(!d&&(c==90000002||c==90002000||c==92000000))
      {
        //G4cout<<"***CHIPStest:***Dibaryon or Eta*** ind="<<ind<<", PDG="<<c<<G4endl;
        G4cout<<"***CHIPStest:***Dibaryon *** ind="<<ind<<", PDG="<<c<<G4endl;
        alarm=true;
      }
      if(!d&&(c==90000003||c==90003000||c==93000000))
      {
        G4cout<<"***CHIPStest:***Tribaryon *** ind="<<ind<<", PDG="<<c<<G4endl;
        alarm=true;
      }
      if(!d) npt++;
      if(d) nDec+=d;
      if(c==223) nOmega++;
      if(c==22) nPhotons++;
      if(c==311||c==321||c==-311||c==-321) nKaons++; // kaons
      if(c==221) nEta++;                             // etas
      if(c==90002002) nAlphas++;                     // Alphas
      if(c==2212 || c==90001000) nProtons++;         // Protons
      if(c==90000001 || c==90001000) dirN++;         // Ditrect nucleons
      if(c==2112 || c==90000001) nNeutrons++;        // Neutrons
      if((c==2112 || c==90000001) && std::fabs(e-1005.)<3.) nSpNeut++;// Dibar-Neutrons
      if(!d && c==90002002 && e-m<7.) nSpAlph++;     // Special Alphas
      if(c==111) nP0++;                              // Neutral  pions
      if(c==-211) nPN++;                             // Negative pions
      if(c==211) nPP++;                              // Positive pions
      if(c==22) nGamma++;                            // Gammas
      if(c==22) EGamma+=e;                           // Energy of gammas
      if(!d) totCharge-=curH->GetCharge();
      if(!d) totBaryN-=curH->GetBaryonNumber();
      if(!d) totSum    -= lorV;
      if(c>80000000 && (c<90000000 || c%1000>500 || c%1000000>500000) ||
         !(c>=0 || c<0))
      {
        G4cout<<"***OUTPUT ERROR*** CHIPStest: bad PDG is found. It is "<<c<<G4endl;
        badPDG=true;
	     }
      rad=rad||c==90002000||c==90003000||c==90004000||c==90000002||c==90000003||c==90000004
             ||c==90002003||c==90003002||c==90004002||c==90002005||c==90005002||c==90004004
             ||c==90006002;
      hyp = hyp || c>90999999;
#ifdef pdebug
      G4cout<<"#"<<ind<<"(d="<<d<<"), PDG="<<c<<",4M="<<lorV<<m<<",T="<<lorV.e()-m<<G4endl;
#endif
    }
#ifdef pdebug
    G4cout<<"CHECK: 4M="<<totSum<<", Charge="<<totCharge<<", BaryN="<<totBaryN<<G4endl;
#endif
    G4double ss=std::fabs(totSum.t())+std::fabs(totSum.x())+std::fabs(totSum.y())+std::fabs(totSum.z());    
	   if (totCharge ||totBaryN || !(ss<.01) || alarm || nGamma&&!EGamma || badPDG)
	   //if (totCharge || ss>.01 || alarm || nSpNeut)
	   //if (totCharge || ss>.01 || alarm || nSpAlph)
    {
      G4cerr<<"***CHIPStest:#"<<ir<<":n="<<tNH<<",4M="<<totSum<<",Charge="<<totCharge
            <<",BaryN="<<totBaryN<<G4endl;
      if(nGamma&&!EGamma)G4cerr<<"***CHIPStest: Egamma=0"<<G4endl;
      totSum    = preSumLV;
      for (int indx=0; indx<tNH; indx++)
      {
        G4QHadron* curH=output->operator[](indx);
        G4double m = curH->GetMass();
        G4LorentzVector lorV = curH->Get4Momentum();
        G4int d=curH->GetNFragments();
        G4int c=curH->GetPDGCode();
        if(!d) totSum    -= lorV;
        G4cerr<<"#"<<indx<<"("<<d<<"), PDG="<<c<<", m/LV="<<m<<lorV<<", T="<<lorV.e()-m
              <<", d4M="<<totSum<<G4endl;
      }
      G4Exception("***CHIPStest: ALARM or charge/energy/momentum is not conserved");
    }
    sum0+=nP0;
    sumP+=nPP;
    sumN+=nPN;
    G4int nPions=nP0+nPP+nPN;
    sumG+=nGamma;
    sumT+=EGamma;
    if(nAlphas)sumAL++;
    if(nProtons)sumPP++;
    if(nNeutrons)sumNN++;
    if(nNeutrons&&!nProtons&&!nAlphas) sum1N++;
    if(nKaons)sumK++;
    if(nEta)sumE++;
    if (nPhotons) nPions=11;
    if (nKaons==2&&npt==2) nPions=1;
    else if (nKaons) nPions=10;
    //histPi.fill(nPions);
    //histNeut.fill(nNeutrons);
    std::for_each(output->begin(), output->end(), DeleteQHadron());
    output->clear();
    delete output;
  }
  time=(clock()/CLOCKS_PER_SEC-time)/fEvt;
  G4cerr<<"CHIPStest::t="<<time<<",Yields:pi-="<<sumN/fEvt<<",pi+="<<sumP/fEvt<<",pi0="
        <<sum0/fEvt<<",K="<<sumK/fEvt<<",eta="<<sumE/fEvt<<",gamma="<<sumG/fEvt<<"(<E>="
        <<sumT/sumG<<"),n="<<sumNN/fEvt<<",p="<<sumPP/fEvt<<",alpha="<<sumAL/fEvt
        <<",onlyN="<<sum1N/fEvt<<G4endl;
  return EXIT_SUCCESS;
}
