///////////////////////////////////////////////////////////////////////////////
// File: CHIPStest.cc
// Author: M.Kossov
// Last modification: 08/99
// Description: Main function for the Test of CHIPS generator for Geant4
///////////////////////////////////////////////////////////////////////////////

//#define debug
//#define pdebug
//#define ppdebug

#include "G4UIterminal.hh"
#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4DynamicParticle.hh"

//#ifdef XX
//  #include "YY.hh"
//#endif

#include "CHBook.h"
#include "G4QNucleus.hh"
#include "G4Quasmon.hh"
#include "G4QHBook.hh"
#include "CHBook.cc"
#include "G4QHBook.cc"

#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>

//int main(int argc, char** argv) 
int main()
{
  //if (argc == 1) 
  //// Interactive mode
  //{
  //  G4UIsession* session = new G4UIterminal;
  //  session->SessionStart();
  //  delete session;
  //}
  ////Batch mode
  //else {
  //  G4String command = "/control/execute ";
  //  G4String fileName = argv[1];
  //  UI->ApplyCommand(command+fileName);    
  //}

  // Histograms
  cout<<"Prepare chips.hbook file for HBOOK histograms"<<endl;
  CHBookFile chipshist("chips.hbook",1,1);

  // Book histograms
  CHBookHisto histPi(3027,"Pion multiplicity",12,0.,12.);
  //
  G4QHBook ntp;
  //
  //-------- set standard output format-------
   G4cout.setf( ios::scientific, ios::floatfield );

  //-------- take parameters from a file -------
  ifstream inFile("chipstest.in", ios::in);
  G4double temperature;
  G4double ssin2g;
  G4double etaetap;
  G4double momb;
  G4double enb;
  G4double fN;
  G4double fD;
  G4double cP;
  G4int    nop;
  G4int    pPDG;
  G4int    tPDG;
  G4int    nEvt;
  G4int    nofdecays;
  G4int    decmask=0;
  inFile>>temperature>>ssin2g>>etaetap>>nop>>momb>>enb>>pPDG>>tPDG>>nEvt>>fN>>fD>>cP;
  cout<<"CHIPS Par's: Temp="<<temperature<<",ss_sea="<<ssin2g<<",P/S="<<etaetap<<",nop="<<nop
      <<",p="<<momb<<",e="<<enb<<",projCode="<<pPDG<<",targCode="<<tPDG<<",nEv="<<nEvt
      <<",fN="<<fN<<",fD="<<fD<<",cP="<<cP<<endl;
  G4QNucleus::SetParameters(fN,fD,cP);
  G4Quasmon::SetParameters(temperature,ssin2g,etaetap);
  //-------- write results onto a file --------
  ofstream outFile( "chipstest.out", ios::out);
  outFile.setf( ios::scientific, ios::floatfield );

  //--------- Example how to write in the outFile -----------
  //outFile << " " << endl;

  // *********** Now momb is a momentum of the incident particle *************
  //G4double p =0.;
  G4double mp=G4QPDGCode(pPDG).GetMass();             // @@ just for the check
  G4double ep=sqrt(mp*mp+momb*momb);                  // @@ just for the check
  if(enb>0.) ep=enb;
  G4double mt=G4QPDGCode(tPDG).GetMass();             // @@ just for the check
  G4QContent pQC=G4QPDGCode(pPDG).GetQuarkContent();
  G4int    cp=pQC.GetCharge();
  G4QContent tQC=G4QPDGCode(tPDG).GetQuarkContent();
  G4int    ct=tQC.GetCharge();
#ifdef debug
  cout<<"Main: pQC"<<pQC<<", pch="<<cp<<", tQC"<<tQC<<", tch="<<ct<<endl;
#endif
  G4int    totC=cp+ct;
  G4LorentzVector proj4Mom(0.,0.,momb,ep);
  G4LorentzVector targ4Mom(0.,0.,0.,mt);  
  G4double sumE=0.;
  G4double sumK=0.;
  G4double sumN=0.;
  G4double sum0=0.;
  G4double sumC=0.;
  G4double fEvt=nEvt;
  // Main LOOP over events ======================================
  for (G4int ir=0; ir<nEvt; ir++)
  {
    // Randomization loop: cycle random generator, using 2 lower digits in nEvt
	G4int    iRandCount = nEvt%100;
    G4double vRandCount = 0.0;
	while (iRandCount>0)
	{
	  vRandCount = G4UniformRand();
	  iRandCount--;
	}
    if(!(ir%10000) && ir) cout<<"CHIPS: "<<ir<<" events are simulated"<<endl;
    G4LorentzVector totSum    = G4LorentzVector(0.,0.,momb,ep+mt);
    G4int           totCharge = totC;
    // Initialization of a quasmon (Temporary wide)
    G4Quasmon* pan= new G4Quasmon(pPDG,tPDG,proj4Mom,targ4Mom,nop);
#ifdef debug
    cout << "===>>> Now call HadronizeQuasmon hadronization function" << endl;
#endif
    G4QHadronVector output = pan->HadronizeQuasmon();
#ifdef debug
    cout << "--->>> Now came out of HadronizeQuasmon hadronization function" << endl;
#endif
	//
    ntp.FillEvt(output);
	//
#ifdef debug
    cout << ">>> Here the histograms are filled" << endl;
#endif
    G4int tNH = output.entries();
    G4int npt=0;
    G4int nPions=0;
    G4int nP0=0;
    G4int nPC=0;
    G4int nKaons=0;
    G4int nEta=0;
    G4int nPhotons=0;
    G4int nOmega=0;
    G4int nDec=0;
#ifdef ppdebug
    cout<<"DONE^^^^^^^^^^^^:ir="<<ir<<": A#ofHadrons="<<tNH<<",p="<<pan<<",o="<<output[0]<<endl;
#endif
    for (G4int ind=0; ind<tNH; ind++)
    {
      G4double m = output[ind]->GetMass();
      G4LorentzVector lorV = output[ind]->Get4Momentum();
      G4int d=output[ind]->GetNFragments();
      G4ThreeVector p = lorV.vect();
      G4int c=output[ind]->GetPDGCode();
      if (!d) npt++;
      if (d) nDec+=output[ind]->GetNFragments();
      if(c==223) nOmega++;
      if(c==22) nPhotons++;
      if(c==311||c==321||c==-311||c==-321) nKaons++; // kaons
      if(c==221) nEta++;                             // etas
      if(c==211 || c==-211) nPC++;                   // Only charged
      if(c==111) nP0++;                              // Only neutral
      if(c==111 || c==211 || c==-211) nPions++;      // All pions
      if(!d) totCharge-=output[ind]->GetCharge();
      if(!d) totSum    -= lorV;
#ifdef pdebug
      cout<<"#"<<ind<<"("<<d<<"), PDG="<<c<<", m="<<m<<", LV="<<lorV<<", Tk="<<lorV.e()-m<<endl;
#endif
	  // Write each particle's type, mass, and 4-momentum to DST
	  //outFile<<ind<<" "<<c<<" "<<m<<" "<<lorV.x()<<" "<<lorV.y()
      //       <<" "<<lorV.z()<<" "<<lorV.t()<<" "<< endl;
    }
#ifdef pdebug
    cout << "CHECK: Lorentz Vector sum:" << totSum << ", Charge="<<totCharge<<endl;
#endif
    G4double ss=abs(totSum.t())+abs(totSum.x())+abs(totSum.y())+abs(totSum.z());    
	if (totCharge || ss>0.005)
    {
      cout<<"***Exception:n="<<ir<<": #ofH="<<tNH<<",res4M="<<totSum<<",c="<<totCharge<<endl;
      for (int indx=0; indx<tNH; indx++)
      {
        G4double mx = output[indx]->GetMass();
        G4LorentzVector lorVx = output[indx]->Get4Momentum();
        G4int dx=output[indx]->GetNFragments();
        G4int cx=output[indx]->GetPDGCode();
        cout<<"#"<<indx<<"("<<dx<<"), PDG="<<cx<<",m="<<mx<<",LV="<<lorVx<<endl;
      }
      G4Exception("********** Charge or energy/momentum is not conserved");
    }
	// WRITE in FILE  ---------------
    sum0+=nP0;
    sumC+=nPC;
    sumN+=nPions;
    if(nKaons)sumK++;
    if(nEta)sumE++;
    if (nPhotons) nPions=11;
    if (nKaons==2&&npt==2) nPions=1;
    else if (nKaons) nPions=10;
    histPi.fill(nPions);
    delete pan; // Destruct theHadronVector (& theCandidateVector) of the Quasmon
  }
  cout<<"The mean number of pions="<<sumN/fEvt<<", pi+-="<<sumC/fEvt<<", pi0="<<sum0/fEvt
      <<", K="<<sumK/fEvt<<", eta="<<sumE/fEvt<<endl;
  // Save all collected histograms
  cout<<"Save all collected histograms to the chips.hbook file"<<endl;
  chipshist.saveAndClose();

  return EXIT_SUCCESS;
}







