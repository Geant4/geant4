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
//       1         2         3         4         5         6         7         
//34567890123456789012345678901234567890123456789012345678901234567890123456789
//
//
// $Id$
//
//      ---------------- G4QEnvironment ----------------
//             by Mikhail Kossov, August 2000.
//  class for Multy Quasmon Environment used by the CHIPS Model
// ------------------------------------------------------------
// Short description: The G4QEnvironment class corresponds to the nuclear
// environment,  containing excited nucleons or nuclear clusters
// (G4Quasmons). In the constructer the nucleus (G4QNucleus) is clusterized
// and then the projectile hadron (or hadrons) can create one (or a few)
// Quasmons, which are fragmented by the Fragment member function. As a
// result a vector of G4QHadrons is created, which can include the residual
// nucleus in the ground state.
//---------------------------------------------------------------------
 
//#define debug
//#define pdebug
//#define qdebug
//#define chdebug
//#define sdebug
//#define ppdebug
//#define cdebug
//#define cldebug
//#define edebug
//#define rdebug
//#define fdebug
//#define ffdebug
//#define pcdebug
//#define mudebug

#include <cmath>
#include <cstdlib>

#include "G4QEnvironment.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

G4QEnvironment::G4QEnvironment(const G4QNucleus& theEnv)
 : theEnvironment(theEnv)
{
  theQFScat = G4QFreeScattering::GetPointer();
  G4int envPDG = theEnv.GetPDG();
  G4QPDGCode envQPDG(envPDG);
  G4int envA = envQPDG.GetBaryNum();
  G4double envM = envQPDG.GetMass();
  theWorld = 0;
  nBarClust = 0;
  f2all = 0;
  totCharge = envQPDG.GetCharge();
  totBaryoN = envA;
  tot4Mom = G4LorentzVector(0.,0.,0.,envM);
  theTargetPDG = 0;
#ifdef debug
  G4cout << "G4QEnviron::Const: t4M=" << tot4Mom << ",tC=" << totCharge
         << ",tB=" << totBaryoN << G4endl;
#endif
}


G4QEnvironment::G4QEnvironment(const G4QHadronVector& projHadrons,
                               const G4int targPDG)
 : theEnvironment(90000000)     // User is responsible for projHadrons(Vector)
{
  //static const G4double mNeut= G4QPDGCode(2112).GetMass();
  //static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QPDGCode pimQPDG(-211);
  theQFScat = G4QFreeScattering::GetPointer();
  theWorld = G4QCHIPSWorld::Get();        // Get a pointer to the CHIPS World
  nBarClust = 0;
  totCharge = 0;
  totBaryoN = 0;
  f2all = 0;
  G4bool fake=false;                      // At present only fake pi-
  theTargetPDG=targPDG;                   // Remenber it for error message
  G4int nHadrons=projHadrons.size();      // A#of hadrons in the input Vector

#ifdef debug
  G4cout<<"--->>G4QE::Const: Called targPDG="<<targPDG<<", nInpHadr="<<nHadrons<<G4endl;
#endif
  if(nHadrons<1 || targPDG==90000000)    // No projectile Hadrons or no target Nucleus
  {
    G4cout << "---Warning---G4QEnv::Const:a#ofINPHadr=" << nHadrons
           << ",tPDG=" << targPDG << G4endl;
    //throw G4QException("***G4QEnvironment: There is no one projectile or vacuum target");
    if(nHadrons)              // The projectiles are copied to the output
    {
      for (G4int ih=0; ih<nHadrons; ih++)
      {
        G4QHadron* curQH    = new G4QHadron(projHadrons[ih]);
#ifdef debug
        G4cout << "*G4QE::Const:iH#" << ih << "," << curQH->GetQC()
               << curQH->Get4Momentum() << G4endl;
#endif

        if (curQH->GetPDGCode() == 10) {
          // Chipolino is found in the input -> Decay

          G4QContent   chQC=curQH->GetQC();
          // Quark content of the Hadron-Chipolino
          G4QChipolino QCh(chQC);
          // Define a Chipolino instance for the Hadron

          G4LorentzVector ch4M=curQH->Get4Momentum(); // 4Mom of the Hadron-Chipolino
          G4QPDGCode h1QPDG=QCh.GetQPDG1();  // QPDG of the first hadron
          G4double   h1M   =h1QPDG.GetMass();// Mass of the first hadron
          G4QPDGCode h2QPDG=QCh.GetQPDG2();  // QPDG of the second hadron
          G4double   h2M   =h2QPDG.GetMass();// Mass of the second hadron
          G4double   chM2  =ch4M.m2();       // Squared Mass of the Chipolino
          if( sqr(h1M+h2M) < chM2 )          // Decay is possible
          {
            G4LorentzVector h14M(0.,0.,0.,h1M);
            G4LorentzVector h24M(0.,0.,0.,h2M);
            if(!G4QHadron(ch4M).DecayIn2(h14M,h24M))
            {
              G4ExceptionDescription ed;
              ed << "QChip DecIn2 error: CM=" << std::sqrt(chM2) << " -> h1="
                 << h1QPDG << "(" << h1M << ") + h2=" << h1QPDG << "(" << h2M
                 << ") = " << h1M+h2M << " **Failed**" << G4endl;
              G4Exception("G4QEnvironment::G4QEnvironment()", "HAD_CHPS_0000",
                          FatalException, ed);
            }
            delete curQH;                    // Kill the primary Chipolino
            G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);
            theQHadrons.push_back(h1H);      // (delete equivalent)
            curQH    = new G4QHadron(h1H);   // ... just to remember independently
            theProjectiles.push_back(curQH); // Remember it for the error message

#ifdef debug
            G4cout<<"G4QE::Constr: QChipolino -> H1="<<h1QPDG<<h14M<<G4endl;
#endif
            curQH = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
            theQHadrons.push_back(curQH);    // (delete equivalent)
#ifdef debug
            G4cout<<"G4QE::Constr: QChipolino -> H2="<<h2QPDG<<h24M<<G4endl;
#endif
          }
          else
          {
            G4ExceptionDescription ed;
            ed << "LowMassChipolino in Input: " << ih << "," << curQH->GetQC()
               << curQH->Get4Momentum() << ", chipoM=" << std::sqrt(chM2) << " < m1="
               << h1M << "(" << h1QPDG << ") + m2=" << h2M << "(" << h2QPDG << ") = "
               << h1M+h2M << G4endl;
            G4Exception("G4QEnvironment::G4QEnvironment()", "HAD_CHPS_0001",
                        FatalException, ed);
          }
        }
        theQHadrons.push_back(curQH);        // (delete equivalent)
        curQH    = new G4QHadron(curQH);     // ... just to remember independently
        theProjectiles.push_back(curQH);     // Remenber it for the error message
      }
    }
    else if(targPDG!=90000000)               // No projHadrons,fill targetNucleus to output
    {
      G4QHadron* curQH    = new G4QHadron(targPDG);
#ifdef debug
      G4cout<<"**G4QE::Const:No iHad,eH="<<curQH->GetQC()<<curQH->Get4Momentum()<<G4endl;
#endif
      theQHadrons.push_back(curQH);          // (delete equivalent)
    }
    if (nHadrons<0) G4cout<<"***Warning****G4QE::Const:NH="<<nHadrons<<" < 0 !"<<G4endl;
    return;
  }
  G4QPDGCode targQPDG(targPDG);
#ifdef debug
  G4cout<<"G4QE::C:targQPDG="<<targQPDG<<G4endl;
#endif
  G4int    targA=targQPDG.GetBaryNum();
  G4double targM=targQPDG.GetMass();
  totCharge=targQPDG.GetCharge();
  totBaryoN=targA;
  tot4Mom=G4LorentzVector(0.,0.,0.,targM);
  // === Print out of the input information at Creation time & tot 4-mom Calculation ===
#ifdef debug
  G4cout<<"G4QE::C:PDG="<<targPDG<<",C="<<totCharge<<",M="<<targM<<",n="<<nHadrons<<G4endl;
#endif
  for(G4int ipr=0; ipr<nHadrons; ipr++)// LOOP is used for the tot4Mom calc. & for printing
  {
    G4QHadron* prHadr = projHadrons[ipr];
    G4QHadron* curQH  = new G4QHadron(prHadr);// Remenber it for _
    theProjectiles.push_back(curQH);          // the error message
    G4LorentzVector h4Mom = prHadr->Get4Momentum();
    tot4Mom      += h4Mom;
    totCharge    += prHadr->GetCharge();
    totBaryoN    += prHadr->GetBaryonNumber();
#ifdef debug
    G4int           hPDG  = prHadr->GetPDGCode();
    G4int           hNFrag= prHadr->GetNFragments();
    G4QContent      hQC   = prHadr->GetQC();
    G4cout<<"G4QE::C:#"<<ipr<<",PDG="<<hPDG<<hQC<<",4M="<<h4Mom<<",hNFr="<<hNFrag<<G4endl;
#endif
  }
#ifdef debug
  G4cout<<"G4QEnv::Const:tC="<<totCharge<<",tB="<<totBaryoN<<",tot4M="<<tot4Mom<<G4endl;
#endif
#ifdef debug
  G4cout<<"G4QEnv::Const: --> tC="<<totCharge<<",tB="<<totBaryoN<<G4endl;
#endif
  G4int nP=theWorld->GetQPEntries();         // A#of init'ed particles in CHIPS World
  //G4int nCl=nP-90;                           // A#of init'ed clusters in CHIPS World
  G4int nCl=nP-53;           // @@ A#ofClusters in CHIPSWorld (53=nQHM in G4QPDGCode.hh)
#ifdef debug
  G4cout<<"G4QEnv:Const:Before NCI:n="<<nP<<",F="<<projHadrons[0]->GetNFragments()<<",tC="
        <<totCharge<<",tB="<<totBaryoN<<", nCl="<<nCl<<G4endl;
#endif
  InitClustersVector(nCl,targA);             // Init Clusters as Particles (to interact)
#ifdef debug
  G4cout<<"G4QEnv::Const:NucClust,n="<<nCl<<",F="<<projHadrons[0]->GetNFragments()<<",tC="
        <<totCharge<<",tB="<<totBaryoN<<G4endl;
#endif
  if(targPDG>80000000)                        // ==> Nuclear target (including NUCPDG)
  {
    theEnvironment.InitByPDG(targPDG);        // Create nuclear environment
#ifdef debug
    G4cout<<"G4QEnv::Const:nH="<<nHadrons<<",PDG="<<projHadrons[0]->GetPDGCode()<<",tC="
          <<totCharge<<",tB="<<totBaryoN<<G4endl;
#endif
    if(nHadrons==1)
    {
      G4QHadron* opHad=projHadrons[0];
      G4int opPDG=opHad->GetPDGCode();
#ifdef debug
      G4cout<<"G4QEnviron::Constructor: *** Only one input hadron*** PDG="<<opPDG<<G4endl;
#endif
      if(opPDG==22)                           // *** Check photon's NuclearSplitThreshold
      {
        G4double exMass=tot4Mom.m();
#ifdef debug
        G4cout<<"G4QEnvironment::Const: exM="<<exMass-targM<<" > mPi0 ?"<<G4endl;
#endif      
        if(exMass<targM+135.977) // Nucleus is below the pion production threshold
        {
          G4QNucleus exEnviron(tot4Mom,targPDG);
          // @@ One can put here the pbpt= (M.K.) @@ What about d,t,alpha splitting?
          if(targM>999.&&!exEnviron.SplitBaryon())//Nucleus is below SplitFragmentThreshold
          {
#ifdef debug
            G4cout<<"G4QEnv::Const:Photon's added to Output, Env="<<theEnvironment<<G4endl;
#endif      
            G4QHadron* photon = new G4QHadron(opHad); // Fill projPhoton to Output
#ifdef debug
            G4cout<<"**G4QE::Const:Phot="<<photon->GetQC()<<photon->Get4Momentum()<<G4endl;
#endif
            theQHadrons.push_back(photon);      // (delete equivalent)
            return;
          }
          else if(targM<=999.)                  // Target is a nucleon
          {
            G4LorentzVector prot4m(0.,0.,0.,targM); // Prototype of secondary proton 4mom
            G4LorentzVector gam4m(0.,0.,0.,0.);     // Prototype for secondary gamma 4mom
            if(!G4QHadron(tot4Mom).DecayIn2(prot4m,gam4m))
            {
#ifdef debug
              G4cout<<"*War*G4QEnv::Const:(P)Photon->Output, Env="<<theEnvironment<<G4endl;
#endif      
              G4QHadron* photon = new G4QHadron(opHad); // Fill projPhoton to Output
#ifdef debug
              G4cout<<"**G4QE::Const:Ph="<<photon->GetQC()<<photon->Get4Momentum()<<G4endl;
#endif
              theQHadrons.push_back(photon);    // (delete equivalent)
              return;
            }
            G4QHadron* proton = new G4QHadron(targPDG,prot4m); // Fill tgProton to Output
            theQHadrons.push_back(proton);      // (delete equivalent)
            G4QHadron* photon = new G4QHadron(22,gam4m);       // Fill prPhoton to Output
            theQHadrons.push_back(photon);      // (delete equivalent)
            theEnvironment.InitByPDG(90000000); // Create nuclear environment
#ifdef debug
            G4cout<<"G4QEnv::Const:Fill gamma and N from gam+N"<<targPDG<<prot4m<<G4endl;
#endif      
            return;
          }
        }
      }
      else if(opPDG==13 || opPDG==15)
      {
        G4int         nuPDG=14;
        if(opPDG==15) nuPDG=16;
        G4LorentzVector mu4m=opHad->Get4Momentum();
        //G4double qpen=-180.*log(G4UniformRand()); // Energy of target-quark-parton(T=180)
        G4double qpen=465.*sqrt(sqrt(G4UniformRand())); // UniformDistr for 3-q nucleon
        G4double qpct=2*G4UniformRand()-1.;         // Cos(thet) of target-quark-parton
        G4double qpst=sqrt(1.-qpct*qpct);           // Sin(theta) of target-quark-parton
        G4double qppt=qpen*qpst;                    // PT of target-quark-parton
        G4double qphi=twopi*G4UniformRand();        // Phi of target-quark-parton
        G4LorentzVector qi4m(qppt*sin(qphi),qppt*cos(qphi),qpen*qpct,qpen); // quark-parton
        G4LorentzVector qt4m=mu4m+qi4m;             // Total 4mom (iniQP+lepton)
        G4LorentzVector nu4m(0.,0.,0.,0.);          // Prototype of secondary neutrino 4mom
        G4LorentzVector qf4m(0.,0.,0.,0.);          // Prototype for secondary quark-parton
        G4QContent targQC=targQPDG.GetQuarkContent(); // QC of the target nucleus (local!)
        targQC+=G4QContent(1,0,0,0,1,0);      // Make iso-shift with fake pi- is added
        G4LorentzVector fn4m=G4LorentzVector(0.,0.,0.,0.); // Prototype of the residual 4M
        G4QNucleus fnN(targQC,fn4m);          // Define the final state nucleus
        G4double   fnm=fnN.GetMZNS();         // GS Mass of the final state nucleus
        //G4QContent resiQC=targQC-neutQC; // QC of resid nucleus (-neutron)
        //G4QNucleus rsN(resiQC,fn4m);          // Define the final state nucleus
        //G4double   rsm=rsN.GetMZNS()+mNeut;   // GS Mass of residual nucleus w/o neutron
        G4double   tm=0.;                     // Prototype of RealMass of the final nucleus
        G4LorentzVector tg4m=G4LorentzVector(0.,0.,0.,targM); // 4mom of all target nucleus
        G4LorentzVector fd4m=tg4m-qi4m;       // 4mom of the residual coloured nuclear sys.
#ifdef debug
        //G4cout<<"-->>G4QEnv::Const:rM="<<rsm<<",fM="<<fnm<<",tM="<<targM<<G4endl;
        G4cout<<"G4QEnvironment::Const:mu4M="<<mu4m<<",t4M="<<qt4m<<",tgQP="<<qi4m<<G4endl;
#endif      
        while (tm<=fnm)
        {
          if(!G4QHadron(qt4m).DecayIn2(nu4m,qf4m))
          {
            G4cout<<"***G4QE::Constr:Muon error (1) 4M="<<mu4m<<". Fill as it is."<<G4endl;
            G4QHadron* lepton = new G4QHadron(opHad); // Fill projMuon to Output
            theQHadrons.push_back(lepton);        // (delete equivalent)
            return;
          }
#ifdef mudebug
          G4cout<<"G4QEnv::Const:i="<<qi4m<<",t="<<qt4m<<"->n="<<nu4m<<"+q="<<qf4m<<G4endl;
#endif
          fn4m=fd4m+qf4m;
          tm=fn4m.m();                 // Real mass of the final nucleus
#ifdef mudebug
          G4cout<<"--G4QEnv::Const:M="<<tm<<",GSM=="<<fnm<<G4endl;
#endif
        }
        fnN.Set4Momentum(fn4m);
        // (mu,q->nu,q) reaction succeded and Neutrino can be pushed to Output
        G4QHadron* neutrino = 0;              // NeutrinoPrototype to be filled to Output
#ifdef mudebug
        G4cout<<"G4QEnv::Const:fM="<<tm<<fn4m<<",GSM="<<fnm<<G4endl;
#endif      
        if(tm<fnm)                            // Final Nucleus is below the GS threshold
        {
          qf4m=G4LorentzVector(0.,0.,0.,fnm); // Final nucleus 4M for the final decay
          qt4m=tg4m+mu4m;
          if(!G4QHadron(qt4m).DecayIn2(nu4m,qf4m)) // Decay in Nucleus+nu_mu
          {
            G4cout<<"***G4QE::Constr:Muon error (2) 4M="<<mu4m<<". Fill as it is."<<G4endl;
            G4QHadron* muon = new G4QHadron(opHad); // Fill projMuon to Output
            theQHadrons.push_back(muon);        // (delete equivalent)
            return;
          }
          G4QHadron* fnuc = new G4QHadron(targQC,qf4m); // Fill Final Nucleus to Output
          //theQHadrons.push_back(fnuc);      // (delete equivalent)
          EvaporateResidual(fnuc);            // Try to evaporate residual (del. equiv.)
          neutrino = new G4QHadron(nuPDG,nu4m);// Fill Neutrino to Output
          theEnvironment.InitByPDG(90000000); // Create nuclear environment
#ifdef debug
          G4cout<<"G4QEnv::Const:Fill neutrino (1) "<<nuPDG<<nu4m<<G4endl;
#endif      
          theQHadrons.push_back(neutrino);    // (delete equivalent)
          return;
        }
        neutrino = new G4QHadron(nuPDG,nu4m); // Fill Neutrino to Output
#ifdef debug
        G4cout<<"G4QEnv::Const:Fill neutrino (2) "<<nuPDG<<nu4m<<G4endl;
#endif      
        theQHadrons.push_back(neutrino);      // (delete equivalent)
        if(tm<fnm+135.98)                     // FinalNucleus is below thePionThreshold(HE)
        {
          if(!fnN.SplitBaryon()) // Final Nucleus is below the splittingFragmentThreshold
          {
#ifdef mudebug
            G4cout<<"G4QEnv::Const: impossible to split nucleon after mu->nu"<<G4endl;
#endif      
            G4LorentzVector ga4m(0.,0.,0.,0.);
            qf4m=G4LorentzVector(0.,0.,0.,fnm);// Final nucleus 4M for the final decay
            if(!G4QHadron(fn4m).DecayIn2(ga4m,qf4m)) // Decay in Nucleus+photon
            {
              G4cout<<"***G4QE::Constr:LepCapError(3),M="<<fn4m.m()<<"<"<<fnm<<G4endl;
              G4QHadron* resid = new G4QHadron(targQC,qt4m); // Fill ResidNucleus to Output
              theQHadrons.push_back(resid);   // (delete equivalent)
              theEnvironment.InitByPDG(90000000);// Create nuclear environment
              return;
            }
            G4QHadron* photon = new G4QHadron(22,ga4m); // Fill projPhoton to Output
#ifdef debug
            G4cout<<"G4QEnv::Const:Fill photon "<<ga4m<<G4endl;
#endif      
            theQHadrons.push_back(photon);    // (delete equivalent)
            G4QHadron* fnuc = new G4QHadron(targQC,qf4m); // Fill Final Nucleus to Output
#ifdef debug
            G4cout<<"G4QEnv::Const:Fill target "<<targQC<<qf4m<<" in any form"<<G4endl;
#endif      
            EvaporateResidual(fnuc);          // Try to evaporate residual (del. equiv.)
            theEnvironment.InitByPDG(90000000);// Create nuclear environment
            return;
          }
        }
        // At this poin it is possible to convert mu- to pi-
        fn4m=qf4m-qi4m;
        opHad->SetQPDG(pimQPDG);              //Convert (mu-)u->d to (virt pi-)u->d capture
        fake=false;                           // normal pi- for q-muon scattering
        //fake=true;                          // fake pi- for q-muon scattering *****
        //if(G4UniformRand()>.5) fake=false;  // normal pi- for q-muon scattering *****
        opHad->Set4Momentum(fn4m);
      }
    }
    for(G4int ih=0; ih<nHadrons; ih++)        // ==> The main LOOP over projQHadrons
    {
      G4QHadron* curHadr=projHadrons[ih];     // Pointer to current projectile Hadron
      G4int hNFrag = curHadr->GetNFragments();// #0 means intermediate (skip)
      G4LorentzVector ch4M=curHadr->Get4Momentum(); // 4-momenyum of the current projectile
#ifdef debug
      G4cout<<"G4QE:C:"<<ih<<",F="<<hNFrag<<",0="<<projHadrons[0]->GetNFragments()<<G4endl;
#endif
      if(!hNFrag&&ch4M.e()>0.)                // => "Final hadron" case
      {
        G4int envPDG=theEnvironment.GetPDG();
        if(envPDG==90000000||(theEnvironment.Get4Momentum().m2())<1.) // ==>"Vacuum"
        {
          G4int hPDG  = curHadr->GetPDGCode();// A PDG Code of the projQHadron
          //if(!hPDG||hPDG==10)        // Check for the validity of the QHadron (@@ 10 OK?)
          if(!hPDG)          // Check for the validity of the QHadron
          {
            //G4cerr<<"--Warning--G4QEnvironment::Constructor: wrong PDG("<<ih<<")="<<hPDG
            //    <<", HQC="<<curHadr->GetQC()<<", HM="<<curHadr->GetMass()<<G4endl;
            //throw G4QException("***G4QEnvironment::Constructor: theInputHadron is Chip");
          }
          else
          {
            G4int hQ = curHadr->GetQCode();  // One more check for valid of the QHadron
            if(hQ<0)
            {
              //G4cerr<<"--Warning--G4QEnv::Constructor:Q<0, PDG=("<<ih<<")"<<hPDG<<G4endl;
              //throw G4QException("***G4QEnvironment::Constructor:theInputHadron is bad");
            }
            else
            {
              G4QHadron* newHadr = new G4QHadron(curHadr);
#ifdef debug
              G4cout<<"*G4QE::Const:H="<<newHadr->GetQC()<<newHadr->Get4Momentum()<<G4endl;
#endif
              theQHadrons.push_back(newHadr); // Fill existing hadron (delete equivalent)
#ifdef debug
              G4cout<<"G4QEnviron::Constructor: Fill h="<<hPDG<<ch4M<<G4endl;
              for(unsigned ipo=0; ipo<theQHadrons.size(); ipo++) // LOOP just for printing
              {
                G4int           hPDG  = theQHadrons[ipo]->GetPDGCode();
                G4LorentzVector h4Mom = theQHadrons[ipo]->Get4Momentum();
                G4int           hNFrag= theQHadrons[ipo]->GetNFragments();
                G4QContent      hQC   = theQHadrons[ipo]->GetQC();
                G4cout<<"h#"<<ipo<<": "<<hPDG<<hQC<<",4M="<<h4Mom<<",nFr="<<hNFrag<<G4endl;
              }
#endif
            } // End of Q-Code check
          } // End of proper PDG for i-th Hadron
        }
        else                                  // Nuclear Environment still exists
        {
          G4QContent      hQC   = curHadr->GetQC();
#ifdef debug
          G4cout<<"G4QE::Const:CreateQuasm, 4M="<<ch4M<<",QC="<<hQC<<",E="<<envPDG<<",tC="
                <<totCharge<<",tB="<<totBaryoN<<G4endl;
#endif
          CreateQuasmon(hQC, ch4M, fake);
        } // End of Existing Nuclear Environment case
      } // End of final hadron case
    } // End of the LOOP over input hadrons
  } // End of nuclear target case (including neutron=90000001 & proton=90001000)
  else                                        // => "Unique hadron" case
  {
    // the nuclEnviron is already init'ed as vacuum + get the first hadron for interaction
    G4QHadron* curHadr=projHadrons[0];        // Pointer to the firstProjecHadron (checked)
    G4int hPDG  = curHadr->GetPDGCode();      // A PDG Code of the projQHadron
    if(!hPDG||hPDG==10)                       // Check for the validity of the QHadron
    {
      G4cout<<"---Warning---G4QEnvironment::Constructor:Vacuum,1st Hadron wrong PDG="<<hPDG
            <<", HQC="<<curHadr->GetQC()<<", HM="<<curHadr->GetMass()<<G4endl;
      //throw G4QException("***G4QEnvironment::Constructor: Fiest input Hadron is wrong");
    }
    else
    {
      G4int hQ = curHadr->GetQCode();         // One more check for valid of the QHadron
      if(hQ<0)
      {
        G4cout<<"---Warning---G4QEnviron::Constructor:Vacuum,Q<0, 1st HPDG="<<hPDG<<G4endl;
        //throw G4QException("***G4QEnvironment::Constructor:theFirstInputHadron's wrong");
      }
      else                                // Now we can get 4Mom &  QC of incedent particle
      {
        G4LorentzVector h4Mom = curHadr->Get4Momentum();
        G4QContent      hQC   = curHadr->GetQC();
        if(!targPDG||targPDG==10) G4cout<<"G4QEnv::CreateQ; (1) PDG="<<targPDG<<G4endl;
        G4QPDGCode      tQPDG(targPDG);
        G4int           tQ    = tQPDG.GetQCode();
        if(tQ<0||targPDG==10)
        {
          G4cout<<"---Warning---G4QEnv::Constructor:TrgQC<0, Chipo?,PDG="<<targPDG<<G4endl;
          //throw G4QException("***G4QEnvironment::Constructor: Target is wrong");
        }
        else                                 // Now we can create a unique Quasmon
        {
          h4Mom+=G4LorentzVector(0.,0.,0.,tQPDG.GetMass()); //Projectile + TargetHadron
          hQC+=tQPDG.GetQuarkContent();
#ifdef debug
          G4cout<<"G4QEnv::Const:VacHadrTarg="<<h4Mom<<hQC<<",E="<<theEnvironment<<G4endl;
#endif
          G4Quasmon* curQuasmon = new G4Quasmon(hQC, h4Mom);
          theQuasmons.push_back(curQuasmon); // Insert Quasmon or hadron/gamma (del. eq.)
        }
      } // End of Q-Code check
    } // End of proper PDG for i-th Hadron
    if(nHadrons>1) for(G4int ih=0; ih<nHadrons; ih++) // fill other Hadrons to Output
    {
      G4QHadron* newHadr = new G4QHadron(curHadr);
#ifdef debug
      G4cout<<"*G4QE::Const:#"<<ih<<","<<curHadr->GetQC()<<curHadr->Get4Momentum()<<G4endl;
#endif
      theQHadrons.push_back(newHadr);        // Fill existing hadron (delete equivalent)
    }
  } // End of Unique Hadron target treatment
#ifdef chdebug
  G4int finCharge=theEnvironment.GetCharge();
  G4int finBaryoN=theEnvironment.GetA();
  G4int nHad=theQHadrons.size();
  if(nHad) for(G4int ih=0; ih<nHad; ih++)
  {
    finCharge+=theQHadrons[ih]->GetCharge();
    finBaryoN+=theQHadrons[ih]->GetBaryonNumber();
  }
  G4int nQuas=theQuasmons.size();
  if(nQuas) for(G4int iq=0; iq<nQuas; iq++)
  {
    finCharge+=theQuasmons[iq]->GetCharge();
    finBaryoN+=theQuasmons[iq]->GetBaryonNumber();
  }
  if(finCharge!=totCharge || finBaryoN!=totBaryoN)
  {
    G4cout<<"*::*G4QEnv::C:(0) tC="<<totCharge<<",C="<<finCharge<<",tB="<<totBaryoN
          <<",B="<<finBaryoN<<",E="<<theEnvironment<<G4endl;
    if(nHad) for(G4int h=0; h<nHad; h++)
    {
      G4QHadron* cH = theQHadrons[h];
      G4cout<<"*:*G4QE::C:h#"<<h<<",QC="<<cH->GetQC()<<",PDG="<<cH->GetPDGCode()<<G4endl;
    }
    if(nQuas) for(G4int q=0; q<nQuas; q++)
    {
      G4Quasmon* cQ = theQuasmons[q];
      G4cout<<"::G4QE::C:q#"<<q<<",C="<<cQ->GetCharge()<<",QuarkCon="<<cQ->GetQC()<<G4endl;
    }
  }
#endif
} // End of the G4QEnvironment constructor


G4QEnvironment::G4QEnvironment(const G4QEnvironment &right)
{
  // theQHadrons (Vector)
  theQFScat = G4QFreeScattering::GetPointer();
  G4int nQH             = right.theQHadrons.size();
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
    G4QHadron* curQH    = new G4QHadron(right.theQHadrons[ih]);
#ifdef debug
    G4cout<<"G4QE::CopyByVal:cH#"<<ih<<","<<curQH->GetQC()<<curQH->Get4Momentum()<<G4endl;
#endif
    theQHadrons.push_back(curQH);            // (delete equivalent)
  }
  // theProjectiles (Vector)
  G4int nP              = right.theProjectiles.size();
  if(nP) for(G4int ip=0; ip<nP; ip++)
  {
    G4QHadron* curP     = new G4QHadron(right.theProjectiles[ip]);
    theProjectiles.push_back(curP);           // (delete equivalent)
  }
  theTargetPDG          = right.theTargetPDG;
  theWorld              = right.theWorld;
  nBarClust             = right.nBarClust;
  f2all                 = right.f2all;
  tot4Mom               = right.tot4Mom;
  totCharge             = right.totCharge;
  totBaryoN             = right.totBaryoN;

  // theQuasmons (Vector)
  G4int nQ              = right.theQuasmons.size();
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
    G4Quasmon* curQ     = new G4Quasmon(right.theQuasmons[iq]);
    theQuasmons.push_back(curQ);             // (delete equivalent)
  }

  // theQCandidates (Vector)
  G4int nQC             = right.theQCandidates.size();
  if(nQC) for(G4int ic=0; ic<nQC; ic++)
  {
    G4QCandidate* curQC = new G4QCandidate(right.theQCandidates[ic]);
    theQCandidates.push_back(curQC);         // (delete equivalent)
  }

  theEnvironment        = right.theEnvironment;
}

G4QEnvironment::G4QEnvironment(G4QEnvironment* right)
{
  // theQHadrons (Vector)
  theQFScat = G4QFreeScattering::GetPointer();
  G4int nQH             = right->theQHadrons.size();
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
    G4QHadron* curQH    = new G4QHadron(right->theQHadrons[ih]);
#ifdef debug
    G4cout<<"G4QE::CopyByPtr:cH#"<<ih<<","<<curQH->GetQC()<<curQH->Get4Momentum()<<G4endl;
#endif
    theQHadrons.push_back(curQH);            // (delete equivalent)
  }

  // theProjectiles (Vector)
  G4int nP              = right->theProjectiles.size();
  if(nP) for(G4int ip=0; ip<nP; ip++)
  {
    G4QHadron* curP     = new G4QHadron(right->theProjectiles[ip]);
    theProjectiles.push_back(curP);           // (delete equivalent)
  }
  theTargetPDG          = right->theTargetPDG;
  theWorld              = right->theWorld;
  nBarClust             = right->nBarClust;
  f2all                 = right->f2all;
  tot4Mom               = right->tot4Mom;
  totCharge             = right->totCharge;
  totBaryoN             = right->totBaryoN;

  // theQuasmons (Vector)
  G4int nQ              = right->theQuasmons.size();
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
    G4Quasmon* curQ     = new G4Quasmon(right->theQuasmons[iq]);
    theQuasmons.push_back(curQ);             // (delete equivalent)
  }

  // theQCandidates (Vector)
  G4int nQC             = right->theQCandidates.size();
  if(nQC) for(G4int ic=0; ic<nQC; ic++)
  {
    G4QCandidate* curQC = new G4QCandidate(right->theQCandidates[ic]);
    theQCandidates.push_back(curQC);         // (delete equivalent)
  }

  theEnvironment        = right->theEnvironment;
}

G4QEnvironment::~G4QEnvironment()
{
#ifdef debug
  G4cout<<"~G4QEnvironment: before theQCandidates nC="<<theQCandidates.size()<<G4endl;
#endif
  for_each(theQCandidates.begin(), theQCandidates.end(), DeleteQCandidate());
#ifdef debug
  G4cout<<"~G4QEnvironment: before theQuasmons nQ="<<theQuasmons.size()<<G4endl;
#endif
  for_each(theQuasmons.begin(), theQuasmons.end(), DeleteQuasmon());
#ifdef debug
  G4cout<<"~G4QEnvironment: before theQHadrons nH="<<theQHadrons.size()<<G4endl;
#endif
  for_each(theQHadrons.begin(), theQHadrons.end(), DeleteQHadron());
#ifdef debug
  G4cout<<"~G4QEnvironment: before theProjectiles nP="<<theProjectiles.size()<<G4endl;
#endif
  for_each(theProjectiles.begin(), theProjectiles.end(), DeleteQHadron());
#ifdef debug
  G4cout<<"~G4QEnvironment: === DONE ==="<<G4endl;
#endif
}

G4double G4QEnvironment::SolidAngle=0.8;     // Part of Solid Angle to capture (@@A-dep.)
G4bool   G4QEnvironment::EnergyFlux=false;   // Flag for Energy Flux use (not MultyQuasmon)
G4bool   G4QEnvironment::WeakDecays=false;   // Flag for WeakDecays(notUsed hadwaredClosed)
G4bool   G4QEnvironment::ElMaDecays=true;    // Flag for Electromagnetic decays of hadrons
G4double G4QEnvironment::PiPrThresh=141.4;   // Pion Production Threshold for gammas
G4double G4QEnvironment::M2ShiftVir=20000.;  // Shift for M2=-Q2=m_pi^2 of the virtualGamma
G4double G4QEnvironment::DiNuclMass=1880.;   // DoubleNucleon Mass for VirtualNormalization

// Open decay of particles with possible electromagnetic channels of decay (gammas)
void G4QEnvironment::OpenElectromagneticDecays(){ElMaDecays=true;}

// Close decay of particles with possible electromagnetic channels of decay (gammas)
void G4QEnvironment::CloseElectromagneticDecays(){ElMaDecays=false;}

// Fill the private static parameters
void G4QEnvironment::SetParameters(G4double solAn, G4bool efFlag, G4double piThresh,
                                   G4double mpisq, G4double dinum)
{
  EnergyFlux=efFlag;       // Flag for Energy Flux use instead of Multy Quasmon
  SolidAngle=solAn;        // Part of Solid Angle to capture secondaries (@@A-dep)
  PiPrThresh=piThresh;     // Pion Production Threshold for gammas
  M2ShiftVir=mpisq;        // Shift for M2=-Q2=m_pi^2 of the virtual gamma
  DiNuclMass=dinum;        // Double Nucleon Mass for virtual normalization
}

const G4QEnvironment& G4QEnvironment::operator=(const G4QEnvironment &right)
{
  if(this != &right)                          // Beware of self assignment
  {
    // theQHadrons (Vector)
    theQFScat = G4QFreeScattering::GetPointer();
    G4int iQH             = theQHadrons.size();
    if(iQH) for(G4int ii=0; ii<iQH; ii++) delete theQHadrons[ii];
    theQHadrons.clear();
    G4int nQH             = right.theQHadrons.size();
    if(nQH) for(G4int ih=0; ih<nQH; ih++)
    {
      G4QHadron* curQH    = new G4QHadron(right.theQHadrons[ih]);
#ifdef debug
      G4cout<<"G4QE::Operator=:c#"<<ih<<","<<curQH->GetQC()<<curQH->Get4Momentum()<<G4endl;
#endif
      theQHadrons.push_back(curQH);            // (delete equivalent)
    }

    theWorld              = right.theWorld;
    nBarClust             = right.nBarClust;
    f2all                 = right.f2all;

    // theQuasmons (Vector)
    G4int iQ              = theQuasmons.size();
    if(iQ) for(G4int jq=0; jq<iQ; jq++) delete theQuasmons[jq];
    theQuasmons.clear();
    G4int nQ              = right.theQuasmons.size();
    if(nQ) for(G4int iq=0; iq<nQ; iq++)
    {
      G4Quasmon* curQ     = new G4Quasmon(right.theQuasmons[iq]);
      theQuasmons.push_back(curQ);             // (delete equivalent)
    }

    // theQCandidates (Vector)
    G4int iQC             = theQCandidates.size();
    if(iQC) for(G4int jc=0; jc<iQC; jc++) delete theQCandidates[jc];
    theQCandidates.clear();
    G4int nQC             = right.theQCandidates.size();
    if(nQC) for(G4int ic=0; ic<nQC; ic++)
    {
      G4QCandidate* curQC = new G4QCandidate(right.theQCandidates[ic]);
      theQCandidates.push_back(curQC);        // (delete equivalent)
    }

    theEnvironment        = right.theEnvironment;
  }
  return *this;
}

// Member function for Quasmon Creation & Environment nucleus modification
void G4QEnvironment::CreateQuasmon(const G4QContent& projQC, const G4LorentzVector& pro4M,
                                   G4bool fake)
{
  //static const G4double third=1./3.;
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mNeu2 = mNeut*mNeut;
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mPro2 = mProt*mProt;
  //static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mPi2 = mPi*mPi;
  //static const G4double mMu  = G4QPDGCode(13).GetMass();
  //static const G4double mMu2 = mMu*mMu;
  //static const G4QContent gamQC(0,0,0,0,0,0);
  //static const G4QContent pimQC(1,0,0,0,1,0);
  //static const G4QContent pipQC(0,1,0,1,0,0);
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent deutQC(3,3,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4QNucleus vacuum(90000000);
  G4QContent valQ(0,0,0,0,0,0);             // Prototype of the Quasmon's Quark Content
  G4LorentzVector q4Mom(0.,0.,0.,0.);       // Prototype of the Quasmon's 4-momentum
  nBarClust = 1;                            // By default only quasi-free nucleons
  G4LorentzVector proj4M=pro4M;             // Fake equivalence to avoid & const
  G4double  projE=proj4M.e();               // energy of the projectile
  G4int projPDG=projQC.GetSPDGCode();       // Minimum hadron for the projectile QC
  if(projE<0.)
  {
    G4cout<<"*Warning*G4QEnvironment::CreateQuasmon:Epr="<<projE<<"<0,QC="<<projQC<<G4endl;
    projE=0.;
    proj4M=G4LorentzVector(0.,0.,0.,0.);
  }
  G4double  projM2=proj4M.m2();             // projectile's squared mass (print & v.gamma)
  G4bool Pr1 = theProjectiles.size() == 1;  // A#ofHadronsInTheInputVector = 1 condition
  if(Pr1 && std::fabs(G4QPDGCode(projPDG).GetMass2()-projM2) > .1 ) Pr1=false; // MassShell
  G4int     targPDG=theEnvironment.GetPDG();// PDG Code of the target nucleus
  if(targPDG>80000000&&targPDG!=90000000&&(theEnvironment.Get4Momentum().m2())>1000.)//Nucl
  {
    G4double  tgMass=theEnvironment.GetMass();// mass of the target (QEnvironment) nucleus
#ifdef debug
    G4cout<<"G4QEnvironment::CreateQ:Interact "<<projQC<<proj4M<<"(m2="<<projM2<<") + A="
          <<targPDG<<",M="<<tgMass<<",tC="<<totCharge<<",tB="<<totBaryoN<<G4endl;
#endif
    G4int envZ=theEnvironment.GetZ();       // A#of protons in the target nucleus
    G4int envN=theEnvironment.GetN();       // A#of neutrons in the target nucleus
    G4int envS=theEnvironment.GetS();       // A#of lambdas in the target nucleus
    G4int envBN=envZ+envN+envS;             // A baryon number of the target nucleus
    G4int nP  =theWorld->GetQPEntries();    // A#of initialized particles in CHIPS World
    //G4int nCl =nP-90;                       // A#of initialized clusters in CHIPS World
    G4int nCl =nP-53;           // @@ A#ofClusters in CHIPSWorld (53=nQHM in G4QPDGCode.hh)
    if(nCl<0) G4cout<<"---Warning---G4QEnv::CreaQ:nP="<<nP<<" for Targ="<<targPDG<<G4endl;
    if     (nCl<3) nBarClust=1;             // Fix the maximum Baryon Number for clusters
    else if(nCl<9) nBarClust=2;
    else
    {
      G4int v=nCl-9;
      G4int d=v/15;
      G4int r=v%15;
      if(r<7) nBarClust=3+d+d;
      else    nBarClust=4+d+d;
    }
#ifdef debug
    G4cout<<"G4QE::CrQ:TNuc:Z="<<envZ<<",N="<<envN<<",nC="<<nBarClust<<",tC="
          <<totCharge<<", tB="<<totBaryoN<<G4endl;
#endif
    G4bool pbpt=projE<PiPrThresh+(M2ShiftVir+projM2)/DiNuclMass;// PhotonBelowPionThreshold
    G4bool din=false;
    G4bool piF=false;
    G4bool gaF=false;
    if((projM2-mPi2<.00001||projE-mPi<2.7)&&projPDG==-211&&!fake) piF=true;//PiAtRestCase
    //if(pbpt&&projPDG==22) din=true; // InCaseOf GammaBelowPiThresh needs DiNucl (?)
    if(pbpt&&projPDG==22) gaF=true; // InCaseOf GammaBelowPiThresh needs DiNucl (?)
    theEnvironment.SetMaxClust(nBarClust);
    nBarClust=theEnvironment.UpdateClusters(din); // Cluster Probabilities upto maxClust
#ifdef debug
    G4cout<<"G4QEnv::CreateQ: Nucleus("<<targPDG<<") is created ("<<nBarClust<<" clast's)";
    for(G4int ic=0;ic<nBarClust;ic++)
      G4cout<<" #"<<ic<<"("<<theEnvironment.GetProbability(ic)<<")";
    G4cout<<G4endl;
#endif
    theEnvironment.PrepareCandidates(theQCandidates,piF,gaF,proj4M);//Calc.Clust's probab's
    G4QNucleus memEnviron=theEnvironment;
#ifdef debug
    G4cout<<"G4QE::CrQ:ClusterProbabCalculation tC="<<totCharge<<",tB="<<totBaryoN<<G4endl;
#endif

    G4bool efFlag = false;     // EnergyFlowFlag=FALSE (@@=DEFOLT=@@ make par)
    // ***** Change if necessary to compare Energy Flux & Multy Quasmon ******

    G4int efCounter=0;                      // Counter of Energy Flux particles
    G4QContent EnFlQC(0,0,0,0,0,0);         // Quark Content of Energy Flux
    G4LorentzVector ef4Mom(0.,0.,0.,0.);    // Summed 4-momentum of Energy Flux
    G4double proj3M=proj4M.rho();
    // --- P-antibar ---  N-antibar  -- LAMBDA-antibar - SIGMA-antibar - SIGMA0-antibar
    if((projPDG==-2212||projPDG==-2112||projPDG==-3122||projPDG==-3112||projPDG==-3212||
        projPDG==-3222) && envBN>1 && proj3M<10.) // OnlyForAtRestReactions(@@to Interface)
    // ^ SIGMA+ antibar
    {
      // @@ Annihilation on one baryon is implemented (no annihilation on clusters! @@?) @@
#ifdef debug
      G4cout<<"G4QE::CreQ:Annihilation on a perif. nucleon, Z="<<envZ<<",N="<<envN<<G4endl;
#endif
      G4double   zpn=envZ+envN;             // a#of nucleons in the nucleus
      G4double   rnd=zpn*G4UniformRand();   // Random number to find a nucleon
      //G4double   envD=.1*envZ*envN/zpn;   // a#of possible quasifree deuterons (@@Param.)
      //G4double   rnd=(zpn+envD)*G4UniformRand(); // Random number to find a baryon
      //G4double   rnd=(zpn+envS)*G4UniformRand(); // Random number to find a baryon
      G4int      targNPDG = 90000000;       // Nucl-Prototype of PDG of Periferal Target
      G4QContent targQC(0,0,0,0,0,0);       // Quark Content of Periferal Target
      if     (rnd<envN)                     // Neutron is a Periferal Target
      {
        targNPDG = 90000001;
        targQC   = neutQC;
      }
      else
      //     if(rnd<=zpn)                      // Proton is a Periferal Target
      {
        targNPDG = 90001000;
        targQC   = protQC;
      }
      //else                                  // Deuteron as a Periferal Target
      //{
      //  targNPDG = 90001001;
      //  targQC   = deutQC;
      //}
      //else                                  // Lambda is a Periferal Target (?)
      //{
      //  targNPDG = 91000000;
      //  targQC   = lambQC;
      //}
      theEnvironment.Reduce(targNPDG);      // Subtract periferal baryon from Nucleus
      G4double resMass=theEnvironment.GetGSMass(); // Nuclear mass after baryon subtraction
      G4double barMass=tgMass-resMass;      // Mass of the bound baryon for annihilation
      tgMass=resMass;                       // New mass of theEnvironment
      q4Mom=G4LorentzVector(0,0,0,barMass)+proj4M;// 4-mom of the intermediate B-Bbar Quasm
      valQ=targQC+projQC;                   // Quark Content of intermediate B-Bbar Quasmon
#ifdef debug
      G4cout<<"G4QEnviron::CQ:"<<targNPDG<<" + Env="<<theEnvironment<<",QC="<<valQ<<G4endl;
#endif
      // Remember the Quasmon parameters, defined by user for recovery after annihilation
      G4Quasmon fakeQ;                      // fake Quasmon to get and restore parameters
      G4double QTemper=fakeQ.GetTemper();   // Temperature defined by user for Quasmons
      G4double QSOverU=fakeQ.GetSOverU();   // S/U defined by user for Quasmons
      G4double QEtaSup=fakeQ.GetEtaSup();   // Eta Suppresion defined by user in Quasmons
      G4Quasmon::SetParameters(180.,QSOverU,.3); // Parameters for N-barN Annihilation
      G4Quasmon::CloseElectromagneticDecays();   // Parameters for N-barN Annihilation
      G4Quasmon* pan = new G4Quasmon(valQ,q4Mom);// N-Nbar Quasm creation (del.at 9th line)
      G4QNucleus vE(90000000);                   // Annihilation in vacuum (in NuclMatter?)
#ifdef debug
      G4cout<<"G4QE::CreQ: before Fragment, vE="<<vE<<",vP="<<vE.GetProbability()<<",QQC="
            <<valQ<<",Q4M="<<q4Mom<<G4endl;
#endif
      G4QHadronVector* output=pan->Fragment(vE,1);//Output of inVacAnnihilation*DESTROY*<-+
#ifdef debug
      G4cout<<"G4QE::CrQ:NucleonAntinucleonAnnihilation's done,N="<<output->size()<<G4endl;
#endif
      G4Quasmon::OpenElectromagneticDecays();  // Parameter for multihadronFragmentatation^
#ifdef debug
      G4cout<<"G4QE::CrQ:>>AnnihilationIsDone,C="<<totCharge<<",B="<<totBaryoN<<G4endl;// ^
#endif
      delete pan;                              // The N-NbarQuasmon is deleted A.S.A.P.   ^
      G4QHadronVector input;                   // Input for MultyQuasmon **DESTROY**<---+ ^
      //G4int trgPDG = theEnvironment.GetPDG();// New PDG Code for the Residual Nucleus ^ ^
      G4LorentzVector trg4M(0.,0.,0.,resMass); // New 4-momentum for the ResidualNucleus^ ^
      G4int tNH = output->size();              // For the selection LOOP                ^ ^
      G4ThreeVector dir = G4RandomDirection(); // For the selection in LOOP (@@ at rest)^ ^
      //G4double ra=std::pow(G4double(totBaryoN),third);  //                            ^ ^
      G4double ra=G4QThd(totBaryoN);  //                                                ^ ^
#ifdef debug
      G4cout<<"G4QE::CQ:N="<<tNH<<",T="<<totCharge<<","<<totBaryoN<<",A="<<ra<<G4endl;//^ ^
#endif
      for (G4int ind=0; ind<tNH; ind++)        // Loop over annihilation  QHadrons      ^ ^
      {
        //G4QHadron* curHadr = output->operator[](ind); // Pointer to theCurrentHadron  ^ ^
        G4QHadron*      curHadr = (*output)[ind];       // Pointer to theCurrentHadron  ^ ^
        G4int           shDFL= curHadr->GetNFragments();// A#of decFragments for proj.  ^ ^
        G4LorentzVector sh4m = curHadr->Get4Momentum(); // 4Mom for the projectile      ^ ^
        G4ThreeVector   shDIR= sh4m.vect().unit();      // unitVector in projMomDirect  ^ ^
        G4int           shPDG= curHadr->GetPDGCode();   // PDG Code of the projectile   ^ ^
        G4int           shCHG= curHadr->GetCharge();    // Charge of the projectile     ^ ^
        G4double        shMOM= sh4m.rho();              // Momentum of the projectile   ^ ^
#ifdef debug
        G4cout<<"G4QE::CrQ:"<<ind<<","<<shDFL<<",PDG="<<shPDG<<",4M="<<sh4m<<G4endl; // ^ ^
#endif
        G4double solAnCut=SolidAngle;                   // Proto ChargeDependantSolAngle^ ^
        if(fabs(ra)<.1) solAnCut=3.;                    // No Nucleus -> no absorption  ^ ^
        else if(shMOM<.1)                               // Meson at rest                ^ ^
        {                                               //                              ^ ^
          if(shCHG<=0.) solAnCut=-3.;                   // Always totally absorbed      ^ ^
          else          solAnCut= 3.;                   // Positive are repelled from A ^ ^
        }                                               //                              ^ ^
        else     solAnCut+=1000*shCHG/shMOM/ra;         // ChargeDepSolAngle(Normal)    ^ ^
        //G4double solAnCut=SolidAngle+20*shCHG*sqrt(1.*envZ)/shMOM;//ChargeDepSolAngle ^ ^
#ifdef debug
        G4cout<<"G4QE::CrQ: PDG="<<shPDG<<", p="<<shMOM<<", r="<<ra<<G4endl; //         ^ ^
#endif
        if(!shDFL)                                      // Final(notDecayed) hadrons    ^ ^
        {
#ifdef debug
          G4cout<<"G4QE::CQ:>H="<<shPDG<<":"<<dir.dot(shDIR)<<">"<<solAnCut<<G4endl; // ^ ^
#endif
          //if((dir.dot(shDIR)>solAnCut||shMOM<120.) && abs(shPDG)>99) // Absorb mesons ^ ^
          if(dir.dot(shDIR)>solAnCut && abs(shPDG)>99)  // Absorb mesons                ^ ^
          {
#ifdef debug
            G4cout<<"G4QE::CQ:>H="<<shPDG<<":"<<dir.dot(shDIR)<<">"<<solAnCut<<", P="// ^ ^
                   << shMOM <<" < 120" << G4endl;                           //          ^ ^
#endif
            //if (efFlag)
            //{ // => Case of "Energy Flux approach"   *** is temporary closed ***      ^ ^
            //  G4QContent shQC = curHadr->GetQC();     //                              ^ ^
            //  // QuarkContent of the Current Hadron                                   ^ ^
            //  ef4Mom+=sh4m;                           //                              ^ ^
            //  EnFlQC+=shQC;                           //                              ^ ^
            //  efCounter++;                            //                              ^ ^
#ifdef debug
		//  G4int hPDG=curHadr->GetPDGCode();       // Only for gebug printing      ^ ^
            //  G4LorentzVector h4M = curHadr->Get4Momentum();  // For debug printing   ^ ^
            //  G4cout<<"G4QE::CrQ:#"<<efCounter<<",PDG="<<hPDG<<",h4M="<<h4M<<G4endl;//^ ^
#endif
            //} //                                                                      ^ ^
            // else                                //=>"MultyQuasFragmentation"(!efFlag)^ ^
            {
              G4QHadron* mqHadron = new G4QHadron(curHadr);
              input.push_back(mqHadron);           // Fill hadron-copy (del equiv)
#ifdef debug
              G4int hPDG=curHadr->GetPDGCode();    // Only for debug printing           ^ ^
              G4LorentzVector h4M = curHadr->Get4Momentum(); // Only for gebug printing ^ ^
              G4cout<<"G4QE::CrQ:Absorb#"<<ind<<", PDG="<<hPDG<<", h4M="<<h4M<<G4endl;//^ ^
#endif
	      }                                      //                                   ^ ^
          }                                        //                                   ^ ^
          else                                     // DirectFilling of the output vector^ ^
          {                                        //                                   ^ ^
#ifdef debug
            G4int hPDG=curHadr->GetPDGCode();      // Only for gebug printing           ^ ^
            G4LorentzVector h4M = curHadr->Get4Momentum(); // Only for gebug printing   ^ ^
            G4cout<<"G4QE::CrQ: Fill OUT #"<<ind<<",PDG="<<hPDG<<",h4M="<<h4M<<G4endl;//^ ^
#endif
            // Just fill a hadron to the output stack (Make EM decays elsewhere)        ^ ^
            G4QHadron* curHadron = new G4QHadron(curHadr); //                           ^ ^
            theQHadrons.push_back(curHadron);      // TheQHadrs are filled as new Hadrs ^ ^
          }
        } // End of the IF over projectiles                                             ^ ^
      } // End of LOOP over "output" of annihilation                                    ^ ^
      for_each(output->begin(), output->end(), DeleteQHadron());     //DESTRUCT output>-^-^
      output->clear();                             //                                   ^ ^
      delete output;                               // ----------------------------------^-*
      if(!efFlag)                                  // =>NotEnergyFlux=MultyQuasmon Case ^
      {
        G4int noh = theQHadrons.size();            // a#oh hadrons in Output UpToNow    ^
        if(noh) for(G4int kh=0; kh<noh; kh++)      // One can escape it but...          ^
        {                                          //                                   ^
#ifdef debug
          G4cout<<"G4QE::CreateQ:H#"<<kh<<", QC="<<theQHadrons[kh]->GetQC() //          ^
                <<", 4M="<<theQHadrons[kh]->Get4Momentum()<<G4endl;         //          ^
#endif
          G4QHadronVector* tmpQHadVec=G4Quasmon().DecayQHadron(theQHadrons[kh]);//d.e<^ ^
          G4int tmpS=tmpQHadVec->size();           //                                 ^ ^
          intQHadrons.resize(tmpS+intQHadrons.size()); // Resize theQHadrons length   ^ ^
          copy(tmpQHadVec->begin(), tmpQHadVec->end(), intQHadrons.end()-tmpS); //    ^ ^
          tmpQHadVec->clear();                     //                                 ^ ^
          delete tmpQHadVec;           // who calls DecayQHadron must clear & delete  ^ ^
        }                                          //                                 ^ ^
        theQHadrons.clear(); // deletedWhenDecayed // Now theQHadrons is EmptyVector->^ ^
#ifdef debug
        G4int nInH=intQHadrons.size();             // Resulting #of hadrons after decay ^
        G4cout<<"G4QE::CrQ:nH="<<nInH<<",C="<<totCharge<<",B="<<totBaryoN<<G4endl;//    ^
#endif
        if(!(input.size()))                        // *RETURN* Without Quasmon creation-^
        {                                          //                                   ^
#ifdef debug
          G4cout<<"*G4QEnv::CrQ:AnnihStack tC="<<totCharge<<",tB="<<totBaryoN<<G4endl;//^
#endif
          return;                                  // Do not clear and delete objects --^
        }                                          //                                   ^
#ifdef debug
        G4cout<<"G4QE::CrQ:fakeQ, restPars tC="<<totCharge<<",tB="<<totBaryoN<<G4endl;//^
#endif
        G4Quasmon::SetParameters(QTemper,QSOverU,QEtaSup);//RecoverQParam's after anihil^
        G4Quasmon::OpenElectromagneticDecays(); // Parameter for multihadron fragmentat.^
        // From this point the new temporary environment is created (multiQuasmon)      ^
        G4QEnvironment* muq = new G4QEnvironment(input,theEnvironment.GetPDG());//<--+  ^
#ifdef debug
        G4cout<<"G4QE::CrQ:befCl&Dest tC="<<totCharge<<", tB="<<totBaryoN<<G4endl; //^  ^
#endif
        for_each(input.begin(), input.end(), DeleteQHadron());     //DESTROING inp >-^--^
        input.clear();                             // =Clearing=>---------->---------^==+
        theEnvironment = muq->GetEnvironment();    // RestoreResidEnv after interact.^  
        G4QHadronVector* outH = muq->GetQHadrons();// Copy of QHadrons *DESTROY* <---^-<--+
        G4QuasmonVector* outQ = muq->GetQuasmons();// Copy of Quasmons *DESTROY* <---^--+ ^
        delete muq;                                //----->----------->--------------^  ^ ^
        noh = outH->size();                        // a#of Not Interacting(Q) Hadrons   ^ ^
#ifdef debug
        G4cout<<"G4QEnv::CreateQ:*** #ofNotInterQH="<<noh<<" is found ***"<<G4endl; //  ^ ^
#endif
        if(noh) for(G4int nh=0; nh<noh; nh++)      // One can escape it but...          ^ ^
        {                                          //                                   ^ ^
#ifdef debug
          G4cout<<"G4QE::CreateQ: NotIntQH#"<<nh<<", QC="<<(*outH)[nh]->GetQC()  //     ^ ^
                <<", 4M="<<(*outH)[nh]->Get4Momentum()<<G4endl;                  //     ^ ^
#endif
          G4QHadronVector* tmpQHadVec=G4Quasmon().DecayQHadron((*outH)[nh]);//del.eq<-+ ^ ^
          G4int tmpS=tmpQHadVec->size();           //                                 ^ ^ ^
          intQHadrons.resize(tmpS+intQHadrons.size()); // Resize theQHadrons length   ^ ^ ^
          copy(tmpQHadVec->begin(), tmpQHadVec->end(), intQHadrons.end()-tmpS); //    ^ ^ ^
          tmpQHadVec->clear();                     //                                 ^ ^ ^
          delete tmpQHadVec;           // who calls DecayQHadron must clear & delete->+ ^ ^
        }                                          //                                   ^ ^
        outH->clear();                             //                                   ^ ^
        delete outH;                               // >---->---->---->---->---->---->---^-+
        G4int nMQ = outQ->size();                  // A#ofQuasmons in MultyQuasmonOutput^
#ifdef debug
        G4LorentzVector eLorV=theEnvironment.Get4Momentum(); //                         ^
        G4cout<<"G4QE::CrQ:nMQ="<<nMQ<<",tC="<<totCharge<<", tB="<<totBaryoN<<G4endl;// ^
        G4cout<<"G4QE::CrQ:Env4M="<<eLorV<<G4endl; //                                   ^
        G4LorentzVector contr4M=eLorV; //                                               ^
#endif
        if(nMQ) for(G4int mh=0; mh<nMQ; mh++)      // Can escape CreationDistruct but...^
        {                                          //                                   ^
          G4Quasmon* curQ = new G4Quasmon((*outQ)[mh]);// Copy to destroy TMP(?)        ^
#ifdef debug
          G4LorentzVector qLorV=curQ->Get4Momentum(); //                                ^
          G4cout<<"G4QE::CrQ:Q#"<<mh<<",4M="<<qLorV<<curQ->GetQC()<<G4endl; //          ^
          contr4M+=qLorV; //                                                            ^
#endif
          theQuasmons.push_back(curQ);             // Fill QuasmonCopies in theQuasmons ^
        }                                          //                                   ^
        for_each(outQ->begin(), outQ->end(), DeleteQuasmon()); // >-------------------->+
        outQ->clear();                             //                                   ^
        delete outQ;                               // >------------>------------------->+
#ifdef debug
        G4int nsHadr  = theQHadrons.size();      // Update the value of OUTPUT entries
        G4cout<<"G4QEnvironment::CreateQ: before return nH="<<nsHadr<<G4endl;
        if(nsHadr) for(G4int jso=0; jso<nsHadr; jso++)// LOOP over output hadrons 
        {
          G4int hsNF  = theQHadrons[jso]->GetNFragments(); // A#of secondary fragments
          if(!hsNF)                                        // Add only final hadrons
          {
            G4LorentzVector hLorV=theQHadrons[jso]->Get4Momentum();
            G4int           hPDGC=theQHadrons[jso]->GetPDGCode();
            G4cout<<"G4QE::CrQ: H#"<<jso<<",4M="<<hLorV<<hPDGC<<G4endl;
            contr4M+=hLorV;
          }
          else
            G4cout<<"G4Q::CrQ:"<<jso<<"NF=0,4M="<<theQHadrons[jso]->Get4Momentum()<<G4endl;
        }
        G4cout<<"G4QEnvironment::CreateQ: before return tot4M="<<contr4M<<G4endl;
#endif
        return;                                    // *** RETURN *** 
      }
      else                                         // ==> Energy Flux case
      {
        if (!efCounter) return;                    // ***RETURN*** Without Quasmon creation
      }
    }                                              // End of Hyperon annihilation case
    else EnFlQC=projQC;                            // For notAntiBar, don't use EnergyFlux
    G4double EnFlP=ef4Mom.rho();                   // Mom. of EnergyFlow forClusterCreation
    PrepareInteractionProbabilities(EnFlQC,EnFlP); // InteractionProbabilities for clusters
    G4int nCandid = theQCandidates.size();
#ifdef debug
    G4cout<<"G4QEnvironment::CrQ: InteractionProbabilities are done, nC="<<nCandid<<G4endl;
#endif
    if(nCandid<=0)
    {
      // G4cout<< "---Warning---G4QEnv::CreaQ:nC=" <<nCandid<< ",E=" <<theEnvironment<<G4endl;
      // throw G4QException("G4QEnvironment::CreateQ: Can not select a cluster");
      G4ExceptionDescription ed;
      ed << "Cannot select a cluster: nC=" << nCandid << ",E="
         << theEnvironment << G4endl;
      G4Exception("G4QEnvironment::CreateQuasmon()", "HAD_CHPS_0000",
                  FatalException, ed);
    }
    G4double maxP = theQCandidates[nCandid-1]->GetIntegProbability();
    G4int i=0;
    G4QContent    curQC;                           // Quark Content of the selected cluster
    if(nCandid==1||maxP==0.)
    {
#ifdef debug
      G4cout<<"***G4QEnv::CrQ:MaxP=0||nCand=1: Use all Env., Env="<<theEnvironment<<G4endl;
#endif
      curQC=theEnvironment.GetQCZNS();
      theEnvironment=vacuum;
    }
    else
    {
      G4double totP = maxP * G4UniformRand();
#ifdef debug
      G4cout<<"G4QEnvironment::CrQ:nC="<<nCandid<<", maxP="<<maxP<<", totP="<<totP<<G4endl;
#endif
      while(theQCandidates[i]->GetIntegProbability()<totP) i++;
      G4QCandidate* curCand = theQCandidates[i];// Pointer to selected cluster to interact
      curQC   = curCand->GetQC();               // Get QuarkContent of the selected cluster
      G4QNucleus targClust(curQC.GetP(),curQC.GetN(),curQC.GetL());//Define Clust as a QNuc
      G4double clMass=targClust.GetGSMass();    // Mass of residual nuclear environment
#ifdef cldebug
      G4cout<<"G4QEnv::CrQ:Cl#"<<i<<"(of "<<nCandid<<"),QC="<<curQC<<",M="<<clMass<<G4endl;
#endif
      G4LorentzVector pq4M=proj4M+G4LorentzVector(0.,0.,0.,clMass); 
      if(pq4M.m()>=clMass)
      {
#ifdef debug
        G4cout<<"G4QEnv::CQ:#"<<i<<"("<<targClust<<curQC<<") Env="<<theEnvironment<<G4endl;
#endif
        theEnvironment.Reduce(targClust.GetPDG());// Subtract selected cluster from Nucleus
      }
      else
      {
        G4double teMass=theEnvironment.GetGSMass(); //Mass of theResidualNuclearEnvironment
        G4LorentzVector te4M=proj4M+G4LorentzVector(0.,0.,0.,teMass);
        if(te4M.m()>=teMass)
        {
#ifdef debug
          G4cout<<"***G4QEnv::CrQ: Deep virtual, use all Env,Env="<<theEnvironment<<G4endl;
#endif
          curQC=theEnvironment.GetQCZNS();
          theEnvironment=vacuum;
        }
        else
        {
          G4QHadron* projH = new G4QHadron(projQC,proj4M);
          theQHadrons.push_back(projH);
          G4cout<<"---Warning---G4QE::CrQ:Fill Proj asItIs QC/4m="<<projQC<<proj4M<<G4endl;
          return;
        }
      }
    }
    G4double envMass=theEnvironment.GetGSMass();   // Mass of residual nuclear environment
    // @@ Pr1 condition (individual particle) can be taken out of brackets for all if's
    if(Pr1&&projPDG==22&&projE<PiPrThresh+(M2ShiftVir+projM2)/DiNuclMass) // ==> Gamma+q
    //if(2>3)                                      //@@ TMP:PhotoAbsorbtion by q is closed
    {
      q4Mom=G4LorentzVector(0.,0.,0.,tgMass-envMass);// PhotoInteracts with BoundedCluster
      valQ=curQC;
#ifdef debug
      G4cout<<"G4QE::CrQ:Q="<<q4Mom<<valQ<<"+vg="<<proj4M<<",Env="<<theEnvironment<<G4endl;
#endif
      G4Quasmon* curQuasmon = new G4Quasmon(valQ, q4Mom, proj4M);//Interaction gam+q inside
      theQuasmons.push_back(curQuasmon);  // Insert Quasmon without incid. gamma (del.eq.)
    }
    else if(Pr1&&(std::fabs(projM2-mPi2)<.00001 && projE-mPi<0.1) && projPDG==-211 &&!fake)
    //if(2>3)                                //@@ ***TMP*** PionAbsorbAtRest by q is closed
    {
      q4Mom=proj4M+G4LorentzVector(0.,0.,0.,tgMass-envMass);// PION + BoundCluster
      valQ=EnFlQC+curQC;
#ifdef debug
      if(projE<mPi)G4cout<<"*VirtualPiM*G4QE::CrQ:Ener(pi-)="<<projE<<"<mPi="<<mPi<<G4endl;
      G4cout<<"G4QEnv::CrQ:Q="<<q4Mom<<valQ<<"+pi="<<proj4M<<",E="<<theEnvironment<<G4endl;
#endif
      G4Quasmon* curQuasmon = new G4Quasmon(valQ, q4Mom, -proj4M);//Interact gam+q inside
      theQuasmons.push_back(curQuasmon);  // Insert Quasmon without incid. gamma (del.eq.)
    }
    //else if(Pr1&&projPDG==2212&&G4UniformRand()<.6+.4*std::exp(-envMass/8192))// keepProj
    else if(2>3)                          // free flying projectile is closed
    {
      q4Mom=proj4M;                       // 4M: QUASMON=Projectile
      valQ=EnFlQC;                        // qc: QUASMON=Projectile
      theEnvironment=memEnviron;
#ifdef debug
      G4cout<<"G4QEnv::CreQAll: Q="<<q4Mom<<valQ<<", QEnv="<<theEnvironment<<G4endl;
#endif
      G4Quasmon* curQuasmon = new G4Quasmon(valQ, q4Mom);
      theQuasmons.push_back(curQuasmon);  // Insert Quasmon (even hadron/gamma) (del.eq.)
    }
    else if(Pr1&&projPDG==2212&&G4UniformRand()>15./(proj4M.e()-mProt))//ExcitatedCluster
    //else if(2>3)                          // No excitation of a cluster by projScattering
    {
      G4double prM=mProt;                 // mass of the projectile (M) (for future gener.)
      G4double prM2=mPro2;                // squared mass of the projectile (M^2)
      G4double scM=mProt;                 // mass of the scattered projectile (M')
      G4double scM2=mPro2;                // squared mass of the scatteredProjectile (M'^2)
      G4QContent scQC=projQC;             // QC of the scattered projectile
      G4QContent chQC(0,0,0,0,0,0);       // Change of the Quasmon QC
      if(G4UniformRand()<.5*envN/envBN)   // (u->u,d, d->d)? in future make it universal
      {
        scM=mNeut;                        // Charge exchange reaction
        scM2=mNeu2;
        scQC=neutQC;
        chQC=projQC-scQC;                 // Charge of the created Quasmon must be changed
      }
      G4double tnM=tgMass-envMass;        // minimal mass of the target cluster (m)
      G4double tnM2=tnM*tnM;              // squared mass of the target cluster (m^2)
      G4double dtnM=tnM+tnM;              // doubled mass of the target cluster (2*m)
      G4double prE=proj4M.e();            // enrgy of the projectile (E)
      G4double mu2=tnM2+dtnM*(prE-scM);   // max squared mass of the excited cluster (mu^2)
      //G4double mu=std::sqrt(mu2);         // max mass of the excited cluster (mu)
      G4double B=.00001;                  // (parameter) slope of the diffraction cone
      G4double rmu2=0.;                   // Chosen sqMass of excitedClust (Prototype)
      G4double tmax=0.;                   // max -t for scattering (Prototype)
      G4double rt=0.;                     // Chosen -t (Prototype)
      G4double om=0.;                     // Energy of the excited cluster (Prototype)
      G4double ep=0.;                     // Energy of the scattered projectile (eps Proto)
      if (prE<prM)G4cout<<"-Warn-G4QEnv::CreQAll:(scat w ex)E="<<prE<<" < M="<<prM<<G4endl;
      G4double Pi=std::sqrt(prE*prE-prM2); // Proj momentum (P)
      G4double po=0.;                      // Scat momentum (p) (Prototype)
      G4double cost=2.;                    // cos(theta) for the scattered proj (Prototype)
      G4int cct=0;                         // Counter of cost attempts (@@ can be limited)
      while ( std::fabs(cost) > 1. )
      {
        cct++;
#ifdef debug
        G4cout<<"-Warning-G4QEnv::CreQAll: c="<<cct<<" (scat w ex) cost="<<cost<<G4endl;
#endif
        rmu2=tnM2*pow(mu2/tnM2,G4UniformRand()); // Chosen SqMass of excitedClust (MMA)
        tmax=mu2-rmu2;                    // max -t for scattering
        rt=-std::log(1.-G4UniformRand()*(1.-std::exp(-B*tmax)))/B; // Chosem -t
        om=(tnM2+rmu2+rt)/dtnM;           // Energy of the excited cluster
        ep=prE+tnM-om;                    // Energy of the scattered projectile (epsilon)
#ifdef debug
        G4cout<<"G4QEnv::CreQAll: m2="<<tnM2<<" < mu2="<<rmu2<<" < "<<mu2<<"=Max2"<<G4endl;
        G4cout<<"G4QEnv::CreQAll: -t="<<rt<<" < "<<tmax<<"=tmax"<<G4endl;
        G4cout<<"G4QEnv::CreQAl: om="<<om<<" > m="<<tnM<<", ep="<<ep<<" > M="<<prM<<G4endl;
#endif
        if(ep<scM)G4cout<<"+Warn-G4QEnv::CreQAl:(scat w ex)Eo="<<prE<<" < M="<<prM<<G4endl;
        po=std::sqrt(ep*ep-scM2);   // Scat momentum (p)
        cost=(prE*ep-0.5*(rt+prM2+scM2))/Pi/po; // cos(theta) for the scattered
      }
      G4double om2=om*om;
      if(om2<rmu2)G4cout<<"-Warn-G4QEnv::CreQA:(scat w ex)e2="<<om<<" < mu2="<<tnM<<G4endl;
#ifdef debug
      G4cout<<"G4QEnv::CreQAll: ct="<<cost<<",pio="<<Pi*po<<",()="<<cost*Pi*po<<G4endl;
      G4double ps=std::sqrt(om2-rmu2);    // Momentum of the excited cluster (p)
#endif
      G4double pfc=po*cost;               // Longitudinal projection of the scattered proj
      G4double pfs=po*std::sqrt(1.-cost*cost); // Transversal projection of scattered proj
      // ---- @@ From here can be a MF for QF nucleon extraction (if used by others)
      G4ThreeVector vdir = proj4M.vect(); // 3-Vector in the projectile direction
      G4ThreeVector vx(0.,0.,1.);         // ProtoOrt in the direction of the projectile
      G4ThreeVector vy(0.,1.,0.);         // First ProtoOrt orthogonal to the direction
      G4ThreeVector vz(1.,0.,0.);         // Second ProtoOrt orthoganal to the direction
      if(vdir.mag2() > 0.)                // the projectile isn't at rest
      {
        vx = vdir.unit();                 // Ort in the direction of the projectile
        G4ThreeVector vv= vx.orthogonal();// Not normed orthogonal vector (!)
        vy = vv.unit();                   // First ort orthogonal to the proj. direction
        vz = vx.cross(vy);                // Second ort orthoganal to the proj. direction
      }
      // ---- @@ End of possible MF (Similar is in G4QCollision)
      G4double phi=twopi*G4UniformRand(); // Phi of the Fermi-Mom
      G4ThreeVector fp=pfc*vx+pfs*(std::sin(phi)*vy+std::cos(phi)*vz);
      G4LorentzVector s4M(fp,ep);
#ifdef debug
      G4cout<<"G4QEnv::CreQA:ps="<<po<<"="<<fp.mag()<<",sM="<<prM<<"="<<s4M.m()<<G4endl;
      G4cout<<"G4QEnv::CreQA:Ee="<<prE*ep<<" =? "<<(prM2+rt/2-Pi*po*cost)<<G4endl;
#endif
      if(std::fabs(s4M.m()-scM)>.001)G4cout<<"-W-G4QE::CQA:M="<<prM<<"#"<<s4M.m()<<G4endl;
      G4LorentzVector c4M=proj4M+G4LorentzVector(0.,0.,0.,tnM)-s4M;
#ifdef debug
      G4cout<<"G4QEnv::CreQA: ec="<<om<<" = "<<c4M.e()<<", pc="<<ps<<" = "
            <<c4M.rho()<<", mc2="<<rmu2<<" = "<<c4M.m2()<<G4endl;
      G4cout<<"G4QEnv::CQA:ht="<<(tnM2+rmu2)/2-tnM*om<<"="<<prM2-prE*ep+Pi*po*cost<<G4endl;
#endif
      if(std::fabs(c4M.m2()-rmu2)>.1)
        G4cout<<"-W-G4QE::CrQ:m2="<<rmu2<<"#"<<c4M.m2()<<",P="<<proj4M<<",M="<<tnM<<G4endl;
      G4QHadron* projH = new G4QHadron(scQC,s4M); // Create scattered Projectile Hadron
      theQHadrons.push_back(projH);       // Fill it to the output vector
      G4Quasmon* curQuasmon = new G4Quasmon(curQC+chQC, c4M); // Q=ExcitedCluster
      theQuasmons.push_back(curQuasmon); // Insert Quasmon (even hadron/gamma) (del.eq.)
    }
    else
    {
      q4Mom=proj4M+G4LorentzVector(0.,0.,0.,tgMass-envMass); // Projectile + BoundCluster
      valQ=EnFlQC+curQC;
#ifdef debug
      G4cout<<"G4QEnv::CreQAll: Q="<<q4Mom<<valQ<<", QEnv="<<theEnvironment<<G4endl;
#endif
      G4Quasmon* curQuasmon = new G4Quasmon(valQ, q4Mom);
      theQuasmons.push_back(curQuasmon); // Insert Quasmon (even hadron/gamma) (del.eq.)
    }
  }
  else
  {
    G4cout<<"---Warning---G4QEnvironment::CreateQuasmon:Strange targPDG="<<targPDG<<G4endl;
    //throw G4QException("***G4QEnvironment::CreateQuasmon: Impossible target");
  }
}

// Calculate a probability to interact with clusters for the givven PDG of the projectile
void G4QEnvironment::PrepareInteractionProbabilities(const G4QContent& projQC, G4double AP)
{
  G4double sum    = 0.;                            // Sum of probabilities of interaction
  G4double probab = 0.;                            // Interaction probability
  G4double denseB = 0.;                            // A#of*prob baryons in dense part
  G4double allB   = 0.;                            // A#of*prob baryons in the nucleus
  G4int pPDG = projQC.GetSPDGCode();               // PDG code of the projectile particle
  if(2>3) {allB=AP; allB=pPDG;}                    // A trick to use not used AP(3M) & pPDG
  for (unsigned index=0; index<theQCandidates.size(); index++)
  {
    G4QCandidate* curCand=theQCandidates[index];   // Intermediate pointer
    G4int cPDG = curCand->GetPDGCode();            // PDG Code of the Candidate
    G4int cST  = curCand->GetStrangeness();        // Strangeness of the candidate
    G4int cBN  = curCand->GetBaryonNumber();       // Baryon Number of the candidate
    G4int cCH  = curCand->GetCharge();             // Charge of the candidate
#ifdef sdebug
    G4cout<<"G4QE::PIP:=---=> #"<<index<<", cPDG="<<cPDG<<",S="<<cST<<G4endl;
#endif
    if(cPDG>80000000&&cPDG!=90000000&&!cST&&cCH>0&&cBN>0&&cCH<=cBN) // ===> Nuclear cluster
    {
      G4int zc = cCH;                              // "Z" of the cluster
      G4int nc = cBN-cCH;                          // "N" of the cluster
      G4double nOfCl=curCand->GetPreProbability(); // A number of clusters of the type
      G4double dOfCl=curCand->GetDenseProbability();// A number of clusters in dense region
#ifdef sdebug
      G4cout<<"G4QE::PIP:Z="<<zc<<",N="<<nc<<",nC="<<nOfCl<<",dC="<<dOfCl<<G4endl;
#endif
      if(cPDG==91000000||cPDG==90001000||cPDG==90000001)
      {
        allB+=nOfCl;
        denseB+=dOfCl;
      }
      G4QContent pQC=curCand->GetQC();             // Quark Content of the candidate
      ////////////G4int pC   = projQC.GetCharge();   // Charge of the projectile
      G4QContent qQC=pQC+projQC;                   // Total Quark content of the Compound
      G4QPDGCode qQPDG(qQC);
      G4int qC   = qQPDG.GetQCode();
      G4double d = abs(zc-nc);
      G4double fact=1./pow(2.,d);
      if (qC<-1) probab=0.;     
      //else if(pPDG==-211&&AP<152.&&cBN<2) probab=0.; // PionCaptureByCluster
      //else if(pPDG==22&&AP<152. || pPDG>90000000)
      else if(pPDG==22 && AP<152.)
      {
        if(cBN<2) probab=nOfCl*cBN*fact; //Gamma Under Pi Threshold (QuarkCapture)
        else probab=0.;
      }
      else if(pPDG==2212)
      {
        //if(cBN<2)probab=nOfCl*cBN*fact; // Moving nucleons hits only nucleons
        //else probab=0.;
        probab=nOfCl*cBN*fact; // Moving nucleons hits composed clusters
        //probab=nOfCl*fact; // Moving nucleons hits compact clusters
      }
      ////////////////////////else if((pPDG==-211&&AP<10.)&&cBN<2) probab=0;//PiCapAtRst(D)
      //else if((pPDG==-211||pPDG==-13)&&AP<27.)probab=dOfCl*cBN*fact;//Pi/Mu-CaptureAtRest
      //else if(pPDG==-211&&AP<10.)            probab=nOfCl*fact;// special PiCaptureAtRest
      //else if(pPDG==-211&&AP<10.)            probab=nOfCl*cBN*(cBN-1)*fact;
      //else                                   probab=nOfCl*fact;
      else                                   probab=nOfCl*cBN*fact;
      //else                                      probab=dOfCl*cBN*fact;
      //if(cBN>1) probab=0.;                       // Suppress clusters
      //if(cBN>2) probab=0.;                       // Suppress heavy clusters
#ifdef sdebug
      G4int pPDG      = projQC.GetSPDGCode();      // PDG code of the projectile particle
      G4int rPDG = qQC.GetSPDGCode();
      G4double baryn = qQC.GetBaryonNumber();
      G4double charge= qQC.GetCharge();
      G4double dq= abs(baryn-charge-charge);
      G4cout<<"G4QE::PIP:P="<<probab<<",ac="<<cBN<<",dq="<<dq<<",f="<<fact<<",qC="
            <<qC<<",rPDG="<<rPDG<<",pPDG="<<pPDG<<",nCP="<<nOfCl<<",dCP="<<dOfCl<<G4endl;
#endif
    }
    else probab=0.;
    sum+=probab;
    curCand->SetIntegProbability(sum);
  }
  if(allB>0.)f2all=(allB-denseB)/allB;
  else       f2all=0.;
} // End of PrepareInteractionProbabilities

//Initialize a Clusters Vector for the Nucleus of the QEnvironment
void G4QEnvironment::InitClustersVector(G4int maxClust, G4int maxA)
{
#ifdef debug
  G4cout<<"G4QEnvironment::InitClustersVector called with nC="<<maxClust<<G4endl;
#endif
  if(maxClust>=0) for (G4int i=0; i<maxClust; i++) 
  {
    //G4int clustQCode = i+90; // Q-code of the cluster in the CHIPS World "IsoNuclei"
    G4int clustQCode = i+53; // @@ cluster Q-code in CHIPS World (53=nQHM in G4QPDGCode.hh)
#ifdef sdebug
    G4cout<<"G4QEnvironment::InitClustersVector: Before Init Q ="<<clustQCode<<G4endl;
#endif
    G4QPDGCode clustQPDG(true,clustQCode);
    //clustQPDG.InitByQCode(clustQCode);
    G4int clusterPDG=clustQPDG.GetPDGCode();
    G4int clustB=clustQPDG.GetBaryNum();
#ifdef sdebug
    G4cout<<"G4QEnvironment::InitClustersVector: Before insert ="<<clusterPDG<<G4endl;
#endif
    //theQCandidates.push_back(new G4QCandidate(clusterPDG)); // (delete equivalent)
    if(clustB<=maxA) theQCandidates.push_back(new G4QCandidate(clusterPDG)); // (del.eq.)
#ifdef sdebug
    G4cout<<"G4QEnvironment::InitClustersVector: Cluster # "<<i<<" with code = "
          <<clusterPDG<<", QC="<<clustQPDG.GetQuarkContent()<<G4endl;
#endif
  }
} // End of InitClastersVector

// Fragmentation of the QEnvironment with MultyQuasmon (the main internal member function)
G4QHadronVector  G4QEnvironment::HadronizeQEnvironment()
{
  static const G4int  NUCPDG = 90000000;
  static const G4QNucleus vacuum(NUCPDG);
  static const G4LorentzVector zeroLV(0.,0.,0.,0.);
  //static const G4QContent zeroQC(0,0,0,0,0,0);
  static const G4QContent PiQC(0,1,0,1,0,0);
  static const G4QContent K0QC(1,0,0,0,0,1);
  static const G4QContent KpQC(0,1,0,0,0,1);
  static const G4QContent SiPQC(0,2,1,0,0,0);
  static const G4QContent SiMQC(2,0,1,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent neutQC(2,1,0,0,0,0);
  //static const G4QContent alQC(6,6,0,0,0,0);
  static const G4QPDGCode nQPDG(2112);
  static const G4QPDGCode pQPDG(2212);
  static const G4QPDGCode lQPDG(3122);
  static const G4QPDGCode s0QPDG(3122);
  static const G4double mPi0 = G4QPDGCode(111).GetMass();
  //static const G4double dPi0 = mPi0+mPi0;
  static const G4double fPi0 = 4*mPi0;
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mK   = G4QPDGCode(321).GetMass();
  static const G4double mK0  = G4QPDGCode(311).GetMass();
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mSigZ= G4QPDGCode(3212).GetMass();
  static const G4double mSigM= G4QPDGCode(3112).GetMass();
  static const G4double mSigP= G4QPDGCode(3222).GetMass();
  //static const G4double mAlph = G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double eps=.003;
  G4int nQuasmons = theQuasmons.size();
#ifdef chdebug
  G4int finCharge=theEnvironment.GetCharge();
  G4int finBaryoN=theEnvironment.GetA();
  G4int nHad=theQHadrons.size();
  if(nHad) for(G4int ih=0; ih<nHad; ih++)
  {
    finCharge+=theQHadrons[ih]->GetCharge();
    finBaryoN+=theQHadrons[ih]->GetBaryonNumber();
  }
  //G4int nQuas=theQuasmons.size();
  if(nQuasmons)for(G4int iq=0; iq<nQuasmons; iq++)
  {
    finCharge+=theQuasmons[iq]->GetCharge();
    finBaryoN+=theQuasmons[iq]->GetBaryonNumber();
  }
  if(finCharge!=totCharge || finBaryoN!=totBaryoN)
  {
    G4cout<<"*::*G4QE::HQ:T(1) tC="<<totCharge<<",C="<<finCharge<<",tB="<<totBaryoN
          <<",B="<<finBaryoN<<",E="<<theEnvironment<<G4endl;
    if(nHad) for(G4int h=0; h<nHad; h++)
    {
      G4QHadron* cH = theQHadrons[h];
      G4cout<<"*::*G4QE::HQ:h#"<<h<<",QC="<<cH->GetQC()<<",PDG="<<cH->GetPDGCode()<<G4endl;
    }
    if(nQuasmons) for(G4int q=0; q<nQuasmons; q++)
    {
      G4Quasmon* cQ = theQuasmons[q];
      G4cout<<"*::*G4QE::HQ:q#"<<q<<",C="<<cQ->GetCharge()<<",QCont="<<cQ->GetQC()<<G4endl;
    }
  }
#endif
#ifdef debug
  G4cout<<"G4QE::HQE:*HADRONIZE Q-ENVIRONMENT="<<theEnvironment<<",nQ="<<nQuasmons<<G4endl;
#endif
  if(nQuasmons<1)                            // "No Quasmons" case -> Fill QEnviron
  {
    G4int nPDG = theEnvironment.GetPDG();    // PDG code of the residual Nucl.Environ.
#ifdef debug
    G4cout<<"G4QE::HQE:***NO QUASMONS***Env="<<nPDG<<theEnvironment.Get4Momentum()<<G4endl;
#endif
    if(nPDG==90000000) return theQHadrons;
    if(nPDG>80000000)
    {
      G4QHadron* rNucleus = new G4QHadron(theEnvironment); // Create HadronEnvironment
      theQHadrons.push_back(rNucleus);       // Fill GS - no further decay (del. equiv.)
#ifdef fdebug
      G4cout<<"G4QEnv::HadrQE: ---->> Fill Environment="<<theEnvironment<<G4endl;
#endif
    }
    return theQHadrons;
  }
  if(theEnvironment.GetPDG()==NUCPDG)        // ==> "Environment is Vacuum" case
  {
#ifdef rdebug
    G4cout<<"G4QEnv::HadrQE: ***Vacuum*** #ofQ="<<nQuasmons<<G4endl;
    G4int totInC=0;
    G4LorentzVector totIn4M(0.,0.,0.,0.);
    for (G4int is=0; is<nQuasmons; is++)     // Sum4mom's of Quasmons for the comparison
    {
      G4Quasmon*      pQ = theQuasmons[is];
      G4LorentzVector Q4M= pQ->Get4Momentum();
      totIn4M           += Q4M;
      totInC            += pQ->GetQC().GetCharge();
    } // End of TotInitial4Momentum summation LOOP over Quasmons
    G4int nsHadr  = theQHadrons.size();      // Update the value of OUTPUT entries
    if(nsHadr) for(G4int jso=0; jso<nsHadr; jso++)// LOOP over output hadrons 
    {
      G4int hsNF  = theQHadrons[jso]->GetNFragments(); // A#of secondary fragments
      if(!hsNF)                                        // Add only final hadrons
      {
        G4LorentzVector hs4Mom = theQHadrons[jso]->Get4Momentum();
        totIn4M          += hs4Mom;
        totInC           += theQHadrons[jso]->GetCharge();
      }
    }
#endif
    G4QNucleus vE(90000000);
    G4int     nlq = 0;                       // Prototype of a#of Living Quasmons
    if(nQuasmons) for(G4int lq=0; lq<nQuasmons; lq++)if(theQuasmons[lq]->GetStatus())nlq++;
    if(nQuasmons) for(G4int iq=0; iq<nQuasmons; iq++)
    {
#ifdef chdebug
      G4int f1Charge=theEnvironment.GetCharge();
      G4int f1BaryoN=theEnvironment.GetA();
      G4int nHad=theQHadrons.size();
      if(nHad) for(G4int ih=0; ih<nHad; ih++)
      {
        f1Charge+=theQHadrons[ih]->GetCharge();
        f1BaryoN+=theQHadrons[ih]->GetBaryonNumber();
      }
      G4int nQuas=theQuasmons.size();
      if(nQuas)for(G4int iqs=0; iqs<nQuas; iqs++)
      {
        f1Charge+=theQuasmons[iqs]->GetCharge();
        f1BaryoN+=theQuasmons[iqs]->GetBaryonNumber();
      }
      if(f1Charge!=totCharge || f1BaryoN!=totBaryoN)
      {
        G4cout<<"*::*G4QE::HQ:(2)q#"<<iq<<",tC="<<totCharge<<",C="<<f1Charge<<",tB="
              <<totBaryoN<<",B="<<f1BaryoN<<",E="<<theEnvironment<<G4endl;
        if(nHad) for(G4int h=0; h<nHad; h++)
        {
          G4QHadron* cH = theQHadrons[h];
          G4cout<<"*:*G4QE::HQ:#"<<h<<",QC="<<cH->GetQC()<<",P="<<cH->GetPDGCode()<<G4endl;
        }
        if(nQuas) for(G4int q=0; q<nQuas; q++)
        {
          G4Quasmon* cQ = theQuasmons[q];
          G4cout<<"*:*G4QE::HQ:q#"<<q<<",C="<<cQ->GetCharge()<<",QC="<<cQ->GetQC()<<G4endl;
        }
      }
#endif
      G4int ist=theQuasmons[iq]->GetStatus();// Status of the Quasmon before fragmentation
      if(ist)
      {
        G4QHadronVector* output=theQuasmons[iq]->Fragment(vE,1);//!!!DESTROY!!! <---------+
        G4int ast=theQuasmons[iq]->GetStatus();  // Quasmon's Status after fragmentation  ^
        if(!ast) nlq--;                          // Reduce nlq if Quasmon decayed         ^
        G4int nHadrons = output->size();         // A#of output Hadrons in the Quasmon    ^
#ifdef debug
        G4cout<<"G4QEnv::HadrQE: ***Vacuum*** Q#"<<iq<<", nHadr="<<nHadrons<<G4endl; //   ^
#endif
        if(nHadrons>0)                           // Copy QHadrons-Quasmon to Output       ^
        {
          for (G4int ih=0; ih<nHadrons; ih++)    // LOOP over QHadrons of the Quasmon     ^
          {
            //G4QHadron* curH=new G4QHadron(output->operator[](ih));// (Del 7 lines below)^
            G4QHadron* curH = new G4QHadron((*output)[ih]); // (Deleted 7 lines below)    ^
#ifdef debug
            G4cout<<"G4QEnv::HadrQE:Vacuum, H#"<<ih<<", QPDG="<<curH->GetQPDG() //        ^
                  <<",4M="<<curH->Get4Momentum()<<G4endl; //                              ^
#endif
            theQHadrons.push_back(curH);         // Fill hadron-copy (delete equivalent)  ^
          }
        }                                        //                                       ^
        else                                     // => "Quasmon can't decay" case         ^
        {                                        //                                       ^
          G4QContent totQC=theQuasmons[iq]->GetQC();//                                    ^
          G4int     tQBN=totQC.GetBaryonNumber();// Baryon Number of not decayed Quasmon  ^
          G4QNucleus     tqN(totQC);             // Define the quasmon as a nucleus       ^
          G4double   gsM=tqN.GetMZNS();          // GS Mass                               ^
          G4LorentzVector tot4M=theQuasmons[iq]->Get4Momentum();
          G4double totQM=tot4M.m();              // Real Mass of Quasmon                  ^
          if(tQBN>0&&totQM>gsM)                  // => "Try Quasmon evaporation" case     ^
          {                                      //                                       ^
            G4QHadron* nuclQ = new G4QHadron(totQC,tot4M); //                             ^
#ifdef fdebug
            G4cout<<"G4QEnv::HadrQE:Vac,tQC"<<totQC<<",t4M="<<tot4M<<G4endl; //           ^
#endif
            EvaporateResidual(nuclQ);            // Evaporate ResNuc (del.equiv)          ^
            theQuasmons[iq]->KillQuasmon();      // Kill evaporated Quasmon               ^
            nlq--;                               //                                       ^
          }
          else if(iq+1<nQuasmons&&nlq>1)         // => "Try to merge with next" case      ^
          {
            G4int s_value=theQuasmons[iq+1]->GetStatus();//Status of the next Quasmon     ^
            theQuasmons[iq+1]->IncreaseBy(theQuasmons[iq]);// Merge with the next Quasmon ^
            theQuasmons[iq]->KillQuasmon();      // Kill the week Quasmon                 ^
            if(s_value) nlq--;                   // Reduce a number of "living Quasmons"  ^
          }
          else if(iq+1==nQuasmons&&iq&&nlq>1)    // => "Quasmon stack is exhosted" case   ^
          {
            G4int s_value=theQuasmons[0]->GetStatus(); // Status of the first Quasmon     ^
            theQuasmons[0]->IncreaseBy(theQuasmons[iq]);// Merge with the first Quasmon   ^
            theQuasmons[iq]->KillQuasmon();      // Kill the week Quasmon                 ^
            if(s_value) nlq--;                   // Reduce a number of "living Quasmons"  ^
          }
          else                                   // "Have a chance to recover" case       ^
          {                                      //                                       ^
#ifdef debug
            G4cout<<"***G4QE::HQE:"<<iq<<",n="<<nHadrons<<",Tot="<<totQC<<totQM<<G4endl;//^
            for (G4int kq=0; kq<nQuasmons; kq++) // LOOP over Quasmons for DEBUG PRINTING ^
              G4cout<<kq<<",St/QC="<<theQuasmons[kq]->GetStatus()<<theQuasmons[kq] //     ^
                    ->GetQC()<<",M="<<theQuasmons[kq]->Get4Momentum().m()<<G4endl; //     ^
#endif
            G4int nOfOUT = theQHadrons.size();   // Total #of QHadrons at this point      ^
            G4double  dM = totQM-gsM;            // Excitation of the Quasmon             ^
            G4bool corrf = true;                 // False when corrected & needs to quit  ^
            while(nOfOUT && corrf)               // LOOP over all existing QHadrons       ^
            {                                    //                                       ^
              G4QHadron*     theLast = theQHadrons[nOfOUT-1];     //  Remember            ^
              G4LorentzVector last4M = theLast->Get4Momentum();   //  all                 ^
              G4QContent      lastQC = theLast->GetQC();          //  content             ^
              G4int           lastS  = lastQC.GetStrangeness();   //  of    // Only       ^
              G4int           totS   = totQC.GetStrangeness();    //  the   // for        ^
              G4int           nFr    = theLast->GetNFragments();  //  Last  // if()       ^
              G4int           gam    = theLast->GetPDGCode();     //        //            ^
              if(gam!=22&&!nFr&&lastS<0&&lastS+totS<0&&nOfOUT>1)//=> Skip K,gam & decayed ^
              {                                  //                                       ^
                G4QHadron* thePrev = theQHadrons[nOfOUT-2];// Kill Prev & make Last->Prev ^
                theQHadrons.pop_back();          // theLastQHadron is excluded from OUTPUT^
                theQHadrons.pop_back();          // thePrevQHadron is excluded from OUTPUT^
                theQHadrons.push_back(thePrev);  // thePrev becomes theLast as an object  ^
                delete     theLast;              // the Last QHadron is destructed        ^
                theLast = thePrev;               // Update parameters (thePrev*->theLast*)^
                last4M = theLast->Get4Momentum();// 4Mom of the previouse Quasmon         ^
                lastQC = theLast->GetQC();       // Quark Content of the previouse Quasmon^
              }                                  //                                       ^
              else                               // Just Clear and destroy theLast        ^
              {                                  //                                       ^
                theQHadrons.pop_back();          // theLastQHadron is excluded from OUTPUT^
                delete         theLast;          // theLastQHadron is deleated as instance^
              }                                  //                                       ^
              totQC+=lastQC;                     // Update (increase) the total QC        ^
              tot4M+=last4M;                     // Update (increase) the total 4-momentum^
              totQM=tot4M.m();                   // Calculate new real total mass         ^
              G4QNucleus nN(totQC);              // Define the Quasmon as a nucleus       ^
              gsM=nN.GetMZNS();                  // Calculate the new GS Mass             ^
              dM = totQM-gsM;                    // Escitation energy for the Quasmon     ^
              if(dM>0)                           // "Mass of Q is big enough" case        ^
              {                                  //                                       ^
                theQuasmons[iq]->InitQuasmon(totQC,tot4M);// Update the week Quasmon      ^
                G4QHadronVector* curout=theQuasmons[iq]->Fragment(vE,1);//!DESTROY! <---+ ^
                ast=theQuasmons[iq]->GetStatus();// Status of the Quasmon               ^ ^
                if(!ast) nlq--;                  // Reduce nlq if Quasmon decayed       ^ ^
                nHadrons=curout->size();         // A#of outputQHadrons in theDecayedQ  ^ ^
#ifdef debug
                G4cout<<"G4QEnv::HadrQE:VacuumRecoverQ#"<<iq<<",n="<<nHadrons<<G4endl;//^ ^
#endif
                if(nHadrons>0)                   // => "QHadrons from Quasmon to Output"^ ^
                {                                //                                     ^ ^
                  for (G4int ih=0; ih<nHadrons; ih++) // LOOP over Hadrons of theQuasmon^ ^
                  {                              //                                     ^ ^
                    //G4QHadron* curH = new G4QHadron(curout->operator[](ih)); //       ^ ^
                    G4QHadron* curH = new G4QHadron((*curout)[ih]); //                  ^ ^
#ifdef debug
                    G4cout<<"G4QEnv::HadrQE:Recovered, H#"<<ih<<", QPDG=" //            ^ ^
                          <<curH->GetQPDG()<<",4M="<<curH->Get4Momentum()<<G4endl;  //  ^ ^
#endif
                    totQC-=curH->GetQC();        // totQC recalculation                 ^ ^
                    tot4M-=curH->Get4Momentum(); // tot4M recalculation                 ^ ^
                    theQHadrons.push_back(curH); // Fill hadron-copy (delete equivalent)^ ^
                    //delete curout->operator[](ih);//>-Necessary to delete instances*>-^ ^
                    delete (*curout)[ih];        // >-*Necessary to delete instances*>--^ ^
                  } // End of LOOP over Hadrons of the Quasmon                          ^ ^
                  curout->clear();               //                                     ^ ^
                  delete curout;                 //>*Necessary to delete VectPointers*>=^ ^
                  corrf = false;                 // Corrected: go out of the while loop ^ ^
                  //break;                       // @@ ??                               ^ ^
                } // End of check for existing output Hadrons in the Quasmon            ^ ^
                else                             //                                     ^ ^
                {                                //                                     ^ ^
                  for_each(curout->begin(), curout->end(), DeleteQHadron()); // >-------^ ^
                  curout->clear();               //                                     ^ ^
                  delete curout;                 //>*Necessary to delete VectPointers>--^ ^
                }                                //                                       ^
              }                                  //                                       ^
              nOfOUT  = theQHadrons.size();      // Update the value of OUTPUT entries    ^
#ifdef rdebug
              G4int tC=totInC;                   // Vacuum: No ResidualEnvironCharge      ^
              G4LorentzVector t4M=totIn4M;       // Vacuum: No ResidualEnvironment 4-Mom  ^
              for (G4int js=0; js<nQuasmons; js++) // Subtract 4mom's of Quasmons from dif^
              {                                  //                                       ^
                G4Quasmon*      pQ = theQuasmons[js]; //                                  ^
                if(pQ->GetStatus())              // Subtract only if Quasmon is alive     ^
                {                                //                                       ^
                  G4LorentzVector Q4M= pQ->Get4Momentum(); //                             ^
                  t4M               -= Q4M;                //                             ^
                  tC                -= pQ->GetQC().GetCharge(); //                        ^
                }                                //                                       ^
                else G4cout<<"G4QE::HQ:SUM-4-Mom s("<<js<<")="<<pQ->GetStatus()<<G4endl;//^
              } // End of Quasmons4Momentum subtractions                                  ^
              if(nOfOUT) for(G4int jpo=0; jpo<nOfOUT; jpo++)// LOOP over output hadrons   ^
              {                                  //                                       ^
                G4int hsNF  = theQHadrons[jpo]->GetNFragments();//A#of secondary fragments^
                if(!hsNF)                                   // Subtract only final hadrons^
                {
                  G4LorentzVector hs4Mom = theQHadrons[jpo]->Get4Momentum(); //           ^
                  t4M                   -= hs4Mom;                           //           ^
                  tC                    -= theQHadrons[jpo]->GetCharge();    //           ^
                }                                //                                       ^
              }                                  //                                       ^
              G4cout<<"G4QE::HQ:|||Vacuum|||4-MomCHECK|||d4M="<<t4M<<",dC="<<tC<<G4endl;//^
#endif
            }                                    // End of the WHILE LOOP                 ^
            //if(!nOfOUT&&nQuasmons==1)          // TRY TO EVAPORATE THE TOTAL SYSTEM     ^
            if((!nOfOUT&&nQuasmons==1)||theEnvironment.GetPDGCode()==NUCPDG)//EvaporTotal ^
            {                                    //                                       ^
              G4int totS=totQC.GetStrangeness(); //  Total Strangeness                    ^
              //G4int totBN=totQC.GetBaryonNumber();// Total Baryon Number                ^
              G4int totPDG=totQC.GetZNSPDGCode();// Convert QC to PDGCOde for the nucleus ^
              if(totS) totPDG-=totS*999999;      // @@ ??                                 ^
#ifdef fdebug
              G4cout<<"G4QE::HQE: totPDG="<<totPDG<<",totM="<<totQM<<G4endl; //           ^
#endif
              G4QHadron* evH = new G4QHadron(totQC,tot4M); // Create a Hadron-ResidualNucl^
              CleanUp();                         //                                       ^
              EvaporateResidual(evH);            // Evaporate ResNuc (del.equiv)          ^
              for_each(output->begin(), output->end(), DeleteQHadron());// >------------->+
              output->clear();                   //                                       ^
              delete output;                     // >--------------->-------------------->+
              return theQHadrons;                //                                       ^
            }                                    //                                       ^
            else if(!nOfOUT)                     // Still remain not used Quasmons        ^
            {                                    //                                       ^
              G4ExceptionDescription ed;         //                                       ^
              ed <<"Can't decay Quasmon: T="<< tot4M << totQC << ",M=" << totQM        // ^
                 <<" < gsM="<< gsM <<", d="<< dM <<",Env="<< theEnvironment << G4endl; // ^
              G4Exception("G4QEnvironment::HadronizeQEnvironment()", "HAD_CHPS_0000",  // ^
                          FatalException, ed);                                         // ^
            }                                    //                                       ^
          } // End of PANIC treatment                                                     ^
        } // End of trouble handling with Quasmon decay in Vacuum                         ^
        for_each(output->begin(), output->end(), DeleteQHadron());  // >----------------->+
        output->clear();                         //                                       ^
        delete output;                           // >----------------->------------------>+
      } // End of check for the already decayed Quasmon
    } // End of the LOOP over Quasmons
  }
  else                                           // ==> "Nuclear environment" case
  {
#ifdef rdebug
    G4cout<<"G4QEnv::HadrQE:FRAGMENTATION IN NUCLEAR ENVIRONMENT nQ="<<nQuasmons<<G4endl;
    G4int totInC=theEnvironment.GetZ();
    G4LorentzVector totIn4M=theEnvironment.Get4Momentum();
    for (G4int is=0; is<nQuasmons; is++) // Sum 4mom's of Quasmons for comparison
    {
      G4Quasmon*      pQ = theQuasmons[is];
      G4LorentzVector Q4M= pQ->Get4Momentum();
      totIn4M           += Q4M;
      totInC            += pQ->GetQC().GetCharge();
    } // End of TotInitial4Momentum summation LOOP over Quasmons
    G4int nsHadr  = theQHadrons.size();        // Update the value of #of OUTPUT entries
    if(nsHadr) for(G4int jso=0; jso<nsHadr; jso++)// LOOP over output hadrons 
    {
      G4int hsNF  = theQHadrons[jso]->GetNFragments(); // A#of secondary fragments
      if(!hsNF)                                        // Add only final hadrons
      {
        G4LorentzVector hs4Mom = theQHadrons[jso]->Get4Momentum();
        totIn4M          += hs4Mom;
        totInC           += theQHadrons[jso]->GetCharge();
      }
    }
#endif
    // @@ Experimental calculations
    G4QContent totInQC=theEnvironment.GetQCZNS();
    G4LorentzVector totIn4M=theEnvironment.Get4Momentum();
    for (G4int is=0; is<nQuasmons; is++) // Sum 4mom's of Quasmons for comparison
    {
      G4Quasmon*      pQ = theQuasmons[is];
      totIn4M           += pQ->Get4Momentum();
      totInQC           += pQ->GetQC();
    } // End of TotInitial4Momentum/QC summation LOOP over Quasmons
    G4double totMass=totIn4M.m();
    G4QNucleus totN(totInQC);
    G4double totM=totN.GetMZNS();
    G4double excE = totMass-totM;
    // @@ End of experimental calculations
    G4int   envA=theEnvironment.GetA();
    //G4int   c3Max = 27;                   // Big number (and any #0) slowes dow a lot
    //G4int   c3Max = 9;                    // Max#of "no hadrons" steps (reduced below?)
    /////////////////G4int   c3Max = 3;
    G4int   c3Max = 1;
    if(excE > fPi0) c3Max=(G4int)(excE/mPi0); // Try more for big excess
    //G4int   c3Max = 1;
    //if(excE > dPi0) c3Max=(G4int)(excE/mPi0); // Try more for big excess
    //G4int   c3Max = 0;
    //if(excE > mPi0) c3Max=(G4int)(excE/mPi0); // Try more for big excess
    //G4int   c3Max = 0;                    // It closes the force decay of Quasmon at all!
    //
    //G4int   premC = 27;
    //G4int   premC = 3;
    G4int   premC = 1;
    //if(envA>1&&envA<13) premC = 24/envA;
    if(envA>1&&envA<19) premC = 36/envA;
    //if(envA>1&&envA<31) premC = 60/envA;
    //if(envA>1&&envA<61) premC = 120/envA;
    G4int  sumstat= 2;                     // Sum of statuses of all Quasmons
    G4bool force  = false;                 // Prototype of the Force Major Flag
    G4int cbR     =0;                      // Counter of the "Stoped by Coulomb Barrier"
    //
    //G4int cbRM    =0;                      // * MaxCounter of "StopedByCoulombBarrier" *
    //G4int cbRM    =1;                      // MaxCounter of the "StopedByCoulombBarrier"
    //G4int cbRM    =3;                      // MaxCounter of the "StopedByCoulombBarrier"
    //G4int cbRM    =9;                      // MaxCounter of the "Stoped byCoulombBarrier"
    G4int cbRM    = c3Max;                 // MaxCounter of the "StopedByCoulombBarrier"
    G4int totC    = 0;                     // Counter to break the "infinit" loop
    //G4int totCM   = 227;                   // Limit for the "infinit" loop counter
    G4int totCM   = envA;                   // Limit for the "infinit" loop counter
    //G4int totCM   = 27;                    // Limit for this counter
    //G4int nCnMax = 1;                      // MaxCounterOfHadrFolts for shortCutSolutions
    //G4int nCnMax = 3;                      // MaxCounterOfHadrFolts for shortCutSolutions
    //G4int nCnMax = 9;                      // MaxCounterOfHadrFolts for shortCutSolutions
    G4int nCnMax = c3Max;                  // MaxCounterOfHadrFolts for shortCutSolutions
    G4bool first=true;                     // Flag of the first interaction (only NucMedia)
    G4int cAN=0;                           // Counter of the nucleon absorptions
    //G4int mcAN=27;                         // Max for the counter of nucleon absorptions
    ///G4int mcAN=9;                          // Max for the counter of nucleon absorptions
    //G4int mcAN=3;                          // Max for the counter of nucleon absorptions
    G4int mcAN=1;                          // Max for the counter of nucleon absorptions
    G4double proNorm=.001;                 // Rescattering parameter mb^-1 
    //G4double proNorm=.0027;                // Rescattering parameter mb^-1
    while (sumstat||totC<totCM)            // ===***=== The MAIN "FOREVER" LOOP ===***===
    {
#ifdef chdebug
      G4int f2Charge=0;
      G4int f2BaryoN=0;
      if(theEnvironment.GetMass()>0.)
      {
        f2Charge=theEnvironment.GetCharge();
        f2BaryoN=theEnvironment.GetA();
      }
      G4int nHad=theQHadrons.size();
      if(nHad) for(G4int ih=0; ih<nHad; ih++)
      {
        f2Charge+=theQHadrons[ih]->GetCharge();
        f2BaryoN+=theQHadrons[ih]->GetBaryonNumber();
      }
      G4int nQuas=theQuasmons.size();
      if(nQuas)for(G4int iqs=0; iqs<nQuas; iqs++)
      {
        f2Charge+=theQuasmons[iqs]->GetCharge();
        f2BaryoN+=theQuasmons[iqs]->GetBaryonNumber();
      }
      if(f2Charge!=totCharge || f2BaryoN!=totBaryoN)
      {
        G4cout<<"*::*G4QE::HQ:(3)(NucEnv)i#"<<totC<<",tC="<<totCharge<<",C="<<f2Charge
              <<",tB="<<totBaryoN<<",B="<<f2BaryoN<<",E="<<theEnvironment<<G4endl;
        if(nHad) for(G4int h=0; h<nHad; h++)
        {
          G4QHadron* cH = theQHadrons[h];
          G4cout<<"*:*G4QE::HQ:#"<<h<<",QC="<<cH->GetQC()<<",P="<<cH->GetPDGCode()<<G4endl;
        }
        if(nQuas) for(G4int q=0; q<nQuas; q++)
        {
          G4Quasmon* cQ = theQuasmons[q];
          G4cout<<"*:*G4QE::HQ:q#"<<q<<",C="<<cQ->GetCharge()<<",QC="<<cQ->GetQC()<<G4endl;
        }
      }
#endif
      totC++;
      if( nQuasmons==1 && sumstat==3 ) cbR++;    // Counter of dead solutions for nQ=1
      else cbR=0;
      G4QContent envQC=theEnvironment.GetQCZNS();// QuarkCont of current NuclearEnvironment
      G4QContent totQC=envQC;                    // Total QuarkContent in the system
      G4double   envM =theEnvironment.GetMass(); // mass of NuclEnvironment (@@GetMZNS())
      G4double   sumM =envM;                     // Sum of all residual masses in theSystem
      G4LorentzVector env4M=theEnvironment.Get4Momentum();      
      G4LorentzVector tot4M=env4M;               // 4-momentum of the Total System
      sumstat         =0;
      G4int     nCount=0;                        // Counter of notsuccessful fragmentations
      G4int     fCount=0;                        // Counter of successful(notFinal) fragm's
      G4int     eCount=0;                        // Counter of not yet decayed Quasmons
      for (G4int iq=0; iq<nQuasmons; iq++)       // Sum up Quasmons for making a decision
      {
        G4Quasmon*      pQ = theQuasmons[iq];
        G4QContent      QQC= pQ->GetQC();
        totQC             += QQC;
        G4LorentzVector Q4M= pQ->Get4Momentum();
        tot4M             += Q4M;
        G4double        QM = Q4M.m();
        sumM              += QM;
        G4int           Qst= pQ->GetStatus();
        sumstat           += Qst;
#ifdef debug
        G4cout<<"G4QEnv::HadrQE:#"<<iq<<", Qst="<<Qst<<", Q="<<Q4M<<Q4M.m()<<QQC<<", Env="
              <<theEnvironment<<",nQ="<<nQuasmons<<G4endl;
#endif
        if(nQuasmons>1 && iq+1==nQuasmons && !Qst && Q4M==zeroLV)
        {
          theQuasmons.pop_back();                // Exclude the zero-Quasmon
          delete pQ;                             // and delet it
          nQuasmons--;
        }
        if(Qst==1||Qst==3||Qst==4)
        {
          fCount++;                              // Incr. counterSuccessfulFragmentations
          nCount=0;                             // SetCounterNotSuccessfulFragmentations
        }
        if(Qst>0)                  eCount++;     // Incr. counter of existing Quasmons 
      } // End of summation LOOP over Quasmons
      G4int      totS  =totQC.GetStrangeness();  // Total Strangeness of the Total System
      G4int      totBN =totQC.GetBaryonNumber(); // Total Baryon Number of the Total System
      G4int      totPDG=0;                       // Total PDG Code for the Current compound
      totM  =0.;                                 // min(GroundSt)Mass of theResidualSystem
      if(totBN < 2)                              // Solution for the light nucleus  
      {
        totPDG=totQC.GetSPDGCode();              // Min totPDGCode for theCurrentCompound
        if(totPDG) totM=G4QPDGCode(totPDG).GetMass(); // min Mass of the Residual System
        else
        {
          // throw G4QException("G4QEnv::HadrQEnv: Impossible PDG for B=1");
          G4Exception("G4QEnvironment::HadronizeQEnvironment()", "HAD_CHPS_0001",
                      FatalException, "Impossible PDG for B=1");
        }
      }
      else
      {
        G4QNucleus totN_temporary(totQC,tot4M);  // Excited nucleus for theResidualSystem
        totN=totN_temporary;
        totM=totN.GetMZNS();                     // min(GroundSt)Mass of theResidualSystem
        totPDG=totN.GetPDG();                    // Total PDG Code for the Current compound
      }
#ifdef fdebug
      G4cout<<"G4QEnv::HadrQE:totC="<<totC<<"<totCM="<<totCM<<",ss="<<sumstat<<G4endl;
#endif
      if ( totC >= totCM || cbR > cbRM)
      {
        CleanUp();
        G4QHadron* evH = new G4QHadron(totQC,tot4M);// Create a Hadron for ResidualNucl
        EvaporateResidual(evH);             // Try to evaporate residual (del. equiv.)
        return theQHadrons;
      }
      // === Now we should be prepared for evaporation ===
      G4int      totChg=totQC.GetCharge();    // Total Electric Charge of the Total System
#ifdef debug
      if(totPDG==90999999||totPDG==90999000||totPDG==90000999||totPDG==89999001)
      G4cout<<"***G4QEnv::HadrQEnv: Meson (1) PDG="<<totPDG<<", M="<<tot4M.m()<<G4endl;
      G4int           nOH=theQHadrons.size(); // A#of output hadrons
      G4LorentzVector s4M=tot4M;              // Total 4-momentum (@@ only for checking)
      if(nOH) for(G4int ih=0; ih<nOH; ih++) s4M+=theQHadrons[ih]->Get4Momentum();     
      G4cout<<"G4QEnv::HadrQE:tBN="<<totBN<<",s="<<sumstat<<",fC="<<fCount<<",eC="<<eCount
            <<",En="<<theEnvironment<<",nH="<<nOH<<",tLV="<<s4M<<totQC<<nCount<<G4endl;
#endif
      if(totBN<2)                             // => "Baryons or Mesons" case (@@antiBaryon)
      {
        totPDG =totQC.GetSPDGCode();
        if(totPDG && totPDG!=10 && totPDG!=1114 && totPDG!=2224)
                                                         totM=G4QPDGCode(totPDG).GetMass();
        else if(totPDG==1114) totM=mNeut+mPi;
        else if(totPDG==2224) totM=mProt+mPi;
        else if(totPDG==10)
        {
          G4QChipolino totChip(totQC);        // define the residual as a Chipolino
          totM  =totChip.GetQPDG1().GetMass()+totChip.GetQPDG2().GetMass();
        }
        else
        {
          G4ExceptionDescription ed;
          ed << "Impossible Hadron in CHIPS: totPDG=" << totPDG << ", totQC="
             << totQC << G4endl;
          G4Exception("G4QEnvironment::HadronizeQEnvironment()", "HAD_CHPS_0002",
                      FatalException, ed); 
        }
      }
      totMass = tot4M.m();                        // Total effective Mass
      G4bool   Premium = eCount && premC && envM; // Premium condition
      G4int    count3  = 0;
      if(sumstat && (fCount||Premium) && !force && count3<c3Max)//=>"Force Quasmons decay"
      {
        if(!fCount) premC--;                  // Reduce premium efforts counter
        if(nQuasmons) for (G4int jq=0; jq<nQuasmons; jq++)//FragmentationLOOP over Quasmons
        {
          G4Quasmon* pQ     = theQuasmons[jq];// Pointer to the CurrentQuasmon <--<--<--+
          G4int      status = pQ->GetStatus();// Old status of the Quasmon              ^
#ifdef debug
          G4cout<<"G4QE::HQE:Status of Q#"<<jq<<" (before Fragment)="<<status<<G4endl;//^
#endif
          if(status)                          // Skip dead Quasmons                     ^
          {                                   //                                        ^
            G4int nQuas=eCount;               //                                        ^
            if(nQuas==1&&first) nQuas=-nQuas; //                                        ^
            G4QHadronVector* output=pQ->Fragment(theEnvironment,nQuas);//<DESTRUCT<--<--^-+
#ifdef debug
            G4cout<<"G4QE::HQE:Q#"<<jq<<",*afterFragm* Env="<<theEnvironment<<G4endl;// ^ ^
#endif
            envM =theEnvironment.GetMass();   // new mass of Nuclear Environment        ^ ^
            status = pQ->GetStatus();         // NewStatus after FragmentationAttempt   ^ ^
            if(!status) eCount--;             // Dec. ExistingQuasmonsCounter for Q=0   ^ ^
            G4int nHadrons = output->size();  //                                        ^ ^
#ifdef rdebug
            G4cout<<"G4QE::HQE:**AfterFragmAttempt**#"<<jq<<",stat="<<status<<", Q4M="//^ ^
                  <<pQ->Get4Momentum()<<", Env="<<theEnvironment<<",nH="<<nHadrons //   ^ ^
                  <<",c3="<<count3<<" < "<<c3Max<<",eC="<<eCount<<G4endl;          //   ^ ^
            G4int tC=totInC-theEnvironment.GetZ(); // Subtract theResidualEnvironCharge ^ ^
            G4LorentzVector t4M=totIn4M;      // Compare with the total                 ^ ^
            G4LorentzVector theEnv4m=theEnvironment.Get4Momentum(); // Environment 4Mom ^ ^
            t4M-=theEnv4m;                    // Subtract the Environment 4-momentum    ^ ^
            G4cout<<"G4QEnv::HadrQE:SUM-4Mom e4M="<<theEnv4m<<theEnvironment<<G4endl;// ^ ^
            for (G4int js=0; js<nQuasmons; ++js)// Subtract 4mom's of Quasmons (compare)^ ^
            {                                 //                                        ^ ^
              G4Quasmon*      prQ = theQuasmons[js];//                                  ^ ^
              if(prQ->GetStatus())            // Subtract only if Quasmon is alive      ^ ^
              {                               //                                        ^ ^
                G4LorentzVector Q4M= prQ->Get4Momentum(); //                            ^ ^
                G4QContent qQC= prQ->GetQC(); //                                        ^ ^
                G4cout<<"G4QE::HQE:SUM-4Mom q("<<js<<")4M="<<Q4M<<",QC="<<qQC<<G4endl;//^ ^
                t4M          -= Q4M;          //                                        ^ ^
                tC           -= prQ->GetQC().GetCharge(); //                            ^ ^
              }                               //                                        ^ ^
              else G4cout<<"G4QE::HQE:SUM-4M,st("<<js<<")="<<prQ->GetStatus()<<G4endl;//^ ^
            } // End of Quasmons4Momentum subtractions                                  ^ ^
            G4int nsbHadr=theQHadrons.size(); // Update the value of OUTPUT entries     ^ ^
            if(nsbHadr) for(G4int jpo=0; jpo<nsbHadr; jpo++)// LOOP over output hadrons ^ ^
            {                                 //                                        ^ ^
              G4int hsNF  = theQHadrons[jpo]->GetNFragments();// A#of out fragments     ^ ^
              if(!hsNF)                       // Subtract only final hadrons            ^ ^
              {                               //                                        ^ ^
                G4LorentzVector hs4Mom = theQHadrons[jpo]->Get4Momentum(); //           ^ ^
                G4int hPDG = theQHadrons[jpo]->GetPDGCode(); //                         ^ ^
                G4cout<<"G4QE::HQE:SUM-4-Mom eh("<<jpo<<")4M="<<hs4Mom<<hPDG<<G4endl;// ^ ^
                t4M          -= hs4Mom;       //                                        ^ ^
                tC           -= theQHadrons[jpo]->GetCharge(); //                       ^ ^
              }                               // End of the "FinalHadron" IF            ^ ^
            }                                 // End of the LOOP over output hadrons    ^ ^
            if(nHadrons) for(G4int kpo=0; kpo<nHadrons; ++kpo)//LOOP over out QHadrons  ^ ^
            {                                 //                                        ^ ^
              //G4QHadron* insH =output->operator[](kpo);// Pointer to theOutputQHadron ^ ^
              G4QHadron* insH = (*output)[kpo];// Pointer to the Output QHadron         ^ ^
              G4int qhsNF  = insH->GetNFragments(); // A#of secondary fragments         ^ ^
              if(!qhsNF) // Should never be here // Subtract -> only final hadrons      ^ ^
              {                               //                                        ^ ^
                G4LorentzVector qhs4Mom = insH->Get4Momentum();// 4M of theOutputHadron ^ ^
                G4int hPDG = insH->GetPDGCode(); //  PDG Code of the Q-output Hadron    ^ ^
                G4cout<<"G4QE::HQE:SUM-4-Mom qh("<<kpo<<")4M="<<qhs4Mom<<hPDG<<G4endl;//^ ^
                t4M          -= qhs4Mom;      //                                        ^ ^
                tC           -= insH->GetCharge(); //                                   ^ ^
              }                               //                                        ^ ^
            }                                 // End of the LOOP over output QHadrons   ^ ^
            G4cout<<"G4QEnv::HadrQE:|||||4-MomCHECK||||d4M="<<t4M<<",dC="<<tC<<G4endl;//^ ^
#endif
            if(!status||status==1||nHadrons)  //OutHadronVector was filled in G4Q::Frag ^ ^
            {                                 //                                        ^ ^
              nCount=0;                       // Reset the NotSuccessfulFragmCounter    ^ ^
              if(nHadrons>0)                  // Transfer QHadrons from Quasm to Output ^ ^
              {                               //                                        ^ ^
                for (G4int ih=0; ih<nHadrons; ++ih) // LOOP over Q-output QHadrons      ^ ^
                {                             //                                        ^ ^
                  G4QHadron* inpH = (*output)[ih]; //                                   ^ ^
                  G4int hC=inpH->GetCharge(); // Charge of the Hadron                   ^ ^
                  G4int hF=inpH->GetNFragments();// Number of fragments                 ^ ^
                  G4double hCB=0.;            // Coulomb Barrier prototype              ^ ^
                  G4double hKE=0.;            // Kinetic Energy of the Hadron prototype ^ ^
                  G4LorentzVector hLV=inpH->Get4Momentum(); //                          ^ ^
#ifdef debug
                  G4cout<<"->G4QEnv::HadrQE:H#"<<ih<<", hC="<<hC<<",hF="<<hF<<",4M=" // ^ ^
                        <<hLV<<inpH->GetPDGCode()<<G4endl;  //                          ^ ^
#endif
                  G4bool can = hC && !hF;     // Charged and not yet decayed hadron     ^ ^
                  if(can)                     //                                        ^ ^
                  {                           //                                        ^ ^
                    G4int hB=inpH->GetBaryonNumber(); //                                ^ ^
                    hCB=theEnvironment.CoulombBarrier(hC,hB); // Coulomb barrier        ^ ^
                    hKE=hLV.e()-hLV.m();      // Kinetic energy of the QHadron          ^ ^
                  }                           //                                        ^ ^
                  if(can && hKE < hCB)        // => "Suck Hadron into Quas or Env" case ^ ^
                  {                           //                                        ^ ^
                    if(status)                // => "Suck into existing Quasmon" case   ^ ^
                    {                         //                                        ^ ^
                      G4QContent tQC=inpH->GetQC()+pQ->GetQC();  //                     ^ ^
                      G4LorentzVector tLV=hLV+pQ->Get4Momentum();//                     ^ ^
                      pQ->InitQuasmon(tQC,tLV); // Reinitialize the current Quasmon     ^ ^
#ifdef debug
                      G4cout<<"G4QE::HQE:Medium, H#"<<ih<<", QPDG="<<inpH->GetQPDG() // ^ ^
                            <<",4M="<<inpH->Get4Momentum()<<" is suckedInQ"<<G4endl; // ^ ^
#endif
                    }                         //                                        ^ ^
                    else                      // => "Suck in the Environment" case      ^ ^
                    {                         //                                        ^ ^
                      G4QContent tQC=inpH->GetQC()+theEnvironment.GetQCZNS();//         ^ ^
                      G4LorentzVector tLV=hLV+theEnvironment.Get4Momentum(); //         ^ ^
                      theEnvironment=G4QNucleus(tQC,tLV); // Reinit currentEnvironment  ^ ^
#ifdef debug
                      G4cout<<"G4QE::HQE:Med,H#"<<ih<<",PDG="<<inpH->GetQPDG()<<",4M="//^ ^
                            <<inpH->Get4Momentum()<<" is suckedInEnvironment"<<G4endl;//^ ^
#endif
                    }                         // End of the STATUS IF                   ^ ^
                  }                           // End of the E/Q suck Hdron IF           ^ ^
                  else if(!hF)                // => "Hadron can go out" case, skip dec  ^ ^
                  {                           //                                        ^ ^
                    G4QHadron* curH = new G4QHadron(inpH); //                           ^ ^
#ifdef debug
                    G4LorentzVector ph4M=curH->Get4Momentum(); // 4-mom of the hadron   ^ ^
                    G4double phX=ph4M.x();    // p_x of the hadron                      ^ ^
                    G4double phY=ph4M.y();    // p_y of the hadron                      ^ ^
                    G4double phZ=ph4M.z();    // p_x of the hadron                      ^ ^
                    G4double phCost=phZ/sqrt(phX*phX+phY*phY+phZ*phZ);// Hadr cos(theta)^ ^
                    G4cout<<"G4QEnv::HadrQE:Medium, H#"<<ih<<",QPDG="<<curH->GetQPDG()//^ ^
                          <<", 4M="<<ph4M<<", ct="<<phCost<<G4endl; //                  ^ ^
#endif
                    G4QContent qhdQC=curH->GetQC(); // QuarkCont of current QHadron     ^ ^
                    G4LorentzVector qhd4M=curH->Get4Momentum(); // 4Mom of currQHadron  ^ ^
                    // ... Here the rescattering starts ...                             ^ ^
                    G4int qhdBN=curH->GetBaryonNumber();// Baryon number of the hadron  ^ ^
                    G4bool scat= false;                 // Scattering happened Flag     ^ ^
                    G4int EnvZ = theEnvironment.GetZ(); // Z of the ResidualEnvironment ^ ^
                    G4int EnvN = theEnvironment.GetN(); // N of the ResidualEnvironment ^ ^
                    G4int EnvS = theEnvironment.GetS(); // S of the ResidualEnvironment ^ ^
                    G4int EnvA = EnvZ + EnvN + EnvS;    // A of the ResidualEnvironment ^ ^
                    G4double EnvM = theEnvironment.GetMZNS(); // ResidualEnvironmentMass   ^ ^
                    G4LorentzVector Env4M= theEnvironment.Get4Momentum(); //            ^ ^
#ifdef debug
                    G4cout<<"=*=>G4QE::HQEnv: Env4M="<<Env4M<<",Z="<<EnvZ<<",N="<<EnvN//^ ^
                          <<",S="<<EnvS<<G4endl;                          //            ^ ^
#endif
                    G4int hPDG = curH->GetPDGCode();      //                            ^ ^
                    G4LorentzVector h4M = curH->Get4Momentum(); //                      ^ ^
                    if(EnvA>1 && qhdBN>-1 && qhdBN<2 && h4M.vect().mag() > 0.000001 && hPDG>111 && //^ ^
                       hPDG!=222 && hPDG!=333)//** Quasi-free interaction is possible **^ ^
                    { //                                                                ^ ^
                      --EnvA;                             //                            ^ ^
                      G4double hM2 = h4M.m2();            //                            ^ ^
                      G4double hP  = h4M.rho();           //                            ^ ^
                      G4int pi0F=0;                       // Pi0 flag                   ^ ^
                      if(hPDG==111 || hPDG==222 || hPDG==333) pi0F=hPDG; //             ^ ^
                      if(pi0F) hPDG=211;                  //                            ^ ^
                      pair<G4double,G4double> Xp=theQFScat->FetchElTot(hP,hPDG,false);//^ ^
                      if(pi0F)                            // Get Pi+ and make the mean//^ ^
                      { //                                                              ^ ^
                        hPDG=-211;                        //                            ^ ^
                        pair<G4double,G4double>Y=theQFScat->FetchElTot(hP,hPDG,false);//^ ^
                        G4double fst=(Xp.first+Y.first)/2;//                            ^ ^
                        G4double snd=(Xp.second+Y.second)/2;//                          ^ ^
                        Xp.first  = fst;                  //                            ^ ^
                        Xp.second = snd;                  //                            ^ ^
                        hPDG=211;                         //                            ^ ^
                      } //                                                              ^ ^
                      G4double XSp = Xp.second;           // TotXS hp                   ^ ^
                      pair<G4double,G4double> Xn=theQFScat->FetchElTot(hP,hPDG,true);// ^ ^
                      if(pi0F)                            // Get Pi+ and make the mean//^ ^
                      { //                                                              ^ ^
                        hPDG=-211;                        //                            ^ ^
                        pair<G4double,G4double> Y=theQFScat->FetchElTot(hP,hPDG,true);//^ ^
                        G4double fst=(Xn.first+Y.first)/2;//                            ^ ^
                        G4double snd=(Xn.second+Y.second)/2;//                          ^ ^
                        Xn.first  = fst;                  //                            ^ ^
                        Xn.second = snd;                  //                            ^ ^
                        hPDG=pi0F;                        //                            ^ ^
                        pi0F=0;                           //                            ^ ^
                      } //                                                              ^ ^
                      G4double XSn = Xn.second;           // TotXS hn                   ^ ^
                      G4double XSZ = XSp * EnvZ;          //                            ^ ^
                      G4double XSN = XSn * EnvN;          //                            ^ ^
                      G4double XSA = XSZ + XSN;           //                            ^ ^
                      if(hM2 > 10000. && XSA > 0.)        // One can try to scatter     ^ ^
                      { //                                                              ^ ^
                        G4double Prob=XSA*sqrt(hM2)*G4QThd(EnvA)*proNorm/h4M.e()/EnvA;//^ ^
                        if(G4UniformRand() < Prob)        // Scattering                 ^ ^
                        { //                                                            ^ ^
                          G4double XEp = Xp.first;        // ElXS hp                    ^ ^
                          G4double XEn = Xn.first;        // ElXS hn                    ^ ^
                          G4double XEZ = XEp * EnvZ;      //                            ^ ^
                          G4double XEN = XEn * EnvN;      //                            ^ ^
                          G4double XEA = XEZ + XEN;       //                            ^ ^
                          G4int NPDG=2112;                // Quasi-Elastic on Neutron   ^ ^
                          G4bool flN=true;                //                            ^ ^
                          G4QNucleus newE(1,0);           // Fake Proto                 ^ ^
                          G4LorentzVector N4M(0.,0.,0.,0.); // Fake Proto               ^ ^
                          if(G4UniformRand() < XEA/XSA)   // Quasi-Elastic              ^ ^
                          { //                                                          ^ ^
                            if(G4UniformRand() < XEZ/XEA) // Quasi-Elastic on Proton    ^ ^
                            { //                                                        ^ ^
                              NPDG= 2212;                 //                            ^ ^
                              flN = false;                //                            ^ ^
                            } //                                                        ^ ^
                            if(flN) newE=G4QNucleus(EnvZ, EnvN-1, EnvS); // QF n        ^ ^
                            else    newE=G4QNucleus(EnvZ-1, EnvN, EnvS); // QF p        ^ ^
                            G4double mT=EnvM - newE.GetMZNS(); // Virtual mass          ^ ^
                            N4M = (mT/EnvM)*Env4M;        //                            ^ ^
                            newE.Set4Momentum(Env4M-N4M); //                            ^ ^
#ifdef debug
                            G4cout<<"==>G4QE::HQE:QEl,NPDG=="<<NPDG<<",N4M="<<N4M    // ^ ^
                                  <<",hPDG="<<hPDG<<",h4M="<<h4M<<G4endl;            // ^ ^
#endif
                            pair<G4LorentzVector,G4LorentzVector> RS =  //              ^ ^
                              theQFScat->Scatter(NPDG, N4M, hPDG, h4M); //              ^ ^
#ifdef debug
                            G4cout<<"**>G4QE::HQE:QEl,N4M="<<RS.first<<",h4M=" //       ^ ^
                               <<RS.second<<",d="<<N4M+h4M-RS.first-RS.second<<G4endl;//^ ^
#endif
                            if((RS.first).e() > 0.) // The Elastic Scattering is made //^ ^
                            { //                                                        ^ ^
				      curH->Set4Momentum(RS.second);// justUpdate, NO scat=true ^ ^
                              G4QHadron* qfN = new G4QHadron(NPDG, RS.first); //        ^ ^
                              theQHadrons.push_back(qfN);// Fill QFN(delete equivalent) ^ ^
                              theEnvironment=newE;       // *** change Environment ***  ^ ^
#ifdef debug
                              G4cout<<"*>G4QE::HQE:QE,PDG="<<NPDG<<",4M=" //            ^ ^
                                    <<RS.first<<",***newEnv***: "<<newE<<G4endl;//      ^ ^
#endif
                              // The totQC & tot4M must be reduced by the qfN           ^ ^
                              tot4M-=RS.first;            // current tot4M is reduced   ^ ^
                              if     (NPDG==2212) totQC-=protQC; // sub n               ^ ^
                              else if(NPDG==2112) totQC-=neutQC; // sub p               ^ ^
                              else G4cout<<"*W*>G4QE::HQE:QE,Bad PDG="<<NPDG<<G4endl;// ^ ^
                            } //                                                        ^ ^
                          } //                                                          ^ ^
                          else                            // Quasi-Inelastic            ^ ^
                          { //                                                          ^ ^
                            if(G4UniformRand() < (XSZ-XEZ)/(XSA-XEA))// QInEl on Proton ^ ^
                            { //                                                        ^ ^
                              NPDG= 2212;                 //                            ^ ^
                              flN = false;                //                            ^ ^
                            } //                                                        ^ ^
                            if(flN) newE=G4QNucleus(EnvZ, EnvN-1, EnvS); // QF n        ^ ^
                            else    newE=G4QNucleus(EnvZ-1, EnvN, EnvS); // QF p        ^ ^
                            G4double mT=EnvM - newE.GetMZNS(); // Virtual mass          ^ ^
                            N4M = (mT/EnvM)*Env4M;        //                            ^ ^
                            newE.Set4Momentum(Env4M-N4M); //                            ^ ^
#ifdef debug
                            G4cout<<"==>G4QE::HQE:QInEl,NPDG=="<<NPDG<<",N4M="<<N4M  // ^ ^
                                  <<",hPDG="<<hPDG<<",h4M="<<h4M<<G4endl;            // ^ ^
#endif
                            G4QHadronVector* Q=theQFScat->InElF(NPDG, N4M, hPDG, h4M);//^ ^
                            if(Q) // Inelastic reaction succeeded                       ^ ^
                            { //                                                        ^ ^
                              theQHadrons.push_back((*Q)[0]); // Fill 1st (del. equiv.) ^ ^
                              theQHadrons.push_back((*Q)[1]); // Fill 2nd (del. equiv.) ^ ^
                              theQHadrons.push_back((*Q)[2]); // Fill 3d  (del. equiv.) ^ ^
                              theEnvironment=newE;            // *** change Environ *** ^ ^
#ifdef debug
                              G4cout<<"*>G4QE::HQE:QIE,PDG1="<<(*Q)[0]->GetPDGCode() // ^ ^
                                    <<",4M1="<<(*Q)[0]->Get4Momentum()<<G4endl;      // ^ ^
                              G4cout<<"*>G4QE::HQE:QIE,PDG2="<<(*Q)[1]->GetPDGCode() // ^ ^
                                    <<",4M1="<<(*Q)[1]->Get4Momentum()<<G4endl;      // ^ ^
                              G4cout<<"*>G4QE::HQE:QIE,PDG3="<<(*Q)[2]->GetPDGCode() // ^ ^
                                    <<",4M1="<<(*Q)[2]->Get4Momentum()<<G4endl;      // ^ ^
                              G4cout<<"*>G4QE::HQE:QIE,***NewEnv***: "<<newE<<G4endl;// ^ ^
#endif
                              scat=true;                      // Don't fill the Primary ^ ^
                              delete Q;
                            } //                                                        ^ ^
                          } //                                                          ^ ^
                        } //                                                            ^ ^
                      } //                                                              ^ ^
                    } //                                                                ^ ^
                    if(!scat) //                                                        ^ ^
                    // ... Rescattering END ........................                    ^ ^
                    theQHadrons.push_back(curH); // Fill hadronCopy (delete equivalent) ^ ^
                    totQC-=qhdQC;             //  Update QC of Env + Quasmons           ^ ^
                    tot4M-=qhd4M;             //  Update 4Mom of Env + Quasmons         ^ ^
                  }                           //                                        ^ ^
                }                             // ==> End of the LOOP over outQHadrons   ^ ^
                pQ->ClearOutput();            // Hadrons are filled, Clear Frag-out <-<-^ ^
                count3=0;                     // Reset counter of empty hadronizations    ^
                //c3Max=1;                    // Reduce repetition Max to accelerate      ^
                first=false;                  // First hadronization is for sure is over  ^
              }                               //                                          ^
              else count3++;                  // Increment counter of empty hadronizations^
            }                                 //                                          ^
            else if(status<0||status==2)      // => "PANIC or NOTHING was done" case      ^
            {                                 //                                          ^
#ifdef debug
              G4cout<<"G4QE::HQE:***PANIC***,status="<<status<<",nC="<<nCount<<G4endl; // ^
#endif
              ++nCount;                       //                                          ^
              if(eCount==1 && status<0 && CheckGroundState(pQ,true))// Correct & Finish   ^
              {                               //                                          ^
                for_each(output->begin(), output->end(), DeleteQHadron()); // ->->->----->^
                output->clear();              //                                          ^
                delete output;                // >----------------->--------------------->^
                pQ->KillQuasmon();            // If BackFusion succeeded, kill the Quasmon^
                delete pQ;
                eCount--;                     // Reduce the number of the living Quasmons ^
                return theQHadrons;           //                                          ^
              }                               //                                          ^
              else if(status<0&&nHadrons)     // This is just a confusion in the status...^
              {                               //                                          ^
                G4cerr<<"***G4QEnv::HadrQE: nH="<<nHadrons<<"< status="<<status<<G4endl;//^
                for_each(output->begin(), output->end(), DeleteQHadron()); // -->-->-->---^
                output->clear();              //                                          ^
                delete output;                // >------------------->------------------->^
                G4Exception("G4QEnvironment::HadronizeQEnvironment()",      //            ^
                            "HAD_CHPS_0003", JustWarning, "Do Nothing Er"); //            ^
              }                               //                                          ^
              else if(status==2 && eCount==1 && cAN<mcAN && envM>500.)// Add N from E to Q^
              {                               //                                          ^
#ifdef debug
		    G4cout<<"G4QE::HQE:E="<<theEnvironment<<",M="<<envM<<",c="<<cAN<<G4endl;//^
#endif
                cAN++;                        // Increment the counter of absorptions     ^
                G4int envPDG = theEnvironment.GetPDG(); // PDGCode of the NuclQEnvironment^
                env4M=theEnvironment.Get4Momentum();    // 4mom of the NucEnv             ^
                G4int envN=theEnvironment.GetN();                    // N of Env          ^
                G4int envZ=theEnvironment.GetZ();                    // Z of Env          ^
                G4int resPDG=envPDG-1;        // Residual for the neutron (prototype)     ^
                G4QContent nucQC=neutQC;      // Nucleon Quark Content                    ^
                if ( envN && (envN+envZ)*G4UniformRand() > envZ ) // Change to a proton   ^
                {                             //                                          ^
                  resPDG=envPDG-1000;         // new PDG for the Environment              ^
                  nucQC=protQC;               // proton QContent                          ^
                }                             //                                          ^
#ifdef debug
		    G4cout<<"G4QE::HQE:P,eZ="<<envZ<<",eN="<<envN<<",rPDG="<<resPDG<<G4endl;//^
#endif
                G4QNucleus resNuc(resPDG);    // Create the residual nucleus (future Env) ^
                G4double resM=resNuc.GetGSMass();       // Mass of the residual nucleus   ^
                G4double eM=theEnvironment.GetGSMass(); // Mass of the current environment^
                G4double nucM=eM-resM;                  // Effective mass of theNucleon   ^
                G4LorentzVector res4M(0.,0.,0.,resM);   // Prototype of newEnv4M (at rest)^
                G4LorentzVector nuc4M(0.,0.,0.,nucM);   // Prototype of newEnv4M (at rest)^
                if(std::fabs(env4M.e()-eM) > 0.001)     // the Environment is not at rest ^
                {                                       //                                ^
                  res4M=(resM/eM)*env4M;                // Proportional 4M for residEnv   ^
                  nuc4M=(nucM/eM)*env4M;                // Proportional 4M for effNucleon ^
                }                                       //                                ^
                theEnvironment=G4QNucleus(res4M,resPDG);// Update the Environment         ^
                theQuasmons[0]->IncreaseBy(nucQC,nuc4M);// Update the Only Quasmon        ^
#ifdef debug
		    G4cout<<"G4QE::HQE:P,Q="<<nucQC<<nuc4M<<",env="<<theEnvironment<<G4endl;//^
#endif
              }                                         //                                ^
              else if(status==2&&nCount>nCnMax)// Treat PANIC for stat=2 (NothingWasDone) ^
              {                               //                                          ^
#ifdef debug
		    G4cout<<"G4QE::HQE:PANIC,nC="<<nCount<<">"<<nCnMax<<G4endl; //            ^
#endif
                G4QContent qQC=pQ->GetQC();   // QuarkContent of the Quasmon              ^
                G4int      pqC=qQC.GetCharge(); // Charge (nP) of the Current Quasmon     ^
                G4int      pqS=qQC.GetStrangeness(); // Strangeness (nL) of theCurrQuasmon^
                G4int      pqB=qQC.GetBaryonNumber(); // BaryNumber of the CurrentQuasmon ^
                G4LorentzVector cq4M=pQ->Get4Momentum(); // 4Mom of the Current Quasmon   ^
                G4double cqMass=cq4M.m();     // Real Mass of the current Quasmon         ^
                G4double fqMass=G4QPDGCode(22).GetNuclMass(pqC,pqB-pqC-pqS,pqS);//CQ FreeM^
#ifdef edebug
                G4cout<<"G4QEnv::HQE:M="<<cqMass<<">fM="<<fqMass<<",S="<<pqS<<",C="<<pqC//^
                      <<",ePDG="<<theEnvironment.GetPDG()<<",qQC="<<qQC<<",eC="<<eCount //^
                      <<G4endl;                                                     //    ^
#endif
                if(pqB>0&&pqS<0&&cqMass>fqMass)// "AntiStrangeNucleus-Chipolino" case     ^
                {                             //                                          ^
                  G4QHadron* nuclQ = new G4QHadron(qQC,cq4M);// Hadron for AntiStrangeNuc.^
                  theEnvironment.DecayAntiStrange(nuclQ,&theQHadrons);// AntiStrange(D.E.)^
                  pQ->KillQuasmon();          // If BackFusion succeeded, kill theQuasmon ^
#ifdef edebug
                  G4cout<<"G4QEnv::HQE:Status after kill (#"<<jq<<")="<<pQ->GetStatus()// ^
                        <<", nH="<<theQHadrons.size()<<G4endl;                         // ^
#endif
                  tot4M-=cq4M;                // Update TotalResidNucleus for hadronizPro.^
                  totQC-=qQC;                 // Update total QC for the HadronizationPro.^
                  eCount--;                   // Reduce the number of the living Quasmons ^
                }                             //                                          ^
                else if(theEnvironment.GetPDG()!=NUCPDG) // ==> "NuclearEnvironment" case ^
                {                             //                                          ^
                  if(eCount>1)                //                                          ^
                  {                           //                                          ^
#ifdef fdebug
                    G4cout<<"G4QE::HQE:TOTEVAP tPDG="<<totPDG<<",t4M="<<tot4M<<G4endl; // ^
#endif
                    G4QHadron* evH = new G4QHadron(totQC,tot4M); // Create Hadron-ResidNuc^
                    CleanUp();                //                                          ^
                    EvaporateResidual(evH);   // Evaporate ResNuc (delete equivalemt)     ^
                    for_each(output->begin(), output->end(), DeleteQHadron());// >--------^
                    output->clear();          //                                          ^
                    delete output;            // >---------->----------->---------------->^
                    return theQHadrons;       //                                          ^
                  }                           //                                          ^
                  G4LorentzVector t4M=cq4M+theEnvironment.Get4Momentum(); // Q+E tot4Mom  ^
                  G4double      tM=t4M.m();   // Real total (Quasmon+Environment) mass    ^
                  envQC=theEnvironment.GetQCZNS(); // QuarkCont of NucEnviron             ^
                  G4QContent curQC=envQC+qQC; // Total Quark Content                      ^
                  G4QNucleus curE(curQC);     // Pseudo nucleus for the Total System      ^
                  G4double   curM=curE.GetGSMass();// min mass of the Total System        ^
#ifdef edebug
                  G4cout<<"G4QEnv::HQE:Q#"<<jq<<",tM="<<tM<<">gsM="<<curM<<curE<<G4endl;//^
#endif
                  if(tM<curM)                 //                                          ^
                  {
                    G4int qPDG=qQC.GetZNSPDGCode();// PDG Code of the Quasmon             ^
                    G4double qMass=G4QPDGCode(qPDG).GetMass(); // GroundStM of theQuasmon ^
#ifdef edebug
                    G4cout<<"G4QE::HQE:nQ="<<nQuasmons<<",eC="<<eCount<<",qPDG="<<qPDG // ^
                          <<",qM="<<qMass<<",eM="<<envM<<",tM="<<tM<<",Q+E="<<qMass+envM//^
                          <<G4endl;           //                                          ^
#endif
                    if(eCount==1&&qPDG&&qMass&&tM>qMass+envM)//==> Q+E decay for one Quasm^
                    //if(nQuasmons==1 && qPDG && qMass && tM>qMass+envM) // ==> Q+E decay ^
                    {                         //                                          ^
                      G4int envPDG = theEnvironment.GetPDG(); // PDGCode of the NuclQEnv. ^
#ifdef edebug
                      G4cout<<"G4QEnv::HadrQEnv: Q+E decay, nQ=1, qPDG=="<<qPDG<<G4endl;//^
#endif
                      // => "Quasmon-Chipolino or Environment-Dibaryon" case              ^
                      if(qPDG==10 || qPDG==92000000 || qPDG==90002000 || qPDG==90000002)//^
                      {                          //                                       ^
                        G4QPDGCode h1QPDG=nQPDG; // QPDG of the first hadron              ^
                        G4double   h1M   =mNeut; // Mass of the first hadron              ^
                        G4QPDGCode h2QPDG=h1QPDG;// QPDG of the second hadron             ^
                        G4double   h2M   =mNeut; // Mass of the second hadron             ^
                        if(qPDG==10)             // CHIPOLINO decay case                  ^
                        {                     //                                          ^
                          G4QChipolino QChip(qQC);// define the Quasmon as a Chipolino    ^
                          h1QPDG=QChip.GetQPDG1();// QPDG of the first hadron             ^
                          h1M   =h1QPDG.GetMass();// Mass of the first hadron             ^
                          h2QPDG=QChip.GetQPDG2();// QPDG of the second hadron            ^
                          h2M   =h2QPDG.GetMass();// Mass of the second hadron            ^
                        }                     //                                          ^
                        else if(qPDG==90002000) // DiProton decay case                    ^
                        {                     //                                          ^
                          h1QPDG=pQPDG;       // QPDG of the first hadron                 ^
                          h1M   =mProt;       // Mass of the first hadron                 ^
                          h2QPDG=h1QPDG;      // QPDG of the second hadron                ^
                          h2M   =mProt;       // Mass of the second hadron                ^
                        }                     //                                          ^
                        else if(qPDG==92000000) // Two lambdas case                       ^
                        {                     //                                          ^
                          h1QPDG=lQPDG;       // QPDG of the first hadron                 ^
                          h1M   =mLamb;       // Mass of the first hadron                 ^
                          h2QPDG=h1QPDG;      // QPDG of the second hadron                ^
                          h2M   =mLamb;       // Mass of the second hadron                ^
                          G4double ddMass=totMass-envM; // Free CM energy                 ^
                          if(ddMass>mSigZ+mSigZ) // Sigma0+Sigma0 is possible             ^
                          {                   // @@ Only two particles PS is used         ^
                            G4double dd2=ddMass*ddMass; // Squared free energy            ^
                            G4double sma=mLamb+mLamb; // Lambda+Lambda sum                ^
                            G4double pr1=0.;          // Prototype to avoid sqrt(-)       ^
                            if(ddMass>sma) pr1=sqrt((dd2-sma*sma)*dd2); // Lamb+Lamb PS   ^
                            sma=mLamb+mSigZ;          // Lambda+Sigma0 sum                ^
                            G4double smi=mSigZ-mLamb; // Sigma0-Lambda difference         ^
                            G4double pr2=pr1;         // Prototype of +L+S0 PS            ^
                            if(ddMass>sma&&ddMass>smi) //                                 ^
                              pr2+=sqrt((dd2-sma*sma)*(dd2-smi*smi)); //                  ^
                            sma=mSigZ+mSigZ;          // Sigma0+Sigma0 sum                ^
                            G4double pr3=pr2;         // Prototype of +Sigma0+Sigma0 PS   ^
                            if(ddMass>sma) pr3+=sqrt((dd2-sma*sma)*dd2); //               ^
                            G4double hhRND=pr3*G4UniformRand(); // Randomize PS           ^
                            if(hhRND>pr2)     // --> "ENnv+Sigma0+Sigma0" case            ^
                            {                 //                                          ^
                              h1QPDG=s0QPDG;  // QPDG of the first hadron                 ^
                              h1M   =mSigZ;   // Mass of the first hadron                 ^
                              h2QPDG=h1QPDG;  // QPDG of the second hadron                ^
                              h2M   =mSigZ;   // Mass of the second hadron                ^
                            }                 //                                          ^
                            else if(hhRND>pr1)// --> "ENnv+Sigma0+Lambda" case            ^
                            {                 //                                          ^
                              h1QPDG=s0QPDG;  // QPDG of the first hadron                 ^
                              h1M   =mSigZ;   // Mass of the first hadron                 ^
                            }                 //                                          ^
                          }                   //                                          ^
                          else if(ddMass>mSigZ+mLamb) // Lambda+Sigma0 is possible        ^
                          {                   // @@ Only two particles PS is used         ^
                            G4double dd2=ddMass*ddMass; // Squared free energy            ^
                            G4double sma=mLamb+mLamb; // Lambda+Lambda sum                ^
                            G4double pr1=0.;          // Prototype to avoid sqrt(-)       ^
                            if(ddMass>sma) pr1=sqrt((dd2-sma*sma)*dd2); // Lamb+Lamb PS   ^
                            sma=mLamb+mSigZ;          // Lambda+Sigma0 sum                ^
                            G4double smi=mSigZ-mLamb; // Sigma0-Lambda difference         ^
                            G4double pr2=pr1;         //+L+S0 PS                          ^
                            if(ddMass>sma && ddMass>smi) //                               ^
                              pr2+=sqrt((dd2-sma*sma)*(dd2-smi*smi)); //                  ^
                            if(pr2*G4UniformRand()>pr1) // --> "ENnv+Sigma0+Lambda" case  ^
                            {                 //                                          ^
                              h1QPDG=s0QPDG;  // QPDG of the first hadron                 ^
                              h1M   =mSigZ;   // Mass of the first hadron                 ^
                            }                 //                                          ^
                          }                   //                                          ^
                        }                     //                                          ^
                        if(h1M+h2M+envM<totMass) // => "Three particles decay" case       ^
                        {                     //                                          ^
                          G4LorentzVector h14M(0.,0.,0.,h1M);           //                ^
                          G4LorentzVector h24M(0.,0.,0.,h2M);           //                ^
                          G4LorentzVector e4M(0.,0.,0.,envM);           //                ^
                          if(!G4QHadron(tot4M).DecayIn3(h14M,h24M,e4M)) //                ^
                          {                   //                                          ^
                            G4ExceptionDescription ed;                  //                ^
                            ed << "QChip+E DecIn3 error: (0)tM=" << tot4M.m() //          ^
                               << "->h1=" << h1QPDG << "(" << h1M << ")+h2="  //          ^
                               << h1QPDG << "(" << h2M << ")+envM=" << envM   //          ^
                               << "==" << h1M+h2M+envM << G4endl;             //          ^
                            G4Exception("G4QEnvironment::HadronizeQEnvironment()",  //    ^
                                        "HAD_CHPS_0004", FatalException, ed);       //    ^
                          }                   //                                          ^
                          G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M); //    ^
                          theQHadrons.push_back(h1H);        // (delete equivalent)       ^
#ifdef debug
                          G4cout<<"G4QE::HQE:(1) H1="<<h1QPDG<<h14M<<G4endl;        //    ^
#endif
                          G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M); //    ^
                          theQHadrons.push_back(h2H);        // (delete equivalent)       ^
#ifdef debug
                          G4cout<<"G4QE::HQE:(1) H2="<<h2QPDG<<h24M<<G4endl;        //    ^
#endif
                          G4QHadron* qeH = new G4QHadron(envPDG,e4M);               //    ^
                          theQHadrons.push_back(qeH);        // (delete equivalent)       ^
#ifdef debug
                          G4cout<<"G4QE::HQE:(1) QEnv="<<envPDG<<e4M<<G4endl;       //    ^
#endif
                        }                     //                                          ^
                        else                  // Try to recover                           ^
                        {                     //                                          ^
                          //if(eCount==1&&CheckGroundState(pQ,true)) // @@ BackFusion     ^
                          if(eCount==1&&CheckGroundState(pQ))// BackFusion attempt        ^
                          {                   //                                          ^
                            pQ->KillQuasmon();// ??                                       ^
                            eCount--;         // Reduce a#of theLivingQuasmons            ^
                            for_each(output->begin(),output->end(),DeleteQHadron()); //   ^
                            output->clear();  // -->-->-->-->-->-->-->-->-->-->-->-->-->--+
                            delete output;    // >---------------->---------------------->+
                            return theQHadrons; //                                        ^
                          }                   //                                          ^
                          for_each(output->begin(),output->end(),DeleteQHadron()); //-->--+
                          output->clear();    // -->-->-->-->-->-->-->-->-->-->-->-->-->--+
                          delete output;      // >---------------->------------->-------->+
#ifdef fdebug
                          G4cout<<"--Warning--G4QE::HQE:tM="<<tot4M.m()<<"< h1="<<h1QPDG//^
                                <<"(M="<<h1M<<")+h2="<<h1QPDG<<"(M="<<h2M<<")+EM="<<envM//^
                                <<"="<<h1M+h2M+envM<<G4endl; //                           ^
                          //throw G4QException("G4QEnv::HQE:(0)Chi+Env mass > totMass");//^
#endif
                          CleanUp();          //                                          ^
                          G4QHadron* evH = new G4QHadron(totQC,tot4M);// ResidualNuclHadr ^
                          EvaporateResidual(evH);  // Evaporate residual (del. equiv.)    ^
                          return theQHadrons; //                                          ^
                        }                     //                                          ^
                      }                       //                                          ^
                      else                    // => "Two particles decay" case            ^
                      {                       //                                          ^
                        G4LorentzVector fq4M(0.,0.,0.,qMass); //                          ^
                        G4LorentzVector qe4M(0.,0.,0.,envM);  //                          ^
                        if(!G4QHadron(tot4M).RelDecayIn2(fq4M,qe4M,cq4M,1.,1.))//Q ch.dir.^
                        {                     //                                          ^
                          G4cerr<<"***G4QEnv::HadQE:(0)tM="<<tot4M.m()<<"-> qPDG="<<qPDG//^
                                <<"(M="<<qMass<<") + envM="<<envM<<")"<<G4endl;         //^
                          for_each(output->begin(),output->end(),DeleteQHadron()); //-->->+
                          output->clear();    // -->-->-->-->-->-->-->-->-->-->-->-->-->->+
                          delete output;      // >----------------->--------------------->+
                          // throw G4QException("***G4QEnv::HadrQEnv: Q+Env DecIn2 error");//^
                          G4Exception("G4QEnvironment::HadronizeQEnvironment()",
                                      "HAD_CHPS_0005", FatalException,
                                      "Q+Env DecIn2 error");
                        }                     //                                          ^
                        G4QHadron* qH = new G4QHadron(qPDG,fq4M);// the out going Quasmon ^
                        theQHadrons.push_back(qH); // (delete equivalent)                 ^
#ifdef debug
                        G4cout<<"G4QE::HQE:QuasmH="<<qPDG<<fq4M<<G4endl;         //       ^
#endif
                        G4QHadron* qeH = new G4QHadron(envPDG,qe4M);//theRecoilEnvironment^
#ifdef debug
                        G4cout<<"G4QE::HQE:EnvironH="<<envPDG<<qe4M<<G4endl;     //       ^
#endif
                        if(envPDG==92000000||envPDG==90002000||envPDG==90000002) //       ^
                          theEnvironment.DecayDibaryon(qeH,&theQHadrons);        //       ^
                        else theQHadrons.push_back(qeH);// (del.equiv.) *** this too ***  ^
                      }                       //                                          ^
                      for_each(output->begin(),output->end(),DeleteQHadron());//--->-->-->+
                      output->clear();        //                                          ^
                      delete output;          // >--------------->----------------->----->+
                      CleanUp();              // Clean up Environ and Quasmon             ^
                      return theQHadrons;     // Finish the hadronization process         ^
                    }                         //                                          ^
                    else status=-1;           // Q+E && totM below MassShell - PANIC      ^
                  }                           //                                          ^
                  else if(eCount>1&&(nCount>nCnMax||theEnvironment.GetA()<2))// 2Quasmons ^
                  {                                  //                                   ^
                    theEnvironment.InitByPDG(NUCPDG);// KillEnvironment(@@ Q's? CleanUp)) ^
#ifdef fdebug
                    G4cout<<"G4QEnv::HQE:Evaporate Env+Quasm Env="<<curE<<G4endl;//       ^
#endif
                    G4QHadron* nucQE = new G4QHadron(curQC,t4M);// Hadron for Quasm+Envir.^
                    EvaporateResidual(nucQE); // Evaporate residual Quasm+Env(del.equiv.) ^
                    pQ->KillQuasmon();        // If BackFusion succeeded, kill theQuasmon ^
#ifdef edebug
                    G4cout<<"G4QEnv::HQE:StatusAfterKill (#"<<jq<<")="<<pQ->GetStatus()// ^
                          <<", nH="<<theQHadrons.size()<<G4endl;                       // ^
#endif
                    tot4M-=t4M;               // Update TotalResidNucleus for hadronizPro.^
                    totQC-=curQC;             // Update total QC for the HadronizationPro.^
                    eCount--;                 // Reduce the number of the living Quasmons ^
                  }                           //                                          ^
                  if(eCount==1 && tM>=curM)   //==>for one Quasmon evaporate ResTotN      ^
                  {                           //                                          ^
                    theEnvironment.InitByPDG(NUCPDG);// Cancele the Environment           ^
                    G4int ttPDG=totQC.GetSPDGCode(); // Total PDG Code (10 - Chipolino)   ^
#ifdef pcdebug
                    G4cout<<"G4QE::HQE:BefEv 4M="<<tot4M<<",QC="<<totQC<<ttPDG<<G4endl;// ^
#endif
                    for_each(output->begin(),output->end(),DeleteQHadron()); //-->-->-->->+
                    output->clear();          //                                          ^
                    delete output;            // >------------>-------------------------->+
                    CleanUp();                // Clean up the Environ and Quasmons        ^
                    G4int ttBN=totQC.GetBaryonNumber(); //                                ^
                    if(ttPDG==10&&ttBN<2)     // Chipolino case                           ^
                    {                         //                                          ^
                      G4QChipolino QCh(totQC);// define the TotalResidual as a Chipolino  ^
                      G4QPDGCode   h1QPDG=QCh.GetQPDG1();  // QPDG of the first hadron    ^
                      G4double     h1M   =h1QPDG.GetMass();// Mass of the first hadron    ^
                      G4QPDGCode   h2QPDG=QCh.GetQPDG2();  // QPDG of the second hadron   ^
                      G4double     h2M   =h2QPDG.GetMass();// Mass of the second hadron   ^
                      G4double ttM=tot4M.m(); // Mass of the Chipolino                    ^
                      if(h1M+h2M<ttM)         // Two particles decay of Chipolino is pos. ^
                      {                       //                                          ^
                        G4LorentzVector h14M(0.,0.,0.,h1M);       //                      ^
                        G4LorentzVector h24M(0.,0.,0.,h2M);       //                      ^
                        if(!G4QHadron(tot4M).DecayIn2(h14M,h24M)) //                      ^
                        {                                         //                      ^
                          G4ExceptionDescription ed;              //                      ^
                          ed << "QChip (1) DecIn2 error: tM=" << ttM << "->h1="      //   ^
                             << h1QPDG << "(" << h1M << ")+h2=" << h1QPDG            //   ^
                             << "(" << h2M << ")=" << h1M+h2M << G4endl;             //   ^
                          G4Exception("G4QEnvironment::HadronizeQEnvironment()",     //   ^
                                      "HAD_CHPS_0006", FatalException, ed);          //   ^
                        }                     //                                          ^
                        G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);    //   ^
                        theQHadrons.push_back(h1H);               // (delete equivalent)  ^
#ifdef debug
                        G4cout<<"G4QE::HQE: QCip-> H1="<<h1QPDG<<h14M<<G4endl;       //   ^
#endif
                        G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);    //   ^
                        theQHadrons.push_back(h2H);               // (delete equivalent)  ^
#ifdef debug
                        G4cout<<"G4QE::HQE: QChip->H2="<<h2QPDG<<h24M<<G4endl;       //   ^
#endif
                      }                       //                                          ^
                      else                    //                                          ^
                      {                       //                                          ^
                        G4ExceptionDescription ed;                              //        ^
                        ed << " QChip (2) DecIn2 error: tM=" << ttM << totQC << "->h1=" //^
                           << h1QPDG << "(" << h2M << "=" << h1M+h2M << G4endl; //        ^
                        G4Exception("G4QEnvironment::HadronizeQEnvironment()",  //        ^
                                    "HAD_CHPS_0007", FatalException, ed);       //        ^
                      }                       //                                          ^
                    }                         //                                          ^
                    else                      //                                          ^
                    {                         //                                          ^
#ifdef edebug
                      if(ttPDG<80000000&&ttBN<1) //                                       ^
                        G4cout<<"---Warning---G4QE::HQE: NotNuc, tPDG="<<ttPDG<<G4endl;// ^
#endif
                      G4QHadron* evH = new G4QHadron(totQC,tot4M);// Hadron for ResidNucl ^
                      EvaporateResidual(evH); // Evaporate residual (del.equiv.)          ^
                    }                         //                                          ^
                    return theQHadrons;       //                                          ^
                  }                           //                                          ^
                  else if(eCount==1 && CheckGroundState(pQ,true)) // Correct and Finish   ^
                  {                           //                                          ^
                    for_each(output->begin(), output->end(), DeleteQHadron()); // -->-->->+
                    output->clear();          //                                          ^
                    delete output;            // >------------------->---------------->-->+
                    pQ->KillQuasmon();        // If BackFusion succeeded, kill Quasm      ^
                    delete pQ;
                    eCount--;                 // Reduce a#of the living Quasmons          ^
                    return theQHadrons;       //                                          ^
                  }                           //                                          ^
                }                             //                                          ^
                else                          // "Vacuum" case                            ^
                {                             //                                          ^
                  G4QPDGCode QPDGQ=pQ->GetQPDG(); // QPDG Code for the Quasmon            ^
                  G4int PDGQ=QPDGQ.GetPDGCode();  // PDG Code of the QUASMON              ^
#ifdef edebug
                  G4cout<<"G4QEnv::HadrQEnv: vacuum PDGQ="<<PDGQ<<G4endl; //              ^
#endif
                  if(!PDGQ) status=-1;        // Unknown Quasmon in Vaquum - PANIC        ^
                  // @@ There still can be a case for 2pSigma+ or 2nSigma- (PDGCode?)     ^
                  else if(PDGQ==3112||PDGQ==3222||PDGQ==90999001||PDGQ==91000999)// S+/S- ^
                  {                           //                                          ^
#ifdef edebug
                    G4cout<<"G4QEnv::HadrQEnv:Sigma Mass="<<cqMass<<G4endl;      //       ^
#endif
                    G4double hyM=mNeut;       // Prototype of the hyperon mass            ^
                    G4int    hyPDG=2112;      // Prototype of the hyperon PDG Code        ^
                    G4double pigM=mPi;        // Prototype of the gamma/pion mass         ^
                    G4int    pigPDG=-211;     // Prototype of the gamma/pion PDG Code     ^
                    if(PDGQ==3112||PDGQ==90999001) // --> "Sigma-" case                   ^
                    { //                                                                  ^
                      if(cqMass>mPi+mLamb)    // "Lambda + Pi- is possible" case          ^
                      { //                                                                ^
                        hyM   = mLamb;        // Lambda mass                              ^
                        hyPDG = 3122;         // Lambda PDG Code                          ^
                      } //                                                                ^
                      else if(cqMass>mSigM)   // "Sigma- gamma decay" case                ^
                      { //                                                                ^
                        hyM=mSigM;            // Sigma- mass                              ^
                        hyPDG=3112;           // Sigma- PDG Code                          ^
                        pigM=0.;              // Gamma mass                               ^
                        pigPDG=22;            // Gamma PDG Code                           ^
                      } //                                                                ^
                    } //                                                                  ^
                    else if(PDGQ==3222||PDGQ==91000999) // --> "Sigma+" case              ^
                    { //                                                                  ^
                      pigPDG= 211;            // Pi+ PDG Code                             ^
                      if(cqMass>mPi+mLamb)    // --- "Lambda + Pi+ is possible" case      ^
                      { //                                                                ^
                        hyM   = mLamb;        // Lambda mass                              ^
                        hyPDG = 3122;         // Lambda PDG Code                          ^
                        pigM  = mPi;          // Pi+ mass                                 ^
                        pigPDG= 211;          // Pi+ PDG Code                             ^
                      } //                                                                ^
                      else if(cqMass>mSigP)   // "Sigma- gamma decay" case                ^
                      { //                                                                ^
                        hyM=mSigP;            // Sigma+ mass                              ^
                        hyPDG=3222;           // Sigma+ PDG Code                          ^
                        pigM=0.;              // Gamma mass                               ^
                        pigPDG=22;            // Gamma PDG Code                           ^
                      } //                                                                ^
                      else if(cqMass>mPi0+mProt&&G4UniformRand()>.5) // "P + Pi0" case    ^
                      { //                                                                ^
                        hyM   = mProt;        // Proton mass                              ^
                        hyPDG = 2212;         // Proton PDG Code                          ^
                        pigM  = mPi0;         // Pi0 mass                                 ^
                        pigPDG= 111;          // Pi0 PDG Code                             ^
                      } //                                                                ^
                      else if(cqMass<mPi+mNeut)// "P+gamma" case as "N+Pi+" is impossible ^
                      { //                                                                ^
                        hyM   = mProt;        // Proton mass                              ^
                        hyPDG = 2212;         // Proton PDG Code                          ^
                        pigM=0.;              // Gamma mass                               ^
                        pigPDG=22;            // Gamma PDG Code                           ^
                      } //                                                                ^
                      // othing should be done for "N + P+" case                          ^
                    } //                                                                  ^
                    G4LorentzVector b4Mom(0.,0.,0.,hyM); // Hyperon mass                  ^
                    G4LorentzVector m4Mom(0.,0.,0.,pigM);// pion/gamma mass               ^
                    if(!G4QHadron(cq4M).DecayIn2(b4Mom, m4Mom)) // "DecayIn2 failed" case ^
                    { //                                                                  ^
                      G4cout<<"---Warning---G4QE::HQE:H="<<hyPDG<<"(m="<<hyM<<")+G/Pi=" //^
                            <<pigPDG<<"(m="<<pigM<<")="<<hyM+pigM<<">"<<cqMass<<G4endl; //^
                      G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC             ^
                      CleanUp();              //                                          ^
                      if(!CheckGroundState(quasH,true)) // Last posibility to correct     ^
                      { //                                                                ^
                        G4QHadron* hadr = new G4QHadron(totQC,tot4M); //                  ^
                        theQHadrons.push_back(hadr);    // Cor or fill as It Is           ^
#ifdef debug
                        G4cout<<"-Warn-G4QE::HQE:Sig,QC="<<totQC<<",4M="<<tot4M<<G4endl;//^
#endif
                        //throw G4QException("G4QEnvironment::HadronizeQEnv:Sig error");//^
                      } //                                                                ^
                      delete quasH;          // Delete the temporary fake Quasmon         ^
                      return theQHadrons;    //                                           ^
                    } //                                                                  ^
#ifdef debug
                    G4cout<<"G4QEnv::HadronizeQEnv: Sigma="<<PDGQ<<cq4M<<" -> Hyperon="// ^
                          <<hyPDG<<b4Mom<<" + Gamma/Pi="<<pigPDG<<m4Mom<<G4endl; //       ^
#endif
                    G4QHadron* curBar = new G4QHadron(hyPDG,b4Mom); //                    ^
                    theQHadrons.push_back(curBar); // Fill the Hyperon (delete equivalent)^
                    G4QHadron* curMes = new G4QHadron(pigPDG,m4Mom);  //                  ^
                    theQHadrons.push_back(curMes); // Fill the gam/pi (delete equivalent) ^
                    pQ->KillQuasmon();        // Make done the current Quasmon            ^
                    tot4M-=cq4M;              // Update theTotalResidNucl of HadrPr.      ^
                    totQC-=qQC;               // Update total residual QC of HadrPr.      ^
                    eCount--;                 // Reduce a#of the living Quasmons          ^
                  }                           //                                          ^
                  else if(PDGQ==90999002||PDGQ==91001999) // pS+/nS-                      ^
                  {                           //                                          ^
#ifdef edebug
                    G4cout<<"G4QEnv::HadrQEnv: Nucleon+Sigma Mass="<<cqMass<<G4endl; //   ^
#endif
                    G4bool dinFlag = false;   // Di-nucleon flag                          ^
                    G4double hyM=mSigM;       // Prototype of the hyperon mass (n+Sigma-) ^
                    G4int    hyPDG=3112;      // Prototype of the hyperon PDG Code        ^
                    G4double pigM=mNeut;      // Prototype of the nucleon mass            ^
                    G4int    pigPDG=2112;     // Prototype of the nucleon PDG Code        ^
                    if     (PDGQ==90999002)   // --> "n+Sigma-" case                      ^
                    {                         //                                          ^
                      if(cqMass<mNeut+mSigM)  // ----> "DiNeutron+Pi-" case               ^
                      { //                                                                ^
                        dinFlag = true;       // For the final decay                      ^
                        hyM=mNeut+mNeut;      // Di-neutron                               ^
                        hyPDG=2112;           // Neutron PDG Code                         ^
                        pigM=mPi;             // Pi- mass                                 ^
                        pigPDG=-211;          // Pi- PDG Code                             ^
                      } //                                                                ^
                    } //                                                                  ^
                    else if(PDGQ==91001999)   // --> "p+Sigma+" case                      ^
                    { //                                                                  ^
                      hyM=mSigP;              // Sigma+                                   ^
                      hyPDG=3222;             // PDG Code of Sigma+                       ^
                      pigM=mProt;             // Proton mass                              ^
                      pigPDG=2212;            // PDG Code of proton                       ^
                      if(cqMass<mProt+mSigP)  // ----> "Proton+Proton" case               ^
                      { //                                                                ^
                        hyM=mProt;            // Proton mass                              ^
                        hyPDG=2212;           // Proton PDG Code                          ^
                        pigM=mProt;           // Proton mass                              ^
                        pigPDG=2212;          // Proton PDG Code                          ^
                      } //                                                                ^
                    } //                                                                  ^
                    G4LorentzVector b4Mom(0.,0.,0.,hyM); // Hyperon (di-nucleon) mass     ^
                    G4LorentzVector m4Mom(0.,0.,0.,pigM);// Nucleon (pion) mass           ^
                    if(!G4QHadron(cq4M).DecayIn2(b4Mom, m4Mom)) // "DecayIn2 failed" case ^
                    { //                                                                  ^
                      G4cout<<"--Warning--G4QE::HQE:S/D="<<hyPDG<<"(m="<<hyM<<")+N/Pi=" //^
                            <<pigPDG<<"(m="<<pigM<<")="<<hyM+pigM<<">"<<cqMass<<G4endl; //^
                      G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC             ^
                      CleanUp();              //                                          ^
                      if(!CheckGroundState(quasH,true)) // Last posibility to correct     ^
                      { //                                                                ^
                        G4QHadron* hadr = new G4QHadron(totQC,tot4M); //                  ^
                        theQHadrons.push_back(hadr);    // Cor or fill as It Is           ^
#ifdef debug
                        G4cout<<"-Warn-G4QE::HQE:Sig,QC="<<totQC<<",4M="<<tot4M<<G4endl;//^
#endif
                        //throw G4QException("G4QEnvironment::HadronizeQEnv:Sig error");//^
                      } //                                                                ^
                      delete quasH;          // Delete the temporary fake Quasmon         ^
                      return theQHadrons;    //                                           ^
                    } //                                                                  ^
#ifdef debug
                    G4cout<<"G4QEnv::HadronizeQEnv: NSigma="<<PDGQ<<cq4M<<"-> Sigma/dN="//^
                          <<hyPDG<<b4Mom<<" + N/Pi="<<pigPDG<<m4Mom<<G4endl; //           ^
#endif
                    if(dinFlag) b4Mom/=2.;    // Split the 4-mom for the dinucleon        ^
                    G4QHadron* curBar = new G4QHadron(hyPDG,b4Mom); //                    ^
                    theQHadrons.push_back(curBar); // Fill the Hyperon (delete equivalent)^
                    if(dinFlag)               //                                          ^
                    {                         //                                          ^
                      G4QHadron* secBar = new G4QHadron(hyPDG,b4Mom);// Cre. 2-nd nucleon ^
                      theQHadrons.push_back(secBar);// Fill 2-nd nucleon (delete equiv.)  ^
                    }                         //                                          ^
                    G4QHadron* curMes = new G4QHadron(pigPDG,m4Mom);  //                  ^
                    theQHadrons.push_back(curMes); // Fill the gam/pi (delete equivalent) ^
                    pQ->KillQuasmon();        // Make done the current Quasmon            ^
                    tot4M-=cq4M;              // Update theTotalResidNucl of HadrPr.      ^
                    totQC-=qQC;               // Update total residual QC of HadrPr.      ^
                    eCount--;                 // Reduce a#of the living Quasmons          ^
                  }                           //                                          ^
                  else if(PDGQ==90999003||PDGQ==91002999) // ppS+/nnS-                    ^
                  {                           //                                          ^
#ifdef edebug
                    G4cout<<"G4QEnv::HadrQEnv: DiNucleon+Sigma Mass="<<cqMass<<G4endl; // ^
#endif
                    G4bool dinFlag = false;   // Di-nucleon flag                          ^
                    G4double hyM=mSigM;       // Prototype of the hyperon mass (n+Sigma-) ^
                    G4int    hyPDG=3112;      // Prototype of the hyperon PDG Code        ^
                    G4double pigM=mNeut+mNeut;// Prototype of the di-nucleon mass         ^
                    G4int    pigPDG=2112;     // Prototype of the nucleon PDG Code        ^
                    if     (PDGQ==90999003)   // --> "n+Sigma-" case                      ^
                    {                         //                                          ^
                      if(cqMass<pigM+mSigM)   // ----> "DiNeutron+Pi-" case               ^
                      { //                                                                ^
                        dinFlag = true;       // For the final decay                      ^
                        pigM=mNeut+mNeut+mNeut;// Tri-neutron                             ^
                        pigPDG=2112;          // Neutron PDG Code                         ^
                        hyM=mPi;              // Pi- mass                                 ^
                        hyPDG=-211;           // Pi- PDG Code                             ^
                      } //                                                                ^
                    } //                                                                  ^
                    else if(PDGQ==91002999)   // --> "p+Sigma+" case                      ^
                    { //                                                                  ^
                      hyM=mSigP;              // Sigma+                                   ^
                      hyPDG=3222;             // PDG Code of Sigma+                       ^
                      pigM=mProt+mProt;       // Di-Proton mass                           ^
                      pigPDG=2212;            // PDG Code of proton                       ^
                      if(cqMass<pigM+mSigP)   // ----> "DiProton+Pi+" case                ^
                      { //                                                                ^
                        dinFlag = true;       // For the final decay                      ^
                        pigM=mProt+mProt+mProt;// Tri-proton                              ^
                        pigPDG=2212;          // Neutron PDG Code                         ^
                        hyM=mPi;             // Pi+ mass                                  ^
                        hyPDG=211;           // Pi+ PDG Code                              ^
                      } //                                                                ^
                    } //                                                                  ^
                    G4LorentzVector b4Mom(0.,0.,0.,hyM); // Hyperon (di-nucleon) mass     ^
                    G4LorentzVector m4Mom(0.,0.,0.,pigM);// Nucleon (pion) mass           ^
                    if(!G4QHadron(cq4M).DecayIn2(b4Mom, m4Mom)) // "DecayIn2 failed" case ^
                    { //                                                                  ^
                      G4cout<<"--Warning--G4QE::HQE:S/Pi="<<hyPDG<<"(m="<<hyM<<")+D/T=" //^
                            <<pigPDG<<"(m="<<pigM<<")="<<hyM+pigM<<">"<<cqMass<<G4endl; //^
                      G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC             ^
                      CleanUp();              //                                          ^
                      if(!CheckGroundState(quasH,true)) // Last posibility to correct     ^
                      { //                                                                ^
                        G4QHadron* hadr = new G4QHadron(totQC,tot4M); //                  ^
                        theQHadrons.push_back(hadr);    // Cor or fill as It Is           ^
#ifdef debug
                        G4cout<<"-Warn-G4QE::HQE:Sig,QC="<<totQC<<",4M="<<tot4M<<G4endl;//^
#endif
                        //throw G4QException("G4QEnvironment::HadronizeQEnv:Sig error");//^
                      } //                                                                ^
                      delete quasH;          // Delete the temporary fake Quasmon         ^
                      return theQHadrons;    //                                           ^
                    } //                                                                  ^
#ifdef debug
                    G4cout<<"G4QEnv::HadronizeQEnv:2NSigma="<<PDGQ<<cq4M<<"-> Sigma/Pi="//^
                          <<hyPDG<<b4Mom<<" + 2N/3N="<<pigPDG<<m4Mom<<dinFlag<<G4endl; // ^
#endif
                    G4QHadron* curBar = new G4QHadron(hyPDG,b4Mom); //                    ^
                    theQHadrons.push_back(curBar); // Fill the Hyperon (delete equivalent)^
                    if(dinFlag) m4Mom/=3.;   // Split the 4-mom for the dinucleon in 3    ^
                    else m4Mom/=2.;          // Split the 4-mom for the dinucleon in 2    ^
                    G4QHadron* curMes = new G4QHadron(pigPDG,m4Mom);  //                  ^
                    theQHadrons.push_back(curMes); // Fill the gam/pi (delete equivalent) ^
                    G4QHadron* secBar = new G4QHadron(pigPDG,m4Mom);  //                  ^
                    theQHadrons.push_back(secBar); // Fill the gam/pi (delete equivalent) ^
                    if(dinFlag)               //                                          ^
                    {                         //                                          ^
                      G4QHadron* triBar = new G4QHadron(pigPDG,m4Mom);// Cre. 3-d nucleon ^
                      theQHadrons.push_back(triBar);// Fill 3-d nucleon (delete equival.) ^
                    }                         //                                          ^
                    pQ->KillQuasmon();        // Make done the current Quasmon            ^
                    tot4M-=cq4M;              // Update theTotalResidNucl of HadrPr.      ^
                    totQC-=qQC;               // Update total residual QC of HadrPr.      ^
                    eCount--;                 // Reduce a#of the living Quasmons          ^
                  }                           //                                          ^
                  else if (PDGQ!=10)          // @@ Chipolino can wait @@                 ^
                  {                           //                                          ^
                    G4double qM =cq4M.m();    // Real mass of the Quasmon                 ^
                    G4double gsM=QPDGQ.GetMass(); // GSmass of the Quasmon                ^
#ifdef edebug
                    G4cout<<"G4QEnv::HadrQEnv:#"<<jq<<",qM="<<qM<<">gsM="<<gsM<<G4endl;// ^
#endif
                    if(fabs(qM-gsM)<0.0001)   // "Fill & Kill" Case                       ^
                    {                         //                                          ^
                      G4QHadron* resQ = new G4QHadron(PDGQ,cq4M); // GSM hadron for CurQ  ^
#ifdef debug
                      G4cout<<"G4QEnv::HadrQEnv:ResQ="<<PDGQ<<cq4M<<G4endl;         //    ^
#endif
                      theQHadrons.push_back(resQ); // @@ Check Dibarions @@ (del.equiv.)  ^
                      pQ->KillQuasmon();      // Make done the current Quasmon            ^
                      tot4M-=cq4M;            // Update theTotalResidNucl of HadrPr.      ^
                      totQC-=qQC;             // Update total residual QC of HadrPr.      ^
                      eCount--;               // Reduce a#of the living Quasmons          ^
                    }                         //                                          ^
                    else if(eCount==1 && qM<gsM && CheckGroundState(pQ,true))// Cor.& Fin.^
                    {                         //                                          ^
#ifdef edebug
                      G4cout<<"G4QEnv::HadrQEnv:**>>** CGS Correction **>>**"<<G4endl; // ^
#endif
                      for_each(output->begin(), output->end(), DeleteQHadron()); //-->-->-^
                      output->clear();        //                                          ^
                      delete output;          // >--------------->--------------->------->+
                      pQ->KillQuasmon();      // If BackFusion -> kill theQuasmon         ^
                      eCount--;               // Reduce a#of the living Quasmons          ^
                      return theQHadrons;     //                                          ^
                    }                         //                                          ^
                    else if(qM<gsM&&(pQ->GetQC().GetSPDGCode()==1114     //               ^
                                     || pQ->GetQC().GetSPDGCode()==2224) //               ^
                          &&qM>theWorld->GetQParticle(QPDGQ)->MinMassOfFragm())//Del&Kill ^
                    {                         //                                          ^
#ifdef edebug
                      G4cout<<"G4QEnv::HadrQEnv:**||** Copy&Decay **||**"<<G4endl; //     ^
#endif
                      G4QHadronVector* decHV=pQ->DecayQuasmon();//Dec.Quasm & fill decHV=*^
                      CopyAndDeleteHadronVector(decHV);// Copy output to QHadrV of G4Env  ^
                      tot4M-=pQ->Get4Momentum(); // tot4M recalculation                   ^
                      totQC-=pQ->GetQC();     // totQC recalculation                      ^
                      pQ->KillQuasmon();      // Make done the current Quasmon            ^
                      eCount--;               // Reduce a#of the living Quasmons          ^
                    }                         //                                          ^
                  }                           //                                          ^
                }                             // End of the Medium/Vacuum IF              ^
              }                               // End of the status ELSE IF                ^
              else if(status==3) count3++;    //                                          ^
              if(status<0)                    // Panic: Quasmon is below theMassShell     ^
              {                               //                                          ^
                //if(eCount==1 && DecayInEnvQ(pQ)) //                                     ^
                //{                             //                                        ^
                //  for_each(output->begin(), output->end(), DeleteQHadron());//--->-->-->+
                //  output->clear();            //                                        ^
                //  delete output;              // >----------------->---------------->-->+
                //  eCount--;                   // Reduce a#of the living Quasmons        ^
                //  pQ->KillQuasmon();          //                                        ^
                //  return theQHadrons;         //                                        ^
                //}                             //                                        ^
                G4int    ppm=jq;                // Initialized by PANIC Quasmon pointer   ^
                G4int    nRQ=0;                 // Prot. of a#of additionalRealQuasmons   ^
#ifdef edebug
                G4cout<<"G4QEnv::HadrQEnv: ***PANIC*** for jq="<<jq<<G4endl; //           ^
#endif
                G4ThreeVector vp= pQ->Get4Momentum().vect(); // PANICQuasmon momentum     ^
                G4double dpm=1.e+30;            // Big number (dot product of momenta)    ^
                if(nQuasmons>1) for(G4int ir=0; ir<nQuasmons; ir++)// Search for partner  ^
                {                               //                                        ^
                  if(ir!=jq)                    // Skip the current (PANIC) Quasmon itself^
                  {                             //                                        ^
                    G4Quasmon* rQ = theQuasmons[ir]; //                                   ^
                    G4int Qst = rQ->GetStatus();// Status of a Quasmon                    ^
#ifdef edebug
                    G4cout<<"G4QEnv::HadrQEnv: ir="<<ir<<",Qstatus="<<Qst<<G4endl; //     ^
#endif
                    if(Qst>0)                   // Skip the dead Quasmon                  ^
                    {
                      nRQ++;                    // Increment real-Quasmon-counter         ^
                      G4double dp=vp.dot(rQ->Get4Momentum().vect()); //                   ^
                      if(dp<dpm)                // Search for the "moving in thesameDir"  ^
                      {                         //                                        ^
                        ppm=ir;                 // Remember the index of MinProj Quasmon  ^
                        dpm=dp;                 // Remember the value of theMinProjection ^
                      }                         //                                        ^
                    }                           //                                        ^
                  } // End of the Quasmon LOOP 
                } // End of the partner-search-for-the-PANIC-Quasmon LOOP                 ^
                if(nRQ)                         // Merge with theBestPartnerQuasmonCandid ^
                {                               //                                        ^
                  G4Quasmon*      rQ = theQuasmons[ppm];   //                             ^
                  G4QContent      rQC= rQ->GetQC();        //                             ^
                  G4LorentzVector r4M= rQ->Get4Momentum(); //                             ^
                  rQC               += pQ->GetQC();        //                             ^
                  r4M               += pQ->Get4Momentum(); //                             ^
                  rQ->InitQuasmon(rQC, r4M);    // Make new Quasmon                       ^
#ifdef edebug
                  G4cout<<"G4QE::HQE:"<<pQ->GetQC()<<"+"<<rQ->GetQC()<<"="<<rQC<<G4endl;//^
#endif
                  pQ->KillQuasmon();            // Delete old Quasmon                     ^
                  eCount--;                     // Decrement counter of living Quasmons   ^
                }                               //                                        ^
                else // No candidate to resolve PANIC was found                           ^
                {                               //                                        ^
#ifdef edebug
                  G4cout<<"G4QEnv::HQE: No Q-cand. nRQ="<<nRQ<<",eC="<<eCount<<G4endl; // ^
#endif
                  //if(eCount==1 && CheckGroundState(pQ,true)) //  BackFusion attempt     ^
                  if(CheckGroundState(pQ,true)) //  The only Q: BackFusion attempt        ^
                  {
                    for_each(output->begin(), output->end(), DeleteQHadron()); //-->-->-->+
                    output->clear();            //                                        ^
                    delete output;              // >-------------->---------------->----->+
                    pQ->KillQuasmon();          //                                        ^
                    eCount--;                   // Reduce a#of the living Quasmons        ^
                    return theQHadrons;         //                                        ^
                  }                             //                                        ^
#ifdef fdebug
                  G4cout<<"G4QEnv::HadrQEnv:NO PANICsolution,t="<<tot4M<<totQC<<G4endl;// ^
#endif
                  totQC=theEnvironment.GetQC(); //                                        ^
                  tot4M=theEnvironment.Get4Momentum(); //                                 ^
                  if(nQuasmons) for(G4int jr=0; jr<nQuasmons; jr++) // Search for partner ^
                  {                             //                                        ^
                    G4Quasmon* rQ = theQuasmons[jr]; // Pointer to the Quasmon            ^
                    G4int Qst = rQ->GetStatus();// Status of a Quasmon                    ^
                    if(jr==jq)                  //                                        ^
                    {                           //                                        ^
                      totQC+=rQ->GetQC();       // QuarkContent of the Quasmon            ^
                      tot4M+=rQ->Get4Momentum();// QuarkContent of the Quasmon            ^
                    }                           //                                        ^
                    else if(Qst)                // Skip dead Quasmons                     ^
                    {                           //                                        ^
                      totQC+=rQ->GetQC();       // QuarkContent of the Quasmon            ^
                      tot4M+=rQ->Get4Momentum();// QuarkContent of the Quasmon            ^
                    }                           //                                        ^
                  } // End of the "No candidate to resolve PANIC" ELSE                    ^
                  pQ->KillQuasmon();            // Kill the only Quasmon                  ^
                  eCount--;                     // Reduce a#of the living Quasmons        ^
                  CleanUp();                    // Clean up THIS Quasmon and Environment  ^
                  G4QHadron* evH = new G4QHadron(totQC,tot4M); // Create ResidNuclHadron  ^
                  EvaporateResidual(evH);       // Try to evaporate residual (del.equiv.) ^
                  for_each(output->begin(), output->end(), DeleteQHadron()); //-->-->-->->+
                  output->clear();              //                                        ^
                  delete output;                // >------------------>--------------->-->+
                  force=true;                   // Make the force decision                ^
                  break;                        // Out of the fragmentation loop >->+     ^
                }                               //                                  |     ^
              }                                 //                                  |     ^
            }                                   //                                  |     ^
            for_each(output->begin(), output->end(), DeleteQHadron()); // ->-->-->--|---->+
            output->clear();                    //                                  |     ^
            delete output;                      // >----------------->--------------|-->->+
          } // End of skip of the dead Quasmons                                     |
#ifdef debug
          G4cout<<"G4QE::HQE:QStat("<<jq<<"="<<status<<pQ->Get4Momentum()<<G4endl;//|
#endif
        } // End of fragmentation LOOP over Quasmons (jq) <--------<----------<-----+
        cAN=0;
      }
      else if(totMass>totM+.001)                // ==> "Try Evaporate or decay" case
      {
#ifdef edebug
        G4cout<<"G4QEnv::HadrQE: NQ="<<nQuasmons<<",tM="<<totMass<<",tPDG="<<totPDG<<",tB="
              <<totBN<<",GSM="<<totM<<",dM="<<totMass-totM<<",totQC="<<totQC<<G4endl;
#endif
        //if(nQuasmons==1)
        if(2>3)                                 // ** closed, because doesn't make a diff
        {
          G4QContent quasQC=totQC-envQC;        // Total QuarkContent of the Only Quasmon
          G4int resQPDG=quasQC.GetSPDGCode();   // GS mass for the Only Quasmon-hadron
          G4int resQB=quasQC.GetBaryonNumber(); // Baryon number of the Only Quasmon
          G4int resQCh=quasQC.GetCharge();      // Charge of the Only Quasmon
          //G4int resQS=quasQC.GetStrangeness();  // Strangeness of the Only Quasmon
          if((resQPDG==0 || resQPDG==10) && resQB>0) resQPDG=quasQC.GetZNSPDGCode();
          G4double resQM=G4QPDGCode(resQPDG).GetMass();// GS Mass of the Only Quasmon
          G4double qCB=theEnvironment.CoulombBarrier(resQCh,resQB); // CoulombBarrier
          G4double de=totMass-envM-resQM-qCB;
#ifdef debug
          G4cout<<"G4QEnv::HadrQE:NQ==1,tM="<<totMass<<",qM="<<resQM<<",eM="<<envM<<",CB="
                <<qCB<<",dE="<<totMass-envM-resQM-qCB<<G4endl;
#endif
          if(de>0.)                             // Make DecayIn2 conserving Q-direction
          {
            G4LorentzVector fq4M=G4LorentzVector(0.,0.,0.,resQM); // Prot. for outQuasmon
            G4LorentzVector fe4M=env4M;         // Prototype for outEnvironment
            G4LorentzVector dir4M=tot4M-env4M;  // Internall quasmon 4-momentum
            if(!G4QHadron(tot4M).RelDecayIn2(fe4M,fq4M,dir4M,1.,1.))
            {
              G4ExceptionDescription ed;
              ed << "Can't decay Q+E: t4M=" << tot4M << ",d=" << de << G4endl;
              G4Exception("G4QEnvironment::HadronizeQEnvironment()", "HAD_CHPS_0008",
                          FatalException, ed); 
            }
              G4QHadron* hQua = new G4QHadron(resQPDG,fq4M);
              theQHadrons.push_back(hQua);      // Fill the hadron-quasmon (delete equiv.)
              G4int envPDG=theEnvironment.GetPDGCode();
              G4QHadron* hEnv = new G4QHadron(envPDG,fe4M);
              theQHadrons.push_back(hEnv);      // Fill the hadron-environ (delete equiv.)
#ifdef debug
              G4cout<<"G4QEnv::HadrQEnv:fQ="<<resQPDG<<fq4M<<", fE="<<envPDG<<fe4M<<G4endl;
#endif
              return theQHadrons;
          }
        }
#ifdef debug
        G4cout<<"G4QEnv::HadrQE: M="<<totMass<<",PDG="<<totPDG<<",B="<<totBN<<",GSM="<<totM
              <<",dM="<<totMass-totM<<",totQC="<<totQC<<G4endl;
#endif
        if(totBN<2)                             // ==> "Baryon/Meson residual Quasmon" case
        {
          if(totPDG==90999999||totPDG==90999000||totPDG==90000999||totPDG==89999001)//"M"ca
          {
            G4cout<<"---Warning---G4QE::HQE:Meson(2) PDG="<<totPDG<<",M="<<totMass<<G4endl;
          }
          else if(totPDG==1114||totPDG==2224)   //==> "DELTA- or DELTA++" case (?antiDELTA)
          {
            G4double   mBar=mProt;
            G4int      bPDG=2212;
            G4double   mMes=mPi;
            G4int      mPDG=211;
            if(totPDG==1114)                    // "DELTA-" case
            {
              mBar=mNeut;
              bPDG=2112;
              mPDG=-211;
            }
            if(totMass<mBar+mMes)
            {
              G4cout<<"--Warning--G4QE::HQE:tM="<<totMass<<"<GSM+mPi0="<<totM+mPi0<<G4endl;
              G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
              CleanUp();
              if(!CheckGroundState(quasH,true))
              {
                G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef debug
                G4cout<<"***G4QE::HQE:FillAsIs(-4),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                //throw G4QException("G4QEnvironment::HadronizeQEnvironment:(1)DecayQEnv");
              }
              delete quasH;  
              return theQHadrons;
            }
            else
            {
              //G4QHadron* delta = new G4QHadron(totQC,tot4M);
              //delta->SetNFragments(2);           // Put a#of Fragments=2
              //theQHadrons.push_back(delta);      // Fill the residual DELTA (del.Eq.)
              // Instead
              //delete delta;
              //
              G4LorentzVector b4Mom(0.,0.,0.,mBar);
              G4LorentzVector m4Mom(0.,0.,0.,mMes);
              if(!G4QHadron(tot4M).DecayIn2(b4Mom, m4Mom))
              {
                G4cout<<"---Warning---G4QEnv::HadronizeQE:B="<<bPDG<<"(m="<<mBar<<") + M="
                      <<mPDG<<"(m="<<mMes<<")="<<mBar+mMes<<" > mDel="<<totMass<<G4endl;
                G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                CleanUp();
                if(!CheckGroundState(quasH,true))
                {
                  G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                  theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef debug
                  G4cout<<"***G4QE::HQE:FillAsIs(-3),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                  //throw G4QException("G4QEnvironment::HadronizeQEnv:Del->Bar+Mes error");
                }
                delete quasH;  
                return theQHadrons;
              }
#ifdef debug
              G4cout<<"G4QEnv::HadronizeQEnv: DELTA="<<totPDG<<tot4M<<" -> Bar="
                    <<bPDG<<b4Mom<<" + Mes="<<mPDG<<m4Mom<<G4endl;
#endif
              G4QHadron* curBar = new G4QHadron(bPDG,b4Mom);
              theQHadrons.push_back(curBar);     // Fill the baryon (delete equivalent)
#ifdef edebug
              G4cout<<"G4QEnv::HadrQEnv:BaryonH="<<bPDG<<b4Mom<<G4endl;
#endif
              G4QHadron* curMes = new G4QHadron(mPDG,m4Mom);
              theQHadrons.push_back(curMes);     // Fill the meson (delete equivalent)
#ifdef edebug
              G4cout<<"G4QEnv::HadrQEnv:MesonH="<<mPDG<<m4Mom<<G4endl;
#endif
              return theQHadrons;
            }
          }
          else if(totPDG==10)                    // ==> "Chipolino" case
          {
            G4QChipolino resChip(totQC);         // define Residual as Chipolino
            G4QPDGCode h1QPDG=resChip.GetQPDG1();// QPDG of the first hadron
            G4int      h1PDG=h1QPDG.GetPDGCode();// PDG code of the first hadron
            G4double   h1M  =h1QPDG.GetMass();   // Mass of the first hadron
            G4QPDGCode h2QPDG=resChip.GetQPDG2();// QPDG of the second hadron
            G4int      h2PDG=h2QPDG.GetPDGCode();// PDG code of the second hadron
            G4double   h2M  =h2QPDG.GetMass();   // Mass of the second hadron
            G4LorentzVector h14Mom(0.,0.,0.,h1M);
            G4LorentzVector h24Mom(0.,0.,0.,h2M);
            if(!G4QHadron(tot4M).DecayIn2(h14Mom, h24Mom))
            {
              G4cout<<"---Warning---G4QEnv::HadronizeQE:h1="<<h1PDG<<"(m="<<h1M<<") + h2="
                    <<h2PDG<<"(m="<<h2M<<")="<<h1M+h2M<<" > mChipo="<<totMass<<G4endl;
              G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
              CleanUp();
              if(!CheckGroundState(quasH,true))
              {
                G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef debug
                G4cout<<"***G4QE::HQE:FillAsIs(-2),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                //throw G4QException("G4QEnvironment::HadrQEnv: Chipo->1+2 decay failed");
              }
              delete quasH;  
              return theQHadrons;
            }
#ifdef debug
            G4cout<<"G4QEnv::HadronizeQEnv: Chipo="<<tot4M<<" -> h1="
                  <<h1PDG<<h14Mom<<" + Mes="<<h2PDG<<h24Mom<<G4endl;
#endif
            G4QHadron* curH1 = new G4QHadron(h1PDG,h14Mom);
            theQHadrons.push_back(curH1);        // Fill the curH1 (delete equivalent)
#ifdef edebug
            G4cout<<"G4QEnv::HadrQEnv:HadronH="<<h1PDG<<h14Mom<<G4endl;
#endif
            G4QHadron* curH2 = new G4QHadron(h2PDG,h24Mom);
            theQHadrons.push_back(curH2);        // Fill the curH2 (delete equivalent)
#ifdef edebug
            G4cout<<"G4QEnv::HadrQEnv:MesAsHadrPartnerH="<<h2PDG<<h24Mom<<G4endl;
#endif
            return theQHadrons;
          }
          else if(totBN<2&&totPDG&&totMass<totM+mPi0+.001)// ==> "Meson/Baryon+gamma" case
          {
            G4LorentzVector h4Mom(0.,0.,0.,totM);
            G4LorentzVector g4Mom(0.,0.,0.,0.);
            if(!G4QHadron(tot4M).DecayIn2(h4Mom, g4Mom))
            {
              G4cout<<"---Warning---G4QEnv::HadronizeQEnv: h="<<totPDG<<"(m="<<totM
                    <<") + gamma > mTot="<<totMass<<G4endl;
              G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
              CleanUp();
              if(!CheckGroundState(quasH,true))
              {
                G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef debug
                G4cout<<"***G4QE::HQE:FillAsIs(-1),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                //throw G4QException("G4QEnvironment::HadronizeQEnv:Gamma Decay failed");
              }
              delete quasH;  
              return theQHadrons;
            }
#ifdef debug
            G4cout<<"G4QE::HQE:"<<tot4M<<"->h="<<totPDG<<h4Mom<<" + gamma="<<g4Mom<<G4endl;
#endif
            G4QHadron* curG = new G4QHadron(22,g4Mom);
            theQHadrons.push_back(curG);         // Fill the gamma (delete equivalent)
#ifdef edebug
            G4cout<<"G4QEnv::HadrQEnv:PhotonH="<<g4Mom<<G4endl;
#endif
            G4QHadron* curH = new G4QHadron(totPDG,h4Mom);
#ifdef edebug
            G4cout<<"G4QEnv::HadrQEnv:GamPartnerH="<<totPDG<<h4Mom<<G4endl;
#endif
            if(totPDG==92000000||totPDG==90002000||totPDG==90000002)
              theEnvironment.DecayDibaryon(curH,&theQHadrons);
            else theQHadrons.push_back(curH);    // Fill the baryon (delete equivalent)
            return theQHadrons;
          }
          else if(totBN<2&&totPDG)               // ==> "Meson/Baryon+pi" case
          {
            G4int piPDG=111;
            G4double mpi=mPi0;
            G4int mbPDG=totPDG;
            G4double mbm=totM;
            if(totPDG==1114)
            {
              piPDG=-211;
              mpi=mPi;
              mbPDG=2112;
              mbm=mNeut;
            }
            else if(totPDG==2224)
            {
              piPDG=211;
              mpi=mPi;
              mbPDG=2212;
              mbm=mProt;
            }
            else if(totPDG==113)
            {
              piPDG=-211;
              mpi=mPi;
              mbPDG=211;
              mbm=mPi;
            }
            G4LorentzVector h4Mom(0.,0.,0.,mbm);
            G4LorentzVector g4Mom(0.,0.,0.,mpi);
            if(!G4QHadron(tot4M).DecayIn2(h4Mom, g4Mom))
            {
              G4cout<<"---Warning---G4QEnv::HadronizeQEnv: h="<<mbPDG<<"(m="<<mbm
                    <<") + pi(m="<<mpi<<")="<<mbm+mpi<<" > mTot="<<totMass<<G4endl;
              G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
              CleanUp();
              if(!CheckGroundState(quasH,true))
              {
                G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef debug
                G4cout<<"***G4QE::HQE:FillAsIs(0),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                //throw G4QException("G4QEnvironment::HadronizeQE: DecIn2 mB+nPi failed");
              }
              delete quasH;  
              return theQHadrons;
            }
#ifdef debug
            G4cout<<"G4QE::HQE:"<<tot4M<<"->h="<<mbPDG<<h4Mom<<"+p="<<piPDG<<g4Mom<<G4endl;
#endif
            G4QHadron* curH = new G4QHadron(mbPDG,h4Mom);
            if(totPDG==92000000||totPDG==90002000||totPDG==90000002)
              theEnvironment.DecayDibaryon(curH,&theQHadrons);
            else theQHadrons.push_back(curH);    // Fill the baryon (delete equivalent)
            G4QHadron* curG = new G4QHadron(piPDG,g4Mom);
#ifdef edebug
            G4cout<<"G4QEnv::HadrQEnv:Gamma/Pi0H="<<piPDG<<g4Mom<<G4endl;
#endif
            theQHadrons.push_back(curG);         // Fill the pi0 (delete equivalent)
            return theQHadrons;
          }
          else                                   // ==> "|B|<2 new Quasmon" case
          {
            G4Quasmon* resid = new G4Quasmon(totQC,tot4M); // delete is 3 lines below <-+
            G4QNucleus vacuum_value(90000000);         //                                     ^
            G4QHadronVector* curout=resid->Fragment(vacuum_value,1);// **!!DESTROY!!**<-<-+   ^
            G4int rest = resid->GetStatus();     // New status after fragm attempt  ^   ^
            if(!rest) eCount--;                  // Dec ExistingQuasmonsCounter     ^   ^
            delete resid;                        //_________________________________^___^
            G4int nHadrons = curout->size();     // a#of Hadrons in the outHV       ^
            if(nHadrons>0)                       // Transfer QHadrons to Output     ^
            {
              for (G4int ih=0; ih<nHadrons; ih++)// LOOP over output QHadrons       ^
              {                                  //                                 ^
#ifdef debug
                G4cout<<"G4QEnv::HadrQE:NewB<2, H#"<<ih //                          ^
                      <<", QPDG="<<(*curout)[ih]->GetQPDG() //                      ^
                      <<", 4M="<<(*curout)[ih]->Get4Momentum()<<G4endl; //          ^
#endif
                //theQHadrons.push_back(curout->operator[](ih));//(delete equ.) <-<-^
                theQHadrons.push_back((*curout)[ih]);          // (delete equ.) <-<-^
              }                                  //                                 ^
            }                                    //                                 ^
            else                                 //                                 ^
            {                                    //                                 ^
              G4ExceptionDescription ed;         //                                 ^
              ed << " Quasmon decay? : MQ=" << tot4M.m() << ",QC=" << totQC //      ^
                 << G4endl;
              G4Exception("G4QEnvironment::HadronizeQEnvironment()", "HAD_CHPS_0009",
                          FatalException, ed);
            }                                    // *** Do not destroy instances ***^
            curout->clear();                     // The instances are filled above  ^
            delete curout;                       // >-------------->--------------->+
            return theQHadrons;
          }
        }
        else
        {
          G4QContent tQC    =totQC;              // Not subtracted copy for error prints
          G4int      NSi    =0;                  // a#of additional Sigma
          G4int      SiPDG  =0;                  // PDG of additional Sigma
          G4double   MSi    =0.;                 // TotalMass of additional Sigma
          G4int      NaK    =0;                  // a#of additional Kaons/anti-Kaons
          G4int      aKPDG  =0;                  // PDG of additional Kaons/anti-Kaons
          G4double   MaK    =0.;                 // TotalMass of additionalKaons/anti-Kaons
          G4int      NPi    =0;                  // a#of additional pions
          G4int      PiPDG  =0;                  // PDG of additional pions
          G4double   MPi    =0.;                 // Total Mass of additional pions
          if    (totBN>0&&totS<0&&totChg+totChg>=totBN)// => "additional K+" case
          {
            aKPDG=321;
            NaK=-totS;
            MaK=mK*NaK;
            tQC+=totS*KpQC;
            totChg+=totS;                        // Charge reduction (totS<0!)
            totS=0;                              // Anti-strangness goes to anti-Kaons
          }
          else if (totBN>0&&totS<0)              // => "additional aK0" case
          {
            aKPDG=311;
            NaK=-totS;
            MaK=mK0*NaK;
            tQC+=totS*K0QC;
            totS=0;                              // Anti-strangness goes to anti-Kaons
          }
          else if (totBN>1&&totS>0&&(totChg<0||totChg>totBN-totS))//=>"additional Sigma"
          {
            NSi=totS;                            // Prototype of a#of Sigmas
            if(totChg<0)                         // Negative Sigmas
            {
              SiPDG=3112;
              if(-totChg<NSi) NSi=-totChg;       // A#of Sigma- is restricted by charge
              MSi=mSigM*NSi;                     // Total mass of Sigma-'s
              tQC-=NSi*SiMQC;                    // Subtract QC of Sigma-'s from totQC
              totChg+=NSi;                       // Increase the TotalResidualCharge
            }
            else
            {
              SiPDG=3222;                        // Positive Sigmas
              G4int exChg=totChg-totBN+totS;     // Excesive positive charge
              if(exChg<NSi) NSi=exChg;           // A#of Sigma+ is restricted by charge
              MSi=mSigP*NSi;                     // Total mass of Sigma+'s
              tQC-=NSi*SiPQC;                    // Subtract QC of Sigma-'s from totQC
              totChg-=NSi;                       // Reduce the TotalResidualCharge
            }
            totS-=NSi;                           // Reduce the TotalResidualStrangeness
            totBN-=NSi;                          // A#of excessive pions is added below
          }
          else if (totBN>0&&totS>totBN&&totBN<totS+totChg)// => "additional K0" case
          {// @@ Here Ksi0 check should be added totS=2>totBN=1&&totBN=1<totS=2+totChg=0
            aKPDG=-311;
            NaK=totS-totBN;
            MaK=mK0*NaK;
            tQC+=NaK*K0QC;
            totS-=NaK;                           // Reduce residualstrangeness
          }
          else if (totBN>0&&totS>totBN&&totChg<0)// => "additional K-" case
          {// @@ Here Ksi- check should be added totS=2>totBN=1&&totChg=-1<0
            aKPDG=-321;
            NaK=totS-totBN;
            MaK=mK0*NaK;
            tQC+=NaK*KpQC;
            totChg+=NaK;                         // Increase residual charge
            totS-=NaK;                           // Reduce residual strangeness
          }
          // === Now residual DELTAS should be subtracted === 
          if      (totBN>0&&totChg>totBN-totS)   // => "additional PI+" case
          {// @@ Here Sigma+ check should be added totChg=1>totBn=1-totS=1
            PiPDG=211;
            NPi=totChg-totBN+totS;
            MPi=mPi*NPi;
            tQC-=NPi*PiQC;
            totChg-=NPi;
          }
          else if (totBN>0&&totChg<0)            // => "additional PI-" case
          {// @@ Here Sigma- check should be added totChg<0
            PiPDG=-211;
            NPi=-totChg;
            MPi=mPi*NPi;
            tQC+=NPi*PiQC;                       // Now anti-Pions must be subtracted
            totChg+=NPi;
          }
          else if (!totBN&&totChg>1-totS)        // => "additional PI+" case
          {// @@ Here Sigma+ check should be added totChg=1>totBn=1-totS=1
            PiPDG=211;
            NPi=totChg+totS-1;
            MPi=mPi*NPi;
            tQC-=NPi*PiQC;
            totChg-=NPi;
          }
          else if (!totBN&&totChg<-1-totS)       // => "additional PI-" case
          {// @@ Here Sigma- check should be added totChg<0
            PiPDG=-211;
            NPi-=totChg+totS+1;
            MPi=mPi*NPi;
            tQC+=NPi*PiQC;                       // Now anti-Pions must be subtracted
            totChg+=NPi;
          }
          G4double      totRM=0.;                // min (GS) Mass of the Residual System
          if(totBN<2)                            // Calculate totPDG & totRM
          {
            totPDG=tQC.GetSPDGCode();            // MinPDGCode for the Residual compound
            if(totPDG==10&&tQC.GetBaryonNumber()>0) totPDG=tQC.GetZNSPDGCode();
            if(totPDG) totRM=G4QPDGCode(totPDG).GetMass(); // minMass of theResidualSystem
            else
            {
              // G4cerr<<"***G4QEnvironment::HadronizeQEnv: totPDG=0"<<G4endl;
              // throw G4QException("G4QEnv::HadrQEnv: Impossible PDG for B=1");
              G4ExceptionDescription ed;
              ed << "Impossible PDG for B=1: totPDG=0" << G4endl;
              G4Exception("G4QEnvironment::HadronizeQEnvironment()",
                          "HAD_CHPS_0010", FatalException, ed);
            }
          }
          else
          {
            G4QNucleus totN_temporary(tQC,tot4M);// Excited nucleus for the Residual System
            totN=totN_temporary;
            totRM=totN.GetMZNS();                // min (GS) Mass of the Residual System
            totPDG=totN.GetPDG();                // Total PDG Code for the Current compound
          }
          if(NaK)                                // ==> "Decay in K0 or K+ + NPi" case
          {//@@ Can (must) be moved to EvaporateResidual ??
            if(!NPi)                             // ==> "Only anti-strange K" case
            {
              G4LorentzVector m4Mom(0.,0.,0.,MaK);
              G4LorentzVector n4Mom(0.,0.,0.,totRM);
              G4double sum=MaK+totRM;
              if(fabs(totMass-sum)<eps)
              {
                m4Mom=tot4M*(MaK/sum);
                n4Mom=tot4M*(totRM/sum);
              }
              else if(totMass<sum || !G4QHadron(tot4M).DecayIn2(m4Mom, n4Mom))
              {
#ifdef edebug
                G4cout<<"***G4QE::HadronizeQE:M="<<aKPDG<<"(m="<<MaK<<")+N="<<totPDG<<"(m="
                      <<totRM<<")="<<sum<<" > mSN="<<totMass<<",d="<<sum-totMass<<G4endl;
#endif
                G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                CleanUp();
                if(!CheckGroundState(quasH,true))
                {
                  G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                  theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef debug
                  G4cout<<"***G4QEnv::HQE:FillAsItIs(1),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                  //throw G4QException("G4QEnvironment::HadronizeQEnv:AntiS-Nuc error");
                }
                delete quasH;  
                return theQHadrons;
              }
#ifdef debug
              G4cout<<"G4QEnv::HadronizeQEnv: SN="<<tot4M<<" -> M="
                    <<aKPDG<<m4Mom<<" + N="<<totPDG<<n4Mom<<totQC<<G4endl;
#endif
              G4LorentzVector oneK=m4Mom;        // 4-mom of only kaon  
              if(NaK>1) oneK = m4Mom/NaK;        // 4-mom of one kaon  
              for (G4int jp=0; jp<NaK; jp++)
              {
                G4QHadron* curK = new G4QHadron(aKPDG,oneK);
                theQHadrons.push_back(curK);     // Fill the curK (delete equivalent)
              }
              G4QHadron* curN = new G4QHadron(totPDG,n4Mom); // @@ Use DecayDib then Evap
              EvaporateResidual(curN);           // Try to evaporate residual (del.eq.)
            }
            else                                 // ==> "Anti-strange K's + DELTA's" case
            {
              G4LorentzVector m4Mom(0.,0.,0.,MPi);
              G4LorentzVector k4Mom(0.,0.,0.,MaK);
              G4LorentzVector n4Mom(0.,0.,0.,totRM);
              if(!G4QHadron(tot4M).DecayIn3(m4Mom, k4Mom, n4Mom))
              {
                G4cout<<"---Warning---G4QE::HadronQE:K="<<aKPDG<<"(m="<<MaK<<")+PI="<<PiPDG
                      <<"(m="<<MPi<<")+N="<<totPDG<<"(m="<<totRM<<")>tM="<<totMass<<G4endl;
                G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                CleanUp();
                if(!CheckGroundState(quasH,true))
                {
                  G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                  theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef debug
                  G4cout<<"***G4QEnv::HQE:FillAsItIs(2),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                  //throw G4QException("G4QEnvironment::HadronizeQE:2AntiS-Nucl(1) error");
                }
                delete quasH;  
                return theQHadrons;
              }
#ifdef fdebug
              G4cout<<"G4QEnv::HadronizeQEnv: SN="<<tot4M<<" -> nK="<<aKPDG<<k4Mom
                    <<" + nPi="<<PiPDG<<m4Mom<<" + N="<<totPDG<<n4Mom<<G4endl;
#endif
              G4LorentzVector onePi=m4Mom;       // 4-mom of only pion  
              if(NPi>1) onePi = m4Mom/NPi;       // 4-mom of one pion  
              for (G4int ip=0; ip<NPi; ip++)
              {
                G4QHadron* curP = new G4QHadron(PiPDG,onePi);
#ifdef debug
                G4cout<<"G4QEnv::HadrQEnv:SPion#"<<ip<<",H="<<PiPDG<<onePi<<G4endl;
#endif
                theQHadrons.push_back(curP);     // Fill the curM (delete equivalent)
              }
              G4LorentzVector oneK=k4Mom;        // 4-mom of one kaon  
              if(NaK>1) oneK = k4Mom/NaK;        // 4-mom of one kaon  
              for (G4int jp=0; jp<NaK; jp++)
              {
                G4QHadron* curP = new G4QHadron(aKPDG,oneK);
#ifdef debug
                G4cout<<"G4QEnv::HadrQEnv:Kaon#"<<jp<<",H="<<aKPDG<<oneK<<G4endl;
#endif
                theQHadrons.push_back(curP);     // Fill the curM (delete equivalent)
              }
              G4QHadron* curN = new G4QHadron(totPDG,n4Mom);
              EvaporateResidual(curN);           // Try to evaporate residual (del.equiv.)
            }
            return theQHadrons;
          }
          else if(NSi)                           // ==> "Decay in Sig+ or Sig- + NPi" case
          {//@@ Can (must) be moved to EvaporateResidual ??
            if(!NPi)                             // ==> "Only Sigma's" case
            {
              G4LorentzVector m4Mom(0.,0.,0.,MSi);
              G4LorentzVector n4Mom(0.,0.,0.,totRM);
              G4double sum=MSi+totRM;
              if(fabs(totMass-sum)<eps)
              {
                m4Mom=tot4M*(MSi/sum);
                n4Mom=tot4M*(totRM/sum);
              }
              else if(totMass<sum || !G4QHadron(tot4M).DecayIn2(m4Mom, n4Mom))
              {
#ifdef edebug
                G4cout<<"***G4QE::HadronizeQE:M="<<aKPDG<<"(s="<<MSi<<")+N="<<totPDG<<"(m="
                      <<totRM<<")="<<sum<<" > mSN="<<totMass<<",d="<<sum-totMass<<G4endl;
#endif
                G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                CleanUp();
                if(!CheckGroundState(quasH,true))
                {
                  G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                  theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef debug
                  G4cout<<"***G4QEnv::HQE:FillAsItIs(2),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                  //throw G4QException("G4QEnvironment::HadronizeQEnv:Sigma-Nuc error");
                }
                delete quasH;  
                return theQHadrons;
              }
#ifdef debug
              G4cout<<"G4QEnv::HadronizeQEnv: SN="<<tot4M<<" -> Sig="
                    <<SiPDG<<m4Mom<<" + N="<<totPDG<<n4Mom<<totQC<<G4endl;
#endif
              G4LorentzVector oneS=m4Mom;        // 4-mom of the only sigma  
              if(NSi>1) oneS = m4Mom/NSi;        // 4-mom of one sigma  
              for (G4int jp=0; jp<NSi; jp++)
              {
                G4QHadron* curS = new G4QHadron(SiPDG,oneS);
                theQHadrons.push_back(curS);     // Fill the curS (delete equivalent)
              }
              G4QHadron* curN = new G4QHadron(totPDG,n4Mom); // @@ Use DecayDib then Evap
              EvaporateResidual(curN);           // Try to evaporate residual (del.eq.)
            }
            else                                 // ==> "Sigma's + DELTA's" case
            {
              G4LorentzVector m4Mom(0.,0.,0.,MPi);
              G4LorentzVector k4Mom(0.,0.,0.,MSi);
              G4LorentzVector n4Mom(0.,0.,0.,totRM);
              if(!G4QHadron(tot4M).DecayIn3(m4Mom, k4Mom, n4Mom))
              {
                G4cout<<"---Warning---G4QE::HadronQE:S="<<SiPDG<<"(m="<<MSi<<")+PI="<<PiPDG
                      <<"(m="<<MPi<<")+N="<<totPDG<<"(m="<<totRM<<")>tM="<<totMass<<G4endl;
                G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                CleanUp();
                if(!CheckGroundState(quasH,true))
                {
                  G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                  theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef debug
                  G4cout<<"***G4QEnv::HQE:FillAsItIs(3),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                  //throw G4QException("G4QEnvironment::HadronizeQE:2Sigma-Nucl(1) error");
                }
                delete quasH;  
                return theQHadrons;
              }
#ifdef fdebug
              G4cout<<"G4QEnv::HadronizeQEnv: SN="<<tot4M<<" -> nS="<<SiPDG<<k4Mom
                    <<" + nPi="<<PiPDG<<m4Mom<<" + N="<<totPDG<<n4Mom<<G4endl;
#endif
              G4LorentzVector onePi=m4Mom;       // 4-mom of the only pion  
              if(NPi>1) onePi = m4Mom/NPi;       // 4-mom of one pion  
              for (G4int ip=0; ip<NPi; ip++)
              {
                G4QHadron* curP = new G4QHadron(PiPDG,onePi);
#ifdef debug
                G4cout<<"G4QEnv::HadrQEnv:SPion#"<<ip<<",H="<<PiPDG<<onePi<<G4endl;
#endif
                theQHadrons.push_back(curP);     // Fill the curM (delete equivalent)
              }
              G4LorentzVector oneS=k4Mom;        // 4-mom of the only kaon  
              if(NSi>1) oneS = k4Mom/NSi;        // 4-mom of one kaon  
              for (G4int jp=0; jp<NSi; jp++)
              {
                G4QHadron* curP = new G4QHadron(SiPDG,oneS);
#ifdef debug
                G4cout<<"G4QEnv::HadrQEnv:Sigma#"<<jp<<",H="<<SiPDG<<oneS<<G4endl;
#endif
                theQHadrons.push_back(curP);     // Fill the curM (delete equivalent)
              }
              G4QHadron* curN = new G4QHadron(totPDG,n4Mom);
              EvaporateResidual(curN);           // Try to evaporate residual (del.equiv.)
            }
            return theQHadrons;
          }
          else if(NPi)                           // ==> "Decay in Pi+ or Pi-" case
          {
            if(NPi==1)                           // ==> "One isobar" case
            {
              G4LorentzVector m4Mom(0.,0.,0.,MPi);
              G4LorentzVector n4Mom(0.,0.,0.,totRM);
              if(!G4QHadron(tot4M).DecayIn2(m4Mom, n4Mom))
              {
                G4cout<<"---Warning---G4QEnv::HadronizeQEnv:M="<<PiPDG<<"(m="<<MPi<<")+N="
                      <<totPDG<<"(m="<<totRM<<")="<<MPi+totRM<<" > mSN="<<totMass<<G4endl;
                G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                CleanUp();
                if(!CheckGroundState(quasH,true))
                {
                  G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                  theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef debug
                  G4cout<<"***G4QEnv::HQE:FillAsItIs(5),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                  //throw G4QException("G4QEnvironment::HadronizeQEnv:Iso-Nucleus error");
                }
                delete quasH;  
                return theQHadrons;
              }
#ifdef debug
              G4cout<<"G4QEnv::HadronizeQEnv: SN="<<tot4M<<" -> M="<<PiPDG<<m4Mom<<" + N="
                    <<totPDG<<n4Mom<<totQC<<G4endl;
#endif
              G4QHadron* curK = new G4QHadron(PiPDG,m4Mom);
              theQHadrons.push_back(curK);       // Fill the curK (delete equivalent)
              G4QHadron* curN = new G4QHadron(totPDG,n4Mom);
              EvaporateResidual(curN);           // Evaporate residual (delete equivalent)
            }
            else                                 // ==> "Many Isobars" case
            {
              G4int N1Pi = NPi/2;                // First pion cluster
              G4int N2Pi = NPi-N1Pi;             // Second pion cluster
              G4double mM  = MPi/NPi;            // Mass of Pi
              G4double m1M = mM*N1Pi;            // Mass of the first Pi-cluster
              G4double m2M = mM*N2Pi;            // Mass of the second Pi-cluster
              G4LorentzVector m4Mom(0.,0.,0.,m1M);
              G4LorentzVector k4Mom(0.,0.,0.,m2M);
              G4LorentzVector n4Mom(0.,0.,0.,totRM);
              if(!G4QHadron(tot4M).DecayIn3(m4Mom, k4Mom, n4Mom))
              {
                G4cout<<"---Warning---G4QEnv::HadronizeQE:N*Pi="<<PiPDG<<"(m="<<mM<<")+N="
                      <<totPDG<<"(m="<<totRM<<") >(?)SN="<<totMass<<G4endl;
                G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                CleanUp();
                if(!CheckGroundState(quasH,true))
                {
                  G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                  theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef debug
                  G4cout<<"***G4QEnv::HQE:FillAsItIs(5),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                  //throw G4QException("G4QEnvironment::HadronizeQE:ManyIsoNucleus error");
                }
                delete quasH;  
                return theQHadrons;
              }
#ifdef debug
              G4cout<<"G4QEnv::HadronizeQEnv: SN="<<tot4M<<" -> N*PI="<<PiPDG
                    <<" (4M1="<<m4Mom<<" + 4M2="<<k4Mom<<") + N="<<totPDG<<n4Mom<<G4endl;
#endif
              G4LorentzVector one1=m4Mom;        // 4-mom of the 1st cluster only pion  
              if(N1Pi>1) one1=m4Mom/N1Pi;        // 4-mom of the 1st cluster one pion  
              for (G4int ip=0; ip<N1Pi; ip++)
              {
                G4QHadron* curP = new G4QHadron(PiPDG,one1);
                theQHadrons.push_back(curP);     // Fill the curP (delete equivalent)
              }
              G4LorentzVector one2=k4Mom;        // 4-mom of the 2nd cluster only pion  
              if(N2Pi>1) one2=k4Mom/N2Pi;        // 4-mom of the 2nd cluster one pion  
              for (G4int jp=0; jp<N2Pi; jp++)
              {
                G4QHadron* curP = new G4QHadron(PiPDG,one2);
                theQHadrons.push_back(curP);     // Fill the curP (delete equivalent)
              }
              G4QHadron* curN = new G4QHadron(totPDG,n4Mom);
              EvaporateResidual(curN);           // Try to evaporate residual (del.equiv.)
            }
            return theQHadrons;
          }
        }
#ifdef fdebug
        G4cout<<"G4QE::HadrQEnv: Try FinalEvaporation t4M="<<tot4M<<",tQC="<<totQC<<G4endl;
#endif
        CleanUp();
        G4QHadron* evH = new G4QHadron(totQC,tot4M); // Create a Hadron for theResidualNucl
        EvaporateResidual(evH);                  // Try to evaporate residual (del.equiv.)
        return theQHadrons;
      }
      else                                       // ==> "Only with GSEnvironment" case
      { 
        if(totPDG==90000000 || fabs(totMass)<0.000001)
        {
          CleanUp();
          return theQHadrons;
        }
        G4double dM=totMass-totM;
#ifdef debug
        G4cout<<"G4QEnv::HadrQEnv:GroundState tM-GSM="<<dM<<",GSM="<<totM<<",tPDG="<<totPDG
              <<",nQ="<<nQuasmons<<G4endl;
#endif
        G4Quasmon*       pQ = theQuasmons[0];    // Pointer to the first Quasmon          
        G4QPDGCode    QQPDG = pQ->GetQPDG();     // QPDG of the first Quasmon
        G4int          QPDG = QQPDG.GetPDGCode();
        G4QNucleus    totRN(totQC,tot4M);        // Nucleus for theTotalResidualNuclearComp
        G4int          spbRN=totRN.SplitBaryon();// PossibilityToSplit baryon from Residual
        if(dM>-0.001)
        {
#ifdef fdebug
          G4cout<<"G4QE::HadrQE:ExcitedNucleus, dM="<<dM<<">0,tBN="<<totBN<<",nQ="<<G4endl;
#endif
          CleanUp();
          G4QHadron* evH = new G4QHadron(totQC,tot4M);// Create a Hadron for ResidualNucl
          EvaporateResidual(evH);                // Try to evaporate residual (del. equiv.)
        }
        else if(nQuasmons==1&&QPDG!=22&&QPDG!=111)// => "Decay Quasmon or Q+Environ" case
        {
          G4int envPDG = theEnvironment.GetPDG();// PDGCode of the NuclQEnvironment
#ifdef debug
          G4cout<<"G4QEnv::HadrQEnv: nQ=1, QPDG=="<<QPDG<<G4endl;
#endif
          if(!QPDG)
          {
            G4ExceptionDescription ed;
            ed <<"(2)Can't Decay QEnv: Quasmon is an unknown QHadron: PDG="<< QPDG<<G4endl;
            G4Exception("G4QEnvironment::HadronizeQEnvironmen()", "HAD_CHPS_0011",
                        FatalException, ed);
          }
          // => "Quasmon-Chipolino or Environment-Dibaryon" case
          else if(QPDG==10||QPDG==92000000||QPDG==90002000||QPDG==90000002)
          {
            G4QContent QQC = pQ->GetQC();        // Quark Content of the Quasmon
            G4QPDGCode h1QPDG=nQPDG;             // QPDG of the first hadron
            G4double   h1M   =mNeut;             // Mass of the first hadron
            G4QPDGCode h2QPDG=h1QPDG;            // QPDG of the second hadron
            G4double   h2M   =mNeut;             // Mass of the second hadron
            if(QPDG==10)
            {
              G4QChipolino QChip(QQC);           // define the Quasmon as a Chipolino
              h1QPDG=QChip.GetQPDG1();           // QPDG of the first hadron
              h1M   =h1QPDG.GetMass();           // Mass of the first hadron
              h2QPDG=QChip.GetQPDG2();           // QPDG of the second hadron
              h2M   =h2QPDG.GetMass();           // Mass of the second hadron
            }
            else if(QPDG==90002000)
            {
              h1QPDG=pQPDG;                      // QPDG of the first hadron
              h1M   =mProt;                      // Mass of the first hadron
              h2QPDG=h1QPDG;                     // QPDG of the second hadron
              h2M   =mProt;                      // Mass of the second hadron
            }
            else if(QPDG==92000000)
            {
              h1QPDG=lQPDG;                      // QPDG of the first hadron
              h1M   =mLamb;                      // Mass of the first hadron
              h2QPDG=h1QPDG;                     // QPDG of the second hadron
              h2M   =mLamb;                      // Mass of the second hadron
              G4double ddMass=totMass-envM;      // Free CM energy
              if(ddMass>mSigZ+mSigZ)             // Sigma0+Sigma0 is possible
              {                                  // @@ Only two particles PS is used
                G4double dd2=ddMass*ddMass;      // Squared free energy
                G4double sma=mLamb+mLamb;        // Lambda+Lambda sum
                G4double pr1=0.;                 // Prototype to avoid sqrt(-)
                if(ddMass>sma) pr1=sqrt((dd2-sma*sma)*dd2); // Lamb+Lamb PS
                sma=mLamb+mSigZ;                 // Lambda+Sigma0 sum
                G4double smi=mSigZ-mLamb;        // Sigma0-Lambda difference
                G4double pr2=pr1;                // Prototype of +L+S0 PS
                if(ddMass>sma && ddMass>smi) pr2+=sqrt((dd2-sma*sma)*(dd2-smi*smi));
                sma=mSigZ+mSigZ;                 // Sigma0+Sigma0 sum
                G4double pr3=pr2;                // Prototype of +Sigma0+Sigma0 PS
                if(ddMass>sma) pr3+=sqrt((dd2-sma*sma)*dd2);
                G4double hhRND=pr3*G4UniformRand(); // Randomize PS
                if(hhRND>pr2)                    // --> "ENnv+Sigma0+Sigma0" case
                {                                //
                  h1QPDG=s0QPDG;                 // QPDG of the first hadron
                  h1M   =mSigZ;                  // Mass of the first hadron
                  h2QPDG=h1QPDG;                 // QPDG of the second hadron
                  h2M   =mSigZ;                  // Mass of the second hadron
                }                                //
                else if(hhRND>pr1)               // --> "ENnv+Sigma0+Lambda" case
                {                                //
                  h1QPDG=s0QPDG;                 // QPDG of the first hadron
                  h1M   =mSigZ;                  // Mass of the first hadron
                }                                //
              }                                  //
              else if(ddMass>mSigZ+mLamb)        // Lambda+Sigma0 is possible
              {                                  // @@ Only two particles PS is used
                G4double dd2=ddMass*ddMass;      // Squared free energy
                G4double sma=mLamb+mLamb;        // Lambda+Lambda sum
                G4double pr1=0.;                 // Prototype to avoid sqrt(-)
                if(ddMass>sma) pr1=sqrt((dd2-sma*sma)*dd2); // Lamb+Lamb PS
                sma=mLamb+mSigZ;                 // Lambda+Sigma0 sum
                G4double smi=mSigZ-mLamb;        // Sigma0-Lambda difference
                G4double pr2=pr1;                // Prototype of +L+S0 PS
                if(ddMass>sma && ddMass>smi) pr2+=sqrt((dd2-sma*sma)*(dd2-smi*smi));
                if(pr2*G4UniformRand()>pr1)      // --> "ENnv+Sigma0+Lambda" case
                {                                //
                  h1QPDG=s0QPDG;                 // QPDG of the first hadron
                  h1M   =mSigZ;                  // Mass of the first hadron
                }                                //
              }                                  //
            }                                    //
            if(h1M+h2M+envM<totMass)             // => "Three parts decay" case
            {
              G4LorentzVector h14M(0.,0.,0.,h1M);
              G4LorentzVector h24M(0.,0.,0.,h2M);
              G4LorentzVector e4M(0.,0.,0.,envM);
              if(!G4QHadron(tot4M).DecayIn3(h14M,h24M,e4M))
              {
                G4cout<<"Warning->G4QE::HQE:M="<<tot4M.m()<<","<<totMass<<"->"<<h1QPDG<<"("
                      <<h1M<<")+"<<h1QPDG<<"("<<h2M<<")+"<<envM<<"="<<h1M+h2M+envM<<G4endl;
                G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                CleanUp();
                if(!CheckGroundState(quasH,true))
                {
                  G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                  theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef debug
                  G4cout<<"***G4QEnv::HQE:FillAsItIs(6),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                  //throw G4QException("G4QEnv::HadrQEnv:QChipo+Environment DecIn3 Error");
                }
                delete quasH;  
                return theQHadrons;

              }
              G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);
              theQHadrons.push_back(h1H);        // (delete equivalent)
#ifdef debug
              G4cout<<"G4QE::HQE:(2) H1="<<h1QPDG<<h14M<<G4endl;
#endif
              G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
              theQHadrons.push_back(h2H);        // (delete equivalent)
#ifdef debug
              G4cout<<"G4QE::HQE:(2) H2-"<<h2QPDG<<h24M<<G4endl;
#endif
              G4QHadron* qeH = new G4QHadron(envPDG,e4M);
              theQHadrons.push_back(qeH);        // (delete equivalent)
#ifdef debug
              G4cout<<"G4QE::HQE:(2) QEenv="<<envPDG<<e4M<<G4endl;
#endif
            }
#ifdef fdebug
            G4cout<<"***G4QEnv::HadQEnv:tM="<<tot4M.m()<<totQC<<"< h1="<<h1QPDG<<"(M="<<h1M
                  <<")+h2="<<h1QPDG<<"(M="<<h2M<<")+eM="<<envM<<"="<<h1M+h2M+envM<<G4endl;
            //throw G4QException("G4QEnv::HadrQEnv:QChipo+Env mass > than decaying mass");
#endif
            CleanUp();
            G4QHadron* evH = new G4QHadron(totQC,tot4M);// Create a Hadron for ResidualNucl
            EvaporateResidual(evH);              // Try to evaporate residual (del. equiv.)
            return theQHadrons;
          }
          else                                   // No environment
          {
            G4int nHadrs=theQHadrons.size();     // #of available hadrons
            for(G4int ih=0; ih<nHadrs; ++ih)
            {
              G4QHadron* ch=theQHadrons[ih];
              G4LorentzVector ch4M=ch->Get4Momentum();
              G4double chM=ch4M.m();
              G4LorentzVector tch4M=ch4M+tot4M;
              if(tch4M.m() > chM + totM)         // Can be corrected
              {
                G4LorentzVector h14M(0.,0.,0.,chM);
                G4LorentzVector h24M(0.,0.,0.,totM);
                if(!G4QHadron(tch4M).DecayIn2(h14M,h24M))
                {
                  G4cout<<"-Warning->G4QE::HQE:M="<<tch4M.m()<<"->"<<chM<<"+"<<totM<<"="
                        <<chM+totM<<G4endl;
                }
                else
                {
                  tot4M=h24M;                    // Change the residual 4M
                  ch->Set4Momentum(h14M);        // Change 4M of the current hadron
                  break;                         // Quit the loop
                }
              }
            }
            G4QHadron* rH = new G4QHadron(totQC,tot4M);// Create a Hadron for ResidualNucl
            theQHadrons.push_back(rH);
          }
        }
        else if(spbRN)// => "Join all quasmons to the residual compound and evaporate" case
        {
#ifdef fdebug
          G4cout<<"***G4QEnv::HadQEnv: Evaporate the total residual tRN="<<totRN<<G4endl;
#endif
          CleanUp();
          G4QHadron* evH = new G4QHadron(totQC,tot4M);// Create Hadron for theResidNucleus
          EvaporateResidual(evH);               // Try to evaporate residual (del.equiv.)
          return theQHadrons;
        }
        //else if(nQuasmons<3||theQHadrons.size()<12)//"Try to correct" case (change cond)
        else if(2>3)  // "Try to correct" case (change condition)
        {
#ifdef debug
          G4cout<<"***G4QEnv::HadrQE: M="<<totMass<<",dM="<<dM<<",nQ="<<nQuasmons<<G4endl;
#endif
          G4int          nOfOUT  = theQHadrons.size();
          while(nOfOUT)
          {
            G4QHadron*     theLast = theQHadrons[nOfOUT-1];
            G4LorentzVector last4M = theLast->Get4Momentum();
            G4QContent      lastQC = theLast->GetQC();
            G4int           lastS  = lastQC.GetStrangeness();
            totS                   = totQC.GetStrangeness();
            G4int           nFr    = theLast->GetNFragments();
            G4int           gam    = theLast->GetPDGCode();
            if(gam!=22&&!nFr&&lastS<0&&lastS+totS<0&&nOfOUT>1) // => "Skip K,gam & decayed"
            {
              G4QHadron* thePrev = theQHadrons[nOfOUT-2];
              theQHadrons.pop_back();         // the last QHadron is excluded from OUTPUT
              theQHadrons.pop_back();         // the prev QHadron is excluded from OUTPUT
              theQHadrons.push_back(thePrev); // thePast becomes theLast as an instance
              delete    theLast;              // theLast QHadron is deleated as an instance
              theLast = thePrev;              // Update parameters(thePrev becomes theLast)
              last4M = theLast->Get4Momentum();
              lastQC = theLast->GetQC();
            }
            else
            {
              theQHadrons.pop_back();         // the last QHadron is excluded from OUTPUT 
              delete         theLast;         // theLastQHadron is deleated as an instance
            }
            totQC+=lastQC;                    // Update (increase) the total QC
            tot4M+=last4M;                    // Update (increase) the total 4-momentum
            totMass=tot4M.m();                // Calculate new real total mass
            G4int bn=totQC.GetBaryonNumber(); // The BaryNum after addition
            totPDG=totQC.GetSPDGCode();
            if(totPDG==10&&totQC.GetBaryonNumber()>0) totPDG=totQC.GetZNSPDGCode();
            if(bn>1)
            {
              totS  =totQC.GetStrangeness();  // Total Strangeness of this System
              if(totS>=0)                     // => "This is a normal nucleus" case
              {
                G4QNucleus newN(totQC,tot4M);
                totPDG=newN.GetPDG();
                totM  =newN.GetMZNS();           // Calculate new minimum (GS) mass
              }
              else if(totS==-1)                  // => "Try to decay in K+/aK0 and finish"
              {
                G4double m1=mK;         
                G4int  PDG1=321;
                G4QNucleus  newNp(totQC-KpQC);
                G4int  PDG2=newNp.GetPDG();
                G4double m2_value=newNp.GetMZNS();
                G4QNucleus  newN0(totQC-K0QC);
                G4double m3_value=newN0.GetMZNS();
                if (m3_value+mK0<m2_value+mK)                // => "aK0+ResA is better" case
                {
                  m1  =mK0;
                  PDG1=311;
                  m2_value  =m3_value;
                  PDG2=newN0.GetPDG();
                }
                if(totMass>m1+m2_value)                // => "can decay" case
                {
                  G4LorentzVector fq4M(0.,0.,0.,m1);
                  G4LorentzVector qe4M(0.,0.,0.,m2_value);
                  if(!G4QHadron(tot4M).DecayIn2(fq4M,qe4M))
                  {
                    G4cout<<"---Warning---G4QE::HadQE:tM="<<tot4M.m()<<"->aK="<<PDG1<<"(M="
                          <<m1<<")+ResA="<<PDG2<<"(M="<<m2_value<<")="<<m1+m2_value<<G4endl;
                    G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                    CleanUp();
                    if(!CheckGroundState(quasH,true))
                    {
                      G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                      theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef debug
                      G4cout<<"***G4QE::HQE:FillAsIs(7),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                      //throw G4QException("G4QEnv::HadrQEnv: aK+ResA DecayIn2 error");
                    }
                    delete quasH;  
                    return theQHadrons;
                  }
                  G4QHadron* H1 = new G4QHadron(PDG1,fq4M);
#ifdef debug
                  G4cout<<"G4QE::HQE:Kaon(Env)="<<PDG1<<fq4M<<G4endl;
#endif
                  theQHadrons.push_back(H1);     // (delete equivalent)
                  G4QHadron* H2 = new G4QHadron(PDG2,qe4M);
#ifdef debug
                  G4cout<<"G4QE::HQE:ResidEnv="<<PDG2<<qe4M<<G4endl;
#endif
                  theQHadrons.push_back(H2);     // (delete equivalent)
                  break;
                }
                else totM=250000.;               // => "Continue reversion" case
              }
              else if(totS==-2)                  //=>"Try to decay in 2(K+/aK0) and finish"
              {
                G4double m1=mK;         
                G4int  PDG1=321;
                G4double m2_value=mK0;         
                G4int  PDG2=311;
                G4QNucleus  newNp0(totQC-KpQC-K0QC);
                G4int  PDG3=newNp0.GetPDG();
                G4double m3_value=newNp0.GetMZNS();    // M-K^0-K^+
                G4QNucleus  newN00(totQC-K0QC-K0QC);
                G4double m4=newN00.GetMZNS();    // M-2*K^0
                G4QNucleus  newNpp(totQC-KpQC-KpQC);
                G4double m5=newNpp.GetMZNS();    // M-2*K^+
                if (m4+mK0+mK0<m3_value+mK+mK0 && m4+mK0+mK0<=m5+mK+mK) //=>"2K0+ResA is theBest"
                {
                  m1  =mK0;
                  PDG1=311;
                  m3_value  =m4;
                  PDG3=newN00.GetPDG();
                }
                else if(m5+mK+mK<m3_value+mK+mK0 && m5+mK+mK<=m4+mK0+mK0)//=>"2Kp+ResA isTheBest"
                {
                  m2_value  =mK;
                  PDG1=321;
                  m3_value  =m5;
                  PDG3=newNpp.GetPDG();
                }
                if(totMass>m1+m2_value+m3_value)             // => "can decay" case
                {
                  G4LorentzVector k14M(0.,0.,0.,m1);
                  G4LorentzVector k24M(0.,0.,0.,m2_value);
                  G4LorentzVector ra4M(0.,0.,0.,m3_value);
                  if(!G4QHadron(tot4M).DecayIn3(k14M,k24M,ra4M))
                  {
                    G4cout<<"--Warning--G4QE::HQE:tM="<<tot4M.m()<<"->aK="<<PDG1<<"(M="<<m1
                          <<")+K2="<<PDG2<<"(M="<<m2_value<<")+A="<<PDG3<<"(M="<<m3_value<<")"<<G4endl;
                    G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                    CleanUp();
                    if(!CheckGroundState(quasH,true))
                    {
                      G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                      theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef debug
                      G4cout<<"***G4QE::HQE:FillAsIs(8),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                      //throw G4QException("G4QEnv::HadrQE:2K+ResidNucleus DecIn3 Error");
                    }
                    delete quasH;  
                    return theQHadrons;
                  }
                  G4QHadron* H1 = new G4QHadron(PDG1,k14M);
                  theQHadrons.push_back(H1);     // (delete equivalent)
#ifdef debug
                  G4cout<<"G4QE::HQE:K1(Env)="<<PDG1<<k14M<<G4endl;
#endif
                  G4QHadron* H2 = new G4QHadron(PDG2,k24M);
                  theQHadrons.push_back(H2);     // (delete equivalent)
#ifdef debug
                  G4cout<<"G4QE::HQE:K2(Env)="<<PDG2<<k24M<<G4endl;
#endif
                  G4QHadron* H3 = new G4QHadron(PDG3,ra4M);
                  theQHadrons.push_back(H3);     // (delete equivalent)
#ifdef debug
                  G4cout<<"G4QE::HQE:ResKKEnv="<<PDG3<<ra4M<<G4endl;
#endif
                  break;
                }
                else totM=270000.;               // => "Continue reversion" case
              }
              else totM=300000.;                 // => "Continue reversion" case
            }
            else
            {
              if     (totPDG==1114||totPDG==2224||totPDG==10) // Decay right now and finish
              {
                G4double m1=mNeut;
                G4int  PDG1=2112;
                G4double m2_value=mPi;
                G4int  PDG2=-211;
                if(totPDG==2224)
                {
                  m1=mProt;
                  PDG1=2212;
                  m2_value=mPi;
                  PDG2=211;
                }
                else if(totPDG==10)              // "Chipolino" case
                {
                  G4QChipolino resChip(totQC);   // define the residual as a Chipolino
                  G4QPDGCode h1=resChip.GetQPDG1();
                  PDG1=h1.GetPDGCode();          // PDG code of the first hadron
                  m1  =h1.GetMass();             // Mass of the first hadron
                  G4QPDGCode h2=resChip.GetQPDG2();
                  PDG2=h2.GetPDGCode();          // PDG code of the second hadron
                  m2_value  =h2.GetMass();       // Mass of the second hadron
                }
                if(totMass>m1+m2)
                {
                  G4LorentzVector fq4M(0.,0.,0.,m1);
                  G4LorentzVector qe4M(0.,0.,0.,m2_value);
                  if(!G4QHadron(tot4M).DecayIn2(fq4M,qe4M))
                  {
                    G4cout<<"---Warning---G4QE::HaQE:tM="<<tot4M.m()<<"-> h1="<<PDG1<<"(M="
                          <<m1<<") + h2="<<PDG2<<"(M="<<m2_value<<")="<<m1+m2_value<<G4endl;
                    G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                    CleanUp();
                    if(!CheckGroundState(quasH,true))
                    {
                      G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                      theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef debug
                      G4cout<<"***G4QE::HQE:FillAsIs(9),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                      //throw G4QException("G4QEnv::HadrQEnv: h1+h2 DecayIn2 Error");
                    }
                    delete quasH;  
                    return theQHadrons;
                  }
                  G4QHadron* H1 = new G4QHadron(PDG1,fq4M);
                  theQHadrons.push_back(H1);     // (delete equivalent)
#ifdef debug
                  G4cout<<"G4QE::HQE:h1="<<PDG1<<fq4M<<G4endl;
#endif
                  G4QHadron* H2 = new G4QHadron(PDG2,qe4M);
#ifdef debug
                  G4cout<<"G4QE::HQE:h2="<<PDG2<<qe4M<<G4endl;
#endif
                  theQHadrons.push_back(H2);     // (delete equivalent)
                  break;
                }
                else totM=350000.;
              }
              else if(totPDG) totM=G4QPDGCode(totPDG).GetMass();
              else totM=400000.;
            }
            totBN=totQC.GetBaryonNumber();      // The BaryNum after addition
            totS=totQC.GetStrangeness();        // The Strangeness after addition
            dM=totMass-totM;
#ifdef fdebug
            G4cout<<"G4QEnv::HadrQE: Add H="<<last4M<<lastQC<<",tM="<<tot4M<<totM<<totQC
                  <<",dM="<<dM<<", tB="<<totBN<<", tS="<<totS<<G4endl;
#endif
            if(dM>-0.001&&totPDG)
            {
              CleanUp();
              G4QHadron* evH = new G4QHadron(totPDG,tot4M);//Create Hadron for ResidNucleus
              EvaporateResidual(evH);           // Evaporate ResNuc (del.equiv)
              break;
            }
            nOfOUT  = theQHadrons.size();       // Update the value of OUTPUT entries
          } // End of WHILE(nOfOUT)
          nOfOUT  = theQHadrons.size();         // Update the value of OUTPUT entries
          if(!nOfOUT)
          {
            G4cout<<"---Warning---G4QEnv::HadrQE:M="<<totMass<<"<gsM="<<totM<<",dM="<<dM
                  <<", tPDG="<<totPDG<<", t4M="<<tot4M<<G4endl;
            // throw G4QException("G4QEnvironment::HadronizeQEnv:Can't decayExhostedQEnv");
            CleanUp();
            G4QHadron* evH = new G4QHadron(totPDG,tot4M);// Create Hadron for ResidNucleus
            EvaporateResidual(evH);             // Evaporate ResidNucl (del.equiv)
          }
        }
        else                                    // "Last decay was fatal" case @@ buggy ?MK
        {
#ifdef debug
          G4cout<<"***G4QEnv::HadrQE: M="<<totMass<<",dM="<<dM<<",nQ="<<nQuasmons<<G4endl;
#endif
          G4Quasmon* quasH = new G4Quasmon(totQC,tot4M);
          CleanUp();
          if(!CheckGroundState(quasH,true))
          {
            G4QHadron* hadr = new G4QHadron(totQC,tot4M);
#ifdef debug
            G4cout<<"G4QE::HQE:CheckGS failed H="<<totQC<<tot4M<<G4endl;
#endif
            theQHadrons.push_back(hadr); // Cor or fill asItIs
          }
          delete quasH;  
        }
        CleanUp();
      }
    } // End of "infinit" WHILE LOOP
  } // End of the "Nuclear Environment" case
  return theQHadrons;
} // End of the main member function HadronizeQEnvironment

// Clean up the QEnvironment to Zero
void G4QEnvironment::CleanUp()
{
  static const G4QNucleus vacuum(90000000);
  theEnvironment=vacuum;
  G4int nQuasmons = theQuasmons.size();
  if (nQuasmons) for (G4int iq=0; iq<nQuasmons; iq++)theQuasmons[iq]->KillQuasmon();
} // End of CleanUp

//Evaporate Residual Nucleus
void G4QEnvironment::EvaporateResidual(G4QHadron* qH,  G4bool fCGS)
{
  static const G4double mAlph = G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double mDeut = G4QPDGCode(2112).GetNuclMass(1,1,0);
  static const G4double mNeut = G4QPDGCode(2112).GetMass();
  static const G4double mProt = G4QPDGCode(2212).GetMass();
  static const G4double mAlPr = mAlph+mProt;
  static const G4double mAlNt = mAlph+mNeut;
  static const G4double dProt = mProt+mProt;
  static const G4double dNeut = mNeut+mNeut;
  static const G4double dAlph = mAlph+mAlph;
  static const G4double eps=.003;
  G4int thePDG = qH->GetPDGCode();           // Get PDG code of the Residual Nucleus
  G4int theBN  = qH->GetBaryonNumber();      // A (Baryon number of the nucleus)
  G4QContent  theQC  = qH->GetQC();          // Quark Content of the hadron
  G4int theS=theQC.GetStrangeness();         // S (Strangeness of the nucleus)
#ifdef debug
  G4cout<<"G4QE::EvaporateRes:Called for PDG="<<thePDG<<",4M="<<qH->Get4Momentum()<<G4endl;
#endif
  if(thePDG==10)
  {
    G4QContent   chQC=qH->GetQC();           // Quark content of the Hadron-Chipolino
    G4QChipolino QCh(chQC);                  // Define a Chipolino instance for the Hadron
    G4LorentzVector ch4M=qH->Get4Momentum(); // 4Mom of the Hadron-Chipolino
    G4QPDGCode h1QPDG=QCh.GetQPDG1();        // QPDG of the first hadron
    G4double   h1M   =h1QPDG.GetMass();      // Mass of the first hadron
    G4QPDGCode h2QPDG=QCh.GetQPDG2();        // QPDG of the second hadron
    G4double   h2M   =h2QPDG.GetMass();      // Mass of the second hadron
    G4double   chM2  =ch4M.m2();             // Squared Mass of the Chipolino
    if( sqr(h1M+h2M) < chM2 )                // Decay is possible
    {
      G4LorentzVector h14M(0.,0.,0.,h1M);
      G4LorentzVector h24M(0.,0.,0.,h2M);
      if(!G4QHadron(ch4M).DecayIn2(h14M,h24M))
      {
        G4ExceptionDescription ed;
        ed << "QChipolino DecIn2 error: CM=" << std::sqrt(chM2) << " -> h1="
           << h1QPDG << "(" << h1M << ") + h2=" << h1QPDG << "(" << h2M
           << ") = " << h1M+h2M << " **Failed**" << G4endl;
        G4Exception("G4QEnvironment::EvaporateResidual()", "HAD_CHPS_0000",
                    FatalException, ed);
      }
      delete qH;                             // Kill the primary Chipolino
      G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);
      theQHadrons.push_back(h1H);            // (delete equivalent)
#ifdef debug
      G4cout<<"G4QE::EvaporateResidual: Chipolino -> H1="<<h1QPDG<<h14M<<G4endl;
#endif
      qH = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
      theQHadrons.push_back(qH);             // (delete equivalent)
#ifdef debug
      G4cout<<"G4QE::EvaporateResidual: Chipolino -> H2="<<h2QPDG<<h24M<<G4endl;
#endif
    }
    else
    {
      // G4cerr<<"***G4QEnv::EvaporateResid: ChipoMeson="<<qH->GetQC()<<qH->Get4Momentum()
      //       <<", chipoM="<<std::sqrt(chM2)<<" < m1="<<h1M<<"("<<h1QPDG<<") + m2="<<h2M
      //       <<"("<<h2QPDG<<") = "<<h1M+h2M<<G4endl;
      // throw G4QException("G4QEnvironment::EvaporateResidual: LowMassChipolino in Input");
      G4ExceptionDescription ed;
      ed << "LowMassChipolino in Input: ChipoMeson=" << qH->GetQC()
         << qH->Get4Momentum() << ", chipoM=" << std::sqrt(chM2) << " < m1="
         << h1M << "(" << h1QPDG << ") + m2=" << h2M << "(" << h2QPDG << ") = "
         << h1M+h2M << G4endl;
      G4Exception("G4QEnvironment::EvaporateResidual()", "HAD_CHPS_0001",
                  FatalException, ed);
    }
    return;
  }
  else if(theS<0)                            // Antistrange nucleus
  {
#ifdef debug
    G4cout<<"G4QE::EvaporateRes: AntistrangeNucleus="<<thePDG<<qH->Get4Momentum()<<G4endl;
#endif
    DecayAntistrange(qH, &theQHadrons);      // (delete equivalent)
    return;
  }
  else if(theBN==1)
  {
#ifdef debug
    G4cout<<"G4QE::EvaporateRes: Baryon="<<thePDG<<qH->Get4Momentum()<<G4endl;
#endif
    DecayBaryon(qH, &theQHadrons);           // (delete equivalent)
    return;
  }
  else if(!theBN) // @@ In future it is usefull to add the MesonExcitationDecay (?!)
  {
#ifdef debug
    G4LorentzVector mesLV=qH->Get4Momentum();
    G4cout<<"G4QE::EvaporateRes:(!)Meson(!) PDG="<<thePDG<<",4M="<<mesLV<<mesLV.m()
          <<",QC="<<qH->GetQC()<<",MPDG="<<G4QPDGCode(thePDG).GetMass()<<G4endl;
#endif
    DecayMeson(qH, &theQHadrons);            // (delete equivalent)
    //theQHadrons.push_back(qH);             // Old solution
    return;
  }
  G4int theC=theQC.GetCharge();              // P
  if(!thePDG) thePDG = theQC.GetSPDGCode();  // If there is no PDG code, get it from QC
  if(thePDG==10 && theBN>0) thePDG=theQC.GetZNSPDGCode();
  if(theS>0) thePDG-=theS*999999;            // @@ May hide hypernuclear problems (G4) ! @@
  G4double totGSM = G4QNucleus(thePDG).GetGSMass();// TheGroundStMass of theTotalResNucleus
  if(theBN==2)
  {
    if(!theC)        totGSM=dNeut;           // nn, nL, LL
    else if(theC==2) totGSM=dProt;           // pp
    else             totGSM=mDeut;           // np, Lp
  }
  else if(theBN==5)
  {
    if     (theC==3) totGSM=mAlPr;           // effective "Alph+p"
    else if(theC==2) totGSM=mAlNt;           // effective "Alph+n"
  }
  else if(theBN==8)   totGSM=dAlph;          // effective "Be8"
  // @@ Should be more (else if) for bigger A=theBN
  G4LorentzVector q4M = qH->Get4Momentum();  // Get 4-momentum of theTotalResidNucleus
  G4double    totMass = q4M.m();             // Get theRealMass of theTotalResidNucleus
  if(fabs(totMass-totGSM)<eps)
  {
    theQHadrons.push_back(qH);               // fill As It Is
  }
  else if(totMass>totGSM)
  {
    theEnvironment.EvaporateNucleus(qH,&theQHadrons);
#ifdef qdebug
    qH=0;
#endif
  }
  else                                       // Correction must be done
  {
#ifdef debug
    G4cout<<"G4QE::EvaRes: *Correct* "<<theQC<<q4M<<totMass<<"<"<<totGSM<<G4endl;
#endif
    G4Quasmon* quasH = new G4Quasmon(theQC,q4M);
    if(fCGS && !CheckGroundState(quasH, true) )
    {
#ifdef debug
      G4cout<<"***G4QE::EvaporResid:GSCorFailed.FillAsItIs,n="<<theQHadrons.size()<<G4endl;
#endif
      theQHadrons.push_back(qH);             // Correction failed: fill as it is
#ifdef qdebug
      qH=0;
#endif
    }
    else
    {
      delete qH;
#ifdef qdebug
      qH=0;
#endif
    }
    delete quasH;
  }
#ifdef qdebug
  if (qH)
  {
    G4cout<<"G4QEnvironment::EvaporateResidual:EndDeleted, PDG="<<qH->GetPDGCode()<<G4endl;
    delete qH;
  }
#endif
  return;
} // End of EvaporateResidual

//Public Hadronisation function with Exception treatment (del is responsibility of User!)
G4QHadronVector* G4QEnvironment::Fragment()
{
#ifdef chdebug
  G4int fCharge=theEnvironment.GetCharge();
  G4int fBaryoN=theEnvironment.GetA();
  G4int nHad=theQHadrons.size();
  if(nHad) for(G4int ih=0; ih<nHad; ih++)
  {
    fCharge+=theQHadrons[ih]->GetCharge();
    fBaryoN+=theQHadrons[ih]->GetBaryonNumber();
  }
  G4int nQuas=theQuasmons.size();
  if(nQuas)for(G4int iqs=0; iqs<nQuas; iqs++)
  {
    fCharge+=theQuasmons[iqs]->GetCharge();
    fBaryoN+=theQuasmons[iqs]->GetBaryonNumber();
  }
  if(fCharge!=totCharge || fBaryoN!=totBaryoN)
  {
    G4cout<<"*::*G4QE::Frag:(4) tC="<<totCharge<<",C="<<fCharge<<",tB="<<totBaryoN
          <<",B="<<fBaryoN<<",E="<<theEnvironment<<G4endl;
    if(nHad) for(G4int h=0; h<nHad; h++)
    {
      G4QHadron* cH = theQHadrons[h];
      G4cout<<"*::*G4QE::HQ:h#"<<h<<",QC="<<cH->GetQC()<<",PDG="<<cH->GetPDGCode()<<G4endl;
    }
    if(nQuas) for(G4int q=0; q<nQuas; q++)
    {
      G4Quasmon* cQ = theQuasmons[q];
      G4cout<<"*::*G4QE::HQ:q#"<<q<<",C="<<cQ->GetCharge()<<",QCont="<<cQ->GetQC()<<G4endl;
    }
  }
#endif
  G4QHadronVector dummy;       // Prototype of the output G4QHadronVector to avoid warnings
  G4QHadronVector* theFragments = &dummy; // Prototype of the output G4QHadronVector
  G4int ExCount =0;                       // Counter of the repetitions
  G4int MaxExCnt=1;                       // A#of of repetitions + 1 (1 for no repetitions)
  G4bool RepFlag=true;                    // To come inside the while
  // For the purpose of the recalculation the Quasmons, Hadrons, Environment must be stored
  G4QuasmonVector* reQuasmons = new G4QuasmonVector; // deleted after the "while LOOP"
  G4int nQ = theQuasmons.size();
  if(nQ)
  {
    for(G4int iq=0; iq<nQ; iq++)
    {
      G4Quasmon* curQ     = new G4Quasmon(theQuasmons[iq]);
      reQuasmons->push_back(curQ);                  // deleted after the "while LOOP"
    }
  }
  G4QHadronVector* reQHadrons = new G4QHadronVector; // deleted after the "while LOOP"
  G4int nH = theQHadrons.size();
  if(nH)
  {
    for(G4int ih=0; ih<nH; ih++)
    {
      G4QHadron* curH     = new G4QHadron(theQHadrons[ih]);
      reQHadrons->push_back(curH);                 // deleted after the "while LOOP"
    }
  }
  G4QNucleus reEnvironment=theEnvironment;
  G4LorentzVector rem4M=tot4Mom;
  while (RepFlag && ExCount<MaxExCnt)
  {
    try
    {
      RepFlag=false;                      // If OK - go out of the while
      theFragments = FSInteraction();     // InterClass creation. User must delet QHadrons.
    }
    catch (G4QException& error)
    {
      G4cout<<"***G4QEnvironment::Fragment: Exception is catched"<<G4endl;
      RepFlag=true;                       // For the Exception - repete
      ExCount++;                          // Increment the repetition counter
      G4cout<<"***G4QEnv::Fragment:Exception #"<<ExCount<<": "<<error.GetMessage()<<G4endl;
      G4LorentzVector dif=rem4M-theEnvironment.Get4Momentum(); // CHECK difference
      G4int nHp=theQHadrons.size();
      G4int nQp = theQuasmons.size();
      G4cout<<"***G4QEnvir::Fragment:nH="<<nHp<<",nQ="<<nQp<<",E="<<theEnvironment<<G4endl;
      for(G4int ph=0; ph<nHp; ph++)
      {
        G4QHadron* cH = theQHadrons[ph];
        dif-=cH->Get4Momentum();
        G4cout<<"***G4QEnvir::Fr:H"<<ph<<"="<<cH->Get4Momentum()<<cH->GetPDGCode()<<G4endl;
      }
      for(G4int pq=0; pq<nQp; pq++)
      {
        G4Quasmon* cQ = theQuasmons[pq];
        dif-=cQ->Get4Momentum();
        G4cout<<"***G4QEnvir::Fr:Quasm#"<<pq<<"="<<cQ->Get4Momentum()<<cQ->GetQC()<<G4endl;
      }
      // *** Cleaning Up of all old output instances for the recalculation purposes ***
      for_each(theFragments->begin(), theFragments->end(), DeleteQHadron()); // old Hadrons
      theFragments->clear();
      for_each(theQHadrons.begin(), theQHadrons.end(), DeleteQHadron()); //internal Hadrons
      theQHadrons.clear();
      for_each(theQuasmons.begin(), theQuasmons.end(), DeleteQuasmon()); // old Quasmons
      theQuasmons.clear();
      G4cout<<"***G4QEnv::Fragment: ----------- End of CleaningUp: 4Mdif="<<dif<<G4endl;
      // **************** Recover all conditions for the recalculation ********************
      theEnvironment=reEnvironment;             // Recover the nuclear environment
      tot4Mom=rem4M;                            // Recover the total 4Momentum of the React
      G4cout<<"***G4QEnv::Fragment:*Recover*Env="<<theEnvironment<<",4M="<<tot4Mom<<G4endl;
      G4int mQ = reQuasmons->size();            // Recover the memorizedQuasmons with print
      for(G4int jq=0; jq<mQ; jq++)
      {
        //G4Quasmon* curQ = new G4Quasmon(reQuasmons->operator[](jq));
        G4Quasmon* curQ = new G4Quasmon((*reQuasmons)[jq]);
        G4cout<<"***G4QE::Fragm:Q("<<jq<<")="<<curQ->Get4Momentum()<<curQ->GetQC()<<G4endl;
        theQuasmons.push_back(curQ);            // (delete equivalent)
      }
      G4int mH = reQHadrons->size();            // Recover the memorizedQHadrons with print
      for(G4int jh=0; jh<mH; jh++)
      {
        //G4QHadron* curH = new G4QHadron(reQHadrons->operator[](jh));
        G4QHadron* curH = new G4QHadron((*reQHadrons)[jh]);
        G4cout<<"***G4QE::Fragm:H("<<jh<<")="<<curH->Get4Momentum()<<curH->GetQC()<<G4endl;
        theQHadrons.push_back(curH);            // (delete equivalent)
      }
    }
  }
  if(reQuasmons->size()) // If something is still in memory then clean it up
  {
    for_each(reQuasmons->begin(),reQuasmons->end(),DeleteQuasmon()); // CleanUp oldQuasmons
    reQuasmons->clear();
  }
  delete reQuasmons;     // All temporary Quasmons memory is wiped out
  if(reQHadrons->size()) // If something is still in memory then clean it up
  {
    for_each(reQHadrons->begin(),reQHadrons->end(),DeleteQHadron()); //CleanUp old QHadrons
    reQHadrons->clear();
  }
  delete reQHadrons;     // All temporary QHadrons memory is wiped out
  if(ExCount>=MaxExCnt)
  {
    G4int nProj=theProjectiles.size();
    G4cerr<<"*G4QEnv::Fragment:Exception.Target="<<theTargetPDG<<". #Proj="<<nProj<<G4endl;
    if(nProj) for(G4int ipr=0; ipr<nProj; ipr++)
    {
      G4QHadron* prH = theProjectiles[ipr];
      G4cerr<<"G4QE::F:#"<<ipr<<",PDG/4M="<<prH->GetPDGCode()<<prH->Get4Momentum()<<G4endl;
    }
    G4Exception("G4QEnvironment::Fragment()", "HAD_CHPS_0000",
                FatalException, "This reaction caused the CHIPSException");
  }
  // Put the postponed hadrons in the begining of theFragments and clean them up
  G4int tmpS=intQHadrons.size();
  if(tmpS)
  {
    //tmpS=theFragments->size();
    //intQHadrons.resize(tmpS+intQHadrons.size());
    //copy(theFragments->begin(), theFragments->end(), intQHadrons.end()-tmpS);
    //tmpS=intQHadrons.size();
    //theFragments->resize(tmpS);  // Resize theFragments
    //copy(intQHadrons.begin(), intQHadrons.end(), theFragments->begin());
    // Can be like this, but by performance it is closer to FNAL (but better than it)
    //copy(theFragments->begin(), theFragments->end(), back_inserter(intQHadrons));
    //theFragments->resize(intQHadrons.size());  // Resize theFragments
    //copy(intQHadrons.begin(), intQHadrons.end(), theFragments->begin());
    // The following (istead of all above) has worse performance !
    theFragments->insert(theFragments->begin(), intQHadrons.begin(), intQHadrons.end() );
    intQHadrons.clear();
  }
  // No we need to check that all hadrons are on the mass shell
  CheckMassShell(theFragments); // @@ the same can be done in the end of G4QFragmentation
  return theFragments;
} // End of the Fragmentation member function

//The Final State Interaction Filter for the resulting output of ::HadronizeQEnvironment()
G4QHadronVector* G4QEnvironment::FSInteraction()
{
  static const G4QPDGCode gQPDG(22);
  static const G4QPDGCode pizQPDG(111);
  static const G4QPDGCode pipQPDG(211);
  static const G4QPDGCode pimQPDG(-211);
  static const G4QPDGCode nQPDG(2112);
  static const G4QPDGCode pQPDG(2212);
  static const G4QPDGCode lQPDG(3122);
  static const G4QPDGCode s0QPDG(3212);
  //static const G4QPDGCode dQPDG(90001001);
  static const G4QPDGCode tQPDG(90001002);
  static const G4QPDGCode he3QPDG(90002001);
  static const G4QPDGCode aQPDG(90002002);
  static const G4QPDGCode a6QPDG(90002004);
  static const G4QPDGCode be6QPDG(90004002);
  //static const G4QPDGCode b7QPDG(90005002);
  //static const G4QPDGCode he7QPDG(90002005);
  static const G4QPDGCode c8QPDG(90006002);
  static const G4QPDGCode a8QPDG(90002006);
  static const G4QPDGCode c10QPDG(90006004);
  static const G4QPDGCode o14QPDG(90008006);
  static const G4QPDGCode o15QPDG(90008007);
  static const G4QContent K0QC(1,0,0,0,0,1);
  static const G4QContent KpQC(0,1,0,0,0,1);
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mPi0 = G4QPDGCode(111).GetMass();
  static const G4double mK   = G4QPDGCode(321).GetMass();
  static const G4double mK0  = G4QPDGCode(311).GetMass();
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mSigZ= G4QPDGCode(3212).GetMass();
  static const G4double mSigM= G4QPDGCode(3112).GetMass();
  static const G4double mSigP= G4QPDGCode(3222).GetMass();
  static const G4double mXiZ = G4QPDGCode(3322).GetMass();
  static const G4double mXiM = G4QPDGCode(3312).GetMass();
  static const G4double mOmM = G4QPDGCode(3334).GetMass();
  static const G4double mDeut= G4QPDGCode(2112).GetNuclMass(1,1,0);
  static const G4double mTrit= G4QPDGCode(2112).GetNuclMass(1,2,0);
  static const G4double mHe3 = G4QPDGCode(2112).GetNuclMass(2,1,0);
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double mHe6 = G4QPDGCode(2112).GetNuclMass(2,4,0);
  static const G4double mBe6 = G4QPDGCode(2112).GetNuclMass(4,2,0);
  //static const G4double mHe7 = G4QPDGCode(2112).GetNuclMass(2,5,0);
  //static const G4double mB7  = G4QPDGCode(2112).GetNuclMass(5,2,0);
  static const G4double mHe8 = G4QPDGCode(2112).GetNuclMass(2,6,0);
  static const G4double mC8  = G4QPDGCode(2112).GetNuclMass(6,2,0);
  static const G4double mC10 = G4QPDGCode(2112).GetNuclMass(6,4,0);
  static const G4double mO14 = G4QPDGCode(2112).GetNuclMass(8,6,0);
  static const G4double mO15 = G4QPDGCode(2112).GetNuclMass(8,7,0);
  static const G4double mKmP = mK+mProt;
  static const G4double mKmN = mK+mNeut;
  static const G4double mK0mP = mK0+mProt;
  static const G4double mK0mN = mK0+mNeut;
  static const G4QNucleus vacuum(90000000);
  static const G4double eps=0.003;
  ///////////////static const G4double third=1./3.;
  ///////////////static const G4double nPDG=90000001;
  G4int envA=theEnvironment.GetBaryonNumber();
  ///////////////G4int envC=theEnvironment.GetCharge();
#ifdef rdebug
  G4int totInC=theEnvironment.GetZ();
  G4LorentzVector totIn4M=theEnvironment.Get4Momentum();
  G4cout<<"G4QEnvironment(G4QE)::FSInter(FSI): ***called*** envA="<<envA<<totIn4M<<G4endl;
  G4int nQuasmons=theQuasmons.size();
  for (G4int is=0; is<nQuasmons; is++) // Sum 4mom's of Quasmons for comparison
  {
    G4Quasmon*      pQ = theQuasmons[is];
    G4LorentzVector Q4M= pQ->Get4Momentum();
    G4cout<<"G4QE::FSI: Quasmon ("<<is<<") is added, 4M="<<Q4M<<G4endl;
    totIn4M           += Q4M;
    totInC            += pQ->GetQC().GetCharge();
  } // End of TotInitial4Momentum summation LOOP over Quasmons
  G4int nsHadr  = theQHadrons.size();        // Update the value of OUTPUT entries
  if(nsHadr) for(G4int jso=0; jso<nsHadr; jso++)// LOOP over output hadrons 
  {
    G4int hsNF  = theQHadrons[jso]->GetNFragments(); // A#of secondary fragments
    if(!hsNF)                                        // Add only final hadrons
    {
      G4LorentzVector hs4Mom = theQHadrons[jso]->Get4Momentum();
      G4cout<<"G4QE::FSI: Hadron ("<<jso<<") is added, 4M="<<hs4Mom<<G4endl;
      totIn4M          += hs4Mom;
      totInC           += theQHadrons[jso]->GetCharge();
    }
  }
  G4cout<<"G4QE::FSI: The resulting 4Momentum="<<totIn4M<<G4endl;
#endif
#ifdef chdebug
  G4int fCharge=theEnvironment.GetCharge();
  G4int fBaryoN=theEnvironment.GetA();
  G4int nHad=theQHadrons.size();
  if(nHad) for(G4int ih=0; ih<nHad; ih++)
  {
    fCharge+=theQHadrons[ih]->GetCharge();
    fBaryoN+=theQHadrons[ih]->GetBaryonNumber();
  }
  G4int nQuas=theQuasmons.size();
  if(nQuas)for(G4int iqs=0; iqs<nQuas; iqs++)
  {
    fCharge+=theQuasmons[iqs]->GetCharge();
    fBaryoN+=theQuasmons[iqs]->GetBaryonNumber();
  }
  if(fCharge!=totCharge || fBaryoN!=totBaryoN)
  {
    G4cout<<"*::*G4QE::FSI:(5) tC="<<totCharge<<",C="<<fCharge<<",tB="<<totBaryoN
          <<",B="<<fBaryoN<<",E="<<theEnvironment<<G4endl;
    if(nHad) for(G4int h=0; h<nHad; h++)
    {
      G4QHadron* cH = theQHadrons[h];
      G4cout<<"*::*G4QE::HQ:h#"<<h<<",QC="<<cH->GetQC()<<",PDG="<<cH->GetPDGCode()<<G4endl;
    }
    if(nQuas) for(G4int q=0; q<nQuas; q++)
    {
      G4Quasmon* cQ = theQuasmons[q];
      G4cout<<"*::*G4QE::HQ:q#"<<q<<",C="<<cQ->GetCharge()<<",QCont="<<cQ->GetQC()<<G4endl;
    }
  }
#endif
  G4QHadronVector* theFragments = new G4QHadronVector;//Internal creation. User must delete
  HadronizeQEnvironment();                // --->> Call the main fragmentation function
#ifdef rdebug
  G4int tC=totInC-theEnvironment.GetZ();  // Subtract theResidualEnvironCharge 
  G4LorentzVector t4M=totIn4M;            // Compare with the total                 
  G4cout<<"G4QEnv::FSI: Initial tot4M="<<t4M<<" to be subtracted"<<G4endl;
  G4LorentzVector theEnv4m=theEnvironment.Get4Momentum(); // Environment 4Mom 
  t4M-=theEnv4m;                          // Subtract the Environment 4-momentum    
  G4cout<<"G4QEnv::FSI: Subtract Environ="<<theEnv4m<<theEnvironment<<G4endl;
  for (G4int js=0; js<nQuasmons; js++)    // Subtract 4mom's of Quasmons (compare)
  {                                       //                                        
    G4Quasmon*      prQ = theQuasmons[js];//                                  
    if(prQ->GetStatus())                  // Subtract only if Quasmon is alive      
    {                                     //                                        
      G4LorentzVector Q4M= prQ->Get4Momentum(); // 4-momentum of the Quasmon
      G4QContent qQC= prQ->GetQC();       //                                        
      G4cout<<"G4QE::FSI: Subtract Quasmon("<<js<<"),4M="<<Q4M<<",QC="<<qQC<<G4endl;
      t4M          -= Q4M;                // Subtract 4-momentum of the Quasmon        
      tC           -= prQ->GetQC().GetCharge(); //                            
    }                                     //                                        
    else G4cout<<"G4QE::FSI:Dead Quasmon("<<js<<")="<<prQ->GetStatus()<<G4endl;
  } // End of Quasmons4Momentum subtractions                                  
  G4int nsbHadr=theQHadrons.size();       // Update the value of OUTPUT entries     
  if(nsbHadr) for(G4int jpo=0; jpo<nsbHadr; jpo++)// LOOP over output hadrons 
  {
    G4int hsNF  = theQHadrons[jpo]->GetNFragments();// A#of out fragments     
    if(!hsNF)                            // Subtract only final hadrons            
    {                                    //                                        
      G4LorentzVector hs4Mom = theQHadrons[jpo]->Get4Momentum(); // Output hadron
      G4int hPDG = theQHadrons[jpo]->GetPDGCode(); // PDG of the Output Hadron
      G4cout<<"G4QE::FSI: Subtract Hadron("<<jpo<<"), 4M="<<hs4Mom<<hPDG<<G4endl; 
      t4M          -= hs4Mom;           //                                        
      tC           -= theQHadrons[jpo]->GetCharge(); // Subtract charge of the OutHadron
    }                                   // End of the "FinalHadron" IF            
  }                                     // End of the LOOP over output hadrons    
  G4cout<<"G4QEnv::FSI:|||||4-MomCHECK||||d4M="<<t4M<<",dCharge="<<tC<<G4endl;
#endif
  unsigned nHadr=theQHadrons.size();
  if(nHadr<=0)
  {
    G4cout<<"---Warning---G4QEnvironment::FSInteraction: nHadrons="<<nHadr<<G4endl;
    //throw G4QException("G4QEnvironment::FSInteraction: No hadrons in the output");
    return theFragments;
  }
  G4int lHadr=theQHadrons[nHadr-1]->GetBaryonNumber();
#ifdef debug
  G4cout<<"G4QE::FSI:after HQE,nH="<<nHadr<<",lHBN="<<lHadr<<",E="<<theEnvironment<<G4endl;
#endif
  if(lHadr>1)                          // TheLastHadron is nucleus:try to decay/evap/cor it
  {
    G4QHadron* theLast = theQHadrons[nHadr-1];
    G4QHadron* curHadr = new G4QHadron(theLast);
    G4LorentzVector lh4M=curHadr->Get4Momentum(); // Actual mass of the last fragment
    G4double lhM=lh4M.m();             // Actual mass of the last fragment
    G4int lhPDG=curHadr->GetPDGCode();            // PDG code of the last fragment
    G4double lhGSM=G4QPDGCode(lhPDG).GetMass();   // GroundStateMass of the last fragment
#ifdef debug
    G4cout<<"G4QE::FSI:lastHadr 4M/M="<<lh4M<<lhM<<",GSM="<<lhGSM<<",PDG="<<lhPDG<<G4endl;
#endif
    if(lhM>lhGSM+eps)                  // ==> Try to evaporate the residual nucleus
    {
      theQHadrons.pop_back();          // the last QHadron-Nucleus is excluded from OUTPUT
      delete theLast;// *!!When kill,DON'T forget to delete theLastQHadron as an instance!*
      EvaporateResidual(curHadr);      // Try to evaporate Hadr-Nucl (@@DecDib)(delete eq.)
      nHadr=theQHadrons.size();
#ifdef debug
      G4cout<<"G4QE::FSI:After nH="<<nHadr<<",PDG="<<curHadr->GetPDGCode()<<G4endl;
#endif
    }
    else if(lhM<lhGSM-eps)             // ==> Try to make the HadronicSteck FSI correction
    {
      theQHadrons.pop_back();          //exclude LastHadronPointer from OUTPUT
      delete theLast;      // *!! When killing, DON'T forget to delete the last QHadron !!*
      G4Quasmon* quasH = new G4Quasmon(curHadr->GetQC(),lh4M); // Fake Quasmon ctreation
      if(!CheckGroundState(quasH,true))// Try to correct with other hadrons
      {
#ifdef debug
        // M.K. Fake complain in the low energy nHe/pHe reactions, while everything is OK
        G4cout<<"---Warning---G4QEnv::FSI:Correction error LeaveAsItIs h4m="<<lh4M<<G4endl;
#endif
        theQHadrons.push_back(curHadr);// Fill theResidualNucleus asItIs(delete equivalent)
        //throw G4QException("G4QEnv::FSI: STOP at Correction Error");
      }
      else
      {
        delete curHadr;              // The intermediate curHadr isn't necessary any more
        nHadr=theQHadrons.size();      // Update nHadr after successful correction
      }
      delete quasH;                    // Delete the temporary fake Quasmon
    }
    else delete curHadr;               // ==> Leave the nucleus as it is (close to the GSM)
  }
#ifdef debug
  G4LorentzVector ccs4M(0.,0.,0.,0.);  // CurrentControlSum of outgoing Hadrons
#endif
  // *** Initial Charge Control Sum Calculation ***
  G4int chContSum=0;                   // ChargeControlSum to keepTrack FSI transformations
  G4int bnContSum=0;                   // BaryoNControlSum to keepTrack FSI transformations
  if(nHadr)for(unsigned ich=0; ich<nHadr; ich++) if(!(theQHadrons[ich]->GetNFragments()))
  {
    chContSum+=theQHadrons[ich]->GetCharge();
    bnContSum+=theQHadrons[ich]->GetBaryonNumber();
  }
#ifdef chdebug
  if(chContSum!=totCharge || bnContSum!=totBaryoN)
  {
    G4cout<<"*::*G4QE::Fr:(6)tC="<<totCharge<<",C="<<chContSum<<",tB="<<totBaryoN
          <<",B="<<bnContSum<<",E="<<theEnvironment<<G4endl;
    if(nHadr) for(unsigned h=0; h<nHadr; h++)
    {
      G4QHadron* cH = theQHadrons[h];
      G4cout<<"*::*G4QE::HQ:h#"<<h<<",QC="<<cH->GetQC()<<",PDG="<<cH->GetPDGCode()<<G4endl;
    }
    if(nQuas) for(G4int q=0; q<nQuas; q++)
    {
      G4Quasmon* cQ = theQuasmons[q];
      G4cout<<"*::*G4QE::HQ:q#"<<q<<",C="<<cQ->GetCharge()<<",QCont="<<cQ->GetQC()<<G4endl;
    }
  }
#endif
  // ***
  if(nHadr)for(unsigned ipo=0; ipo<theQHadrons.size(); ipo++)//FindBigestNuclFragm & DecayA
  {
    unsigned jpo=ipo;
    nHadr=theQHadrons.size();
    lHadr=theQHadrons[nHadr-1]->GetBaryonNumber();
    G4QHadron* theCurr = theQHadrons[ipo];    // Pointer to the Current Hadron
    G4int hBN  = theCurr->GetBaryonNumber();
    G4int sBN  = theCurr->GetStrangeness();
    G4int cBN  = theCurr->GetCharge();
    G4int hPDG = theCurr->GetPDGCode();
    G4LorentzVector h4Mom = theCurr->Get4Momentum();
#ifdef debug
    G4int hNF  = theCurr->GetNFragments();
    G4cout<<"G4QE::FSI:h#"<<ipo<<",PDG="<<hPDG<<h4Mom<<",mGS="<<G4QPDGCode(hPDG).GetMass()
          <<",F="<<hNF<<",nH="<<theQHadrons.size()<<G4endl;
#endif
    if(hBN>lHadr && ipo+1<theQHadrons.size()) // CurrHadron=BiggestFragment -> Swap w/ Last
    {
      G4QHadron* curHadr = new G4QHadron(theCurr);// Remember CurHadron for evaporation
      G4QHadron* theLast = theQHadrons[theQHadrons.size()-1]; // theLastHadron (Cur<Last)
      G4QPDGCode lQP=theLast->GetQPDG();      // The QPDG of the last
      if(lQP.GetPDGCode()!=10) theCurr->SetQPDG(lQP); //CurHadr instead of LastHadr
      else theCurr->SetQC(theLast->GetQC());  // CurHadrPDG instead of LastHadrPDG
      theCurr->Set4Momentum(theLast->Get4Momentum()); // ... continue substitution
      h4Mom = theCurr->Get4Momentum();
      hBN  = theCurr->GetBaryonNumber();
      cBN  = theCurr->GetCharge();
      sBN  = theCurr->GetStrangeness();
      hPDG = theCurr->GetPDGCode();
      theQHadrons.pop_back();         // pointer to theLastHadron is excluded from OUTPUT
      delete theLast;// *!!When kill,DON'T forget to delete theLastQHadron asAnInstance !!*
      theQHadrons.push_back(curHadr);
      nHadr=theQHadrons.size();
    }
    if(hPDG==89002000||hPDG==89001001||hPDG==89000002)// 2pt dec. of anti-strange (3pt dec)
    {
#ifdef debug
      G4cout<<"G4QE::FSI:***ANTISTRANGE*** i="<<ipo<<",PDG="<<hPDG<<",BaryN="<<hBN<<G4endl;
#endif
      G4double hM=h4Mom.m(); // 89002000
      G4double hMi=hM+eps;
      G4QPDGCode fQPDG = pQPDG;
      G4double fM = mProt;
      G4int  sPDG = 321;
      G4double sM = mK;
      G4int  tPDG = 0;
      G4double tM = 0.;
      if(hPDG==89002000)                     // Use the prototypes above
      {
        if(hMi<mKmP)
        {
          if(hMi>mProt+mPi+mPi0)
          {
            sPDG=211;
            sM  =mPi;
            tPDG=111;
            tM  =mPi0;
          }
          else if(hMi>mProt+mPi) // @@ Does not conserve strangeness (Week decay)
          {
#ifdef debug
            G4cout<<"**G4QE::FSI:ANTISTRANGE*++*STRANGENESS,PDG="<<hPDG<<",M="<<hM<<G4endl;
#endif
            sPDG=211;
            sM  =mPi;
          }
          else sPDG=0;
        }
      }
      else if(hPDG==89001001)
      {
        fQPDG= nQPDG;
        fM   = mNeut;
        sPDG = 321;
        sM   = mK;
        if(hMi>mK0mP&&G4UniformRand()>.5)
        {
          fQPDG= pQPDG;
          fM   = mProt;
          sPDG = 311;
          sM   = mK0;
        }
        else if(hMi<mKmN)
        {
          if(hMi>mProt+mPi0+mPi0)
          {
            fQPDG= pQPDG;
            fM   = mProt;
            sPDG = 111;
            sM   = mPi0;
            tPDG = 111;
            tM   = mPi0;
            if(hMi>mNeut+mPi+mPi0&&G4UniformRand()>.67)
            {
              fQPDG= nQPDG;
              fM   = mNeut;
              tPDG = 211;
              tM   = mPi;
            }
            if(hMi>mProt+mPi+mPi&&G4UniformRand()>.5)
            {
              sPDG = 211;
              sM   = mPi;
              tPDG =-211;
              tM   = mPi;
            }
          }
          else if(hMi>mProt+mPi0) // @@ Does not conserve strangeness (Week decay)
          {
#ifdef debug
            G4cout<<"**G4QE::FSI:*ANTISTRANGE*+*STRANGENESS*PDG="<<hPDG<<",M="<<hM<<G4endl;
#endif
            fQPDG= pQPDG;
            fM   = mProt;
            sPDG = 111;
            sM   = mPi0;
          }
          else sPDG=0;      // @@ Still can try to decay in gamma+neutron (electromagnetic)
        }
      }
      else if(hPDG==89000002)
      {
        fQPDG= nQPDG;
        fM   = mNeut;
        sPDG = 311;
        sM   = mK0;
        if(hMi<mK0mN)
        {
          if(hMi>mNeut+mPi+mPi)
          {
            sPDG = 211;
            sM   = mPi;
            tPDG =-211;
            tM   = mPi;
          }
          if(hMi>mProt+mPi+mPi0)
          {
            fQPDG= pQPDG;
            fM   = mProt;
            sPDG = 111;
            sM   = mPi0;
            tPDG =-211;
            tM   = mPi;
          }
          else if(hMi>mProt+mPi) // @@ Does not conserve strangeness (Week decay)
          {
#ifdef debug
            G4cout<<"**G4QE::FSI:**ANTISTRANGE*0*STRANGENE**PDG="<<hPDG<<",M="<<hM<<G4endl;
#endif
            fQPDG= pQPDG;
            fM   = mProt;
            sPDG =-211;
            sM   = mPi;
          }
          else sPDG=0;      // @@ Still can try to decay in gamma+neutron (electromagnetic)
        }
      }
      if(!sPDG)
      {
#ifdef debug
        G4cout<<"***G4QE::FSI:***ANTISTRANGE***CANN'T DECAY,PDG="<<hPDG<<",M="<<hM<<G4endl;
#endif
      }
      else if(!tPDG)           // 2 particle decay
      {
        G4bool fOK=true;
        G4LorentzVector f4M(0.,0.,0.,fM);
        G4LorentzVector s4M(0.,0.,0.,sM);
        G4double sum=fM+sM;
        if(fabs(hM-sum)<=eps)
        {
          f4M=h4Mom*(fM/sum);
          s4M=h4Mom*(sM/sum);
        }
        else if(hM<sum || !G4QHadron(h4Mom).DecayIn2(f4M,s4M))
        {
          G4cout<<"---Warning---G4QE::FSI: Still try(2),M="<<hM<<"->"<<fM<<"("<<fQPDG<<")+"
                <<sM<<"("<<sPDG<<")="<<sum<<G4endl;
          // Tipical scenario of recovery:
          //             1. Check that the Environment is vacuum (must be), 
          //if(theEnvironment==vacuum)
          if(!theEnvironment.GetA())
          {
            //           2. Extract and put in qH, substitute by the Last and make quasH,
            G4QHadron* theLast = theCurr;       // Prototype of thePointer to theLastHadron
            G4QHadron* qH = new G4QHadron(theCurr); // Copy of the Current Hadron
            if(ipo+1<theQHadrons.size())            // If ipo<Last, swap CurHadr & LastHadr
            {
              theLast = theQHadrons[theQHadrons.size()-1];//PointerTo theLastHadr(ipo<Last)
              G4QPDGCode lQP=theLast->GetQPDG();    // The QPDG of the last
              if(lQP.GetPDGCode()!=10) theCurr->SetQPDG(lQP); //CurHadr instead of LastHadr
              else theCurr->SetQC(theLast->GetQC());// CurHadrPDG instead of LastHadrPDG
              theCurr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
            }                                       //ELSE: it's already theLast -> no swap
            theQHadrons.pop_back();                 //exclude LastHadronPointer from OUTPUT
            delete theLast;// *!! When killing, DON'T forget to delete the last QHadron !!*
            G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());//Fake Quasmon
            //           3. Try to use other hadrons to recover this one (under Mass Shell)
            if(!CheckGroundState(quasH,true))       // Try to correct with other hadrons
            {
              G4cout<<"---Warning---G4QEnv::FSI:Failed (2) LeaveAsItIs 4m="<<h4Mom<<G4endl;
              theQHadrons.push_back(qH);            // Fill as it is (delete equivalent)
            }
            else
            {
              delete qH;
              //         4. Decrement jpo index (copy of ipo) to find theBiggestNuclearFrag
              jpo--;
              //         5. Recheck the probuct, which replaced the Last and check others
              nHadr=theQHadrons.size();
            }
            //           6. Delete the intermediate quasmon  
            delete quasH;
            //           7. Forbid the decay
            fOK=false;
          }
          else
          {
            // G4cerr<<"***G4QEnv::FSI: No recovery (1) Env="<<theEnvironment<<G4endl;
            // throw G4QException("G4QEnv::FSI:ANTISTRANGE DecayIn2 did not succeed");
            G4ExceptionDescription ed;
            ed << "FSI:ANTISTRANGE DecayIn2 did not succeed: FSI: No recovery (1) Env="
               << theEnvironment << G4endl;
            G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0000",
                        FatalException, ed);
          }
        }
        if(fOK)
        {
          theQHadrons[ipo]->SetQPDG(fQPDG);
          theQHadrons[ipo]->Set4Momentum(f4M);
          G4QHadron* sH = new G4QHadron(sPDG,s4M);
          theQHadrons.push_back(sH);               // (delete equivalent)
        }
      }
      else
      {
        G4bool fOK=true;
        G4LorentzVector f4M(0.,0.,0.,fM);
        G4LorentzVector s4M(0.,0.,0.,sM);
        G4LorentzVector t4M(0.,0.,0.,tM);
        G4double sum=fM+sM+tM;
        if(fabs(hM-sum)<=eps)
        {
          f4M=h4Mom*(fM/sum);
          s4M=h4Mom*(sM/sum);
          t4M=h4Mom*(tM/sum);
        }
        else if(hM<sum || !G4QHadron(h4Mom).DecayIn3(f4M,s4M,t4M))
        {
          G4cout<<"---Warning---G4QE::FSI: Still try(3), M"<<hM<<"->"<<fM<<"("<<fQPDG<<")+"
                <<sM<<"("<<sPDG<<")+"<<tM<<"("<<tPDG<<")="<<sum<<G4endl;
          //if(theEnvironment==vacuum)
          if(!theEnvironment.GetA())
          {
            G4QHadron* theLast = theCurr;    // Prototype of the pointer to the Last Hadron
            G4QHadron* qH = new G4QHadron(theCurr); // Copy of the Current Hadron
            if(ipo+1<theQHadrons.size())       // If ipo<Last, swap CurHadr and theLastHadr
            {
              theLast = theQHadrons[theQHadrons.size()-1];//Pointer to LastHadron(ipo<Last)
              G4QPDGCode lQP=theLast->GetQPDG();    // The QPDG of the last
              if(lQP.GetPDGCode()!=10) theCurr->SetQPDG(lQP); //CurHadr instead of LastHadr
              else theCurr->SetQC(theLast->GetQC());// CurHadrPDG instead of LastHadrPDG
              theCurr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
            }
            theQHadrons.pop_back();        // exclude theLastHadron pointer from the OUTPUT
            delete theLast;// *!! When killing, DON'T forget to delete the last QHadron !!*
            G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());//Fake Quasmon
            if(!CheckGroundState(quasH,true))         // Try to correct by other hadrons
            {
              G4cout<<"---Warning---G4QEnv::FSI:Failed (3) LeaveAsItIs,4M="<<h4Mom<<G4endl;
              theQHadrons.push_back(qH);              // Fill as it is (delete equivalent)
            }
            else
            {
              delete qH;
              jpo--;
              nHadr=theQHadrons.size();
            }
            delete quasH;
            fOK=false;
          }
          else
          {
            // G4cout<<"***G4QEnv::FSI: No recovery (2) Env="<<theEnvironment<<G4endl;
            // throw G4QException("G4QEnv::FSI:ANTISTRANGE DecayIn3 did not succeed");
            G4ExceptionDescription ed;
            ed << "ANTISTRANGE DecayIn3 did not succeed: No recovery (2) Env="
               << theEnvironment << G4endl;
            G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0001",
                        FatalException, ed);
          }
        }
        if(fOK)
        {
          theQHadrons[ipo]->SetQPDG(fQPDG);
          theQHadrons[ipo]->Set4Momentum(f4M);
          G4QHadron* sH = new G4QHadron(sPDG,s4M);
          theQHadrons.push_back(sH);               // (delete equivalent)
          G4QHadron* tH = new G4QHadron(tPDG,t4M);
          theQHadrons.push_back(tH);               // (delete equivalent)
        }
      }
      nHadr=theQHadrons.size();
    }
    else if(hPDG==89999003||hPDG==90002999||hPDG==90000003||hPDG==90003000||
            hPDG==90999002||hPDG==91001999) // "3-particles decays of dibaryons and 3N"
    {
#ifdef debug
      G4cout<<"G4QE::FSI:***nD-/pD++/nnn/ppp***i="<<ipo<<",PDG="<<hPDG<<",A="<<hBN<<G4endl;
#endif
      G4double hM=h4Mom.m();
      G4QPDGCode nuQPDG=nQPDG; // n,n,pi-
      G4double nucM = mNeut;
      G4int  barPDG = 2112;
      G4double barM = mNeut;
      G4int   tPDG  = -211;
      G4double tM   = mPi;
      if(hPDG==90002999)       // p,p,pi+
      {
        nuQPDG = pQPDG;        // Substitute p for the first n
        nucM   = mProt;
        barPDG = 2212;         // Substitute p for the second n
        barM   = mProt;
        tPDG   = 211;          // Substitute pi+ for the first pi-
      }
      else if(hPDG==90003000)  // 3p
      {
        nuQPDG = pQPDG;        // Substitute p for the first n
        nucM   = mProt;
        barPDG = 2212;         // Substitute p for the second n
        barM   = mProt;
        tPDG   = 2212;         // Substitute p for pi-
        tM     = mProt;
      }
      else if(hPDG==90999002)  // n,Lambda,pi-/n,Sigma0,pi-/n,Sigma-,gamma(@@)
      {
        if(hM>mSigZ+mNeut+mPi)
        {
          G4double ddMass=hM-mPi;          // Free CM energy
          G4double dd2=ddMass*ddMass;      // Squared free energy
          G4double sma=mLamb+mNeut;        // Neutron+Lambda sum
          G4double pr1=0.;                 // Prototype to avoid sqrt(-)
          if(ddMass>sma) pr1=sqrt((dd2-sma*sma)*dd2); // Neut+Lamb PS
          sma=mNeut+mSigZ;                 // Neutron+Sigma0 sum
          G4double smi=mSigZ-mNeut;        // Sigma0-Neutron difference
          G4double pr2=pr1;                // Prototype of +N+S0 PS
          if(ddMass>sma && ddMass>smi) pr2+=sqrt((dd2-sma*sma)*(dd2-smi*smi));
          if(pr2*G4UniformRand()>pr1)      // --> "ENnv+Sigma0+Lambda" case
          {
            barPDG = 3212;     // Substitute Sigma0 for the second n
            barM   = mSigZ;
          }
          else
          {
            barPDG = 3122;     // Substitute Lambda for the second n
            barM   = mLamb;
          }
        }
        else if(hM>mLamb+mNeut+mPi)
        {
          barPDG = 3122;       // Substitute Lambda for the second n
          barM   = mLamb;
        }
        else if(hM>mSigM+mNeut)// @@ Decay in 2
        {
          barPDG = 3112;       // Substitute Sigma- for the second n
          barM   = mSigM;
          tPDG   = 22;
          tM     = 0;
        }
      }
      else if(hPDG==91001999)  // p,Lambda,pi+/p,Sigma0,pi+/p,Sigma+,gamma(@@)
      {
        nuQPDG = pQPDG;        // Substitute p for the first n
        nucM   = mProt;
        tPDG   = 211;          // Substitute pi+ for the first pi-
        if(hM>mSigZ+mProt+mPi)
        {
          G4double ddMass=hM-mPi;          // Free CM energy
          G4double dd2=ddMass*ddMass;      // Squared free energy
          G4double sma=mLamb+mProt;        // Lambda+Proton sum
          G4double pr1=0.;                 // Prototype to avoid sqrt(-)
          if(ddMass>sma) pr1=sqrt((dd2-sma*sma)*dd2); // Lamb+Prot PS
          sma=mProt+mSigZ;                 // Proton+Sigma0 sum
          G4double smi=mSigZ-mProt;        // Sigma0-Proton difference
          G4double pr2=pr1;                // Prototype of +P+S0 PS
          if(ddMass>sma && ddMass>smi) pr2+=sqrt((dd2-sma*sma)*(dd2-smi*smi));
          if(pr2*G4UniformRand()>pr1)      // --> "ENnv+Sigma0+Lambda" case
          {
            barPDG = 3212;     // Substitute Sigma0 for the second n
            barM   = mSigZ;
          }
          else
          {
            barPDG = 3122;     // Substitute Lambda for the second n
            barM   = mLamb;
          }
        }
        if(hM>mLamb+mProt+mPi)
        {
          barPDG = 3122;         // Substitute Lambda for the second n
          barM   = mLamb;
        }
        else if(hM>mSigP+mProt)  // @@ Decay in 2
        {
          barPDG = 3222;         // Substitute Sigma- for the second n
          barM   = mSigP;
          tPDG   = 22;
          tM     = 0;
        }
      }
      else if(hPDG==90000003)  // 3n
      {
        tPDG   = 2112;         // Substitute n for pi-
        tM     = mNeut;
      }
      G4bool fOK=true;
      G4LorentzVector nu4M(0.,0.,0.,nucM);
      G4LorentzVector ba4M(0.,0.,0.,barM);
      G4LorentzVector pi4M(0.,0.,0.,tM);
      G4double sum=nucM+barM+tM;
      if(fabs(hM-sum)<=eps)
      {
        nu4M=h4Mom*(nucM/sum);
        ba4M=h4Mom*(barM/sum);
        pi4M=h4Mom*(tM/sum);
      }
      if(hM<sum || !G4QHadron(h4Mom).DecayIn3(nu4M,ba4M,pi4M))
      {
#ifdef debug
        G4int eA=theEnvironment.GetA();
        G4cout<<"***G4QEnv::FSI:T="<<hPDG<<"("<<hM<<")-> N="<<nuQPDG<<"(M="<<nucM<<") + B="
              <<barPDG<<"("<<barM<<")+N/pi="<<tPDG<<"("<<tM<<")="<<sum<<", A="<<eA<<G4endl;
#endif
        //if(!eA)
        //{
          G4QHadron* theLast = theCurr;        // Prototype of the pointer to theLastHadron
          G4QHadron* qH = new G4QHadron(theCurr); // Copy of the Current Hadron
#ifdef debug
          G4cout<<"***G4QE::FSI:#"<<ipo<<",4MQC="<<qH->Get4Momentum()<<qH->GetQC()<<G4endl;
#endif
          if(ipo+1<theQHadrons.size())         // If ipo<Last, swap CurHadr and theLastHadr
          {
            G4int nhd1=theQHadrons.size()-1;
            theLast = theQHadrons[nhd1];// Pointer to theLastHadron (ipo<L)
            G4LorentzVector l4M=theLast->Get4Momentum();
            G4QPDGCode lQP=theLast->GetQPDG();    // The QPDG of the last
            if(lQP.GetPDGCode()!=10) theCurr->SetQPDG(lQP); //CurHadr instead of LastHadr
            else theCurr->SetQC(theLast->GetQC());// CurHadrPDG instead of LastHadrPDG
#ifdef debug
            G4cout<<"---Warning---G4QE::FSI:l#"<<nhd1<<",4M="<<l4M<<",PDG="<<lQP<<G4endl;
#endif
            theCurr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
          }
          theQHadrons.pop_back();              // exclude theLastHadron pointer from OUTPUT
          delete theLast;  // *!! When killing, DON'T forget to delete the last QHadron !!*
          G4QContent cqQC=qH->GetQC();
          G4LorentzVector tqLV=qH->Get4Momentum();
          G4Quasmon* quasH = new G4Quasmon(cqQC,tqLV);//Create fakeQuasm
          if(!CheckGroundState(quasH,true))         // Try to correct by other hadrons
          {
#ifdef chdebug
            G4cout<<":W:G4QE::FSI:E="<<theEnvironment<<",Q="<<cqQC<<tqLV<<tqLV.m()<<G4endl;
#endif
            theQHadrons.push_back(qH);              // Fill as it is (delete equivalent)
          }
          else
          {
            delete qH;
            jpo--;
            nHadr=theQHadrons.size();
          }
          delete quasH;
          fOK=false;
        //}
        //else
        //{
        //  G4cerr<<"***G4QEnv::FSI:NoRec(3)E="<<theEnvironment<<eA<<",PDG="<<hPDG<<G4endl;
        //  throw G4QException("G4QEnv::FSI:ISO-dibaryon or 3n/3p DecayIn3 error");
      //}
      }
      if(fOK)
      {
        theQHadrons[ipo]->SetQPDG(nuQPDG);
        theQHadrons[ipo]->Set4Momentum(nu4M);
        G4QHadron* baH = new G4QHadron(barPDG,ba4M);
        theQHadrons.push_back(baH);               // (delete equivalent)
        G4QHadron* piH = new G4QHadron(tPDG,pi4M);
        theQHadrons.push_back(piH);               // (delete equivalent)
        nHadr=theQHadrons.size();
      }
    }
    else if(hBN>1 && !sBN && (cBN<0 || cBN>hBN)) // "nN+mD- or nP+mD++ decay"
    {
#ifdef debug
      G4cout<<"G4QE::FSI:nNmD-/nPmD++ #"<<ipo<<",P="<<hPDG<<",B="<<hBN<<",C="<<cBN<<G4endl;
#endif
      G4double hM=h4Mom.m();
      G4QPDGCode nuQPDG=nQPDG; // "(n+m)*N,m*pi-" case === Default
      G4double nucM = mNeut;
      G4int  barPDG = 2112;
      G4double barM = mNeut;
      G4int    nN   = hBN-1;   // a#of baryons - 1
      G4int   tPDG  = -211;
      G4double tM   = mPi;
      G4int    nPi  = -cBN;    // a#of Pi-'s
      if( cBN>hBN)             // reinitialization for the "(n+m)*P,m*pi+" case
      {
        nuQPDG = pQPDG;        // Substitute p for the first n
        nucM   = mProt;
        barPDG = 2212;         // Substitute p for the second n
        barM   = mProt;
        tPDG   = 211;          // Substitute pi+ for the first pi-
        nPi    = cBN-hBN;      // a#0f Pi+'s
      }
      if(nPi>1)   tM*=nPi;
      if(nN >1) barM*=nN;
      G4bool fOK=true;
      G4LorentzVector nu4M(0.,0.,0.,nucM);
      G4LorentzVector ba4M(0.,0.,0.,barM);
      G4LorentzVector pi4M(0.,0.,0.,tM);
      G4double sum=nucM+barM+tM;
      if(fabs(hM-sum)<=eps)
      {
        nu4M=h4Mom*(nucM/sum);
        ba4M=h4Mom*(barM/sum);
        pi4M=h4Mom*(tM/sum);
      }
      if(hM<sum || !G4QHadron(h4Mom).DecayIn3(nu4M,ba4M,pi4M))
      {
#ifdef debug
        G4cout<<"***G4QEnv::FSI:IsN M="<<hM<<","<<hPDG<<"->N="<<nuQPDG<<"(M="<<nucM<<")+"
              <<nN<<"*B="<<barPDG<<"(M="<<barM<<")+"<<nPi<<"*pi="<<tPDG<<"(M="<<tM<<")="
              <<nucM+barM+tM<<G4endl;
        G4cout<<"***G4QEnv::FSI:IsN BaryN="<<hBN<<",Charge="<<cBN<<",Stran="<<sBN<<G4endl;
#endif
        if(!theEnvironment.GetA())           // Emergency recovery
        {
          G4QHadron* theLast = theCurr;      // Prototype of the pointer to the Last Hadron
          G4QHadron* qH = new G4QHadron(theCurr); // Copy of the Current Hadron
          if(ipo+1<theQHadrons.size())       // If ipo<Last, swap CurHadr and theLastHadr
          {
            theLast = theQHadrons[theQHadrons.size()-1];// Ptr to theLastHadron (ipo<Last)
            G4QPDGCode lQP=theLast->GetQPDG();    // The QPDG of the last
            if(lQP.GetPDGCode()!=10) theCurr->SetQPDG(lQP); //CurHadr instead of LastHadr
            else theCurr->SetQC(theLast->GetQC());// CurHadrPDG instead of LastHadrPDG
            theCurr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
          }
          theQHadrons.pop_back();            // exclude theLastHadron pointer from OUTPUT
          delete theLast;  // *!! When killing, DON'T forget to delete the last QHadron !!*
          G4QContent cqQC=qH->GetQC();
          G4LorentzVector cq4M=qH->Get4Momentum();
          G4Quasmon* quasH = new G4Quasmon(cqQC,cq4M);// Create fakeQuasmon for the Last
          if(!CheckGroundState(quasH,true))         // Try to correct by other hadrons
          {
            G4cout<<"---Warning---G4QEnv::FSI:IN Failed, FillAsItIs: "<<cqQC<<cq4M<<G4endl;
            theQHadrons.push_back(qH);              // Fill as it is (delete equivalent)
          }
          else
          {
            delete qH;
            jpo--;
            nHadr=theQHadrons.size();
          }
          delete quasH;
          fOK=false;
        }
        else
        {
          G4ExceptionDescription ed;
          ed << "IN Multy ISO-dibaryon DecayIn3 did not succeed: IN,NoRec(4) Env="
             << theEnvironment << ",PDG=" << hPDG << G4endl;
          G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0002",
                      FatalException, ed);
        }
      }
      if(fOK)
      {
        theQHadrons[ipo]->SetQPDG(nuQPDG);
        theQHadrons[ipo]->Set4Momentum(nu4M);
        if(nN>1) ba4M=ba4M/nN;
        for(G4int ib=0; ib<nN; ib++)
        {
          G4QHadron* baH = new G4QHadron(barPDG,ba4M);
          theQHadrons.push_back(baH);               // (delete equivalent)
        }
        if(nPi>1) pi4M=pi4M/nPi;
        for(G4int im=0; im<nPi; im++)
        {
          G4QHadron* piH = new G4QHadron(tPDG,pi4M);
          theQHadrons.push_back(piH);               // (delete equivalent)
        }
        nHadr=theQHadrons.size();
      }
    }
#ifdef debug
    G4int           hNFrag= theQHadrons[ipo]->GetNFragments(); //Recover after swapping
    G4QContent      hQC   = theQHadrons[ipo]->GetQC();         // ...
    hPDG                  = theQHadrons[ipo]->GetPDGCode();    // ...
    h4Mom                 = theQHadrons[ipo]->Get4Momentum();  // ...
    ccs4M+=h4Mom;                                              // Calculate CurSum of Hadrs
    G4cout<<"G4QE::FSI:#"<<ipo<<": h="<<hPDG<<hQC<<",h4M="<<h4Mom<<h4Mom.m()<<",hNF="
          <<hNFrag<<G4endl;
#endif
    ipo=jpo;            // Take into account the roll back in case of the Last substitution
  }
#ifdef debug
  G4cout<<"G4QE::FSI: --->>CurrentControlSumOf4MomOfHadrons="<<ccs4M<<G4endl;
#endif
  nHadr=theQHadrons.size();
#ifdef chdebug
  // *** (1) Charge Control Sum Calculation for the Charge Conservation Check ***
  G4int ccContSum=0;                   // Intermediate ChargeControlSum 
  G4int cbContSum=0;                   // Intermediate BaryonNumberControlSum 
  if(nHadr)for(unsigned ic1=0; ic1<nHadr; ic1++) if(!(theQHadrons[ic1]->GetNFragments()))
  {
    ccContSum+=theQHadrons[ic1]->GetCharge();
    cbContSum+=theQHadrons[ic1]->GetBaryonNumber();
  }
  if(ccContSum-chContSum || cbContSum-bnContSum)
  {
    G4cout<<"*::*G4QE::FSI:(7)dC="<<ccContSum-chContSum<<",dB="<<cbContSum-bnContSum
          <<G4endl;
    //throw G4QException("G4QEnvironment::FSInteract: (1) Charge is not conserved");
  }
  // ***
#endif
  G4double p2cut=250000.;        // 250000=(2*p_Ferm)**2
  if(envA>0) p2cut/=envA*envA;
  //G4double p2cut2=0.;          //cut for the alpha creation
  //
  G4int bfCountM=3;
  if(envA>10) bfCountM*=(envA-1)/3;
  G4bool bfAct = true;
  G4int bfCount= 0;
  G4LorentzVector tmp4Mom=tot4Mom;
  G4LorentzVector postp4M(0.,0.,0.,0.);
  G4int nPost=intQHadrons.size();
  if(nPost) for(G4int psp=0; psp<nPost; psp++)
    if(!(intQHadrons[psp]->GetNFragments())) postp4M+=intQHadrons[psp]->Get4Momentum();
  while(bfAct&&bfCount<bfCountM) // "Infinite" LOOP of the ThermoNuclearBackFusion
  {
    tot4Mom=tmp4Mom-postp4M;     // Prepare tot4Mom for the "En/Mom conservation" reduction
    bfAct=false;
    bfCount++;
    nHadr=theQHadrons.size();
    if(nHadr) for(unsigned hadron=0; hadron<theQHadrons.size(); hadron++)// BackFusion LOOP
    {
      G4QHadron* curHadr = theQHadrons[hadron]; // Get a pointer to the current Hadron
      G4int         hPDG = curHadr->GetPDGCode();
      G4QPDGCode    hQPDG(hPDG);
      G4double      hGSM = hQPDG.GetMass();     // Ground State Mass of the first fragment
#ifdef debug
      G4cout<<"G4QE::FSI:LOOP START,h#"<<hadron<<curHadr->Get4Momentum()<<hPDG<<G4endl;
#endif
      if(hPDG==89999003||hPDG==90002999)
      {
        G4cout<<"---WARNING---G4QE::FSI:**nD-/pD++**(3),PDG="<<hPDG<<" CORRECTION"<<G4endl;
        G4LorentzVector h4Mom=curHadr->Get4Momentum();
        G4double      hM=h4Mom.m();
        G4QPDGCode fQPDG=nQPDG;
        G4double     fM =mNeut;
        G4int      sPDG =2112;
        G4double     sM =mNeut;
        G4int      tPDG =-211;
        G4double     tM =mPi;
        if(hPDG==90002999)
        {
          fQPDG=pQPDG;
          fM   =mProt;
          sPDG =2212;
          sM   =mProt;
          tPDG =211;
        }
        G4bool fOK=true;
        G4LorentzVector f4M(0.,0.,0.,fM);
        G4LorentzVector s4M(0.,0.,0.,sM);
        G4LorentzVector t4M(0.,0.,0.,tM);
        G4double sum=fM+sM+tM;
        if(fabs(hM-sum)<=eps)
        {
          f4M=h4Mom*(fM/sum);
          s4M=h4Mom*(sM/sum);
          t4M=h4Mom*(tM/sum);
        }
        else if(hM<sum || !G4QHadron(h4Mom).DecayIn3(f4M,s4M,t4M))
        {
          G4cout<<"---WARNING---G4QE::FSI: Still trying, NDM="<<hM<<"->"<<fM<<"("<<fQPDG
                <<")+"<<sM<<"("<<sPDG<<")+"<<tM<<"("<<tPDG<<")="<<sum<<G4endl;
          if(!theEnvironment.GetA())
          {
            G4QHadron* theLast = curHadr;          // Prototype of Pointer to theLastHadron
            G4QHadron* qH = new G4QHadron(curHadr);// Copy of the Current Hadron
            if(hadron+1<theQHadrons.size())        // If hadr<Last,swap CurHadr & LastHadr
            {
              theLast = theQHadrons[theQHadrons.size()-1]; // Pointer to LastHadr (nh<Last)
              G4QPDGCode lQP=theLast->GetQPDG();    // The QPDG of the last
              if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP); //CurHadr instead of LastHadr
              else curHadr->SetQC(theLast->GetQC());// CurHadrPDG instead of LastHadrPDG
              curHadr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
            }
            theQHadrons.pop_back();        // exclude theLastHadron pointer from the OUTPUT
            delete theLast;// *!! When killing, DON'T forget to delete the last QHadron !!*
            G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());//Fake Quasmon
            if(!CheckGroundState(quasH,true))      // Try to correct by other hadrons
            {
              G4cout<<"---Warning---G4QE::FSI:NDel Failed LeaveAsItIs, 4m="<<h4Mom<<G4endl;
              theQHadrons.push_back(qH);           // Leave as it is (delete equivalent)
            }
            else
            {
              delete qH;
              nHadr=theQHadrons.size();
            }
            delete quasH;
            fOK=false;
          }
          else
          {
            // G4cout<<"***G4QEnv::FSI: No ND recovery Env="<<theEnvironment<<G4endl;
            // throw G4QException("G4QEnv::FSI:ND DecayIn3 did not succeed");
            G4ExceptionDescription ed;
            ed << "ND DecayIn3 did not succeed: No ND recovery Env="
               << theEnvironment << G4endl;
            G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0003",
                        FatalException, ed);
          }
        }
        if(fOK)
        {
          curHadr->SetQPDG(fQPDG);
          curHadr->Set4Momentum(f4M);
          G4QHadron* sH = new G4QHadron(sPDG,s4M);
          theQHadrons.push_back(sH);               // (delete equivalent)
          G4QHadron* tH = new G4QHadron(tPDG,t4M);
          theQHadrons.push_back(tH);               // (delete equivalent)
        }
        hPDG = curHadr->GetPDGCode();            // Change PDG Code of theFirstFragment
        hQPDG= G4QPDGCode(hPDG);
        hGSM = hQPDG.GetMass();                  // Change GroundStateMass of theFirstFragm
      }
      nHadr=theQHadrons.size();
      if(hPDG==89001001||hPDG==89002000||hPDG==89000002)
      {
        G4cout<<"---WARNING---G4QE::FSI:***(K+N)*** (2),PDG="<<hPDG<<" CORRECTION"<<G4endl;
        G4LorentzVector h4Mom=curHadr->Get4Momentum();
        G4double      hM=h4Mom.m();
        G4QPDGCode fQPDG=nQPDG;
        G4double     fM =mNeut;
        G4int      sPDG =311;
        G4double     sM =mK0;
        if(hPDG==89000002)
        {
          fQPDG=pQPDG;
          fM   =mProt;
          sPDG =321;
          sM   =mK;
        }
        if(hPDG==89001001)
        {
          if(hM<mK0+mProt || G4UniformRand()>.5)
          {
            sPDG =321;
            sM   =mK;
          }
          else
          {
            fQPDG=pQPDG;
            fM   =mProt;
          }
        }
        G4bool fOK=true;
        G4LorentzVector f4M(0.,0.,0.,fM);
        G4LorentzVector s4M(0.,0.,0.,sM);
        G4double sum=fM+sM;
        if(fabs(hM-sum)<=eps)
        {
          f4M=h4Mom*(fM/sum);
          s4M=h4Mom*(sM/sum);
        }
        else if(hM<sum || !G4QHadron(h4Mom).DecayIn2(f4M,s4M))
        {
          G4cout<<"---WARNING---G4QE::FSI: Still trying (2),NDM="<<hM<<"->"<<fM<<"("<<fQPDG
                <<")+"<<sM<<"("<<sPDG<<")="<<sum<<G4endl;
          if(!theEnvironment.GetA())
          {
            G4QHadron* theLast = curHadr;          // Prototype of Pointer to theLastHadron
            G4QHadron* qH = new G4QHadron(curHadr);// Copy of the Current Hadron
            if(hadron+1<theQHadrons.size())        // If hadr<Last, swap CurHadr & LastHadr
            {
              theLast = theQHadrons[theQHadrons.size()-1]; // Pointer to LastHadr (nh<Last)
              G4QPDGCode lQP=theLast->GetQPDG();    // The QPDG of the last
              if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP); //CurHadr instead of LastHadr
              else curHadr->SetQC(theLast->GetQC());// CurHadrPDG instead of LastHadrPDG
              curHadr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
            }
            theQHadrons.pop_back();        // exclude theLastHadron pointer from the OUTPUT
            delete theLast;// *!! When killing, DON'T forget to delete the last QHadron !!*
            G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());//Fake Quasmon
            if(!CheckGroundState(quasH,true))      // Try to correct by other hadrons
            {
              G4cout<<"---Warning---G4QE::FSI:KN Failed LeaveAsItIs 4m="<<h4Mom<<G4endl;
              theQHadrons.push_back(qH);           // Leave as it is (delete equivalent)
            }
            else
            {
              delete qH;
              nHadr=theQHadrons.size();
            }
            delete quasH;
            fOK=false;
          }
          else
          {
            // G4cerr<<"***G4QEnv::FSI: No KN recovery Env="<<theEnvironment<<G4endl;
            // throw G4QException("G4QEnv::FSI:KN DecayIn2 did not succeed");
            G4ExceptionDescription ed;
            ed << "KN DecayIn2 did not succeed: No KN recovery Env="
               << theEnvironment << G4endl;
            G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0004",
                        FatalException, ed);
          }
        }
        if(fOK)
        {
          curHadr->SetQPDG(fQPDG);
          curHadr->Set4Momentum(f4M);
          G4QHadron* sH = new G4QHadron(sPDG,s4M);
          theQHadrons.push_back(sH);               // (delete equivalent)
        }
        hPDG = curHadr->GetPDGCode();            // Change PDG Code of theFirstFragment
        hQPDG= G4QPDGCode(hPDG);
        hGSM = hQPDG.GetMass();                  // Change GroundStateMass of theFirstFragm
      }
      nHadr=theQHadrons.size();
      G4int           hS = curHadr->GetStrangeness();
      G4int           hF = curHadr->GetNFragments();
      G4LorentzVector h4m= curHadr->Get4Momentum();
      G4double hM        = h4m.m();             // Real Mass of the first fragment
      G4int hB           = curHadr->GetBaryonNumber();
      //////////////////////G4int hC           = curHadr->GetCharge();
#ifdef debug
      if(!hF && ( (hPDG>80000000 && hPDG<90000000) || hPDG==90000000 ||
                  (hPDG>90000000 && (hPDG%1000000>200000 || hPDG%1000>300) ) ) )
        G4cout<<"**G4QEnv::FSInteraction: PDG("<<hadron<<")="<<hPDG<<", M="<<hM<<G4endl;
#endif
#ifdef debug
      G4cout<<"G4QE::FSI:h="<<hPDG<<",S="<<hS<<",B="<<hB<<",#"<<hadron<<"<"<<nHadr<<G4endl;
#endif
      //if(hadron&&!hF&&hB>0&&!hS&&(nHadr>3||hB<2)) // ThermoBackFus (VIMP for gamA TotCS)
      //if(hadron&&!hF&&hB>0&&!hS&&(nHadr>2||hB<4)) // ThermoBackFus (VIMP for gamA TotCS)
      if(hadron&&!hF&&hB>0&&!hS) // ThermoBackFusion cond. (VIMP for gamA TotCS)
      //if(hadron&&!hF&&hB>0&&hB<4&&!hS) // ThermoBackFusion cond. (VIMP for gamA TotCS)
      //if(hadron&&!hF&&hB>0&&!hS&&nHadr>2)//ThermoBackFusion MAX condition (VIMP for gamA)
      //if(2>3)                         // Close the ThermoBackFusion (VIMP for gamA TotCS)
      {
#ifdef debug
        //if(nHadr==3)
          G4cout<<"G4QE::FSI: h="<<hPDG<<",B="<<hB<<",h#"<<hadron<<" < nH="<<nHadr<<G4endl;
#endif
        G4QContent hQC = curHadr->GetQC();
        if(hadron&&!hF&&hB>0) for(unsigned pt=0; pt<hadron; pt++)
        {
          G4QHadron* backH = theQHadrons[pt];   // Get pointer to one of thePreviousHadrons
          G4int   bF = backH->GetNFragments();
          G4LorentzVector b4m= backH->Get4Momentum();
          G4double bM= b4m.m();                 // Real Mass of the second fragment
          G4QContent bQC = backH->GetQC();
          G4int bPDG=bQC.GetZNSPDGCode();
          G4QPDGCode bQPDG(bPDG);
          G4double bGSM=bQPDG.GetMass();        // Ground State Mass of the second fragment
          G4int   bB = backH->GetBaryonNumber();

          //////////////////G4int   bC = backH->GetCharge();
          G4QContent sQC=bQC+hQC;
          G4int sPDG=sQC.GetZNSPDGCode();
          G4QPDGCode sQPDG(sPDG);
          G4double tM=sQPDG.GetMass();
          G4LorentzVector s4M=h4m+b4m;
          G4double sM2=s4M.m2();
          G4double sM=sqrt(sM2);
          G4double dsM2=sM2+sM2;
          G4double rm=bM-hM;
          G4double sm=bM+hM;
          G4double pCM2=(sM2-rm*rm)*(sM2-sm*sm)/(dsM2+dsM2);
          G4int   bS = backH->GetStrangeness();
#ifdef debug
          //if(nHadr==3)
          G4cout<<"G4QE::FSI:"<<pt<<",B="<<bB<<",S="<<bS<<",p="<<pCM2<<"<"<<p2cut<<",hB="
                <<hB<<",bM+hM="<<bM+hM<<">tM="<<tM<<",tQC="<<sQC<<G4endl;
#endif
          //if(!bF&&(bB==1||hB==1)&&bM+hM>tM+.001&&pCM2<p2cut)      // Only baryons == pcut
          //if(!bF&&!bS&&(bB==1&&hB>0||hB==1&&bB>0)&&bM+hM>tM+.00001
          //   && (pCM2<p2cut&&nHadr>3||pCM2<p2cut2&&nHadr==3))
          //if(!bF&&(bB==1||hB==1)&&bM+hM>tM+.001&&(pCM2<p2cut&&nHadr>3 ||
          //   pCM2<p2cut2&&nHadr==3&&bPDG>90000000))
          //if(!bF&&bB<4&&bM+hM>tM+.001&&pCM2<p2cut)
          if(!bF&&!bS&&bB>0&&bM+hM>tM+.001&&pCM2<p2cut)
          //if(!bF&&bB<4&&bM+hM>tM+.001&&(pCM2<p2cut || bB+hB==4&&pCM2<p2cut2))
          //if(!bF&&(bB==1||hB==1)&&(nHadr>3||bPDG>90000000)&&bM+hM>tM+.001&&pCM2<p2cut)
          //if(!bF&&(bB==1&&!bC||hB==1&&!hC)&&bM+hM>tM+.001&&pCM2<p2cut)// Only n == pcut
          //if(!bF&&(bB==1||hB==1)&&bM+hM>tM+.001&&sM-bM-hM<cut)  // Only baryons == ecut
          //if(!bF&&bB&&bB<fL&&bM+hM>tM+.001&&sM-bM-hM<cut)    // Light fragments == ecut
          {
#ifdef fdebug
            G4int bPDG = backH->GetPDGCode();
            if(sPDG==89999003||sPDG==90002999||sPDG==89999002||sPDG==90001999)
              G4cout<<"G4QE::FSI:**nD-/pD++**,h="<<hPDG<<",hB="<<hB<<",b="<<bPDG<<",bB="
                    <<bB<<G4endl;
            //if(nHadr==3)
            G4cout<<"G4QE::FSI:*FUSION*#"<<hadron<<"["<<hPDG<<"]"<<hM<<"+#"<<pt<<"["<<bPDG
                  <<"]"<<bM<<"="<<bM+hM<<", sM="<<sM<<">["<<sPDG<<"]"<<tM<<",p2="<<pCM2
                  <<"<"<<p2cut<<G4endl;
#endif
            bfAct=true;
            //@@Temporary decay in gamma
            G4bool three=false;
            G4QPDGCode fQPDG=sQPDG;
            G4QPDGCode rQPDG=gQPDG;
            hQPDG=gQPDG;
            G4LorentzVector f4Mom(0.,0.,0.,tM);
            G4LorentzVector g4Mom(0.,0.,0.,0.);
            G4LorentzVector t4Mom(0.,0.,0.,0.);
            if(sPDG==89999002)                               // A=1
            {
              fQPDG=nQPDG;
              rQPDG=pimQPDG;
              f4Mom=G4LorentzVector(0.,0.,0.,mNeut);
              g4Mom=G4LorentzVector(0.,0.,0.,mPi);
            }
            else if(sPDG==90001999)
            {
              fQPDG=pQPDG;
              rQPDG=pipQPDG;
              f4Mom=G4LorentzVector(0.,0.,0.,mProt);
              g4Mom=G4LorentzVector(0.,0.,0.,mPi);
            }
            else if(sPDG==90000002)                        // A=2
            {
              fQPDG=nQPDG;
              rQPDG=nQPDG;
              f4Mom=G4LorentzVector(0.,0.,0.,mNeut);
              g4Mom=f4Mom;
            }
            else if(sPDG==90002000)
            {
              fQPDG=pQPDG;
              rQPDG=pQPDG;
              f4Mom=G4LorentzVector(0.,0.,0.,mProt);
              g4Mom=f4Mom;
            }
            else if(sPDG==92000000)
            {
              fQPDG=lQPDG;
              rQPDG=lQPDG;
              f4Mom=G4LorentzVector(0.,0.,0.,mLamb);
              g4Mom=f4Mom;
              if(sM>mSigZ+mSigZ)             // Sigma0+Sigma0 is possible
              {                                  // @@ Only two particles PS is used
                G4double sma=mLamb+mLamb;        // Lambda+Lambda sum
                G4double pr1=0.;                 // Prototype to avoid sqrt(-)
                if(sM>sma) pr1=sqrt((sM2-sma*sma)*sM2); // Lamb+Lamb PS
                sma=mLamb+mSigZ;                 // Lambda+Sigma0 sum
                G4double smi=mSigZ-mLamb;        // Sigma0-Lambda difference
                G4double pr2=pr1;                // Prototype of +L +S0 PS
                if(sM>sma && sM>smi) pr2+=sqrt((sM2-sma*sma)*(sM2-smi*smi));
                sma=mSigZ+mSigZ;                 // Sigma0+Sigma0 sum
                G4double pr3=pr2;                // Prototype of +Sigma0+Sigma0 PS
                if(sM>sma) pr3+=sqrt((sM2-sma*sma)*sM2);
                G4double hhRND=pr3*G4UniformRand(); // Randomize PS
                if(hhRND>pr2)                    // --> "ENnv+Sigma0+Sigma0" case
                {                                //
                  fQPDG=s0QPDG;
                  f4Mom=G4LorentzVector(0.,0.,0.,mSigZ);
                  rQPDG=s0QPDG;
                  g4Mom=f4Mom;
                }                                //
                else if(hhRND>pr1)               // --> "ENnv+Sigma0+Lambda" case
                {                                //
                  fQPDG=s0QPDG;
                  f4Mom=G4LorentzVector(0.,0.,0.,mSigZ);
                }                                //
              }
              else if(sM>mSigZ+mLamb)            // Lambda+Sigma0 is possible
              {                                  // @@ Only two particles PS is used
                G4double sma=mLamb+mLamb;        // Lambda+Lambda sum
                G4double pr1=0.;                 // Prototype to avoid sqrt(-)
                if(sM>sma) pr1=sqrt((sM2-sma*sma)*sM2); // Lamb+Lamb PS
                sma=mLamb+mSigZ;                 // Lambda+Sigma0 sum
                G4double smi=mSigZ-mLamb;        // Sigma0-Lambda difference
                G4double pr2=pr1;                // Prototype of +L +S0 PS
                if(sM>sma && sM>smi) pr2+=sqrt((sM2-sma*sma)*(sM2-smi*smi));
                if(pr2*G4UniformRand()>pr1)      // --> "ENnv+Sigma0+Lambda" case
                {                                //
                  fQPDG=s0QPDG;
                  f4Mom=G4LorentzVector(0.,0.,0.,mSigZ);
                }                                //
              }                                  //
            }
            else if(sPDG==89999003)                       // A=2
            {
              hQPDG=nQPDG;
              rQPDG=nQPDG;
              fQPDG=pimQPDG;
              t4Mom=G4LorentzVector(0.,0.,0.,mNeut);
              g4Mom=G4LorentzVector(0.,0.,0.,mNeut);
              f4Mom=G4LorentzVector(0.,0.,0.,mPi);
              three=true;
            }
            else if(sPDG==90002999)
            {
              hQPDG=pQPDG;
              rQPDG=pQPDG;
              fQPDG=pipQPDG;
              t4Mom=G4LorentzVector(0.,0.,0.,mProt);
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mPi);
              three=true;
            }
            else if(sPDG==90000003)                        // A=3
            {
              hQPDG=nQPDG;
              rQPDG=nQPDG;
              fQPDG=nQPDG;
              t4Mom=G4LorentzVector(0.,0.,0.,mNeut);
              g4Mom=G4LorentzVector(0.,0.,0.,mNeut);
              f4Mom=G4LorentzVector(0.,0.,0.,mNeut);
              three=true;
            }
            else if(sPDG==90003000)
            {
              hQPDG=pQPDG;
              rQPDG=pQPDG;
              fQPDG=pQPDG;
              t4Mom=G4LorentzVector(0.,0.,0.,mProt);
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mProt);
              three=true;
            }
            else if(sPDG==90001003)                     // A=4
            {
              rQPDG=nQPDG;
              fQPDG=tQPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mNeut);
              f4Mom=G4LorentzVector(0.,0.,0.,mTrit);
            }
            else if(sPDG==90003001)
            {
              rQPDG=pQPDG;
              fQPDG=he3QPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mHe3);
            }
            else if(sPDG==90002003)                         // A=5
            {
              rQPDG=nQPDG;
              fQPDG=aQPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mNeut);
              f4Mom=G4LorentzVector(0.,0.,0.,mAlph);
            }
            else if(sPDG==90003002)
            {
              rQPDG=pQPDG;
              fQPDG=aQPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mAlph);
            }
            else if(sPDG==90004002)                          // A=6
            {
              hQPDG=pQPDG;
              rQPDG=pQPDG;
              fQPDG=aQPDG;
              t4Mom=G4LorentzVector(0.,0.,0.,mProt);
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mAlph);
              three=true;
            }
            else if(sPDG==90002005)                        // A=7
            {
              rQPDG=nQPDG;
              fQPDG=a6QPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mNeut);
              f4Mom=G4LorentzVector(0.,0.,0.,mHe6);
            }
            else if(sPDG==90005002)
            {
              rQPDG=pQPDG;
              fQPDG=be6QPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mBe6);
            }
            else if(sPDG==90004004)                        // A=8
            {
              fQPDG=aQPDG;
              rQPDG=aQPDG;
              f4Mom=G4LorentzVector(0.,0.,0.,mAlph);
              g4Mom=f4Mom;
            }
            //else if(sPDG==90006002)
            //{
            //  hQPDG=pQPDG;
            //  rQPDG=pQPDG;
            //  fQPDG=be6QPDG;
            //  t4Mom=G4LorentzVector(0.,0.,0.,mProt);
            //  g4Mom=G4LorentzVector(0.,0.,0.,mProt);
            //  f4Mom=G4LorentzVector(0.,0.,0.,mBe6);
            //  three=true;
            //}
            //else if(sPDG==90002006)
            //{
            //  hQPDG=nQPDG;
            //  rQPDG=nQPDG;
            //  fQPDG=a6QPDG;
            //  t4Mom=G4LorentzVector(0.,0.,0.,mNeut);
            //  g4Mom=G4LorentzVector(0.,0.,0.,mNeut);
            //  f4Mom=G4LorentzVector(0.,0.,0.,mHe6);
            //  three=true;
            //}
            else if(sPDG==90002007)                      // A=9
            {
              rQPDG=nQPDG;
              fQPDG=a8QPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mNeut);
              f4Mom=G4LorentzVector(0.,0.,0.,mHe8);
            }
            else if(sPDG==90005004)                      // A=9
            {
              rQPDG=pQPDG;
              fQPDG=aQPDG;
              hQPDG=aQPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mAlph);
              t4Mom=G4LorentzVector(0.,0.,0.,mAlph);
              three=true;
            }
            else if(sPDG==90007002)                      // A=9
            {
              rQPDG=pQPDG;
              fQPDG=c8QPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mC8);
            }
            else if(sPDG==90008004)                      // A=12
            {
              hQPDG=pQPDG;
              rQPDG=pQPDG;
              fQPDG=c10QPDG;
              t4Mom=G4LorentzVector(0.,0.,0.,mProt);
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mC10);
              three=true;
            }
            else if(sPDG==90009006)                     // A=15
            {
              rQPDG=pQPDG;
              fQPDG=o14QPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mO14);
            }
            else if(sPDG==90009007)                     // A=16
            {
              rQPDG=pQPDG;
              fQPDG=o15QPDG;
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mO15);
            }
            else if(sPDG==90010006)                     // A=16
            {
              hQPDG=pQPDG;
              rQPDG=pQPDG;
              fQPDG=o14QPDG;
              t4Mom=G4LorentzVector(0.,0.,0.,mProt);
              g4Mom=G4LorentzVector(0.,0.,0.,mProt);
              f4Mom=G4LorentzVector(0.,0.,0.,mO14);
              three=true;
            }
#ifdef debug
            G4cout<<"G4QE::FSI: "<<three<<",r="<<rQPDG<<",f="<<fQPDG<<",t="<<hQPDG<<G4endl;
#endif
            if(!three)
            {
              if(!G4QHadron(s4M).DecayIn2(f4Mom,g4Mom))
              {
                G4ExceptionDescription ed;
                ed << "Fusion (1) DecIn2 error: (2)*FUSION*,tM[" << sPDG << "]="
                   << tM << ">sM=" << sM << " of " << h4m << hM << hQC << hGSM
                   << " & " << b4m << bM << bQC << bGSM << G4endl;
                G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0005",
                            FatalException, ed);
              }
              else
              {
#ifdef debug
                G4cout<<"G4QE::FSI:*FUSION IS DONE*,fPDG="<<sPDG<<",PDG1="<<hPDG<<",PDG2="
                      <<bPDG<<G4endl;
#endif
                curHadr->SetQPDG(fQPDG);
                curHadr->Set4Momentum(f4Mom);
                backH->SetQPDG(rQPDG);
                backH->Set4Momentum(g4Mom);
#ifdef debug
                G4cout<<"G4QE::FSI:h="<<h4m<<",b="<<b4m<<",s="<<s4M<<G4endl;
                G4cout<<"G4QE::FSI:f="<<f4Mom<<",g="<<g4Mom<<",s="<<f4Mom+g4Mom<<G4endl;
#endif
              }
            }
            else
            {
              if(!G4QHadron(s4M).DecayIn3(f4Mom,g4Mom,t4Mom))
              {
                G4ExceptionDescription ed;
                ed << "Fusion(2) DecayIn3 error: (3)*FUSION*,tM[" << sPDG
                   << "]=" << tM << ">sM=" << sM << " of " << h4m << hM << hQC
                   << hGSM << " & " << b4m << bM << bQC << bGSM << G4endl;
                G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0006",
                            FatalException, ed);
              }
              else
              {
#ifdef debug
                G4cout<<"G4QE::FSI:DONE,n="<<nHadr<<",PDG="<<sPDG<<",1="<<hPDG<<",2="<<bPDG
                      <<G4endl;
#endif
                curHadr->SetQPDG(fQPDG);
                curHadr->Set4Momentum(f4Mom);
                backH->SetQPDG(rQPDG);
                backH->Set4Momentum(g4Mom);
                G4QHadron* newH = new G4QHadron(hQPDG.GetPDGCode(),t4Mom);
                theQHadrons.push_back(newH);      // (delete equivalent for newH)
                nHadr=theQHadrons.size();
#ifdef debug
                G4cout<<"G4QE::FSI:h="<<h4m<<",b="<<b4m<<G4endl;
                G4cout<<"G4QE::FSI:s="<<s4M<<" = Sum"<<f4Mom+g4Mom+t4Mom<<G4endl;
                G4cout<<"G4QE::FSI:*Products*,nH="<<nHadr<<",f="<<fQPDG<<f4Mom<<",b="
                      <<rQPDG<<g4Mom<<",new="<<hQPDG<<t4Mom<<",nH="<<nHadr<<",nD="
                      <<theQHadrons.size()<<G4endl;
#endif
              }
            }
            tot4Mom+=b4m;                       // Instead of the fused hadron
            tot4Mom-=g4Mom;                     // subtract from the total the new hadron
            /////////////////////////////break; // Make fusion only for one (?)
            // Instead the curHadr parameters should be updated ______
            hQPDG=fQPDG;
            hPDG=hQPDG.GetPDGCode();
            hQC=fQPDG.GetQuarkContent();
            hS=hQC.GetStrangeness();
            hB=hQC.GetBaryonNumber();
            hGSM = hQPDG.GetMass();
            h4m=f4Mom;
            hM=h4m.m();
            // End of Instead ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#ifdef debug
            G4cout<<"G4QE::FSI:cH4M="<<curHadr->Get4Momentum()<<G4endl;
#endif
          } // End of the fusion check
        } // End of the LOOP over previous hadrons
      } // End of the FUSION check
#ifdef chdebug
      // *** (2) Charge Control Sum Calculation for the Charge Conservation Check ***
      ccContSum=0;                   // Intermediate ChargeControlSum 
      cbContSum=0;                   // Intermediate BaryonNumberControlSum 
      if(nHadr)for(unsigned ic2=0;ic2<nHadr;ic2++) if(!(theQHadrons[ic2]->GetNFragments()))
      {
        ccContSum+=theQHadrons[ic2]->GetCharge();
        cbContSum+=theQHadrons[ic2]->GetBaryonNumber();
      }
      unsigned pHadr=intQHadrons.size();
      if(pHadr)for(unsigned ic3=0;ic3<pHadr;ic3++) if(!(intQHadrons[ic3]->GetNFragments()))
      {
        ccContSum+=intQHadrons[ic3]->GetCharge();
        cbContSum+=intQHadrons[ic3]->GetBaryonNumber();
      }
      if(ccContSum-chContSum || cbContSum-bnContSum)
      {
        G4cout<<"*::*G4QE::FSI:(8) dC="<<ccContSum-chContSum<<",dB="<<cbContSum-bnContSum
              <<G4endl;
        //throw G4QException("G4QEnvironment::FSInteract: (2) Charge is not conserved");
      }
      // ***
#endif
      G4LorentzVector cH4Mom = curHadr->Get4Momentum(); // 4-mom of the current hadron
      tot4Mom-=cH4Mom;                          // Reduce theTotal 4mom by theCurrent 4mom
      totCharge-=curHadr->GetCharge();          // @!@
      totBaryoN-=curHadr->GetBaryonNumber();    // @!@
#ifdef pdebug
        G4cout<<"G4QE::FSI:Cur4M="<<tot4Mom<<",tC="<<totCharge<<",tB="<<totBaryoN<<G4endl;
#endif
      if(hadron+1==nHadr)                       // The Last Hadron in the set
      {
#ifdef pdebug
        G4cout<<"G4QE::FSI:Last4M="<<tot4Mom<<",tC="<<totCharge<<",tB="<<totBaryoN<<G4endl;
#endif
        G4double misM=0.;                       // @!@
        G4int mPDG=0;                           // @!@
        if(totCharge>=0 && totBaryoN > 0)       // @!@
        {                                       // @!@
          mPDG=90000000+999*totCharge+totBaryoN;// @!@
          misM=G4QPDGCode(mPDG).GetMass();      // @!@
        }                                       // @!@
        G4double re =tot4Mom.e();
        G4double rpx=tot4Mom.px();
        G4double rpy=tot4Mom.py();
        G4double rpz=tot4Mom.pz();
        G4double re2=re*re;                     // @!@
        G4double dmo=rpx*rpx+rpy*rpy+rpz*rpz;
        G4double dem=re2+dmo;
        G4double dm2=re2-dmo;                   // @!@
        G4double sdm=0.;                        // @!@
        if(dm2>0.) sdm=std::sqrt(dm2);          // @!@
#ifdef debug
        G4cout<<"G4QE::FSI: Is En&Mom conserved? t4M="<<tot4Mom<<",dM="<<sdm<<", mM="<<misM
              <<",mPDG="<<mPDG<<",dCH="<<totCharge<<",dBN="<<totBaryoN<<G4endl;
#endif
        G4LorentzVector cor4M(0.,0.,0.,0.);     // Prototype for the missing particle
        if(dem>0.1)                             // Energy or momentum is not conserved
        {
          G4bool corf=false;
#ifdef pdebug
          G4cout<<"--Warning--G4QE::FSI:dE/Mc4M="<<tot4Mom<<sdm<<". Correct it!"<<G4endl;
#endif
	    if(sdm < .01 || (re2 > 0. && !totCharge && !totBaryoN && sdm/re2 < .0001)) // @!@
          {
#ifdef pdebug
            G4cout<<"...G4QE::FSI:E/M conservation is corrected by a photon"<<G4endl;
#endif
            cor4M=tot4Mom;                                 // Complete correction
            G4QHadron* theH = new G4QHadron(22,tot4Mom);
            theQHadrons.push_back(theH);    // (delete equivalent for the proton)
            corf=true;
          }
          else                                             // @!@
          {
            if(dmo<0.0001 && re>900.)               // MomentumIsConserved - recoverMissing
            {
              if(fabs(re-mNeut)<.01)
              {
#ifdef pdebug
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by neutron"<<G4endl;
#endif
                cor4M=G4LorentzVector(0.,0.,0.,mNeut);
                G4QHadron* theH = new G4QHadron(90000001,G4LorentzVector(0.,0.,0.,mNeut));
                theQHadrons.push_back(theH);  // (delete equivalent for the proton)
                corf=true;
              }
              else if(fabs(re-mProt)<.01)
              {
#ifdef pdebug
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by proton"<<G4endl;
#endif
                cor4M=G4LorentzVector(0.,0.,0.,mProt);
                G4QHadron* theH = new G4QHadron(90001000,G4LorentzVector(0.,0.,0.,mProt));
                theQHadrons.push_back(theH);  // (delete equivalent for the proton)
                corf=true;
              }
              else if(fabs(re-mDeut)<.01)
              {
#ifdef pdebug
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by deuteron"<<G4endl;
#endif
                cor4M=G4LorentzVector(0.,0.,0.,mDeut);
                G4QHadron* theH = new G4QHadron(90001001,G4LorentzVector(0.,0.,0.,mDeut));
                theQHadrons.push_back(theH);  // (delete equivalent for the proton)
                corf=true;
              }
              else if(fabs(re-mTrit)<.01)
              {
#ifdef pdebug
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by tritium"<<G4endl;
#endif
                cor4M=G4LorentzVector(0.,0.,0.,mTrit);
                G4QHadron* theH = new G4QHadron(90001002,G4LorentzVector(0.,0.,0.,mTrit));
                theQHadrons.push_back(theH);  // (delete equivalent for the proton)
                corf=true;
              }
              else if(fabs(re-mHe3)<.01)
              {
#ifdef pdebug
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by He3"<<G4endl;
#endif
                cor4M=G4LorentzVector(0.,0.,0.,mHe3);
                G4QHadron* theH = new G4QHadron(90002001,G4LorentzVector(0.,0.,0.,mHe3));
                theQHadrons.push_back(theH);  // (delete equivalent for the proton)
                corf=true;
              }
              else if(fabs(re-mAlph)<.01)
              {
#ifdef pdebug
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by alpha"<<G4endl;
#endif
                cor4M=G4LorentzVector(0.,0.,0.,mAlph);
                G4QHadron* theH = new G4QHadron(90002002,G4LorentzVector(0.,0.,0.,mAlph));
                theQHadrons.push_back(theH);  // (delete equivalent for the proton)
                corf=true;
              }
              else if(fabs(re-mNeut-mNeut)<.01)
              {
                cor4M=G4LorentzVector(0.,0.,0.,mNeut+mNeut);
#ifdef pdebug
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by 2 neutrons"<<G4endl;
#endif
                G4QHadron* theH1 = new G4QHadron(90000001,G4LorentzVector(0.,0.,0.,mNeut));
                theQHadrons.push_back(theH1); // (delete equivalent for the proton)
                G4QHadron* theH2 = new G4QHadron(90000001,G4LorentzVector(0.,0.,0.,mNeut));
                theQHadrons.push_back(theH2); // (delete equivalent for the proton)
                corf=true;
              }
              else if(fabs(re-mProt-mProt)<.01)
              {
#ifdef pdebug
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by 2 protons"<<G4endl;
#endif
                cor4M=G4LorentzVector(0.,0.,0.,mProt+mProt);
                G4QHadron* theH1 = new G4QHadron(90001000,G4LorentzVector(0.,0.,0.,mProt));
                theQHadrons.push_back(theH1); // (delete equivalent for the proton)
                G4QHadron* theH2 = new G4QHadron(90001000,G4LorentzVector(0.,0.,0.,mProt));
                theQHadrons.push_back(theH2); // (delete equivalent for the proton)
                corf=true;
              }
              else G4Exception("G4QEnvironment::FSInteract()", "HAD_CHPS_0007",
                               JustWarning, "Try heavier nuclei at rest");
            }
            else if(std::abs(sdm-misM) < 0.01)        // on flight correction @!@
            {
#ifdef pdebug
              G4cout<<"...G4QE::FSI:E/M conservation is corrected by ResidualNucl"<<G4endl;
#endif
              if(!misM) mPDG=22;
              G4QHadron* theH = new G4QHadron(mPDG,tot4Mom); // Create Residual Nucleus
              cor4M=tot4Mom;                  // Complete correction
              if(std::fabs(sdm-misM) <= 0.01) theQHadrons.push_back(theH); // As is
              else EvaporateResidual(theH);   // Evaporate Residual Nucleus
              corf=true;
            }
            else if(tot4Mom.e() > 0 && cH4Mom.e() > 0 && nHadr > 1) // TeV error check
            {
              G4QHadron* prevHadr = theQHadrons[nHadr-2]; // GetPointer to Prev to theLast
              G4LorentzVector pH4Mom = prevHadr->Get4Momentum(); // 4mom of thePrevHadron
              G4double cHM = curHadr->GetMass();  // Mass of the current hadron
              G4double pHM = prevHadr->GetMass(); // Mass of the previous hadron
#ifdef pdebug
		  G4cout<<"G4QE::FSI:Bt4M="<<tot4Mom<<",c4M="<<cH4Mom<<",p4M="<<pH4Mom<<G4endl;
#endif
              G4LorentzVector tt4Mom=tot4Mom+cH4Mom+pH4Mom;
              G4double totRM=tt4Mom.m();
#ifdef pdebug
		  G4cout<<"G4QE::FS:"<<tt4Mom<<",tM="<<totRM<<",cM="<<cHM<<",pM="<<pHM<<G4endl;
#endif
              if(cHM+pHM<=totRM)                  // *** Make the final correction ***
	        {
                if(!G4QHadron(tt4Mom).DecayIn2(pH4Mom,cH4Mom))
                {
                  G4cout<<"***G4QE::FSI:**Correction**tot4M="<<tt4Mom<<totRM<<">sM="
                        <<cHM+cHM<<G4endl;
#ifdef pdebug
                  G4ExceptionDescription ed;
                  ed << "CORRECT DecIn2 error: **Correction**tot4M=" << tt4Mom
                     << totRM << ">sM=" << cHM+cHM << G4endl;
                  G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0007",
                              JustWarning, ed);
#endif
                }
#ifdef chdebug
                G4cout<<"-:-Warning-:-G4QE::FSI:***CORRECTION IS DONE*** d="<<dem<<G4endl;
#endif
                cor4M=tot4Mom;                                 // Complete correction
                curHadr->Set4Momentum(cH4Mom);
                prevHadr->Set4Momentum(pH4Mom);
                corf=true;
              }
              else
              {
#ifdef pdebug
                G4cerr<<"*!*G4QE::FSI: "<<cHM<<"+"<<pHM<<"="<<cHM+pHM<<">"<<totRM<<G4endl;
                G4ExceptionDescription ed;
                ed <<"TEMPORARY EXCEPTION: "<<cHM<<"+"<<pHM<<" = "<<cHM+pHM<<" > "<<totRM
                   <<", tot4M="<<tot4Mom<<", c4M="<<cH4Mom<<", p4M="<<pH4Mom<< G4endl;
                G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0008", 
                            JustWarning, ed);
#endif 
              }
            }
            else
            {
              G4cerr<<"*!*G4QE::FSI: tE="<<tot4Mom.e()<<", nHadr="<<nHadr<<G4endl;
              G4ExceptionDescription ed;
              ed << "TEMPORARY EXCEPTION: *check energy!* tot4M=" << tot4Mom << ", c4M="
                 << cH4Mom << ", nHadr="<< nHadr << " > 1 ?" << G4endl;
              G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0009", 
                          JustWarning, ed);
            }
            tot4Mom=tot4Mom-cor4M;
#ifdef pdebug
            G4cout<<"---Warning---G4QE::FSI:En/MomCons.Error is corrected:"<<cor4M<<G4endl;
#endif
          }
          if(nHadr>2 && !corf)
          {
            G4double cHM = curHadr->GetMass();  // Mass of the current hadron
            G4int ch=0;
            for(ch=nHadr-3; ch>-1; --ch)
            {
              G4QHadron* prevHadr = theQHadrons[ch]; // GetPointer to Hadr prev to theLast
              G4LorentzVector pH4Mom = prevHadr->Get4Momentum();// 4M of thePreviousHadron
              G4double pHM = prevHadr->GetMass(); // Mass of the current hadron
              tot4Mom+=cH4Mom+pH4Mom;
              G4double totRM=tot4Mom.m();
              if(cHM+pHM<=totRM)                  // *** Make the final correction ***
              {
                if(!G4QHadron(tot4Mom).DecayIn2(pH4Mom,cH4Mom))
                {
                  G4cout<<"***G4QEnv::FSI:**Correction**,tot4M="<<tot4Mom<<totRM<<" > sM="
                        <<cHM+cHM<<G4endl;
#ifdef debug
                  G4ExceptionDescription ed;
                  ed << "CORRECTION DecIn2Error: **Correction**,tot4M=" << tot4Mom
                     << totRM << " > sM=" << cHM+cHM << G4endl;
                  G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0010",
                              FatalException, ed);
#endif
                }
#ifdef chdebug
                G4cout<<"-:-!!!-:-G4QE::FSI:***CORRECTION IS DONE*** d="<<dem<<G4endl;
#endif
                curHadr->Set4Momentum(cH4Mom);
                prevHadr->Set4Momentum(pH4Mom);
                break;                            // Get out of the correction LOOP
              }
              else tot4Mom-=cH4Mom+pH4Mom;
            }
#ifdef ppdebug
            if(ch<0)
            {
              G4ExceptionDescription ed;
              ed << "EnMomCorrectionFailed: EnergyMomentumCorrection FAILED "
                 << G4endl;
              G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0011",
                          FatalException, ed);
            }
#endif
          } // End of additional attempt to correct EM by other hadrons
        }
#ifdef debug
        else  G4cout<<"...G4QE::FSI:E/M conservation is good enough"<<G4endl;
        G4cout<<"G4QE::FSI:EMCorrection by "<<theQHadrons.size()-nHadr<<" hadrons"<<G4endl;
#endif
        break;
      }
    } // End of the Back fusion LOOP
    // >| 2     | 2  | 2     | 2     | 2      | 2 - old    | 1. If gamma: add to sum4Mom
    //  |>0->sum| 3  | 3     | 3     | 3      | 3 - old    | 2. Compare BN with the Last
    //  | 5     |>5  | 4     | 4     | 4      | 4 - old    | 3. Swap if larger, del theLast
    //  | 0     | 0  |>0->sum| 5<-sum| 5->evap| 2 - new    | 4. If the Last: add the sum
    //  | 4     | 4  | 5     | ex    |        |(0 - gamma?)| 5. Decay/Eporate the Last
    //  | 3     | ex |                        | 3 - new
#ifdef chdebug
    // *** (3) Charge Control Sum Calculation for the Charge Conservation Check ***
    ccContSum=0;                              // Intermediate ChargeControlSum 
    cbContSum=0;                              // Intermediate BaryonNumberControlSum 
    if(nHadr)for(unsigned ic3=0; ic3<nHadr; ic3++) if(!(theQHadrons[ic3]->GetNFragments()))
    {
      ccContSum+=theQHadrons[ic3]->GetCharge();
      cbContSum+=theQHadrons[ic3]->GetBaryonNumber();
    }
    if(ccContSum-chContSum || cbContSum-bnContSum)
    {
      G4cout<<"*::*G4QE::FSI:(9) dC="<<ccContSum-chContSum<<",dB="<<cbContSum-bnContSum
            <<G4endl;
      //throw G4QException("G4QEnvironment::FSInteract: (3) Charge is not conserved");
    }
    // ***
#endif
    G4LorentzVector sum(0.,0.,0.,0.);
    G4int gamCount=0;
    nHadr=theQHadrons.size();
    G4bool frag=false;
    if(nHadr>2)for(unsigned f=0; f<theQHadrons.size(); f++) //Check that there's a fragment
    {
      G4int fBN=theQHadrons[f]->GetBaryonNumber(); // Baryon number of the fragment
#ifdef debug
      G4int fPDG=theQHadrons[f]->GetPDGCode();     // PDG code of the fragment
      G4LorentzVector fLV=theQHadrons[f]->Get4Momentum(); // 4Mom of the fragment
      G4cout<<"G4QE::FSI:"<<f<<",PDG="<<fPDG<<",fBN="<<fBN<<",f4M="<<fLV<<G4endl;
#endif
      if(fBN>1)                                    // At least one fragment (A>1) is found
      {
        frag=true;
        break;
      }
    }
#ifdef debug
    G4cout<<"G4QE::FSI:===Before Gamma Compression===, nH="<<nHadr<<",frag="<<frag<<G4endl;
#endif
    if(nHadr>2 && frag) for(G4int h=nHadr-1; h>=0; h--)//Collect gammas & kill DecayedHadrs
    {
      G4QHadron* curHadr = theQHadrons[h];    // Get a pointer to the current Hadron
      G4int   hF = curHadr->GetNFragments();
      G4int hPDG = curHadr->GetPDGCode();
      if(hPDG==89999003||hPDG==90002999)
        G4cout<<"---Warning---G4QEnv::FSI:nD-/pD++(1)="<<hPDG<<G4endl;
#ifdef debug
      G4cout<<"G4QE::FSI: h#"<<h<<", hPDG="<<hPDG<<", hNFrag="<<hF<<G4endl;
#endif
      if(hF||hPDG==22)                        // It should be compressed
      {
        G4QHadron* theLast = theQHadrons[theQHadrons.size()-1];//Get Ptr to the Last Hadron
        if(hPDG==22)
        {
          G4LorentzVector g4M=curHadr->Get4Momentum();
          sum+=g4M;                           // Add 4Mom of gamma to the "sum"
          gamCount++;
#ifdef debug
          G4cout<<"G4QE::FSI: gam4M="<<g4M<<" is added to s4M="<<sum<<G4endl;
#endif
        }
        nHadr = theQHadrons.size()-1;
        if(h < static_cast<G4int>(nHadr))     // Need swap with the Last
        {
          curHadr->SetNFragments(0);
          curHadr->Set4Momentum(theLast->Get4Momentum());
          G4QPDGCode lQP=theLast->GetQPDG();    // The QPDG of the last
          if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP); //CurHadr instead of LastHadr
          else curHadr->SetQC(theLast->GetQC());// CurHadrPDG instead of LastHadrPDG
#ifdef debug
          G4cout<<"G4QE::FSI: Exchange with the last is done"<<G4endl;
#endif
        }
        theQHadrons.pop_back();               // theLastQHadron is excluded from QHadrons
        delete theLast;//!!When kill,DON'T forget to delete theLastQHadron as an instance!!
#ifdef debug
        G4cout<<"G4QE::FSI: The last is compessed"<<G4endl;
#endif
      }
    }
#ifdef debug
    G4cout<<"G4QE::FSI: nH="<<nHadr<<"="<<theQHadrons.size()<<", sum="<<sum<<G4endl;
#endif
#ifdef chdebug
    // *** (4) Charge Control Sum Calculation for the Charge Conservation Check ***
    ccContSum=0;                   // Intermediate ChargeControlSum 
    cbContSum=0;                   // Intermediate BaryonNumberControlSum 
    if(nHadr)for(unsigned ic4=0; ic4<nHadr; ic4++) if(!(theQHadrons[ic4]->GetNFragments()))
    {
        ccContSum+=theQHadrons[ic4]->GetCharge();
        cbContSum+=theQHadrons[ic4]->GetBaryonNumber();
    }
    if(ccContSum-chContSum || cbContSum-bnContSum)
    {
      G4cout<<"*::*G4QE::FSI:(A) dC="<<ccContSum-chContSum<<",dB="<<cbContSum-bnContSum
            <<G4endl;
      //throw G4QException("G4QEnvironment::FSInteract: (4) Charge is not conserved");
    }
    // ***
#endif
    if(nHadr>1)for(unsigned hdr=0; hdr<theQHadrons.size()-1; hdr++)//Ord:theBigestIsTheLast
    {
      G4QHadron* curHadr = theQHadrons[hdr];  // Get a pointer to the current Hadron
#ifdef debug
      G4cout<<"G4QE::FSI:ORD,h="<<hdr<<"<"<<nHadr<<",hPDG="<<curHadr->GetPDGCode()<<G4endl;
#endif
      G4QHadron* theLast = theQHadrons[theQHadrons.size()-1]; //Get Ptr to the Last Hadron
      G4int hB           = curHadr->GetBaryonNumber();
      G4int lB           = theLast->GetBaryonNumber();
#ifdef debug
      G4cout<<"G4QE::FSI:hBN="<<hB<<"<lBN="<<lB<<",lstPDG="<<theLast->GetPDGCode()<<G4endl;
#endif
      if(lB<hB)                               // Must be swapped
      {
        G4QPDGCode   hQPDG = curHadr->GetQPDG();
        G4LorentzVector h4m= curHadr->Get4Momentum();
        curHadr->Set4Momentum(theLast->Get4Momentum());
        G4QPDGCode lQP=theLast->GetQPDG();    // The QPDG of the last
        if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP); //CurHadr instead of LastHadr
        else curHadr->SetQC(theLast->GetQC());// CurHadrPDG instead of LastHadrPDG
        theLast->Set4Momentum(h4m);
        theLast->SetQPDG(hQPDG);
      }
    }
    nHadr=theQHadrons.size();
#ifdef chdebug
    // *** (5) Charge Control Sum Calculation for the Charge Conservation Check ***
    ccContSum=0;                              // Intermediate ChargeControlSum 
    cbContSum=0;                              // Intermediate BaryonNumberControlSum 
    if(nHadr)for(unsigned ic5=0; ic5<nHadr; ic5++) if(!(theQHadrons[ic5]->GetNFragments()))
    {
      ccContSum+=theQHadrons[ic5]->GetCharge();
      cbContSum+=theQHadrons[ic5]->GetBaryonNumber();
    }
    if(ccContSum-chContSum || cbContSum-bnContSum)
    {
      G4cout<<"*::*G4QE::FSI:(B) dC="<<ccContSum-chContSum<<",dB="<<cbContSum-bnContSum
            <<G4endl;
      //throw G4QException("G4QEnvironment::FSInteract: (5) Charge is not conserved");
    }
    // ***
#endif
    if(gamCount)
    {
      G4QHadron* theLast = theQHadrons[nHadr-1];// Get a pointer to the Last Hadron
      if(theLast->GetBaryonNumber()>0)        // "Absorb photons & evaporate/decay" case
      {
        G4QHadron* theNew  = new G4QHadron(theLast); // Make New Hadron of the Last Hadron
#ifdef ffdebug
        G4cout<<"G4QE::FSI:BeforeLastSub,n="<<nHadr<<",PDG="<<theNew->GetPDGCode()<<G4endl;
#endif
        theQHadrons.pop_back();               // the last QHadron is excluded from OUTPUT
        delete theLast;// *!When kill,DON'T forget to delete theLastQHadron as anInstance!*
        nHadr--;                              // TheLastHadron is only virtually exists now
        G4int newPDG=theNew->GetPDGCode();
        G4LorentzVector new4M=theNew->Get4Momentum(); // 4-mom of the fragment
#ifdef debug
        G4cout<<"G4QE::FSI:gSum4M="<<sum<<" is added to "<<new4M<<", PDG="<<newPDG<<G4endl;
#endif
        G4LorentzVector exRes4M=new4M+sum;    //Icrease 4Mom of theLast by "sum of gammas"
        G4QNucleus exResidN(exRes4M,newPDG);
        //G4double mGamEva=2700.;             // @@Threshold for the evaporation
        G4double mGamEva=1700.;               // @@Threshold for the evaporation
        if(exResidN.SplitBaryon())
        //if(2>3)                             //CloseTheFirstPriorityResN+gamSumEvaporation
        {
          theNew->Set4Momentum(exRes4M);      // Icrease 4Mom of theLast by "sum" to Evapor
#ifdef ffdebug
          G4cout<<"G4QE::FSI:BeforeE(1),n="<<nHadr<<",nPDG="<<theNew->GetPDGCode()<<G4endl;
#endif
          EvaporateResidual(theNew);          // Try to evaporate the Nucl.(@@DecDib)(d.e.)
        }
        else if(theNew->GetPDGCode()==90002002&&exRes4M.m()>mHe3+mNeut&&G4UniformRand()>.5)
        {
          theNew->Set4Momentum(exRes4M);      // Icrease 4Mom of theLast by "sum" to Evapor
          G4LorentzVector n4M(0.,0.,0.,mNeut);
          G4LorentzVector h4M(0.,0.,0.,mHe3);
          if(!theNew->DecayIn2(n4M,h4M))
          {
            G4ExceptionDescription ed;
            ed << "GamSUPPRES DecIn2(n+He3)error: GamSup, tM=" << exRes4M.m()
               << "<n+He3=" << mNeut+mHe3 << G4endl;
            G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0012",
                        FatalException, ed);
          }
#ifdef ffdebug
          G4cout<<"G4QE::FSI:Gamma Suppression succided, n="<<n4M<<", He3="<<h4M<<G4endl;
#endif
          theNew->Set4Momentum(n4M);
          theNew->SetQPDG(nQPDG);             // convert the alpha to the neutron
          theQHadrons.push_back(theNew);      // (delete equivalent for theHad=neutron)
          G4QHadron* theHe3 = new G4QHadron(90002001,h4M);// Make a New Hadr for the He3
          theQHadrons.push_back(theHe3);      // (delete equivalent for the proton)
        }
        else if(nHadr)                        // Get LastHadrBefResNuc, absorb gam & decay
        //else if(2>3)                        // Close the pair absorbtion of gamma
        {
          if(nHadr>1)for(unsigned sh=0; sh<theQHadrons.size()-1; sh++)//Ord:MinE is TheLast
          {
            G4QHadron* curHadr = theQHadrons[sh];// Get a pointer to the current Hadron
            G4QHadron* thePrev = theQHadrons[theQHadrons.size()-1]; //GetPtr to theLastHadr
            G4LorentzVector h4M= curHadr->Get4Momentum();
            G4LorentzVector l4M= thePrev->Get4Momentum();
#ifdef ffdebug
            G4cout<<"G4QE::FSI:SO,h="<<sh<<"<"<<nHadr<<",PDG/LV="<<curHadr->GetPDGCode()
                  <<h4M<<G4endl;
#endif
            G4double hM=h4M.m();
            G4double hT=h4M.e()-hM;
            G4double lT=l4M.e()-l4M.m();
#ifdef ffdebug
            G4cout<<"G4QE::FSI:hT="<<hT<<"<T="<<lT<<".PDG="<<thePrev->GetPDGCode()<<G4endl;
#endif
            if(hM>mGamEva&&lT>hT)             // Must be swapped as the current is smaller
            {
              G4QPDGCode   hQPDG = curHadr->GetQPDG();
              curHadr->Set4Momentum(l4M);
              G4QPDGCode lQP=thePrev->GetQPDG();    // The QPDG of the previous
              if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP); //CurHadr instead of PrevHadr
              else curHadr->SetQC(thePrev->GetQC());// CurHadrPDG instead of PrevHadrPDG
              thePrev->Set4Momentum(h4M);
              thePrev->SetQPDG(hQPDG);
            }
          }
          nHadr=theQHadrons.size();
          G4QHadron* thePrev = theQHadrons[nHadr-1]; // GetPtr to the BeforeResidNuclHadron
          if(thePrev->Get4Momentum().m()>mGamEva)
          {
            G4QHadron* theHad  = new G4QHadron(thePrev);// MakeNewHadr of theBeforeResNuclH
#ifdef ffdebug
            G4cout<<"G4QE::FSI:BeforeResidNucHadr nH="<<nHadr<<",hPDG="
                  <<theHad->GetPDGCode()<<G4endl;
#endif
            theQHadrons.pop_back();           // theLastQHadron excluded from OUTPUT
            delete thePrev;// *!When kill,DON'T forget to delete theLastQHadrAsAnInstance!*
            G4LorentzVector n4M=theNew->Get4Momentum();// 4Mom of theLast (biggest nucleus)
            G4LorentzVector h4M=theHad->Get4Momentum();// 4Mom of the previous Hadron in HV
            G4LorentzVector dh4M=exRes4M+h4M; // 4Mom of LH+PH+sum(gam) for theDecay
            G4double dhM=dh4M.m();            // M of LH+PH+sum(gammas) for theDecay
            if(theHad->GetPDGCode()==90001001&&dhM>n4M.m()+mProt+mNeut&&G4UniformRand()>.5)
            //if(2>3)                           // Close Possibility toSplitDeuteron
            {
              G4double nuM=n4M.m();
              h4M=G4LorentzVector(0.,0.,0.,mNeut);
              G4LorentzVector p4M(0.,0.,0.,mProt);
              G4double sum_value=nuM+mNeut+mProt;
              if(fabs(dhM-sum_value)<eps)
              {
                n4M=dh4M*(nuM/sum_value);
                h4M=dh4M*(mNeut/sum_value);
                p4M=dh4M*(mProt/sum_value);
              }
              else if(dhM<sum_value || !G4QHadron(dh4M).DecayIn3(n4M,h4M,p4M))
              {
                G4ExceptionDescription ed;
                ed << "Gamma SUPPRESSION by D DecIn3error: GamSupByD,M="
                   << dhM << "<A+p+n=" << sum_value << G4endl;
                G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0013",
                            FatalException, ed);
              }
#ifdef ffdebug
              G4cout<<"G4QE::FSI:GamSuppression by d succided,h="<<h4M<<",A="<<n4M<<G4endl;
#endif
              theHad->Set4Momentum(h4M);
              theHad->SetQPDG(nQPDG);         // convert the deuteron to the neutron
              theQHadrons.push_back(theHad);  // (delete equivalent for theHad=neutron)
              G4QHadron* theProt = new G4QHadron(90001000,p4M);// Make NewHadr for Proton
              theQHadrons.push_back(theProt); // (delete equivalent for the proton)
              theNew->Set4Momentum(n4M);
              EvaporateResidual(theNew);      // TryToEvaporate theResNuc onceMore(del.eq.)
            }
            else
            {
              if(!G4QHadron(dh4M).DecayIn2(n4M,h4M))
              {
                G4ExceptionDescription ed;
                ed << "GamSUPPRESSION (3) DecIn2 error: GamSup,M=" << dh4M.m()
                   << "<A+h=" << n4M.m()+h4M.m() << G4endl;
                G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0014",
                            FatalException, ed);
              }
#ifdef ffdebug
              G4cout<<"G4QE::FSI:Gamma Suppression succided, h="<<h4M<<", A="<<n4M<<G4endl;
#endif
              theHad->Set4Momentum(h4M);
              theQHadrons.push_back(theHad);  // (delete equivalent for theHad)
              theNew->Set4Momentum(n4M);
              EvaporateResidual(theNew);      // Try to evaporate theResNuc (del.eq.)
            }
          }
          else
          {
            theNew->Set4Momentum(exRes4M);    // Icrease 4MomOfTheLast by "sum" for Evapor
#ifdef ffdebug
            G4cout<<"G4QE::FSI:BeforE(2),n="<<nHadr<<",PDG="<<theNew->GetPDGCode()<<G4endl;
#endif
            EvaporateResidual(theNew);  // Try to evaporate the Nucl.(@@DecDib)(delete eq.)
          }
        }
        else                            // Absorb gammas to theResidNucleus and evaporateIt
        {
          theNew->Set4Momentum(exRes4M);// Icrease 4Mom of the Last by the "sum" for Evap
          EvaporateResidual(theNew);    // Try to evaporate the Nucl.(@@DecDib)(delete eq.)
#ifdef ffdebug
          G4cout<<"G4QE::FSI:Bef.E(3),n="<<nHadr<<",PDG="<<newPDG<<",4M="<<exRes4M<<G4endl;
          unsigned nHN=theQHadrons.size();
          G4cout<<"G4QE::FSI:AfterEvaporation: nNew="<<nHN<<G4endl;
          if(nHN>nHadr)for(unsigned idp=nHadr; idp<nHN; idp++)
          G4cout<<"G4QE::FSI: h#"<<idp<<", PDG="<<theQHadrons[idp]->GetPDGCode()<<G4endl;
#endif
        }
        //G4int onH=nHadr;
        nHadr=theQHadrons.size();
        //if(nHadr>onH) bfAct=true;
      } // End of "the last is the nucleus" case
    } // End of "There are gammas to assimilate"
  } // End of the While-LOOOP for the Back Fusion
  tot4Mom=tmp4Mom; // Recover tot4Mom after the "En/Mom conservation" reduction
  // Final attempt to alpha-decay the residual nucleus, suppressing the gamma ===
  G4int gamcnt=0; // Counter of the residual gammas at this level
  nHadr=theQHadrons.size();
  unsigned maxB=nHadr-1;
#ifdef chdebug
  // *** (6) Charge Control Sum Calculation for the Charge Conservation Check ***
  ccContSum=0;                   // Intermediate ChargeControlSum 
  cbContSum=0;                   // Intermediate BaryonNumberControlSum 
  if(nHadr)for(unsigned ic6=0; ic6<nHadr; ic6++) if(!(theQHadrons[ic6]->GetNFragments()))
  {
    ccContSum+=theQHadrons[ic6]->GetCharge();
    cbContSum+=theQHadrons[ic6]->GetBaryonNumber();
  }
  if(ccContSum-chContSum || cbContSum-bnContSum)
  {
    G4cout<<"*::*G4QE::FSI:(C) dC="<<ccContSum-chContSum<<",dB="<<cbContSum-bnContSum
          <<G4endl;
    //throw G4QException("G4QEnvironment::FSInteract: (6) Charge is not conserved");
  }
  // ***
#endif
  lHadr=theQHadrons[maxB]->GetBaryonNumber();
  G4int tHadr=lHadr;                             // Total Baryon number
  if(nHadr>1)for(unsigned ipo=0; ipo<theQHadrons.size()-1; ipo++) // Find BiggestNuclFragm
  {
    G4int hPDG = theQHadrons[ipo]->GetPDGCode();
    if(hPDG==22) gamcnt++;
    else
    {
      G4int hBN  = theQHadrons[ipo]->GetBaryonNumber();
      tHadr+=hBN;
#ifdef debug
      G4cout<<"G4QE::FSI:h#"<<ipo<<":hPDG="<<hPDG<<",hBN="<<hBN<<",nH="<<theQHadrons.size()
            <<G4endl;
#endif
      if(hBN>lHadr)
      {
        lHadr=hBN;
        maxB=ipo;
      }                                           // the current biggest nuclear fragment
    }
  }
#ifdef debug
  G4cout<<"G4QE::FSI:max#"<<maxB<<",lB="<<lHadr<<",tBN="<<tHadr<<",gam="<<gamcnt<<G4endl;
#endif
  nHadr=theQHadrons.size();
#ifdef chdebug
  // *** (7) Charge Control Sum Calculation for the Charge Conservation Check ***
  ccContSum=0;                   // Intermediate ChargeControlSum
  cbContSum=0;                   // Intermediate BaryonNumberControlSum
  if(nHadr)for(unsigned ic7=0; ic7<nHadr; ic7++) if(!(theQHadrons[ic7]->GetNFragments()))
  {
    ccContSum+=theQHadrons[ic7]->GetCharge();
    cbContSum+=theQHadrons[ic7]->GetBaryonNumber();
  }
  if(ccContSum-chContSum || cbContSum-bnContSum)
  {
    G4cout<<"*::*G4QE::FSI:(D) dC="<<ccContSum-chContSum<<",dB="<<cbContSum-bnContSum
          <<G4endl;
    //throw G4QException("G4QEnvironment::FSInteract: (7) Charge is not conserved");
  }
  // ***
#endif
  if(gamcnt&&tHadr>1)                           // Only if there are gammas one should act
  {
    if(maxB+1<nHadr)                            // If maxB<Last, swap theCurH and theLastH
    {
      G4QHadron* theCurr = theQHadrons[maxB];   // Pointer to the Current Hadron
      G4QHadron* theLast = theQHadrons[nHadr-1];// Pointer to the Last Hadron
      G4QHadron* curHadr = new G4QHadron(theCurr);//Remember theCurrentHadron to put on top
      G4QPDGCode lQP=theLast->GetQPDG();        // The QPDG of the last
      if(lQP.GetPDGCode()!=10) theCurr->SetQPDG(lQP); //CurHadr instead of PrevHadr
      else theCurr->SetQC(theLast->GetQC());    // CurHadrPDG instead of LastHadrPDG
      theCurr->Set4Momentum(theLast->Get4Momentum()); // ... continue substitution
      theQHadrons.pop_back();                   // Rnt to theLastHadron is excluded from HV
      delete theLast;// *!!When kill,DON'T forget to delete the last QHadron as anInst. !!*
      theQHadrons.push_back(curHadr);           // The current Hadron, which is the Biggest
    }
    nHadr=theQHadrons.size();                   // Must be the same
    // Now it is necessary to absorb the photon (photons) and try to radiate alpha or evap.
    G4LorentzVector gamSum(0.,0.,0.,0.);
    if(nHadr>1)for(unsigned gp=0; gp<nHadr-1; gp++)// Find Gumma, remember and kill
    {
      G4QHadron* theCurr = theQHadrons[gp];       // Pointer to the Current Hadron
      G4int hPDG=theCurr->GetPDGCode();
#ifdef debug
      G4cout<<"G4QE::FSI:gp#"<<gp<<", PDG="<<hPDG<<", is found"<<G4endl;
#endif
      if(hPDG==22)                                // Photon is foun ond the "gp" position
      {
        gamSum=gamSum+theCurr->Get4Momentum();    // Accumulate the 4Momenta of the photon
#ifdef debug
        G4cout<<"G4QE::FSI:Photon gp#"<<gp<<",nH="<<nHadr<<", update gS="<<gamSum<<G4endl;
#endif
        unsigned nLast=nHadr-1;                   // position of theLastHadron (gp<nHadr-1)
        G4QHadron* theLast = theQHadrons[nLast];  // Pointer to the Last Hadron
#ifdef debug
        G4int wcn=0;
#endif
        while(nLast>=gp && theLast->GetPDGCode()==22) // "TheLast is a photon too" LOOP
        {
#ifdef debug
          ++wcn;
#endif
          if(nLast>gp)
          {
            gamSum=gamSum+theLast->Get4Momentum();// Accumulate 4-momentum of theLastPhoton
#ifdef debug
            G4cout<<"G4QE::FSI:TheLastPhotonIsFound #"<<wcn<<",update gS="<<gamSum<<G4endl;
#endif
          }
          theQHadrons.pop_back();                 // Pnt to theLastHadr.is excluded from HV
          delete theLast;// *!!When kill,DON'T forget to delete theLastQHadron as anInst!!*
          nHadr=theQHadrons.size();
          nLast=nHadr-1;
          theLast = theQHadrons[nLast];
        }
        if(nLast>gp)                              // -> swapping with the last
        {
          G4QPDGCode lQP=theLast->GetQPDG();      // The QPDG of the last
          if(lQP.GetPDGCode()!=10) theCurr->SetQPDG(lQP); //CurHadr instead of PrevHadr
          else theCurr->SetQC(theLast->GetQC());  // CurHadrPDG instead of LastHadrPDG
          theCurr->Set4Momentum(theLast->Get4Momentum()); // ... continue substitution
          theQHadrons.pop_back();                 // Pnt to theLastHadr.is excluded from HV
          delete theLast;// *!|When kill,DON'T forget to delete theLastQHadron as anInst!!*
          nHadr=theQHadrons.size();
#ifdef debug
          G4cout<<"G4QE::FSI:RepBy lPDG="<<lQP<<", nH="<<nHadr<<", gS="<<gamSum<<G4endl;
#endif
        }
      }
    }
    // @@ Now it is necessary to try to emit alpha or evaporate the residual nucleus
    G4QHadron* theLast = theQHadrons[nHadr-1];   // Pointer to the Last Hadron
    if(theLast->GetPDGCode()==22)
    {
      gamSum=gamSum+theLast->Get4Momentum();     // Accumulate 4-momentum of the LastPhoton
      theQHadrons.pop_back();                    // Pnt to theLastHadr.is excluded from HV
      delete theLast; // *!!When kill,DON'T forget to delete theLastQHadron as an inst.!!*
      nHadr=theQHadrons.size();
#ifdef debug
      G4cout<<"-Warning-G4QE::FSI: LastPhotonIsKilled, nH="<<nHadr<<",gS="<<gamSum<<G4endl;
#endif
      theLast = theQHadrons[nHadr-1];
    }
    G4int nEx=nHadr-2;                           // position to be exchanged with theLast
    while(theLast->GetBaryonNumber()<1 && nEx>=0)// theLastHadron must be a nucleus (A>0)
    {
      G4QHadron* theEx=theQHadrons[nEx];         // A hadron to be exchanged with theLast
      G4LorentzVector ex4Mom=theEx->Get4Momentum();
      G4QPDGCode exQPDG=theEx->GetQPDG();
      G4QContent exQC=theEx->GetQC();
      G4QPDGCode lQP=theLast->GetQPDG();         // The QPDG of the last
      if(lQP.GetPDGCode()!=10) theEx->SetQPDG(lQP); //CurHadr instead of PrevHadr
      else theEx->SetQC(theLast->GetQC());       // CurHadrPDG instead of LastHadrPDG
      theEx->Set4Momentum(theLast->Get4Momentum());
      if(exQPDG.GetPDGCode()!=10) theLast->SetQPDG(exQPDG);
      else theLast->SetQC(exQC);                 // CurHadrPDG instead of LastHadrPDG
      theLast->Set4Momentum(ex4Mom);
      nEx--;
    }
    G4QHadron* curHadr = new G4QHadron(theLast); // Pnt to theCurrentHadron is theLastCopy
    theQHadrons.pop_back();                 // Pnt to theLastHadron is excluded from OUTPUT
    delete theLast;// *!!When kill,DON'T forget to delete the LastQHadron as an instance!!*
    G4int theLB= curHadr->GetBaryonNumber();
    G4LorentzVector tR4M=curHadr->Get4Momentum()+gamSum;
    G4double tRM=tR4M.m();                       // TotMass of theResidualNucleus to decay
    if(theLB>4)
    {
      G4QContent lrQC=curHadr->GetQC()-G4QContent(6,6,0,0,0,0);
      G4QNucleus lrN(lrQC);
      G4double lrM=lrN.GetMZNS();
      if(tRM>lrM+mAlph)
      {
        G4LorentzVector lr4M(0.,0.,0.,lrM);
        G4LorentzVector al4M(0.,0.,0.,mAlph);
        if(!G4QHadron(tR4M).DecayIn2(lr4M,al4M))
        {
          curHadr->Set4Momentum(tR4M);
          EvaporateResidual(curHadr); // delete equivalent
#ifdef fdebug
          G4cout<<"G4QE::FSI: After Evap (1) nH="<<theQHadrons.size()<<G4endl;
#endif
        }
        else
        {
          delete curHadr;
          G4int APDG=lrN.GetPDG();
#ifdef debug
          G4cout<<"G4QE::FSI: Final A+alpha, A="<<APDG<<lr4M<<", a="<<al4M<<G4endl;
#endif
          G4QHadron* lrH = new G4QHadron(APDG,lr4M);
          theQHadrons.push_back(lrH);      // (delete equivalent for lrH)
          G4QHadron* alH = new G4QHadron(90002002,al4M);
          theQHadrons.push_back(alH);      // (delete equivalent for alH)
        }
      }
      else
      {
        curHadr->Set4Momentum(tR4M);
        EvaporateResidual(curHadr); // delete equivalent
#ifdef fdebug
        G4cout<<"G4QE::FSI: After Evap (2) nH="<<theQHadrons.size()<<G4endl;
#endif
      }
    }
    else
    {
      curHadr->Set4Momentum(tR4M);
      EvaporateResidual(curHadr); // delete equivalent
#ifdef fdebug
      G4cout<<"G4QE::FSI: After Evap (5) nH="<<theQHadrons.size()<<G4endl;
#endif
    }
  }
  //Now just fill the output theFravment vector (User is responsible to ClearAndDestroy it)
  nHadr=theQHadrons.size();
#ifdef chdebug
  // *** (8) Charge Control Sum Calculation for the Charge Conservation Check ***
  ccContSum=0;                   // Intermediate ChargeControlSum 
  cbContSum=0;                   // Intermediate BaryonNumberControlSum
  if(nHadr)for(unsigned ic8=0; ic8<nHadr; ic8++) if(!(theQHadrons[ic8]->GetNFragments()))
  {
    ccContSum+=theQHadrons[ic8]->GetCharge();
    cbContSum+=theQHadrons[ic8]->GetBaryonNumber();
  }
  if(ccContSum-chContSum || cbContSum-bnContSum)
  {
    G4cout<<"*::*G4QE::FSI:(E) dC="<<ccContSum-chContSum<<",dB="<<cbContSum-bnContSum
          <<G4endl;
    //throw G4QException("G4QEnvironment::FSInteract: (8) Charge is not conserved");
  }
  // ***
#endif
  if(nHadr) for(unsigned hd=0; hd<theQHadrons.size(); hd++)
  {
    //G4QHadron* curHadr = new G4QHadron(theQHadrons[hd]);
    G4QHadron* curHadr = theQHadrons[hd];
    G4int hPDG=curHadr->GetPDGCode();
    if(hPDG==22 && fabs(curHadr->Get4Momentum().e())<.00001) // E=0, gamma in the OUTPUT
    {
      unsigned lin=theQHadrons.size()-1;
      G4QHadron* theLast = theQHadrons[lin];// Pointer to theLastHadron in theQHadrVector
      if(lin>hd)
      {
        G4QPDGCode lQP=theLast->GetQPDG();         // The QPDG of the last
        if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP); //CurHadr instead of PrevHadr
        else curHadr->SetQC(theLast->GetQC());       // CurHadrPDG instead of LastHadrPDG
        curHadr->Set4Momentum(theLast->Get4Momentum()); // ... continue substitution (4Mom)
      }
      //else break; // It's not necessary to delete: not copy to theFragments is enough
      delete theLast;// *!!When kill, DON'T forget to delete theLastQHadron asAnInstance!!*
      theQHadrons.pop_back();              //pointer to theLastHadron is excluded from HV
      hPDG=curHadr->GetPDGCode();          // Redefine thePDGCode of theCurrentHadron
      if(lin==hd) break;
    }
#ifdef fdebug
    G4cout<<"G4QE::FSI: LOOP starts nH="<<nHadr<<", h#"<<hd<<", PDG="<<hPDG<<G4endl;
#endif
    if(hPDG==89999003||hPDG==90002999) G4cout<<"G4QEnv::FSI:nD-/pD++(0)="<<hPDG<<G4endl;
    if(hPDG==89999004||hPDG==90003999) G4cout<<"G4QEnv::FSI:nnD-/ppD++(0)="<<hPDG<<G4endl;
    //if(hPDG==89998003||hPDG==90002998)G4cout<<"G4QE::FSI:D-Pi-/D++Pi+PDG="<<hPDG<<G4endl;
    if(hPDG>100000000) G4cout<<"***G4QE::FSI: h#"<<hd<<", wrong PDGCode="<<hPDG<<G4endl;
#ifdef fdebug
    G4cout<<"G4QE::FSI:Copy is made with PDG="<<hPDG<<G4endl;//To compare with Original
#endif
    G4int hBN=curHadr->GetBaryonNumber();
    G4int hCG=curHadr->GetCharge();
    G4int hST=curHadr->GetStrangeness();
    G4int hNF=curHadr->GetNFragments();
    G4bool fOK=true;
    if(hNF==-1) curHadr->SetNFragments(0);
    else if(hPDG==91000000) curHadr->SetQPDG(G4QPDGCode(3122));// Move to theNextLoop
    else if(hPDG==90998003||hPDG==91002998)    // n,pi-,Sigma- OR p,pi+,Sigma+ (StrangeIso)
    {                                          // @@ can be converted here to B>2
#ifdef fdebug
      G4cout<<"G4QE::FSI: Pi+Nuc+Sigma state decay is found PDG="<<hPDG<<G4endl;
#endif
      G4LorentzVector r4M=curHadr->Get4Momentum(); // Real 4-mom of theIsoNucleus
      G4double reM=r4M.m();                   // Real mass of the IsoNucleus
      G4bool dub=false;
      G4int PDGnu=2112;
      G4double mNucl=mNeut;
      G4int PDGpi=-211;
      G4double mPion=mPi;
      G4int PDGsi=3112;
      G4double mSi=mSigM;
      G4double sum=mNucl+mPion+mSi;
      if(hPDG==90998003&&reM<sum)                      // Default -- values
      {
        PDGsi=2112;
        mSi=mNeut;
        mPion=mPi+mPi;
        sum=mNucl+mPion+mSi;
        dub=true;
      }
      if(hPDG==91002998)                      // Change -- default values to +++ values
      {
        mNucl=mProt;
        PDGnu=2212;
        PDGpi=211;
        PDGsi=3222;
        mSi  =mSigP;
        sum=mNucl+mPion+mSi;
        if(reM<sum)
        {
          PDGsi=2212;
          mSi=mProt;
          sum=mNucl+mPion+mSi;
        }
      }
      G4LorentzVector n4M(0.,0.,0.,mNucl);
      G4LorentzVector p4M(0.,0.,0.,mPion);
      G4LorentzVector k4M(0.,0.,0.,mSi);
      if(fabs(reM-sum)<eps)
      {
        //G4cout<<"G4QE::FSI:*TMP* PiDelta split PDG="<<hPDG<<G4endl;//To find out why?
        n4M=r4M*(mNucl/sum);
        p4M=r4M*(mPion/sum);
        k4M=r4M*(mSi/sum);
      }
      else if(reM<sum || !G4QHadron(r4M).DecayIn3(n4M,p4M,k4M))
      {
        G4cout<<"---Warning---G4QE::FSI:Pi+N+Sigma recovery INPDG="<<hPDG<<","<<reM<<" < "
              <<mNucl<<"(PDG="<<PDGnu<<")+Pi="<<mPion<<")+Sigma="<<mSi<<"="<<sum<<G4endl;
        if(!theEnvironment.GetA())
        {
          G4QHadron* theLast = curHadr;          // Prototype of Pointer to theLastHadron
          G4QHadron* qH = new G4QHadron(curHadr);// Copy of the Current Hadron
          if(hd+1<theQHadrons.size())            // If ipo<Last, swap CurHadr & LastHadr
          {
            theLast = theQHadrons[theQHadrons.size()-1]; // Pointer to LastHadr(ipo<Last)
            G4QPDGCode lQP=theLast->GetQPDG();         // The QPDG of the last
            if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP); //CurHadr instead of PrevHadr
            else curHadr->SetQC(theLast->GetQC());     // CurHadrPDG instead of LastHadrPDG
            curHadr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
          }
          theQHadrons.pop_back();        // exclude theLastHadron pointer from the OUTPUT
          delete theLast; // *!! When killing, DON'T forget to delete the last QHadron !!*
          G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());//Fake Quasmon
          if(!CheckGroundState(quasH,true))      // Try to correct by other hadrons
          {
            qH->SetNFragments(-1);              //@@ To avoid looping
            G4cout<<"---Warning---G4QE::FSI:PiNSig Failed, LeaveAsItIs 4m="<<r4M<<G4endl;
            theQHadrons.push_back(qH);           // Leave as it is (delete equivalent)
          }
          else
          {
#ifdef fdebug
            for(unsigned hp=0; hp<theQHadrons.size(); hp++)
            {
              G4QHadron* cpHadr = new G4QHadron(theQHadrons[hp]);
              G4int hpPDG=cpHadr->GetPDGCode();
              G4LorentzVector hpLV=cpHadr->Get4Momentum();
        G4cout<<"G4QE::FSI:h#"<<hp<<": hPDG="<<hpPDG<<", h4M="<<hpLV<<G4endl;
            }
#endif
            delete qH;
            nHadr=theQHadrons.size();
          }
          delete quasH;
          fOK=false;
        }
        else
        {
          G4ExceptionDescription ed;
          ed << "PiNucSigma Final Decay Error: No Final PiNSig recovery, Env="
             << theEnvironment << G4endl; 
          G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0015",
                      FatalException, ed);
        }
#ifdef fdebug
        G4cout<<"G4QEnv::FSI: PiNSi recover #"<<hd<<",PDG="<<curHadr->GetPDGCode()<<G4endl;
#endif
      }
      if(fOK)
      {
#ifdef fdebug
        G4cout<<"G4QE::FSI:PiNSigma==>"<<r4M<<"->N="<<PDGnu<<n4M<<"+Pi="<<PDGpi<<p4M
              <<"+ Sigma="<<PDGsi<<k4M<<G4endl;
#endif
        curHadr->Set4Momentum(n4M);
        curHadr->SetQPDG(G4QPDGCode(PDGnu));// At least one nucleun always exists
        if(dub)
        {
          p4M/=2.;
          G4QHadron* Pid = new G4QHadron(PDGpi,p4M); // Create a Hadron for the pion
          theQHadrons.push_back(Pid);       // Fill pion (delete equivalent)
        }
        G4QHadron* Pi = new G4QHadron(PDGpi,p4M); // Create a Hadron for the pion
        theQHadrons.push_back(Pi);          // Fill pion (delete equivalent)
        G4QHadron* Si = new G4QHadron(PDGsi,k4M); // Create a Hadron for the Sigma
        theQHadrons.push_back(Si);          // Fill Sigma (delete equivalent)
      }
#ifdef fdebug
      G4cout<<"G4QE::FSI:*TMP* PiNSigma end up PDG="<<hPDG<<G4endl;//To find out why?
#endif
    }
    else if(hPDG==89998003||hPDG==90002998)   // Isonucleus (pi- + DEL- or pi+ + DEL++)
    {                                         // @@ can be converted here to B>1
      //G4cout<<"G4QE::FSI:*TMP* PiDelta decay PDG="<<hPDG<<G4endl;//To find out why?
      G4double mNucl=mNeut;
      G4int PDGnu=2112;
      G4int PDGpi=-211;
      if(hPDG==90002998)                      // Change DEL- default values to DEL++ values
      {
        mNucl=mProt;
        PDGnu=2212;
        PDGpi=211;
      }
      //G4double nucM=mNucl*hBN;
      G4LorentzVector r4M=curHadr->Get4Momentum(); // Real 4-mom of theIsoNucleus
      G4double reM=r4M.m();                   // Real mass of the IsoNucleus
      G4LorentzVector n4M(0.,0.,0.,mNucl);
      G4LorentzVector p4M(0.,0.,0.,mPi);
      G4LorentzVector k4M(0.,0.,0.,mPi);
      G4double sum=mNucl+mPi+mPi;
      if(fabs(reM-sum)<eps)
      {
        //G4cout<<"G4QE::FSI:*TMP* PiDelta split PDG="<<hPDG<<G4endl;//To find out why?
        n4M=r4M*(mNucl/sum);
        p4M=r4M*(mPi/sum);
        k4M=r4M*(mPi/sum);
      }
      else if(reM<sum || !G4QHadron(r4M).DecayIn3(n4M,p4M,k4M))
      {
        G4cout<<"---Warning---G4QE::FSI: Isonuc+Pi recovery INPDG="<<hPDG<<","<<reM<<" < "
              <<mNucl<<"(PDG="<<PDGnu<<") + 2*"<<mPi<<"="<<sum<<G4endl;
        if(!theEnvironment.GetA())
        {
          G4QHadron* theLast = curHadr;          // Prototype of Pointer to theLastHadron
          G4QHadron* qH = new G4QHadron(curHadr);// Copy of the Current Hadron
          if(hd+1<theQHadrons.size())            // If ipo<Last, swap CurHadr & LastHadr
          {
            theLast = theQHadrons[theQHadrons.size()-1]; // Pointer to LastHadr(ipo<Last)
            G4QPDGCode lQP=theLast->GetQPDG();         // The QPDG of the last
            if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP); //CurHadr instead of PrevHadr
            else curHadr->SetQC(theLast->GetQC());     // CurHadrPDG instead of LastHadrPDG
            curHadr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
          }
          theQHadrons.pop_back();        // exclude theLastHadron pointer from the OUTPUT
          delete theLast; // *!! When killing, DON'T forget to delete the last QHadron !!*
          G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());//Fake Quasmon
          if(!CheckGroundState(quasH,true))      // Try to correct by other hadrons
          {
            qH->SetNFragments(-1);              //@@ To avoid looping
            G4cout<<"---Warning---G4QE::FSI:IsoN+Pi Failed, LeaveAsItIs 4m="<<r4M<<G4endl;
            theQHadrons.push_back(qH);           // Leave as it is (delete equivalent)
          }
          else
          {
            delete qH;
            nHadr=theQHadrons.size();
          }
          delete quasH;
          fOK=false;
        }
        else
        {
          G4ExceptionDescription ed;
          ed << "IsoNucl+Pi FinalDecayError: No Final IsoN+Pi recovery, Env="
             << theEnvironment << G4endl;
          G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0016",
                      FatalException, ed);
        }
#ifdef fdebug
        G4cout<<"G4QEnv::FSI: PiDel recover #"<<hd<<",PDG="<<curHadr->GetPDGCode()<<G4endl;
#endif
      }
      if(fOK)
      {
#ifdef fdebug
        G4cout<<"G4QE::FSI:IsoNuc+Pi==>"<<r4M<<"->N="<<PDGnu<<n4M<<"+Pi="<<PDGpi<<p4M
              <<"+ Pi="<<PDGpi<<k4M<<G4endl;
#endif
        curHadr->Set4Momentum(n4M);
        curHadr->SetQPDG(G4QPDGCode(PDGnu)); // At least one nucleun always exists
        G4QHadron* Pi1 = new G4QHadron(PDGpi,p4M); // Create a Hadron for the pion
        theQHadrons.push_back(Pi1);          // Fill pion (delete equivalent)
        G4QHadron* Pi2 = new G4QHadron(PDGpi,k4M); // Create a Hadron for the pion
        theQHadrons.push_back(Pi2);          // Fill pion (delete equivalent)
      }
#ifdef fdebug
      G4cout<<"G4QE::FSI:*TMP* PiDelta end up PDG="<<hPDG<<G4endl;//To find out why?
#endif
    }
    else if(hBN>0 && !hST && (hCG<0||hCG>hBN)) // Isonucleus (n*N+k*DEL- or n*P+k*DEL++)
    {
      G4double mNucl=mNeut;
      G4int PDGnu=2112;
      G4int PDGpi=-211;
      G4int nPi=-hCG;                         // Prototype of the minimum number of pions
      if(hCG>0)                               // Change DEL- default values to DEL++ values
      {
        mNucl=mProt;
        PDGnu=2212;
        PDGpi=211;
        nPi=hCG-hBN;
      }
      G4double nucM=mNucl*hBN;
      G4double pioM=mPi*nPi;
      G4LorentzVector r4M=curHadr->Get4Momentum(); // Real 4-mom of theIsoNucleus
      G4double reM=r4M.m();                   // Real mass of the IsoNucleus
      G4LorentzVector n4M(0.,0.,0.,nucM);
      G4LorentzVector k4M(0.,0.,0.,pioM);
      G4double sum=nucM+pioM;
      if(fabs(reM-sum)<eps)
      {
        n4M=r4M*(nucM/sum);
        k4M=r4M*(pioM/sum);
      }
      else if(reM<sum || !G4QHadron(r4M).DecayIn2(n4M,k4M))
      {
#ifdef fdebug
        G4cout<<"---Warning---G4QE::FSI: Isonucleus recovery INPDG="<<hPDG<<", M="<<reM
              <<" < "<<nucM<<"+"<<pioM<<"="<<sum<<G4endl;
#endif
        if(!theEnvironment.GetA())
        {
          G4QHadron* theLast = curHadr;          // Prototype of Pointer to theLastHadron
          G4QHadron* qH = new G4QHadron(curHadr);// Copy of the Current Hadron
          G4QContent tQC=qH->GetQC();            // Quark content of the hadron
          G4LorentzVector t4M=qH->Get4Momentum();// 4Momentum of the hadron
          if(hd+1<theQHadrons.size())            // If ipo<Last, swap CurHadr & LastHadr
          {
            theLast = theQHadrons[theQHadrons.size()-1]; // Pointer to LastHadr(ipo<Last)
            G4QPDGCode lQP=theLast->GetQPDG();         // The QPDG of the last
            if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP); //CurHadr instead of PrevHadr
            else curHadr->SetQC(theLast->GetQC());     // CurHadrPDG instead of LastHadrPDG
            curHadr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
          }
          theQHadrons.pop_back();        // exclude theLastHadron pointer from the OUTPUT
          delete theLast; // *!! When killing, DON'T forget to delete the last QHadron !!*
          G4Quasmon* quasH = new G4Quasmon(tQC,t4M); // Fake Quasmon for the Recovery
          if(!CheckGroundState(quasH,true))      // Try to correct by other hadrons
          {
            G4int tPDG=qH->GetPDGCode();
            qH->SetNFragments(-1);              // @@ To avoid looping
            G4cout<<"---Warning---G4QE::FSI:IsoN="<<tPDG<<tQC<<" FAsIs 4m="<<t4M<<G4endl;
            theQHadrons.push_back(qH);           // Leave as it is (delete equivalent)
          }
          else
          {
#ifdef fdebug
            for(unsigned hp=0; hp<theQHadrons.size(); hp++)
            {
              G4QHadron* cpHadr = new G4QHadron(theQHadrons[hp]);
              G4int hpPDG=cpHadr->GetPDGCode();
              G4LorentzVector hpLV=cpHadr->Get4Momentum();
              G4cout<<"G4QE::FSI:h#"<<hp<<": hPDG="<<hpPDG<<", h4M="<<hpLV<<G4endl;
            }
#endif
            delete qH;                           // Temporary Hadron is used for recovery
            nHadr=theQHadrons.size();
          }
          delete quasH;                          // TemporaryQuasmon is used for recovery
          fOK=false;
        }
        else
        {
          G4ExceptionDescription ed;
          ed << "IsoNucleus FinalDecayError: No FinalIsoNucRecovery, Env="
             << theEnvironment << G4endl;
          G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0017",
                      FatalException, ed);
        }
#ifdef fdebug
        G4cout<<"G4QEnv::FSI: Isonucleus recovery outPDG="<<curHadr->GetPDGCode()<<G4endl;
#endif
      }
      if(fOK)
      {
#ifdef fdebug
        G4cout<<"G4QE::FSI:IsoN==>"<<r4M<<"->N="<<PDGnu<<n4M<<"+Pi="<<PDGpi<<k4M<<G4endl;
#endif
        if(hBN>1) n4M/=hBN;                    // Calculate the 4mom per one nucleon
        curHadr->Set4Momentum(n4M);
        curHadr->SetQPDG(G4QPDGCode(PDGnu));   // At least one nucleun always exists
        if(hBN>1) for(G4int ih=1; ih<hBN; ih++)
        {
          G4QHadron* Hi = new G4QHadron(PDGnu,n4M); // Create a Hadron for the baryon
          theQHadrons.push_back(Hi);    // (del.eq.- user is responsible for del)
          //theFragments->push_back(Hi);        // Fill nucleons (delete equivalent)
        }
        if(nPi>1) k4M/=nPi;
        for(G4int ip=0; ip<nPi; ip++)
        {
          G4QHadron* Hj = new G4QHadron(PDGpi,k4M); // Create a Hadron for the pion
          theQHadrons.push_back(Hj);          // Fill pion (delete equivalent)
        }
      }
#ifdef fdebug
      G4cout<<"G4QEnv::FSI: Isonucleus decay result h#="<<hd<<", outPDG="<<hPDG<<G4endl;
#endif
    }
      else if ( hBN > 1 && 
               ( (hBN == hCG && !hST) || 
                 (!hCG && !hST) || 
                 (!hCG && hST==hBN) ) ) //(n*P, n*N, or n*L)
      {
      // *** Temporary Correction *** (To find out where is the error)
      if(hPDG==90000003 && fabs(curHadr->Get4Momentum().m()-mNeut-mNeut)<.1) {
        hPDG=90000002;
        hBN=2;
        G4cout<<"--Warning--G4QEnv::FSI:3->2 neutrons conversion (***Check it***)"<<G4endl;
      }
      // End of the Temporary Correction
      if     (!hCG&&!hST)     hPDG=90000001;
      else if(hBN==hCG&&!hST) hPDG=90001000;
      else if(!hCG&&hST==hBN) hPDG=91000000;
      else {
        // throw G4QException("***G4QEnvironment::FSInteract: MultyDibaryon cant be here");
        G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0018",
                    FatalException, "MultyDibaryon cant be here"); 
      }
      G4LorentzVector newLV=(curHadr->Get4Momentum())/hBN;
      curHadr->Set4Momentum(newLV);
      curHadr->SetQPDG(G4QPDGCode(hPDG));
      for(G4int bd=1; bd<hBN; bd++)
      {
        G4QHadron* secHadr = new G4QHadron(curHadr);
        theQHadrons.push_back(secHadr);    // (del.eq.- user is responsible for del)
        //theFragments->push_back(secHadr);// (del.equiv. - user is responsible for that)
      }
      }
    else if(hST<0 && hBN>0) // AntistrangeNucleus (@@ see above, already done)
    {
      G4LorentzVector r4M=curHadr->Get4Momentum(); // Real 4-mom of theAntiStrangeNucleus
      G4double reM=r4M.m();              // Real mass of the antiStrange nucleus
      G4QContent nQC=curHadr->GetQC();   // QCont of the Antistrange nucleus
      G4QNucleus  newN0(nQC-K0QC);
      G4int  RPDG=newN0.GetPDG();
      G4double mR=newN0.GetMZNS();       // Residual mass for K0
      G4double mKaon=mK0;                // Prototype is mass of K0
      G4int    kPDG =311;                // Prototype is PDG of K0
      G4QNucleus  newNp(nQC-KpQC);
      G4double mp=newNp.GetMZNS();       // Residual mass for K+
      if(mp+mK<mR+mK0)                   // Select K+ as the minimum mass of products
      {
        mR=mp;
        RPDG=newNp.GetPDG();
        mKaon=mK;
        kPDG=321;
      }
      G4double sum=mR+mKaon;
      if(sum>reM)                        // for GEANT4 (Can not decay in kaon and residual)
      {
        if(kPDG==321)                    // *** Very seldom *** "K+->pi+ conversion"
        {
          kPDG=211;
          mKaon=mPi;
        }
        else                             // *** Very seldom *** "K0(S=-1)->pi- conversion"
        {
          kPDG=111;
          mKaon=mPi0;
        }
        sum=mR+mKaon;
      }
      G4LorentzVector n4M(0.,0.,0.,mR);
      G4LorentzVector k4M(0.,0.,0.,mKaon);
      if(fabs(reM-sum)<eps)
      {
        n4M=r4M*(mR/sum);
        k4M=r4M*(mKaon/sum);
      }
      else if(reM<sum || !G4QHadron(r4M).DecayIn2(n4M,k4M))
      {
#ifdef debug

        G4cout<<"---Warning---G4QE::FSI: Try to recover ASN="<<hPDG<<","<<reM<<"<"<<mR<<"+"
              <<mKaon<<"="<<sum<<G4endl;
#endif
        if(!theEnvironment.GetA())
        {
          G4QHadron* theLast = curHadr;          // Prototype of Pointer to theLastHadron
          G4QHadron* qH = new G4QHadron(curHadr);// Copy of the Current Hadron
          if(hd+1<theQHadrons.size())           // If ipo<Last, swap CurHadr & LastHadr
          {
            theLast = theQHadrons[theQHadrons.size()-1]; // Pointer to LastHadr(ipo<Last)
            G4QPDGCode lQP=theLast->GetQPDG();         // The QPDG of the last
            if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP); //CurHadr instead of PrevHadr
            else curHadr->SetQC(theLast->GetQC());     // CurHadrPDG instead of LastHadrPDG
            curHadr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
          }
          theQHadrons.pop_back();        // exclude theLastHadron pointer from the OUTPUT
          delete theLast; // *!! When killing, DON'T forget to delete the last QHadron !!*
          G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());//Fake Quasmon
          if(!CheckGroundState(quasH,true))      // Try to correct by other hadrons
          {
            qH->SetNFragments(-1);              //@@ To avoid looping
#ifdef debug
            G4cout<<"---Warning---G4QE::FSI:AntiStraN Failed LeaveAsItIs 4m="<<r4M<<G4endl;
#endif
            theQHadrons.push_back(qH);           // Leave as it is (delete equivalent)
          }
          else
          {
            delete qH;
            nHadr=theQHadrons.size();
          }
          delete quasH;
          fOK=false;
        }
        else
        {
          G4ExceptionDescription ed;
          ed << "AntistrangeNucleus decayError: No Final AntiSN recovery, E="
             << theEnvironment << G4endl;
          G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0019",
                      FatalException, ed);
        }
      }
      if(fOK)
      {
#ifdef fdebug
        G4cout<<"G4QE::FSI:AntiSN==>"<<r4M<<"->A="<<RPDG<<n4M<<"+K="<<kPDG<<k4M<<G4endl;
#endif
        curHadr->Set4Momentum(k4M);                 // Kaon
        curHadr->SetQPDG(G4QPDGCode(kPDG));
        G4QHadron* theRes = new G4QHadron(RPDG,n4M);// Make a New Hadr for the residual
        EvaporateResidual(theRes);                  // Try to evaporate theResNuc (del.eq.)
      }
    }
    //                      =Lambda(Sigma0)=      = Sigma+ =        = Sigma- =
    else if(hPDG>90500000 && hPDG!=91000000 && hPDG!=91000999 && hPDG!=90999001 &&
                             hPDG!=91999000 && hPDG!=91999999 && hPDG!=92998999 )
    //                        ==  Xi-  ==       ==  Xi0  ==      === Omega- ===
    {
#ifdef pdebug
      G4cout<<"***G4QEnvironment::FSI:*G4* Hypernucleus PDG="<<hPDG<<" must decay"<<G4endl;
#endif
      G4int nL=curHadr->GetStrangeness();      // A number of Lambdas in the Hypernucleus
      G4int nB=curHadr->GetBaryonNumber();     // Total Baryon Number of the Hypernucleus
      G4int nC=curHadr->GetCharge();           // Charge of the Hypernucleus
      G4int nSM=0;                             // A#0f unavoidable Sigma-
      G4int nSP=0;                             // A#0f unavoidable Sigma+
      G4int nXZ=0;                             // A#0f unavoidable Xi0
      G4int nXM=0;                             // A#0f unavoidable Xi-
      G4int nOM=0;                             // A#0f unavoidable Omega-
      // @@ Now it does not include more complicated clusters as Omega-,Xi- (Improve)
      if(nC < 0)                               // Negative hypernucleus
      {
        if(-nC <= nL)                          // Negative Charge is smaller than Strange
        {
          if(nL <= nB)                         // ==> Big Hyper-nucleus
          {
            nSM=-nC;                           // Can be compensated by Sigma-
            nL+=nC;                            // Reduce the residual strangeness
          }
          else                                 // ==> Lack of BaryonNumber (better use Xi-)
          {
            G4int nX=nL-nB;
            if(-nC <= nX && nL >= nX+nX)       // Part or all Xi are Xi- (e.g. Xi-,Lambda)
            {
              nXM=-nC;
              if(nXM != nX) nXZ=nX-nXM;        // The rest are Xi0
              nL-=nX+nX;                       // Xi reduce srangeness twice faster
            }
            else if(nL >= nX-nC)               // Sigma- should be used in addition to Xi-
            {                                  // e.g. Xi-, Sigma-
              nXM=nX;
              nSM=-nC-nX;
              nL-=nX-nC;
            }
            // @@ Don't close this warning as it costs days to find out this growing point!
            else G4cout<<"-Warning-G4QEn::FSI:*Improve*,-nC<=nL,nL>nB, PDG="<<hPDG<<G4endl;
          }
        }
        else                                   // -nC can not be totally compensated by nL
        {
          nSM=nL;                              // The maximumNumber of Sigma- (theRest pi-)
          nL=0;                                // Kill the residualStrangeness (@notEnough)
        }
      }
      else if(nC>nB-nL)                        // Extra positive hypernucleus (nC>=0 here)
	{                                        // Can't compensate theCharge by onlyProtons
        if(nC<=nB)
        {
          if(nL <= nB)
          {
            G4int dH=nB-nC;
            nSP=nL-dH;                         // Can be compensated by Sigma+
            nL=dH;                             // Reduce the residual strangeness
          }
          else if(nL<nB+nB)                    // Xi0 can be used
          {
            nXZ=nL-nB;                         // This is how many Xi0
            nL=nB-nL+nB;                       // a#of S=1 particles
            if(nC <= nL)
            {
              nSP=nC;
              nL-=nSP;
            }
            else
            {
              nSP=nL;
              nL=0;
            }
          }
          else if(nL > nB && ((nL-nB)/2)*2 == nL-nB && nC+nL <= nB+nB)// Omega- can be used
          {
            nOM=(nL-nB)/2;                     // This is how many Omega- can be made
            nL=nB+nB-nL-nC;                    // a#of Lambdas should be 0 !
            nSP=nOM-nC;
          }
          // @@ Do not close this warning as it costs days to find out this growing point !
          else G4cout<<"-Warning-G4QEn::FSI:*Improve*,nC>nB-nL,nC<=nB, PDG="<<hPDG<<G4endl;
        }
        else
        {
          nSP=nL;                              // The maximumNumber of Sigma+ (theRest pi+)
          nL=0;                                // Kill the residual strangeness
        }
      }
      // else = the charge can be compensated by nucleons
      G4LorentzVector r4M=curHadr->Get4Momentum(); // Real 4-momentum of the hypernucleus
      G4double reM=r4M.m();                    // Real mass of the hypernucleus
      // Subtract Lamb/Sig/Xi/Omega from the Nucleus and decay
      G4int SS=nL+nSP+nSM+nXZ+nXM;
      if(SS<nB && !nOM)                        // Should not be Xi's or Omeg in the nucleus
      {
        G4int rlPDG = hPDG-nL*1000000-nSP*1000999-nSM*999001;
        G4int    sPDG=3122;                    // Prototype for the Hyperon PDG (Lambda)
        G4double MLa=mLamb;                    // Prototype for one Hyperon decay
#ifdef pdebug
        G4cout<<"G4QEnvironment::FSI:*G4* nS+="<<nSP<<", nS-="<<nSM<<", nL="<<nL<<G4endl;
#endif
        if(nSP || nSM)
        {
          if(nL)
          {
            // Hopefully it never happens
            G4cout<<"-Warning-G4QE::FSI:*Improve*,PDG="<<hPDG<<": both Sigma&Lamb"<<G4endl;
            // @@ Bad Correction, which does not conserve the charge !! (-> add decay in 3)
            if(nSP) nL+=nSP;
            else    nL+=nSM;
          }
          if(nSP)
          {
            nL=nSP;
            sPDG=3222;
            MLa=mSigP;
          }
          else
          {
            nL=nSM;
            sPDG=3112;
            MLa=mSigM;
          }
        }
#ifdef pdebug
        G4cout<<"G4QEnvironment::FSI:*G4* mS="<<MLa<<", sPDG="<<sPDG<<", nL="<<nL<<G4endl;
#endif
        if(nL>1) MLa*=nL;
        G4double rlM=G4QNucleus(rlPDG).GetMZNS();// Mass of theNewResidualNonstrangeNucleus
        if(!nSP&&!nSM&&nL==1&&reM>rlM+mSigZ&&G4UniformRand()>.5)// Sigma0 @@ AddToHadroniz?
        {
          sPDG=3212;
          MLa=mSigZ;
        }
        G4int rnPDG = hPDG-nL*999999;          // Convert Lambdas to neutrons (for convInN)
        G4QNucleus rnN(rnPDG);                 // New nonstrange nucleus
        G4double rnM=rnN.GetMZNS();            // Mass of the new nonstrange nucleus
        // @@ Take into account that can be Iso-Hypernucleus (Add PI+,R & Pi-,R decays)
        if(rlPDG==90000000)                    // Multy Hyperon (HyperNuc of only hyperons)
        {
          if(nL>1) r4M=r4M/nL;                 // split the 4-mom for the MultyLambda
          for(G4int il=0; il<nL; il++)         // loop over Lambdas
          {
            G4QHadron* theLam = new G4QHadron(sPDG,r4M);// Make a New Hadr for the Hyperon
            theQHadrons.push_back(theLam);     // (del.eq.- user is responsible for del)
          }
        }
        else if(reM>rlM+MLa-eps)               // Lambda(or Sigma) can be split
        {
          G4LorentzVector n4M(0.,0.,0.,rlM);
          G4LorentzVector h4M(0.,0.,0.,MLa);
          G4double sum=rlM+MLa;
          if(fabs(reM-sum)<eps)
          {
            n4M=r4M*(rlM/sum);
            h4M=r4M*(MLa/sum);
          }
          else if(reM<sum || !G4QHadron(r4M).DecayIn2(n4M,h4M))
          {
            G4ExceptionDescription ed;
            ed << "Hypernuclear L-decay error: Hypern, M=" << reM << "<A+n*L="
               << sum << ",d=" << sum-reM << G4endl;
            G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0020",
                        FatalException, ed);
          }
#ifdef pdebug
          G4cout<<"*G4QE::FSI:HypN="<<r4M<<"->A="<<rlPDG<<n4M<<",n*Lamb="<<nL<<h4M<<G4endl;
#endif
          curHadr->Set4Momentum(n4M);
          curHadr->SetQPDG(G4QPDGCode(rlPDG)); // converted hypernucleus
          if(rlPDG==90000002)                  // Additional action with curHadr change
          {
            G4LorentzVector newLV=n4M/2.;
            curHadr->Set4Momentum(newLV);
            curHadr->SetQPDG(G4QPDGCode(90000001));
            G4QHadron* secHadr = new G4QHadron(curHadr);
            theQHadrons.push_back(secHadr);    // (del.eq.- user is responsible for del)
          }
          else if(rlPDG==90002000)             // Additional action with curHadr change
          {
            G4LorentzVector newLV=n4M/2.;
            curHadr->Set4Momentum(newLV);
            curHadr->SetQPDG(G4QPDGCode(90001000));
            G4QHadron* secHadr = new G4QHadron(curHadr);
            theQHadrons.push_back(secHadr);    // (del.eq.- user is responsible for del)
          }
          // @@(?) Add multybaryon decays if necessary (Now it anyhow is made later)
#ifdef pdebug
          G4cout<<"*G4QE::FSI: resNucPDG="<<curHadr->GetPDGCode()<<G4endl;
#endif
          if(nL>1) h4M=h4M/nL;                 // split the lambda's 4-mom if necessary
          for(G4int il=0; il<nL; il++)         // loop over excessive
          {
            G4QHadron* theLamb = new G4QHadron(sPDG,h4M);// Make a New Hadr for the Hyperon
            theQHadrons.push_back(theLamb);    // (del.eq.- user is responsible for del)
          }
        }
        else if(reM>rnM+mPi0-eps&&!nSP&&!nSM)  // Lambda->N only if Sigmas are absent
        {
          G4int nPi=static_cast<G4int>((reM-rnM)/mPi0);
          if (nPi>nL) nPi=nL;
          G4double npiM=nPi*mPi0;
          G4LorentzVector n4M(0.,0.,0.,rnM);
          G4LorentzVector h4M(0.,0.,0.,npiM);
          G4double sum=rnM+npiM;
          if(fabs(reM-sum)<eps)
          {
            n4M=r4M*(rnM/sum);
            h4M=r4M*(npiM/sum);
          }
          else if(reM<sum || !G4QHadron(r4M).DecayIn2(n4M,h4M))
          {
            G4ExceptionDescription ed;
            ed << "Hypernuclear decay error: HyperN,M=" << reM << "<A+n*Pi0="
               << sum << ",d=" << sum-reM << G4endl;
            G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0021",
                        FatalException, ed);
          }
          curHadr->Set4Momentum(n4M);
          curHadr->SetQPDG(G4QPDGCode(rnPDG)); // converted hypernucleus
#ifdef pdebug
          G4cout<<"*G4QE::FSI:HypN "<<r4M<<"->A="<<rnPDG<<n4M<<",n*Pi0="<<nPi<<h4M<<G4endl;
#endif
          if(nPi>1) h4M=h4M/nPi;               // split the 4-mom if necessary
          for(G4int ihn=0; ihn<nPi; ihn++)     // loop over additional pions
          {
            G4QHadron* thePion = new G4QHadron(111,h4M);// Make a New Hadr for the pion
            theQHadrons.push_back(thePion);    // (del.eq.- user is responsible for del)
            //theFragments->push_back(thePion);  // (delete equivalent for the pion)
          }
          if(rnPDG==90000002)                  // Additional action with curHadr change
          {
            G4LorentzVector newLV=n4M/2.;
            curHadr->Set4Momentum(newLV);
            curHadr->SetQPDG(G4QPDGCode(90000001));
            G4QHadron* secHadr = new G4QHadron(curHadr);
            theQHadrons.push_back(secHadr);    // (del.eq.- user is responsible for del)
            //theFragments->push_back(secHadr);  // (del.eq.- user is responsible for del)
          }
          else if(rnPDG==90002000)             // Additional action with curHadr change
          {
            G4LorentzVector newLV=n4M/2.;
            curHadr->Set4Momentum(newLV);
            curHadr->SetQPDG(G4QPDGCode(90001000));
            G4QHadron* secHadr = new G4QHadron(curHadr);
            theQHadrons.push_back(secHadr);    // (del.eq.- user is responsible for del)
            //theFragments->push_back(secHadr);  // (del.eq.- user is responsible for del)
          }
          // @@ Add multybaryon decays if necessary
        }
        else if(reM>rnM-eps)    // Decay in nonstrange and gamma
        {
          G4LorentzVector n4M(0.,0.,0.,rnM);
          G4LorentzVector h4M(0.,0.,0.,0.);
          G4double sum=rnM;
          if(fabs(reM-sum)<eps) n4M=r4M;
          else if(reM<sum || !G4QHadron(r4M).DecayIn2(n4M,h4M))
          {
            G4ExceptionDescription ed;
            ed << "Hypernuclear GammaDecay error: Hypern,M=" << reM << "<A+n*Pi0="
               << sum << ",d=" << sum-reM << G4endl;
            G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0022",
                        FatalException, ed);
          }
          curHadr->Set4Momentum(n4M);
          curHadr->SetQPDG(G4QPDGCode(rnPDG)); // converted hypernucleus
#ifdef pdebug
          G4cout<<"*G4QE::FSI:HypNg "<<r4M<<"->A="<<rnPDG<<n4M<<",gamma="<<h4M<<G4endl;
#endif
          G4QHadron* theGamma = new G4QHadron(22,h4M);// Make a New Hadr for the gamma
          theQHadrons.push_back(theGamma);     // (del.eq.- user is responsible for del)
          if(rnPDG==90000002)                  // Additional action with curHadr change
          {
            G4LorentzVector newLV=n4M/2.;
            curHadr->Set4Momentum(newLV);
            curHadr->SetQPDG(G4QPDGCode(90000001));
            G4QHadron* secHadr = new G4QHadron(curHadr);
            theQHadrons.push_back(secHadr);    // (del.eq.- user is responsible for del)
            //theFragments->push_back(secHadr);  // (del.eq.- user is responsible for del)
          }
          else if(rnPDG==90002000)             // Additional action with curHadr change
          {
            G4LorentzVector newLV=n4M/2.;
            curHadr->Set4Momentum(newLV);
            curHadr->SetQPDG(G4QPDGCode(90001000));
            G4QHadron* secHadr = new G4QHadron(curHadr);
            theQHadrons.push_back(secHadr);    // (del.eq.- user is responsible for del)
            //theFragments->push_back(secHadr);  // (del.eq.- user is responsible for del)
          }
          // @@ Add multybaryon decays if necessary
        }
#ifdef pdebug
        else // If this Error shows up (lowProbable appearance) => now it is left as is
        {
          G4cout<<"-Warning-G4QEnv::F:R="<<rlM<<",S+="<<nSP<<",S-="<<nSM<<",L="<<nL<<",d1="
                <<rlM+MLa-reM<<G4endl;
          G4cout<<"-Warning-G4QEnv::FSI:HyperN="<<hPDG<<", M="<<reM<<"<"<<rnM+mPi0<<", d2="
                <<rnM+mPi0-reM<<G4endl;
          //throw G4QException("G4QEnvironment::FSInteract: Hypernuclear conversion");
        }
#endif
      }
      else if(SS==nB && !nOM) // Decay with Xi, but without residual nucleus
      {
#ifdef pdebug
        G4cout<<"G4QEnvironment::FSI:OnlyHyp,B="<<nB<<",SP="<<nSP<<",SM="<<nSM<<",L="<<nL
              <<",XM="<<nXM<<",XZ="<<nXZ<<G4endl;
#endif
        G4int SP=0;
        if(nL) ++SP;
        if(nSP) ++SP;
        if(nSM) ++SP;
        if(nXZ) ++SP;
        if(nXM) ++SP;
        if(SP>1 && SP<4)
        {
          G4int fPDG=3122;
          G4double fM=mLamb;
          G4int fS=nL;
          G4int sPDG=3222;
          G4double sM=mSigP;
          G4int sS=nSP;
          G4int uPDG=3322;
          G4double uM=mXiZ;
          G4int uS=nSP;
          if(nL)
          {
            if(nSM)
            {
              sPDG=3112;
              sM  =mSigM;
              sS  =nSM;
              if(SP==3)
              {
                if(nXM)
                {
                  uPDG=3312;
                  uM  =mXiM;
                  uS  =nXM;
                }
              }
            }
            else if(nXZ)
            {
              sPDG=3322;
              sM  =mXiZ;
              sS  =nXZ;
              if(SP==3)
              {
                if(nSP)
                {
                  uPDG=3222;
                  uM  =mSigP;
                  uS  =nSP;
                }
                else if(nXM)
                {
                  uPDG=3312;
                  uM  =mXiM;
                  uS  =nXM;
                }
              }
            }
            else if(nXM)
            {
              sPDG=3312;
              sM  =mXiM;
              sS  =nXM;
              //if(SP==3) uXiZ is a default
            }
            else if(SP==3)                          // Lambda,Sigma+ (only Xi0 is possible)
            {
              if(nXZ)
              {
                uPDG=3322;
                uM  =mXiZ;
                uS  =nXZ;
              }
              else G4cout<<"-Warning-G4QE::FSI: *Improve* StrangeComb, PDG="<<hPDG<<G4endl;
            }
          }
          else if(nSM)
          {
            fPDG=3112;
            fM  =mSigM;
            fS  =nSM;
            if(nXZ)
            {
              sPDG=3322;
              sM  =mXiZ;
              sS  =nXZ;
              if(SP==3)
              {
                if(nXM)
                {
                  uPDG=3312;
                  uM  =mXiM;
                  uS  =nXM;
                }
                else if(nXZ)
                {
                  uPDG=3322;
                  uM  =mXiZ;
                  uS  =nXZ;
                }
              }
            }
            else if(nXM)
            {
              sPDG=3312;
              sM  =mXiM;
              sS  =nXM;
            }
          }
          else if(nSP)
          {
            fPDG=3222;
            fM  =mSigP;
            fS  =nSP;
            if(nXZ)
            {
              sPDG=3322;
              sM  =mXiZ;
              sS  =nXZ;
              if(SP==3)
              {
                if(nXZ)
                {
                  uPDG=3322;
                  uM  =mXiZ;
                  uS  =nXZ;
                }
              }
            }
            else if(nXM)
            {
              sPDG=3312;
              sM  =mXiM;
              sS  =nXM;
            }
          }
          else if(nXZ)
          {
            fPDG=3322;
            fM  =mXiZ;
            fS  =nXZ;
            if(nXM)
            {
              sPDG=3312;
              sM  =mXiM;
              sS  =nXM;
            }
          }
          else G4cout<<"-Warning-G4QE::FSI: *Improve* StrangeFragment, PDG="<<hPDG<<G4endl;
          // Now make two or three particles decay
          if(SP==2)                           // @@ Only 2BodyDecay is implemented >Improve
          {
#ifdef pdebug
            G4cout<<"G4QEnvironment::FSI:2HypD,fS="<<fS<<",fPDG="<<sPDG<<",fM="<<fM<<",sS="
                  <<sS<<",sPDG="<<sPDG<<",sM="<<sM<<",rM="<<reM<<G4endl;
#endif
            fM*=fS;
            sM*=sS;
            G4double mm_value=fM+sM;
            G4double MM=reM+eps;
            if(MM<=mm_value && fPDG==3122 && sPDG==3222) // Lamba,Sigma+ => Xi0,p (can be 50%)
            {
              fPDG= 2212;
              fM  = mProt;
              sPDG= 3322;
              sM  = mXiZ;
              mm_value  = fM+sM;
            }
            if(MM>mm_value)                    // can be split or decayed
            {
              G4LorentzVector f4M(0.,0.,0.,fM);
              G4LorentzVector s4M(0.,0.,0.,sM);
              G4double sum=fM+sM;
              if(fabs(reM-sum)<=eps)           // splitting
              {
                f4M=r4M*(fM/sum);
                s4M=r4M*(sM/sum);
              }
              else if(reM<sum || !G4QHadron(r4M).DecayIn2(f4M,s4M))
              {
                // G4cerr<<"***G4QE::FSI:Hyp2,M="<<reM<<"< sum="<<sum<<",d="<<sum-reM<<G4endl;
                // throw G4QException("***G4QEnvironment::FSInter: HypernucOnlyStran2 error");
                G4ExceptionDescription ed;
                ed << "HypernucOnlyStran2 error: Hyp2,M=" << reM << "< sum=" << sum
                   << ",d=" << sum-reM << G4endl;
                G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0023",
                            FatalException, ed);
              }
#ifdef pdebug
              G4cout<<"G4QEnv::FSI:2,HN="<<r4M<<hPDG<<"->First="<<fS<<f4M<<fPDG<<",Second="
                    <<sS<<s4M<<sPDG<<G4endl;
#endif
              if(fS>1)
              {
                f4M/=fS;                           // split the Hyperons 4-mom if necessary
                for(G4int il=1; il<fS; il++)       // loop over excessive Hyperons
                {
                  G4QHadron* theHyp = new G4QHadron(fPDG,f4M);// Make NewHadr for theHyper
                  theQHadrons.push_back(theHyp);    // (del.eq.- user is respons for del)
                }
              }
              curHadr->Set4Momentum(f4M);
              curHadr->SetQPDG(G4QPDGCode(fPDG));  // converted hypernucleus
              if(sS>1) s4M/=sS;                    // split the Hyperon 4-mom if necessary
              for(G4int il=0; il<sS; il++)         // loop over excessive
              {
                G4QHadron* theHyp = new G4QHadron(sPDG,s4M);// Make NewHadr for theHyperon
                theQHadrons.push_back(theHyp);     // (del.eq.- user is respons for del)
              }
            }
            // @@ Don't close this warning as it costs days to find out this growing point!
            else G4cout<<"-Warning-G4QE::FSI:*Improve* S2, PDG="<<hPDG<<",M="<<reM<<",fS="
                       <<fS<<",sS="<<sS<<",fM="<<fM<<",sM="<<sM<<G4endl;
          }
          else  // --> 3Body @@ Only 3BodyDecay is implemented >Improve
          {
#ifdef pdebug
            G4cout<<"G4QEnvironment::FSI:3HypD,fS="<<fS<<",fPDG="<<sPDG<<",fM="<<fM<<",sS="
                  <<sS<<",sPDG="<<sPDG<<",fM="<<fM<<",uS="<<sS<<",uPDG="<<sPDG<<",uM="<<sM
                  <<",rM="<<reM<<G4endl;
#endif
            fM*=fS;
            sM*=sS;
            uM*=sS;
            G4double mm_value=fM+sM+uM;
            G4double MM=reM+eps;
            // @@ There is much more alternative 3B states...
            //if(MM<=mm_value && fPDG==3122 && sPDG==3222) // Lamba,Sigma+ => Xi0,p (can be 50%)
            //{
            //  fPDG= 2212;
            //  fM  = mProt;
            //  sPDG= 3322;
            //  sM  = mXiZ;
            //  mm_value  = fM+sM;
            //}
            if(MM>mm_value)                        // can be split or decayed
            {
              G4LorentzVector f4M(0.,0.,0.,fM);
              G4LorentzVector s4M(0.,0.,0.,sM);
              G4LorentzVector u4M(0.,0.,0.,uM);
              G4double sum=fM+sM+uM;
              if(fabs(reM-sum)<=eps)               // splitting
              {
                f4M=r4M*(fM/sum);
                s4M=r4M*(sM/sum);
                u4M=r4M*(uM/sum);
              }
              else if(reM<sum || !G4QHadron(r4M).DecayIn2(f4M,s4M))
              {
                G4ExceptionDescription ed;
                ed << "HypernucOnlyStran3 error: Hyp3,M=" << reM << "< sum="
                   << sum << ",d=" << sum-reM << G4endl;
                G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0024",
                            FatalException, ed);
              }
#ifdef pdebug
              G4cout<<"G4QEnv::FSI:3,HN="<<r4M<<hPDG<<"->First="<<fS<<f4M<<fPDG<<",Second="
                    <<sS<<s4M<<sPDG<<",Third="<<fS<<f4M<<fPDG<<G4endl;
#endif
              if(fS>1)
              {
                f4M/=fS;                           // split the Hyperons 4-mom if necessary
                for(G4int il=1; il<fS; ++il)       // loop over excessive Hyperons
                {
                  G4QHadron* theHyp = new G4QHadron(fPDG,f4M);// Make NewHadr for theHyper
                  theQHadrons.push_back(theHyp);    // (del.eq.- user is respons for del)
                }
              }
              curHadr->Set4Momentum(f4M);          // The hypernucleus itself is modified
              curHadr->SetQPDG(G4QPDGCode(fPDG));  // converted hypernucleus
              if(sS>1) s4M/=sS;                    // split the Hyperon 4-mom if necessary
              for(G4int il=0; il<sS; ++il)         // loop over excessive
              {
                G4QHadron* theHyp = new G4QHadron(sPDG,s4M);// Make NewHadr for theHyperon
                theQHadrons.push_back(theHyp);     // (del.eq.- user is respons for del)
              }
              if(uS>1) u4M/=uS;                    // split the Hyperon 4-mom if necessary
              for(G4int il=0; il<uS; ++il)         // loop over excessive
              {
                G4QHadron* theHyp = new G4QHadron(uPDG,u4M);// Make NewHadr for theHyperon
                theQHadrons.push_back(theHyp);     // (del.eq.- user is respons for del)
              }
            }
            // @@ Don't close this warning as it costs days to find out this growing point!
            else G4cout<<"-Warn-G4QE::FSI:*Improve* S3, PDG="<<hPDG<<",M="<<reM<<",fS="<<fS
                      <<",sS="<<sS<<",uS="<<uS<<",fM="<<fM<<",sM="<<sM<<",uM="<<uM<<G4endl;
          }
        }
      }
      else if(nOM) // Decay with Omega-, but without residual nucleus
      {
#ifdef pdebug
        G4cout<<"G4QEnvir::FSI:Omega-,B="<<nB<<",SP="<<nSP<<",OM="<<nOM<<",L="<<nL<<G4endl;
#endif
        G4int SP=0;
        if(nL) ++SP;
        if(nSP) ++SP;
        if(nOM) ++SP;
        if(!nL)
        {
          G4int    fPDG= 3334;
          G4double fM  = mOmM;
          G4int    fS  = nOM;
          G4int    sPDG= 3222;
          G4double sM  = mSigP;
          G4int    sS  = nSP;
#ifdef pdebug
          G4cout<<"G4QEnvironment::FSI:2HypD,fS="<<fS<<",fPDG="<<sPDG<<",fM="<<fM<<",sS="
                <<sS<<",sPDG="<<sPDG<<",sM="<<sM<<",rM="<<reM<<G4endl;
#endif
          fM*=fS;
          sM+=sS;
          if(reM>fM+sM-eps)                 // can be split or decayed
          {
            G4LorentzVector f4M(0.,0.,0.,fM);
            G4LorentzVector s4M(0.,0.,0.,sM);
            G4double sum=fM+sM;
            if(fabs(reM-sum)<eps)           // splitting
            {
              f4M=r4M*(fM/sum);
              s4M=r4M*(sM/sum);
            }
            else if(reM<sum || !G4QHadron(r4M).DecayIn2(f4M,s4M))
            {
              G4ExceptionDescription ed;
              ed << "HypernucOnlyStran3 error: Hyp,M=" << reM << "<OM+SP="
                 << sum << ",d=" << sum-reM << G4endl;
              G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0025",
                          FatalException, ed);
            }
#ifdef pdebug
            G4cout<<"G4QEnv::FSI:OHNuc="<<r4M<<hPDG<<"->First="<<fS<<f4M<<fPDG<<",Second="
                  <<sS<<s4M<<sPDG<<G4endl;
#endif
            if(fS>1)
            {
              f4M/=fS;                           // split the Hyperons 4-mom if necessary
              for(G4int il=1; il<fS; il++)       // loop over excessive Hyperons
              {
                G4QHadron* theHyp = new G4QHadron(fPDG,f4M);// Make NewHadr for theHyper
                theQHadrons.push_back(theHyp);    // (del.eq.- user is respons for del)
              }
            }
            curHadr->Set4Momentum(f4M);
            curHadr->SetQPDG(G4QPDGCode(fPDG));  // converted hypernucleus
            if(sS>1) s4M/=sS;                    // split the Hyperon 4-mom if necessary
            for(G4int il=0; il<sS; il++)         // loop over excessive
            {
              G4QHadron* theHyp = new G4QHadron(sPDG,s4M);// Make NewHadr for theHyperon
              theQHadrons.push_back(theHyp);     // (del.eq.- user is respons for del)
            }
          }
          // @@ Don't close this warning as it costs days to find out this growing point!
          else G4cout<<"-Warning-G4QE::FSI:*Improve* OMSP, PDG="<<hPDG<<",M="<<reM<<",nSP="
                     <<nSP<<",nOM="<<nOM<<",mOM="<<fM<<",mSP="<<sM<<G4endl;
        }
        // @@ Do not close this warning as it costs days to find out this growing point !
        else G4cout<<"-Warning-G4QE::FSI:*Improve* nLamb="<<nL<<" # 0, PDG="<<hPDG<<G4endl;
      }
      // @@ Do not close this warning as it costs days to find out this growing point !
      else G4cout<<"-Warning-G4QE::FSI:*Improve*,S="<<SS<<">B="<<nB<<",PDG="<<hPDG<<G4endl;
    }
    if(!fOK) hd--;
  } // 
#ifdef pdebug
  G4cout<<"G4QE::FSI:==>OUT,nH="<<theQHadrons.size()<<",nF="<<theFragments->size()<<G4endl;
#endif
  // *** (Final) Charge Control Sum Calculation for the Charge Conservation Check ***
  //nHadr=theFragments->size();
  nHadr=theQHadrons.size();
  G4int cfContSum=0;                   // Final ChargeControlSum 
  G4int bfContSum=0;                   // Final BaryonNumberControlSum 
  if(nHadr)for(unsigned f=0; f<nHadr; f++)
  {
    G4QHadron* curHadr = new G4QHadron(theQHadrons[f]);
#ifdef debug
    G4cout<<"G4QE::FSI:#"<<f<<curHadr->Get4Momentum()<<curHadr->GetPDGCode()<<G4endl;
#endif
    if(!(curHadr->GetNFragments()))
    {
      G4int chg=curHadr->GetCharge();
      cfContSum+=chg;
      G4int brn=curHadr->GetBaryonNumber();
      bfContSum+=brn;
      G4int str=curHadr->GetStrangeness();
      if ( brn > 1 && 
           ( (!str && (chg == brn || !chg)) || 
             (!chg && str == brn) ) ) // Check for multibaryon(split)
      {
        G4int             bPDG=90000001; // Prototype: multineutron
        if     (chg==brn) bPDG=90001000; // Multyproton
        else if(str==brn) bPDG=91000000; // Multilambda
        G4LorentzVector sp4m=curHadr->Get4Momentum()/brn; // Split 4-mom of the Multibaryon
        curHadr->Set4Momentum(sp4m);
        curHadr->SetQPDG(G4QPDGCode(bPDG)); // Substitute Multibaryon by a baryon
        for (G4int ib=1; ib<brn; ib++)
        {
          G4QHadron* bH = new G4QHadron(bPDG, sp4m);
         theFragments->push_back(bH);     //(del eq. - user is responsible for deleting)
        }
      }
      theFragments->push_back(curHadr);  //(del eq. - user is responsible for deleting)
    }
  }
#ifdef chdebug
  if(cfContSum-chContSum || bfContSum-bnContSum)
  {
    G4ExceptionDescription ed;
    ed << "::(F) Charge is not conserved: Ch=" << cfContSum-chContSum << ",Bn="
       << bfContSum-bnContSum << G4endl;
    G4Exception("G4QEnvironment::FSInteraction()", "HAD_CHPS_0026",
                FatalException, ed);
  }
#endif
  // ***
  return theFragments;
} // End of "FSInteraction"

//The public Quasmons duplication with delete responsibility of User (!)
G4QuasmonVector* G4QEnvironment::GetQuasmons()
{
  G4int nQ=theQuasmons.size();
#ifdef debug
  G4cout<<"G4QEnvironment::GetQuasmons is called nQ="<<nQ<<G4endl;
#endif
  G4QuasmonVector* quasmons = new G4QuasmonVector;   // !! User is responsible to delet it
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
#ifdef debug
    G4cout<<"G4QEnv::GetQuasm:Q#"<<iq<<",QQPDG="<<theQuasmons[iq]->GetQPDG()<<",QQC="
          <<theQuasmons[iq]->GetQC()<<",M="<<theQuasmons[iq]->Get4Momentum().m()<<G4endl;
#endif
    G4Quasmon* curQ = new G4Quasmon(theQuasmons[iq]);
    quasmons->push_back(curQ);                 // (delete equivalent - user is responsible)
  }
#ifdef debug
  G4cout<<"G4QEnvironment::GetQuasmons ===OUT==="<<G4endl;
#endif
  return quasmons;
} // End of GetQuasmons

//The public QHadrons duplication with delete responsibility of User (!)
G4QHadronVector* G4QEnvironment::GetQHadrons()
{
  G4int nH=theQHadrons.size();
#ifdef debug
  G4cout<<"G4QEnvironment::GetQHadrons is called nH="<<nH<<G4endl;
#endif
  G4QHadronVector* hadrons = new G4QHadronVector;  // !! User is responsible to delet it
  if(nH) for(G4int ih=0; ih<nH; ih++)
  {
#ifdef debug
    G4cout<<"G4QEnv::GetQHadrons:H#"<<ih<<",HQPDG="<<theQHadrons[ih]->GetQPDG()<<",HQC="
          <<theQHadrons[ih]->GetQC()<<",HM="<<theQHadrons[ih]->GetMass()<<G4endl;
#endif
    G4QHadron* curH = new G4QHadron(theQHadrons[ih]);
    hadrons->push_back(curH);                       // (del. equiv. - user is responsibile)
  }
#ifdef debug
  G4cout<<"G4QEnvironment::GetQHadrons ===OUT=== Copied nQH="<<hadrons->size()<<G4endl;
#endif
  return hadrons;
} // End of GetQHadrons

//Public QHadrons cleaning up after extraction (GetQHadrons) between Construct and Fragment
void G4QEnvironment::CleanUpQHadrons()
{
  for_each(theQHadrons.begin(), theQHadrons.end(), DeleteQHadron());
  theQHadrons.clear();
} // End of CleanUpQHadrons

//The public FillQHadrons filling. It must be used only internally in CHIPS
void G4QEnvironment::FillQHadrons(G4QHadronVector* input)
{
  G4int nH=input->size();
#ifdef debug
  G4cout<<"G4QEnvironment::FillQHadrons is called nH="<<nH<<G4endl;
#endif
  if(nH) for(G4int ih=0; ih<nH; ih++)
  {
#ifdef debug
    G4cout<<"G4QEnv::FillQHadrons:H#"<<ih<<",HQPDG="<<(*input)[ih]->GetQPDG()<<",HQC="
          <<(*input)[ih]->GetQC()<<",HM="<<(*input)[ih]->GetMass()<<G4endl;
#endif
    G4QHadron* curH = new G4QHadron((*input)[ih]);
    theQHadrons.push_back(curH);                        // (del. equiv. 
  }
#ifdef debug
  G4cout<<"G4QEnvironment::FillQHadrons ===OUT=== Filled nQH="<<theQHadrons.size()<<G4endl;
#endif
} // End of FillQHadrons

//Decay of the excited Baryon in baryon & meson (gamma)
void G4QEnvironment::DecayBaryon(G4QHadron* qH, G4QHadronVector* HV)
{
  static const G4QPDGCode gQPDG(22);
  static const G4QPDGCode pizQPDG(111);
  static const G4QPDGCode pipQPDG(211);
  static const G4QPDGCode pimQPDG(-211);
  static const G4QPDGCode kmQPDG(-321);
  static const G4QPDGCode kzQPDG(-311);
  static const G4QPDGCode nQPDG(2112);
  static const G4QPDGCode pQPDG(2212);
  static const G4QPDGCode lQPDG(3122);
  static const G4QPDGCode laQPDG(3122);
  static const G4QPDGCode smQPDG(3112);
  static const G4QPDGCode szQPDG(3212);
  static const G4QPDGCode spQPDG(3222);
  static const G4QPDGCode kszQPDG(3322);
  static const G4QPDGCode ksmQPDG(3312);
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mPi0 = G4QPDGCode(111).GetMass();
  static const G4double mK   = G4QPDGCode(321).GetMass();
  static const G4double mK0  = G4QPDGCode(311).GetMass();
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mSigM= G4QPDGCode(3112).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mSigZ= G4QPDGCode(3212).GetMass();
  static const G4double mSigP= G4QPDGCode(3222).GetMass();
  //static const G4double mKsiM= G4QPDGCode(3312).GetMass();
  //static const G4double mKsi0= G4QPDGCode(3322).GetMass();
  //static const G4double mOmeg= G4QPDGCode(3334).GetMass();
  static const G4double mNPi0 = mPi0+ mNeut;
  static const G4double mNPi  = mPi + mNeut;
  static const G4double mPPi0 = mPi0+ mProt;
  static const G4double mPPi  = mPi + mProt;
  static const G4double mLPi0 = mPi0+ mLamb;
  static const G4double mLPi  = mPi + mLamb;
  static const G4double mSpPi = mPi + mSigP;
  static const G4double mSmPi = mPi + mSigM;
  static const G4double mPK   = mK  + mProt;
  static const G4double mPKZ  = mK0 + mProt;
  static const G4double mNKZ  = mK0 + mNeut;
  static const G4double mSpPi0= mPi0+ mSigP;
  static const G4double mSzPi0= mPi0+ mSigZ;
  static const G4double mSzPi = mPi + mSigZ;
  static const G4double eps  = 0.003;
  //static const G4QNucleus vacuum(90000000);
  G4int theLB= qH->GetBaryonNumber();          // Baryon number of the Baryon
  if(theLB!=1)
  {
    G4cerr<<"***G4QEnvironment::DecayBaryon: A!=1 -> fill as it is"<<G4endl;
#ifdef ppdebug
    // throw G4QException("G4QEnv::DecayBaryon: Unknown Baryon with A!=1");
    G4Exception("G4QEnvironment::DecayBaryon()", "HAD_CHPS_0000",
                FatalException, "Unknown Baryon with A!=1");
#endif
    HV->push_back(qH);                 // Fill AsIs (delete equivalent)
    return;
  }
  G4int theLC= qH->GetCharge();                // Chsrge of the Baryon
  G4int theLS= qH->GetStrangeness();           // Strangness of the Baryon
  //G4int          qPDG = qH->GetPDGCode();      // PDG Code of the decaying baryon
  G4LorentzVector q4M = qH->Get4Momentum();    // Get 4-momentum of the Baryon
  G4double         qM = q4M.m();               // Mass of the Baryon
#ifdef debug
  G4cout<<"G4QEnv::DecayBaryon: *Called* S="<<theLS<<",C="<<theLC<<",4M="<<q4M<<qM<<G4endl;
#endif
  // Select a chanel of the baryon decay
  G4QPDGCode     fQPDG = pQPDG;                 // Prototype for Proton
  G4double       fMass= mProt;
  G4QPDGCode     sQPDG = pizQPDG;                  // Prototype for Pi0
  G4double       sMass= mPi0;
  if(!theLS)                    // This is only for not strange baryons
  {
    if(!theLC)                  // Neutron like: n+gam, n+Pi0, p+Pi- are possible
    {
      if(qM<mNPi0)              // Only n+gamma decay is possible
      {
        if(qM < mNeut+eps)
        {
#ifdef debug
          G4cout<<"G4QEnv::DecayBaryon: Fill Neutron AsIs"<<G4endl;
#endif
          qH->SetQPDG(nQPDG);
          HV->push_back(qH); // Fill AsIs (delete equivalent)
          if(qM+eps<mNeut) G4cout<<"-Warning-G4QEnv::DecayBaryon: N0 AsIs, M="<<qM<<G4endl;
          return;
        }
        fQPDG=nQPDG;            // Baryon is neutron
        fMass=mNeut;
        sQPDG=gQPDG;            // Meson is gamma
        sMass=0.;
#ifdef debug
        G4cout<<"-Warning-G4QEnv::DecayBaryon: Gamma+n, M="<<qM<<",d="<<qM-mNeut<<G4endl;
#endif
      }
      else if(qM>mPPi)          // Both n+pi0 (p=2/3) & p+Pi- (p=1/3) are possible
      {
        if(G4UniformRand()>.333333333) // Clebsh value for the Isobar decay
        {
          fQPDG=nQPDG;          // Baryon is neutron (meson is Pi0)
          fMass=mNeut;
        }
        else
        {
          sQPDG=pimQPDG;        // Meson is Pi- (baryon is proton)
          sMass=mPi;
        }
      }
      else                      // Only n+Pi0 decay is possible
      {
        fQPDG=nQPDG;            // Baryon is neutron
        fMass=mNeut;
      }
    }
    else if(theLC==1)           // Proton like: p+gam, p+Pi0, n+Pi+ are possible
    {
      if(qM<mPPi0)              // Only p+gamma decay is possible
      {
        if(qM < mProt+eps)
        {
#ifdef debug
          G4cout<<"G4QEnv::DecayBaryon: Fill Proton AsIs"<<G4endl;
#endif
          qH->SetQPDG(pQPDG);
          HV->push_back(qH); // Fill AsIs (delete equivalent)
#ifdef debug
          if(qM+eps<mProt)G4cout<<"-Warning-G4QEnv::DecayBaryon: Pr+ AsIs, M="<<qM<<G4endl;
#endif
          return;
        }
        sQPDG=gQPDG;            // Meson is gamma (baryon is proton)
        sMass=0.;
      }
      else if(qM>mNPi)          // Both p+pi0 (p=2/3) & n+Pi+ (p=1/3) are possible
      {
        if(G4UniformRand()<.333333333) // Clebsh value for the Isobar decay
        {
          fQPDG=nQPDG;          // Baryon is neutron
          fMass=mNeut;
          sQPDG=pipQPDG;        // Meson is Pi+
          sMass=mPi;
        }
        // p+Pi0 is a default case
      }
      // p+Pi0 is a default case
    }
    else if(theLC==2)           // Delta++ like: only p+PiP is possible
    {
      if(qM>mPPi)               // Only p+gamma decay is possible
      {
        sQPDG=pipQPDG;          // Meson is Pi+ (baryon is proton)
        sMass=mPi;
      }
      else                      // @@ Can be aReason to search for anError in Fragmentation
      {
#ifdef debug
        G4cout<<"-Warning-G4QE::DecBary:*AsIs* DEL++ M="<<qM<<"<"<<mPPi<<G4endl;
#endif
        HV->push_back(qH);               // Fill AsIs (delete equivalent)
        return;
      }
    }
    else if(theLC==-1)          // Delta- like: only n+PiM is possible
    {
      if(qM+eps>mNPi)           // Only p+gamma decay is possible
      {
        fQPDG=nQPDG;            // Baryon is neutron
        fMass=mNeut;
        sQPDG=pimQPDG;          // Meson is Pi-
        sMass=mPi;
      }
      else                      // @@ Can be aReason to search for anError in Fragmentation
      {
#ifdef debug
        G4cout<<"-Warning-G4QE::DecBary:*AsIs* DEL++ M="<<qM<<"<"<<mNPi<<G4endl;
#endif
        HV->push_back(qH);               // Fill AsIs (delete equivalent)
        return;
      }
    }
    else 
    {
#ifdef debug
      G4cout<<"-Warning-G4QE::DecBary:*AsIs* UnknBaryon (S=0) QC="<<qH->GetQC()<<G4endl;
#endif
      HV->push_back(qH);                 // Fill AsIs (delete equivalent)
      return;
    }
  }
  else if(theLS==1)             // ==>->->-> S=1 <-<-<-<==
  {
    if(!theLC)                  // -->> Lambda/Sz: L+g,L+Pi0,Sz+Pi0,Sm+Pip,Sp+Pim,p+Km,n+Kz
    {
      if(qM<mLPi0)              // Only Lambda+gamma decay is possible
      {
        if(qM < mLamb+eps)
        {
#ifdef debug
          G4cout<<"G4QEnv::DecayBaryon: Fill Lambda AsIs"<<G4endl;
#endif
          qH->SetQPDG(lQPDG);
          HV->push_back(qH); // Fill AsIs (delete equivalent)
          if(qM+eps<mLamb)G4cout<<"-Warning-G4QEnv::DecayBaryon: La0 AsIs, M="<<qM<<G4endl;
          return;
        }
        fQPDG=lQPDG;            // Baryon is Lambda
        fMass=mLamb;
        sQPDG=gQPDG;            // Meson is gamma
        sMass=0.;
      }
      else if(qM<mSzPi0)        // Only Lambda+Pi0 is possible
      {
        fQPDG=lQPDG;            // Baryon is Lambda
        fMass=mLamb;
      }
      else if(qM<mSpPi)         // Both Lambda+Pi0 & Sigma0+Pi0 are possible
      {
        if(G4UniformRand()>.6)  // @@ Relative probability (take into account Phase Space)
        {
          fQPDG=szQPDG;         // Baryon is Sigma0
          fMass=mSigZ;
        }
        else
        {
          fQPDG=lQPDG;          // Baryon is Lambda
          fMass=mLamb;
        }
      }
      else if(qM<mSmPi)         // Lambda+Pi0, Sigma0+Pi0, & SigmaP+PiM are possible
      {
        G4double ra=G4UniformRand();
        if(ra<.4)               // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=lQPDG;          // Baryon is Lambda
          fMass=mLamb;
        }
        else if(ra<.7)          // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=szQPDG;         // Baryon is Sigma0
          fMass=mSigZ;
        }
        else
        {
          fQPDG=spQPDG;         // Baryon is SigmaP
          fMass=mSigP;
          sQPDG=pimQPDG;        // Meson is Pi-
          sMass=mPi;
        }
      }
      else if(qM<mPK)           // Lambda+Pi0, Sig0+Pi0, SigP+PiM, SigM+PiP are possible
      {
        G4double ra=G4UniformRand();
        if(ra<.35)              // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=lQPDG;          // Baryon is Lambda
          fMass=mLamb;
        }
        else if(ra<.6)          // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=szQPDG;         // Baryon is Sigma0
          fMass=mSigZ;
        }
        else if(ra<.8)          // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=smQPDG;         // Baryon is SigmaM
          fMass=mSigM;
          sQPDG=pipQPDG;        // Meson is Pi+
          sMass=mPi;
        }
        else
        {
          fQPDG=spQPDG;         // Baryon is SigmaP
          fMass=mSigP;
          sQPDG=pimQPDG;        // Meson is Pi-
          sMass=mPi;
        }
      }
      else if(qM<mNKZ)          // Lambda+Pi0, Sig0+Pi0, SigP+PiM, SigM+PiP are possible
      {
        G4double ra=G4UniformRand();
        if(ra<.3)               // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=lQPDG;          // Baryon is Lambda
          fMass=mLamb;
        }
        else if(ra<.5)          // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=szQPDG;         // Baryon is Sigma0
          fMass=mSigZ;
        }
        else if(ra<.7)          // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=smQPDG;         // Baryon is SigmaM
          fMass=mSigM;
          sQPDG=pipQPDG;        // Meson is Pi+
          sMass=mPi;
        }
        else if(ra<.9)          // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=spQPDG;         // Baryon is SigmaP
          fMass=mSigP;
          sQPDG=pimQPDG;        // Meson is Pi-
          sMass=mPi;
        }
        else
        {
          fQPDG=pQPDG;          // Baryon is Proton
          fMass=mProt;
          sQPDG=kmQPDG;         // Meson is K-
          sMass=mK;
        }
      }
      else                      // Only n+Pi0 decay is possible
      {
        G4double ra=G4UniformRand();
        if(ra<.3)               // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=lQPDG;          // Baryon is Lambda
          fMass=mLamb;
        }
        else if(ra<.5)          // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=szQPDG;         // Baryon is Sigma0
          fMass=mSigZ;
        }
        else if(ra<.65)         // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=smQPDG;         // Baryon is SigmaM
          fMass=mSigM;
          sQPDG=pipQPDG;        // Meson is Pi+
          sMass=mPi;
        }
        else if(ra<.8)          // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=spQPDG;         // Baryon is SigmaP
          fMass=mSigP;
          sQPDG=pimQPDG;        // Meson is Pi-
          sMass=mPi;
        }
        else if(ra<.9)          // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=pQPDG;          // Baryon is Proton
          fMass=mProt;
          sQPDG=kmQPDG;         // Meson is K-
          sMass=mK;
        }
        else
        {
          fQPDG=nQPDG;          // Baryon is Neutron
          fMass=mNeut;
          sQPDG=kzQPDG;         // Meson is K0
          sMass=mK0;
        }
      }
    }
    else if(theLC==1)           // SigmaP: SigP+gam,SigP+Pi0,Sp+PiP,L+Pi0,p+K0 are possible
    {
      if(qM<mLPi)               // Only SigmaPlus+gamma decay is possible
      {
        if(qM < mSigP+eps)
        {
#ifdef debug
          G4cout<<"G4QEnv::DecayBaryon: Fill SigmaPlus AsIs"<<G4endl;
#endif
          qH->SetQPDG(spQPDG);
          HV->push_back(qH); // Fill AsIs (delete equivalent)
          if(qM+eps<mSigP)G4cout<<"-Warning-G4QEnv::DecayBaryon: Si+ AsIs, M="<<qM<<G4endl;
          return;
        }
        fQPDG=spQPDG;           // Baryon is SigmaP
        fMass=mSigP;
        sQPDG=gQPDG;            // Meson is gamma
        sMass=0.;
      }
      else if(qM<mSpPi0)        // Only Lambda+PiP is possible
      {
        fQPDG=lQPDG;            // Baryon is Lambda
        fMass=mLamb;
        sQPDG=pipQPDG;          // Meson is Pi+
        sMass=mPi;
      }
      else if(qM<mSzPi)         // Both Lambda+PiP & Sigma0+Pi0 are possible
      {
        if(G4UniformRand()<.6)  // @@ Relative probability (take into account Phase Space)
        {
          fQPDG=lQPDG;          // Baryon is Lambda
          fMass=mLamb;
          sQPDG=pipQPDG;        // Meson is Pi+
          sMass=mPi;
        }
        else
        {
          fQPDG=spQPDG;         // Baryon is SigmaP
          fMass=mSigP;
        }
      }
      else if(qM<mPKZ)          // Lambda+PiP, SigmaP+Pi0, & Sigma0+PiP are possible
      {
        G4double ra=G4UniformRand();
        if(ra<.4)               // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=lQPDG;          // Baryon is Lambda
          fMass=mLamb;
          sQPDG=pipQPDG;        // Meson is Pi+
          sMass=mPi;
        }
        else if(ra<.7)          // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=spQPDG;         // Baryon is SigmaP
          fMass=mSigP;
        }
        else
        {
          fQPDG=szQPDG;         // Baryon is SigmaZ
          fMass=mSigZ;
          sQPDG=pipQPDG;        // Meson is Pi+
          sMass=mPi;
        }
      }
      else                      // Lambda+PiP, SigmaP+Pi0, Sigma0+PiP, p+K0 are possible
      {
        G4double ra=G4UniformRand();
        if(ra<.35)              // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=lQPDG;          // Baryon is Lambda
          fMass=mLamb;
          sQPDG=pipQPDG;        // Meson is Pi+
          sMass=mPi;
        }
        else if(ra<.6)          // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=spQPDG;         // Baryon is SigmaP
          fMass=mSigP;
        }
        else if(ra<.8)          // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=szQPDG;         // Baryon is SigmaZ
          fMass=mSigZ;
          sQPDG=pipQPDG;        // Meson is Pi+
          sMass=mPi;
        }
        else
        {
          fQPDG=pQPDG;          // Baryon is Proton
          fMass=mProt;
          sQPDG=kzQPDG;         // Meson is K0
          sMass=mK0;
        }
      }
    }
    else if(theLC==-1)          // SigmaM: SigM+gam,SigM+Pi0,S0+PiM,L+Pi-,n+KM are possible
    {
      if(qM<mLPi)               // Only SigmaMinus + gamma decay is possible
      {
        if(qM < mSigM+eps)
        {
#ifdef debug
          G4cout<<"G4QEnv::DecayBaryon: Fill SigmaMinus AsIs"<<G4endl;
#endif
          qH->SetQPDG(smQPDG);
          HV->push_back(qH); // Fill AsIs (delete equivalent)
          if(qM+eps<mSigM)G4cout<<"-Warning-G4QEnv::DecayBaryon: Si- AsIs, M="<<qM<<G4endl;
          return;
        }
        fQPDG=smQPDG;           // Baryon is SigmaP
        fMass=mSigM;
        sQPDG=gQPDG;            // Meson is gamma
        sMass=0.;
      }
      else if(qM<mSzPi)         // Only Lambda+PiM is possible
      {
        fQPDG=lQPDG;            // Baryon is Lambda
        fMass=mLamb;
        sQPDG=pimQPDG;          // Meson is Pi-
        sMass=mPi;
      }
      else if(qM<mSzPi)         // Both Lambda+PiM & Sigma0+PiM are possible
      {
        if(G4UniformRand()<.6)  // @@ Relative probability (take into account Phase Space)
        {
          fQPDG=lQPDG;          // Baryon is Lambda
          fMass=mLamb;
          sQPDG=pimQPDG;        // Meson is Pi-
          sMass=mPi;
        }
        else
        {
          fQPDG=szQPDG;         // Baryon is Sigma0
          fMass=mSigZ;
          sQPDG=pimQPDG;        // Meson is Pi-
          sMass=mPi;
        }
      }
      else if(qM<mPKZ)          // Lambda+PiM, Sigma0+PiM, & SigmaM+Pi0 are possible
      {
        G4double ra=G4UniformRand();
        if(ra<.4)               // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=lQPDG;          // Baryon is Lambda
          fMass=mLamb;
          sQPDG=pimQPDG;        // Meson is Pi-
          sMass=mPi;
        }
        else if(ra<.7)          // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=szQPDG;         // Baryon is Sigma0
          fMass=mSigZ;
          sQPDG=pimQPDG;        // Meson is Pi-
          sMass=mPi;
        }
        else
        {
          fQPDG=smQPDG;         // Baryon is SigmaM
          fMass=mSigM;
        }
      }
      else                      // Lambda+PiM, Sig0+PiM, SigM+Pi0, n+KM are possible
      {
        G4double ra=G4UniformRand();
        if(ra<.35)              // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=lQPDG;          // Baryon is Lambda
          fMass=mLamb;
          sQPDG=pimQPDG;        // Meson is Pi-
          sMass=mPi;
        }
        else if(ra<.6)          // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=szQPDG;         // Baryon is Sigma0
          fMass=mSigZ;
          sQPDG=pimQPDG;        // Meson is Pi-
          sMass=mPi;
        }
        else if(ra<.8)          // @@ Rel. probab. (take into account Phase Space)
        {
          fQPDG=smQPDG;         // Baryon is SigmaM
          fMass=mSigM;
        }
        else
        {
          fQPDG=nQPDG;          // Baryon is Proton
          fMass=mNeut;
          sQPDG=kmQPDG;         // Meson is K-
          sMass=mK;
        }
      }
    }
    else if(theLC==2 || theLC==-2) // SigmaPlus+PiPlus or SigmaMinus+PiMinus
    {
      if(theLC==2 && qM>=mSpPi) // SigmaPlus+PiPlus decay is possible
      {
        fQPDG=spQPDG;           // Baryon is SigmaP
        fMass=mSigP;
        sQPDG=pipQPDG;          // Pi+ Meson
        sMass=mPi;
      }
      if(theLC==-2 && qM>=mSmPi)// SigmaPlus+PiPlus decay is possible
      {
        fQPDG=smQPDG;           // Baryon is SigmaP
        fMass=mSigM;
        sQPDG=pimQPDG;          // Pi- Meson
        sMass=mPi;
      }
      else
      {
#ifdef debug
        G4cout<<"-Warning-G4QE::DecBary:*AsIs* Baryon(S=1,|C|=2),QC="<<qH->GetQC()<<G4endl;
#endif
        HV->push_back(qH);                 // Fill AsIs (delete equivalent)
        return;
      }
    }
    else 
    {
      //KsiM: KsiM+Pi0=1456.29, Ksi0+Pi=1454.4, L+K=1609.36, Sig0+K=1686.32, SigM+K0=1695.1
      //KsiZ: Ksi0+Pi0=1449.81, KsiM+Pi=1460.9, L+K0=1613.3, Sig0+K0=1690.3, SigP+K=1683.05
      //Omeg: Omeg+Pi0=1807.43, Ksi0+K=1808.5, KsiM+K0=1818.96
      G4cout<<"-Warning-G4QE::DecBary:*AsIs* UnknownBaryon(S=1)QC="<<qH->GetQC()<<G4endl;
      HV->push_back(qH);                 // Fill AsIs (delete equivalent)
      return;
    }
  }
  else 
  {
#ifdef debug
    G4cout<<"---Warning---G4QE::DecBary:*AsIso*UnknBaryon(AntiS),QC="<<qH->GetQC()<<G4endl;
#endif
    HV->push_back(qH);                 // Fill AsIs (delete equivalent)
    return;
  }
  G4LorentzVector f4Mom(0.,0.,0.,fMass);
  G4LorentzVector s4Mom(0.,0.,0.,sMass);
  G4double sum=fMass+sMass;
  if(fabs(qM-sum)<eps)
  {
    f4Mom=q4M*(fMass/sum);
    s4Mom=q4M*(sMass/sum);
  }
  else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
  {
#ifdef debug
    G4cout<<"---Warning---G4QE::DecBar:fPDG="<<fQPDG<<"(M="<<fMass<<")+sPDG="<<sQPDG<<"(M="
          <<sMass<<") > TotM="<<q4M.m()<<G4endl;
#endif
    if(!theEnvironment.GetA())
    {
      G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());
      if(!CheckGroundState(quasH,true)) HV->push_back(qH); // Cor or fill asItIs
      else delete qH;  
      delete quasH;
      return;
    }
    else
    {
      delete qH;
      G4ExceptionDescription ed;
      ed << "Baryon DecayIn2 error: Can't Correct, *EmptyEnv*="
         << theEnvironment << G4endl;
      G4Exception("G4QEnvironment::DecayBaryon()", "HAD_CHPS_0001",
                  FatalException, ed);
    }
  }
#ifdef debug
  G4cout<<"G4QEnv::DecayBaryon: *DONE* f4M="<<f4Mom<<",fPDG="<<fQPDG<<", s4M="<<s4Mom
        <<",sPDG="<<sQPDG<<G4endl;
#endif
  delete qH;
  //
  G4QHadron* H1 = new G4QHadron(fQPDG,f4Mom); // Create a Hadron for the 1-st baryon
  HV->push_back(H1);                  // Fill "H1" (delete equivalent)
  G4QHadron* H2 = new G4QHadron(sQPDG,s4Mom); // Create a Hadron for the 2-nd baryon
  HV->push_back(H2);                  // Fill "H2" (delete equivalent)
} // End of DecayBaryon

//Decay of the excited Meson in meson & meson (gamma)
void G4QEnvironment::DecayMeson(G4QHadron* qH, G4QHadronVector* HV)
{
  static const G4QPDGCode gQPDG(22);
  static const G4QPDGCode pizQPDG(111);
  static const G4QPDGCode pipQPDG(211);
  static const G4QPDGCode pimQPDG(-211);
  static const G4QPDGCode kmQPDG(-321);
  static const G4QPDGCode kzQPDG(-311);
  static const G4QPDGCode kpQPDG(321);
  static const G4QPDGCode akzQPDG(311);
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mPi0 = G4QPDGCode(111).GetMass();
  static const G4double mK   = G4QPDGCode(321).GetMass();
  static const G4double mK0  = G4QPDGCode(311).GetMass();
  static const G4double m2Pi0 = mPi0+ mPi0;
  static const G4double mPiPi0= mPi + mPi0;
  static const G4double mPiPi = mPi + mPi;
  static const G4double mKPi0 = mPi0+ mK;
  static const G4double mK0Pi0= mPi0+ mK0;
  static const G4double mKPi  = mPi + mK;
  static const G4double mK0Pi = mPi + mK0;
  static const G4double mK0K  = mK0 + mK;
  static const G4double mK0K0 = mK0 + mK0;
  static const G4double mKK   = mK  + mK;
  static const G4double eps  = 0.003;
  //static const G4QNucleus vacuum(90000000);
  G4int theLB= qH->GetBaryonNumber();          // Baryon number of the Meson
  if(theLB)
  {
    G4cerr<<"*Warning*G4QEnvironment::DecayMeson:A="<<theLB<<" != 0 -> Fill asIs"<<G4endl;
#ifdef ppdebug
    // throw G4QException("G4QEnv::DecayMeson: Unknown Meson with A!=0");
    G4Exception("G4QEnvironment::DecayMeson()", "HAD_CHPS_0000",
                FatalException, "Unknown Meson with A!=0");
#endif
    HV->push_back(qH);                 // Fill AsIs (delete equivalent)
    return;
  }
  G4LorentzVector q4M = qH->Get4Momentum();    // Get 4-momentum of the Meson
  G4double         qM = q4M.m();               // Excited Mass of the Meson
  G4int theLC= qH->GetCharge();                // Chsrge of the Meson
  G4int theLS= qH->GetStrangeness();           // Strangness of the Meson
  if(qM < eps)
  {
    HV->push_back(qH);    // Fill AsIs (delete equivalent)
    if(theLC || theLS) G4cout<<"-Warning-G4QEnv::DecMes: S="<<theLS<<", C="<<theLC<<G4endl;
    return;
  }
#ifdef debug
  G4cout<<"G4QEnv::DecayMeson: *Called* S="<<theLS<<",C="<<theLC<<",4M="<<q4M<<qM<<G4endl;
#endif
  // Select a chanel of the Excited Meson decay
  G4QPDGCode     fQPDG = pipQPDG;              // Prototype for Meson1 = Pi+
  G4double       fMass= mPi;
  G4QPDGCode     sQPDG = pizQPDG;              // Prototype for Meson2 = Pi0
  G4double       sMass= mPi0;
  if(!theLS)                    // This is only for not strange Excited Mesons
  {
    if(!theLC)                  // Rho0 like
    {
      if(qM < m2Pi0)            // Only Pi0+gamma decay is possible
      {
        if(qM < mPi0+eps)
        {
#ifdef debug
          G4cout<<"G4QEnv::DecayMeson: Fill Pi0 AsIs"<<G4endl;
#endif
          qH->SetQPDG(pizQPDG);
          HV->push_back(qH);    // Fill AsIs (delete equivalent)
          if(qM+eps < mPi0)G4cout<<"-Warning-G4QEnv::DecayMeson: Pi0 AsIs, M="<<qM<<G4endl;
          return;
        }
        fQPDG=gQPDG;            // Meson1 is gamma
        fMass=0.;
#ifdef debug
        G4cout<<"-Warning-G4QEnv::DecayMeson: Gamma+Pi0, M="<<qM<<",d="<<qM-mPi0<<G4endl;
#endif
      }
      else if(qM < mPiPi)       // only Pi0+Pi0 is possible
      {
        fQPDG=pizQPDG;          // Meson1 is Pi0
        fMass=mPi0;
      }
      else
      {
        if(G4UniformRand() < .333333333) // @@ can be smaller or dec \w qM !
        {
          fQPDG=pizQPDG;        // Meson1 is Pi0
          fMass=mPi0;
        }
        else
        {
          sQPDG=pimQPDG;        // Meson2 is Pi- (Pi+/Pi- decay)
          sMass=mPi;
        }
      }
    }
    else if(theLC==1)           // Rho+ like decay
    {
      if(qM < mPiPi0)           // Only gamma+Pi+ decay is possible
      {
        if(qM < mPi+eps)
        {
#ifdef debug
          G4cout<<"G4QEnv::DecayMeson: Fill Pi+ AsIs"<<G4endl;
#endif
          qH->SetQPDG(pipQPDG);
          HV->push_back(qH);    // Fill AsIs (delete equivalent)
          if(qM+eps < mPi) G4cout<<"-Warning-G4QEnv::DecayMeson: Pi+ AsIs, M="<<qM<<G4endl;
          return;
        }
        sQPDG=gQPDG;            // Meson is gamma (gamma + Pi+ decay)
        sMass=0.;
#ifdef debug
        G4cout<<"-Warning-G4QEnv::DecayMeson: Gamma+Pi+, M="<<qM<<",d="<<qM-mPi0<<G4endl;
#endif
      }
      //else // Pi0+Pi+ is a default case
    }
    else if(theLC==2)           // Meson++ like: only PiP+PiP is possible
    {
      if(qM > mPiPi)            // Only Pi+ + Pi+ decay is possible
      {
        sQPDG=pipQPDG;          // Meson2 is Pi+
        sMass=mPi;
      }
      else                      // @@ Can be aReason to search for anError in Fragmentation
      {
#ifdef debug
        G4cout<<"-Warning-G4QE::DecayMeson:*AsIs* Meson++ M="<<qM<<",d="<<qM-mPiPi<<G4endl;
#endif
        HV->push_back(qH);      // Fill AsIs (delete equivalent)
        return;
      }
    }
    else if(theLC==-1)          // Rho- like decay
    {
      if(qM < mPiPi0)           // Only gamma + Pi- decay is possible
      {
        if(qM < mPi+eps)
        {
#ifdef debug
          G4cout<<"G4QEnv::DecayMeson: Fill Pi- AsIs"<<G4endl;
#endif
          qH->SetQPDG(pimQPDG);
          HV->push_back(qH);    // Fill AsIs (delete equivalent)
          if(qM+eps < mPi) G4cout<<"-Warning-G4QEnv::DecayMeson: Pi- AsIs, M="<<qM<<G4endl;
          return;
        }
        fQPDG=pimQPDG;          // Meson1 is Pi-
        sQPDG=gQPDG;            // Meson2 is gamma (gamma + Pi- decay)
        sMass=0.;
#ifdef debug
        G4cout<<"-Warning-G4QEnv::DecayMeson: Gamma+Pi-, M="<<qM<<",d="<<qM-mPi<<G4endl;
#endif
      }
      else fQPDG=pimQPDG;       // Meson1 is Pi- instead of Pi+
    }
    else if(theLC==-2)          // Meson-- like: only p+PiP is possible
    {
      if(qM > mPiPi)            // Only Pi- + Pi- decay is possible
      {
        fQPDG=pimQPDG;          // Meson1 is Pi-
        sQPDG=pimQPDG;          // Meson2 is Pi-
        sMass=mPi;
      }
      else                      // @@ Can be aReason to search for anError in Fragmentation
      {
#ifdef debug
        G4cout<<"-Warning-G4QE::DecayMeson:*AsIs* Meson-- M="<<qM<<",d="<<qM-mPiPi<<G4endl;
#endif
        HV->push_back(qH);      // Fill AsIs (delete equivalent)
        return;
      }
    }
    else                        // @@ Growing point: No solution for |C| > 2 so far
    {
#ifdef debug
      G4cout<<"-Warning-G4QE::DecayMeson:*AsIs* UnknMeson (S=0) QC="<<qH->GetQC()<<G4endl;
#endif 
      HV->push_back(qH);        // Fill AsIs (delete equivalent)
      return;
    }
  }
  else if( theLS == 1 )         // ==>->->-> Strange mesons (K-* like) S = 1 <-<-<-<==
  {
    if(!theLC)                  // -->> K0* like decay
    {
      if(qM < mKPi)             // Only K0+gamma decay is possible KPi < K0Pi0
      {
        if(qM < mK0+eps)        // Can not decay, hopefully it is close to K0
        {
#ifdef debug
          G4cout << "G4QEnv::DecayMeson: Fill K0 AsIs" << G4endl;
#endif
          qH->SetQPDG(kzQPDG);
          HV->push_back(qH);    // Fill AsIs (delete equivalent)
          if(qM+eps < mK0) G4cout<<"-Warning-G4QEnv::DecayMeson: K0 AsIs, M="<<qM<<G4endl;
          return;
        }
        fQPDG=kzQPDG;           // Meson is K0
        fMass=mK0;
        sQPDG=gQPDG;            // Meson is gamma
        sMass=0.;
      }
      else if(qM < mK0Pi0)      // Only K- + Pi+ is possible
      {
        sQPDG=kmQPDG;           // Meson2 is K-
        sMass=mK; 
      }
      else                      // Both K0 + Pi0 & K- + Pi+ are possible
      {
        if(G4UniformRand()>.5)  // @@ Relative probability (take into account Phase Space)
        {
          sQPDG=kmQPDG;         // Meson2 is K-
          sMass=mK; 
        }
        else
        {
          fQPDG=kzQPDG;         // Meson2 is K0
          fMass=mK0; 
        }
      }
    }
    else if(theLC==-1)          //  -->> K-* like decay
    {
      if(qM < mKPi0)            // Only K- + gamma decay is possible
      {
        if(qM < mK+eps)         // Can not decay, hopefully it is close to K-
        {
#ifdef debug
          G4cout << "G4QEnv::DecayMeson: Fill K- AsIs" << G4endl;
#endif
          qH->SetQPDG(kmQPDG);
          HV->push_back(qH);    // Fill AsIs (delete equivalent)
          if(qM+eps < mK) G4cout<<"-Warning-G4QEnv::DecayMeson: K- AsIs, M="<<qM<<G4endl;
          return;
        }
        fQPDG=kmQPDG;           // Meson is K-
        fMass=mK;
        sQPDG=gQPDG;            // Meson is gamma
        sMass=0.;
      }
      else if(qM < mKPi0)       // Only K- + Pi0 is possible
      {
        fQPDG=kmQPDG;           // Meson1 is K-
        fMass=mK; 
      }
      else                      // Both K0 + Pi0 & K- + Pi+ are possible
      {
        if(G4UniformRand()>.5)  // @@ Relative probability (take into account Phase Space)
        {
          fQPDG=kmQPDG;         // Meson1 is K-
          fMass=mK; 
        }
        else
        {
          fQPDG=pimQPDG;        // Meson1 is Pi-
          sQPDG=kzQPDG;         // Meson2 is K0
          sMass=mK0; 
        }
      }
    }
    else if(theLC== 1)          //  -->> K0 + Pi+  decay only
    {
      if(qM+eps < mK0Pi)        // Nothing can be done for this bad combination (Recover!)
      {
        G4cout<<"-Warniong-G4QEnv::DecayMeson: LowMassPositive Strange Meson AsIs"<<G4endl;
        HV->push_back(qH);      // Fill AsIs (delete equivalent)
        return;
      }
      else                      // Only K- + Pi0 is possible
      {
        sQPDG=kzQPDG;           // Meson2 is K0
        sMass=mK0; 
      }
    }
    else if(theLC== -2)          //  -->> K- + Pi-  decay only
    {
      if(qM+eps < mKPi)         // Nothing can be done for this bad combination (Recover!)
      {
        G4cout<<"-Warniong-G4QEnv::DecayMeson: LowMassDNegativeStrange Meson AsIs"<<G4endl;
        HV->push_back(qH);      // Fill AsIs (delete equivalent)
        return;
      }
      else                      // Only K- + Pi- is possible
      {
        fQPDG=pimQPDG;          // Meson1 is Pi-
        sQPDG=kmQPDG;           // Meson2 is K-
        sMass=mK;
      }
    }
    else 
    {
      G4cout<<"-Warning-G4QE::DecMeson:*AsIs*UnknownMeson.(S=1), QC="<<qH->GetQC()<<G4endl;
      HV->push_back(qH);        // Fill AsIs (delete equivalent)
      return;
    }
  }
  else if(theLS==2)             // ==>->->-> Double srange Meson: S = 2 <-<-<-<==
  {
    if(theLC== 0)               //  -->> K0 + K0  decay only
    {
      if(qM+eps < mK0K0)        // Nothing can be done for this bad combination (Recover!)
      {
        G4cout<<"-Warniong-G4QEnv::DecayMeson: LowMassNeutral DStrange Meson AsIs"<<G4endl;
        HV->push_back(qH);      // Fill AsIs (delete equivalent)
        return;
      }
      else                      // Only K- + Pi0 is possible
      {
        fQPDG=kzQPDG;           // Meson1 is K0
        fMass=mK0;
        sQPDG=kzQPDG;           // Meson2 is K0
        sMass=mK0; 
      }
    }
    else if(theLC== -1)         //  -->> K- + K0  decay only
    {
      if(qM+eps < mK0K)         // Nothing can be done for this bad combination (Recover!)
      {
        G4cout<<"-Warniong-G4QEnv::DecayMeson: LowMassNegativeDStrange Meson AsIs"<<G4endl;
        HV->push_back(qH);      // Fill AsIs (delete equivalent)
        return;
      }
      else                      // Only K- + Pi- is possible
      {
        fQPDG=kmQPDG;           // Meson1 is K-
        fMass=mK;
        sQPDG=kzQPDG;           // Meson2 is K0
        sMass=mK0;
      }
    }
    else if(theLC==-2)          //  -->> K- + K-  decay only
    {
      if(qM+eps < mKK)          // Nothing can be done for this bad combination (Recover!)
      {
        G4cout<<"-Warniong-G4QEnv::DecayMeson:LowMassDNegativeADStrangeMeson AsIs"<<G4endl;
        HV->push_back(qH);      // Fill AsIs (delete equivalent)
        return;
      }
      else                      // Only K+ + K+ is possible
      {
        fQPDG=kmQPDG;           // Meson1 is K-
        fMass=mK;
        sQPDG=kmQPDG;           // Meson2 is K-
        sMass=mK;
      }
    }
    else 
    {
      G4cout<<"-Warning-G4QE::DecMeson:*AsIs*UnknownMeson.(S=2), QC="<<qH->GetQC()<<G4endl;
      HV->push_back(qH);        // Fill AsIs (delete equivalent)
      return;
    }
  }
  else if( theLS ==-1 )         // ==>->->-> AntiStrange mesons (K+* like) S =-1 <-<-<-<==
  {
    if(!theLC)                  // -->> aK0* like decay
    {
      if(qM < mKPi)             // Only aK0+gamma decay is possible
      {
        if(qM < mK0+eps)        // Can not decay, hopefully it is close to K0
        {
#ifdef debug
          G4cout << "G4QEnv::DecayMeson: Fill K0 AsIs" << G4endl;
#endif
          qH->SetQPDG(akzQPDG);
          HV->push_back(qH);    // Fill AsIs (delete equivalent)
          if(qM+eps < mK0) G4cout<<"-Warning-G4QEnv::DecayMeson: aK0 AsIs, M="<<qM<<G4endl;
          return;
        }
        fQPDG=akzQPDG;          // Meson is aK0
        fMass=mK0;
        sQPDG=gQPDG;            // Meson is gamma
        sMass=0.;
      }
      else if(qM < mK0Pi0)      // Only K+ + Pi- is possible
      {
        fQPDG=pimQPDG;          // Meson1 is Pi-
        sQPDG=kpQPDG;           // Meson2 is K+
        sMass=mK; 
      }
      else                      // Both aK0 + Pi0 & K+ + Pi- are possible
      {
        if(G4UniformRand()>.5)  // @@ Relative probability (take into account Phase Space)
        {
          fQPDG=pimQPDG;        // Meson1 is Pi-
          sQPDG=kpQPDG;         // Meson2 is K+
          sMass=mK; 
        }
        else
        {
          fQPDG=akzQPDG;        // Meson1 is aK0
          fMass=mK0; 
        }
      }
    }
    else if(theLC == 1)         //  -->> aK+* like decay
    {
      if(qM < mKPi0)            // Only aK+ + gamma decay is possible
      {
        if(qM < mK+eps)         // Can not decay, hopefully it is close to K+
        {
#ifdef debug
          G4cout << "G4QEnv::DecayMeson: Fill K+ AsIs" << G4endl;
#endif
          HV->push_back(qH);    // Fill AsIs (delete equivalent)
          if(qM+eps < mK) G4cout<<"-Warning-G4QEnv::DecayMeson: K+ AsIs, M="<<qM<<G4endl;
          return;
        }
        fQPDG=kpQPDG;           // Meson is K+
        fMass=mK;
        sQPDG=gQPDG;            // Meson is gamma
        sMass=0.;
      }
      else if(qM < mKPi0)       // Only K+ + Pi0 is possible
      {
        fQPDG=kpQPDG;           // Meson1 is K+
        fMass=mK; 
      }
      else                      // Both K+ + Pi0 & aK0 + Pi+ are possible
      {
        if(G4UniformRand()>.5)  // @@ Relative probability (take into account Phase Space)
        {
          fQPDG=kpQPDG;         // Meson1 is K+
          fMass=mK; 
        }
        else
        {
          sQPDG=akzQPDG;        // Meson2 is aK0
          sMass=mK0; 
        }
      }
    }
    else if(theLC==-1)          //  -->> aK0 + Pi-  decay only
    {
      if(qM+eps < mK0Pi)        // Nothing can be done for this bad combination (Recover!)
      {
        G4cout<<"-Warniong-G4QEnv::DecayMeson: LowMassNegativeAStrange Meson AsIs"<<G4endl;
        HV->push_back(qH);      // Fill AsIs (delete equivalent)
        return;
      }
      else                      // Only aK0 + Pi- is possible
      {
        fQPDG=pimQPDG;          // Meson1 is Pi-
        sQPDG=akzQPDG;          // Meson2 is K0
        sMass=mK0; 
      }
    }
    else if(theLC== 2)          //  -->> K+ + Pi+  decay only
    {
      if(qM+eps < mKPi)         // Nothing can be done for this bad combination (Recover!)
      {
        G4cout<<"-Warniong-G4QEnv::DecayMeson: LowMassDPositiveStrange Meson AsIs"<<G4endl;
        HV->push_back(qH);      // Fill AsIs (delete equivalent)
        return;
      }
      else                      // Only K- + Pi- is possible
      {
        sQPDG=kpQPDG;           // Meson2 is K+
        sMass=mK;
      }
    }
    else 
    {
      G4cout<<"-Warning-G4QE::DecMeson:*AsIs*UnknownMeson.(S=-1),QC="<<qH->GetQC()<<G4endl;
      HV->push_back(qH);        // Fill AsIs (delete equivalent)
      return;
    }
  }
  else if(theLS==-2)            // ==>->->-> Double srange Meson: S =-2 <-<-<-<==
  {
    if(theLC== 0)               //  -->> aK0 + aK0  decay only
    {
      if(qM+eps < mK0K0)        // Nothing can be done for this bad combination (Recover!)
      {
        G4cout<<"-Warniong-G4QEnv::DecayMeson: LowMassNeutralADStrange Meson AsIs"<<G4endl;
        HV->push_back(qH);      // Fill AsIs (delete equivalent)
        return;
      }
      else                      // Only K- + Pi0 is possible
      {
        fQPDG=akzQPDG;          // Meson1 is aK0
        fMass=mK0;
        sQPDG=akzQPDG;          // Meson2 is aK0
        sMass=mK0; 
      }
    }
    else if(theLC== 1)          //  -->> K+ + K0  decay only
    {
      if(qM+eps < mK0K)         // Nothing can be done for this bad combination (Recover!)
      {
        G4cout<<"-Warniong-G4QEnv::DecayMeson: LowMassPositiveADStrangeMeson AsIs"<<G4endl;
        HV->push_back(qH);      // Fill AsIs (delete equivalent)
        return;
      }
      else                      // Only K+ + K0 is possible
      {
        fQPDG=kpQPDG;           // Meson1 is K+
        fMass=mK;
        sQPDG=akzQPDG;          // Meson2 is K0
        sMass=mK0;
      }
    }
    else if(theLC== 2)          //  -->> K+ + K+  decay only
    {
      if(qM+eps < mKK)          // Nothing can be done for this bad combination (Recover!)
      {
        G4cout<<"-Warniong-G4QEnv::DecayMeson:LowMassDPositiveADStrangeMeson AsIs"<<G4endl;
        HV->push_back(qH);      // Fill AsIs (delete equivalent)
        return;
      }
      else                      // Only K+ + K+ is possible
      {
        fQPDG=kpQPDG;           // Meson1 is K+
        fMass=mK;
        sQPDG=kpQPDG;           // Meson2 is K+
        sMass=mK;
      }
    }
    else 
    {
      G4cout<<"-Warning-G4QE::DecMeson:*AsIs*UnknownMeson.(S=-2),QC="<<qH->GetQC()<<G4endl;
      HV->push_back(qH);        // Fill AsIs (delete equivalent)
      return;
    }
  }
  else 
  {
    //#ifdef debug
    G4cout<<"---Warning---G4QE::DecMeso: *Fill AsIso* UnknMeson, QC="<<qH->GetQC()<<G4endl;
    //#endif
    HV->push_back(qH);                 // Fill AsIs (delete equivalent)
    return;
  }
  G4LorentzVector f4Mom(0.,0.,0.,fMass);
  G4LorentzVector s4Mom(0.,0.,0.,sMass);
  G4double sum=fMass+sMass;
  if(fabs(qM-sum)<eps)
  {
    f4Mom=q4M*(fMass/sum);
    s4Mom=q4M*(sMass/sum);
  }
  else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
  {
#ifdef debug
    G4cout<<"---Warning---G4QE::DecMes:fPDG="<<fQPDG<<"(M="<<fMass<<")+sPDG="<<sQPDG<<"(M="
          <<sMass<<") > TotM="<<q4M.m()<<G4endl;
#endif
    if(!theEnvironment.GetA())
    {
      G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());
      if(!CheckGroundState(quasH,true)) HV->push_back(qH); // Cor or fill asItIs
      else delete qH;  
      delete quasH;
      return;
    }
    else
    {
      delete qH;
      G4ExceptionDescription ed;
      ed << "Meson DecayIn2 error: Can't Correct, *EmptyEnv*="
         << theEnvironment << G4endl;
      G4Exception("G4QEnvironment::DecayMeson()", "HAD_CHPS_0001",
                  FatalException, ed);
    }
  }
#ifdef debug
  G4cout<<"G4QEnv::DecayMeson: *DONE* f4M="<<f4Mom<<",fPDG="<<fQPDG<<", s4M="<<s4Mom
        <<",sPDG="<<sQPDG<<G4endl;
#endif
  delete qH;
  //
  G4QHadron* H1 = new G4QHadron(fQPDG,f4Mom); // Create a Hadron for the 1-st baryon
  HV->push_back(H1);                  // Fill "H1" (delete equivalent)
  G4QHadron* H2 = new G4QHadron(sQPDG,s4Mom); // Create a Hadron for the 2-nd baryon
  HV->push_back(H2);                  // Fill "H2" (delete equivalent)
} // End of DecayMeson

//Decay of the Antistrange Nucleus in nucleus & K+/antiK0 (at present only S=-1)
void G4QEnvironment::DecayAntistrange(G4QHadron* qH, G4QHadronVector* HV)
{
  static const G4QPDGCode kpQPDG(321);          // QPDG for Antistrange positive kaon
  static const G4QPDGCode kzQPDG(311);          // QPDG for Antistrange neutral anti-kaon
  static const G4double mK =G4QPDGCode(321).GetMass(); // Mass of antistrange positive Kaon
  static const G4double mK0=G4QPDGCode(311).GetMass(); // Mass of antistrange neutral aK0
  static const G4double eps=0.003;
  //static const G4QNucleus vacuum(90000000);
  G4int theLS= qH->GetStrangeness();           // Strangness of the Nucleus
  if(theLS>=0)
  {
    G4cerr<<"***G4QEnvironment::DecayAntistrange: S="<<theLS<<", but must be <0"<<G4endl;
    //#ifdef ppdebug
    // throw G4QException("G4QEnv::DecayAntistrange: Not Antistrange nucleus");
    G4Exception("G4QEnvironment::DecayAntistrange()", "HAD_CHPS_0000",
                FatalException, "Not Antistrange nucleus"); 
    //#endif
    //HV->push_back(qH);                 // Fill AsIs (delete equivalent)
    return;
  }
  //else if(theLS<-1)
  //{
  //  G4cout<<"*Warning*G4QEnviron::DecayAntistrange: S="<<theLS<<",AsIs->Improve"<<G4endl;
  //  HV->push_back(qH);                 // Fill AsIs (delete equivalent)
  //  return;
  //}
  G4int astr=-theLS;                           // Number of K+ (or anti-K0)
  G4int theLB= qH->GetBaryonNumber();          // Baryon number of the Nucleus
  G4int theLC= qH->GetCharge();                // Chsrge of the Nucleus
  G4int qPDG = qH->GetPDGCode();               // PDG Code of the decaying Nucleus
  G4int K0PDG= qPDG+astr*999999;               // Residual nonStrange nucleus for S*antiK0
  G4QPDGCode K0QPDG(K0PDG);                    // QPDG of the nuclear residual for S*antiK0
  G4double rK0M=K0QPDG.GetMass();              // Mass of the nuclear residual for S*antiK0
  G4int KpPDG= qPDG+astr*999000;               // Residual nonStrange nucleus for S*K+
  G4QPDGCode KpQPDG(KpPDG);                    // QPDG of the nuclear residual for S*K+
  G4double rKpM=KpQPDG.GetMass();              // Mass of the nuclear residual for S*K+
  G4LorentzVector q4M = qH->Get4Momentum();    // Get 4-momentum of the Nucleus
  G4double         qM = q4M.m();               // Mass of the Nucleus
#ifdef debug
  G4cout<<"G4QE::DecayAntistrang:S="<<theLS<<",C="<<theLC<<",B="<<theLB<<",M="<<qM<<G4endl;
#endif
  // Select a chanel of the decay: @@ The Kaon binding energy is not taken into account !!
  G4QPDGCode     fQPDG = kzQPDG;               // Prototype for Kaon (anti-K0)
  G4double       fMass= mK0;
  G4QPDGCode     sQPDG = K0QPDG;               // Prototype for residual nucleus to Kaon
  G4double       sMass= rK0M;
  if(astr*mK0+rK0M>qM)                              // Can not be K0
  {
    if(astr*mK+rKpM>qM)                             // Can not be K+ too
    {
#ifdef debug
      // @@ Survices, but...
      G4cout<<"*Warning*G4QEnvironment::DecayAntistrange: Too low mass, keep AsIs"<<G4endl;
#endif
      HV->push_back(qH);               // Fill AsIs (delete equivalent)
      return;
    }
    else                                       // Switch to K+
    {
      fQPDG = kpQPDG;                          // Positive Kaon
      fMass= mK;
      sQPDG = KpQPDG;                          // Residual nucleus to K+
      sMass= rKpM;
    }
  }
  else if(astr*mK+rKpM<qM && theLC>theLB-theLC)     // Switch to K+ if Z>N
  {
    fQPDG = kpQPDG;                            // Positive Kaon
    fMass= mK;
    sQPDG = KpQPDG;                            // Residual nucleus to K+
    sMass= rKpM;
  }
  G4double afMass=fMass;
  if(astr>1) afMass*=astr;
  G4LorentzVector f4Mom(0.,0.,0.,afMass);
  G4LorentzVector s4Mom(0.,0.,0.,sMass);
  G4double sum=afMass+sMass;
  if(fabs(qM-sum)<eps)
  {
    f4Mom=q4M*(afMass/sum);
    s4Mom=q4M*(sMass/sum);
  }
  else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
  {
#ifdef debug
    G4cout<<"--Warning--G4QE::DecAntistrange: fPDG="<<fQPDG<<"(M="<<fMass<<")+sPDG="<<sQPDG
          <<"(M="<<sMass<<") > TotM="<<q4M.m()<<G4endl;
#endif
    // G4cerr<<"***G4QEnv::DecayAntistrange: M="<<qM<<", sum="<<sum<<G4endl;
    // throw G4QException("G4QEnv::DecayAntistrange: Nucleus DecayIn2 error");
    G4ExceptionDescription ed;
    ed << "Nucleus DecayIn2 error: DecayAntistrange: M=" << qM << ", sum="
       << sum << G4endl;
    G4Exception("G4QEnvironment::DecayAntistrange()", "HAD_CHPS_0001",
                FatalException, ed);
  }
#ifdef debug
  G4cout<<"G4QEnv::DecayAntistrange: *Done* f4M="<<f4Mom<<",fPDG="<<fQPDG<<", s4M="<<s4Mom
        <<",sPDG="<<sQPDG<<G4endl;
#endif
  delete qH;
  //
  if(astr>1) f4Mom/=astr;
  G4QHadron* H1 = new G4QHadron(fQPDG,f4Mom); // Create a Hadron for the 1-st kaon
  HV->push_back(H1);                  // Fill "H1" (delete equivalent)
  for(G4int ia=1; ia < astr; ++ia)
  {
    H1 = new G4QHadron(fQPDG,f4Mom);          // Create a Hadron for other kaons
    HV->push_back(H1);                // Fill "H1" (delete equivalent)
  }
  G4QHadron* H2 = new G4QHadron(sQPDG,s4Mom); // Create a Hadron for the Residual Nucleus
  HV->push_back(H2);                  // Fill "H2" (delete equivalent)
} // End of DecayAntistrange

// Check that it's possible to decay the TotalResidualNucleus in Quasmon+Environ & correct
// In case of correction the "quasm" is never deleted! If corFlag==true - correction.
G4bool G4QEnvironment::CheckGroundState(G4Quasmon* quasm, G4bool corFlag)
{
  static const G4QPDGCode gamQPDG(22);
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4QContent pimQC(1,0,0,0,1,0);
  static const G4QContent pipQC(0,1,0,1,0,0);
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  //static const G4double dmPi  = mPi+mPi;
  ///@@@///
  ///////////////corFlag=true;
  ///@@@///
  G4QContent valQ=quasm->GetQC();               // Quark content of the Quasmon
  G4int    resQPDG=valQ.GetSPDGCode();          // Reachable in a member function
  G4int    resB=valQ.GetBaryonNumber();         // Baryon number of the Quasmon
  G4int    resC=valQ.GetCharge();               // Charge of the Quasmon
  G4int    resS=valQ.GetStrangeness();          // Strangeness of the Quasmon
  if(resQPDG==10 && resB>0) resQPDG=valQ.GetZNSPDGCode();
  G4double resQMa=G4QPDGCode(resQPDG).GetMass();// GS Mass of the Residual Quasmon
  G4double resEMa=0.;                           // GS Mass of the EmptyResidualEnvironment
  G4bool   bsCond=false;                        // FragSeparatCondition for QuasmonInVacuum
  G4LorentzVector enva4M=G4LorentzVector(0.,0.,0.,0.); // Prototype of the environment Mass
  G4QContent reTQC=valQ;                        // Prototype QuarkContent of the ResidNucl
  G4LorentzVector reTLV=quasm->Get4Momentum();  // Prototyoe 4-Mom of the Residual Nucleus
  G4double reTM=reTLV.m();                      // Real mass of the Quasmon
  G4int envPDG=theEnvironment.GetPDG();

  if ( resB > 1 && 
       ( ( !resS && 
           ( (resC == resB && reTM > resC*mProt) || (!resC && reTM > resB*mNeut) )
         ) ||
         (resS == resB && reTM > resS*mLamb)
       )
     ) // Immediate Split(@@Decay) MultiBaryon
  {
#ifdef chdebug
    G4cout<<"::G4QE::CGS:*MultyB*E="<<envPDG<<",B="<<resB<<",C="<<resC<<",S"<<resS<<G4endl;
#endif
    if(envPDG!=90000000)
    {
      resEMa=theEnvironment.GetMZNS();          // GSMass of the Residual Environment
      enva4M=theEnvironment.Get4Momentum();     // 4-Mom of the Residual Environment
      G4LorentzVector toLV=reTLV+enva4M;        // Total 4-mom for decay
      enva4M=G4LorentzVector(0.,0.,0.,resEMa);  // GSM of the Environment
      reTLV=G4LorentzVector(0.,0.,0.,resQMa);   // GSM of the Quasmon
      if(G4QHadron(toLV).DecayIn2(reTLV,enva4M))
      {
#ifdef debug
        G4cout<<"G4QE::CGS: fill EnvPDG="<<envPDG<<",4M="<<enva4M<<" and continue"<<G4endl;
#endif
        theQHadrons.push_back(new G4QHadron(envPDG,enva4M)); // (delete equivalent)
        theEnvironment=G4QNucleus(G4QContent(0,0,0,0,0,0), G4LorentzVector(0.,0.,0.,0.));
      }
      else G4cout<<"*G4QE::CGS:tM="<<toLV.m()<<toLV<<"<q="<<resQMa<<"+EM="<<resEMa<<G4endl;
    }
    G4int  baPDG=3122;                          // Prototype for MultiLambda
    if(!resS) baPDG = (!resC) ? 2112 : 2212;    // MultiNeutron or MultiProton
#ifdef debug
    G4cout<<"G4QE::CGS: fill "<<resB<<" of "<<baPDG<<" with t4M="<<reTLV<<G4endl;
#endif
    reTLV/=resB;                                // Recalculate! Anyway go out...
    for(G4int ib=0; ib<resB; ib++) theQHadrons.push_back(new G4QHadron(baPDG,reTLV));
    return true;
  }
#ifdef cdebug
  if(resQPDG==89998005)
   G4cout<<"G4QE::CGS:Q="<<valQ<<resQPDG<<",GM="<<resQMa<<",4M="<<reTLV<<reTLV.m()<<G4endl;
#endif
  G4double resSMa=resQMa;                       // Prototype MinimalSplitMass of ResidNucl
  enva4M=theEnvironment.Get4Momentum();         // 4-Mom of the Residual Environment
  if(envPDG!=90000000 && fabs(enva4M.m2())>1.)  // => "Environment is not vacuum" case
  { // @@@@@@@@@@@@@@@@@@@ CALL SUBROUTINE ? @@@@@@@@@
    resEMa=theEnvironment.GetMZNS();            // GSMass of the Residual Environment
#ifdef cdebug
    G4cout<<"G4QE::CGS:EnvironmentExists,gsM="<<resEMa<<",4M="<<enva4M<<enva4M.m()<<G4endl;
#endif
    reTQC+=theEnvironment.GetQCZNS();           // Quark content of the Residual Nucleus
    reTLV+=enva4M;                              // 4-Mom of Residual Nucleus
    //resSMa+=resEMa;                           // Minimal Split Mass of Residual Nucleus
    resSMa=G4QPDGCode(reTQC).GetMass();         // GS Mass of the Residual Quasmon+Environ
  }
  else                                          // Calculate BaryonSeparatCond for vacQuasm
  {
    G4double resQM=reTLV.m();                   // CM Mass of the Residual vacQuasmon
    G4int    baryn=valQ.GetBaryonNumber();      // Baryon Number of the Residual vacQuasmon
    if(baryn>1)                                 // => "Can split baryon?"
    {
      if(valQ.GetN())                           // ===> "Can split neutron?"
      {
        G4QContent resQC=valQ-neutQC;           // QC of Residual for the Neutron
        G4int    resPDG=resQC.GetSPDGCode();    // PDG of Residual for the Neutron
        if(resPDG==10&&resQC.GetBaryonNumber()>0) resPDG=resQC.GetZNSPDGCode();
        G4double resMas=G4QPDGCode(resPDG).GetMass(); // GS Mass of the Residual
        if(resMas+mNeut<resQM) bsCond=true;
      }
      else if(valQ.GetP())                      // ===> "Can split proton?"
      {
        G4QContent resQC=valQ-protQC;           // QC of Residual for the Proton
        G4int    resPDG=resQC.GetSPDGCode();    // PDG of Residual for the Proton
        if(resPDG==10&&resQC.GetBaryonNumber()>0) resPDG=resQC.GetZNSPDGCode();
        G4double resMas=G4QPDGCode(resPDG).GetMass(); // GS Mass of the Residual
        if(resMas+mProt<resQM) bsCond=true;
      }
      else if(valQ.GetL())                      // ===> "Can split Lambda?"
      {
        G4QContent resQC=valQ-lambQC;           // QC of Residual for the Lambda
        G4int    resPDG=resQC.GetSPDGCode();    // PDG of Residual for the Lambda
        if(resPDG==10&&resQC.GetBaryonNumber()>0) resPDG=resQC.GetZNSPDGCode();
        G4double resMas=G4QPDGCode(resPDG).GetMass(); // GS Mass of the Residual
        if(resMas+mLamb<resQM) bsCond=true;
      }
    }
  }
  G4double resTMa=reTLV.m();                    // CM Mass of the ResidualNucleus (Q+Env)
  //if(resTMa>resSMa && (resEMa || bsCond)) return true;// Why not ?? @@ (See G4Q the same)
  G4int nOfOUT = theQHadrons.size();            // Total #of QHadrons at this point
#ifdef cdebug
  G4int    reTPDG=reTQC.GetSPDGCode();
  if(reTPDG==10&&reTQC.GetBaryonNumber()>0) reTPDG=reTQC.GetZNSPDGCode();
  G4cout<<"G4QEnv::CheckGS:(tM="<<resTMa<<"<rQM+rEM="<<resSMa<<",d="<<resSMa-resTMa
        <<" || rEM="<<resEMa<<"=0 & "<<!bsCond<<"=1) & n="<<nOfOUT<<">0 & F="<<corFlag
        <<" then the correction must be done for PDG="<<reTPDG<<G4endl;
#endif
  if ( (resTMa < resSMa || (!resEMa && !bsCond) ) && nOfOUT > 0 && corFlag) // *CORRECTION*
  {
    G4QHadron*  theLast = theQHadrons[nOfOUT-1];
    G4int cNf=theLast->GetNFragments();
    G4int crPDG=theLast->GetPDGCode();
#ifdef cdebug
    G4cout<<"G4QE::CGS: **Correction** lNF="<<cNf<<",lPDG="<<crPDG<<",Q4M="<<reTLV<<G4endl;
#endif
    G4LorentzVector hadr4M = theLast->Get4Momentum();
    G4double  hadrMa=hadr4M.m();
    G4LorentzVector tmpTLV=reTLV+hadr4M;        // Tot (ResidNucl+LastHadron) 4-Mom
#ifdef cdebug
    G4cout<<"G4QE::CGS:YES, s4M/M="<<tmpTLV<<tmpTLV.m()<<">rSM+hM="<<resSMa+hadrMa<<G4endl;
#endif
    G4double tmpTM=tmpTLV.m();
    if(tmpTM>resSMa+hadrMa &&!cNf && crPDG!=22) // Q(E)L contain QM(+EM)+lM ***Last CORR***
    {
      if(resEMa)                                // => "NonVacuumEnvironment exists" case
      {
        G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,resQMa); // GS Mass of Quasmon
        enva4M = G4LorentzVector(0.,0.,0.,resEMa);                 // GS Mass of ResidEnvir
        if(tmpTM>=resQMa+resEMa+hadrMa && G4QHadron(tmpTLV).DecayIn3(hadr4M,quas4M,enva4M))
        {
          //@@CHECK CoulBar (only for ResQuasmon in respect to ResEnv) and may be evaporate
#ifdef debug
          G4cout<<"G4QE::CGS: Modify the Last 4-momentum to "<<hadr4M<<G4endl;
#endif
          theLast->Set4Momentum(hadr4M);
          G4QHadron* quasH = new G4QHadron(valQ, quas4M);
          G4QContent theEQC=theEnvironment.GetQCZNS();
          G4QHadron* envaH = new G4QHadron(theEQC,enva4M);
#ifdef debug
          G4cout<<"G4QE::CGS: Fill Quasm "<<valQ<<quas4M<<" in any form"<<G4endl;
#endif
          EvaporateResidual(quasH, false);   // Try to evaporate residual quasmon (del.eq.)
#ifdef debug
          G4cout<<"G4QE::CGS: Fill envir "<<theEQC<<enva4M<<" in any form"<<G4endl;
#endif
          // @@ Substitute by EvaporateResidual (if it is not used in the evaporateResid)
          envaH->Set4Momentum(enva4M);
          EvaporateResidual(envaH);     // Try to evaporate residual environment(del.eq.)
          // Kill environment as it is already included in the decay
          theEnvironment=G4QNucleus(G4QContent(0,0,0,0,0,0), G4LorentzVector(0.,0.,0.,0.));
        }
        else
        {
#ifdef cdebug
          G4cout<<"***G4QEnv::CheckGroundState: Decay in Frag+ResQ+ResE error"<<G4endl;
#endif
          return false;
        }
      }
      else                                      // => "Environ is vacuum" case (DecayIn2)
      {
        G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,resQMa); // GS Mass of Quasmon
        G4QHadron* quasH = new G4QHadron(valQ, quas4M);
        if(tmpTM>=resQMa+hadrMa && G4QHadron(tmpTLV).DecayIn2(hadr4M,quas4M))//DeIn2(noEnv)
        {
          //@@CHECK CoulBar (only for ResQuasmon in respect to ResEnv) & may be evaporate
          theLast->Set4Momentum(hadr4M);
#ifdef debug
          G4cout<<"G4QE::CGS: modify theLast 4M="<<hadr4M<<hadr4M.m()<<G4endl;
#endif
          quasH->Set4Momentum(quas4M);
#ifdef debug
          G4cout<<"G4QE::CGS: fill newQH "<<valQ<<quas4M<<quas4M.m()<<" inAnyForm"<<G4endl;
#endif
          EvaporateResidual(quasH, false);   // Try to evaporate residual quasmon (del.eq.)
        }
        else
        {
          delete quasH;                         // Delete "Quasmon Hadron"
#ifdef cdebug
          G4cout<<"***G4QEnv::CheckGS: Decay in Fragm+ResQ did not succeeded"<<G4endl;
#endif
          return false;
        }
      }
    }
    else                                        // *** Try Last+Previous CORRECTION ***
    {
#ifdef cdebug
      G4cout<<"G4QEnv::CheckGS: the Last did not help, nH="<<nOfOUT<<G4endl;
#endif
      if(nOfOUT>1)                              // Cor with Last&Previous can be tryed
      {
        G4QHadron*  thePrev = theQHadrons[nOfOUT-2]; // Get pointer to Prev before Last
        if(thePrev->GetNFragments()||thePrev->GetNFragments()) return false;// Dec H/g
        G4LorentzVector prev4M = thePrev->Get4Momentum();
        G4double  prevMa=prev4M.m();            // Mass of previous hadron
        tmpTLV+=prev4M;                         // Increase Total 4-Mom of TotalResidNucl
        G4int      totPDG=reTQC.GetSPDGCode();  // PDG Code of Total Residual Nucleus 
        if(totPDG==10&&reTQC.GetBaryonNumber()>0) totPDG=reTQC.GetZNSPDGCode();
        G4double   tQMa=G4QPDGCode(totPDG).GetMass(); // GS Mass of the Residual Nucleus
#ifdef cdebug
        G4cout<<"G4QE::CGS:T3="<<tmpTLV<<tmpTLV.m()<<">t+h+p="<<tQMa+hadrMa+prevMa<<G4endl;
#endif
        if(tmpTLV.m()>tQMa+hadrMa+prevMa)
        {
          G4LorentzVector nuc4M = G4LorentzVector(0.,0.,0.,tQMa);// 4-Mom of ResidNucAtRest
          G4QHadron* nucH = new G4QHadron(reTQC, nuc4M);
          if(!G4QHadron(tmpTLV).DecayIn3(hadr4M,prev4M,nuc4M))
          {
            delete nucH;                        // Delete "Residual Nucleus Hadron"
            G4cout<<"---Warning---G4QE::CGS:DecayIn ResNuc+LastH+PrevH Error"<<G4endl;
            return false;
          }
          else
          {
            theLast->Set4Momentum(hadr4M);
            thePrev->Set4Momentum(prev4M);
            nucH->Set4Momentum(nuc4M);
#ifdef cdebug
            G4cout<<"G4QE::CGS:*SUCCESS**>CHECK, D4M="<<tmpTLV-hadr4M-prev4M-nuc4M<<G4endl;
#endif
#ifdef debug
            G4cout<<"G4QE::CGS: Fill nucleus "<<reTQC<<nuc4M<<" in any form"<<G4endl;
#endif
            EvaporateResidual(nucH, false);   // Try to evaporate residual nucleus(del.eq.)
          }
        }
        else                                  // ===> Try to use any hadron from the output
        {
#ifdef cdebug
          G4cout<<"G4QE::CGS:Prev&Last didn't help,nH="<<nOfOUT<<">2, MQ="<<resQMa<<G4endl;
#endif
          G4int nphot=-1;                     // #of photons
          G4int npip=-1;                      // #of Pi+
          G4int npiz=-1;                      // #of Pi0
          G4int npim=-1;                      // #of Pi-
          G4double emaz=0.;                   // The maximum energy for pi0 (to sel high E)
          for(G4int id=nOfOUT-1; id>=0; id--) // Search for photons and pi+, and pi-
          {
            G4QHadron* curHadr = theQHadrons[id];
            G4int hPDG=curHadr->GetPDGCode();
            if(hPDG == 22) nphot=id;              // Get the earliest gamma
            else if(hPDG==211 && npip<1) npip=id; // Get the latest Pi+
            else if(hPDG==111)
            {
              G4double piE=curHadr->Get4Momentum().e();
#ifdef chdebug
              G4cout<<"::G4QE::CheckGroundState:"<<id<<",Epi0="<<piE<<">"<<emaz<<G4endl;
#endif
              if(piE>emaz)
              {
                npiz=id;
                emaz=piE;
              }
            }
            else if(hPDG==-211 && npim<1) npim=id; // Get the latest Pi-
          }
          if(nphot>=0)                 // Photon is found, try to use it to resolve PANIC
          {
            G4QHadron* curHadr = theQHadrons[nphot];      // Pointer to the photon
            G4LorentzVector ch4M=curHadr->Get4Momentum(); // 4-Mom of the Photon
            G4LorentzVector tt4M=ch4M+reTLV;// (resQMa=GSMass of the ResidQuasmon(+Env.))
            G4double ttM=tt4M.m();          // Mass of the Phot+ResidQm compaund system
            if(resQMa<ttM)                  // PANIC can be resolved with this Photon
            {
              G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,resSMa);//GSMass of Quasm
              G4QHadron* rqH = new G4QHadron(reTQC,quas4M); // Prototype of outputResidQ
              if(!G4QHadron(tt4M).DecayIn2(ch4M,quas4M))
              {
                delete rqH;                 // Delete tmp "Residual Quasmon Hadron"
#ifdef cdebug
                G4cout<<"***G4QEnv::CheckGS: Decay in Photon+ResQ tM="<<ttM<<G4endl;
#endif
              }
              else
              {
                curHadr->Set4Momentum(ch4M);// Change 4M of Photon (reduced by decay)
                rqH->Set4Momentum(quas4M);  // Fill 4M of the GS Residual Quasmon
#ifdef cdebug
                G4cout<<"G4QE::CGS:nPhot="<<nphot<<",ph4M="<<ch4M<<"+r4M="<<quas4M<<G4endl;
#endif
#ifdef debug
                G4cout<<"G4QE::CGS: Fill Resid "<<reTQC<<quas4M<<" in any form"<<G4endl;
#endif
                EvaporateResidual(rqH, false);//Try to evaporate residual quasmon (del.eq.)
                if(envPDG!=90000000) theEnvironment=G4QNucleus(G4QContent(0,0,0,0,0,0),
                                                          G4LorentzVector(0.,0.,0.,0.));
                return true;
              }
            } // End of the KINEMATIC CHECK FOR THE PHOTON if
          } // End of nphot IF
#ifdef chdebug
          if(resQPDG==89998004) G4cout<<"::G4QE::CGS:S="<<resS<<",B="<<resB<<",C="<<resC
                            <<",+#"<<npip<<",-#"<<npim<<",0#"<<npiz<<",E="<<envPDG<<G4endl;
#endif
          //if(npip>=0&&resQPDG==89998004 || npim>=0&&resQPDG==90003998)// D+D+pi->N+N+pi
          if (envPDG == 90000000 && !resS && resB > 1 && 
              ((npip >= 0 && resC == -2) || (npim >= 0 && resC-resB == 2)) )
   {
            G4int npi=npip;               // (Delta-)+(Delta-)+k*n+(pi+)->(k+2)*n+(pi-)
            G4int piPD=-211;
            G4int nuPD=2112;
            G4double nuM=mNeut;
            if(resC!=-2)                  // (Delta++)+(Delta++)+k*p+(pi-)->(k+2)*p+(pi-)
            {
              npi=npim;
              piPD=211;
              nuPD=2212;
              nuM=mProt;
            }
            G4QPDGCode piQPDG(piPD);
            G4int rB=resB-1;
            G4double suB=rB*nuM;
            G4double suM=suB+nuM+mPi;
            G4QHadron* curHadr = theQHadrons[npi]; // Pointer to the pion of oposit sign
            G4LorentzVector ch4M=curHadr->Get4Momentum(); // 4-Mom of the Pion
            G4LorentzVector tt4M=ch4M+reTLV;// (resQMa=GSMass of the ResidQuasmon(+Env.))
            G4double ttM=tt4M.m();          // Mass of the Pion+ResidQm compaund system
#ifdef cdebug
            if(resQPDG==89998005)
              G4cout<<"G4QE::CGS:Sm="<<suM<<"<Tot="<<ttM<<tt4M
                    <<",pi="<<ch4M<<",Q="<<reTLV.m()<<reTLV<<G4endl;
#endif
            if(suM<ttM)                    // PANIC can be resolved with this Pion
            {
              G4LorentzVector fn4M = G4LorentzVector(0.,0.,0.,suB);//First nucleon(s)
              G4LorentzVector sn4M = G4LorentzVector(0.,0.,0.,nuM);//Second nucleon
              G4LorentzVector pi4M = G4LorentzVector(0.,0.,0.,mPi);//Pion
              if(!G4QHadron(tt4M).DecayIn3(fn4M,sn4M,pi4M))
              {
#ifdef cdebug
                if(resQPDG==89998005)
                  G4cout<<"***G4QEnv::CheckGS:DecayIn3 2N+Pi,tM="<<ttM<<","<<suM<<G4endl;
#endif
              }
              else
              {
                if(rB>1) fn4M/=rB;
                for(G4int ib=0; ib<rB; ib++)
                {
                  G4QHadron* fnH = new G4QHadron(nuPD,fn4M);// First Nucleon(s)
#ifdef debug
                  G4cout<<"G4QE::CGS: fill Nucleon #"<<ib<<", "<<nuPD<<fn4M<<G4endl;
#endif
                  theQHadrons.push_back(fnH); // Fill First Nucleon(s) (del. equivalent)
                }
                G4QHadron* snH = new G4QHadron(nuPD,sn4M);// Second Nucleon
#ifdef debug
                G4cout<<"G4QE::CGS: fill the Last Nucleon, "<<nuPD<<sn4M<<G4endl;
#endif
                theQHadrons.push_back(snH); // Fill Second Nucleon (delete equivalent)
                curHadr->Set4Momentum(pi4M);// Change 4M of the Pion (reduced by decay)
                curHadr->SetQPDG(piQPDG);   // Change Charge of thePion
#ifdef cdebug
                if(resQPDG==89998005)
                  G4cout<<"G4QE::CGS:1="<<nuPD<<fn4M<<rB<<",2="<<sn4M<<",p="<<pi4M<<G4endl;
#endif
                return true;
              }
            } // End of the KINEMATIC CHECK FOR THE PI+/PI- if
          } // End of npip/pin Isonucleus IF
          if(envPDG==90000000&&!resS&&resB>1&&npiz>=0&&(resC<-1||resC-resB>1))
          {
#ifdef chdebug
            G4cout<<"*::*G4QE::CGS:Pi0,rPDG="<<resQPDG<<",rC="<<resC<<",rB="<<resB<<G4endl;
#endif
            G4int npi=-resC;                // k*(Delta-)+m*n+pi0->(k+m)*k+(pi-)
            G4int piPD=-211;
            G4int nuPD=2112;
            G4double nuM=mNeut;
            if(resC!=-2)                    // k*(Delta++)+m*p+pi0->(k+m)*p+k*(pi+)
            {
              npi=resC-resB;
              piPD=211;
              nuPD=2212;
              nuM=mProt;
            }
            G4QPDGCode piQPDG(piPD);
            G4double suB=resB*nuM;          // Total mass of nucleons
            G4double suM=npi*mPi;           // Total mass of pions
            G4double sum=suB+suM;           // Total mass of secondaries
            G4QHadron* curHadr = theQHadrons[npiz]; // Pointer to pi0
            G4LorentzVector ch4M=curHadr->Get4Momentum(); // 4-Mom of the Pion
            G4LorentzVector tt4M=ch4M+reTLV;// (resQMa=GSMass of the ResidQuasmon(+Env.))
            G4double ttM=tt4M.m();          // Mass of the Pion+ResidQm compaund system
#ifdef chdebug
            G4cout<<"::G4QE::CGS:sum="<<sum<<"<"<<ttM<<tt4M<<",pi0="<<ch4M<<",Q="
                  <<reTLV.m()<<reTLV<<G4endl;
#endif
            if(sum<ttM)                     // PANIC can be resolved with this Pi0
            {
              G4LorentzVector fn4M = G4LorentzVector(0.,0.,0.,suB); // Nucleon(s)
              G4LorentzVector pi4M = G4LorentzVector(0.,0.,0.,suM); // Pion(s)
              if(G4QHadron(tt4M).DecayIn2(fn4M,pi4M))
              {
                if(npi>1) pi4M/=npi;
                curHadr->Set4Momentum(pi4M);// Change 4M of the Pion (reduced by decay)
                curHadr->SetQPDG(piQPDG);   // Change Charge of thePion
                if(npi>1) for(G4int ip=1; ip<npi; ip++)
                {
                  G4QHadron* piH = new G4QHadron(piPD,pi4M);// Pion(s)
#ifdef debug
                  G4cout<<"G4QE::CGS: fill Pion #"<<ip<<", "<<piPD<<pi4M<<G4endl;
#endif
                  theQHadrons.push_back(piH); // Fill Pion(s) (delete equivalent)
                }
                if(resB>1) fn4M/=resB;
                for(G4int ib=0; ib<resB; ib++)
                {
                  G4QHadron* fnH = new G4QHadron(nuPD,fn4M);// Nucleon(s)
#ifdef debug
                  G4cout<<"G4QE::CGS: fill IsoNucleon #"<<ib<<", "<<nuPD<<fn4M<<G4endl;
#endif
                  theQHadrons.push_back(fnH); // Fill Nucleon(s) (delete equivalent)
                }
#ifdef chdebug
                G4cout<<"::G4QE::CGS:nucl="<<nuPD<<fn4M<<resB<<",pion="<<pi4M<<npi<<G4endl;
#endif

                return true;
              }
#ifdef chdebug
              else G4cout<<"*::*G4QEnv::CheckGS:DecayIn3:*Pi0* tM="<<ttM<<","<<sum<<G4endl;
#endif
            } // End of the KINEMATIC CHECK FOR THE PI0 if
          } // End of npiz Isonucleus IF
#ifdef cdebug
          G4cout<<"G4QE::CGS: Photons can't help nP="<<nphot<<". TryChangeCharge."<<G4endl;
#endif
          // > Photons did not help, try to find an appropriate partner to join and decay
          G4int    reTBN=reTQC.GetBaryonNumber(); // Baryon number of theHadronicState
          G4int    reTCH=reTQC.GetCharge();       // Charge of theHadronicState
          G4bool isoN = reTCH-reTBN>0 || reTCH<0; // UnavoidableIsonucleus (Delta cond.)
          G4bool norN = reTCH<=reTBN || reTCH>=0; // "Regular nucleus" condition
          G4double nnM=resSMa;               // Fake prototype of the NormalNucleusMass
          G4QContent ipiQC=pipQC;            // Prototype of QCont for the Residual Pion+
          G4QContent nnQC=reTQC+pimQC;       // Prototype of theIsoReduceNucleus(Delta++)
          G4int nnPDG=nnQC.GetSPDGCode();    // Prot. PDGCode of the ResidNormalNucleus
          if((!nnPDG||nnPDG==10)&&nnQC.GetBaryonNumber()>0) nnPDG=nnQC.GetZNSPDGCode();
#ifdef cdebug
          G4cout<<"G4QE::CGS: nnPDR="<<nnPDG<<". TryChangeCharge nOUT="<<nOfOUT<<",Iso="
                <<isoN<<",Nor="<<norN<<",C="<<reTCH<<",B="<<reTBN<<G4endl;
#endif
          if(isoN)                           // Calculations for the Isonuclear Residual
          {
            if(reTCH<0)                      // "at least one Delta-" isostate (chngPort)
            {
              ipiQC=pimQC;                   // Change QCont for the Residual Pion-
              nnQC=reTQC+pipQC;              // Change QCont for theNormalNucleus(Delta-)
              nnPDG=nnQC.GetSPDGCode();      // Change PDGCode of theResidNormalNucleus
              if(nnPDG==10&&nnQC.GetBaryonNumber()>0) nnPDG=nnQC.GetZNSPDGCode();
            }
            G4QPDGCode nnQPDG(nnPDG);        // Now can even have Q-code !
            if(nnPDG<80000000) nnM=nnQPDG.GetMass(); // Mass for the Fundamental Hadron
            else               nnM=nnQPDG.GetNuclMass(nnPDG); // Mass for the Nucleus
          }
          G4bool chx2g=true;
          G4bool force=false; // Force-flag to initiate gamma decays (ChEx=>"Pi0"->2gamma)
          while(chx2g)
          {
            if(force) chx2g=false;
            for(G4int hd=nOfOUT-1; hd>=0; hd--)// Try to use any hadron to resolve PANIC
            {
              G4QHadron* curHadr = theQHadrons[hd];
              G4int chNF=curHadr->GetNFragments();
              G4int chCH=curHadr->GetCharge();
              G4int chBN=curHadr->GetBaryonNumber();
              //G4int chS=curHadr->GetStrangeness();
              G4LorentzVector ch4M=curHadr->Get4Momentum(); // 4Mom of the Current Hadron
#ifdef cdebug
              G4cout<<"G4QE::CGS:#"<<hd<<",ch="<<chCH<<",b="<<chBN<<",4M="<<ch4M<<G4endl;
#endif
              if(!chNF)
              {
                G4LorentzVector tt4M=ch4M+reTLV;// resSMa=GSMass of the ResidQuasmon(+Env)
                G4double chM=ch4M.m();          // Mass of the CurrentHadron from theOUTPUT
                G4double ttM=tt4M.m();          // TotalMass of CurHadr+Residual compaund
                if(isoN)                        // "1 Delta Isonucleus" case
                {
                  if(nnM+mPi+chM<ttM)           // PANIC can be resolved with thisCurHadron
                  {
#ifdef cdebug
                      G4cout<<"G4QE::CGS:CurH+ResQ+Pion t="<<tt4M<<ttM<<",cM="<<chM<<",rM="
                            <<nnM<<", d="<<ttM-chM-nnM-mPi<<G4endl;
#endif
                    ch4M = G4LorentzVector(0.,0.,0.,chM); // Mass of current Hadron
                    G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,nnM); // GSMass of RQ
                    G4LorentzVector ipi4M = G4LorentzVector(0.,0.,0.,mPi);// GSMass of Pion
                    if(G4QHadron(tt4M).DecayIn3(ch4M,ipi4M,quas4M))
                    {
                      curHadr->Set4Momentum(ch4M);// Change 4M of the Current Hadron
                      G4QHadron* rpH = new G4QHadron(ipiQC,ipi4M);// Prototype of ResidPion
#ifdef debug
                      G4cout<<"G4QE::CGS: fill Pion "<<ipiQC<<ipi4M<<G4endl;
#endif
                      theQHadrons.push_back(rpH); // Fill Resid Pion (delete equivalent)
                      G4QHadron* rqH = new G4QHadron(nnQC,quas4M);// Prototype of OutResidQ
#ifdef debug
                      G4cout<<"G4QE::CGS:Fill isoRes "<<nnQC<<quas4M<<" inAnyForm"<<G4endl;
#endif
#ifdef cdebug
                      //if(resQPDG==89998004)
                        G4cout<<"G4QE::CGS:#"<<hd<<"is h="<<curHadr->GetPDGCode()<<ch4M
                              <<curHadr->Get4Momentum()<<" + rq="<<nnPDG<<quas4M<<" + pi="
                              <<ipiQC<<ipi4M<<G4endl;
#endif
                      EvaporateResidual(rqH, false); // Try to evaporate residual (del.eq.)
                      if(envPDG!=90000000)theEnvironment=G4QNucleus(G4QContent(0,0,0,0,0,0)
                                                        ,G4LorentzVector(0.,0.,0.,0.));
                      return true;
                    }
#ifdef cdebug
                    else G4cout<<"***G4QE::CGS:DecIn3 CurH+ResQ+Pion dM="<<ttM-chM<<G4endl;
#endif
                  }
                  if ( (reTCH < 0 && chCH > 0) || (reTCH > reTBN && chCH < chBN) ) //IsoEx?
                  {
                    G4QContent chQC=curHadr->GetQC(); // QuarkCont of the CurrentHadron
                    if(reTCH<0)chQC+=pimQC;           // Add the negativPion QC to CurHadr
                    else       chQC+=pipQC;           // Add the positivePion QC to CurHadr
                    G4QPDGCode nnQPDG=G4QPDGCode(nnPDG);// New QPDG of the Residual
                    nnM=nnQPDG.GetMass();             // New Mass of the Residual
                    G4QPDGCode chQPDG=G4QPDGCode(chQC.GetSPDGCode());// New QPDG of CurHadr
                    chM=chQPDG.GetMass();             // New Mass of the CurHadron
                    if(force && nnPDG==111) nnM=0.;   // Decay of Pi0->2 gammas is possible
#ifdef cdebug
                    G4cout<<"G4QE::CGS:ChargeExchange,cx="<<chx2g<<",rC="<<reTCH<<",rB="
                          <<reTBN<<",rM="<<nnM<<",hC="<<chCH<<",hB="<<chBN<<",hM="<<chM
                          <<",rM+hB="<<nnM+chM<<" < "<<ttM<<G4endl;
#endif
                    if(nnM+chM<ttM)
                    {
                      G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,nnM);//GSMass of RQ
                      G4LorentzVector gam4M = G4LorentzVector(0.,0.,0.,0.);//4Mom of gamma1
                      ch4M = G4LorentzVector(0.,0.,0.,chM);//GSMass of ChrgExchanged Hadron
                      G4QHadron* rqH = new G4QHadron(nnQPDG,quas4M);//ChrgExResidualQuasmon
                      if(!nnM)                      // Decay ResidualVirtualQ: Pi0->2 gamma
                      {
                        if(!G4QHadron(tt4M).DecayIn3(ch4M,quas4M,gam4M))
                        {
                          delete rqH;               // Delete tmp "Residual Quasmon Hadron"
#ifdef cdebug
                          G4cout<<"***G4QE::CGS:DecayIn3 CurH+2Gammas,d="<<ttM-chM<<G4endl;
#endif
                        }
                        else
                        {
                          if(chCH+reTCH-chQC.GetCharge()-nnQC.GetCharge())
                            G4cout<<"**G4QE::CGS:ChEx CH+2G i="<<reTCH<<"+h="<<chCH<<", f="
                                  <<nnQC.GetCharge()<<"+o="<<chQC.GetCharge()<<G4endl;
                          curHadr->Set4Momentum(ch4M);// Change 4M of the Current Hadron
                          curHadr->SetQPDG(chQPDG);   // Change QPDG of the Current Hadron
#ifdef cdebug
                          G4cout<<"G4QE::CGS:SubstituteH#"<<hd<<"->"<<chQPDG<<ch4M<<G4endl;
#endif
                          rqH->Set4Momentum(quas4M);  // Fill 4M of the GS Residual Quasmon
                          rqH->SetQPDG(gamQPDG);      // Change QPDG of the ResidualQuasmon
                          theQHadrons.push_back(rqH); // Fill Gamma 1 as QHadron (del. eq.)
#ifdef debug
                          G4cout<<"G4QE::CGS:Fill (SubRQ) Gamma 1,(22)4M="<<quas4M<<G4endl;
#endif
                          G4QHadron* gamH = new G4QHadron(gamQPDG, gam4M);
                          theQHadrons.push_back(gamH);// Fill Gamma 2 as QHadron (del. eq.)
#ifdef debug
                          G4cout<<"G4QE::CGS:Fill newMadeGamma 2, (22) 4M="<<gam4M<<G4endl;
#endif
                          if(envPDG!=90000000) theEnvironment=
                          G4QNucleus(G4QContent(0,0,0,0,0,0),G4LorentzVector(0.,0.,0.,0.));
                          return true;
                        }
                      } // End of "the Normal decay without 2 gammas"
                      else                        // Normal decay (without "Pi0"->2 gammas)
                      {
                        if(!G4QHadron(tt4M).DecayIn2(ch4M,quas4M))
                        {
                          delete rqH;               // Delete tmp "Residual Quasmon Hadron"
#ifdef cdebug
                          G4cout<<"**G4QE::CGS:DecayIn2 CurH+ResQ d="<<ttM-chM-nnM<<G4endl;
#endif
                        }
                        else
                        {
                          if(chCH+reTCH-chQC.GetCharge()-nnQC.GetCharge())
                            G4cout<<"**G4QE::CGS:ChEx CH+RQ i="<<reTCH<<"+h="<<chCH<<", f="
                                  <<nnQC.GetCharge()<<"+o="<<chQC.GetCharge()<<G4endl;
                          curHadr->Set4Momentum(ch4M);// Change 4M of the Current Hadron
                          curHadr->SetQPDG(chQPDG);   // Change QPDG of the Current Hadron
                          rqH->Set4Momentum(quas4M);  // Fill 4M of the GS Residual Quasmon
#ifdef cdebug
                          G4cout<<"G4QE::CGS:#"<<hd<<",h="<<ch4M<<"+rq="<<quas4M<<G4endl;
#endif
#ifdef debug
                          G4cout<<"G4QE::CGS:FilFr "<<nnQPDG<<quas4M<<" inAnyForm"<<G4endl;
#endif
                          EvaporateResidual(rqH, false); // Try to evaporate resQ (del.eq.)
                          if(envPDG!=90000000) theEnvironment=
                          G4QNucleus(G4QContent(0,0,0,0,0,0),G4LorentzVector(0.,0.,0.,0.));
                          return true;
                        }
                      } // End of "the Normal decay without 2 gammas"
                    }
                    else
                    {
#ifdef cdebug
                      G4cout<<"**G4QE::CGS:rM+hB="<<nnM+chM<<">"<<ttM<<",B="<<chBN<<G4endl;
#endif
                      if(chBN>1)
                      {
                        G4QContent tcQC=chQC+nnQC; //QuarkCont for theTotalCompound nucleus
                        G4QPDGCode tcQPDG(tcQC);         // QPDG for the Total Compound
                        G4double   tcM=tcQPDG.GetMass(); // GS Mass of the TotalCompound
                        G4QHadron* tcH = new G4QHadron(tcQPDG,tt4M);// Hadron=TotalCompound
                        G4int tcS=tcQC.GetStrangeness();            // Total Strangeness
                        G4int tcC=tcQC.GetCharge();                 // Total Charge
                        G4int tcBN=tcQC.GetBaryonNumber();          // Total Baryon Number
                        //Try to decay in DiBaryon or MultyBaryon
                        if ( tcBN == 2 || (!tcS && !tcC) || tcS == tcBN || tcC == tcBN )
                        {
                          if(tcBN==2) theEnvironment.DecayDibaryon(tcH,&theQHadrons); // DB
                          else     theEnvironment.DecayMultyBaryon(tcH,&theQHadrons); // MB
                          G4QHadron* theLast_value = theQHadrons[theQHadrons.size()-1];
                          curHadr->Set4Momentum(theLast_value->Get4Momentum());//4-Mom of CurHadr
                          G4QPDGCode lQP=theLast_value->GetQPDG();
                          if(lQP.GetPDGCode()!=10) curHadr->SetQPDG(lQP);
                          else curHadr->SetQC(theLast_value->GetQC());
                          theQHadrons.pop_back(); // theLastQHadron is excluded from OUTPUT
                          delete theLast_value;// *!!When kill,delete theLastQHadr asAnInstance!*
                        }
                        else if(tcM<ttM)// @@ Can't evaporate here, can only radiate gamma
                        {
                          G4LorentzVector tc4M=G4LorentzVector(0.,0.,0.,tcM);//4M of TotCom
                          G4LorentzVector gc4M=G4LorentzVector(0.,0.,0.,0.); //4M of gamma
                          if(!G4QHadron(tt4M).DecayIn2(tc4M,gc4M))
                          {
                            delete tcH;                 // Delete tmp TotalCompoundHadron
                            curHadr->Set4Momentum(tt4M);// Change 4M of the Current Hadron
                            curHadr->SetQPDG(tcQPDG);   // Change QPDG of theCurrentHadron
#ifdef cdebug
                            G4cout<<"**G4QE::CGS:DecayIn2 TotComp+gam d="<<ttM-tcM<<G4endl;
#endif
                          }
                          else
                          {
                            curHadr->Set4Momentum(gc4M);//Change 4Mom of the Current Hadron
                            curHadr->SetQPDG(gamQPDG);//Change PDG of theCurHadron to gamma
                            tcH->Set4Momentum(tc4M); // Fill 4-Mom of the GS Total Compound
#ifdef cdebug
                            G4cout<<"G4QE::CGS:#"<<hd<<",ch="<<ch4M<<"+t4M="<<tc4M<<G4endl;
#endif
#ifdef debug
                            G4cout<<"G4QE::CGS:FilTC "<<tcQPDG<<tc4M<<" inAnyForm"<<G4endl;
#endif
                            EvaporateResidual(tcH, false);// Try to evaporate hadron (d.e.)
                          }
                        }
                        else // @@ Fill the TotalCompound instead of the CurrentHadron @@
                        {
#ifdef cdebug
                          G4cout<<"G4QE::CGS:**CEF,M="<<tcM<<">"<<ttM<<",B="<<chBN<<G4endl;
#endif
                          delete tcH;                  // Delete tmp "TotalCompound"
                          curHadr->Set4Momentum(tt4M); // Change 4Mom of the Current Hadron
                          curHadr->SetQPDG(tcQPDG);    // Change QPDG of the Current Hadron
                        }
                        return true;
                      } // Check that the residual is a nucleus, not just a nucleon
                    } // Check if pion can still be radiated or must be absorbed
                  } // Check that the iso-exchange could help
                } // End of the IF: KINEMATIC CHECK FOR THE CURRENT HADRON(Isonuclear Case)
                else if(norN)
                {
                  if(resSMa+chM<ttM)          // PANIC can be resolved with this CurHadron
                  {
                    G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,resSMa);// GSMass ofQ
                    G4QHadron* rqH = new G4QHadron(reTQC,quas4M);// Prototype OutputResidQ
                    if(!G4QHadron(tt4M).DecayIn2(ch4M,quas4M))
                    {
                      delete rqH;                   // Delete tmp "Residual Quasmon Hadron"
#ifdef cdebug
                      G4cout<<"***G4QE::CheckGS:Decay in CurH+ResQ dM="<<ttM-chM<<G4endl;
#endif
                    }
                    else
                    {
                      curHadr->Set4Momentum(ch4M);  // Change 4M of the Current Hadron
                      rqH->Set4Momentum(quas4M);    // Fill 4M of the GS Residual Quasmon
#ifdef cdebug
                      G4cout<<"G4QEnv::CheckGS:#"<<hd<<",ch4M="<<curHadr->GetPDGCode()
                            <<ch4M<<" + ResQ4M="<<totPDG<<quas4M<<G4endl;
#endif
#ifdef debug
                      G4cout<<"G4QE::CGS:Fill GSRes "<<reTQC<<quas4M<<" inAnyForm"<<G4endl;
#endif
                      EvaporateResidual(rqH, false); // Try to evaporate residHadron (d.e.)
                      if(envPDG!=90000000)theEnvironment=G4QNucleus(G4QContent(0,0,0,0,0,0)
                                                            ,G4LorentzVector(0.,0.,0.,0.));
                      return true;
                    }
                  } // End of the KINEMATIC CHECK FOR THE CURRENT HADRON if (NormNuclCase)
                } // End of IsoNucleus/NormalNucleus choice
              } // End of the NumberOfFragments=0 (NOT DECAYED PARTICLE) if
            } // End of the LOOP over hadrons and all attempts to resolve PANIC
#ifdef cdebug
            G4cout<<"G4QEnv::CheckGS:***Any hadron from the OUTPUT did not help"<<G4endl;
#endif
            force=true;
          } // End of while for chx2g
          if(resB>1)
          {
            G4int hind=-1;
            if(npiz>-1) hind=npiz;
            else
            {
              if(resC+resC>resB && npim>-1) hind=npim;
              else                          hind=npip;
            }
            if(hind>-1)
            {
              G4QHadron* curHadr = theQHadrons[hind];
              G4int chNF=curHadr->GetNFragments();
              if(!chNF)
              {
                G4LorentzVector ch4M=curHadr->Get4Momentum(); // 4Mom of the Current Hadron
                G4LorentzVector tt4M=ch4M+reTLV; // resSMa=GSMass of the ResidQuasmon(+Env)
                G4double        ttM=tt4M.m();    // TotalMass of CurHadr+Residual compaund
                G4QContent      ttQC=valQ+curHadr->GetQC(); // total Quark Content
                G4int           ttPDG=ttQC.GetZNSPDGCode();
                G4QPDGCode      ttQPDG(ttPDG);
                G4double        ttGSM=ttQPDG.GetMass();
                if(ttM>ttGSM && ttQC.GetStrangeness()>0)//Hypernucleus can be L->N degraded
                {
                  ttPDG-=ttQC.GetStrangeness()*999999; // S Neutrons instead of S Lambdas
                  ttQPDG=G4QPDGCode(ttPDG);            // Update QPDGcode defining fragment
                  ttGSM=ttQPDG.GetMass();              // Update the degraded mass value
#ifdef debug
                  G4cout<<"G4QEnv::CheckGS:Hypernucleus degraded to QPDG="<<ttQPDG<<G4endl;
#endif
                }
                if(ttM>ttGSM)   // Decay of Pi0 in 2 gammas with the residual nucleus
                {
#ifdef cdebug
                  G4cout<<"G4QEnv::CheckGS: Decay in ResQ="<<ttQPDG<<" & 2 gammas"<<G4endl;
#endif
                  G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,ttGSM); // GSMass of RQ
                  G4LorentzVector fgam4M = G4LorentzVector(0.,0.,0.,0.); // 4Mom of gamma 1
                  G4LorentzVector sgam4M = G4LorentzVector(0.,0.,0.,0.); // 4Mom of gamma 2
                  if(!G4QHadron(tt4M).DecayIn3(quas4M,fgam4M,sgam4M))G4cout<<"*3*"<<G4endl;
                  else
                  {
                    curHadr->Set4Momentum(fgam4M); // Change 4M of the Pion to Gamma
                    curHadr->SetQPDG(gamQPDG);      // Change QPDG of the Pion to Gamma
                    G4QHadron* sgH = new G4QHadron(gamQPDG,sgam4M); // Gamma 2
                    theQHadrons.push_back(sgH); // Fill Gamma 2 as QHadron (del. eq.)
                    G4QHadron* rqH = new G4QHadron(ttQPDG,quas4M); // GSResidualQuasmon
                    theQHadrons.push_back(rqH); // Fill GSResidQuasmon as QHadron (del.eq.)
                    return true;
                  }
                }
#ifdef debug
                else G4cout<<"-W-G4QEn::CheckGS:M="<<ttM<<" < GSM="<<ttGSM<<ttQPDG<<G4endl;
#endif
              }
            }
          }
          return false;
        } // End of the KINEMATIC LIMIT FOR THE L&P CORRECTION if/else (try any hadron)
      } // End of the POSSIBILITY OF PREV+LAST (OR MORE) CORRECTION if
      else return false;
    } // End of the CORRECTION WITH THE LAST if/else
  } // End of the CORRECTION IS POSSIBLE if
  else return false;                     // Correction can not be done
  return true;                           // If correction was done successfully
} // End of "CheckGroundState"

// Try to decay the Total Residual Nucleus in Environ+Quasmon
void G4QEnvironment::CopyAndDeleteHadronVector(G4QHadronVector* HV)
{
  G4int nHadrons = HV->size();
  if(nHadrons)
  {
    for(G4int ih=0; ih<nHadrons; ++ih) // LOOP over output QHadrons
    {
      //G4QHadron* inH = HV->operator[](ih);         // Pointer to the i-th QHadron
      G4QHadron* inH = (*HV)[ih];                  // Pointer to the i-th QHadron
      G4int hNF  = inH->GetNFragments();           // A#of secondary fragments
      if(!hNF)                                     // Fill only final hadrons
      {
#ifdef debug
        G4cout<<"G4QEnv::Copy&DeleteHV:#"<<ih<<", hPDG="<<inH->GetPDGCode()<<G4endl;
#endif
        G4QHadron* curH = new G4QHadron(inH);      // will be deleted with allQHadronVector
        theQHadrons.push_back(curH);               // Fill hadron-copy (delete equivalent)
      }
    }
    for_each(HV->begin(), HV->end(), DeleteQHadron()); // Delete instances
    HV->clear();                                       // Delete pointers
  }
#ifdef debug
  else G4cout<<"***G4QEnv::Kopy&DelHV: No hadrons in the QHadronVector"<<G4endl;
#endif
  delete HV;                                       // Delete the inputQHadronVector
} // End of "CopyAndDeleteHadronVector"

// Try to decay the Total Residual Nucleus in Environ+Quasmon
G4bool G4QEnvironment::DecayInEnvQ(G4Quasmon* quasm)
{
  G4QContent valQ=quasm->GetQC();                 // Quark content of the Quasmon
  G4int    resQPDG=valQ.GetSPDGCode();            // Reachable in a member function
  if(resQPDG==10&&valQ.GetBaryonNumber()>0) resQPDG=valQ.GetZNSPDGCode();
  G4double resQMa=G4QPDGCode(resQPDG).GetMass();  // GS Mass of the Residual Quasmon
  G4LorentzVector enva4M=G4LorentzVector(0.,0.,0.,0.);
  G4LorentzVector reTLV=quasm->Get4Momentum();    // Prototyoe of the 4-Mom of the ResidNuc
  G4double resSMa=resQMa;                         // Prototype of MinSplitMass of ResidNucl
  G4int envPDG=theEnvironment.GetPDG();           // PDG Code of the Environment
  if(envPDG!=90000000)                            // => "Environment is not vacuum" case
  {
    G4double resEMa=theEnvironment.GetMZNS();     // GSMass of the Residual Environment
    enva4M=G4LorentzVector(0.,0.,0.,resEMa);      // 4-Mom of the Residual Environment
    reTLV+=enva4M;                                // 4-Mom of Residual Nucleus
    resSMa+=resEMa;                               // Minimal Split Mass of Residual Nucleus
    G4double resTMa=reTLV.m();                    // CM Mass of theResidNucleus (Quasm+Env)
 //#ifdef debug
    G4cout<<"G4QEnv::DecayInEnvQ: totM="<<reTLV<<resTMa<<" > rQM+rEM="<<resSMa<<G4endl;
 //#endif
    if(resTMa>resSMa)
    {
      G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,resQMa); // GS Mass of Quasmon
      G4QHadron* quasH = new G4QHadron(valQ, quas4M);
      G4QHadron* envaH = new G4QHadron(envPDG, enva4M);
      if(!G4QHadron(reTLV).DecayIn2(enva4M,quas4M))
      {
        delete quasH;                             // Delete "Quasmon Hadron"
        delete envaH;                             // Delete "Environment Hadron"
        G4cout<<"---Warning---G4Q::DecInEnvQ:Decay in Environment+ResidualQuasmon"<<G4endl;
        return false;
      }
      else
      {
        quasH->Set4Momentum(quas4M);
        EvaporateResidual(quasH);                // Try to evaporate quasmon (del. equiv.)
        envaH->Set4Momentum(enva4M);
        EvaporateResidual(envaH);                // Try to evaporate residual (del. equiv.)
      }
    }
    else return false;
  }
  else return false;                              // => "Environment is vacuum" case
  return true;
} // End of "DecayInEnvQ"

// Add a Quasmon to the Environment
void G4QEnvironment::AddQuasmon(G4Quasmon* Q)
{
  theQuasmons.push_back(Q);
  totCharge+=Q->GetCharge();
  totBaryoN+=Q->GetBaryonNumber();
  tot4Mom  +=Q->Get4Momentum();
#ifdef debug
  G4cout<<"G4QEnv::AddQuasmon:t4M="<<tot4Mom<<",tC="<<totCharge<<",tB="<<totBaryoN<<G4endl;
#endif
}

// Put all hadrons in the vector to the mass shell
void G4QEnvironment::CheckMassShell(G4QHadronVector* HV)
{
  static const G4double eps=.003;
  G4int nHadrons = HV->size();
  if(nHadrons)
  {
    for(G4int ih=0; ih<nHadrons; ++ih)                  // LOOP over output QHadrons
    {
      G4QHadron* inH = (*HV)[ih];                       // Pointer to the i-th QHadron
      G4int hNF  = inH->GetNFragments();                // A#of secondary fragments
      if(!hNF)                                          // Fill only final hadrons
      {
        G4LorentzVector q4M = inH->Get4Momentum();      // Get 4-momentum of the Hadron
        G4double        qM2 = q4M.m2();                 // Squared Mass of the Hadron
        G4double        qGM = inH->GetQPDG().GetMass(); // Get Ground State Mass
        G4double        qGM2= qGM*qGM;                  // Squared Ground State Mass
#ifdef debug
	  G4cout << "G4QEnv::ChkMassShell:#" << ih <<", hPDG="<< inH->GetPDGCode() <<", dM2="
               << qM2 - qGM*qGM << G4endl;
#endif
        if( fabs(qM2 - qGM*qGM) > eps)                    // Correct the mass shell
        {
          G4double mins= 10000000000.;                    // Minimum excitation found
          G4int    nj  = -1;                               // Minimum j-index
          for(G4int jh=0; jh<nHadrons; ++jh)              // LOOP over output QHadrons
          {
            if(jh != ih)
            {
              G4QHadron* jnH = (*HV)[jh];                 // Pointer to the j-th QHadron
              hNF  = inH->GetNFragments();                // A#of secondary fragments
              if(!hNF)                                    // Fill only final hadrons
              {
                G4LorentzVector j4M = jnH->Get4Momentum();      // Get 4-mom of the Hadron
                G4double        jGM = jnH->GetQPDG().GetMass(); // Get Ground State Mass
                G4double        jGM2= jGM*jGM;                  // Doubled GS Mass
                G4LorentzVector s4M = j4M+q4M;                  // Compound 4-mom
                G4double        jqM = qGM * jGM;                // Worling var.
                G4double        s2M = s4M.m2()-qGM2-jGM2-jqM-jqM;// Closeness parameter
                if(s2M > eps && s2M < mins)
                {
                  mins = s2M;
                  nj   = jh;
                }
		  }
            }
          } // End of the searching LOOP for the rest of hadrons
          //if(nj<0) G4cout<<"-W-G4QE::ChkMShell:NotCorr,M2="<<qM2<<",GSM2="<<qGM2<<G4endl;
          //else                                                // It's possible to correct
          if(nj >= 0)                                           // It's possible to correct
          {
            G4QHadron* jnH = (*HV)[nj];                         // Pointer to j-th QHadron
            G4LorentzVector j4M = jnH->Get4Momentum();          // Get 4-mom of the Hadron
            G4double        jGM = jnH->GetQPDG().GetMass();     // Get Ground State Mass
            G4LorentzVector c4M = q4M + j4M;                    // Get 4-mom of Compound
            G4LorentzVector i4Mom(0.,0.,0.,qGM);
            G4LorentzVector j4Mom(0.,0.,0.,jGM);
            /*
            // DHW 16 June 2011: variable set but not used.  Comment out to fix compiler
            //                   warning.    
            G4bool done = true;
            */
            if(!G4QHadron(c4M).DecayIn2(i4Mom, j4Mom))
            {
              G4cout<<"-Warning-G4QEnv::ChkMShell: tM="<< c4M.m() <<", iM="<< qGM <<", jM="
                    << jGM <<", d="<< c4M.m()-qGM-jGM << G4endl;
              /*
              // DHW 16 June 2011: set but not used.
              done = false;
              */
            }
            else
            {
              inH->Set4Momentum(i4Mom);
              jnH->Set4Momentum(j4Mom);
            }
          }
        }
      }
    }
  }
}
