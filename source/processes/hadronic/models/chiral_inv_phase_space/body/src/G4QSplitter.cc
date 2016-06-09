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
//       1         2         3         4         5         6         7         8         9
//34567890123456789012345678901234567890123456789012345678901234567890123456789012345678901
//
//
// $Id: G4QSplitter.cc,v 1.6 2006/11/27 10:44:55 mkossov Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//      ---------------- G4QSplitter ----------------
//             by Mikhail Kossov, August 2005.
//  class for Hadron-Hadron String Interaction used by the CHIPS Model
// -------------------------------------------------------------------
 
//#define chdebug
//#define debug
//#define sdebug
//#define ppdebug
//#define cdebug
//#define cldebug
//#define edebug
//#define fdebug
//#define pdebug
//#define rdebug
//#define ffdebug
//#define pcdebug
//#define mudebug

#include "G4QSplitter.hh"
#include <cmath>
using namespace std;
 
G4QSplitter::G4QSplitter(G4QHadron projHadron, const G4bool projEnvFlag,
                     const G4int targPDG, const G4bool targEnvFlag) :
  theProjEnvFlag(projEnvFlag), theTargEnvFlag(targEnvFlag), theWeight(1.)
{
  static const G4double mPi0 = G4QPDGCode(111).GetMass();
  //static const G4double Temperature = 180.; // Temperature as a parameter
  //static const G4double Temp2 = Temperature*Temperature; // Squared Temperature
  theWorld= G4QCHIPSWorld::Get();             // Get a pointer to the CHIPS World
#ifdef debug
  G4cout<<"-->G4QString::Constr:pPDG="<<projHadron.GetPDGCode()<<",tPDG="<<targPDG<<G4endl;
#endif
  // Target initialization
  G4QPDGCode targQPDG(targPDG);
  totBaryNum= targQPDG.GetBaryNum();          // Only a Prototype
  totCharge = targQPDG.GetCharge();           // Only a Prototype
  theTargQC = targQPDG.GetQuarkContent();
  G4double tM = targQPDG.GetMass();
  theTarg4Mom = G4LorentzVector(0.,0.,0.,tM);
  G4QHadron targHadron(targQPDG,theTarg4Mom); // Target Hadron

  // Projectile initialization
  theProj4Mom = projHadron.Get4Momentum();
		G4ThreeVector projBoost = theProj4Mom.boostVector(); // Projectile BoostVector to LS
  G4ThreeVector projRBoost= -projBoost;       // Projevtile Boost vector to projectile CMS
  G4QPDGCode projQPDG(projHadron.GetPDGCode());
  theProjQC  = projQPDG.GetQuarkContent();
  totBaryNum+= projQPDG.GetBaryNum();         // Final control value
  totCharge += projQPDG.GetCharge();          // Final control value
  tot4Mom  = theProj4Mom+theTarg4Mom;         // Final control value

  G4double projM = projQPDG.GetMass();
  G4double projM2 = projM*projM;
  G4double pM2 = theProj4Mom.m2();
  G4double tM2 = tM*tM;
  G4double dM2 = fabs(pM2-projM2);
  if(dM2>1.) G4cout<<"-Warn-G4QS::Constr:dM2="<<dM2<<",M2="<<projM2<<",LVM2="<<pM2<<G4endl;
  G4double pM=sqrt(pM2);                     // @@ do we need pM ? @@ (in print)

  // === Print out of the input information at Creation time & tot 4-mom Calculation ======
#ifdef pdebug
  G4cout<<"G4QS::Cons:PQC="<<theProjQC<<",TQC="<<theTargQC<<",P4Mom="<<theProj4Mom<<
				theProj4Mom.m2()<<theProj4Mom.m()<<G4endl;
  G4cout<<"G4QSpl::Constr: tC="<<totCharge<<",tB="<<totBaryNum<<",tot4M="<<tot4Mom<<G4endl;
#endif
  //G4int nP=theWorld->GetQPEntries();       // A#of init'ed particles in CHIPS World (@@?)
  //G4int nCl=nP-90;                          // A#of init'ed clusters in CHIPS World (@@?)
		//#ifdef pdebug
  //G4cout<<"G4QS:Const:Before QEX:n="<<nP<<G4endl;
		//#endif
  // @@@@@@ ===> Here the Quark Exchange Quasmon Creation must be added <=== @@@@@@
  // @@ --- old ---
  //G4int nQTarg=theTargQC.GetTot();           // a#of Quark-partons in the Target
  //G4int nQProj=theProjQC.GetTot();           // a#of Quark-partons in the Progectile
  // @@ --- new ---
  G4LorentzVector tSum=theProj4Mom+theTarg4Mom;
  G4double sumM=tM+pM;
  G4double dn=2.17*(tSum.m2()/sumM/sumM-1.);   // additional partons (gluons) - real
  G4int np=static_cast<int>(dn);               // additional partons (gluons) - basic int
  if(G4UniformRand()>1.-dn+np) np++;           // randomized int additional partons(gluons)
  G4int nQTarg=theTargQC.GetTot()+np;          // a#of Quark-partons in the Target
  G4int nQProj=theProjQC.GetTot()+np;          // a#of Quark-partons in the Progectile
  // @@ --- new ---
  //G4int mQTarg=theTargQC.GetTot();           // a#of Quark-partons in the Target
  //G4int mQProj=theProjQC.GetTot();           // a#of Quark-partons in the Progectile
  ////G4double tothM=tot4Mom.m()/2.;              // C.M. mass of the protons
  ////G4int nQCM=classG4Quasmon.CalculateNumberOfQPartons(tothM);
  //G4double tothM2=tot4Mom.m2()/4.;              // C.M. mass of the protons
  //G4int nQCM=static_cast<int>((1.+sqrt(tothM2/Temp2+1.))/2.);
  //G4int nQTarg=nQCM;           // a#of Quark-partons in the Target
  //if(nQTarg<mQTarg) nQTarg=mQTarg;
  //G4int nQProj=nQCM;           // a#of Quark-partons in the Progectile
  //if(nQProj<mQProj) nQProj=mQProj;
  // @@ --- end ---
  G4double interc=1.;                        // @@ A parameter ?!
#ifdef pdebug
  G4cout<<"G4QSpl::Constr: nP="<<nQProj<<", nT="<<nQTarg<<", intercept="<<interc<<G4endl;
#endif
  // @@ Now projectile can be only meson or baryon @@ -- @@ Improve for clusters @@ --
  G4LorentzVector pq4Mom(0.,0.,0.,0.);       // Prototype of LV of quark of progectile  
  G4double rPMass=0.;                        // Prototype of the residProjMass (Meson case)
  G4bool FreeFraF=false;                     // Prototype of the free exchange Flag
  //if(G4UniformRand()<0.) FreeFraF=true;      // @@@@@@@@ Confirm the free exchange
  if(nQProj<2) G4cout<<"***G4QSpl::Constr: nQProj="<<nQProj<<"<2 ***FatalError***"<<G4endl;
		else if(nQProj>2)                          // ---> Baryon case (clusters are not implem.)
		{
    //if(nQProj>3)G4cout<<"-Wor-G4QS::Const:nQProj="<<nQProj<<">3 is not implem'd"<<G4endl;
    G4double xP=0.;
    if(FreeFraF) xP=RandomizeMomFractionFree(nQProj);
    else
    {
      xP=RandomizeMomFractionString(nQProj);
      theWeight*=pow(xP,interc);
      G4cout<<"************G4QS::C: string xP="<<xP<<G4endl;
    }
    rPMass = sqrt(pM2*(1.-xP));              // Residual Projectile mass
#ifdef pdebug
    G4cout<<"G4QS::C: nQProj="<<nQProj<<", xProj="<<xP<<", rPMass="<<rPMass<<G4endl;
#endif
  }
  G4LorentzVector pr4Mom(0.,0.,0.,rPMass);   // Prototype of LV of the residual projectile
  G4bool projFl=false;
  G4bool targFl=false;
  G4bool tmpBl=projHadron.DecayIn2(pq4Mom,pr4Mom);
  if(tmpBl) projFl=true;                     // Proj decay is OK
  //if(projHadron.DecayIn2(pq4Mom,pr4Mom)) projFl=true;
  else G4cout<<"*Warning*G4QSplitter::Constr:ProjDecIn2 rpM="<<rPMass<<", pM="<<pM<<G4endl;
#ifdef pdebug
  G4cout<<"G4QSpl::Con:"<<projFl<<" split PROJ in R4M="<<pr4Mom<<" & Q4M="<<pq4Mom<<G4endl;
#endif
  G4LorentzVector tq4Mom(0.,0.,0.,0.);       // Prototype of LV of quark of the target  
  //if(nQTarg<3)G4cout<<"***G4QStr::Const: nQTarg="<<nQTarg<<"<3 ***FatalError***"<<G4endl;
  //if(nQTarg>3)G4cout<<"-W-G4QS::Const: nQTarg="<<nQTarg<<">3 is not implemented"<<G4endl;
  G4double xT=0.;
  if(FreeFraF) xT=RandomizeMomFractionFree(nQTarg);
  else
  {
    xT=RandomizeMomFractionString(nQTarg);
    theWeight*=pow(xT,interc);
    G4cout<<"************G4QS::C: string xT="<<xT<<G4endl;
  }
  G4double rTMass = sqrt(tM2*(1.-xT));       // Residual Target mass
#ifdef pdebug
  G4cout<<"G4QS::C: nQTarg="<<nQTarg<<", xProj="<<xT<<", rTMass="<<rTMass<<G4endl;
#endif
  G4LorentzVector tr4Mom(0.,0.,0.,rTMass);   // Prototype of LV of the residual projectile
  if(targHadron.DecayIn2(tq4Mom,tr4Mom)) targFl=true; // Targ decay is OK
  else G4cout<<"**Worning**G4QStr::Constr:TargDecIn2 rtM="<<rTMass<<", tM="<<tM<<G4endl;
#ifdef pdebug
  G4cout<<"G4QStr::Con:"<<targFl<<" split TARG in R4M="<<tr4Mom<<" & Q4M="<<tq4Mom<<G4endl;
#endif
  G4bool elasFl=false;                       // ByDefault avoid the elastic scattering
  if (targFl && projFl)                      // --- @@ Now only for pp case @@ ---
  {
    G4LorentzVector newProj4M=pr4Mom+tq4Mom;
    G4double nPM2=newProj4M.m2();
    G4double npM=sqrt(nPM2);
    G4bool pcorFl=false;                     // ByDefault no MassSh correction for newProj.
    // @@ Just to try ---- Start
    //G4LorentzVector newTarg4M=tr4Mom+pq4Mom;
    //G4double nTM2=newTarg4M.m2();
    //G4double ntM=sqrt(nTM2);
#ifdef pdebug
    //G4cout<<"G4QStr::C:ntM="<<ntM<<" <? tM="<<tM<<"+mPi0="<<mPi0<<" = "<<tM+mPi0<<G4endl;
#endif
    //if(ntM<tM+mPi0 && npM<pM+mPi0) elasFl=true; // Target&Project are underMinMass => El
    // @@ Just to try ---- End
#ifdef pdebug
    G4cout<<"G4QS::Cons: npM="<<npM<<" <? pM="<<pM<<"+mPi0="<<mPi0<<" = "<<pM+mPi0<<G4endl;
#endif
    if(npM<pM+mPi0) // The projectile is under min mass (@@ In Future put other cut @@)
    {
      G4double valk=tq4Mom.e();         // Momentum of the target quark
						//#ifdef pdebug
      G4double dvalk=valk-tq4Mom.rho();
      if(fabs(dvalk)>.00001) G4cout<<"***kp***G4QS::C:vk="<<valk<<",dvk="<<dvalk<<G4endl;
						//#endif
      G4double dmas=pM2-pr4Mom.m2();    // Difference of squared masses
      G4double vale=pr4Mom.e();
      G4ThreeVector pr=pr4Mom.v();
      G4double valp=pr.mag();
      G4ThreeVector upr=pr/valp;             // Unit vector in the pr direction
      G4double cost=-(dmas/(valk+valk)-vale)/valp; // Initial cos(theta)
      if(fabs(cost)>1.)
      {
#ifdef pdebug
        G4cout<<"***p***>>>G4QS::C: cost="<<cost<<G4endl;
#endif
        // Get max value of |cos(theta)| and change tq value to get the pM on mass shell
        if(cost>0.)                        // --->>> cost>0.
        {
          valk=dmas/(vale-valp)/2;         // Momentum is too big (must be reduced)
          G4ThreeVector tqf=upr*valk;      // final targetQuark 3-vector
          tq4Mom.set(valk,tqf);            // Fill new 4-mom for targetQuark
          tr4Mom.set(tM-valk,-tqf);        // Fill new 4-mom for targetResidual
          newProj4M=tq4Mom+pr4Mom;         // Final Projectile 4-mom in LS
        }
        else                               // --->>> cost<0.
								{
          G4double hvalk=dmas/(vale+valp); // Momentum's too small (must be increased, LIM)
          if(hvalk>tM)                     // Momentum can not be increased to this value
          {
#ifdef pdebug
            G4cout<<"**p-Cor**>G4QS::C: hvalk="<<hvalk<<" > tM="<<tM<<", dm="<<dmas<<", e="
                  <<vale<<", p="<<valp<<", ct="<<cost<<G4endl;
#endif
            // Put the scatteredProjectile on the massShell, and rescatter the target quark
            G4LorentzVector tmpProj4M=newProj4M+pq4Mom; // TempCriticalCompound for Project
            newProj4M=G4LorentzVector(0.,0.,0.,pM);     // New Project is on the mass shell
            G4QHadron tmpHadron(tmpProj4M);             // Create a Hadron for decay
            pq4Mom=G4LorentzVector(0.,0.,0.,0.);        // NewProjParton on lightMassShell
            tmpBl=tmpHadron.DecayIn2(newProj4M,pq4Mom); // Decay the Compound in newProg+pQ
            if(!tmpBl)G4cout<<"G4QS::C:DecIn2 err "<<sqrt(tmpProj4M.m2())<<"<"<<pM<<G4endl;
          }
          else
										{
#ifdef pdebug
            G4cout<<"***---***>G4QS::C: hvalk="<<hvalk<<" < tM="<<tM<<G4endl;
#endif
            valk=hvalk/2;                  // Momentum is too small (must be reduced)
            G4ThreeVector tqf=upr*(-valk); // Final targetQuark 3-vector
            tq4Mom.set(valk,tqf);          // Fill new 4-mom for targetQuark
            tr4Mom.set(tM-valk,-tqf);      // Fill new 4-mom for targetResidual
            newProj4M=tq4Mom+pr4Mom;
          }
        }
        if(fabs(newProj4M.m()-pM)>.0001)G4cout<<"G4QS::C:"<<newProj4M.m()<<","<<pM<<G4endl;
      }
      else
      {
        // Turn the target quark-parton to the projectile residual in LS = CMSofTheTarget
#ifdef pdebug
        G4cout<<"---p--->>>G4QS::C: cost="<<cost<<G4endl;
#endif
        G4ThreeVector tq=tq4Mom.v();           // 3-mom of the targetQ (LS)
        G4double      tqlp=upr.dot(tq);        // Projection of tq on the projectile
        G4ThreeVector tql=upr*tqlp;            // tq part along pr
        G4ThreeVector tqt=tq-tql;              // tq part perpendicular to pr
        G4double      cosi=tqlp/valk;          // Initial cos(theta)
        G4ThreeVector tqlf=tql*(cost/cosi);    // final tq part along pr
        G4double      sini=sqrt(1.-cosi*cosi); // Initial sin(theta)
        G4double      sint=sqrt(1.-cost*cost); // Desired sin(theta)
        G4ThreeVector tqpf=tqt*(sint/sini);    // final tq part perpendicular pr
        G4ThreeVector tqf=tqlf+tqpf;           // final tq
        tq4Mom.setV(tqf);                      // Change the momentum direction of targetQ
        tr4Mom.setV(-tqf);                     // Change the momentum direction of targetR
        if(fabs(tqf.mag()-valk)>.0001)G4cout<<"*G4QS::C:||="<<tqf.mag()<<","<<valk<<G4endl;
        newProj4M=tq4Mom+pr4Mom;
        if(fabs(newProj4M.m()-pM)>.001)G4cout<<"*G4QS::C:"<<newProj4M.m()<<","<<pM<<G4endl;
      }
#ifdef pdebug
      G4cout<<"G4QStr::C: Proj under GS, newP4M="<<newProj4M<<", pq4M="<<pq4Mom<<G4endl;
#endif
      pcorFl=true;                           // Projectile is on the GS mass shell
    }
    G4bool tcorFl=false;                     // ByDefault no MassSh correction for newTarg.
    G4LorentzVector newTarg4M=tr4Mom+pq4Mom;
    G4double nTM2=newTarg4M.m2();
    G4double ntM=sqrt(nTM2);
    //newTarg4M=tr4Mom+pq4Mom;
    //nTM2=newTarg4M.m2();
    //ntM=sqrt(nTM2);
#ifdef pdebug
    G4cout<<"G4QStr::C: ntM="<<ntM<<" <? tM="<<tM<<"+mPi0="<<mPi0<<" = "<<tM+mPi0<<G4endl;
#endif
    if(ntM<tM+mPi0 && !pcorFl) // The target is under min mass (@@ InFut put another cut@@)
    {
      // Turn the projectile quark-parton to the target rsidual in CMS of the Projectile
      G4LorentzVector pqc4M=pq4Mom;        // projectileQuark => projCM system <step1>
      pqc4M.boost(projRBoost);             // projectileQuark => projCM system <step2>
      G4double valk=pqc4M.e();             // Momentum of the target quark in projCM
						//#ifdef pdebug
      G4double dvalk=valk-pqc4M.rho();
      if(fabs(dvalk)>.00001) G4cout<<"***kt***G4QS::C:vk="<<valk<<",dvk="<<dvalk<<G4endl;
						//#endif
      G4double dmas=tM2-tr4Mom.m2();       // Difference of squared masses (targ - targRes)
      G4LorentzVector trc4M=tr4Mom;        // targetResidual => projCM system <step1>
      trc4M.boost(projRBoost);             // targetResidual => projCM system <step2>
      G4double vale=trc4M.e();             // Energy of the targetResidual in projCM
      G4ThreeVector tr=trc4M.v();          // 3-mom of the targetResidual in projCM
      G4double valp=tr.mag();              // momentum of the targetResidual in projCM
      if(fabs(dmas-tM2+trc4M.m2())>.1) G4cout<<"**t**G4QS::C: trM2="<<tr4Mom.m2()<<"="
                                            <<trc4M.m2()<<","<<vale*vale-valp*valp<<G4endl;
      G4ThreeVector utr=tr/valp;           // Unit vector in the tr direction in projCM
      G4double cost=-(dmas/(valk+valk)-vale)/valp; // Initial cos(theta)
      if(fabs(cost)>1.)
      {
#ifdef pdebug
        G4cout<<"***t***>>>G4QS::C: cost="<<cost<<G4endl;
#endif
        // Get max value of |cos(theta)| and change pq value to get the tM on mass shell
        if(cost>0.)                            // --->>> cost>0.
        {
          valk=dmas/(vale-valp)/2;             // Momentum is too big (must be reduced)
          G4ThreeVector pqf=utr*valk;          // final projectileQuark 3-vector
          G4LorentzVector pqc4M(valk,pqf);     // Fill new 4-mom for projectileQuark in CM
          pq4Mom=pqc4M;             // <step1> // Fill new 4-mom for projectileQuark in LS
          pq4Mom.boost(projBoost);  // <step2> // Fill new 4-mom for projectileQuark in LS
          G4LorentzVector prc4M(pM-valk,-pqf); // Fill new 4-mom for projectResidual in CM
          pr4Mom=prc4M;             // <step1> // Fill new 4-mom for projectResidual in LS
          pr4Mom.boost(projBoost);  // <step2> // Fill new 4-mom for projectResidual in LS
          newTarg4M=pq4Mom+tr4Mom;             // Final Target 4-mom in LS
        }
        else                                   // --->>> cost<-1 => cost=-1
								{
          G4double hvalk=dmas/(vale+valp); // Momentum's too small (must be increased, LIM)
          if(hvalk>pM)                     // Momentum can not be increased to this value
          {
#ifdef pdebug
            G4cout<<"**t-Cor**>G4QS::C: hvalk="<<hvalk<<" > pM="<<pM<<", dm="<<dmas<<", e="
                  <<vale<<", p="<<valp<<", ct="<<cost<<G4endl;
#endif
            // Put the scatteredProjectile on the massShell, rescatter the targetQuark (LS)
            G4LorentzVector tmpTarg4M=newTarg4M+tq4Mom; // TempCriticalCompound for Target
            newTarg4M=G4LorentzVector(0.,0.,0.,tM);     // New Target is on the mass shell
            G4QHadron tmpHadron(tmpTarg4M);             // Create a CompHadron for decay
            tq4Mom=G4LorentzVector(0.,0.,0.,0.);        // NewTargParton on lightMassShell
            tmpBl=tmpHadron.DecayIn2(newTarg4M,tq4Mom); // Decay the Compound in newTarg+tQ
            if(!tmpBl)G4cout<<"G4QS::C:DecIn2-err "<<sqrt(tmpTarg4M.m2())<<"<"<<tM<<G4endl;
          }
          else
										{
#ifdef pdebug
            G4cout<<"***---***>G4QS::C: hvalk="<<hvalk<<" < pM="<<pM<<G4endl;
#endif
            valk=hvalk/2;                        // Momentum is too small (mustBeIncreased)
            G4ThreeVector pqf=utr*(-valk);       // Final projectileQuark 3-vector in CM
            G4LorentzVector pqc4M(valk,pqf);     // Fill new 4-mom for projectQuark in CM
            pq4Mom=pqc4M;             // <step1> // Fill new 4-mom for projectQuark in LS
            pq4Mom.boost(projBoost);  // <step1> // Fill new 4-mom for projectQuark in LS
            G4LorentzVector nTc4M=pqc4M+trc4M;
#ifdef pdebug
            G4cout<<"G4QS::C:After Boost="<<projBoost<<G4endl;
            G4cout<<"G4QS::C:Aft: q="<<pqc4M<<pqc4M.m2()<<",t="<<nTc4M<<nTc4M.m2()<<G4endl;
            G4cout<<"G4QS::C:E="<<nTc4M.e()<<" = q="<<pqc4M.e()<<"="<<valk<<" + r="<<
                  trc4M.e()<<"="<<vale<<" = "<<pqc4M.e()+trc4M.e()<<",c="<<pq4Mom<<G4endl;
            if(fabs(nTc4M.m()-tM)>.0001) G4cout<<"*G4QS::C:"<<nTc4M.m()<<"!="<<tM<<G4endl;
            G4double pp=pqf.dot(tr);             // trc3M*pqc=-valk*valp
            G4double ee=vale*valk;
            G4cout<<"G4QS::C:tM2="<<tM2<<"="<<tM*tM<<",pp="<<pp*2<<"="<<-valk*valp*2<<",d="
                  <<ee-pp<<"="<<dmas<<"="<<(tM2-trc4M.m2())<<",u="<<utr.dot(utr)<<G4endl;
            G4double sen=nTc4M.e();
            G4double smo=nTc4M.rho();
            G4cout<<"G4QS::C:qM2="<<pqc4M.m2()<<",rM2="<<trc4M.m2()<<",E="<<sen<<"="<<
														vale+valk<<",P="<<smo<<"="<<valp-valk<<",M="<<sqrt(sen*sen-smo*smo)<<G4endl;
#endif
            G4LorentzVector prc4M(pM-valk,-pqf); // Fill new 4-mom for projectResid in CM
            pr4Mom=prc4M;             // <step1> // Fill new 4-mom for projectResid in LS
            pr4Mom.boost(projBoost);  // <step2> // Fill new 4-mom for projectResid in LS
            newTarg4M=pq4Mom+tr4Mom;
          }
        }
        if(fabs(newTarg4M.m()-tM)>.0001)
          G4cout<<"******************************G4QS::C:"<<newTarg4M.m()<<"="<<pM<<G4endl;
      }
      else // -1<cos(theta)<1 - only rotate the projQ
      {
        // Turn the projectile quark-parton to the target residual in CMS of the Projectile
#ifdef pdebug
        G4cout<<"---t--->>>G4QS::C: cost="<<cost<<G4endl;
#endif
        G4LorentzVector prc4M=pr4Mom;          // projectileQuark => projCM system <step1>
        prc4M.boost(projRBoost);               // projectileQuark => projCM system <step2>
        G4ThreeVector pq=pqc4M.v();            // 3-vector of the projQ in projCM
        G4double      pqlp=utr.dot(pq);        // Projection of pq on the target in projCM
        G4ThreeVector pql=utr*pqlp;            // pq part along pr
        G4ThreeVector pqt=pq-pql;              // pq part perpendicular to tr in projCM
        G4double      cosi=pqlp/valk;          // Initial cos(theta)
        G4ThreeVector pqlf=pql*(cost/cosi);    // final pq part along tr
        G4double      sini=sqrt(1.-cosi*cosi); // Initial sin(theta)
        G4double      sint=sqrt(1.-cost*cost); // Desired sin(theta)
        G4ThreeVector pqpf=pqt*(sint/sini);    // final pq part perpendicular tr
        G4ThreeVector pqf=pqlf+pqpf;           // final pq
        pqc4M.setV(pqf);                       // Change the momentumDirection of projQ(CM)
        pq4Mom=pqc4M;                          // Fill new 4-mom for projQuark in LS<step1>
        pq4Mom.boost(projBoost);               // Fill new 4-mom for projQuark in LS<step2>
        prc4M.setV(-pqf);                      // Change the momentumDirection of projR(CM)
#ifdef pdebug
        G4LorentzVector nT4M=pqc4M+trc4M;
        if(fabs(nT4M.m()-tM)>.001)G4cout<<"****t****G4QS::C:M="<<nT4M.m()<<"="<<tM<<G4endl;
#endif
        pr4Mom=prc4M;               // <step1> // Fill new 4-mom for projResidual in LS
        pr4Mom.boost(projBoost);    // <step2> // Fill new 4-mom for projResidual in LS
        if(fabs(pqf.mag()-valk)>.0001)G4cout<<"*G4QS::C:||="<<pqf.mag()<<"="<<valk<<G4endl;
        newTarg4M=pq4Mom+tr4Mom;
        if(fabs(newTarg4M.m()-tM)>.001)G4cout<<"*G4QS::C:"<<newTarg4M.m()<<"="<<tM<<G4endl;
      }
#ifdef pdebug
      G4cout<<"G4QStr::C: Targ under GS, newT4M="<<newTarg4M<<", tq4M="<<tq4Mom<<G4endl;
#endif
      tcorFl=true;                  // Target is on the GS mass shell
      newProj4M=pr4Mom+tq4Mom;      // Recalculate the Projectile (!)
      nPM2=newProj4M.m2();
      npM=sqrt(nPM2);
#ifdef pdebug
      G4cout<<"G4QStr::C:npM="<<npM<<" <? pM="<<pM<<"+mPi0="<<mPi0<<" = "<<pM+mPi0<<G4endl;
#endif
      if(npM<pM+mPi0) elasFl=true;    // Force elastic scattering (@@InFut putAnotherCut@@)
    }
    else if(ntM<tM+mPi0) elasFl=true; // Force elastScattering (@@InFut put another cut@@)
    if(elasFl)                        // Ellastic scattering happened
    {
      // **** Put the hadrons on the mass shell conserving the CMS scattering angle ****
      G4LorentzVector theTot4M=theProj4Mom+theTarg4Mom;// 4-momentum of CMS of "targ+proj"
		    G4ThreeVector cmsBoost = theTot4M.boostVector(); // CMS Boost Vector "CMS to LS"
      G4ThreeVector cmsRBoost= -cmsBoost;              // CMS Boost Vector "LS to CMS"
      G4LorentzVector cmsProj4M=theProj4Mom; // LS projectile => CMS <step1>
      cmsProj4M.boost(cmsRBoost);            // LS projectile => CMS <step2>
      G4LorentzVector cmsTarg4M=theTarg4Mom; // LS target => CMS <step1>
      cmsTarg4M.boost(cmsRBoost);            // LS target => CMS <step2>
      G4double pcm=cmsProj4M.rho();          // CMS momentum for the elastic scattering
      //#ifdef pdebug
      if(fabs(cmsTarg4M.rho()-pcm) > 0.0001)
        G4cout<<"-Worning-G4QSplitter::Constr: P="<<cmsTarg4M.rho()<<"#"<<pcm<<G4endl;
						//#endif
      G4LorentzVector cmsNewPr4M=newProj4M;  // LS finalProj => CMS <step1>
      cmsNewPr4M.boost(cmsRBoost);           // LS finalProj => CMS <step2>
      G4ThreeVector puV=cmsNewPr4M.v()/cmsNewPr4M.rho(); // Direction of the projectile
      //#ifdef pdebug
      G4LorentzVector cmsNewTg4M=newTarg4M;  // LS finalTarg => CMS <step1> @@ TMP
      cmsNewTg4M.boost(cmsRBoost);           // LS finalTarg => CMS <step2> @@ TMP
      G4ThreeVector tuV=cmsNewTg4M.v()/cmsNewTg4M.rho(); // Direction of the projectile
      if(1.+puV.dot(tuV) > 0.001)
        G4cout<<"-Worning-G4QSplitter::Constr: ct="<<puV.dot(tuV)<<G4endl;
						//#endif
      cmsProj4M.setV(puV*pcm);
      newProj4M=cmsProj4M;
      newProj4M.boost(cmsBoost);             // CMS FinalProjectile => LS <step2>
      G4QHadron* projH = new G4QHadron(projHadron); // Prototype of the Projectile Hadron
      projH->Set4Momentum(newProj4M);
      theQHadrons.push_back(projH);          // Fill theProjectileHadron(delete equivalent)
#ifdef pdebug
      G4cout<<"G4QSplitter::Constr: Fill ElastProjH="<<projQPDG<<newProj4M<<G4endl;
#endif
      cmsTarg4M.setV(puV*(-pcm));
      newTarg4M=cmsTarg4M;
      newTarg4M.boost(cmsBoost);             // CMS FinalTarget => LS <step2>
      G4QHadron* targH = new G4QHadron(targQPDG,newTarg4M); // Prototype of theTargetHadron
      theQHadrons.push_back(targH);          // Fill the Target Hadron (delete equivalent)
#ifdef pdebug
      G4cout<<"G4QSplitter::Constr: Fill ElastTargH="<<targQPDG<<newTarg4M<<G4endl;
#endif
    }
    else
    {
      // Inelastic scattering: one or both hadrons are excited (charge exchange is not in).
      if(pcorFl)                             // Projectile is on the mass shell
      {
        G4QHadron* projH = new G4QHadron(projHadron); // Prototype of the Projectile Hadron
        projH->Set4Momentum(newProj4M);
        theQHadrons.push_back(projH);        // Fill theProjectileHadron(delete equivalent)
        G4Quasmon* targQ = new G4Quasmon(targQPDG.GetQuarkContent(),newTarg4M);
        theQuasmons.push_back(targQ);        // Insert Projectile Quasmon (delete equival.)
      }
      if(tcorFl)                             // Target is on the mass shell
      {
        G4Quasmon* projQ = new G4Quasmon(projHadron.GetQC(),newProj4M);
        theQuasmons.push_back(projQ);        // Insert Projectile Quasmon (delete equival.)
        G4QHadron* targH = new G4QHadron(targQPDG,newTarg4M);//Prototype of theTargetHadron
        theQHadrons.push_back(targH);        // Fill the Target Hadron (delete equivalent)
      }
						else                                   // Both are excited (two Quasmons only) 
      {
        G4Quasmon* projQ = new G4Quasmon(projHadron.GetQC(),newProj4M);
        theQuasmons.push_back(projQ);        // Insert Projectile Quasmon (delete equival.)
        G4Quasmon* targQ = new G4Quasmon(targQPDG.GetQuarkContent(),newTarg4M);
        theQuasmons.push_back(targQ);        // Insert Projectile Quasmon (delete equival.)
      }
    }
  }
  else G4cout<<"-Wor-G4QSplitter::Constr:T="<<targFl<<" or P="<<projFl<<" err"<<G4endl;
#ifdef pdebug
  G4cout<<"G4QSplitter::Constructor: *** End of Constructor ***, elF="<<elasFl<<G4endl;
#endif
  // This is the first Check of the control values -- @@ Must be remade @@ --
#ifdef chdebug
  G4int finCharge=0;
  G4int finBaryoN=0;
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
  G4cout<<"G4QStr::C:nH="<<nHad<<",nQ="<<nQuas<<",C="<<finCharge<<",B="<<finBaryoN<<G4endl;
  if(finCharge!=totCharge || finBaryoN!=totBaryNum)
  {
    G4cerr<<"***G4QStr::C:tC="<<totCharge<<",C="<<finCharge<<",tB="<<totBaryNum
          <<",B="<<finBaryoN<<G4endl;
    if(nHad) for(G4int h=0; h<nHad; h++)
    {
      G4QHadron* cH = theQHadrons[h];
      G4cerr<<"::G4QS::C:h#"<<h<<",QC="<<cH->GetQC()<<",PDG="<<cH->GetPDGCode()<<G4endl;
    }
    if(nQuas) for(G4int q=0; q<nQuas; q++)
    {
      G4Quasmon* cQ = theQuasmons[q];
      G4cerr<<"::G4QS::C:q#"<<q<<",C="<<cQ->GetCharge()<<",QuarkCon="<<cQ->GetQC()<<G4endl;
    }
  }
#endif
} // End of the G4QSplitter constructor

G4QSplitter::G4QSplitter(const G4QSplitter &right)
{
  // theQuasmons (Vector)
  G4int nQ             = right.theQuasmons.size();
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
    G4Quasmon* curQ    = new G4Quasmon(right.theQuasmons[iq]);
#ifdef fdebug
    G4cout<<"G4QS::CopyByVal:Q#"<<iq<<","<<curQ->GetQC()<<curQ->Get4Momentum()<<G4endl;
#endif
    theQuasmons.push_back(curQ);             // (delete equivalent)
  }
  theProjEnvFlag  = right.theProjEnvFlag;
  theTargEnvFlag  = right.theTargEnvFlag;
  theWeight       = right.theWeight;
  theProjQC       = right.theProjQC;
  theTargQC       = right.theTargQC;
  theProj4Mom     = right.theProj4Mom;
  theTarg4Mom     = right.theTarg4Mom;
  
  theWorld        =  right.theWorld; 
		tot4Mom         =	 right.tot4Mom;
		totCharge       =	 right.totCharge;
		totBaryNum      =	 right.totBaryNum;
}

const G4QSplitter& G4QSplitter::operator=(const G4QSplitter &right)
{// ========================================================================
  if(this != &right)                          // Beware of self assignment
  {
    // theQuasmons (Vector)
    G4int iQ             = theQuasmons.size();
    if(iQ) for(G4int jq=0; jq<iQ; jq++) delete theQuasmons[jq];
    theQuasmons.clear();
    G4int nQ             = right.theQuasmons.size();
    if(nQ) for(G4int iq=0; iq<nQ; iq++)
    {
      G4Quasmon* curQ    = new G4Quasmon(right.theQuasmons[iq]);
#ifdef fdebug
      G4cout<<"G4QS::CopyByVal:Q#"<<iq<<","<<curQ->GetQC()<<curQ->Get4Momentum()<<G4endl;
#endif
      theQuasmons.push_back(curQ);             // (delete equivalent)
    }
    theProjEnvFlag  = right.theProjEnvFlag;
    theTargEnvFlag  = right.theTargEnvFlag;
    theWeight       = right.theWeight;
    theProjQC       = right.theProjQC;
    theTargQC       = right.theTargQC;
    theProj4Mom     = right.theProj4Mom;
    theTarg4Mom     = right.theTarg4Mom;
  
    theWorld        =  right.theWorld; 
		  tot4Mom         =	 right.tot4Mom;
		  totCharge       =	 right.totCharge;
		  totBaryNum      =	 right.totBaryNum;
  }
  return *this;
}

G4QSplitter::G4QSplitter(G4QSplitter* right)
{
  // theQuasmons (Vector)
  G4int nQ             = right->theQuasmons.size();
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
    G4Quasmon* curQ    = new G4Quasmon(right->theQuasmons[iq]);
#ifdef fdebug
    G4cout<<"G4QS::CopyByPoint:Q#"<<iq<<","<<curQ->GetQC()<<curQ->Get4Momentum()<<G4endl;
#endif
    theQuasmons.push_back(curQ);             // (delete equivalent)
  }
  theProjEnvFlag  = right->theProjEnvFlag;
  theTargEnvFlag  = right->theTargEnvFlag;
  theWeight       = right->theWeight;
  theProjQC       = right->theProjQC;
  theTargQC       = right->theTargQC;
  theProj4Mom     = right->theProj4Mom;
  theTarg4Mom     = right->theTarg4Mom;
  
  theWorld        =  right->theWorld; 
		tot4Mom         =	 right->tot4Mom;
		totCharge       =	 right->totCharge;
		totBaryNum      =	 right->totBaryNum;
}

G4QSplitter::~G4QSplitter()
{
#ifdef debug
  G4cout<<"~G4QSplitter: before theQHadrons nH="<<theQHadrons.size()<<G4endl;
#endif
  for_each(theQHadrons.begin(), theQHadrons.end(), DeleteQHadron());
#ifdef debug
  G4cout<<"~G4QSplitter: before theQuasmons nQ="<<theQuasmons.size()<<G4endl;
#endif
  for_each(theQuasmons.begin(), theQuasmons.end(), DeleteQuasmon());
#ifdef debug
  G4cout<<"~G4QSplitter: === DONE ==="<<G4endl;
#endif
}

// ================== Initialization of Static Parameters ============
//G4double G4QSplitter::StParName=0.;     // Static Parameter (Example)
//G4bool   G4QSplitter::stFlag=false;     // Static Flag (Example)

// Fill the private static parameters
//void G4QSplitter::SetParameters(G4double aStPar, G4bool aStFlag)
//{//  ====================================================================================
//  StParName=aStPar;       // (Example)
//  stFlag=aStFlag;         // (Example)
//}


//The public Hadronisation function with the Exception treatment (del respons. of User !)
G4QHadronVector* G4QSplitter::Fragment()
{//              ========================== -- @@ Must be changed @@ --
		// Make the final check before filling the output -- @@ Must be changed @@ --
#ifdef chdebug
  G4int fCharge=0;
  G4int fBaryoN=0;
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
  if(fCharge!=totCharge || fBaryoN!=totBaryNum)
  {
    G4cerr<<"***G4QS::Frag:tC="<<totCharge<<",C="<<fCharge<<",tB="<<totBaryNum
          <<",B="<<fBaryoN<<G4endl;
    if(nHad) for(G4int h=0; h<nHad; h++)
    {
      G4QHadron* cH = theQHadrons[h];
      G4cerr<<":G4QS::HQ:h#"<<h<<",QC="<<cH->GetQC()<<",PDG="<<cH->GetPDGCode()<<G4endl;
    }
    if(nQuas) for(G4int q=0; q<nQuas; q++)
    {
      G4Quasmon* cQ = theQuasmons[q];
      G4cerr<<":G4QS::HQ:q#"<<q<<",C="<<cQ->GetCharge()<<",QCont="<<cQ->GetQC()<<G4endl;
    }
  }
#endif
  G4QHadronVector dummy;       // Prototype of the output G4QuasmonVector to avoid wornings
  G4QHadronVector* finalQHadrons = &dummy; // Prototype of the output G4QuasmonVector
  G4int nH = theQHadrons.size();
  if(nH)
  {
    for(G4int ih=0; ih<nH; ih++)
    {
      G4QHadron* curH     = new G4QHadron(theQHadrons[ih]);
      finalQHadrons->push_back(curH);                  // deleted after the "while LOOP"
    }
  }
  for_each(theQuasmons.begin(),theQuasmons.end(),DeleteQuasmon()); //CleanUp Quasm's
  theQuasmons.clear();
  return finalQHadrons;
} // End of the Fragmentation member function

//The public extraction of the number of the created (in Constructor) G4QHadrons
G4int G4QSplitter::GetNOfHadrons() {return theQHadrons.size();}

//The public extraction of the number of the created (in Constructor) G4Quasmons
G4int G4QSplitter::GetNOfQuasmons() {return theQuasmons.size();}

//The public extraction of the projectile environment flag ("true" - exists)
G4bool G4QSplitter::GetProjEnvFlag() {return theProjEnvFlag;}

//The public extraction of the target environment flag ("true" - exists)
G4bool G4QSplitter::GetTargEnvFlag() {return theTargEnvFlag;}

//The public Quasmons duplication with delete responsibility of User (!)
G4QuasmonVector* G4QSplitter::GetQuasmons()
{//              ==============================
  G4int nQ=theQuasmons.size();
#ifdef pdebug
  G4cout<<"G4QSplitter::GetQuasmons is called nQ="<<nQ<<G4endl;
#endif
  G4QuasmonVector* quasmons = new G4QuasmonVector; // Intermediate
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
#ifdef pdebug
    G4cout<<"G4QStr::GetQuasmons:Q#"<<iq<<",QPDG="<<theQuasmons[iq]->GetQPDG()<<",QQC="
          <<theQuasmons[iq]->GetQC()<<",Q4M="<<theQuasmons[iq]->Get4Momentum()<<G4endl;
#endif
    if(theQuasmons[iq]->GetStatus())
    {
      G4Quasmon* curQ = new G4Quasmon(theQuasmons[iq]);
      quasmons->push_back(curQ);                   // (del. equiv. - user is responsibile)
    }
  }
#ifdef pdebug
  G4cout<<"G4QSplitter::GetQuasmons ===OUT==="<<G4endl;
#endif
  return quasmons;
} // End of GetQuasmons

//The public Quasmons duplication with delete responsibility of User (!)
G4QHadronVector* G4QSplitter::GetHadrons()
{//              =============================
  G4int nH=theQHadrons.size();
#ifdef pdebug
  G4cout<<"G4QSplitter::GetHadrons is called nH="<<nH<<G4endl;
#endif
  G4QHadronVector* hadrons = new G4QHadronVector; // Intermediate
  if(nH) for(G4int ih=0; ih<nH; ih++)
  {
#ifdef pdebug
    G4cout<<"G4QStr::GetHadrons: H#"<<ih<<",QPDG="<<theQHadrons[ih]->GetQPDG()<<",QQC="
          <<theQHadrons[ih]->GetQC()<<",H_Mass="<<theQHadrons[ih]->GetMass()<<G4endl;
#endif
    if(!theQHadrons[ih]->GetNFragments())
    {
      G4QHadron* curH = new G4QHadron(theQHadrons[ih]);
      hadrons->push_back(curH);                   // (del. equiv. - user is responsibile)
    }
  }
#ifdef pdebug
  G4cout<<"G4QSplitter::GetHadrons ===OUT==="<<G4endl;
#endif
  return hadrons;
} // End of GetHadrons

// Randomize the momentum fraction for the CHIPS of nPart free partons [x*(1-x)^(n-2)]
G4double G4QSplitter::RandomizeMomFractionFree(G4int nPart)
{//              ==============================================
  // @@ TMP --- Begin ---
  if(2>1)
		{
    if(nPart<2)
    {
      G4cerr<<"**G4QSplitter::RandMomFractionString: n="<<nPart<<" < 2, retun 0"<<G4endl;
      return 0.;
    }
    G4double x=0.5;
    if(nPart==2) return x;          // GS meson
    G4double r=G4UniformRand();
    if(r==0.) return 0.;
    if(r==1.) return 1.;
		  if     (nPart==3) x=r;          // GS baryon
		  else if(nPart==4) x=1.-sqrt(r); // GS quaternion
    else x=1.-pow(r,1./(nPart-2.)); // nPart>4
    return x;
  }
  //if(nPart!=3) G4cerr<<"*>>>>>*G4QSplitter::RandMomFractionFree: n="<<nPart<<G4endl;
  // @@ TMP --- End ---
  if(nPart<2)
  {
    G4cerr<<"**G4QSplitter::RandMomFractionFree: n="<<nPart<<" < 2, retun 0."<<G4endl;
    return 0.;
  }
  G4double x=0.5;
  if(nPart==2) return x;       // GS meson
  G4double r=G4UniformRand();
  if(r==0.) return 0.;
  if(r==1.) return 1.;
		if(nPart==3) x=sqrt(r);      // GS baryon
		else if(nPart==4)            // GS quaternion
		{
    if    (r==0.5) x=0.5;
    else if(r<0.5) x=sqrt(r+r)*(.5+.1579*(r-.5));
    else           x=1.-sqrt(2.-r-r)*(.5+.1579*(.5-r));
  }
  else
		{
    G4int n1=nPart-2;
    G4double r1=n1;
    G4double r2=r1-1.;
    G4double rr=r2/r1;
    G4double rp=pow(rr,n1);
    G4double p2=rp+rp;
    if  (r==rr)  x=p2;
    else
				{
      if(r<rr)
      {
								G4double pr=0.;
								G4double pra=0.;
        if(nPart>8)
								{
          if(nPart>10)
								  {
            if(nPart>11)                      // >11
            {
              pr=.614/pow((nPart+1.25),.75);
              pra=.915/pow((nPart+6.7),1.75);
            }
												else                              // 11
            {
              pr=.09945;
              pra=.00667;
            }
          }
          else
								  {
            if(nPart>9)                       // 10
            {
              pr=.1064;
              pra=.00741;
            }
												else                              // 9
            {
              pr=.11425;
              pra=.00828;
            }
          }
        }
        else
								{
          if(nPart>6)
								  {
            if(nPart>7)                       // 8
            {
              pr=.12347;
              pra=.00926;
            }
												else                              // 7
            {
              pr=.13405;
              pra=.01027;
            }
          }
          else
								  {
            if(nPart>5)                       // 6
            {
              pr=.1454;
              pra=.01112;
            }
												else                              // 5
            {
              pr=.15765;
              pra=.00965;
            }
          }
        }
        x=pow((r/p2),(1.-rr+pra))*(rr+pr*(r-p2));
      }
      else
      {
								G4double sr=0.;
        if(nPart>8)
								{
          if(nPart>10)
								  {
            if(nPart>11) sr=.86/(nPart+1.05); // >11
												else         sr=.0774;            // 11
          }
          else
								  {
            if(nPart>9) sr=.0849;             // 10
												else        sr=.0938;             // 9
          }
        }
        else
								{
          if(nPart>6)
								  {
            if(nPart>7) sr=.1047;             // 8
												else        sr=.1179;             // 7
          }
          else
								  {
            if(nPart>5) sr=.1339;             // 6
												else        sr=.15135;            // 5
          }
        }
        x=1.-sqrt((1.-r)/(1.-p2))*(1.-rr+sr*(p2-r));
      }
    }
  }
  return x;
} // End of RandomizeMomFractionFree

// Randomize MomentumFraction for CHIPS of nPart partons for Pomeron exchange [(1-x)^(n-2)]
G4double G4QSplitter::RandomizeMomFractionString(G4int nPart)
{//              ================================================
  if(nPart<2)
  {
    G4cerr<<"**G4QSplitter::RandMomFractionString: n="<<nPart<<" < 2, retun 0"<<G4endl;
    return 0.;
  }
  G4double x=0.5;
  if(nPart==2) return x;          // GS meson
  G4double r=G4UniformRand();
  if(r==0.) return 0.;
  if(r==1.) return 1.;
		if     (nPart==3) x=r;          // GS baryon
		else if(nPart==4) x=1.-sqrt(r); // GS quaternion
  else x=1.-pow(r,1./(nPart-2.)); // nPart>4
  return x;
} // End of RandomizeMomFractionString
