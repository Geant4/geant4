//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//       1         2         3         4         5         6         7         8         9
//34567890123456789012345678901234567890123456789012345678901234567890123456789012345678901
//
//
// $Id: G4QuasmonString.cc,v 1.1 2004-12-14 16:01:20 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QuasmonString ----------------
//             by Mikhail Kossov, August 2000.
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

#include "G4QuasmonString.hh"
#include <cmath>
using namespace std;
 

G4QuasmonString::G4QuasmonString(G4QHadron projHadron, const G4bool projEnvFlag,
                     const G4int targPDG, const G4bool targEnvFlag) :
  theProjEnvFlag(projEnvFlag), theTargEnvFlag(targEnvFlag)
{
  static const G4double mPi0 = G4QPDGCode(111).GetMass(); // Pi0 mass
  //static const G4QPDGCode pimQPDG(-211);
  theWorld= G4QCHIPSWorld::Get();             // Get a pointer to the CHIPS World
#ifdef debug
  G4cout<<">>>>>>G4QuasmonString::Constructor:Called pH="<<projHadron<<",tPDG="<<targPDG<<G4endl;
#endif
  // Target initialization
  G4QPDGCode targQPDG(targPDG);
  totBaryNum= targQPDG.GetBaryNum();          // Prototype
  totCharge = targQPDG.GetCharge();           // Prototype
  theTargQC = targQPDG.GetQuarkContent();
  G4double tM = targQPDG.GetMass();
  theTarg4Mom = G4LorentzVector(0.,0.,0.,tM);
  G4QHadron targHadron(targQPDG,theTarg4Mom); // Target Hadron

  // Projectile initialization
  theProj4Mom = projHadron.Get4Momentum();
  G4QPDGCode projQPDG(projHadron.GetPDGCode());
  theTargQC  = projQPDG.GetQuarkContent();
  totBaryNum+= projQPDG.GetBaryNum();         // Final control value
  totCharge += projQPDG.GetCharge();          // Final control value
  tot4Mom  = theProj4Mom+theTarg4Mom;         // Final control value

  // === Print out of the input information at Creation time & tot 4-mom Calculation ======
#ifdef pdebug
  G4cout<<"G4QS::Cons:PQC="<<theProjQC<<",TQC="<<theTargQC<<",P4Mom="<<theProj4Mom<<G4endl;
  G4cout<<"G4QStr::Constr: tC="<<totCharge<<",tB="<<totBaryNum<<",tot4M="<<tot4Mom<<G4endl;
#endif
  //G4int nP=theWorld->GetQPEntries();       // A#of init'ed particles in CHIPS World (@@?)
  //G4int nCl=nP-90;                          // A#of init'ed clusters in CHIPS World (@@?)
		//#ifdef pdebug
  //G4cout<<"G4QS:Const:Before QEX:n="<<nP<<G4endl;
		//#endif
  // @@@@@@ ===> Here the Quark Exchange Quasmon Creation must be added <=== @@@@@@
  G4int nQTarg=theTargQC.GetTot();           // a#of Quark-partons in the Target
  G4int nQProj=theProjQC.GetTot();           // a#of Quark-partons in the Progectile
  G4double projM = projQPDG.GetMass();
  G4double projM2 = projM*projM;
  G4double pM2 = theProj4Mom.m2();
  G4double tM2 = tM*tM;
  G4double dM2 = fabs(pM2-projM2);
  if(dM2>.01) G4cout<<"-Wor-G4QS::Constr:dM2="<<dM2<<",M2="<<projM2<<",LVM2="<<pM2<<G4endl;
  G4double pM=sqrt(pM2);                     // @@ do we need pM ? @@ (in print)
  // @@ Now projectile can be only meson or baryon @@ -- @@ Improve for clusters @@ --
  G4LorentzVector pq4Mom(0.,0.,0.,0.);       // Prototype of LV of quark of progectile  
  G4double rPMass=0.;                        // Prototype of the residProjMass (Meson case)
  if(nQProj<2)  G4cout<<"***G4QSt::Constr: nQProj="<<nQProj<<"<2 ***Fata error***"<<G4endl;
		else if(nQProj>2)                          // ---> Baryon case (clusters are not implem.)
		{
    if(nQProj>3) G4cout<<"-Wor-G4QS::Constr:nQProj="<<nQProj<<">3 is not implem'd"<<G4endl;
    G4double x=RandomizeMomentumFraction(3);
    rPMass = sqrt(pM2*(1.-x));               // Residual Projectile mass
  }
  G4LorentzVector pr4Mom(0.,0.,0.,rPMass);   // Prototype of LV of the residual projectile
  G4bool projFl=false;
  G4bool targFl=false;
  G4bool tmpBl=projHadron.DecayIn2(pq4Mom,pr4Mom);
  if(tmpBl) projFl=true;
  //if(projHadron.DecayIn2(pq4Mom,pr4Mom)) projFl=true;
  else G4cout<<"**Worning**G4QuasmonString::Constr:ProjDecIn2 rpM="<<rPMass<<", pM="<<pM<<G4endl;
  G4LorentzVector tq4Mom(0.,0.,0.,0.);       // Prototype of LV of quark of the target  
  if(nQTarg<3) G4cout<<"***G4QStr::Constr: nQTarg="<<nQTarg<<"<3 ***Fata error***"<<G4endl;
  if(nQTarg>3) G4cout<<"-W-G4QS::Constr: nQTarg="<<nQProj<<">3 is not implemented"<<G4endl;
  G4double rTMass = sqrt(tM2*(1.-RandomizeMomentumFraction(3))); // Residual Target mass
  G4LorentzVector tr4Mom(0.,0.,0.,rTMass);   // Prototype of LV of the residual projectile
  if(targHadron.DecayIn2(tq4Mom,tr4Mom)) targFl=true;
  else G4cout<<"**Worning**G4QStr::Constr:TargDecIn2 rtM="<<rTMass<<", tM="<<tM<<G4endl;
  G4bool elasFl=false;                          // ByDefault avoid the elastic scattering
  if (targFl && projFl)                      // --- @@ Now only for pp case @@ ---
  {
    G4LorentzVector newProj4M=pr4Mom+tq4Mom;
    G4double nPM2=newProj4M.m2();
    G4double npM=sqrt(nPM2);
    G4bool pcorFl=false;                          // ByDefault no correction for newProj.
    if(npM<pM+mPi0) // The projectile is under min mass (@@ In Future put other cut @@)
    {
      // Put the scattered projectile on the mass shell, and rescatter the target quark
      G4LorentzVector tmpProj4M=newProj4M+pq4Mom; // Temporary critical compound for proj.
      newProj4M=G4LorentzVector(0,0,0,pM);        // New proj is on the mass shell
      G4QHadron tmpHadron(tmpProj4M);
      tmpBl=tmpHadron.DecayIn2(newProj4M,pq4Mom);
      if(!tmpBl)G4cout<<"*W*G4QSt::C:DecIn2 M="<<sqrt(tmpProj4M.m2())<<", pM="<<pM<<G4endl;
      pcorFl=true;
    }
    G4LorentzVector newTarg4M=tr4Mom+pq4Mom;
    G4double nTM2=newTarg4M.m2();
    G4double ntM=sqrt(nTM2);
    if(ntM<tM+mPi0 && !pcorFl) // The target is under min mass (@@ InFut put another cut@@)
    {
      // Put the scattered projectile on the mass shell, and rescatter the target quark
      G4LorentzVector tmpTarg4M=newTarg4M+tq4Mom; // Temporary critical compound for target
      newTarg4M=G4LorentzVector(0,0,0,tM);        // New target is on the mass shell
      G4QHadron tmpHadron(tmpTarg4M);
      tmpBl=tmpHadron.DecayIn2(newTarg4M,pq4Mom);
      if(!tmpBl)G4cout<<"*W*G4QSt::C:DecIn2 M="<<sqrt(tmpTarg4M.m2())<<", tM="<<tM<<G4endl;
      newProj4M=pr4Mom+tq4Mom;
      nPM2=newProj4M.m2();
      npM=sqrt(nPM2);
      if(npM<pM+mPi0) elasFl=true; // Force elastic scattering (@@InFut put another cut@@)
    }
    else if(ntM<tM+mPi0) elasFl=true; // Force elastScattering (@@InFut put another cut@@)

  }
  else G4cout<<"**Worning**G4QuasmonString::Constr:T="<<targFl<<" or P="<<projFl<<" err"<<G4endl;
  // This is the first Check of the control values -- @@ Must be remade @@ --
  // @@@@@@@@@@@@@@ Temporary for the testing purposes --- Begin
  theProjEnvFlag=elasFl;
  // @@@@@@@@@@@@@@ Temporary for the testing purposes --- End
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
    G4cerr<<"***G4QStr::C:tC="<<totCharge<<",C="<<finCharge<<",tB="<<totBaryoN
          <<",B="<<finBaryoN<<",E="<<theEnvironment<<G4endl;
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
} // End of the G4QuasmonString constructor

G4QuasmonString::G4QuasmonString(const G4QuasmonString &right)
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
  theProjQC       = right.theProjQC;
  theTargQC       = right.theTargQC;
  theProj4Mom     = right.theProj4Mom;
  theTarg4Mom     = right.theTarg4Mom;
  
  theWorld        =  right.theWorld; 
		tot4Mom         =	 right.tot4Mom;
		totCharge       =	 right.totCharge;
		totBaryNum      =	 right.totBaryNum;
}

const G4QuasmonString& G4QuasmonString::operator=(const G4QuasmonString &right)
{// ========================================================================
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
  theProjQC       = right.theProjQC;
  theTargQC       = right.theTargQC;
  theProj4Mom     = right.theProj4Mom;
  theTarg4Mom     = right.theTarg4Mom;
  
  theWorld        =  right.theWorld; 
		tot4Mom         =	 right.tot4Mom;
		totCharge       =	 right.totCharge;
		totBaryNum      =	 right.totBaryNum;

  return *this;
}

G4QuasmonString::G4QuasmonString(G4QuasmonString* right)
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
  theProjQC       = right->theProjQC;
  theTargQC       = right->theTargQC;
  theProj4Mom     = right->theProj4Mom;
  theTarg4Mom     = right->theTarg4Mom;
  
  theWorld        =  right->theWorld; 
		tot4Mom         =	 right->tot4Mom;
		totCharge       =	 right->totCharge;
		totBaryNum      =	 right->totBaryNum;
}

G4QuasmonString::~G4QuasmonString()
{
#ifdef debug
  G4cout<<"~G4QuasmonString: before theQHadrons nH="<<theQHadrons.size()<<G4endl;
#endif
  for_each(theQHadrons.begin(), theQHadrons.end(), DeleteQHadron());
#ifdef debug
  G4cout<<"~G4QuasmonString: before theQuasmons nQ="<<theQuasmons.size()<<G4endl;
#endif
  for_each(theQuasmons.begin(), theQuasmons.end(), DeleteQuasmon());
#ifdef debug
  G4cout<<"~G4QuasmonString: === DONE ==="<<G4endl;
#endif
}

// ================== Initialization of Static Parameters ============
//G4double G4QuasmonString::StParName=0.;     // Static Parameter (Example)
//G4bool   G4QuasmonString::stFlag=false;     // Static Flag (Example)

// Fill the private static parameters
//void G4QuasmonString::SetParameters(G4double aStPar, G4bool aStFlag)
//{//  ====================================================================================
//  StParName=aStPar;       // (Example)
//  stFlag=aStFlag;         // (Example)
//}


//The public Hadronisation function with the Exception treatment (del respons. of User !)
G4QHadronVector* G4QuasmonString::Fragment()
{//              ========================== -- @@ Must be changed @@ --
		// Make the final check before filling the output -- @@ Must be changed @@ --
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
    G4cerr<<"***G4QS::Frag:tC="<<totCharge<<",C="<<fCharge<<",tB="<<totBaryoN
          <<",B="<<fBaryoN<<",E="<<theEnvironment<<G4endl;
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

//The public Quasmons duplication with delete responsibility of User (!)
G4QuasmonVector* G4QuasmonString::GetQuasmons()
{//              =========================
  G4int nQ=theQuasmons.size();
#ifdef pdebug
  G4cout<<"G4QuasmonString::GetQuasmons is called nQ="<<nQ<<G4endl;
#endif
  G4QuasmonVector* quasmons = new G4QuasmonVector; // Intermediate
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
#ifdef pdebug
    G4cout<<"G4QStr::GetQuasmons:Q#"<<iq<<",QPDG="<<theQuasmons[iq]->GetQPDG()<<",QQC="
          <<theQuasmons[iq]->GetQC()<<",QM="<<theQuasmons[iq]->GetMass()<<G4endl;
#endif
    G4Quasmon* curQ = new G4Quasmon(theQuasmons[iq]);
    quasmons->push_back(curQ);                     // (del. equiv. - user is responsibile)
  }
#ifdef pdebug
  G4cout<<"G4QuasmonString::GetQuasmons ===OUT==="<<G4endl;
#endif
  return quasmons;
} // End of GetQuasmons

//The public Quasmons duplication with delete responsibility of User (!)
G4QHadronVector* G4QuasmonString::GetHadrons()
{//              =========================
  G4int nH=theQHadrons.size();
#ifdef pdebug
  G4cout<<"G4QuasmonString::GetHadrons is called nH="<<nH<<G4endl;
#endif
  G4QHadronVector* hadrons = new G4QHadronVector; // Intermediate
  if(nH) for(G4int ih=0; ih<nH; ih++)
  {
#ifdef pdebug
    G4cout<<"G4QStr::GetHadrons: H#"<<ih<<",QPDG="<<theQHadrons[ih]->GetQPDG()<<",QQC="
          <<theQHadrons[ih]->GetQC()<<",H_Mass="<<theQHadrons[ih]->GetMass()<<G4endl;
#endif
    G4QHadron* curH = new G4QHadron(theQHadrons[ih]);
    hadrons->push_back(curH);                       // (del. equiv. - user is responsibile)
  }
#ifdef pdebug
  G4cout<<"G4QuasmonString::GetHadrons ===OUT==="<<G4endl;
#endif
  return hadrons;
} // End of GetHadrons

//The public Quasmons duplication with delete responsibility of User (!)
G4double G4QuasmonString::RandomizeMomentumFraction(G4int nPart)
{//              =========================================
  if(nPart<2)
  {
    G4cerr<<"**G4QuasmonString::RandomizeMomFraction: n="<<nPart<<" < 2, retun 0."<<G4endl;
    return 0.;
  }
  G4double x=0.5;
  if(nPart==2) return x;       // GS meson
  G4double r=G4UniformRand();
  if(r==0.) return 0.;
  if(r==1.) return 1.;
		if(nPart==3)                 // GS baryon
		{
    if    (r==0.5) x=0.5;
    else if(r<0.5) x=sqrt(r+r)*(.5+.1579*(r-.5));
    else           x=1.-sqrt(2.-r-r)*(.5+.1579*(.5-r));
  }
  else
		{
    G4int n1=nPart-1;
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
        if(nPart>7)
								{
          if(nPart>9)
								  {
            if(nPart>10)                      // >10
            {
              pr=.614/pow((nPart+1.25),.75);
              pra=.915/pow((nPart+6.7),1.75);
            }
												else                              // 10
            {
              pr=.09945;
              pra=.00667;
            }
          }
          else
								  {
            if(nPart>8)                       // 9
            {
              pr=.1064;
              pra=.00741;
            }
												else                              // 8
            {
              pr=.11425;
              pra=.00828;
            }
          }
        }
        else
								{
          if(nPart>5)
								  {
            if(nPart>6)                       // 7
            {
              pr=.12347;
              pra=.00926;
            }
												else                              // 6
            {
              pr=.13405;
              pra=.01027;
            }
          }
          else
								  {
            if(nPart>4)                       // 5
            {
              pr=.1454;
              pra=.01112;
            }
												else                              // 4
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
        if(nPart>7)
								{
          if(nPart>9)
								  {
            if(nPart>10) sr=.86/(nPart+1.05); // >10
												else         sr=.0774;            // 10
          }
          else
								  {
            if(nPart>8) sr=.0849;             // 9
												else        sr=.0938;             // 8
          }
        }
        else
								{
          if(nPart>5)
								  {
            if(nPart>6) sr=.1047;             // 7
												else        sr=.1179;             // 6
          }
          else
								  {
            if(nPart>4) sr=.1339;             // 5
												else        sr=.15135;            // 4
          }
        }
        x=1.-sqrt((1.-r)/(1.-p2))*(1.-rr+sr*(p2-r));
      }
    }
  }
  return x;
} // End of GetQuasmons
