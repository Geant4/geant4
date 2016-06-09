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
// $Id: G4QEnvironment.cc,v 1.79.2.1 2004/02/08 10:53:26 hpw Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
//      ---------------- G4QEnvironment ----------------
//             by Mikhail Kossov, August 2000.
//  class for Multy Quasmon Environment used by the CHIPS Model
// ------------------------------------------------------------
 
//#define debug
//#define pdebug
//#define rdebug
//#define edebug
//#define fdebug
//#define ppdebug
//#define sdebug

#include "G4QEnvironment.hh" 
#include "G4HadronicException.hh"

G4QEnvironment::G4QEnvironment(const G4QHadronVector& projHadrons, const G4int targPDG) :
  theEnvironment(90000000)                   // User is responsible for projHadrons(Vector)
{
  G4int nHadrons=projHadrons.size();         // A#of hadrons in the input Vector
  if(!nHadrons||targPDG==90000000)           // No projectile Hadrons or no target Nucleus
  {
    G4cerr<<"---Worning---G4QEnv::Const:a#ofINPHadr="<<nHadrons<<",tPDG="<<targPDG<<G4endl;
    //throw G4QException("***G4QEnvironment: There is no one projectile or vacuum target");
    if(nHadrons)                             // The projectiles are copied to the output
	{
      for(G4int ih=0; ih<nHadrons; ih++)
      {
        G4QHadron* curQH    = new G4QHadron(projHadrons[ih]);
#ifdef fdebug
        G4cout<<"*G4QE::Const:iH#"<<ih<<","<<curQH->GetQC()<<curQH->Get4Momentum()<<G4endl;
#endif
        theQHadrons.push_back(curQH);        // (delete equivalent)
      }
	}
    else if(targPDG!=90000000)               // No projHadrons,fill targetNucleus to output
    {
      G4QHadron* curQH    = new G4QHadron(targPDG);
#ifdef fdebug
        G4cout<<"**G4QE::Const:No iHad,eH="<<curQH->GetQC()<<curQH->Get4Momentum()<<G4endl;
#endif
      theQHadrons.push_back(curQH);          // (delete equivalent)
    }
    return;
  }
  G4int    targA=G4QPDGCode(targPDG).GetBaryNum();
  G4double targM=G4QPDGCode(targPDG).GetMass();
  tot4Mom=G4LorentzVector(0.,0.,0.,targM);
  // === Print out of the input information at Creation time & tot 4-mom Calculation ======
#ifdef pdebug
  G4cout<<"G4QEnvironment::Const:tPDG="<<targPDG<<",tM="<<targM<<",nPr="<<nHadrons<<G4endl;
#endif
  for(G4int ipr=0; ipr<nHadrons; ipr++)// LOOP is used for the tot4Mom calc. & for printing
  {
    G4LorentzVector h4Mom = projHadrons[ipr]->Get4Momentum();
    tot4Mom      += h4Mom;
#ifdef pdebug
    G4int           hPDG  = projHadrons[ipr]->GetPDGCode();
    G4int           hNFrag= projHadrons[ipr]->GetNFragments();
    G4QContent      hQC   = projHadrons[ipr]->GetQC();
    G4cout<<"G4E::C:#"<<ipr<<",PDG="<<hPDG<<hQC<<",4M="<<h4Mom<<",hNFrag="<<hNFrag<<G4endl;
#endif
  }
#ifdef pdebug
  G4cout<<"G4QEnv::Const:nF="<<projHadrons[0]->GetNFragments()<<",tot4M="<<tot4Mom<<G4endl;
#endif
  theWorld= G4QCHIPSWorld::Get();            // Get a pointer to the CHIPS World
  G4int nP=theWorld->GetQPEntries();         // A#of init'ed particles in CHIPS World
  G4int nCl=nP-80;                           // A#of init'ed clusters in CHIPS World
#ifdef pdebug
  G4cout<<"G4QEnv:Const:Before NCI:n="<<nP<<",F="<<projHadrons[0]->GetNFragments()<<G4endl;
#endif
  InitClustersVector(nCl,targA);             // Init Clusters as Particles (to interact)
#ifdef pdebug
  G4cout<<"G4QEnv::Const:NucClust,n="<<nCl<<",F="<<projHadrons[0]->GetNFragments()<<G4endl;
#endif
  if(targPDG>80000000)                        // ==> Nuclear target (including NUCPDG)
  {
    theEnvironment.InitByPDG(targPDG);        // Create nuclear environment
#ifdef pdebug
      G4cout<<"G4QEnv::Const:nH="<<nHadrons<<",PDG="<<projHadrons[0]->GetPDGCode()<<G4endl;
#endif
    if(nHadrons==1 && projHadrons[0]->GetPDGCode()==22)
	{
      G4double exMass=tot4Mom.m();
#ifdef pdebug
      G4cout<<"G4QEnvironment::Const: exM="<<exMass-targM<<" > mPi0 ?"<<G4endl;
#endif      
      if(exMass<targM+135.977) // Nucleus is below the pion threshold
	  {
        G4QNucleus exEnviron(tot4Mom,targPDG);
        // One can put here the pbpt= (M.K.)
        if(!exEnviron.SplitBaryon())   // Nucleus is below the splitting fragment threshold
		{
#ifdef pdebug
          G4cout<<"G4QEnv::Const:Photon is added to the Output,E="<<theEnvironment<<G4endl;
#endif      
          G4QHadron* photon = new G4QHadron(projHadrons[0]); // Fill prog photon to output
#ifdef fdebug
          G4cout<<"**G4QE::Const:photon="<<photon->GetQC()<<photon->Get4Momentum()<<G4endl;
#endif
          theQHadrons.push_back(photon);      // (delete equivalent)
          return;
		}
	  }
    }
    for(G4int ih=0; ih<nHadrons; ih++)        // ==> The main LOOP over projQHadrons
    {
      G4QHadron* curHadr=projHadrons[ih];     // Pointer to current projectile Hadron
      G4int hNFrag = curHadr->GetNFragments();// #0 means intermediate (skip)
      G4LorentzVector ch4M=curHadr->Get4Momentum(); // 4-momenyum of the current projectile
#ifdef pdebug
      G4cout<<"G4QE:C:"<<ih<<",F="<<hNFrag<<",0="<<projHadrons[0]->GetNFragments()<<G4endl;
#endif
      if(!hNFrag&&ch4M.e()>0.)                // => "Final hadron" case
	  {
        G4int envPDG=theEnvironment.GetPDG();
        if(envPDG==90000000)                  // ==> "Interaction with vacuum" case
        {
          G4int hPDG  = curHadr->GetPDGCode();// A PDG Code of the projQHadron
          if(!hPDG||hPDG==10)          // Check for the validity of the QHadron (@@ 10 OK?)
          {
            G4cerr<<"---Worning---G4QEnvironment::Constructor: wrong PDG("<<ih<<")="<<hPDG
                  <<", HQC="<<curHadr->GetQC()<<", HM="<<curHadr->GetMass()<<G4endl;
            //throw G4QException("***G4QEnvironment::Constructor: theInputHadron is Chip");
          }
          else
          {
            G4int hQ = curHadr->GetQCode();  // One more check for valid of the QHadron
            if(hQ<0)
	        {
              G4cerr<<"---Worning---G4QEnv::Constructor:Q<0, PDG=("<<ih<<")"<<hPDG<<G4endl;
              //throw G4QException("***G4QEnvironment::Constructor:theInputHadron is bad");
	        }
            else
            {
              G4QHadron* newHadr = new G4QHadron(curHadr);
#ifdef fdebug
              G4cout<<"*G4QE::Const:H="<<newHadr->GetQC()<<newHadr->Get4Momentum()<<G4endl;
#endif
              theQHadrons.push_back(newHadr); // Fill existing hadron (delete equivalent)
#ifdef pdebug
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
#ifdef pdebug
          G4cout<<"G4QE::Const:CreateQuasm, 4M="<<ch4M<<",QC="<<hQC<<",E="<<envPDG<<G4endl;
#endif
          CreateQuasmon(hQC, ch4M);
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
      G4cerr<<"---Worning---G4QEnvironment::Constructor:Vacuum,1st Hadron wrong PDG="<<hPDG
            <<", HQC="<<curHadr->GetQC()<<", HM="<<curHadr->GetMass()<<G4endl;
      //throw G4QException("***G4QEnvironment::Constructor: Fiest input Hadron is wrong");
    }
    else
    {
      G4int hQ = curHadr->GetQCode();         // One more check for valid of the QHadron
      if(hQ<0)
	  {
        G4cerr<<"---Worning---G4QEnviron::Constructor:Vacuum,Q<0, 1st HPDG="<<hPDG<<G4endl;
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
          G4cerr<<"---Worning---G4QEnv::Constructor:TrgQC<0, Chipo?,PDG="<<targPDG<<G4endl;
          //throw G4QException("***G4QEnvironment::Constructor: Target is wrong");
		}
        else                                 // Now we can create a unique Quasmon
		{
          h4Mom+=G4LorentzVector(0.,0.,0.,tQPDG.GetMass()); //Projectile + TargetHadron
          hQC+=tQPDG.GetQuarkContent();
#ifdef pdebug
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
#ifdef fdebug
      G4cout<<"*G4QE::Const:#"<<ih<<","<<curHadr->GetQC()<<curHadr->Get4Momentum()<<G4endl;
#endif
      theQHadrons.push_back(newHadr);        // Fill existing hadron (delete equivalent)
	}
  } // End of Unique Hadron target treatment
} // End of the G4QEnvironment constructor

G4QEnvironment::G4QEnvironment(const G4QEnvironment &right)
{
  // theQHadrons (Vector)
  G4int nQH             = right.theQHadrons.size();
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
    G4QHadron* curQH    = new G4QHadron(right.theQHadrons[ih]);
#ifdef fdebug
    G4cout<<"G4QE::CopyByVal:cH#"<<ih<<","<<curQH->GetQC()<<curQH->Get4Momentum()<<G4endl;
#endif
    theQHadrons.push_back(curQH);            // (delete equivalent)
  }

  theWorld              = right.theWorld;
  nBarClust             = right.nBarClust;
  f2all                 = right.f2all;
  tot4Mom               = right.tot4Mom;

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
  G4int nQH             = right->theQHadrons.size();
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
    G4QHadron* curQH    = new G4QHadron(right->theQHadrons[ih]);
#ifdef fdebug
    G4cout<<"G4QE::CopyByPtr:cH#"<<ih<<","<<curQH->GetQC()<<curQH->Get4Momentum()<<G4endl;
#endif
    theQHadrons.push_back(curQH);            // (delete equivalent)
  }

  theWorld              = right->theWorld;
  nBarClust             = right->nBarClust;
  f2all                 = right->f2all;
  tot4Mom               = right->tot4Mom;

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
  std::for_each(theQCandidates.begin(), theQCandidates.end(), DeleteQCandidate());
#ifdef debug
  G4cout<<"~G4QEnvironment: before theQuasmons nQ="<<theQuasmons.size()<<G4endl;
#endif
  std::for_each(theQuasmons.begin(), theQuasmons.end(), DeleteQuasmon());
#ifdef debug
  G4cout<<"~G4QEnvironment: before theQHadrons nH="<<theQHadrons.size()<<G4endl;
#endif
  std::for_each(theQHadrons.begin(), theQHadrons.end(), DeleteQHadron());
#ifdef debug
  G4cout<<"~G4QEnvironment: === DONE ==="<<G4endl;
#endif
}

G4double G4QEnvironment::SolidAngle=0.8;     // Part of Solid Angle to capture (@@A-dep.)
G4bool   G4QEnvironment::EnergyFlux=false;   // Flag for Energy Flux use (not MultyQuasmon)
G4double G4QEnvironment::PiPrThresh=141.4;   // Pion Production Threshold for gammas
G4double G4QEnvironment::M2ShiftVir=20000.;  // Shift for M2=-Q2=m_pi^2 of the virtualGamma
G4double G4QEnvironment::DiNuclMass=1880.;   // DoubleNucleon Mass for VirtualNormalization
// Fill the private static parameters
void G4QEnvironment::SetParameters(G4double solAn, G4bool efFlag, G4double piThresh,
                                   G4double mpisq, G4double dinum)
{//  ======================================================================================
  EnergyFlux=efFlag;       // Flag for Energy Flux use instead of Multy Quasmon
  SolidAngle=solAn;        // Part of Solid Angle to capture secondaries (@@A-dep)
  PiPrThresh=piThresh;     // Pion Production Threshold for gammas
  M2ShiftVir=mpisq;        // Shift for M2=-Q2=m_pi^2 of the virtual gamma
  DiNuclMass=dinum;        // Double Nucleon Mass for virtual normalization
}

const G4QEnvironment& G4QEnvironment::operator=(const G4QEnvironment &right)
{// ========================================================================
  // theQHadrons (Vector)
  G4int nQH             = right.theQHadrons.size();
  if(nQH) for(G4int ih=0; ih<nQH; ih++)
  {
    G4QHadron* curQH    = new G4QHadron(right.theQHadrons[ih]);
#ifdef fdebug
    G4cout<<"G4QE::CopyOper=:cH#"<<ih<<","<<curQH->GetQC()<<curQH->Get4Momentum()<<G4endl;
#endif
    theQHadrons.push_back(curQH);            // (delete equivalent)
  }

  theWorld              = right.theWorld;
  nBarClust             = right.nBarClust;
  f2all                 = right.f2all;

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
    theQCandidates.push_back(curQC);        // (delete equivalent)
  }

  theEnvironment        = right.theEnvironment;
  return *this;
}

// Member function for Quasmon Creation & Environment nucleus modification
void G4QEnvironment::CreateQuasmon(const G4QContent& projQC, const G4LorentzVector& proj4M)
{//========================================================================================
  //static const G4double mNeut= G4QPDGCode(2112).GetMass();
  //static const G4double mProt= G4QPDGCode(2212).GetMass();
  //static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mPi2 = mPi*mPi;
  //static const G4QContent gamQC(0,0,0,0,0,0);
  //static const G4QContent pimQC(1,0,0,0,1,0);
  //static const G4QContent pipQC(0,1,0,1,0,0);
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4QNucleus vacuum(90000000);
  G4QContent valQ(0,0,0,0,0,0);             // Prototype of the Quasmon's Quark Content
  G4LorentzVector q4Mom(0.,0.,0.,0.);       // Prototype of the Quasmon's 4-momentum
  nBarClust = 1;                            // By default only quasi-free nucleons
  G4double  projE=proj4M.e();               // energy of the projectile
  if(projE<0.)
  {
    G4cout<<"***G4QEnvironment::CreateQuasmon: projE="<<projE<<"<=0, QC="<<projQC<<G4endl;
    //throw G4QException("G4QEnvironment::CreateQuasmon: Energy of the Projectile <0.");
  }
  G4double  projM2=proj4M.m2();             // projectile's squared mass (print & v.gamma)
  G4int     targPDG=theEnvironment.GetPDG();// PDG Code of the target nucleus
  if(targPDG>80000000&&targPDG!=90000000)   // Interaction with a nuclear target
  {
    G4double  tgMass=theEnvironment.GetMass();// mass of the target (QEnvironment) nucleus
#ifdef pdebug
    G4cout<<"G4QEnvironment::CreateQ:Interact "<<projQC<<proj4M<<"(m2="
          <<projM2<<") + A="<<targPDG<<",M="<<tgMass<<G4endl;
#endif
    G4int envZ=theEnvironment.GetZ();       // A#of protons in the nucleus
    G4int envN=theEnvironment.GetN();       // A#of neutrons in the nucleus
    G4int envS=theEnvironment.GetS();       // A#of lambdas in the nucleus
    G4int nP  =theWorld->GetQPEntries();    // A#of initialized particles in CHIPS World
    G4int nCl =nP-80;                       // A#of initialized clusters in CHIPS World
    if(nCl<0) G4cerr<<"---Worning---G4QEnv::CreaQ:nP="<<nP<<" for Targ="<<targPDG<<G4endl;
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
#ifdef pdebug
	G4cout<<"G4QE::CrQ:TNuc:Z="<<envZ<<",N="<<envN<<",S="<<envS<<",nC="<<nBarClust<<G4endl;
#endif
    G4int projPDG=projQC.GetSPDGCode();     // Minimum hadron for the projectile QC
    G4bool pbpt=projE<PiPrThresh+(M2ShiftVir+projM2)/DiNuclMass;// PhotonBelowPionThreshold
    G4bool din=false;
    G4bool piF=false;
    G4bool gaF=false;
    //if(abs(projM2-mPi2)<.00001&&projE-mPi<0.1&&projPDG==-211) din=true;//PiAtRest case
    if(abs(projM2-mPi2)<.00001&&projE-mPi<0.1&&projPDG==-211) piF=true;// PiMiAtRest case
    //if(pbpt&&projPDG==22) din=true; // InCaseOf GammaBelowPiThresh needs DiNucl (?)
    if(pbpt&&projPDG==22) gaF=true; // InCaseOf GammaBelowPiThresh needs DiNucl (?)
    theEnvironment.SetMaxClust(nBarClust);
    nBarClust=theEnvironment.UpdateClusters(din); // Cluster Probabilities upto maxClust
#ifdef pdebug
	G4cout<<"G4QEnv::CreateQ: Nucleus("<<targPDG<<") is created ("<<nBarClust<<" clast's)";
    for(G4int ic=0;ic<nBarClust;ic++)
      G4cout<<" #"<<ic<<"("<<theEnvironment.GetProbability(ic)<<")";
    G4cout<<G4endl;
#endif
    theEnvironment.PrepareCandidates(theQCandidates,piF,gaF,proj4M);//Calc.Clust's probab's
#ifdef pdebug
	G4cout<<"G4QEnvironment::CreateQ: Cluster probab is calculated."<<G4endl;
#endif
    G4bool efFlag=false;                    // EnergyFlowFlag=FALSE (@@=DEFOLT=@@ make par)
    // ********** Change if necessary to compare Energy Flux & Multy Quasmon **************
    G4int efCounter=0;                      // Counter of Energy Flux particles
    G4QContent EnFlQC(0,0,0,0,0,0);         // Quark Content of Energy Flux
    G4LorentzVector ef4Mom(0.,0.,0.,0.);    // Summed 4-momentum of Energy Flux
    G4double proj3M=proj4M.rho();
    //   ---   Pbar     ---    Nbar  -------  LAMBDAbar  ----  SIGMA-bar  ---  SIGMA0bar
    if((projPDG==-2212||projPDG==-2112||projPDG==-3122||projPDG==-3112||projPDG==-3212||
        projPDG==-3222) && proj3M<10.)      // Only for AtRest reactions(move to Interface)
	//  ^^^^^  SIGMA+bar
	{
      // @@ Annihilation on one baryon is implemented (no annihilation on clusters! @@?) @@
#ifdef pdebug
      G4cout<<"G4QE::CreQ:Annihilation on a perif. nucleon, Z="<<envZ<<",N="<<envN<<G4endl;
#endif
      G4double   zpn=envZ+envN;             // a#of nucleons in the nucleus
      G4double   rnd=(zpn+envS)*G4UniformRand(); // Random number to find a baryon
      G4int      targBPDG = 0;              // Bary-Prototype of PDG of Periferal Target
      G4int      targNPDG = 90000000;       // Nucl-Prototype of PDG of Periferal Target
      G4QContent targQC(0,0,0,0,0,0);       // Quark Content of Periferal Target
      if     (rnd<envN)                     // Neutron is a Periferal Target
      {
        targBPDG = 2112;
        targNPDG = 90000001;
        targQC   = neutQC;
	  }
      else if(rnd<zpn)                      // Proton is a Periferal Target
      {
        targBPDG = 2212;
        targNPDG = 90001000;
        targQC   = protQC;
	  }
      else                                  // Lambda is a Periferal Target
      {
        targBPDG = 3122;
        targNPDG = 91000000;
        targQC   = lambQC;
	  }
      theEnvironment.Reduce(targNPDG);      // Subtract periferal baryon from Nucleus
#ifdef pdebug
      G4cout<<"G4QEnvironment::CQ:"<<targNPDG<<" is selected Env="<<theEnvironment<<G4endl;
#endif
      G4double resMass=theEnvironment.GetGSMass(); // Nuclear mass after baryon subtraction
      G4double barMass=tgMass-resMass;      // Mass of the bound baryon for annihilation
      tgMass=resMass;                       // New mass of theEnvironment
      q4Mom=G4LorentzVector(0,0,0,barMass)+proj4M;// 4-mom of the intermediate B-Bbar Quasm
      valQ=targQC+projQC;                   // Quark Content of intermediate B-Bbar Quasmon
      G4Quasmon* pan = new G4Quasmon(valQ,q4Mom);// N-Nbar Quasm creation (del.at 9th line)
      G4QNucleus vE(90000000);               // The annihilation in vacuum (in NuclMatter?)
#ifdef pdebug
      G4cout<<"G4QE::CreQ: before Fragment vE="<<vE<<",QQC="<<valQ<<",Q4M="<<q4Mom<<G4endl;
#endif
      G4QHadronVector* output=pan->Fragment(vE,1);//Output of inVacAnnihilation*DESTROY*<-+
#ifdef pdebug
	  G4cout<<"G4QEnvironment::CreateQ: Before delet pan."<<G4endl;  //                   ^
#endif
      delete pan;                              // The N-NbarQuasmon is deleted A.S.A.P.   ^
      G4QHadronVector input;                   // Input for MultyQuasmon **DESTROY**<---+ ^
      //G4int trgPDG = theEnvironment.GetPDG();// New PDG Code for the Residual Nucleus ^ ^
      G4LorentzVector trg4M(0.,0.,0.,resMass); // New 4-momentum for the ResidualNucleus^ ^
      G4int tNH = output->size();              // For the selection LOOP                ^ ^
      G4ThreeVector dir = RndmDir();           // For the selection in LOOP (@@ at rest)^ ^
#ifdef pdebug
	  G4cout<<"G4QEnvir::CreateQ: Loop over "<<tNH<<" annihilation hadrons"<<G4endl;//  ^ ^
#endif
      for (G4int ind=0; ind<tNH; ind++)        // Loop over annihilation  QHadrons      ^ ^
      {
        G4QHadron*   curHadr = output->operator[](ind); // Pointer to theCurrentHadron  ^ ^
        G4int           shDFL= curHadr->GetNFragments();// A#of decFragments for proj.  ^ ^
        G4LorentzVector sh4m = curHadr->Get4Momentum(); // 4Mom for the projectile      ^ ^
        G4ThreeVector   shDIR= sh4m.vect().unit(); // unitVector in the hadron mom. dir.^ ^
#ifdef pdebug
		G4cout<<"G4QE::CrQ:"<<ind<<","<<shDFL<<",PDG="<<curHadr->GetPDGCode()<<G4endl;//^ ^
#endif
        if(!shDFL)                               // Final (not decayed) QHadron (d==0)  ^ ^
		{
#ifdef pdebug
		  G4cout<<"G4QE::CQ:f="<<efFlag<<","<<dir.dot(shDIR)<<","<<SolidAngle<<G4endl;//^ ^
#endif
		  if(dir.dot(shDIR)>SolidAngle)          // Sum Up Hadrons and make Energy Flow ^ ^
		  {
            if(efFlag)                           // => Case of "Energy Flux approach"   ^ ^
		    {
              G4QContent shQC = curHadr->GetQC();// QuarkContent of the Current Hadron  ^ ^
              ef4Mom+=sh4m;
              EnFlQC+=shQC;
              efCounter++;
#ifdef pdebug
              G4int hPDG=curHadr->GetPDGCode();    // Only for gebug printing           ^ ^
              G4LorentzVector h4M = curHadr->Get4Momentum();  // Only for gebug printing^ ^
			  G4cout<<"G4QE::CrQ:#"<<efCounter<<", PDG="<<hPDG<<", h4M="<<h4M<<G4endl;//^ ^
#endif
		    }
		    else                                   //=>"MultyQuasFragmentation"(!efFlag)^ ^
		    {
              G4QHadron* mqHadron = new G4QHadron(curHadr);
              input.push_back(mqHadron);           // Fill hadron-copy (del equiv)  <...^ ^
#ifdef pdebug
              G4int hPDG=curHadr->GetPDGCode();    // Only for debug printing           ^ ^
              G4LorentzVector h4M = curHadr->Get4Momentum(); // Only for gebug printing ^ ^
		      G4cout<<"G4QE::CrQ: Fill #"<<ind<<", PDG="<<hPDG<<", h4M="<<h4M<<G4endl;//^ ^
#endif
		    }                                      //                                   ^ ^
		  }                                        //                                   ^ ^
		  else                                     // DirectFilling of the output vector^ ^
		  {                                        //                                   ^ ^
#ifdef fdebug
            G4int hPDG=curHadr->GetPDGCode();      // Only for gebug printing           ^ ^
            G4LorentzVector h4M = curHadr->Get4Momentum(); // Only for gebug printing   ^ ^
		    G4cout<<"G4QE::CrQ: Fill OUT #"<<ind<<",PDG="<<hPDG<<",h4M="<<h4M<<G4endl;//^ ^
#endif
            G4QHadron* curHadron = new G4QHadron(curHadr); //                           ^ ^
            theQHadrons.push_back(curHadron);      // TheQHadrs are filled as new Hadrs ^ ^
          }
		} // End of the LOOP over projectiles                                           ^ ^
	  } // End of LOOP over "output" of annihilation                                    ^ ^
      std::for_each(output->begin(), output->end(), DeleteQHadron());//DESTRUCT output>-^-^
      output->clear();                             //                                   ^ ^
      delete output;                               // ==================================^=^
      if(!efFlag)                                  // =>NotEnergyFlux=MultyQuasmon Case ^
	  {
        if(!(input.size()))                        // *RETURN* Without Quasmon creation-^
		{                                          //                                   ^
#ifdef pdebug
	      G4cout<<"*G4QEnv::CreateQ: Peripheral annihilation stack is empty"<<G4endl; //^
#endif
		  return;                                  // Do not clear and delete objects --^
        }                                          //                                   ^
#ifdef pdebug
	    G4cout<<"G4QEnvironment::CreateQ: fakeQuasmon & restoring parameters"<<G4endl;// ^
#endif
        G4Quasmon fakeQ;                   // fake Quasmon to get and restore parameters ^
        G4double QTemper=fakeQ.GetTemper();// Temperature defined by user for Quasmons   ^
        G4double QSOverU=fakeQ.GetSOverU();// S/U defined by user for Quasmons           ^
        G4double QEtaSup=fakeQ.GetEtaSup();// Eta Suppresion defined by user in Quasmons ^
        G4Quasmon::SetParameters(180.,.1,.3);// Fixed Parameters for N-barN Annihilation ^
        // From this point the new temporary environment is created (recursive)          ^
        G4QEnvironment* muq = new G4QEnvironment(input,theEnvironment.GetPDG());//---+   ^
#ifdef pdebug
	    G4cout<<"G4QEnvironment::CreateQ: before input.clearAndDestroy()"<<G4endl; //^   ^
#endif
        std::for_each(input.begin(), input.end(), DeleteQHadron());//DESTROING inp >-^---^
        input.clear();                             // ===============================^===^
        theEnvironment = muq->GetEnvironment();    // RestoreResidEnv after interact.^  
        G4QuasmonVector* outQ = muq->GetQuasmons();// Copy of quasmons *DESTROY* <---^---+
        G4QHadronVector* outH = muq->GetQHadrons();// Copy of hadrons *DESTROY* <----^-+ ^
        G4int noh = outH->size();                  // a#oh hadrons in TmpEnv UpToNow ^ ^ ^
        if(noh) for(G4int kh=0; kh<noh; kh++)      // One can escape it but...       ^ ^ ^
        {
          G4QHadron* curH = new G4QHadron(outH->operator[](kh));//CopyOfDestroiedTMP ^ ^ ^
#ifdef fdebug
          G4cout<<"G4QE::CreateQ:#"<<kh<<","<<curH->GetQC()<<curH->Get4Momentum()<<G4endl;
#endif
          theQHadrons.push_back(curH);             // Fill new Hadrons in theQHadrons^ ^ ^
		}
        std::for_each(outH->begin(), outH->end(), DeleteQHadron());// >--------------^=^ ^
        outH->clear();                             //                                ^ ^ ^
        delete outH;                               // >==============================^=^ ^
#ifdef pdebug
	    G4cout<<"G4QEnvironment::CreateQ: before delete muq"<<G4endl; //             ^   ^
#endif
        delete muq;                                //================================^   ^
        G4Quasmon::SetParameters(QTemper,QSOverU,QEtaSup); // Recover QuasmonParameters  ^
	    G4int nMQ = outQ->size();                  // A#of Quasmons in MultyQuasmonOutput^
#ifdef pdebug
	    G4cout<<"G4QEnvironment::CreateQ: after GetQuasmon nMQ="<<nMQ<<G4endl; //        ^
#endif
        if(nMQ) for(G4int mh=0; mh<nMQ; mh++)      // Can escape creation/distruct but...^
        {
          G4Quasmon* curQ = new G4Quasmon(outQ->operator[](mh));// Copy to destroy TMP(?)^
          theQuasmons.push_back(curQ);             // Fill Quasmon-copies in theQuasmons ^
	    }
        std::for_each(outQ->begin(), outQ->end(), DeleteQuasmon()); // >-----------------^
        outQ->clear();                             //                                    ^
        delete outQ;                               // >==================================^
#ifdef pdebug
	    G4cout<<"G4QEnvironment::CreateQ: befor return"<<G4endl;
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
#ifdef pdebug
	G4cout<<"G4QEnvironment::CreateQ: Interaction Probabilities are calculated"<<G4endl;
#endif
    G4int nCandid = theQCandidates.size();
    G4double maxP = theQCandidates[nCandid-1]->GetIntegProbability();
    G4int i=0;
    G4QContent    curQC;                           // Quark Content of the selected cluster
    if(nCandid<=0)
	{
	  G4cerr<<"---Worning---G4QEnv::CreaQ:nC="<<nCandid<<",maxP="<<maxP<<",E="
            <<theEnvironment<<G4endl;
      //throw G4QException("G4QEnvironment::CreateQ: Can not select a cluster");
	}
    if(nCandid==1||maxP==0.)
	{
#ifdef pdebug
	  G4cout<<"***G4QEnv::CrQ:MaxP=0||nCand=1: Use all Env., Env="<<theEnvironment<<G4endl;
#endif
      curQC=theEnvironment.GetQCZNS();
      theEnvironment=vacuum;
    }
    else
    {
      G4double totP = maxP * G4UniformRand();
#ifdef pdebug
	  G4cout<<"G4QEnvironment::CrQ:nC="<<nCandid<<", maxP="<<maxP<<", totP="<<totP<<G4endl;
#endif
      while(theQCandidates[i]->GetIntegProbability()<totP) i++;
      G4QCandidate* curCand = theQCandidates[i];// Pointer to selected cluster to interact
      curQC   = curCand->GetQC();               // Get QuarkContent of the selected cluster
      G4QNucleus targClust(curQC.GetP(),curQC.GetN(),curQC.GetL());//Define Clust as a QNuc
      G4double clMass=targClust.GetGSMass();    // Mass of residual nuclear environment
      G4LorentzVector pq4M=proj4M+G4LorentzVector(0.,0.,0.,clMass); 
      if(pq4M.m()>=clMass)
	  {
#ifdef pdebug
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
#ifdef pdebug
	      G4cout<<"***G4QEnv::CrQ: Deep virtual, use all Env,Env="<<theEnvironment<<G4endl;
#endif
          curQC=theEnvironment.GetQCZNS();
          theEnvironment=vacuum;
        }
        else
		{
          G4QHadron* projH = new G4QHadron(projQC,proj4M);
          theQHadrons.push_back(projH);
	      G4cerr<<"---Worning---G4QE::CrQ:Fill Proj asItIs QC/4m="<<projQC<<proj4M<<G4endl;
          return;
        }
	  }
	}
    G4double envMass=theEnvironment.GetGSMass();   // Mass of residual nuclear environment
    if(projPDG==22&&projE<PiPrThresh+(M2ShiftVir+projM2)/DiNuclMass) // Gam+quark Interact
    //if(2>3)                                      //@@ TMP:PhotoAbsorbtion by q is closed
	{
      q4Mom=G4LorentzVector(0.,0.,0.,tgMass-envMass);// PhotoInteracts with BoundedCluster
      valQ=curQC;
#ifdef pdebug
      G4cout<<"G4QE::CrQ:Q="<<q4Mom<<valQ<<"+vg="<<proj4M<<",Env="<<theEnvironment<<G4endl;
#endif
      G4Quasmon* curQuasmon = new G4Quasmon(valQ, q4Mom, proj4M);//Interaction gam+q inside
      theQuasmons.push_back(curQuasmon);  // Insert Quasmon without incid. gamma (del.eq.)
	}
    else if(abs(projM2-mPi2)<.00001&&projE-mPi<0.1&&projPDG==-211)
    //if(2>3)                                //@@ ***TMP*** PionAbsorbAtRest by q is closed
	{
      q4Mom=proj4M+G4LorentzVector(0.,0.,0.,tgMass-envMass);// PION + BoundCluster
      valQ=EnFlQC+curQC;
      if(projE<mPi)G4cout<<"*InputError*G4QE::CrQ:Ener(pi-)="<<projE<<"<mPi="<<mPi<<G4endl;
#ifdef pdebug
      G4cout<<"G4QEnv::CrQ:Q="<<q4Mom<<valQ<<"+pi="<<proj4M<<",E="<<theEnvironment<<G4endl;
#endif
      G4Quasmon* curQuasmon = new G4Quasmon(valQ, q4Mom, -proj4M);//Interact gam+q inside
      theQuasmons.push_back(curQuasmon);  // Insert Quasmon without incid. gamma (del.eq.)
	}
    else
	{
      q4Mom=proj4M+G4LorentzVector(0.,0.,0.,tgMass-envMass);//Projectile + BoundCluster
      valQ=EnFlQC+curQC;
#ifdef pdebug
      G4cout<<"G4QEnv::CreQAll: Q="<<q4Mom<<valQ<<", QEnv="<<theEnvironment<<G4endl;
#endif
      G4Quasmon* curQuasmon = new G4Quasmon(valQ, q4Mom);
      theQuasmons.push_back(curQuasmon); // Insert Quasmon (even hadron/gamma) (del.eq.)
	}
  }
  else
  {
    G4cerr<<"---Worning---G4QEnvironment::CreateQuasmon:Strange targPDG="<<targPDG<<G4endl;
    //throw G4QException("***G4QEnvironment::CreateQuasmon: Impossible target");
  }
}

// Calculate a probability to interact with clusters for the givven PDG of the projectile
void G4QEnvironment::PrepareInteractionProbabilities(const G4QContent& projQC, G4double AP)
//   ===========================================================================Proj.3Mom.=
{
  G4double sum    = 0.;                            // Sum of probabilities of interaction
  G4double probab = 0.;                            // Interaction probability
  G4double denseB = 0.;                            // A#of*prob baryons in dense part
  G4double allB   = 0.;                            // A#of*prob baryons in the nucleus
  if(2>3) allB=AP;                                 // Trick to use not used AP
  //in PR//G4int pPDG = projQC.GetSPDGCode();        // PDG code of the projectile particle
  for (unsigned index=0; index<theQCandidates.size(); index++)
  {
    G4QCandidate* curCand=theQCandidates[index];   // Intermediate pointer
    G4int cPDG  = curCand->GetPDGCode();
    if(cPDG>80000000&&cPDG!=90000000)              // ===> Nuclear cluster case
	{
      G4QNucleus cN(cPDG);
      G4int zc = cN.GetZ();                        // "Z" of the cluster
      G4int nc = cN.GetN();                        // "N" of the cluster
      ///////////G4int sc = cN.GetS();                         // "S" of the cluster
      G4int ac = cN.GetA();                        // "A" of the cluster
      G4double nOfCl=curCand->GetPreProbability(); // A number of clusters of the type
      G4double dOfCl=curCand->GetDenseProbability();// A number of clusters in dense region
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
      //else if((pPDG==-211&&AP<10.||pPDG==22&&AP<150.)&&ac<2)probab=0.;//PiAtRest/GamUndPi
      //else if(pPDG==22&&AP<152.&&ac<2) probab=nOfCl*ac*fact*.5;//GamUnderPi(QuarkCapture)
      ////////////////////////else if((pPDG==-211&&AP<10.)&&ac<2) probab=0;//PiCapAtRest(D)
      //else if(pPDG==-211&&AP<10.)            probab=nOfCl*fact;// special PiCaptureAtRest
      //else if(pPDG==-211&&AP<10.)            probab=nOfCl*ac*(ac-1)*fact;
      else                                   probab=nOfCl*ac*fact;
      //else                                      probab=dOfCl*ac*fact;
      //if(ac>1) probab=0.;                       // Suppress clusters
      //if(ac>2) probab=0.;                       // Suppress heavy clusters
#ifdef sdebug
      G4int pPDG      = projQC.GetSPDGCode();      // PDG code of the projectile particle
      G4int rPDG = qQC.GetSPDGCode();
      G4double baryn = qQC.GetBaryonNumber();
      G4double charge= qQC.GetCharge();
      G4double dq= abs(baryn-charge-charge);
	  G4cout<<"G4QE::PrepIntProb:C="<<cPDG<<",P="<<probab<<",ac="<<ac<<",dq="<<dq<<",qC="
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
//   ==============================================================
{
#ifdef pdebug
  G4cout<<"G4QEnvironment::InitClustersVector called with nC="<<maxClust<<G4endl;
#endif
  if(maxClust>=0) for (G4int i=0; i<maxClust; i++) 
  {
    G4int clustQCode = i+80; // Q-code of the cluster in the CHIPS World "IsoNuclei"
    G4QPDGCode clustQPDG;
#ifdef sdebug
	G4cout<<"G4QEnvironment::InitClustersVector: Before Init Q ="<<clustQCode<<G4endl;
#endif
    clustQPDG.InitByQCode(clustQCode);
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

// Fragmentation of the QEnvironment with MultyQuasmon (the main member function)
G4QHadronVector  G4QEnvironment::HadronizeQEnvironment()
//               ================***********************
{
  static const G4int  NUCPDG = 90000000;
  static const G4QNucleus vacuum(NUCPDG);
  static const G4QContent PiQC(0,1,0,1,0,0);
  static const G4QContent K0QC(1,0,0,0,0,1);
  static const G4QContent KpQC(0,1,0,0,0,1);
  //static const G4QContent alQC(6,6,0,0,0,0);
  static const G4QPDGCode nQPDG(2112);
  static const G4QPDGCode pQPDG(2212);
  static const G4QPDGCode lQPDG(3122);
  static const G4double mPi0 = G4QPDGCode(111).GetMass();
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mK   = G4QPDGCode(321).GetMass();
  static const G4double mK0  = G4QPDGCode(311).GetMass();
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mSigM= G4QPDGCode(3112).GetMass();
  static const G4double mSigP= G4QPDGCode(3222).GetMass();
  //static const G4double mAlph = G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double eps=.003;
  G4int nQuasmons = theQuasmons.size();
#ifdef pdebug
  G4cout<<"G4QE::HQE:*HADRONIZE Q-ENVIRONMENT="<<theEnvironment<<",nQ="<<nQuasmons<<G4endl;
#endif
  if(nQuasmons<1)                            // "No Quasmons" case -> Fill QEnviron
  {
    G4int nPDG = theEnvironment.GetPDG();    // PDG code of the residual Nucl.Environ.
#ifdef pdebug
	G4cout<<"G4QE::HQE:***NO QUASMONS***Env="<<nPDG<<theEnvironment.Get4Momentum()<<G4endl;
#endif
    if(nPDG==90000000) return theQHadrons;
    if(nPDG>80000000)
	{
      G4QHadron* rNucleus = new G4QHadron(theEnvironment); // Create HadronEnvironment
      theQHadrons.push_back(rNucleus);       // Fill GS - no further decay (del. equiv.)
#ifdef fdebug
	  G4cout<<"G4QEnv::HadrQE: >>>> Fill Environment="<<theEnvironment<<G4endl;
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
      G4int ist=theQuasmons[iq]->GetStatus();// Status of the Quasmon before fragmentation
      if(ist)
	  {
        G4QHadronVector* output=theQuasmons[iq]->Fragment(vE,1);//!!!DESTROY!!! <---------+
        G4int ast=theQuasmons[iq]->GetStatus();  // Quasmon's Status after fragmentation  ^
        if(!ast) nlq--;                          // Reduce nlq is Quasmon decayed         ^
        G4int nHadrons = output->size();         // A#of output Hadrons in the Quasmon    ^
#ifdef pdebug
        G4cout<<"G4QEnv::HadrQE: ***Vacuum*** Q#"<<iq<<", nHadr="<<nHadrons<<G4endl; //   ^
#endif
        if(nHadrons>0)                           // Copy QHadrons-Quasmon to Output       ^
	    {
    	  for (G4int ih=0; ih<nHadrons; ih++)    // LOOP over QHadrons of the Quasmon     ^
          {
            G4QHadron* curH = new G4QHadron(output->operator[](ih));// (Del 7 lines below)^
#ifdef fdebug
            G4cout<<"G4QEnv::HadrQE:Vacuum, H#"<<ih<<", QPDG="<<curH->GetQPDG() //        ^
                  <<",4M="<<curH->Get4Momentum()<<G4endl; //                              ^
#endif
            theQHadrons.push_back(curH);         // Fill hadron-copy (delete equivalent)  ^
          }
	    }
        else                                     // => "Quasmon can't decay" case         ^
	    {
          G4QContent totQC=theQuasmons[iq]->GetQC();//                                    ^
          G4int     tQBN=totQC.GetBaryonNumber();// Baryon Number of not decayed Quasmon  ^
          G4QNucleus     tqN(totQC);             // Define the quasmon as a nucleus       ^
          G4double   gsM=tqN.GetMZNS();          // GS Mass                               ^
          G4LorentzVector tot4M=theQuasmons[iq]->Get4Momentum();
          G4double totQM=tot4M.m();              // Real Mass of Quasmon                  ^
          if(tQBN>0&&totQM>gsM)                  // => "Try Quasmon evaporation" case     ^
		  {
            G4QHadron* nuclQ = new G4QHadron(totQC,tot4M); //                             ^
            EvaporateResidual(nuclQ);            // Evaporate ResNuc (del.equiv)          ^
            theQuasmons[iq]->KillQuasmon();      // Kill evaporated Quasmon               ^
            nlq--;                               //                                       ^
		  }
          else if(iq+1<nQuasmons&&nlq>1)         // => "Try to merge with next" case      ^
		  {
            G4int s=theQuasmons[iq+1]->GetStatus();//Status of the next Quasmon           ^
            theQuasmons[iq+1]->IncreaseBy(theQuasmons[iq]);// Merge with the next Quasmon ^
            theQuasmons[iq]->KillQuasmon();      // Kill the week Quasmon                 ^
            if(s) nlq--;                         // Reduce a number of "living Quasmons"  ^
		  }
          else if(iq+1==nQuasmons&&iq&&nlq>1)    // => "Quasmon stack is exhosted" case   ^
		  {
            G4int s=theQuasmons[0]->GetStatus(); // Status of the first Quasmon           ^
            theQuasmons[0]->IncreaseBy(theQuasmons[iq]);// Merge with the first Quasmon   ^
            theQuasmons[iq]->KillQuasmon();      // Kill the week Quasmon                 ^
            if(s) nlq--;                         // Reduce a number of "living Quasmons"  ^
		  }
		  else                                   // "Have a chance to recover" case       ^
		  {                                      //                                       ^
#ifdef pdebug
		    G4cout<<"***G4QE::HQE:"<<iq<<",n="<<nHadrons<<",Tot="<<totQC<<totQM<<G4endl;//^
		    for (G4int kq=0; kq<nQuasmons; kq++) // LOOP over Quasmons for DEBUG PRINTING ^
		      G4cout<<kq<<",St/QC="<<theQuasmons[kq]->GetStatus()<<theQuasmons[kq]->GetQC()
                    <<",M="<<theQuasmons[kq]->Get4Momentum().m()<<G4endl; //              ^
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
                G4int ast=theQuasmons[iq]->GetStatus();  // Status of the Quasmon       ^ ^
                if(!ast) nlq--;                  // Reduce nlq if Quasmon decayed       ^ ^
                G4int nHadrons=curout->size();   // A#of outputQHadrons in theDecayedQ  ^ ^
#ifdef pdebug
                G4cout<<"G4QEnv::HadrQE:VacuumRecoverQ#"<<iq<<",n="<<nHadrons<<G4endl;//^ ^
#endif
                if(nHadrons>0)                   // => "QHadrons from Quasmon to Output"^ ^
				{                                //                                     ^ ^
    	          for (G4int ih=0; ih<nHadrons; ih++) // LOOP over Hadrons of theQuasmon^ ^
				  {                              //                                     ^ ^
                    G4QHadron* curH = new G4QHadron(curout->operator[](ih)); //         ^ ^
#ifdef fdebug
                    G4cout<<"G4QEnv::HadrQE:Recovered, H#"<<ih<<", QPDG="<<curH->GetQPDG()
                          <<",4M="<<curH->Get4Momentum()<<G4endl;    //                 ^ ^
#endif
                    totQC-=curH->GetQC();        // totQC recalculation                 ^ ^
					tot4M-=curH->Get4Momentum(); // tot4M recalculation                 ^ ^
                    theQHadrons.push_back(curH); // Fill hadron-copy (delete equivalent)^ ^
                    delete curout->operator[](ih); //>-*Necessary to delete instances*>-^ ^
                  } // End of LOOP over Hadrons of the Quasmon                          ^ ^
                  curout->clear();               //                                     ^ ^
                  delete curout;                 //>*Necessary to delete VectPointers*>=^ ^
                  corrf = false;                 // Corrected: go out of the while loop ^ ^
                  //break;                       // @@ ??                               ^ ^
	            } // End of check for existing output Hadrons in the Quasmon            ^ ^
                else                             //                                     ^ ^
				{                                //                                     ^ ^
                  std::for_each(curout->begin(), curout->end(), DeleteQHadron()); // >--^ ^
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
              G4int totBN=totQC.GetBaryonNumber();// Total Baryon Number                  ^
              G4int totPDG=totQC.GetZNSPDGCode();// Convert QC to PDGCOde for the nucleus ^
              if(totS) totPDG-=totS*999999;      // @@ ??                                 ^
#ifdef pdebug
		      G4cout<<"G4QE::HQE: totPDG="<<totPDG<<",totM="<<totQM<<G4endl; //           ^
#endif
              G4QHadron* evH = new G4QHadron(totQC,tot4M); // Create a Hadron-ResidualNucl^
              if(totS<0&&totBN>0) DecayAntiStrange(evH);// Decay anti-strange nucleus (K+)^
              else if(totBN==2) DecayDibaryon(evH); // Decay dibaryon (delete equivalent) ^
              else EvaporateResidual(evH);      // Evaporate ResNuc (del.equiv)           ^
              std::for_each(output->begin(), output->end(), DeleteQHadron());// >---------^
              output->clear();                   //                                       ^
              delete output;                     // >=====================================^
              CleanUp();                         //                                       ^
              return theQHadrons;                //                                       ^
            }                                    //                                       ^
		    else if(!nOfOUT)                     // Still remain not used Quasmons        ^
			{                                    //                                       ^
		      G4cerr<<"***G4QE::HQ:T="<<tot4M<<totQC<<",M="<<totQM<<"< gsM=" //           ^
                    <<gsM<<",d="<<dM<<",Env="<<theEnvironment<<G4endl; //                 ^
              throw G4QException("G4QEnvironment::HadronizeQEnv: Can't decay Quasmon");// ^
	        }                                    //                                       ^
		  } // End of PANIC treatment                                                     ^
		} // End of trouble handling with Quasmon decay in Vacuum                         ^
        std::for_each(output->begin(), output->end(), DeleteQHadron());  // >-------------^
        output->clear();                         //                                       ^
        delete output;                           // >=====================================^
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
    G4int nsHadr  = theQHadrons.size();        // Update the value of OUTPUT entries
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
    //G4int   c3Max = 27;
    //G4int   c3Max = 3;
    G4int   c3Max = 1;
    //G4int   premC = 27;
    //G4int   premC = 3;
    G4int   premC = 1;
    G4int   envA=theEnvironment.GetA();
    //if(envA>1&&envA<13) premC = 24/envA;
    if(envA>1&&envA<19) premC = 36/envA;
    //if(envA>1&&envA<25) premC = 48/envA;
    //if(envA>1&&envA<31) premC = 60/envA;
    G4int  sumstat= 2;                     // Sum of statuses of all Quasmons
    G4bool force  = false;                 // Prototype of the Force Major Flag
    G4int cbR     =0;                      // Counter of the "Stoped by Coulomb Barrier"
    G4int cbRM    =3;                      // MaxCounter of the "Stoped by Coulomb Barrier"
    G4int totC    = 0;                     // Counter to break the "infinit" loop
    G4int totCM   = 227;                   // Limit for this counter
    //G4int totCM   = 27;                    // Limit for this counter
    while (sumstat||totC<totCM)            // ===***=== The MAIN "FOREVER" LOOP ===***===
	{
      totC++;
      if(nQuasmons==1&&sumstat==3) cbR++;
      else cbR=0;
      G4QContent envQC=theEnvironment.GetQCZNS();// QuarkCont of current NuclearEnvironment
      G4QContent totQC=envQC;                    // Total QuarkContent in the system
      G4double   envM =theEnvironment.GetMass(); // mass of NuclEnvironment (@@GetMZNS())
      G4double   sumM =envM;                     // Sum of all residual masses in theSystem
      G4LorentzVector env4M=theEnvironment.Get4Momentum();      
      G4LorentzVector tot4M=env4M;               // 4-momentum of the Total System
      sumstat         =0;
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
#ifdef pdebug
	    G4cout<<"G4QEnv::HadrQE:#"<<iq<<", Qst="<<Qst<<", Q="<<Q4M<<Q4M.m()<<QQC<<", Env="
              <<theEnvironment<<G4endl;
#endif
        if(Qst==1||Qst==3||Qst==4) fCount++;     // Incr. counterSuccessfulFragmentations
        if(Qst>0)                  eCount++;     // Incr. counter of existing Quasmons 
	  } // End of summation LOOP over Quasmons
      G4int      totS  =totQC.GetStrangeness();  // Total Strangeness of the Total System
      G4int      totBN =totQC.GetBaryonNumber(); // Total Baryon Number of the Total System
      G4int      totPDG=0;                       // Total PDG Code for the Current compound
      G4double   totM  =0.;                      // min(GroundSt)Mass of theResidualSystem
      if(totBN<2)  
	  {
        totPDG=totQC.GetSPDGCode();              // Min totPDGCode for theCurrentCompound
        if(totPDG) totM=G4QPDGCode(totPDG).GetMass(); // min Mass of the Residual System
        else throw G4QException("G4QEnv::HadrQEnv: Impossible PDG for B=1");
      }
      else
	  {
        G4QNucleus totN(totQC,tot4M);            // Excited nucleus for theResidualSystem
        totM=totN.GetMZNS();                     // min(GroundSt)Mass of theResidualSystem
        totPDG=totN.GetPDG();                    // Total PDG Code for the Current compound
	  }
#ifdef pdebug
	  G4cout<<"G4QEnv::HadrQE:totC="<<totC<<"<totCM="<<totCM<<",ss="<<sumstat<<G4endl;
#endif
      if(totC>=totCM||cbR>cbRM)
      {
        CleanUp();
        G4QHadron* evH = new G4QHadron(totQC,tot4M);// Create a Hadron for ResidualNucl
        if(totS<0&&totBN>0)DecayAntiStrange(evH);
        else if(totBN==2) DecayDibaryon(evH); // Try to decay unstable Dibaryon
        else EvaporateResidual(evH);       // Try to evaporate residual (delete equivalent)
        return theQHadrons;
	  }
      // === Now we should be prepared for evaporation ===
      G4int      totChg=totQC.GetCharge();    // Total Electric Charge of the Total System
#ifdef pdebug
      if(totPDG==90999999||totPDG==90999000||totPDG==90000999||totPDG==89999001)
		 G4cout<<"***G4QEnv::HadrQEnv: Meson (1) PDG="<<totPDG<<", M="<<tot4M.m()<<G4endl;
      G4int           nOH=theQHadrons.size(); // A#of output hadrons
      G4LorentzVector s4M=tot4M;              // Total 4-momentum (@@ only for checking)
      if(nOH) for(G4int ih=0; ih<nOH; ih++) s4M+=theQHadrons[ih]->Get4Momentum();     
	  G4cout<<"G4QEnv::HadrQE:tBN="<<totBN<<",s="<<sumstat<<",fC="<<fCount<<",eC="<<eCount
            <<",En="<<theEnvironment<<",nH="<<nOH<<",tLV="<<s4M<<totQC<<G4endl;
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
          G4cerr<<"***G4QEnv::HadrQE:totPDG="<<totPDG<<", totQC="<<totQC<<G4endl;
          throw G4QException("G4QEnvironment::HadronizeQEnvir: ImpossibleHadron in CHIPS");
		}
	  }
      G4double totMass = tot4M.m();           // Total effective Mass
      G4bool   Premium = eCount&&premC&&envM; // Premium condition
      G4int    count3  =0;
      if(sumstat&&(fCount||Premium)&&!force&&count3<c3Max)//=>"Still try to decay Quasmons"
	  {
        if(!fCount)premC--;                   // Reduce premium efforts counter
	    if(nQuasmons) for (G4int jq=0; jq<nQuasmons; jq++)//FragmentationLOOP over Quasmons
	    {
	      G4Quasmon* pQ     = theQuasmons[jq];// Pointer to the CurrentQuasmon <--<--<--+
          G4int      status = pQ->GetStatus();// Old status of the Quasmon              ^
          if(status)                          // Skip dead Quasmons                     ^
          {                                   //                                        ^
 	        G4QHadronVector* output=pQ->Fragment(theEnvironment,eCount);//*!DEST* <--<--^-+
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
	        for (G4int js=0; js<nQuasmons; js++)// Subtract 4mom's of Quasmons (compare)^ ^
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
            {
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
            if(nHadrons) for(G4int kpo=0; kpo<nHadrons; kpo++)//LOOP over out QHadrons  ^ ^
            {                                 //                                        ^ ^
	  		  G4QHadron* insH =output->operator[](kpo);// Pointer to theOutputQHadron   ^ ^
              G4int qhsNF  = insH->GetNFragments(); // A#of secondary fragments         ^ ^
              if(!qhsNF)                      // Subtract only final hadrons            ^ ^
			  {                               //                                        ^ ^
                G4LorentzVector qhs4Mom = insH->Get4Momentum();// 4M of theOutputQHadron^ ^
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
              if(nHadrons>0)                  // Transfer QHadrons from Quasm to Output ^ ^
              {                               //                                        ^ ^
    	        for (G4int ih=0; ih<nHadrons; ih++) // LOOP over Q-output QHadrons      ^ ^
                {                             //                                        ^ ^
				  G4QHadron* inpH =output->operator[](ih); //                           ^ ^
                  G4int hC=inpH->GetCharge(); // Charge of the Hadron                   ^ ^
                  G4int hF=inpH->GetNFragments();// Number of fragments                 ^ ^
                  G4double hCB=0.;            // Coulomb Barrier                        ^ ^
                  G4double hKE=0.;            // Kinetic Energy of the Hadron           ^ ^
                  G4LorentzVector hLV=inpH->Get4Momentum(); //                          ^ ^
#ifdef pdebug
                  G4cout<<"G4QEnv::HadrQE: hC="<<hC<<",hF="<<hF<<",4M="<<hLV<<G4endl;// ^ ^
#endif
                  G4bool can=hC&&!hF;         // Charged and not yet decayed hadron     ^ ^
                  if(can)                     //                                        ^ ^
                  {                           //                                        ^ ^
                    G4int hB=inpH->GetBaryonNumber(); //                                ^ ^
                    hCB=theEnvironment.CoulombBarrier(hC,hB); //                        ^ ^
                    hKE=hLV.e()-hLV.m();      //                                        ^ ^
				  }                           //                                        ^ ^
                  if(can&&hKE<hCB)            // => "Suck Hadron in Quasm or Env" case  ^ ^
				  {                           //                                        ^ ^
                    if(status)                // => "Suck in the existing Quasmon" case ^ ^
                    {                         //                                        ^ ^
                      G4QContent tQC=inpH->GetQC()+pQ->GetQC();  //                     ^ ^
                      G4LorentzVector tLV=hLV+pQ->Get4Momentum();//                     ^ ^
                      pQ->InitQuasmon(tQC,tLV); // Reinitialize the current Quasmon     ^ ^
#ifdef pdebug
	                  G4cout<<"G4QE::HQE:Medium, H#"<<ih<<", QPDG="<<inpH->GetQPDG() // ^ ^
                            <<",4M="<<inpH->Get4Momentum()<<" is suckedInQ"<<G4endl; // ^ ^
#endif
				    }                         //                                        ^ ^
                    else                      // => "Suck in the existing Quasmon" case ^ ^
                    {                         //                                        ^ ^
                      G4QContent tQC=inpH->GetQC()+theEnvironment.GetQCZNS();//         ^ ^
                      G4LorentzVector tLV=hLV+theEnvironment.Get4Momentum(); //         ^ ^
                      theEnvironment=G4QNucleus(tQC,tLV); // Reinit currentEnvironment  ^ ^
#ifdef pdebug
	                  G4cout<<"G4QE::HQE:Med,H#"<<ih<<",PDG="<<inpH->GetQPDG()<<",4M="//^ ^
                            <<inpH->Get4Momentum()<<" is suckedInEnvironment"<<G4endl;//^ ^
#endif
				    }                         // End of the STATUS IF                   ^ ^
				  }                           // End of the E/Q suck Hdron IF           ^ ^
                  else if(!hF)                // => "Hadron can go out" case            ^ ^
                  {                           //                                        ^ ^
                    G4QHadron* curH = new G4QHadron(inpH); //                           ^ ^
#ifdef fdebug
                    G4LorentzVector ph4M=curH->Get4Momentum(); // 4-mom of the hadron   ^ ^
                    G4double phX=ph4M.x();    // p_x of the hadron                      ^ ^
                    G4double phY=ph4M.y();    // p_y of the hadron                      ^ ^
                    G4double phZ=ph4M.z();    // p_x of the hadron                      ^ ^
                    G4double phCost=phZ/sqrt(phX*phX+phY*phY+phZ*phZ);// Hadr cos(theta)^ ^
	                G4cout<<"G4QEnv::HadrQE:Medium, H#"<<ih<<",QPDG="<<curH->GetQPDG()//^ ^
                          <<", 4M="<<ph4M<<", ct="<<phCost<<G4endl; //                  ^ ^
#endif
                    theQHadrons.push_back(curH); // Fill hadron-copy (delete equivalent)^ ^
                    totQC-=curH->GetQC();     //                                        ^ ^
                    tot4M-=curH->Get4Momentum(); //                                     ^ ^
				  }                           //                                        ^ ^
				}                             //                                        ^ ^
                pQ->ClearOutput();            // Hadrons are filled, Clear Frag-out <-<-^ ^
                count3=0;                     // Reset counter of empty hadronizations    ^
	          }                               //                                          ^
              else count3++;                  // Increment counter of empty hadronizations^
			}                                 //                                          ^
            else if(status<0||status==2)      // => "PANIC or NOTHING was done" case      ^
            {                                 //                                          ^
              if(eCount==1 && status<0 && CheckGroundState(pQ,true))// Correct & Finish   ^
              {                               //                                          ^
                std::for_each(output->begin(), output->end(), DeleteQHadron()); // ->->->-^
                output->clear();              //                                          ^
                delete output;                // >========================================^
                pQ->KillQuasmon();            // If BackFusion succeeded, kill the Quasmon^
                eCount--;                     // Reduce the number of the living Quasmons ^
                return theQHadrons;           //                                          ^
              }                               //                                          ^
              else if(status<0&&nHadrons)     // This is just a confusion in the status...^
              {                               //                                          ^
		        G4cerr<<"***G4QEnv::HadrQE: nH="<<nHadrons<<"< status="<<status<<G4endl;//^
                std::for_each(output->begin(), output->end(), DeleteQHadron()); // ->->->-^
                output->clear();              //                                          ^
                delete output;                // >========================================^
                throw G4QException("G4QEnvironment::HadronizeQEnvironment:DoNothPANIC");//^
			  }                               //                                          ^
              else if(status==2)              // Treat PANIC for stat=2 (NothingWasDone)  ^
              {                               //                                          ^
                G4QContent qQC=pQ->GetQC();   // QuarkContent of the Quasmon              ^
                G4int      pqC=qQC.GetCharge(); // Charge (nP) of the Current Quasmon     ^
                G4int      pqS=qQC.GetStrangeness(); // Strangeness (nL) of theCurrQuasmon^
                G4int      pqB=qQC.GetBaryonNumber(); // BaryNumber of the CurrentQuasmon ^
                G4LorentzVector cq4M=pQ->Get4Momentum(); // 4Mom of the Current Quasmon   ^
                G4double cqMass=cq4M.m();     // Real Mass of the current Quasmon         ^
                G4double fqMass=G4QPDGCode(22).GetNuclMass(pqC,pqB-pqC-pqS,pqS);//CQ FreeM^
#ifdef edebug
    		    G4cout<<"G4QEnv::HQE:M="<<cqMass<<">fM="<<fqMass<<",S="<<pqS<<",C="<<pqC//^
                      <<", envPDG="<<theEnvironment.GetPDG()<<",qQC="<<qQC<<G4endl; //    ^
#endif
                if(pqB>0&&pqS<0&&cqMass>fqMass)// "AntiStrangeNucleus-Chipolino" case     ^
				{                             //                                          ^
                  G4QHadron* nuclQ = new G4QHadron(qQC,cq4M);// Hadron for AntiStrangeNuc.^
                  DecayAntiStrange(nuclQ);    // Decay the AntyStrangeNucl (Del.Equiv.)   ^
                  pQ->KillQuasmon();          // If BackFusion succeeded, kill theQuasmon ^
#ifdef edebug
     		      G4cout<<"G4QEnv::HQE:Status after kill (#"<<jq<<")="<<pQ->GetStatus()// ^
                        <<", nH="<<theQHadrons.size()<<G4endl;                         // ^
#endif
                  tot4M-=cq4M;                // Update TotalResidNucleus for hadronizPro.^
                  totQC-=qQC;                 // Update total QC for the HadronizationPro.^
                  eCount--;                   // Reduce the number of the living Quasmons ^
                }                             //                                          ^
                else if(theEnvironment.GetPDG()!=NUCPDG) // "Nuclear Environment" case    ^
				{                             //                                          ^
                  G4LorentzVector t4M=cq4M+theEnvironment.Get4Momentum(); // Q+E tot4Mom  ^
                  G4double      tM=t4M.m();   // Real total (Quasmon+Environment) mass    ^
                  G4QContent envQC=theEnvironment.GetQCZNS(); // QuarkCont of NucEnviron  ^
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
		              {
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
			            }                     //                                          ^
                        if(h1M+h2M+envM<totMass) // => "Three parts decay" case           ^
                        {                     //                                          ^
                          G4LorentzVector h14M(0.,0.,0.,h1M);           //                ^
                          G4LorentzVector h24M(0.,0.,0.,h2M);           //                ^
                          G4LorentzVector e4M(0.,0.,0.,envM);           //                ^
                          if(!G4QHadron(tot4M).DecayIn3(h14M,h24M,e4M)) //                ^
                          {                   //                                          ^
				            G4cerr<<"***G4QE::HQE:(0)tM="<<tot4M.m()<<"->h1="<<h1QPDG //  ^
                                  <<"("<<h1M<<")+h2="<<h1QPDG<<"("<<h2M<<")+envM="<<envM//^
                                  <<"=="<<h1M+h2M+envM<<G4endl; //                        ^
				            throw G4QException("G4QE::HadrQEnv:QChip+E DecIn3 error"); // ^
			              }                   //                                          ^
                          G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M); //    ^
                          theQHadrons.push_back(h1H);        // (delete equivalent)       ^
#ifdef edebug
						  G4cout<<"G4QE::HQE:(1) H1="<<h1QPDG<<h14M<<G4endl;        //    ^
#endif
                          G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M); //    ^
                          theQHadrons.push_back(h2H);        // (delete equivalent)       ^
#ifdef edebug
                          G4cout<<"G4QE::HQE:(1) H2="<<h2QPDG<<h24M<<G4endl;        //    ^
#endif
                          G4QHadron* qeH = new G4QHadron(envPDG,e4M);               //    ^
                          theQHadrons.push_back(qeH);        // (delete equivalent)       ^
#ifdef edebug
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
                            std::for_each(output->begin(),output->end(),DeleteQHadron());
                            output->clear();  // -->-->-->-->-->-->-->-->-->-->-->-->-->--^
                            delete output;    // >========================================^
                            return theQHadrons; //                                        ^
                          }                   //                                          ^
                          std::for_each(output->begin(),output->end(),DeleteQHadron());//>^
                          output->clear();    // -->-->-->-->-->-->-->-->-->-->-->-->-->--^
                          delete output;      // >========================================^
#ifdef edebug
                          G4cerr<<"--Worning--G4QE::HQE:tM="<<tot4M.m()<<"< h1="<<h1QPDG//^
                                <<"(M="<<h1M<<")+h2="<<h1QPDG<<"(M="<<h2M<<")+EM="<<envM//^
                                <<"="<<h1M+h2M+envM<<G4endl; //                           ^
			              //throw G4QException("G4QEnv::HQE:(0)Chi+Env mass > totMass");//^
#endif
                          CleanUp();               //                                     ^
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
                          std::for_each(output->begin(),output->end(),DeleteQHadron());//>^
                          output->clear();    // -->-->-->-->-->-->-->-->-->-->-->-->-->--^
                          delete output;      // >========================================^
				          throw G4QException("G4QEnv::HadrQEnv: Q+Env DecIn2 error"); //  ^
			            }                     //                                          ^
                        G4QHadron* qH = new G4QHadron(qPDG,fq4M);// the out going Quasmon ^
                        theQHadrons.push_back(qH); // (delete equivalent)                 ^
#ifdef edebug
						  G4cout<<"G4QE::HQE:QuasmH="<<qPDG<<fq4M<<G4endl;          //    ^
#endif
                        G4QHadron* qeH = new G4QHadron(envPDG,qe4M);//theRecoilEnvironment^
#ifdef edebug
						  G4cout<<"G4QE::HQE:EnvironH="<<envPDG<<qe4M<<G4endl;      //    ^
#endif
                        if(envPDG==92000000||envPDG==90002000||envPDG==90000002) //       ^
						                                           DecayDibaryon(qeH); // ^
                        else theQHadrons.push_back(qeH);// (del.equiv.) *** this too ***  ^
		              }                       //                                          ^
                      std::for_each(output->begin(),output->end(),DeleteQHadron());//->->-^
                      output->clear();        //                                          ^
                      delete output;          // >========================================^
                      CleanUp();              // Clean up Environ and Quasmon             ^
                      return theQHadrons;     // Finish the hadronization process         ^
                    }                         //                                          ^
                    else status=-1;           // Q+E && totM below MassShell - PANIC      ^
                  }                           //                                          ^
                  if(eCount==1 && tM>=curM)   //==>for one Quasmon evaporate ResTotN      ^
                  {                           //                                          ^
                    theEnvironment.InitByPDG(NUCPDG); // Cancele the Environment          ^
#ifdef edebug
                    G4cout<<"G4QEnv::HadQEnv: Before evaporation t4M="<<tot4M<<G4endl; // ^
#endif
                    std::for_each(output->begin(),output->end(),DeleteQHadron());//->->->-^
                    output->clear();          //                                          ^
                    delete output;            // >========================================^
                    CleanUp();                // Clean up the Environ and Quasmons        ^
                    G4QHadron* evH = new G4QHadron(totQC,tot4M);// Hadron for ResidNucl   ^
                    EvaporateResidual(evH);   // Evaporate residual (del.equiv.)          ^
                    return theQHadrons;       //                                          ^
                  }                           //                                          ^
                  else if(eCount==1 && CheckGroundState(pQ,true)) // Correct and Finish   ^
                  {                           //                                          ^
                    std::for_each(output->begin(), output->end(), DeleteQHadron()); // ->-^
                    output->clear();          //                                          ^
                    delete output;            // >========================================^
                    pQ->KillQuasmon();        // If BackFusion succeeded, kill Quasm      ^
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
                  // @@ There still can be a case for p,Sigma+ or n,Sigma- (PDGCode?)     ^
                  else if(PDGQ==3112||PDGQ==3222) // "Sigma+/Sigma-" case                 ^
                  {                           //                                          ^
#ifdef edebug
				    G4cout<<"G4QEnv::HadrQEnv:Sigma dm="<<cqMass-fpMass<<G4endl; //       ^
#endif
                    G4double hyM=mNeut;       // Prototype of the hyperon mass            ^
                    G4int    hyPDG=2112;      // Prototype of the hyperon PDG Code        ^
                    G4double pigM=mPi;        // Prototype of the gamma/pion mass         ^
                    G4int    pigPDG=-211;     // Prototype of the gamma/pion PDG Code     ^
                    if(PDGQ==3112)            // "Sigma-" case                            ^
					{ //                                                                  ^
                      if(cqMass>mPi+mLamb)    // "Lambda + Pii is possible" case          ^
					  { //                                                                ^
                        hyM   = mLamb;        // Lambda mass                              ^
                        hyPDG = 3122;         // Lambda PDG Code                          ^
					  } //                                                                ^
                      else if(cqMass>mSigM)   // "Sigma- gamma decay" case                ^
					  { //                                                                ^
                        hyM=mSigM;            // Sigma- mass                              ^
                        hyPDG=PDGQ;           // Sigma- PDG Code                          ^
                        pigM=0.;              // Gamma mass                               ^
                        pigPDG=22;            // Gamma PDG Code                           ^
					  } //                                                                ^
                    } //                                                                  ^
                    else if(PDGQ==3222)       // "Sigma+" case                            ^
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
                        hyPDG=PDGQ;           // Sigma+ PDG Code                          ^
                        pigM=0.;              // Gamma mass                               ^
                        pigPDG=22;            // Gamma PDG Code                           ^
					  } //                                                                ^
                      if(cqMass>mPi0+mProt&&G4UniformRand()>.5) // --- "P + Pi0" case     ^
					  { //                                                                ^
                        hyM   = mProt;        // Proton mass                              ^
                        hyPDG = 2212;         // Proton PDG Code                          ^
                        pigM  = mPi0;         // Pi0 mass                                 ^
                        pigPDG= 111;          // Pi0 PDG Code                             ^
					  } //                                                                ^   
                    } //                                                                  ^
                    G4LorentzVector b4Mom(0.,0.,0.,hyM); // Hyperon mass                  ^
                    G4LorentzVector m4Mom(0.,0.,0.,pigM);// pion/gamma mass               ^
                    if(!G4QHadron(cq4M).DecayIn2(b4Mom, m4Mom)) // "DecayIn2 failed" case ^
					{ //                                                                  ^
                      G4cerr<<"---Worning---G4QE::HQE:H="<<hyPDG<<"(m="<<hyM<<")+G/Pi=" //^
							<<pigPDG<<"(m="<<pigM<<")="<<hyM+pigM<<">"<<cqMass<<G4endl; //^
                      G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC             ^
                      CleanUp();              //                                          ^
		              if(!CheckGroundState(quasH,true)) // Last posibility to correct     ^
					  { //                                                                ^
                        G4QHadron* hadr = new G4QHadron(totQC,tot4M); //                  ^
                        theQHadrons.push_back(hadr);    // Cor or fill as It Is           ^
#ifdef edebug
                        G4cerr<<"-Worn-G4QE::HQE:Sig,QC="<<totQC<<",4M="<<tot4M<<G4endl;//^
#endif
    	                //throw G4QException("G4QEnvironment::HadronizeQEnv:Sig error");//^
                      } //                                                                ^
                      delete quasH;          // Delete the temporary fake Quasmon         ^
                      return theQHadrons;    //                                           ^
                    } //                                                                  ^
#ifdef edebug
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
                  else if (PDGQ!=10)          // @@ Chipolino can wait @@                 ^
                  {                           //                                          ^
                    G4double qM =cq4M.m();    // Real mass of the Quasmon                 ^
                    G4double gsM=QPDGQ.GetMass(); // GSmass of the Quasmon                ^
#ifdef edebug
    		        G4cout<<"G4QEnv::HadrQEnv:#"<<jq<<",qM="<<qM<<">gsM="<<gsM<<G4endl;// ^
#endif
					if(abs(qM-gsM)<0.0001)    // "Fill & Kill" Case                       ^
                    {                         //                                          ^
                      G4QHadron* resQ = new G4QHadron(PDGQ,cq4M); // GSM hadron for CurQ  ^
#ifdef edebug
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
                      std::for_each(output->begin(), output->end(), DeleteQHadron());//->-^
                      output->clear();        //                                          ^
                      delete output;          // >========================================^
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
			  {                   
                //if(eCount==1 && DecayInEnvQ(pQ)) //                                     ^
                //{                             //                                        ^
                //  std::for_each(output->begin(), output->end(), DeleteQHadron());//->---^
                //  output->clear();            //                                        ^
                //  delete output;              // >======================================^
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
                    std::for_each(output->begin(), output->end(), DeleteQHadron()); //-->-^
                    output->clear();            //                                        ^
                    delete output;              // >======================================^
                    pQ->KillQuasmon();          //                                        ^
                    eCount--;                   // Reduce a#of the living Quasmons        ^
					return theQHadrons;         //                                        ^
                  }                             //                                        ^
#ifdef edebug
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
                  CleanUp();                    // Clean up THIS Quasmon                  ^
                  G4QHadron* evH = new G4QHadron(totQC,tot4M); // Create ResidNuclHadron  ^
                  EvaporateResidual(evH);       // Try to evaporate residual (del.equiv.) ^
                  std::for_each(output->begin(), output->end(), DeleteQHadron());//->-->--^
                  output->clear();              //                                        ^
                  delete output;                // >======================================^
                  force=true;                   // Make the force decision                ^
                  break;                        // Out of the fragmentation loop >->+     ^
			    }                               //                                  |     ^
			  }                                 //                                  |     ^
			}                                   //                                  |     ^
            std::for_each(output->begin(), output->end(), DeleteQHadron()); // >----|-----^
            output->clear();                    //                                  |     ^
            delete output;                      // >================================|=====^
		  }                                     //                                  |
	    } // End of fragmentation LOOP over Quasmons (jq) <--------<----------<-----+
      }
      // @@@@ For 1 livingQuasmon & totMass>minRQ+mREnv -> decay in RE+RQ (in RQ direction)
      else if(totMass>totM+.001)                // ==> "Try Evaporadecayin2tion or decay" case
	  {
#ifdef edebug
		G4cout<<"G4QEnv::HadrQE: M="<<totMass<<",PDG="<<totPDG<<",B="<<totBN<<",GSM="<<totM
              <<",dM="<<totMass-totM<<",totQC="<<totQC<<G4endl;
#endif
        if(totBN<2)                             // ==> "Baryon/Meson residual Quasmon" case
		{
          if(totPDG==90999999||totPDG==90999000||totPDG==90000999||totPDG==89999001)//"M"ca
		  {
		    G4cerr<<"---Worning---G4QE::HQE:Meson(2) PDG="<<totPDG<<",M="<<totMass<<G4endl;
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
		      G4cerr<<"--Worning--G4QE::HQE:tM="<<totMass<<"<GSM+mPi0="<<totM+mPi0<<G4endl;
              G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
              CleanUp();
		      if(!CheckGroundState(quasH,true))
              {
                G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef edebug
                G4cerr<<"***G4QE::HQE:FillAsIs(-4),QC="<<totQC<<",4M="<<tot4M<<G4endl;
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
                G4cerr<<"---Worning---G4QEnv::HadronizeQE:B="<<bPDG<<"(m="<<mBar<<") + M="
					  <<mPDG<<"(m="<<mMes<<")="<<mBar+mMes<<" > mDel="<<totMass<<G4endl;
                G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                CleanUp();
		        if(!CheckGroundState(quasH,true))
                {
                  G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                  theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef edebug
                  G4cerr<<"***G4QE::HQE:FillAsIs(-3),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
    	          //throw G4QException("G4QEnvironment::HadronizeQEnv:Del->Bar+Mes error");
                }
                delete quasH;  
                return theQHadrons;
              }
#ifdef edebug
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
              G4cerr<<"---Worning---G4QEnv::HadronizeQE:h1="<<h1PDG<<"(m="<<h1M<<") + h2="
				    <<h2PDG<<"(m="<<h2M<<")="<<h1M+h2M<<" > mChipo="<<totMass<<G4endl;
              G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
              CleanUp();
		      if(!CheckGroundState(quasH,true))
              {
                G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef edebug
                G4cerr<<"***G4QE::HQE:FillAsIs(-2),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
    	        //throw G4QException("G4QEnvironment::HadrQEnv: Chipo->1+2 decay failed");
              }
              delete quasH;  
              return theQHadrons;
            }
#ifdef edebug
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
              G4cerr<<"---Worning---G4QEnv::HadronizeQEnv: h="<<totPDG<<"(m="<<totM
                    <<") + gamma > mTot="<<totMass<<G4endl;
              G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
              CleanUp();
		      if(!CheckGroundState(quasH,true))
              {
                G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef edebug
                G4cerr<<"***G4QE::HQE:FillAsIs(-1),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
    	        //throw G4QException("G4QEnvironment::HadronizeQEnv:Gamma Decay failed");
              }
              delete quasH;  
              return theQHadrons;
            }
#ifdef edebug
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
            if(totPDG==92000000||totPDG==90002000||totPDG==90000002) DecayDibaryon(curH);
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
              G4cerr<<"---Worning---G4QEnv::HadronizeQEnv: h="<<mbPDG<<"(m="<<mbm
                    <<") + pi(m="<<mpi<<")="<<mbm+mpi<<" > mTot="<<totMass<<G4endl;
              G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
              CleanUp();
		      if(!CheckGroundState(quasH,true))
              {
                G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef edebug
                G4cerr<<"***G4QE::HQE:FillAsIs(0),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
    	        //throw G4QException("G4QEnvironment::HadronizeQE: DecIn2 mB+nPi failed");
              }
              delete quasH;  
              return theQHadrons;
            }
#ifdef edebug
	        G4cout<<"G4QE::HQE:"<<tot4M<<"->h="<<mbPDG<<h4Mom<<"+p="<<piPDG<<g4Mom<<G4endl;
#endif
            G4QHadron* curH = new G4QHadron(mbPDG,h4Mom);
            if(totPDG==92000000||totPDG==90002000||totPDG==90000002) DecayDibaryon(curH);
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
            G4QNucleus vacuum(90000000);         //                                     ^
 	        G4QHadronVector* curout=resid->Fragment(vacuum,1);//**!!DESTROY!!** <-<-+   ^
            G4int rest = resid->GetStatus();     // New status after fragm attempt  ^   ^
            if(!rest) eCount--;                  // Dec ExistingQuasmonsCounter     ^   ^
            delete resid;                        //_________________________________^___^
            G4int nHadrons = curout->size();     // a#of Hadrons in the outHV       ^
            if(nHadrons>0)                       // Transfer QHadrons to Output     ^
	        {
    	      for (G4int ih=0; ih<nHadrons; ih++)// LOOP over output QHadrons       ^
			  {                                  //                                 ^
#ifdef edebug
				G4cout<<"G4QEnv::HadrQE:NewB<2, H#"<<ih
                      <<", QPDG="<<curout->operator[](ih)->GetQPDG()
                      <<", 4M="<<curout->operator[](ih)->Get4Momentum()<<G4endl; // ^
#endif
                theQHadrons.push_back(curout->operator[](ih)); // (delete equ.) <-<-^
			  }                                  //                                 ^
	        }                                    //                                 ^
			else                                 //                                 ^
			{
              G4cerr<<"***G4QEnv::HadrQEnv:MQ="<<tot4M.m()<<",QC="<<totQC<<G4endl;//^
			  throw G4QException("G4QEnvironment::HadronizeQEnv: Quasmon decay?");//^
			}                                    // *** Do not destroy instances ***^
            curout->clear();                     // The instances are filled above  ^
            delete curout;                       // >===============================^
            return theQHadrons;
		  }
		}
        else
		{
          G4QContent tQC    =totQC;              // Not subtracted copy for error prints
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
              G4cerr<<"***G4QEnvironment::HadronizeQEnv: totPDG=0"<<G4endl;
              throw G4QException("G4QEnv::HadrQEnv: Impossible PDG for B=1");
            }
          }
          else
	      {
            G4QNucleus totN(tQC,tot4M);          // Excited nucleus for the Residual System
            totRM=totN.GetMZNS();                // min (GS) Mass of the Residual System
            totPDG=totN.GetPDG();                // Total PDG Code for the Current compound
	      }
          if(NaK)                                // ==> "Decay in K0 or K+ + NPi" case
	      {//@@ Can (must) be moved to EvaporateResidual ??
            if(NaK==1&&!NPi)                     // ==> "One anti-strange K" case
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
                G4cerr<<"***G4QE::HadronizeQE:M="<<aKPDG<<"(m="<<MaK<<")+N="<<totPDG<<"(m="
                      <<totRM<<")="<<sum<<" > mSN="<<totMass<<",d="<<sum-totMass<<G4endl;
#endif
                G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                CleanUp();
		        if(!CheckGroundState(quasH,true))
                {
                  G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                  theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef edebug
                  G4cerr<<"***G4QEnv::HQE:FillAsItIs(1),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
                  //throw G4QException("G4QEnvironment::HadronizeQEnv:AntiS-Nuc error");
                }
                delete quasH;  
                return theQHadrons;
              }
#ifdef edebug
	          G4cout<<"G4QEnv::HadronizeQEnv: SN="<<tot4M<<" -> M="
                    <<aKPDG<<m4Mom<<" + N="<<totPDG<<n4Mom<<totQC<<G4endl;
#endif
              G4QHadron* curK = new G4QHadron(aKPDG,m4Mom);
              theQHadrons.push_back(curK);       // Fill the curK (delete equivalent)
              G4QHadron* curN = new G4QHadron(totPDG,n4Mom); // @@ Use DecayDib then Evap
              EvaporateResidual(curN);           // Try to evaporate residual (del.eq.)
		    }
            else if(NaK&&NPi)                    // ==> "Anti-strange K's + DELTA's" case
		    {
              G4LorentzVector m4Mom(0.,0.,0.,MPi);
              G4LorentzVector k4Mom(0.,0.,0.,MaK);
              G4LorentzVector n4Mom(0.,0.,0.,totRM);
              if(!G4QHadron(tot4M).DecayIn3(m4Mom, k4Mom, n4Mom))
              {
                G4cerr<<"---Worning---G4QE::HadronQE:K="<<aKPDG<<"(m="<<MaK<<")+PI="<<PiPDG
                      <<"(m="<<MPi<<"+N="<<totPDG<<"(m="<<totRM<<")>toM="<<totMass<<G4endl;
                G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                CleanUp();
		        if(!CheckGroundState(quasH,true))
                {
                  G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                  theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef edebug
                  G4cerr<<"***G4QEnv::HQE:FillAsItIs(2),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
    	          //throw G4QException("G4QEnvironment::HadronizeQE:2AntiS-Nucl(1) error");
                }
                delete quasH;  
                return theQHadrons;
              }
#ifdef pdebug
	          G4cout<<"G4QEnv::HadronizeQEnv: SN="<<tot4M<<" -> nK="<<aKPDG<<k4Mom
                    <<" + nPi="<<PiPDG<<m4Mom<<" + N="<<totPDG<<n4Mom<<G4endl;
#endif
              G4LorentzVector onePi=(1./NPi)*m4Mom;// 4-mom of one pion  
              for (G4int ip=0; ip<NPi; ip++)
			  {
                G4QHadron* curP = new G4QHadron(PiPDG,onePi);
#ifdef pdebug
			    G4cout<<"G4QEnv::HadrQEnv:Pion#"<<ip<<",H="<<PiPDG<<onePi<<G4endl;
#endif
                theQHadrons.push_back(curP);     // Fill the curM (delete equivalent)
	 		  }
              G4LorentzVector oneK=(1./NaK)*k4Mom; // 4-mom of one kaon  
              for (G4int jp=0; jp<NaK; jp++)
			  {
                G4QHadron* curP = new G4QHadron(aKPDG,oneK);
#ifdef pdebug
			    G4cout<<"G4QEnv::HadrQEnv:Kaon#"<<jp<<",H="<<PiPDG<<oneK<<G4endl;
#endif
                theQHadrons.push_back(curP);     // Fill the curM (delete equivalent)
	 		  }
              G4QHadron* curN = new G4QHadron(totPDG,n4Mom);
              EvaporateResidual(curN);           // Try to evaporate residual (del.equiv.)
		    }
		    else                                 // ==> "Two anti-strange Kaons" case
		    {
              G4int N1K = NaK/2;                 // First kaon cluster
              G4int N2K = NaK-N1K;               // Second kaon cluster
              G4double mM  = MaK/NaK;            // Mass of Pi
              G4double m1M = mM*N1K;             // Mass of the first Pi-cluster
              G4double m2M = mM*N2K;             // Mass of the second Pi-cluster
              G4LorentzVector m4Mom(0.,0.,0.,m1M);
              G4LorentzVector k4Mom(0.,0.,0.,m2M);
              G4LorentzVector n4Mom(0.,0.,0.,totRM);
              if(!G4QHadron(tot4M).DecayIn3(m4Mom, k4Mom, n4Mom))
              {
                G4cerr<<"---Worning---G4QEnv::HadronizeQEnv:N*K="<<aKPDG<<"(m="<<mM<<")+N="
                      <<totPDG<<"(m="<<totRM<<") >(?)SN="<<totMass<<G4endl;
                G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                CleanUp();
		        if(!CheckGroundState(quasH,true))
                {
                  G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                  theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef edebug
                  G4cerr<<"***G4QEnv::HQE:FillAsItIs(3),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
    	          //throw G4QException("G4QEnvironment::HadronizeQE:2antiS-Nucl(2) error");
                }
                delete quasH;  
                return theQHadrons;
              }
#ifdef pdebug
	          G4cout<<"G4QEnv::HadronizeQEnv: SN="<<tot4M<<" -> N*K="<<aKPDG
                    <<" (4M1="<<m4Mom<<" + 4M2="<<k4Mom<<") + N="<<totPDG<<n4Mom<<G4endl;
#endif
              G4LorentzVector one1=(1./N1K)*m4Mom;// 4-mom of one kaon  
              for (G4int ip=0; ip<N1K; ip++)
			  {
                G4QHadron* curP = new G4QHadron(aKPDG,one1);
#ifdef pdebug
			    G4cout<<"G4QEnv::HadrQEnv:FirstKaon#"<<ip<<",H="<<aKPDG<<one1<<G4endl;
#endif
                theQHadrons.push_back(curP);     // Fill the curP (delete equivalent)
	 		  }
              G4LorentzVector one2=(1./N2K)*k4Mom;// 4-mom of one kaon  
              for (G4int jp=0; jp<N2K; jp++)
			  {
                G4QHadron* curP = new G4QHadron(aKPDG,one2);
#ifdef pdebug
			    G4cout<<"G4QEnv::HadrQEnv:SecondKaon#"<<jp<<",H="<<aKPDG<<one1<<G4endl;
#endif
                theQHadrons.push_back(curP);     // Fill the curP (delete equivalent)
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
                G4cerr<<"---Worning---G4QEnv::HadronizeQEnv:M="<<PiPDG<<"(m="<<MPi<<")+N="
				      <<totPDG<<"(m="<<totRM<<")="<<MPi+totRM<<" > mSN="<<totMass<<G4endl;
                G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                CleanUp();
		        if(!CheckGroundState(quasH,true))
                {
                  G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                  theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef edebug
                  G4cerr<<"***G4QEnv::HQE:FillAsItIs(4),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
    	          //throw G4QException("G4QEnvironment::HadronizeQEnv:Iso-Nucleus error");
                }
                delete quasH;  
                return theQHadrons;
              }
#ifdef edebug
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
                G4cerr<<"---Worning---G4QEnv::HadronizeQE:N*Pi="<<PiPDG<<"(m="<<mM<<")+N="
                      <<totPDG<<"(m="<<totRM<<") >(?)SN="<<totMass<<G4endl;
                G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                CleanUp();
		        if(!CheckGroundState(quasH,true))
                {
                  G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                  theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef edebug
                  G4cerr<<"***G4QEnv::HQE:FillAsItIs(5),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
    	          //throw G4QException("G4QEnvironment::HadronizeQE:ManyIsoNucleus error");
                }
                delete quasH;  
                return theQHadrons;
              }
#ifdef edebug
	          G4cout<<"G4QEnv::HadronizeQEnv: SN="<<tot4M<<" -> N*PI="<<PiPDG
                    <<" (4M1="<<m4Mom<<" + 4M2="<<k4Mom<<") + N="<<totPDG<<n4Mom<<G4endl;
#endif
              G4LorentzVector one1=(1./N1Pi)*m4Mom;// 4-mom of one pion  
              for (G4int ip=0; ip<N1Pi; ip++)
			  {
                G4QHadron* curP = new G4QHadron(PiPDG,one1);
                theQHadrons.push_back(curP);     // Fill the curP (delete equivalent)
	 		  }
              G4LorentzVector one2=(1./N2Pi)*k4Mom;// 4-mom of one pion  
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
#ifdef edebug
		G4cout<<"G4QE::HadrQEnv: Try FinalEvaporation t4M="<<tot4M<<",tQC="<<totQC<<G4endl;
#endif
        CleanUp();
        G4QHadron* evH = new G4QHadron(totQC,tot4M); // Create a Hadron for theResidualNucl
        EvaporateResidual(evH);                  // Try to evaporate residual (del.equiv.)
        return theQHadrons;
	  }
      else                                       // ==> "Only GSEnvironment exists" case
      { 
        if(totPDG==90000000||abs(totMass)<0.000001)
		{
          CleanUp();
          return theQHadrons;
		}
        G4double dM=totMass-totM;
#ifdef pdebug
		G4cout<<"G4QEnv::HadrQEnv:GroundState tM-GSM="<<dM<<",GSM="<<totM<<",tPDG="<<totPDG
              <<",nQ="<<nQuasmons<<G4endl;
#endif
	    G4Quasmon*       pQ = theQuasmons[0];    // Pointer to the first Quasmon          
        G4QPDGCode    QQPDG = pQ->GetQPDG();     // QPDG of the Quasmon
        G4int          QPDG = QQPDG.GetPDGCode();
        G4QNucleus    totRN(totQC,tot4M);        // Nucleus for theTotalResidualNuclearComp
        G4int          spbRN=totRN.SplitBaryon();// PossibilityToSplit baryon from Residual
        if(dM>-0.001)
		{
#ifdef pdebug
		  G4cout<<"G4QE::HadrQE:ExcitedNucleus, dM="<<dM<<">0,tBN="<<totBN<<",nQ="<<G4endl;
#endif
          G4QHadron* evH = new G4QHadron(totQC,tot4M);// Create a Hadron for ResidualNucl
          EvaporateResidual(evH);                // Try to evaporate residual (del. equiv.)
		}
        else if(nQuasmons==1&&QPDG!=22&&QPDG!=111)// => "Decay Quasmon or Q+Environ" case
		{
          G4int envPDG = theEnvironment.GetPDG();// PDGCode of the NuclQEnvironment
#ifdef pdebug
		  G4cout<<"G4QEnv::HadrQEnv: nQ=1, QPDG=="<<QPDG<<G4endl;
#endif
          if(!QPDG)
		  {
			G4cerr<<"***G4QEnv::HadrQE: Quasmon is an unknown QHadron: PDG="<<QPDG<<G4endl;
			throw G4QException("G4QEnvironment::HadronizeQEnvironment:(2)Cann'tDecayQEnv");
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
			}
            if(h1M+h2M+envM<totMass)             // => "Three parts decay" case
			{
              G4LorentzVector h14M(0.,0.,0.,h1M);
              G4LorentzVector h24M(0.,0.,0.,h2M);
              G4LorentzVector e4M(0.,0.,0.,envM);
              if(!G4QHadron(tot4M).DecayIn3(h14M,h24M,e4M))
			  {
				G4cerr<<"Worning->G4QE::HQE:M="<<tot4M.m()<<","<<totMass<<"->"<<h1QPDG<<"("
					  <<h1M<<")+"<<h1QPDG<<"("<<h2M<<")+"<<envM<<"="<<h1M+h2M+envM<<G4endl;
                G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                CleanUp();
		        if(!CheckGroundState(quasH,true))
                {
                  G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                  theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef pdebug
                  G4cerr<<"***G4QEnv::HQE:FillAsItIs(6),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
				  //throw G4QException("G4QEnv::HadrQEnv:QChipo+Environment DecIn3 Error");
                }
                delete quasH;  
                return theQHadrons;

			  }
              G4QHadron* h1H = new G4QHadron(h1QPDG.GetPDGCode(),h14M);
              theQHadrons.push_back(h1H);        // (delete equivalent)
#ifdef pdebug
			  G4cout<<"G4QE::HQE:(2) H1="<<h1QPDG<<h14M<<G4endl;
#endif
              G4QHadron* h2H = new G4QHadron(h2QPDG.GetPDGCode(),h24M);
              theQHadrons.push_back(h2H);        // (delete equivalent)
#ifdef pdebug
              G4cout<<"G4QE::HQE:(2) H2-"<<h2QPDG<<h24M<<G4endl;
#endif
              G4QHadron* qeH = new G4QHadron(envPDG,e4M);
              theQHadrons.push_back(qeH);        // (delete equivalent)
#ifdef pdebug
              G4cout<<"G4QE::HQE:(2) QEenv="<<envPDG<<e4M<<G4endl;
#endif
			}
#ifdef pdebug
            G4cerr<<"***G4QEnv::HadQEnv:tM="<<tot4M.m()<<totQC<<"< h1="<<h1QPDG<<"(M="<<h1M
                  <<")+h2="<<h1QPDG<<"(M="<<h2M<<")+eM="<<envM<<"="<<h1M+h2M+envM<<G4endl;
			//throw G4QException("G4QEnv::HadrQEnv:QChipo+Env mass > than decaying mass");
#endif
            CleanUp();
            G4QHadron* evH = new G4QHadron(totQC,tot4M);// Create a Hadron for ResidualNucl
            EvaporateResidual(evH);            // Try to evaporate residual (del. equiv.)
            return theQHadrons;
		  }
        }
        else if(spbRN)// => "Join all quasmons to the residual compound and evaporate" case
		{
#ifdef pdebug
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
#ifdef pdebug
		  G4cout<<"***G4QEnv::HadrQE: M="<<totMass<<",dM="<<dM<<",nQ="<<nQuasmons<<G4endl;
#endif
          G4int          nOfOUT  = theQHadrons.size();
          while(nOfOUT)
          {
            G4QHadron*     theLast = theQHadrons[nOfOUT-1];
            G4LorentzVector last4M = theLast->Get4Momentum();
            G4QContent      lastQC = theLast->GetQC();
            G4int           lastS  = lastQC.GetStrangeness();
            G4int           totS   = totQC.GetStrangeness();
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
                G4double m2=newNp.GetMZNS();
                G4QNucleus  newN0(totQC-K0QC);
                G4double m3=newN0.GetMZNS();
                if (m3+mK0<m2+mK)                // => "aK0+ResA is better" case
	            {
                  m1  =mK0;
                  PDG1=311;
                  m2  =m3;
                  PDG2=newN0.GetPDG();
	            }
                if(totMass>m1+m2)                // => "can decay" case
                {
                  G4LorentzVector fq4M(0.,0.,0.,m1);
                  G4LorentzVector qe4M(0.,0.,0.,m2);
                  if(!G4QHadron(tot4M).DecayIn2(fq4M,qe4M))
			      {
                    G4cerr<<"---Worning---G4QE::HadQE:tM="<<tot4M.m()<<"->aK="<<PDG1<<"(M="
                          <<m1<<")+ResA="<<PDG2<<"(M="<<m2<<")="<<m1+m2<<G4endl;
                    G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                    CleanUp();
		            if(!CheckGroundState(quasH,true))
                    {
                      G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                      theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef pdebug
                      G4cerr<<"***G4QE::HQE:FillAsIs(7),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
		  		      //throw G4QException("G4QEnv::HadrQEnv: aK+ResA DecayIn2 error");
                    }
                    delete quasH;  
                    return theQHadrons;
			      }
                  G4QHadron* H1 = new G4QHadron(PDG1,fq4M);
#ifdef pdebug
                  G4cout<<"G4QE::HQE:Kaon(Env)="<<PDG1<<fq4M<<G4endl;
#endif
                  theQHadrons.push_back(H1);     // (delete equivalent)
                  G4QHadron* H2 = new G4QHadron(PDG2,qe4M);
#ifdef pdebug
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
                G4double m2=mK0;         
                G4int  PDG2=311;
                G4QNucleus  newNp0(totQC-KpQC-K0QC);
                G4int  PDG3=newNp0.GetPDG();
                G4double m3=newNp0.GetMZNS();    // M-K^0-K^+
                G4QNucleus  newN00(totQC-K0QC-K0QC);
                G4double m4=newN00.GetMZNS();    // M-2*K^0
                G4QNucleus  newNpp(totQC-KpQC-KpQC);
                G4double m5=newNpp.GetMZNS();    // M-2*K^+
                if (m4+mK0+mK0<m3+mK+mK0 && m4+mK0+mK0<=m5+mK+mK) //=>"2K0+ResA is theBest"
	            {
                  m1  =mK0;
                  PDG1=311;
                  m3  =m4;
                  PDG3=newN00.GetPDG();
	            }
                else if(m5+mK+mK<m3+mK+mK0 && m5+mK+mK<=m4+mK0+mK0)//=>"2Kp+ResA isTheBest"
	            {
                  m2  =mK;
                  PDG1=321;
                  m3  =m5;
                  PDG3=newNpp.GetPDG();
	            }
                if(totMass>m1+m2+m3)             // => "can decay" case
                {
                  G4LorentzVector k14M(0.,0.,0.,m1);
                  G4LorentzVector k24M(0.,0.,0.,m2);
                  G4LorentzVector ra4M(0.,0.,0.,m3);
                  if(!G4QHadron(tot4M).DecayIn3(k14M,k24M,ra4M))
			      {
                    G4cerr<<"--Worning--G4QE::HQE:tM="<<tot4M.m()<<"->aK="<<PDG1<<"(M="<<m1
                          <<")+K2="<<PDG2<<"(M="<<m2<<")+A="<<PDG3<<"(M="<<m3<<")"<<G4endl;
                    G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                    CleanUp();
		            if(!CheckGroundState(quasH,true))
                    {
                      G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                      theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef pdebug
                      G4cerr<<"***G4QE::HQE:FillAsIs(8),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
				      //throw G4QException("G4QEnv::HadrQE:2K+ResidNucleus DecIn3 Error");
                    }
                    delete quasH;  
                    return theQHadrons;
			      }
                  G4QHadron* H1 = new G4QHadron(PDG1,k14M);
                  theQHadrons.push_back(H1);     // (delete equivalent)
#ifdef pdebug
                  G4cout<<"G4QE::HQE:K1(Env)="<<PDG1<<k14M<<G4endl;
#endif
                  G4QHadron* H2 = new G4QHadron(PDG2,k24M);
                  theQHadrons.push_back(H2);     // (delete equivalent)
#ifdef pdebug
                  G4cout<<"G4QE::HQE:K2(Env)="<<PDG2<<k24M<<G4endl;
#endif
                  G4QHadron* H3 = new G4QHadron(PDG3,ra4M);
                  theQHadrons.push_back(H3);     // (delete equivalent)
#ifdef pdebug
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
                G4double m2=mPi;
                G4int  PDG2=-211;
				if(totPDG==2224)
				{
                  m1=mProt;
                  PDG1=2212;
                  m2=mPi;
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
                  m2  =h2.GetMass();             // Mass of the second hadron
				}
                if(totMass>m1+m2)
				{
                  G4LorentzVector fq4M(0.,0.,0.,m1);
                  G4LorentzVector qe4M(0.,0.,0.,m2);
                  if(!G4QHadron(tot4M).DecayIn2(fq4M,qe4M))
			      {
                    G4cerr<<"---Worning---G4QE::HaQE:tM="<<tot4M.m()<<"-> h1="<<PDG1<<"(M="
					      <<m1<<") + h2="<<PDG2<<"(M="<<m2<<")="<<m1+m2<<G4endl;
                    G4Quasmon* quasH = new G4Quasmon(totQC,tot4M); // totQC not tQC!
                    CleanUp();
		            if(!CheckGroundState(quasH,true))
                    {
                      G4QHadron* hadr = new G4QHadron(totQC,tot4M); // totQC not tQC!
                      theQHadrons.push_back(hadr);   // Cor or fill as It Is
#ifdef pdebug
                      G4cerr<<"***G4QE::HQE:FillAsIs(9),QC="<<totQC<<",4M="<<tot4M<<G4endl;
#endif
				      //throw G4QException("G4QEnv::HadrQEnv: h1+h2 DecayIn2 Error");
                    }
                    delete quasH;  
                    return theQHadrons;
			      }
                  G4QHadron* H1 = new G4QHadron(PDG1,fq4M);
                  theQHadrons.push_back(H1);     // (delete equivalent)
#ifdef pdebug
                  G4cout<<"G4QE::HQE:h1="<<PDG1<<fq4M<<G4endl;
#endif
                  G4QHadron* H2 = new G4QHadron(PDG2,qe4M);
#ifdef pdebug
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
            G4double dM=totMass-totM;
#ifdef pdebug
		    G4cout<<"G4QEnv::HadrQE: Add H="<<last4M<<lastQC<<",tM="<<tot4M<<totM<<totQC
                  <<",dM="<<dM<<", tB="<<totBN<<", tS="<<totS<<G4endl;
#endif
            if(dM>-0.001&&totPDG)
		    {
              G4QHadron* evH = new G4QHadron(totPDG,tot4M);//Create Hadron for ResidNucleus
              EvaporateResidual(evH);           // Evaporate ResNuc (del.equiv)
              break;
		    }
            nOfOUT  = theQHadrons.size();       // Update the value of OUTPUT entries
		  } // End of WHILE(nOfOUT)
          nOfOUT  = theQHadrons.size();         // Update the value of OUTPUT entries
		  if(!nOfOUT)
		  {
		    G4cerr<<"---Worning---G4QEnv::HadrQE:M="<<totMass<<"<gsM="<<totM<<",dM="<<dM
                  <<", tPDG="<<totPDG<<", t4M="<<tot4M<<G4endl;
			// throw G4QException("G4QEnvironment::HadronizeQEnv:Can't decayExhostedQEnv");
            G4QHadron* evH = new G4QHadron(totPDG,tot4M);// Create Hadron for ResidNucleus
            EvaporateResidual(evH);             // Evaporate ResidNucl (del.equiv)
	      }
		}
        else                                    // "Last decay was fatal" case @@ buggy ?MK
		{
#ifdef pdebug
		  G4cout<<"***G4QEnv::HadrQE: M="<<totMass<<",dM="<<dM<<",nQ="<<nQuasmons<<G4endl;
#endif
          G4Quasmon* quasH = new G4Quasmon(totQC,tot4M);
          CleanUp();
		  if(!CheckGroundState(quasH,true))
          {
            G4QHadron* hadr = new G4QHadron(totQC,tot4M);
#ifdef pdebug
            G4cout<<"G4QE::HQE:CheckGS failed H="<<totQC<<tot4M<<G4endl;
#endif
            theQHadrons.push_back(hadr); // Cor or fill asItIs
          }
          delete quasH;  
        }
        CleanUp();
	  }
	}
  } // End of the "Nuclear Environment" case
  return theQHadrons;
} // End of the main member function HadronizeQEnvironment

// Clean up the QEnvironment to Zero
void G4QEnvironment::CleanUp()
//   =========================
{
  static const G4QNucleus vacuum(90000000);
  theEnvironment=vacuum;
  G4int nQuasmons = theQuasmons.size();
  if (nQuasmons) for (G4int iq=0; iq<nQuasmons; iq++)theQuasmons[iq]->KillQuasmon();
} // End of CleanUp

//Evaporate Residual Nucleus
void G4QEnvironment::EvaporateResidual(G4QHadron* qH, G4bool corFlag)
{//  ================================================================
  static const G4double mHel6 = G4QPDGCode(2112).GetNuclMass(2,4,0);
  static const G4double mAlph = G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double mDeut = G4QPDGCode(2112).GetNuclMass(1,1,0);
  static const G4double mNeut = G4QPDGCode(2112).GetMass();
  static const G4double mProt = G4QPDGCode(2212).GetMass();
  static const G4double mLamb = G4QPDGCode(3122).GetMass();
  static const G4double mPi   = G4QPDGCode(211).GetMass();
  static const G4double mPi0  = G4QPDGCode(111).GetMass();
  static const G4double mK    = G4QPDGCode(321).GetMass();
  static const G4double mK0   = G4QPDGCode(311).GetMass();
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4QContent deutQC(3,3,0,0,0,0);
  static const G4QContent alphQC(6,6,0,0,0,0);
  G4int       thePDG = qH->GetPDGCode();         // Get PDG code of the Residual Nucleus
  /// @@@@@@@ *** TEMPORARY TO AVOID HYPERMUCLEI FOR GEANT4 *** @@@@@@@
  if(thePDG!=91000000 && thePDG>90999999)
  {
    G4int S=(thePDG-90000000)/1000000;
    thePDG-=S*999999;
    qH->SetQPDG(G4QPDGCode(thePDG));
  }
  /// @@@ *** ^^^ END OF TEMPORARY ^^^ *** @@@
  if(thePDG<80000000)
  {
#ifdef fdebug
    G4cout<<"G4QE::EvaporateRes: FundamentalParticle="<<thePDG<<qH->Get4Momentum()<<G4endl;
#endif
    theQHadrons.push_back(qH);// TheFundamentalParticles must be FilledAsTheyAre(del.eq)
    return;
  }
  G4QContent  theQC  = qH->GetQC();          // Quark Content of the hadron
  G4int theBN=theQC.GetBaryonNumber();       // A
  G4int theC=theQC.GetCharge();              // P
  G4int theS=theQC.GetStrangeness();         // S
  G4int theN=theBN-theC-theS;                // N
  if(!thePDG) thePDG = theQC.GetSPDGCode();  // If there is no PDG code, get it from QC
  if(thePDG==10&&theQC.GetBaryonNumber()>0) thePDG=theQC.GetZNSPDGCode();
  G4double totGSM = G4QNucleus(thePDG).GetGSMass();// TheGroundStMass of theTotalResNucleus
  G4LorentzVector q4M = qH->Get4Momentum();  // Get 4-momentum of theTotalResidNucleus
  G4double    totMass = q4M.m();             // Get theRealMass of theTotalResidNucleus
#ifdef pdebug
  G4cout<<"G4QEnvironment::EvapResid(EvaRes):===IN==> PDG="<<thePDG<<",4Mom="<<q4M<<G4endl;
  G4cout<<"G4QEnviron::EvaRes: A="<<theBN<<",Z="<<theC<<", N="<<theN<<",S="<<theS<<G4endl;
#endif
  if     (thePDG==90000000)                  // ==> "Nothing in the INPUT Hadron" KEEP IT!
  {
    delete qH;
#ifdef debug
    G4cerr<<"---Worning---G4QEnvir::EvaRes: Residual Nucleus is vacuum,4Mom="<<q4M<<G4endl;
#endif
    return;
  }
  else if(thePDG==91000000||thePDG==90001000||thePDG==90000001)
  //else if(2>3)// One can easily close this decay as it will be done later (time of calc?)
  {
    G4double gsM=mNeut;
    if(thePDG==90001000)      gsM=mProt;
	else if(thePDG==91000000) gsM=mLamb;
    if(fabs(totMass-gsM)<.001)
    {
#ifdef fdebug
	  G4cout<<"G4QE::EvaporateR:GSM="<<gsM<<", H="<<thePDG<<qH->Get4Momentum()<<G4endl;
#endif
      theQHadrons.push_back(qH); // (delete equivalent)
    }
    else if(thePDG==90000001&&totMass<gsM)
	{
       G4cerr<<"***G4QE::EvaRes:Baryon "<<thePDG<<" is belowMassShell M="<<totMass<<G4endl;
       throw G4QException("G4QEnvironment::EvaporateResidual: Baryon is below Mass Shell");
    }
    else
    {
#ifdef debug
	  G4cout<<"G4QE::ER:DB="<<thePDG<<",M="<<totMass<<">"<<gsM<<",d="<<totMass-gsM<<G4endl;
#endif
      G4LorentzVector h4Mom(0.,0.,0.,gsM); // GSMass must be updated in previous while-LOOP
      G4LorentzVector g4Mom(0.,0.,0.,0.);
      if(!G4QHadron(q4M).DecayIn2(h4Mom, g4Mom))
      {
        G4cerr<<"***G4QEnv::EvaRes:h="<<thePDG<<"(GSM="<<gsM<<")+gam>tM="<<totMass<<G4endl;
        throw G4QException("G4QEnvironment::EvaporateResidual:BaryonDecInBaryon+GamError");
      }
#ifdef debug
	  G4cout<<"G4QE::EvaRes:"<<totMass<<q4M<<"->"<<thePDG<<h4Mom<<"+g="<<g4Mom<<",n="
            <<theQHadrons.size()<<G4endl;
#endif
      G4QHadron* curH = new G4QHadron(thePDG,h4Mom);
#ifdef fdebug
    G4cout<<"G4QE::EvaporateR:Hadr="<<thePDG<<h4Mom<<G4endl;
#endif
      theQHadrons.push_back(curH);         // Fill the TotalResidualNucleus (delete equiv.)
      G4QHadron* curG = new G4QHadron(22,g4Mom);
#ifdef fdebug
    G4cout<<"G4QE::EvaporateR:Gamma="<<g4Mom<<G4endl;
#endif
      theQHadrons.push_back(curG);         // Fill the gamma (delete equivalent)
      delete qH;
    }
  }
  else if(theBN==2) DecayDibaryon(qH);           //=> "Dibaryon" case (del eq.)
  else if(theBN>0&&theS<0) DecayAntiStrange(qH); // "AntyStrange Nucleus" case (del eq.)
  else if(theBN>0&&!theS&&(theC<0||theC>theBN)) DecayIsonucleus(qH);// "UnavoidableIsonucl"
  else if(thePDG==89999003||thePDG==90002999)    //=> "ISO-dibarion" (@@Look at prev.line!)
  {
    // @@ Check that it never comes here !!
    G4int  nucPDG = 2112;
    G4double nucM = mNeut;
    G4int   piPDG = -211;
    if(thePDG==90002999)
	{
      nucPDG = 2212;
      nucM   = mProt;
      piPDG  = 211;
	}
    if(totMass>mPi+nucM+nucM)
	{
      G4LorentzVector n14M(0.,0.,0.,nucM);
      G4LorentzVector n24M(0.,0.,0.,nucM);
      G4LorentzVector pi4M(0.,0.,0.,mPi);
      if(!G4QHadron(q4M).DecayIn3(n14M,n24M,pi4M))
	  {
        G4cerr<<"***G4QEnv::HadQEnv: tM="<<totMass<<"-> 2N="<<nucPDG<<"(M="
		      <<nucM<<") + pi="<<piPDG<<"(M="<<mPi<<")"<<G4endl;
		throw G4QException("G4QEnv::EvaporateResidual:ISO-dibaryon DecayIn3 error");
	  }
      delete qH;
      G4QHadron* h1H = new G4QHadron(nucPDG,n14M);
#ifdef fdebug
      G4cout<<"G4QE::EvaporateR:Bar1="<<nucPDG<<n14M<<G4endl;
#endif
      theQHadrons.push_back(h1H);                // (delete equivalent)
      G4QHadron* h2H = new G4QHadron(nucPDG,n24M);
#ifdef fdebug
      G4cout<<"G4QE::EvaporateR:Bar2="<<nucPDG<<n24M<<G4endl;
#endif
      theQHadrons.push_back(h2H);                // (delete equivalent)
      G4QHadron* piH = new G4QHadron(piPDG,pi4M);
#ifdef fdebug
      G4cout<<"G4QE::EvaporateR:Pi="<<piPDG<<pi4M<<G4endl;
#endif
      theQHadrons.push_back(piH);                // (delete equivalent)
	}
	else
	{
      delete qH;
      G4cerr<<"***G4QEnv::EvaRes: IdPDG="<<thePDG<<", q4M="<<q4M<<", M="<<totMass
            <<" < M_2N+Pi, d="<<totMass-2*nucM-mPi<<G4endl;
      throw G4QException("G4QEnvironment::EvaporateResidual:ISO-DiBar is under MassShell");
	}
  }
  else if(abs(totMass-totGSM)<.001)   // Fill as it is or decay Be8, He5, Li5 (@@ add more)
  {

    if(thePDG==90004004) DecayAlphaAlpha(qH);    // "Alpha+Alpha Decay" case (del eq.)
    else if(thePDG==90004002) DecayAlphaDiN(qH);  // Decay alpha+2p (alpha+2n is stable)
    else if((theC==theBN||theN==theBN||theS==theBN)&&theBN>1) DecayMultyBaryon(qH); 
    else if(theBN==5)    DecayAlphaBar(qH);   // Try to decay unstable A5 system (del eq.)
    else                 theQHadrons.push_back(qH); // Fill as it is (del eq.)
  }
  else if(theBN>0&&thePDG>88000000&&thePDG<89000000)//==> 2antiK in the nucleus (!Comment!)
  {
    G4cerr<<"---Worning---G4QEnv::EvaRes:MustNotBeHere.PDG="<<thePDG<<",S="<<theS<<G4endl;
    G4int bZ=theQC.GetCharge();
    G4int bN=theBN-bZ;
    G4int k1PDG = 321;
    G4double mK1= mK;
    G4int k2PDG = 321;
    G4double mK2= mK;
    G4int  nucPDG = thePDG;
    if(bZ>=bN) nucPDG+=999000;
    else
    {
      nucPDG+=999999;
      k1PDG = 311;
      mK1= mK0;
    }
    if(bZ>bN) nucPDG+=999000;
    else
    {
      nucPDG+=999999;
      k2PDG = 311;
      mK2= mK0;
    }
    G4double nucM = G4QNucleus(nucPDG).GetGSMass();
    G4cerr<<"-Worning-G4QE::EvR:M="<<nucM<<","<<nucPDG<<",1="<<k1PDG<<",2="<<k2PDG<<G4endl;
    G4LorentzVector n4M(0.,0.,0.,nucM);
    G4LorentzVector k14M(0.,0.,0.,mK1);
    G4LorentzVector k24M(0.,0.,0.,mK2);
    if(!G4QHadron(q4M).DecayIn3(n4M,k14M,k24M))
	{
      G4cerr<<"***G4QEnv::EvaRes: tM="<<totMass<<"-> N="<<nucPDG<<"(M="<<nucM<<") + k1="
            <<k1PDG<<"(M="<<mK1<<") + k2="<<k2PDG<<"(M="<<mK2<<")"<<G4endl;
	  throw G4QException("G4QEnv::EvaporateResidual: Nucleus+2antiK DecayIn3 error");
	}
    delete qH;
    G4QHadron* k1H = new G4QHadron(k1PDG,k14M);
#ifdef fdebug
    G4cout<<"G4QE::EvaRes:k1="<<k1PDG<<k14M<<G4endl;
#endif
    theQHadrons.push_back(k1H);                // (delete equivalent)
    G4QHadron* k2H = new G4QHadron(k2PDG,k24M);
#ifdef fdebug
    G4cout<<"G4QE::EvaRes:k2="<<k2PDG<<k24M<<G4endl;
#endif
    theQHadrons.push_back(k2H);                // (delete equivalent)
    G4QHadron* nH = new G4QHadron(nucPDG,n4M);
#ifdef fdebug
    G4cout<<"G4QE::EvaRes:resNuc="<<nucPDG<<n4M<<G4endl;
#endif
    theQHadrons.push_back(nH);                 // (delete equivalent)
  }
  // From here the EVA code starts
  else if(thePDG>80000000&&thePDG!=90000000||thePDG==2112||thePDG==2212||thePDG==3122)
  { // @@ Improve for Sigma+, Sigma-, Ksi0 & Ksi- content in the Total Residual Nucleus
    if(thePDG<80000000)                        // Switch from QHadronCode to QNuclearCode
    {
      if     (thePDG==2112) thePDG=90000001;   // n
      else if(thePDG==2212) thePDG=90001000;   // p
      else if(thePDG==3122) thePDG=91000000;   // lambda
	}
    G4QNucleus qNuc(q4M,thePDG);               // Make a Nucleus for theTotalResidNucleus
    G4double GSMass =qNuc.GetGSMass();         // GSMass of the Total Residual Nucleus
    G4QContent totQC=qNuc.GetQCZNS();          // QuarkCont of theTotalResidNucleus (theQC)
    G4int    bA     =qNuc.GetA();              // A#of baryons in theTotal Residual Nucleus
    G4int    bZ     =qNuc.GetZ();              // A#of protons in theTotal Residual Nucleus
    G4int    bN     =qNuc.GetN();              // A#of neutrons in theTotal ResidualNucleus
    G4int    bS     =qNuc.GetS();              // A#of lambdas in theTotal Residual Nucleus
#ifdef debug
    if(bZ==2&&bN==5)G4cout<<"G4QE::EvR: (2,5) GSM="<<GSMass<<" > "
						  <<G4QPDGCode(2112).GetNuclMass(2,4,0)+mNeut<<G4endl;
    if(bZ==1&&bN==3)G4cout<<"G4QE::EvR: (1,3) GSM="<<GSMass<<" > "
						  <<G4QPDGCode(2112).GetNuclMass(1,2,0)+mNeut<<G4endl;
    G4double dM=totMass-GSMass;
	G4cout<<"G4QEnv::EvaRes:"<<qNuc<<",PDG="<<thePDG<<",M="<<totMass<<",dM="<<dM<<G4endl;
    ////////if(dM>7.) throw G4QException("G4QEnvironment::EvaporateResidual: CALLED");
#endif
    G4int   bsCond =qNuc.SplitBaryon();        // (Bary/Deut/Alph)SeparCond for TotResNucl
    G4bool  dbsCond=qNuc.Split2Baryons();      // (Two Baryons)SeparCond for TotResidNucl
#ifdef debug
	G4cout<<"G4QEnv::EvaRes:bs="<<bsCond<<",dbs="<<dbsCond<<G4endl;
#endif
    if(abs(totMass-GSMass)<.003&&!bsCond&&!dbsCond) theQHadrons.push_back(qH);//FillAsItIs
    else if((bA==1||!bsCond&&!dbsCond)&&totMass>GSMass+.003)//==>Fuse&DecayTech (no gamDec)
	//if(2>3)                                  // Close "Fuse&Decay Technology" ***@@@***
	{
#ifdef debug
	  G4cout<<"G4QE::EvaR:SplitBar, s="<<bsCond<<",M="<<totMass<<" > GSM="<<GSMass<<G4endl;
#endif
      G4int nOfOUT = theQHadrons.size();       // Total#of QHadrons in Vector at this point
      while(nOfOUT && corFlag)                 // Try BackFusionDecays till something is in
	  {
        G4QHadron*     theLast = theQHadrons[nOfOUT-1];
        G4int          lastBN = theLast->GetBaryonNumber();
        G4int          nFragm = theLast->GetNFragments();
        //////////////////G4int          gam    = theLast->GetPDGCode();
#ifdef debug
		G4cout<<"G4QE::EvaRes:*BackFus,BN="<<lastBN<<",nF="<<nFragm<<",n="<<nOfOUT<<G4endl;
#endif
		while(nFragm)                            // => "Delete Decayed Hadrons" case
		{
          G4QHadron* thePrev = theQHadrons[nOfOUT-2];
          nFragm = thePrev->GetNFragments();
#ifdef debug
          G4int          prevPDG = thePrev->GetPDGCode();
		  G4cout<<"G4QEnv::EvaRes:DelTheLast, nFr="<<nFragm<<", pPDG="<<prevPDG<<G4endl;
#endif
          theQHadrons.pop_back();          // the prev QHadron is excluded from OUTPUT
          delete theLast;//!!When kill,DON'T forget to del. theLastQHadron as an instance!!
          theLast = thePrev;               // Update theLastPntr(Prev instead of Last)
		  nOfOUT--;
		}
        if(nOfOUT)
		{
          if(lastBN<1&&nOfOUT>1)           // => "Skip Meson/Gams & Antibaryons" case
		  {
            G4QHadron* thePrev = theQHadrons[nOfOUT-2];//***Exchange between Last & Prev***
            theQHadrons.pop_back();        // the last QHadron is excluded from OUTPUT
            theQHadrons.pop_back();        // the prev QHadron is excluded from OUTPUT
            theQHadrons.push_back(theLast);// the Last becomes the Prev (1st part of exch.)
            theQHadrons.push_back(thePrev);// the Prev becomes the Last (2nd part of exch.)
            theLast = thePrev;             // Update theLastPointer (Prev instead of Last)
		  }
          G4LorentzVector last4M = theLast->Get4Momentum();
          G4QContent  lastQC = theLast->GetQC();
          G4double lastM  = last4M.m();    // Mass of the BackFused Fragment
          totQC+=lastQC;                   // Update (increase) the total QC
          q4M+=last4M;                     // Update (increase) the total 4-momentum
          totMass=q4M.m();                 // Calculate new real total mass
          G4int totPDG=totQC.GetSPDGCode();// The updated PDG for the TotalResidualNucleus
          if(totPDG==10&&totQC.GetBaryonNumber()>0) totPDG=totQC.GetZNSPDGCode();
          G4int totBN=totQC.GetBaryonNumber();// BaryonNumber of the Total Residual Nucleus
          G4double dM=totMass-GSMass -lastM;
#ifdef debug
		  G4cout<<"G4QE::EvaRes: tM="<<totMass<<"-LM="<<lastM<<lastQC<<"-GSM="<<GSMass<<"="
                <<dM<<G4endl;
#endif
          if(dM>-0.001)
		  {
            G4QHadron* evH = new G4QHadron(totPDG,q4M);// Create QHadron for TotalResidNucl
            if(dM<=0.)
            {
              theQHadrons.pop_back();    // lastQHadron is excluded from QHadrV asIs in TRN
              delete theLast; //When kill, DON'T forget to delete lastQHadron asAnInstance!
              if(totBN==2)DecayDibaryon(evH); // Fill dibaryon (with decay products)
              else theQHadrons.push_back(evH);// Fill TRN to HVect asIs (delete equivalent)
		    }
            else                        // Decay TotalResidualNucleus in GSM+Last and Break
		    {
              G4LorentzVector r4Mom(0.,0.,0.,GSMass);
              if(!G4QHadron(q4M).DecayIn2(last4M,r4Mom))
              {
                theQHadrons.pop_back(); // lastQHadron is excluded from QHadrV as is in TRN
                delete theLast; //When kill,DON'T forget to delete lastQHadron asAnInstance
                theQHadrons.push_back(evH);// Fill TRN to Vect as it is (delete equivalent)
#ifdef debug
                G4cout<<"***G4QE::EvaR:DecayIn L"<<lastQC<<"+RN"<<totQC<<" failed"<<G4endl;
#endif
	          }
              else
              {
                delete evH;                // Delete theHadron of the TotalResidualNucleus
                theLast->Set4Momentum(last4M);// Already exists:don't create&fill,->set4Mom
                G4QHadron* nuclH = new G4QHadron(thePDG,r4Mom);
                if(thePDG==92000000||thePDG==90002000||thePDG==90000002)
                                                                  DecayDibaryon(nuclH);//DE
                else theQHadrons.push_back(nuclH);// Fill the Residual Nucleus (del.eq.)
              }
              break;
		    }
		  }
          thePDG=totPDG;                   // Make ResidualNucleus outOf theTotResidualNucl
		  GSMass=G4QPDGCode(thePDG).GetMass();// Update the Total Residual Nucleus mass
          theQHadrons.pop_back();          // the last QHadron is excluded from OUTPUT
          delete theLast;//!!When kill,DON'T forget to delete theLastQHadron asAnInstance!!
          nOfOUT--;                        // Update the value of OUTPUT entries
		}
	  }
      if(!nOfOUT || !corFlag)
      {
        G4LorentzVector h4Mom(0.,0.,0.,GSMass);//GSMass must be updated inPreviousWhileLOOP
        G4LorentzVector g4Mom(0.,0.,0.,0.);
        if(!G4QHadron(q4M).DecayIn2(h4Mom, g4Mom))
        {
          G4cerr<<"***G4QE::EvaR: h="<<thePDG<<"(GSM="<<GSMass<<")+g>tM="<<totMass<<G4endl;
          throw G4QException("G4QEnvironment::EvaporateResidual: Decay in Gamma failed");
        }
#ifdef debug
	    G4cout<<"G4QE::EvaRes: "<<q4M<<"->totResN="<<thePDG<<h4Mom<<"+gam="<<g4Mom<<G4endl;
#endif
        G4QHadron* curH = new G4QHadron(thePDG,h4Mom);
#ifdef fdebug
        G4cout<<"G4QE::EvaRes:Fragment="<<thePDG<<h4Mom<<G4endl;
#endif
        if(thePDG==92000000||thePDG==90002000||thePDG==90000002) DecayDibaryon(curH);//(DE)
        else theQHadrons.push_back(curH);  // Fill the TotalResidualNucleus (del.equiv.)
        G4QHadron* curG = new G4QHadron(22,g4Mom);
#ifdef fdebug
        G4cout<<"G4QE::EvaRes:Gamma="<<g4Mom<<G4endl;
#endif
        theQHadrons.push_back(curG);       // Fill the gamma (delete equivalent)
	  }
      delete qH;
	}
    else if(bA==2) DecayDibaryon(qH);      // Decay the residual dibaryon (del.equivalent)
    else if(bA>0&&bS<0) DecayAntiStrange(qH); // Decay the state with antistrangeness
    else if(totMass<GSMass+.003&&(bsCond||dbsCond))//==>" M<GSM but decay is possible" case
    {
#ifdef debug
	  G4cout<<"G4QE::EvR:2B="<<dbsCond<<",c="<<bsCond<<",M="<<totMass<<"<"<<GSMass<<G4endl;
#endif
      G4double gResM  =1000000.;           // Prototype of mass of residual for a neutron
      G4int    gResPDG=0;                  // Prototype of residualPDGCode for a neutron
      if(bN==4&&bZ==2&&!bS)                // It's He6 nucleus
	  {
        G4QContent resQC=totQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        gResPDG= thePDG;                   // PDG of the Residual Nucleus
        gResM  = mHel6;                    // min mass of the Residual Nucleus
	  }
      G4double nResM  =1000000.;           // Prototype of mass of residual for a neutron
      G4int    nResPDG=0;                  // Prototype of ResidualPDGCode for a neutron
      if(bsCond==112&&bN>0&&bA>1)          // There's aNeutr in theNucl, which can be split
	  {
        G4QContent resQC=totQC-neutQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        nResPDG=resN.GetPDG();             // PDG of the Residual Nucleus
        if     (nResPDG==90000001) nResM=mNeut;
        else if(nResPDG==90001000) nResM=mProt;
        else if(nResPDG==91000000) nResM=mLamb;
        else nResM=resN.GetMZNS();         // min mass of the Residual Nucleus
	  }
      G4double pResM  =1000000.;           // Prototype of mass of residual for a proton
      G4int    pResPDG=0;                  // Prototype of PDGCode of residual for a proton
      if(bsCond==2212&&bZ>0&&bA>1)         // There's aProton in Nucl, which can be split
	  {
        G4QContent resQC=totQC-protQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        pResPDG=resN.GetPDG();             // PDG of the Residual Nucleus
        if     (pResPDG==90000001) pResM=mNeut;
        else if(pResPDG==90001000) pResM=mProt;
        else if(pResPDG==91000000) pResM=mLamb;
        else pResM  =resN.GetMZNS();       // min mass of the Residual Nucleus
	  }
      G4double lResM  =1000000.;           // Prototype of mass of residual for a Lambda
      G4int    lResPDG=0;                  // Prototype of PDGCode of residual for a Lambda
      if(bsCond==3122&&bS>0&&bA>1)         // There's aLambd in theNucl, which can be split
	  {
        G4QContent resQC=totQC-lambQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        lResPDG=resN.GetPDG();             // PDG of the Residual Nucleus
        if     (lResPDG==90000001) lResM=mNeut;
        else if(lResPDG==90001000) lResM=mProt;
        else if(lResPDG==91000000) lResM=mLamb;
        else lResM  =resN.GetMZNS();       // min mass of the Residual Nucleus
	  }
      G4double dResM  =1000000.;           // Prototype of mass of residual for a Alpha
      G4int    dResPDG=0;                  // Prototype of PDGCode of residual for a Alpha
      if(bsCond==90001001&&bN>0&&bZ>0&&bA>2)// There's aDeuter in Nucl, which canBeRadiated
	  {
        G4QContent resQC=totQC-deutQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        dResPDG=resN.GetPDG();             // PDG of the Residual Nucleus
        if     (dResPDG==90000001) dResM=mNeut;
        else if(dResPDG==90001000) dResM=mProt;
        else if(dResPDG==91000000) dResM=mLamb;
        else dResM  =resN.GetMZNS();       // minMass of the Residual Nucleus
	  }
      G4double aResM  =1000000.;           // Prototype of mass of residual for a Alpha
      G4int    aResPDG=0;                  // Prototype of PDGCode of residual for a Alpha
      if(bsCond==90002002&&bN>1&&bZ>1&&bA>4)// There's Alpha in theNucl, which can be split
	  {
        G4QContent resQC=totQC-alphQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        aResPDG=resN.GetPDG();             // PDG of the Residual Nucleus
        if     (aResPDG==90000001) aResM=mNeut;
        else if(aResPDG==90001000) aResM=mProt;
        else if(aResPDG==91000000) aResM=mLamb;
        else aResM  =resN.GetMZNS();       // minMass of the Residual Nucleus
	  }
      G4double nnResM  =1000000.;          // Prototype of mass of residual for a dineutron
      G4int    nnResPDG=0;                 // Prototype of ResidualPDGCode for a dineutron
      if(dbsCond&&bN>1&&bA>2)              // It's nucleus and there is a dineutron
	  {
        G4QContent resQC=totQC-neutQC-neutQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        nnResPDG=resN.GetPDG();            // PDG of the Residual Nucleus
        if     (nnResPDG==90000001) nnResM=mNeut;
        else if(nnResPDG==90001000) nnResM=mProt;
        else if(nnResPDG==91000000) nnResM=mLamb;
        else nnResM  =resN.GetMZNS();      // min mass of the Residual Nucleus
	  }
      G4double ppResM  =1000000.;          // Prototype of mass of residual for a diproton
      G4int    ppResPDG=0;                 // Prototype of ResidualPDGCode for a diproton
      if(dbsCond&&bZ>1&&bA>2)              // It's nucleus and there is a diproton
	  {
        G4QContent resQC=totQC-protQC-protQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        ppResPDG=resN.GetPDG();            // PDG of the Residual Nucleus
        if     (ppResPDG==90000001) ppResM=mNeut;
        else if(ppResPDG==90001000) ppResM=mProt;
        else if(ppResPDG==91000000) ppResM=mLamb;
        else ppResM  =resN.GetMZNS();      // min mass of the Residual Nucleus
	  }
      G4double npResM  =1000000.;          // Prototype of ResidualMass for proton+neutron
      G4int    npResPDG=0;                 // Prototype of ResidualPDGCode for a prot+neut
      if(dbsCond&&bN>0&&bZ>0&&bA>2)        // It's nucleus and there is aProton & aNeutron
	  {
        G4QContent resQC=totQC-neutQC-protQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        npResPDG=resN.GetPDG();            // PDG of the Residual Nucleus
        if     (npResPDG==90000001) npResM=mNeut;
        else if(npResPDG==90001000) npResM=mProt;
        else if(npResPDG==91000000) npResM=mLamb;
        else npResM  =resN.GetMZNS();      // min mass of the Residual Nucleus
	  }
      G4double lnResM  =1000000.;          // Prototype of residualMass for lambda+neutron
      G4int    lnResPDG=0;                 // Prototype of ResidualPDGCode for aLambda+Neut
      if(dbsCond&&bN>0&&bS>0&&bA>2)        // It's nucleus and there is aLambda & aNeutron
	  {
        G4QContent resQC=totQC-lambQC-protQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        lnResPDG=resN.GetPDG();            // PDG of the Residual Nucleus
        if     (lnResPDG==90000001) lnResM=mNeut;
        else if(lnResPDG==90001000) lnResM=mProt;
        else if(lnResPDG==91000000) lnResM=mLamb;
        else lnResM  =resN.GetMZNS();      // min mass of the Residual Nucleus
	  }
      G4double lpResM  =1000000.;          // Prototype of residualMass for a proton+lambda
      G4int    lpResPDG=0;                 // Prototype of ResidualPDGCode for theProt+lamb
      if(dbsCond&&bS>0&&bZ>0&&bA>2)        // It's nucleus and there is aProton and aLambda
	  {
        G4QContent resQC=totQC-neutQC-protQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        lpResPDG=resN.GetPDG();            // PDG of the Residual Nucleus
        if     (lpResPDG==90000001) lpResM=mNeut;
        else if(lpResPDG==90001000) lpResM=mProt;
        else if(lpResPDG==91000000) lpResM=mLamb;
        else lpResM  =resN.GetMZNS();      // minMass of the Residual Nucleus
	  }
      G4double llResM  =1000000.;          // Prototype of mass of residual for a di-lambda
      G4int    llResPDG=0;                 // Prototype of ResidPDGCode for the di-lambda
      if(dbsCond&&bS>1&&bA>2)              // It's nucleus and there is a di-lambda
	  {
        G4QContent resQC=totQC-neutQC-protQC;
        G4QNucleus resN(resQC);            // Pseudo nucleus for the Residual Nucleus
        llResPDG=resN.GetPDG();            // PDG of the Residual Nucleus
        if     (llResPDG==90000001) llResM=mNeut;
        else if(llResPDG==90001000) llResM=mProt;
        else if(llResPDG==91000000) llResM=mLamb;
        else llResM  =resN.GetMZNS();      // min mass of the Residual Nucleus
	  }
#ifdef debug
      G4cout<<"G4QEnv::EvaRes: rP="<<pResPDG<<",rN="<<nResPDG<<",rL="<<lResPDG<<",N="
            <<bN<<",Z="<<bZ<<",nL="<<bS<<",totM="<<totMass<<",n="<<totMass-nResM-mNeut
            <<",p="<<totMass-pResM-mProt<<",l="<<totMass-lResM-mLamb<<G4endl;
#endif
      if(   thePDG==90004004 || thePDG==90002004 && totMass>mHel6+.003
         || bA>4 && bsCond && bN>1 && bZ>1 && totMass>aResM+mAlph
         || bA>1 && bsCond && (   bN>0&&totMass>nResM+mNeut
                               || bZ>0&&totMass>pResM+mProt
							   || bS>0&&totMass>lResM+mLamb)
         ||bA>2&&(bN>0&&bZ>0&&(bsCond&&totMass>dResM+mDeut||dbsCond && totMass>dResM+mDeut)
				               ||dbsCond && (   bN>1&&totMass>nnResM+mNeut+mNeut
                                             || bZ>1&&totMass>ppResM+mProt+mProt
                                             || bS>1&&totMass>llResM+mLamb+mLamb
                                             || bN&&bS&&totMass>lnResM+mLamb+mNeut
                                             || bZ&&bS&&totMass>lpResM+mLamb+mProt)))
	  {
        G4int barPDG = 90002002;           // Just for the default case of Be8->alpha+alpha
        G4int resPDG = 90002002;
        G4int thdPDG = 0;
        G4double barM= mAlph;
        G4double resM= mAlph;
        G4double thdM= mNeut;              // This default value is used in the IF
        G4double tMC=totMass+.0002;
		if(gResPDG&&tMC>mHel6+.003)        // Can make radiative decay of He6 (priority 0)
		{
          barPDG=90002004;
          resPDG=22;
          barM  =mHel6;
          resM  =0.;
		}
		else if(nResPDG&&tMC>nResM+mNeut)  // Can radiate a neutron (priority 1)
		{
          barPDG=90000001;
          resPDG=nResPDG;
          barM  =mNeut;
          resM  =nResM;
		}
		else if(pResPDG&&totMass+.001>pResM+mProt)   // Can radiate a proton (priority 2)
		{
          barPDG=90001000;
          resPDG=pResPDG;
          barM  =mProt;
          resM  =pResM;
		}
		else if(lResPDG&&tMC>lResM+mLamb)  // Can radiate a Lambda (priority 3)
		{
          barPDG=91000000;
          resPDG=lResPDG;
          barM  =mLamb;
          resM  =lResM;
		}
        else if(thePDG!=90004004&&bN>1&&bZ>1&&bA>4&&tMC>aResM+mAlph)// Decay in alpha (p4)
		{
          barPDG=90002002;
          resPDG=aResPDG;
          barM  =mAlph;
          resM  =aResM;
		}
		else if(dResPDG&&tMC>dResM+mDeut)  // Can radiate a Deuteron (priority 5)
		{
          barPDG=90001001;
          resPDG=dResPDG;
          barM  =mDeut;
          resM  =dResM;
		}
		else if(ppResPDG&&tMC>ppResM+mProt+mProt)// Can radiate a DiProton (priority 6)
		{
          barPDG=90001000;
          resPDG=ppResPDG;
          thdPDG=90001000;
          barM  =mProt;
          resM  =ppResM;
          thdM  =mProt;
		}
		else if(nnResPDG&&tMC>nnResM+mNeut+mNeut)// Can radiate a DiNeutron (priority 7)
		{
          barPDG=90000001;
          resPDG=nnResPDG;
          thdPDG=90000001;
          barM  =mNeut;
          resM  =nnResM;
		}
		else if(npResPDG&&tMC>npResM+mProt+mNeut)// Can radiate a neutron+proton (prior 8)
		{
          barPDG=90001000;
          resPDG=npResPDG;
          thdPDG=90000001;
          barM  =mProt;
          resM  =npResM;
		}
		else if(lnResPDG&&tMC>lnResM+mLamb+mNeut)// Can radiate a Lambda+neutron (prior 9)
		{
          barPDG=91000000;
          resPDG=lnResPDG;
          thdPDG=90000001;
          barM  =mLamb;
          resM  =lnResM;
		}
		else if(lpResPDG&&tMC>lpResM+mLamb+mProt)// Can radiate a Lambda+proton (prior 10)
		{
          barPDG=91000000;
          resPDG=lpResPDG;
          thdPDG=90001000;
          barM  =mLamb;
          resM  =lpResM;
          thdM  =mProt;
		}
		else if(llResPDG&&tMC>llResM+mLamb+mLamb)// Can radiate a DiLambda (priority 11)
		{
          barPDG=91000000;
          resPDG=llResPDG;
          thdPDG=91000000;
          barM  =mLamb;
          resM  =llResM;
          thdM  =mLamb;
		}
        else if(thePDG!=90004004&&tMC>GSMass)// If it's not Be8 decay in gamma & GSM
		{
          barPDG=thePDG;
          resPDG=22;
          barM  =GSMass;
          resM  =0.;
		}
        else if(thePDG!=90004004)
		{
          G4cerr<<"***G4QEnv::EvaR:PDG="<<thePDG<<",M="<<totMass<<"< GSM="<<GSMass<<G4endl;
          throw G4QException("G4QEnvironment::EvaporateResidual:M<GSM & can't decayInPNL");
		}
        G4LorentzVector a4Mom(0.,0.,0.,barM);
        G4LorentzVector b4Mom(0.,0.,0.,resM);
        if(!thdPDG)
        {
          if(!qH->DecayIn2(a4Mom,b4Mom))
          {
            theQHadrons.push_back(qH);          // Fill as it is (delete equivalent)
            G4cerr<<"---Worning---G4QEnv::EvaRes:rP="<<pResPDG<<",rN="<<nResPDG<<",rL="
                  <<lResPDG<<",N="<<bN<<",Z="<<bZ<<",L="<<bS<<",totM="<<totMass<<",n="
                  <<totMass-nResM-mNeut<<",p="<<totMass-pResM-mProt<<",l="
                  <<totMass-lResM-mLamb<<G4endl;
            G4cerr<<"---Worning---G4QE::EvR:DecIn2Error b="<<barPDG<<",r="<<resPDG<<G4endl;
	      }
          else
          {
            delete qH;
            G4QHadron* HadrB = new G4QHadron(barPDG,a4Mom);
#ifdef fdebug
            G4cout<<"G4QE::EvaRes:(1) Baryon="<<barPDG<<a4Mom<<G4endl;
#endif
            theQHadrons.push_back(HadrB);       // Fill the baryon (delete equivalent)
            G4QHadron* HadrR = new G4QHadron(resPDG,b4Mom);
#ifdef fdebug
            G4cout<<"G4QE::EvaRes:(1) Residual="<<resPDG<<b4Mom<<G4endl;
#endif
            // @@ Self-call !!
            if(HadrR->GetBaryonNumber()>1) EvaporateResidual(HadrR);//Continue decay(del.e)
            else theQHadrons.push_back(HadrR);  // Fill ResidNucl=Baryon to OutHadronVector
          }
        }
        else
        {
          G4LorentzVector c4Mom(0.,0.,0.,thdM);
          if(!qH->DecayIn3(a4Mom,b4Mom,c4Mom))
          {
            theQHadrons.push_back(qH);          // Fill as it is (delete equivalent)
            G4cout<<"---Worning---G4QE::EvR:rNN="<<nnResPDG<<",rNP="<<npResPDG<<",rPP="
                  <<ppResPDG<<",N="<<bN<<",Z="<<bZ<<",L="<<bS<<",tM="<<totMass<<",nn="
                  <<totMass-nnResM-mNeut-mNeut<<",np="<<totMass-npResM-mProt-mNeut<<",pp="
                  <<totMass-ppResM-mProt-mProt<<G4endl;
            G4cerr<<"---Worning---G4QE::EvR:DecIn2Error,b="<<barPDG<<",r="<<resPDG<<G4endl;
	      }
          else
          {
            delete qH;
            G4QHadron* HadrB = new G4QHadron(barPDG,a4Mom);
#ifdef fdebug
            G4cout<<"G4QE::EvaRes:(2) Baryon1="<<barPDG<<a4Mom<<G4endl;
#endif
            theQHadrons.push_back(HadrB);       // Fill the first baryon (del.equiv.)
            G4QHadron* HadrC = new G4QHadron(thdPDG,c4Mom);
#ifdef fdebug
            G4cout<<"G4QE::EvaRes:(2) Baryon2="<<thdPDG<<c4Mom<<G4endl;
#endif
            theQHadrons.push_back(HadrC);       // Fill the second baryon (del.equiv.)
            G4QHadron* HadrR = new G4QHadron(resPDG,b4Mom);
#ifdef fdebug
            G4cout<<"G4QE::EvaRes:(2) Residual="<<resPDG<<b4Mom<<G4endl;
#endif
            // @@ Self-call !!
            if(HadrR->GetBaryonNumber()>1) EvaporateResidual(HadrR);//Continue decay(DelEq)
            else theQHadrons.push_back(HadrR); // Fill ResidNucl=Baryon to OutputHadrVector
          }
        }
	  }
      else if (abs(totMass-GSMass)<.003) theQHadrons.push_back(qH);// FillAsItIs (del.eq.)
      else                             // "System is below mass shell and can't decay" case
	  {
#ifdef fdebug
        G4cout<<"***G4QEnv::EvaRes: tM="<<totMass<<"("<<thePDG<<") < GSM="<<GSMass<<", d="
              <<totMass-GSMass<<", QC="<<qH->GetQC()<<qH->Get4Momentum()<<G4endl;
#endif
        //@@ Why this does not work? - Wait for the close message
        //G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());
        //theEnvironment=G4QNucleus(90000000,G4LorentzVector(0.,0.,0.,0.));
		//if(!CheckGroundState(quasH,true)) theQHadrons.push_back(qH);//Correct or FillAsIs
        //else delete qH;  
        //delete quasH;
        // @@ Temporary
        theQHadrons.push_back(qH); // Correct or fill as it is
      }
    }
    //else if(bA==5) DecayAlphaBar(qH);// Decay alpha-nucleon state (delete equivalent)
    //else if(bZ==4&&bN==2&&!bS) DecayAlphaDiN(qH); // Decay alpha+2protons state (del.eq.)
    //else if(bZ==4&&bN==4&&!bS) DecayAlphaAlpha(qH); // Decay alpha+alpha state (del.eq.)
    else                                        // ===> Evaporation of excited system
	{
#ifdef pdebug
      G4cout<<"G4QE::EvaRes:***EVA***tPDG="<<thePDG<<",M="<<totMass<<">GSM="<<GSMass<<",d="
            <<totMass-GSMass<<", N="<<qNuc.Get4Momentum()<<qNuc.Get4Momentum().m()<<G4endl;
#endif
      G4LorentzVector b4M;
      G4LorentzVector r4M;
      G4bool evC=true;
      G4bool dcC=false;
      G4int bPDG=0;
      G4int rPDG=0;
      G4double bM   = 0.;               // Prototype of Real Mass of the EvaporatedDibaryon
      G4double rM   = 0.;               // Prototype of Real Mass of the residual nucleus
      G4int bB=0;                       // Proto of Baryon Number of the evaporated baryon
      G4int rB=0;                       // Proto of Baryon Number of the residual nucleus
      G4QHadron* bHadron = new G4QHadron;
      G4QHadron* rHadron = new G4QHadron;
      G4int evcn=0;
      //G4int evcm=27;
      G4int evcm=9;
      while(evC&&evcn<evcm)
      {
        evC=true;
        dcC=false;
        evcn++;
        if(!qNuc.EvaporateBaryon(bHadron,rHadron))
	    {
#ifdef pdebug
          G4cout<<"***G4QEnv::EvaRes: ***EVA*** PDG="<<thePDG<<",tM="<<totMass<<G4endl;
#endif
          delete bHadron;
          delete rHadron;
          G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());
#ifdef fdebug
          G4cout<<"***G4QE::EvaRes: Residual="<<qH->GetQC()<<qH->Get4Momentum()<<G4endl;
#endif
          theEnvironment=G4QNucleus(90000000,G4LorentzVector(0.,0.,0.,0.));
		  if(!CheckGroundState(quasH,true)) theQHadrons.push_back(qH); // Cor or fillAsItIs
          else delete qH;  
          delete quasH;
          return;
          //evC=false;
          //dcC=true;
#ifdef ppdebug
          //@@Temporary
          //throw G4QException("G4QEnvironment::EvaporateResidual:Failed toEvap.bar/alph");
#endif
	    }
        if(!dcC)
        {
          evC=false;
          b4M=bHadron->Get4Momentum();
          r4M=rHadron->Get4Momentum();
          bM   = b4M.m();                       // Real mass of the evaporated dibaryon
          rM   = r4M.m();                       // Real mass of the residual nucleus
          bB=bHadron->GetBaryonNumber();        // Baryon number of the evaporated baryon
          rB=rHadron->GetBaryonNumber();        // Baryon number of the residual nucleus
          bPDG=bHadron->GetPDGCode();
          rPDG=rHadron->GetPDGCode();
#ifdef debug
          G4int bC=bHadron->GetCharge();        // Baryon number of the evaporated baryon
          //G4int rC=rHadron->GetCharge();       // Baryon number of the residual nucleus
          G4double bCB=qNuc.CoulombBarrier(bC,bB);
          //G4double rCB=qNuc.CoulombBarrier(rC,rB);
          G4cout<<"G4QEnv::EvaRes:Attempt #"<<evcn<<" > "<<evcm<<", rPDG="<<rPDG<<", bPDG="
                <<bPDG<<", bE="<<b4M.e()-b4M.m()<<" > bCB="<<bCB<<G4endl;
#endif
          //if(b4M.e()-b4M.m()<bCB&&evcn<evcm) evC=true;
		}
	  }  // End of while
      if(!dcC)
      {
#ifdef debug
        G4cout<<"G4QEnv::EvaRes:*** EVA IS DONE *** F="<<bPDG<<b4M<<",bB="<<bB<<", ResNuc="
              <<rPDG<<r4M<<",rB="<<rB<<G4endl;
#endif
        delete qH;
        if(bB<2)theQHadrons.push_back(bHadron); // Fill EvaporatedBaryon (del.equivalent)
        else if(bB==2) DecayDibaryon(bHadron);  // => "Dibaryon" case needs decay
        else if(bB==4) theQHadrons.push_back(bHadron); // "Alpha radiation" case (del.eq.)
        else if(bB==5) DecayAlphaBar(bHadron);  // "Alpha+Baryon Decay" case (del.equiv.)
        else if(bPDG==90004002) DecayAlphaDiN(bHadron); // Decay alph+2p(alph+2n is stable)
        else if(bPDG==90004004) DecayAlphaAlpha(bHadron);// "Alph+Alph Decay" case(del.eq.)
        else
	    {
          delete bHadron;
          G4cerr<<"***G4QE::EvaRes:bB="<<bB<<">2 - unexpected evaporated fragment"<<G4endl;
          throw G4QException("G4QEnvironment::EvaporateResidual: Wrong evaporation act");
	    }
        if(rB>2) EvaporateResidual(rHadron);    // Continue evaporation (@@ Self-call)
        else if(rB==2)                   // => "Dibaryon" case needs decay @@ DecayDibaryon
	    {
          G4double rGSM = rHadron->GetQPDG().GetMass();// Ground State mass of the dibaryon
#ifdef debug
		  G4cout<<"G4QE::EvaRes:ResidDibM="<<rM<<",GSM="<<rGSM<<",M-GSM="<<rM-rGSM<<G4endl;
#endif

          if(rM<=rGSM-0.001)
		  {
            delete rHadron;
            G4cerr<<"***G4QEnv::EvaRes: <residual> M="<<rM<<" < GSM="<<rGSM<<G4endl;
            throw G4QException("G4QEnvironment::EvaporateResidual:Evapor below MassShell");
		  }
          else if(abs(rM-rGSM)<0.001&&rPDG==90001001)theQHadrons.push_back(rHadron);//(D.E)
          else DecayDibaryon(rHadron);         // => "Dibaryon Decay" case (del.equivalent)
	    }
        else if(rB==5) DecayAlphaBar(rHadron); // "Alpha+Baryon Decay" case (del.equival.)
        else if(rPDG==90004002) DecayAlphaDiN(rHadron);//Decay alph+2p (alph+2n is stable)
        else if(rPDG==90004004) DecayAlphaAlpha(rHadron);// "Alpha+Alpha Decay" case(delEq)
        else theQHadrons.push_back(rHadron); // Fill ResidNucl=Baryon to OutputHadronVector
	  } // End of fail in decay check
      //else delete qH;
      else
      {
        G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());
        theEnvironment=G4QNucleus(90000000,G4LorentzVector(0.,0.,0.,0.));
		if(!CheckGroundState(quasH,true)) theQHadrons.push_back(qH); // Cor or fillAsItIs
        else delete qH;  
        delete quasH;
      }
	} // End of Evaporation of excited system
#ifdef debug
    G4cout<<"G4QEnv::EvaRes: === End of the evaporation attempt"<<G4endl;
#endif
  }
  else                                          // => "Decay if impossible evaporate" case
  {
#ifdef debug
    G4cout<<"G4QEnv::EvaRes: InputHadron4M="<<q4M<<", PDG="<<thePDG<<G4endl;
#endif
    if(thePDG)
    {
      if(thePDG==10)                            // "Chipolino decay" case 
	  {
        G4QContent totQC = qH->GetQC();         // Quark content of the hadron
        G4QChipolino resChip(totQC);            // define the Residual as a Chipolino
        G4QPDGCode h1=resChip.GetQPDG1();
        G4double m1  =h1.GetMass();             // Mass of the first hadron
        G4QPDGCode h2=resChip.GetQPDG2();
        G4double m2  =h2.GetMass();             // Mass of the second hadron
        if(totMass+.0001>m1+m2)
        {
          delete qH;                            // Chipolino should not be in a sequence
          G4LorentzVector fq4M(0.,0.,0.,m1);
          G4LorentzVector qe4M(0.,0.,0.,m2);
          if(!G4QHadron(q4M).DecayIn2(fq4M,qe4M))
		  {
            G4cerr<<"***G4QEnv::EvaRes:tM="<<totMass<<"-> h1M="<<m1<<" + h2M="<<m2<<G4endl;
		    throw G4QException("G4QEnvironment::EvaporateResid:Chip->h1+h2 DecIn2 error");
	      }
          G4QHadron* H2 = new G4QHadron(h2.GetPDGCode(),qe4M);
#ifdef fdebug
          G4cout<<"G4QE::EvaRes:(2) h2="<<h2.GetPDGCode()<<qe4M<<G4endl;
#endif
          theQHadrons.push_back(H2);            // (delete equivalent)
          G4QHadron* H1 = new G4QHadron(h1.GetPDGCode(),fq4M);
#ifdef fdebug
          G4cout<<"G4QE::EvaRes:(2) h1="<<h1.GetPDGCode()<<fq4M<<G4endl;
#endif
          theQHadrons.push_back(H1);            // (delete equivalent)
		}
        else
	    {
          delete qH;
          G4cerr<<"**G4QE::ER:M="<<totMass<<"<"<<m1<<"+"<<m2<<",d="<<m1+m2-totMass<<G4endl;
          throw G4QException("G4QEnvironment::EvaporateResidual:ChipolinoUnder MassShell");
	    }
	  }
      else                                      // "Hadron" case
	  {
        G4double totM=G4QPDGCode(thePDG).GetMass();
        if(abs(totMass-totM)<0.001||abs(thePDG)-10*(abs(thePDG)/10)>2)
                                                                 theQHadrons.push_back(qH);
        else if ((thePDG==221||thePDG==331)&&totMass>mPi+mPi) // "Decay in pipi" case
	    {
          G4LorentzVector fq4M(0.,0.,0.,mPi);
          G4LorentzVector qe4M(0.,0.,0.,mPi);
          if(!G4QHadron(q4M).DecayIn2(fq4M,qe4M))
		  {
            G4cerr<<"***G4QEnv::EvaporateResidual: tM="<<totMass<<"-> pi+ + pi-"<<G4endl;
		    throw G4QException("G4QEnv::EvaporateResidual: H->Pi+Pi DecayIn2 error");
	      }
          delete qH;
          G4QHadron* H1 = new G4QHadron(211,fq4M);
#ifdef fdebug
          G4cout<<"G4QE::EvaRes:(3) PiPlus="<<fq4M<<G4endl;
#endif
          theQHadrons.push_back(H1);            // (delete equivalent)
          G4QHadron* H2 = new G4QHadron(-211,qe4M);
#ifdef fdebug
          G4cout<<"G4QE::EvaRes:(3) PiMinus="<<qe4M<<G4endl;
#endif
          theQHadrons.push_back(H2);            // (delete equivalent)
	    }
        else if ((thePDG==221||thePDG==331)&&totMass>mPi0+mPi0) // "Decay in 2pi0" case
	    {
          G4LorentzVector fq4M(0.,0.,0.,mPi0);
          G4LorentzVector qe4M(0.,0.,0.,mPi0);
          if(!G4QHadron(q4M).DecayIn2(fq4M,qe4M))
		  {
            G4cerr<<"***G4QEnv::EvaporateResidual: tM="<<totMass<<"-> pi0 + pi0"<<G4endl;
		    throw G4QException("G4QEnv::EvaporateResidual: H->Pi+Pi DecayIn2 error");
	      }
          delete qH;
          G4QHadron* H1 = new G4QHadron(111,fq4M);
#ifdef fdebug
          G4cout<<"G4QE::EvaRes:(4) Pi01="<<fq4M<<G4endl;
#endif
          theQHadrons.push_back(H1);            // (delete equivalent)
          G4QHadron* H2 = new G4QHadron(111,qe4M);
#ifdef fdebug
          G4cout<<"G4QE::EvaRes:(4) Pi02="<<qe4M<<G4endl;
#endif
          theQHadrons.push_back(H2);            // (delete equivalent)
	    }
        else if (totMass>totM)                  // "Radiative Hadron decay" case
	    {
          G4LorentzVector fq4M(0.,0.,0.,0.);
          G4LorentzVector qe4M(0.,0.,0.,totM);
          if(!G4QHadron(q4M).DecayIn2(fq4M,qe4M))
		  {
            G4cerr<<"***G4QEnv::EvaporateRes:tM="<<totMass<<"->h1M="<<totM<<"+gam"<<G4endl;
		    throw G4QException("G4QEnv::EvaporateResidual:H*->H+gamma DecIn2 error");
	      }
          delete qH;
          G4QHadron* H2 = new G4QHadron(thePDG,qe4M);
#ifdef fdebug
          G4cout<<"G4QE::EvaRes:(5) tot="<<thePDG<<qe4M<<G4endl;
#endif
          theQHadrons.push_back(H2);            // (delete equivalent)
          G4QHadron* H1 = new G4QHadron(22,fq4M);
#ifdef fdebug
          G4cout<<"G4QE::EvaRes:(5) GamFortot="<<fq4M<<G4endl;
#endif
          theQHadrons.push_back(H1);            // (delete equivalent)
	    }
        else if (thePDG==111||thePDG==221||thePDG==331) // "Gamma+Gamma decay" case
	    {
          G4LorentzVector fq4M(0.,0.,0.,0.);
          G4LorentzVector qe4M(0.,0.,0.,0.);
          if(!G4QHadron(q4M).DecayIn2(fq4M,qe4M))
		  {
            G4cerr<<"***G4QEnv::EvaporateResidual:tM="<<totMass<<"->gamma + gamma"<<G4endl;
		    throw G4QException("G4QEnv::EvaporateResidual:pi/eta->g+g DecIn2 error");
	      }
          delete qH;
          G4QHadron* H2 = new G4QHadron(22,qe4M);
#ifdef fdebug
          G4cout<<"G4QE::EvaRes:(6) gam1="<<qe4M<<G4endl;
#endif
          theQHadrons.push_back(H2);            // (delete equivalent)
          G4QHadron* H1 = new G4QHadron(22,fq4M);
#ifdef fdebug
          G4cout<<"G4QE::EvaRes:(6) gam2="<<fq4M<<G4endl;
#endif
          theQHadrons.push_back(H1);            // (delete equivalent)
	    }
        else
	    {
          G4cerr<<"***G4QEnv::EvaRes: RNuc="<<thePDG<<theQC<<", q4M="<<q4M<<", M="<<totMass
                <<" < GSM="<<totM<<", 2Pi="<<mPi+mPi<<", 2Pi0="<<mPi0+mPi0<<G4endl;
          throw G4QException("G4QEnvironment::EvaporateResidual:Hadron's under MassShell");
	    }
	  }
	}
    else
    {
      G4cerr<<"***G4QE::EvaRes: RN="<<thePDG<<theQC<<",q4M="<<q4M<<",qM="<<totMass<<G4endl;
      throw G4QException("G4QEnv::EvaporateResidual: This is not a nucleus nor a hadron");
    }
  }
#ifdef debug
  G4cout<<"G4QEnv::EvaRes: ====>>>> End of the EvaporateResidual function"<<G4endl;
#endif
} // End of EvaporateResidual

//Make Random Unit 3D-Vector
G4ThreeVector G4QEnvironment::RndmDir()
{//  ==================================
  G4double x = G4UniformRand();
  G4double y = G4UniformRand();
  G4double z = G4UniformRand();
  G4double r2= x*x+y*y+z*z;
  while(r2>1.||r2<.01)
  {
    x = G4UniformRand();
    y = G4UniformRand();
    z = G4UniformRand();
    r2=x*x+y*y+z*z;
  }
  G4double r=sqrt(r2);
  return G4ThreeVector(x/r,y/r,z/r);
} // End of RndmDir

//The public Hadronisation function with the Exception treatment (del respons. of User !)
G4QHadronVector* G4QEnvironment::Fragment()
{//              ==========================
  G4QHadronVector dummy;       // Prototype of the output G4QHadronVector to avoid wornings
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
      theFragments = FSInteraction();     //InterClass creation. User must delet QHadrons.
    }
    catch (G4QException& error)
    {
      G4cerr<<"***G4QEnvironment::Fragment: Exception is catched"<<G4endl;
      RepFlag=true;                       // For the Exception - repete
      ExCount++;                          // Increment the repetition counter
      G4cerr<<"***G4QEnv::Fragment:Exception #"<<ExCount<<": "<<error.GetMessage()<<G4endl;
      G4LorentzVector dif=rem4M-theEnvironment.Get4Momentum(); // CHECK difference
      G4int nHp=theQHadrons.size();
      G4int nQp = theQuasmons.size();
      G4cerr<<"***G4QEnvir::Fragment:nH="<<nHp<<",nQ="<<nQp<<",E="<<theEnvironment<<G4endl;
      for(G4int ph=0; ph<nHp; ph++)
      {
        G4QHadron* cH = theQHadrons[ph];
        dif-=cH->Get4Momentum();
		G4cerr<<"***G4QEnvir::Fr:H"<<ph<<"="<<cH->Get4Momentum()<<cH->GetPDGCode()<<G4endl;
      }
      for(G4int pq=0; pq<nQp; pq++)
      {
        G4Quasmon* cQ = theQuasmons[pq];
        dif-=cQ->Get4Momentum();
        G4cerr<<"***G4QEnvir::Fr:Quasm#"<<pq<<"="<<cQ->Get4Momentum()<<cQ->GetQC()<<G4endl;
      }
      // *** Cleaning Up of all old output instances for the recalculation purposes ***
      std::for_each(theFragments->begin(), theFragments->end(), DeleteQHadron());//oldHadrs
      theFragments->clear();
      std::for_each(theQHadrons.begin(), theQHadrons.end(), DeleteQHadron());//intern.Hadrs
      theQHadrons.clear();
      std::for_each(theQuasmons.begin(), theQuasmons.end(), DeleteQuasmon());//old Quasmons
      theQuasmons.clear();
      G4cerr<<"***G4QEnv::Fragment: ----------- End of CleaningUp: 4Mdif="<<dif<<G4endl;
      // **************** Recover all conditions for the recalculation ********************
      theEnvironment=reEnvironment;             // Recover the nuclear environment
      tot4Mom=rem4M;                            // Recover the total 4Momentum of the React
      G4cerr<<"***G4QEnv::Fragment:*Recover*Env="<<theEnvironment<<",4M="<<tot4Mom<<G4endl;
      G4int mQ = reQuasmons->size();            // Recover the memorizedQuasmons with print
      for(G4int jq=0; jq<mQ; jq++)
      {
        G4Quasmon* curQ = new G4Quasmon(reQuasmons->operator[](jq));
        G4cerr<<"***G4QE::Fragm:Q("<<jq<<")="<<curQ->Get4Momentum()<<curQ->GetQC()<<G4endl;
        theQuasmons.push_back(curQ);            // (delete equivalent)
      }
      G4int mH = reQHadrons->size();            // Recover the memorizedQHadrons with print
      for(G4int jh=0; jh<mH; jh++)
      {
        G4QHadron* curH = new G4QHadron(reQHadrons->operator[](jh));
        G4cerr<<"***G4QE::Fragm:H("<<jh<<")="<<curH->Get4Momentum()<<curH->GetQC()<<G4endl;
        theQHadrons.push_back(curH);            // (delete equivalent)
      }
    }
  }
  if(reQuasmons->size()) // If something is still in memory then clean it up
  {
    std::for_each(reQuasmons->begin(),reQuasmons->end(),DeleteQuasmon());//CleanUp oldQuasm
    reQuasmons->clear();
  }
  delete reQuasmons;     // All temporary Quasmons memory is wiped out
  if(reQHadrons->size()) // If something is still in memory then clean it up
  {
    std::for_each(reQHadrons->begin(),reQHadrons->end(),DeleteQHadron());//CleanUp oldQHadr
    reQHadrons->clear();
  }
  delete reQHadrons;     // All temporary QHadrons memory is wiped out
  if(ExCount>=MaxExCnt)
  {
    G4cerr<<"***G4QEnv::Fragment: Exception is pushed Up to Hadronics *** ^^^ ***"<<G4endl;
    // G4Exception("G4QEnvironment::Fragment","027",FatalException,"*General CHIPSException");
    throw G4HadronicException(__FILE__, __LINE__, "*Fragment:Exception Up to Hadronics");
  }
  return theFragments;
} // End of the Fragmentation member function

//The Final State Interaction Filter
G4QHadronVector* G4QEnvironment::FSInteraction()
{//              ===============================
  static const G4QPDGCode gQPDG(22);
  static const G4QPDGCode pipQPDG(211);
  static const G4QPDGCode pimQPDG(-211);
  static const G4QPDGCode nQPDG(2112);
  static const G4QPDGCode pQPDG(2212);
  static const G4QPDGCode lQPDG(3122);
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
  static const G4double mSigM= G4QPDGCode(3112).GetMass();
  static const G4double mSigP= G4QPDGCode(3222).GetMass();
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
#ifdef pdebug
  G4cout<<"G4QEnvironment(G4QE)::FSInteraction(FSI): ***called*** envA="<<envA<<G4endl;
#endif
  G4QHadronVector* theFragments = new G4QHadronVector;//Internal creation. User must delete
  HadronizeQEnvironment();            // >>>>>>>>> Call the main fragmentation function
  unsigned nHadr=theQHadrons.size();
  if(nHadr<=0)
  {
    G4cerr<<"---Worning---G4QEnvironment::FSInteraction: nHadrons="<<nHadr<<G4endl;
	//throw G4QException("G4QEnvironment::FSInteraction: No hadrons in the output");
    return theFragments;
  }
  G4int lHadr=theQHadrons[nHadr-1]->GetBaryonNumber();
#ifdef pdebug
  G4cout<<"G4QE::FSI:after HQE,nH="<<nHadr<<",lHBN="<<lHadr<<",E="<<theEnvironment<<G4endl;
#endif
  if(lHadr>1)                          // The last hadron is nucleus: try to decay/evap. it
  {
    G4QHadron* theLast = theQHadrons[nHadr-1];
    G4QHadron* curHadr = new G4QHadron(theLast);
#ifdef fdebug
    G4cout<<"G4QE::FSI:B,nH="<<nHadr<<",PDG="<<curHadr->GetPDGCode()<<",E="<<theEnvironment
          <<G4endl;
#endif
    theQHadrons.pop_back();            // the last QHadron-Nucleus is excluded from OUTPUT
    delete theLast;//**!!When kill,DON'T forget to delete theLastQHadron as an instance!!**
    EvaporateResidual(curHadr);        // Try to evaporate Hadr-Nucl (@@DecDib)(delete eq.)
    nHadr=theQHadrons.size();
#ifdef pdebug
    G4cout<<"G4QE::FSI:After,n="<<nHadr<<",P="<<theQHadrons[nHadr-1]->GetPDGCode()<<G4endl;
#endif
  }
#ifdef fdebug
  G4LorentzVector ccs4M(0.,0.,0.,0.);           // CurrentControlSum of outgoing Hadrons
#endif
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
#ifdef fdebug
    G4int hNF  = theCurr->GetNFragments();
    G4cout<<"G4QE::FSI:h#"<<ipo<<":h="<<hPDG<<h4Mom<<",F="<<hNF<<",nH="<<theQHadrons.size()
          <<G4endl;
#endif
    if(hBN>lHadr) // Current Hadron is the Biggest fragment -> Swap with theLast Hadron
	{
      G4QHadron* theLast = theCurr;          // Prototype of the pointer to the Last Hadron
      G4QHadron* curHadr = new G4QHadron(theCurr);// Remember CurrentHadron for evaporation
      if(ipo+1<theQHadrons.size())          // If ipo<Last, swap theCurHadr and theLastHadr
      {
        theLast = theQHadrons[theQHadrons.size()-1];// Pointer to theLastHadron (ipo<Last)
        theCurr->SetQPDG(theLast->GetQPDG());// theCurHadron is substituted by the LastHadr
        theCurr->Set4Momentum(theLast->Get4Momentum()); // ... continue substitution
      }
      theQHadrons.pop_back();           // pointer to theLastHadron is excluded from OUTPUT
      delete theLast; //*!!When kill,DON'T forget to delete theLastQHadron asAnInstance !!*
      theQHadrons.push_back(curHadr);
      nHadr=theQHadrons.size();
      h4Mom = theCurr->Get4Momentum();
      hBN  = theCurr->GetBaryonNumber();
      cBN  = theCurr->GetCharge();
      sBN  = theCurr->GetStrangeness();
      hPDG = theCurr->GetPDGCode();
	}
    if(hPDG==89002000||hPDG==89001001||hPDG==89000002)// 2pt dec. of anti-strange (3pt dec)
    {
#ifdef fdebug
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
#ifdef fdebug
            G4cerr<<"**G4QE::FSI:ANTISTRANGE*++*STRANGENESS,PDG="<<hPDG<<",M="<<hM<<G4endl;
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
#ifdef fdebug
            G4cerr<<"**G4QE::FSI:*ANTISTRANGE*+*STRANGENESS*PDG="<<hPDG<<",M="<<hM<<G4endl;
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
#ifdef fdebug
            G4cerr<<"**G4QE::FSI:**ANTISTRANGE*0*STRANGENE**PDG="<<hPDG<<",M="<<hM<<G4endl;
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
#ifdef fdebug
        G4cerr<<"***G4QE::FSI:***ANTISTRANGE***CANN'T DECAY,PDG="<<hPDG<<",M="<<hM<<G4endl;
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
          G4cerr<<"---Worning---G4QE::FSI: Still try(2),M="<<hM<<"->"<<fM<<"("<<fQPDG<<")+"
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
              theCurr->SetQPDG(theLast->GetQPDG()); // CurHadrPDG instead of LastHadrPDG
              theCurr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
            }                                       //ELSE: it's already theLast -> no swap
            theQHadrons.pop_back();                 //exclude LastHadronPointer from OUTPUT
            delete theLast; //*!! When killing, DON'T forget to delete the last QHadron !!*
            G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());//Fake Quasmon
            //           3. Try to use other hadrons to recover this one (under Mass Shell)
            if(!CheckGroundState(quasH,true))       // Try to correct with other hadrons
            {
              G4cerr<<"---Worning---G4QEnv::FSI:Failed (2) LeaveAsItIs 4m="<<h4Mom<<G4endl;
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
            G4cerr<<"***G4QEnv::FSI: No recovery (1) Env="<<theEnvironment<<G4endl;
	        throw G4QException("G4QEnv::FSI:ANTISTRANGE DecayIn2 did not succeed");
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
          G4cerr<<"---Worning---G4QE::FSI: Still try(3), M"<<hM<<"->"<<fM<<"("<<fQPDG<<")+"
                <<sM<<"("<<sPDG<<")+"<<tM<<"("<<tPDG<<")="<<sum<<G4endl;
          //if(theEnvironment==vacuum)
         if(!theEnvironment.GetA())
          {
            G4QHadron* theLast = theCurr;    // Prototype of the pointer to the Last Hadron
            G4QHadron* qH = new G4QHadron(theCurr); // Copy of the Current Hadron
            if(ipo+1<theQHadrons.size())       // If ipo<Last, swap CurHadr and theLastHadr
            {
              theLast = theQHadrons[theQHadrons.size()-1];//Pointer to LastHadron(ipo<Last)
              theCurr->SetQPDG(theLast->GetQPDG()); // CurHadrPDG instead of theLastHadrPDG
              theCurr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
            }
            theQHadrons.pop_back();        // exclude theLastHadron pointer from the OUTPUT
            delete theLast; //*!! When killing, DON'T forget to delete the last QHadron !!*
            G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());//Fake Quasmon
            if(!CheckGroundState(quasH,true))         // Try to correct by other hadrons
            {
              G4cerr<<"---Worning---G4QEnv::FSI:Failed (3) LeaveAsItIs,4M="<<h4Mom<<G4endl;
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
            G4cerr<<"***G4QEnv::FSI: No recovery (2) Env="<<theEnvironment<<G4endl;
	        throw G4QException("G4QEnv::FSI:ANTISTRANGE DecayIn3 did not succeed");
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
#ifdef fdebug
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
      else if(hPDG==90999002)  // n,Lambda,pi-
	  {
        if(hM>mLamb+mNeut+mPi)
        {
          barPDG = 3122;         // Substitute Lambda for the second n
          barM   = mSigM;
        }
        else if(hM>mSigM+mNeut)  // @@ Decay in 2
        {
          barPDG = 3112;         // Substitute Sigma- for the second n
          barM   = mLamb;
          tPDG   = 22;
          tM     = 0;
        }
	  }
      else if(hPDG==91001999)  // p,Lambda,pi+
	  {
        nuQPDG = pQPDG;        // Substitute p for the first n
        nucM   = mProt;
        tPDG   = 211;          // Substitute pi+ for the first pi-
        if(hM>mLamb+mProt+mPi)
        {
          barPDG = 3122;         // Substitute Lambda for the second n
          barM   = mLamb;
        }
        else if(hM>mSigP+mProt)  // @@ Decay in 2
        {
          barPDG = 3222;         // Substitute Sigma- for the second n
          barM   = mLamb;
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
        G4int eA=theEnvironment.GetA();
#ifdef fdebug
        G4cerr<<"***G4QEnv::FSI:T="<<hPDG<<"("<<hM<<")-> N="<<nuQPDG<<"(M="<<nucM<<") + B="
              <<barPDG<<"("<<barM<<")+N/pi="<<tPDG<<"("<<tM<<")="<<sum<<", A="<<eA<<G4endl;
#endif
        if(!eA)
        {
          G4QHadron* theLast = theCurr;        // Prototype of the pointer to theLastHadron
          G4QHadron* qH = new G4QHadron(theCurr); // Copy of the Current Hadron
#ifdef fdebug
          G4cerr<<"***G4QE::FSI:#"<<ipo<<",4MQC="<<qH->Get4Momentum()<<qH->GetQC()<<G4endl;
#endif
          if(ipo+1<theQHadrons.size())         // If ipo<Last, swap CurHadr and theLastHadr
          {
            G4int nhd1=theQHadrons.size()-1;
            theLast = theQHadrons[nhd1];// Pointer to theLastHadron (ipo<L)
            G4LorentzVector l4M=theLast->Get4Momentum();
            G4QPDGCode lQPDG=theLast->GetQPDG();
#ifdef fdebug
			G4cerr<<"---Worning---G4QE::FSI:l#"<<nhd1<<",4M="<<l4M<<",PDG="<<lQPDG<<G4endl;
#endif
            theCurr->SetQPDG(theLast->GetQPDG());// CurHadrPDG is substituted by LastHadPDG
            theCurr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
          }
          theQHadrons.pop_back();              // exclude theLastHadron pointer from OUTPUT
          delete theLast; //*!! When killing, DON'T forget to delete the last QHadron !!*
          G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());//Crea fakeQuasm
          if(!CheckGroundState(quasH,true))         // Try to correct by other hadrons
          {
            G4cerr<<"---Worning---G4QEnv::FSI: Failed (3) FillAsItIs 4M="<<h4Mom<<G4endl;
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
          G4cerr<<"***G4QEnv::FSI: No recovery (3) Env="<<theEnvironment<<G4endl;
	      throw G4QException("G4QEnv::FSI:ISO-dibaryon or 3n/3p DecayIn3 did not succeed");
        }
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
#ifdef fdebug
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
#ifdef fdebug
        G4cerr<<"***G4QEnv::FSI:IsN M="<<hM<<","<<hPDG<<"->N="<<nuQPDG<<"(M="<<nucM<<")+"
              <<nN<<"*B="<<barPDG<<"(M="<<barM<<")+"<<nPi<<"*pi="<<tPDG<<"(M="<<tM<<")="
              <<nucM+barM+tM<<G4endl;
        G4cerr<<"***G4QEnv::FSI:IsN BaryN="<<hBN<<",Charge="<<cBN<<",Stran="<<sBN<<G4endl;
#endif
        if(!theEnvironment.GetA())           // Emergency recovery
        {
          G4QHadron* theLast = theCurr;      // Prototype of the pointer to the Last Hadron
          G4QHadron* qH = new G4QHadron(theCurr); // Copy of the Current Hadron
          if(ipo+1<theQHadrons.size())       // If ipo<Last, swap CurHadr and theLastHadr
          {
            theLast = theQHadrons[theQHadrons.size()-1];// Ptr to theLastHadron (ipo<Last)
            theCurr->SetQPDG(theLast->GetQPDG());// CurHadPDG is substituted by LastHadrPDG
            theCurr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
          }
          theQHadrons.pop_back();            // exclude theLastHadron pointer from OUTPUT
          delete theLast;//*!! When killing, DON'T forget to delete the last QHadron !!*
          G4QContent cqQC=qH->GetQC();
          G4LorentzVector cq4M=qH->Get4Momentum();
          G4Quasmon* quasH = new G4Quasmon(cqQC,cq4M);// Create fakeQuasmon for the Last
          if(!CheckGroundState(quasH,true))         // Try to correct by other hadrons
          {
            G4cerr<<"---Worning---G4QEnv::FSI:IN Failed, FillAsItIs: "<<cqQC<<cq4M<<G4endl;
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
          G4cerr<<"***G4QEnv::FSI:IN No recovery (3) Env="<<theEnvironment<<G4endl;
	      throw G4QException("G4QEnv::FSI:IN Multy ISO-dibaryon DecayIn3 did not succeed");
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
#ifdef fdebug
    G4int           hNFrag= theQHadrons[ipo]->GetNFragments(); //Recover after swapping
    G4QContent      hQC   = theQHadrons[ipo]->GetQC();         // ...
    hPDG                  = theQHadrons[ipo]->GetPDGCode();    // ...
    h4Mom                 = theQHadrons[ipo]->Get4Momentum();  // ...
    ccs4M+=h4Mom;                                              // Calculate CurSum of Hadrs
    G4cout<<"G4QE::FSI:#"<<ipo<<": h="<<hPDG<<hQC<<",h4M="<<h4Mom<<",hNF="<<hNFrag<<G4endl;
#endif
    ipo=jpo;            // Take into account the roll back in case of the Last substitution
  }
#ifdef fdebug
  G4cout<<"G4QE::FSI: >>>CurrentControlSumOf4MomOfHadrons="<<ccs4M<<G4endl;
#endif
  G4double p2cut=250000./envA/envA; // 250000=(2*p_Ferm)**2
  //G4double p2cut2=0.; //cut for the alpha creation
  //
  G4int bfCountM=3;
  if(envA>10) bfCountM*=(envA-1)/3;
  G4bool bfAct = true;
  G4int bfCount= 0;
  G4LorentzVector tmp4Mom=tot4Mom;
  while(bfAct&&bfCount<bfCountM) // "Infinite" LOOP of the ThermoNuclearBackFusion
  {
    tot4Mom=tmp4Mom;
    bfAct=false;
    bfCount++;
    nHadr=theQHadrons.size();
    if(nHadr) for(unsigned hadron=0; hadron<theQHadrons.size(); hadron++)// BackFusion LOOP
    {
      G4QHadron* curHadr = theQHadrons[hadron]; // Get a pointer to the current Hadron
      G4int         hPDG = curHadr->GetPDGCode();
      G4QPDGCode    hQPDG(hPDG);
      G4double      hGSM = hQPDG.GetMass();     // Ground State Mass of the first fragment
#ifdef fdebug
      G4cout<<"G4QE::FSI:LOOP START,h#"<<hadron<<curHadr->Get4Momentum()<<hPDG<<G4endl;
#endif
      if(hPDG==89999003||hPDG==90002999)
      {
        G4cerr<<"---WORNING---G4QE::FSI:**nD-/pD++**(3),PDG="<<hPDG<<" CORRECTION"<<G4endl;
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
          G4cerr<<"---WORNING---G4QE::FSI: Still trying, NDM="<<hM<<"->"<<fM<<"("<<fQPDG
                <<")+"<<sM<<"("<<sPDG<<")+"<<tM<<"("<<tPDG<<")="<<sum<<G4endl;
          if(!theEnvironment.GetA())
          {
            G4QHadron* theLast = curHadr;          // Prototype of Pointer to theLastHadron
            G4QHadron* qH = new G4QHadron(curHadr);// Copy of the Current Hadron
            if(hadron+1<theQHadrons.size())        // If hadr<Last,swap CurHadr & LastHadr
            {
              theLast = theQHadrons[theQHadrons.size()-1]; // Pointer to LastHadr (nh<Last)
              curHadr->SetQPDG(theLast->GetQPDG());// CurHadrPDG's substituted by LastHaPDG
              curHadr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
            }
            theQHadrons.pop_back();        // exclude theLastHadron pointer from the OUTPUT
            delete theLast; //*!! When killing, DON'T forget to delete the last QHadron !!*
            G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());//Fake Quasmon
            if(!CheckGroundState(quasH,true))      // Try to correct by other hadrons
            {
              G4cerr<<"---Worning---G4QE::FSI:NDel Failed LeaveAsItIs, 4m="<<h4Mom<<G4endl;
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
            G4cerr<<"***G4QEnv::FSI: No ND recovery Env="<<theEnvironment<<G4endl;
	        throw G4QException("G4QEnv::FSI:ND DecayIn3 did not succeed");
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
        G4cerr<<"---WORNING---G4QE::FSI:***(K+N)*** (2),PDG="<<hPDG<<" CORRECTION"<<G4endl;
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
          G4cerr<<"---WORNING---G4QE::FSI: Still trying (2),NDM="<<hM<<"->"<<fM<<"("<<fQPDG
                <<")+"<<sM<<"("<<sPDG<<")="<<sum<<G4endl;
          if(!theEnvironment.GetA())
          {
            G4QHadron* theLast = curHadr;          // Prototype of Pointer to theLastHadron
            G4QHadron* qH = new G4QHadron(curHadr);// Copy of the Current Hadron
            if(hadron+1<theQHadrons.size())        // If hadr<Last, swap CurHadr & LastHadr
            {
              theLast = theQHadrons[theQHadrons.size()-1]; // Pointer to LastHadr (nh<Last)
              curHadr->SetQPDG(theLast->GetQPDG());// CurHadrPDG's substituted by LastHaPDG
              curHadr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
            }
            theQHadrons.pop_back();        // exclude theLastHadron pointer from the OUTPUT
            delete theLast; //*!! When killing, DON'T forget to delete the last QHadron !!*
            G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());//Fake Quasmon
            if(!CheckGroundState(quasH,true))      // Try to correct by other hadrons
            {
              G4cerr<<"---Worning---G4QE::FSI:KN Failed LeaveAsItIs 4m="<<h4Mom<<G4endl;
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
            G4cerr<<"***G4QEnv::FSI: No KN recovery Env="<<theEnvironment<<G4endl;
	        throw G4QException("G4QEnv::FSI:KN DecayIn2 did not succeed");
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
#ifdef fdebug
      if(!hF&&(hPDG>80000000&&hPDG<90000000||hPDG==90000000||
               hPDG>90000000&&(hPDG%1000000>200000||hPDG%1000>300)))
        G4cerr<<"**G4QEnv::FSInteraction: PDG("<<hadron<<")="<<hPDG<<", M="<<hM<<G4endl;
#endif
#ifdef fdebug
      G4cout<<"G4QE::FSI:h="<<hPDG<<",S="<<hS<<",B="<<hB<<",#"<<hadron<<"<"<<nHadr<<G4endl;
#endif
	  //if(hadron&&!hF&&hB>0&&!hS&&(nHadr>3||hB<2)) // ThermoBackFus (VIMP for gamA TotCS)
	  //if(hadron&&!hF&&hB>0&&!hS&&(nHadr>2||hB<4)) // ThermoBackFus (VIMP for gamA TotCS)
	  if(hadron&&!hF&&hB>0&&!hS) // ThermoBackFusion cond. (VIMP for gamA TotCS)
	  //if(hadron&&!hF&&hB>0&&hB<4&&!hS) // ThermoBackFusion cond. (VIMP for gamA TotCS)
	  //if(hadron&&!hF&&hB>0&&!hS&&nHadr>2)//ThermoBackFusion MAX condition (VIMP for gamA)
	  //if(2>3)                         // Close the ThermoBackFusion (VIMP for gamA TotCS)
      {
#ifdef fdebug
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
#ifdef fdebug
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
            G4QPDGCode hQPDG=gQPDG;
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
            else if(sPDG==90002003)                     // A=5
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
#ifdef fdebug
            G4cout<<"G4QE::FSI: "<<three<<",r="<<rQPDG<<",f="<<fQPDG<<",t="<<hQPDG<<G4endl;
#endif
            if(!three)
            {
              if(!G4QHadron(s4M).DecayIn2(f4Mom,g4Mom))
              {
                G4cerr<<"***G4QE::FSI:(2)*FUSION*,tM["<<sPDG<<"]="<<tM<<">sM="<<sM<<" of "
                      <<h4m<<hM<<hQC<<hGSM<<" & "<<b4m<<bM<<bQC<<bGSM<<G4endl;
                throw G4QException("***G4QEnvironment::FSInter:Fusion (1) DecIn2 error");
              }
              else
			  {
#ifdef fdebug
                G4cout<<"G4QE::FSI:*FUSION IS DONE*,fPDG="<<sPDG<<",PDG1="<<hPDG<<",PDG2="
                      <<bPDG<<G4endl;
#endif
                curHadr->SetQPDG(fQPDG);
                curHadr->Set4Momentum(f4Mom);
                backH->SetQPDG(rQPDG);
                backH->Set4Momentum(g4Mom);
#ifdef fdebug
                G4cout<<"G4QE::FSI:h="<<h4m<<",b="<<b4m<<",s="<<s4M<<G4endl;
                G4cout<<"G4QE::FSI:f="<<f4Mom<<",g="<<g4Mom<<",s="<<f4Mom+g4Mom<<G4endl;
#endif
              }
		    }
            else
            {
              if(!G4QHadron(s4M).DecayIn3(f4Mom,g4Mom,t4Mom))
              {
                G4cerr<<"***G4QE::FSI:(3)*FUSION*,tM["<<sPDG<<"]="<<tM<<">sM="<<sM<<" of "
                      <<h4m<<hM<<hQC<<hGSM<<" & "<<b4m<<bM<<bQC<<bGSM<<G4endl;
                throw G4QException("G4QEnvironment::FSInteract:Fusion(2) DecayIn3 error");
              }
              else
			  {
#ifdef fdebug
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
#ifdef fdebug
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
#ifdef fdebug
            G4cout<<"G4QE::FSI:cH4M="<<curHadr->Get4Momentum()<<G4endl;
#endif
		  } // End of the fusion check
	    } // End of the LOOP over previous hadrons
	  } // End of the FUSION check
      G4LorentzVector cH4Mom = curHadr->Get4Momentum(); // 4-mom of the current hadron
      tot4Mom-=cH4Mom;                          // Reduce theTotal 4mom by theCurrent 4mom
      if(hadron+1==nHadr)                       // The Last Hadron in the set
	  {
        G4double re=tot4Mom.e();
        G4double rpx=tot4Mom.px();
        G4double rpy=tot4Mom.py();
        G4double rpz=tot4Mom.pz();
        G4double dmo=rpx*rpx+rpy*rpy+rpz*rpz;
        G4double dem=re*re+dmo;
#ifdef fdebug
        G4cout<<"G4QE::FSI: Is Energy&Mom conserved? t4M="<<tot4Mom<<",d2="<<dem<<G4endl;
#endif
		//if(2>3)                                 //@@En/Mom conservation check is closed
        if(dem>0.0001)                          // Energy or momentum is not conserved
	    {
          //#ifdef pdebug
          if(dem>.1)
          {
            G4cerr<<"---Worning---G4QE::FSI:E/Mcons4M="<<tot4Mom<<dem<<".Correct!"<<G4endl;
            G4LorentzVector cor4M(0.,0.,0.,mNeut);  // Prototype for the neutron
            if(dmo<0.0001 && re>900.)               // MomentumIsConserved - recoverMissing
			{
              if(fabs(re-mNeut)<.01)
			  {
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by neutron"<<G4endl;
                //cor4M=G4LorentzVector(0.,0.,0.,mNeut); // Made as a prototype
                G4QHadron* theH = new G4QHadron(90000001,G4LorentzVector(0.,0.,0.,mNeut));
                theQHadrons.push_back(theH);        // (delete equivalent for the proton)
              }
			  else if(fabs(re-mProt)<.01)
			  {
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by proton"<<G4endl;
                cor4M=G4LorentzVector(0.,0.,0.,mProt);
                G4QHadron* theH = new G4QHadron(90001000,G4LorentzVector(0.,0.,0.,mProt));
                theQHadrons.push_back(theH);        // (delete equivalent for the proton)
              }
			  else if(fabs(re-mDeut)<.01)
			  {
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by deuteron"<<G4endl;
                cor4M=G4LorentzVector(0.,0.,0.,mDeut);
                G4QHadron* theH = new G4QHadron(90001001,G4LorentzVector(0.,0.,0.,mDeut));
                theQHadrons.push_back(theH);        // (delete equivalent for the proton)
              }
			  else if(fabs(re-mTrit)<.01)
			  {
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by tritium"<<G4endl;
                cor4M=G4LorentzVector(0.,0.,0.,mTrit);
                G4QHadron* theH = new G4QHadron(90001002,G4LorentzVector(0.,0.,0.,mTrit));
                theQHadrons.push_back(theH);        // (delete equivalent for the proton)
              }
			  else if(fabs(re-mHe3)<.01)
			  {
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by He3"<<G4endl;
                cor4M=G4LorentzVector(0.,0.,0.,mHe3);
                G4QHadron* theH = new G4QHadron(90002001,G4LorentzVector(0.,0.,0.,mHe3));
                theQHadrons.push_back(theH);        // (delete equivalent for the proton)
              }
			  else if(fabs(re-mAlph)<.01)
			  {
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by alpha"<<G4endl;
                cor4M=G4LorentzVector(0.,0.,0.,mAlph);
                G4QHadron* theH = new G4QHadron(90002002,G4LorentzVector(0.,0.,0.,mAlph));
                theQHadrons.push_back(theH);        // (delete equivalent for the proton)
              }
			  else if(fabs(re-mNeut-mNeut)<.01)
			  {
                cor4M=G4LorentzVector(0.,0.,0.,mNeut+mNeut);
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by 2 neutrons"<<G4endl;
                G4QHadron* theH1 = new G4QHadron(90000001,G4LorentzVector(0.,0.,0.,mNeut));
                theQHadrons.push_back(theH1);        // (delete equivalent for the proton)
                G4QHadron* theH2 = new G4QHadron(90000001,G4LorentzVector(0.,0.,0.,mNeut));
                theQHadrons.push_back(theH2);        // (delete equivalent for the proton)
              }
			  else if(fabs(re-mProt-mProt)<.01)
			  {
                G4cout<<"...G4QE::FSI:E/M conservation is corrected by 2 protons"<<G4endl;
                cor4M=G4LorentzVector(0.,0.,0.,mProt+mProt);
                G4QHadron* theH1 = new G4QHadron(90001000,G4LorentzVector(0.,0.,0.,mProt));
                theQHadrons.push_back(theH1);        // (delete equivalent for the proton)
                G4QHadron* theH2 = new G4QHadron(90001000,G4LorentzVector(0.,0.,0.,mProt));
                theQHadrons.push_back(theH2);        // (delete equivalent for the proton)
              }
			  else throw G4QException("***G4QEnv::FSInteract: Increase theCorrectionEps");
			}
            else
            {
              // @@ Temporary very bad correction: just to pass through @@
              curHadr->Set4Momentum(cH4Mom+tot4Mom); // modify 4-mom of the Last hadron
              G4cout<<"*!* G4QE::FSI:E/M dif is added to theLast="<<tot4Mom+cH4Mom<<G4endl;
              if(curHadr->GetBaryonNumber()>1)
              {
                G4QHadron* theNew  = new G4QHadron(curHadr); // Make NewHadr of theLastHadr
                theQHadrons.pop_back();      // the last QHadron is excluded from theOUTPUT
                delete curHadr;//!!When kill,DON'T forget to delete theLastQH as an inst.!!
                EvaporateResidual(theNew); // Try to evaporate Residual
              }
              //throw G4QException("***G4QEnv::FSInteract: Energy/Momentum conservation");
            }
            tot4Mom=tot4Mom-cor4M;
            G4cerr<<"---Worning---G4QE::FSI:En/MomCons.Error is corrected:"<<cor4M<<G4endl;
          }
		  //#endif
          G4QHadron* prevHadr = theQHadrons[nHadr-2]; // GetPointer to Hadr prev to theLast
          G4LorentzVector pH4Mom = prevHadr->Get4Momentum(); // 4-mom of thePreviousHadron
          G4double cHM = curHadr->GetMass();    // Mass of the current hadron
          G4double pHM = prevHadr->GetMass();   // Mass of the current hadron
          tot4Mom+=cH4Mom+pH4Mom;
          G4double totRM=tot4Mom.m();
          if(cHM+pHM<=totRM+.03)                // *** Make the final correction ***
		  {
            if(!G4QHadron(tot4Mom).DecayIn2(pH4Mom,cH4Mom))
            {
              G4cerr<<"***G4QE::FSI:**Correction**tot4M="<<tot4Mom<<totRM<<">sM="<<cHM+cHM
                    <<G4endl;
              throw G4QException("***G4QEnvironment::FSInteract:CORRECTION DecIn2 error");
            }
#ifdef fdebug
            G4cout<<"G4QE::FSI:***CORRECTION IS DONE*** "<<G4endl;
#endif
            curHadr->Set4Momentum(cH4Mom);
            prevHadr->Set4Momentum(pH4Mom);
          }
          else
          {
            G4cerr<<"*!*G4QE::FSI:NoCor "<<cHM<<"+"<<pHM<<"="<<cHM+pHM<<">"<<totRM<<G4endl;
		    //#ifdef pdebug
            throw G4QException("***G4QEnvironment::FSInteraction: TEMPORARY EXCEPTION");
            //#endif
          }
        }
#ifdef fdebug
        else G4cout<<"G4QE::FSI: Yes, it is. d="<<dem<<G4endl;
#endif
      }
    } // End of the Back fusion LOOP
    // >| 2     | 2  | 2     | 2     | 2      | 2 - old    | 1. If gamma: add to sum4Mom
    //  |>0->sum| 3  | 3     | 3     | 3      | 3 - old    | 2. Compare BN with the Last
    //  | 5     |>5  | 4     | 4     | 4      | 4 - old    | 3. Swap if larger, del theLast
    //  | 0     | 0  |>0->sum| 5<-sum| 5->evap| 2 - new    | 4. If the Last: add the sum
    //  | 4     | 4  | 5     | ex    |        |(0 - gamma?)| 5. Decay/Eporate the Last
    //  | 3     | ex |                        | 3 - new
    G4LorentzVector sum(0.,0.,0.,0.);
    nHadr=theQHadrons.size();
#ifdef fdebug
    G4cout<<"G4QE::FSI: ===Before the gamma compression===, nH="<<nHadr<<G4endl;
#endif
    G4bool frag=false;
    if(nHadr>2)for(unsigned f=0; f<theQHadrons.size()-1; f++)//Check that there is a fragm
    {
      if(theQHadrons[f]->GetBaryonNumber()>1) frag=true; // The fragment is found
    }
	if(nHadr>2 && frag) for(G4int h=nHadr-1; h>=0; h--)//Collect gammas & kill DecayedHadrs
    {
      G4QHadron* curHadr = theQHadrons[h];      // Get a pointer to the current Hadron
      G4int   hF = curHadr->GetNFragments();
      G4int hPDG = curHadr->GetPDGCode();
      if(hPDG==89999003||hPDG==90002999)
        G4cerr<<"---Worning---G4QEnv::FSI:nD-/pD++(1)="<<hPDG<<G4endl;
#ifdef fdebug
	  G4cout<<"G4QE::FSI: h#"<<h<<", hPDG="<<hPDG<<", hNFrag="<<hF<<G4endl;
#endif
      //if(hPDG==22) gamCount++;
      //if(hF)       decCount++;
      //if(hF||hPDG==22&&gamCount>1)             // It should be compressed
	  if(hF||hPDG==22)                         // It should be compressed
	  {
        G4QHadron* theLast = theQHadrons[theQHadrons.size()-1];//Get Ptr to the Last Hadron
        if(hPDG==22) sum+=curHadr->Get4Momentum(); // Add 4Mom of gamma to the "sum"
        if(h<static_cast<G4int>(theQHadrons.size())-1) // Need swap with the Last
	    {
          curHadr->SetNFragments(0);
          curHadr->Set4Momentum(theLast->Get4Momentum());
          //curHadr->SetQC(theLast->GetQC());
          curHadr->SetQPDG(theLast->GetQPDG());
#ifdef fdebug
	      G4cout<<"G4QE::FSI: Exchange with the last is done"<<G4endl;
#endif
	    }
        theQHadrons.pop_back();                // theLastQHadron is excluded from QHadrons
        delete theLast;//!!When kill,DON'T forget to delete theLastQHadron as an instance!!
        nHadr--;
#ifdef fdebug
	    G4cout<<"G4QE::FSI: The last is compessed"<<G4endl;
#endif
      }
    }
#ifdef fdebug
    G4cout<<"G4QE::FSI: nH="<<nHadr<<"="<<theQHadrons.size()<<", sum="<<sum<<G4endl;
#endif
    if(nHadr>1)for(unsigned hdr=0; hdr<theQHadrons.size()-1; hdr++)//Ord:theBigestIstheLast
    {
#ifdef fdebug
      G4cout<<"G4QE::FSI:ORD,h="<<hdr<<"<"<<nHadr<<",hPDG="<<theQHadrons[hdr]->GetPDGCode()
            <<G4endl;
#endif
      G4QHadron* curHadr = theQHadrons[hdr];   // Get a pointer to the current Hadron
      G4QHadron* theLast = theQHadrons[theQHadrons.size()-1]; //Get Ptr to the Last Hadron
      G4int hB           = curHadr->GetBaryonNumber();
      G4int lB           = theLast->GetBaryonNumber();
#ifdef fdebug
      G4cout<<"G4QE::FSI:hBN="<<hB<<"<lBN="<<lB<<",lstPDG="<<theLast->GetPDGCode()<<G4endl;
#endif
      if(lB<hB)                                // Must be swapped
	  {
        G4QPDGCode   hQPDG = curHadr->GetQPDG();
        G4LorentzVector h4m= curHadr->Get4Momentum();
        curHadr->Set4Momentum(theLast->Get4Momentum());
        curHadr->SetQPDG(theLast->GetQPDG());
        theLast->Set4Momentum(h4m);
        theLast->SetQPDG(hQPDG);
	  }
    }
    nHadr=theQHadrons.size();
    G4QHadron* theLast = theQHadrons[nHadr-1]; // Get a pointer to the Last Hadron
    if(theLast->GetBaryonNumber()>1) // "Absorb photons & try to evaporate/decay" case
    {
      G4QHadron* theNew  = new G4QHadron(theLast); // Make a New Hadron of the Last Hadron
#ifdef fdebug
      G4cout<<"G4QE::FSI:BeforeLastSub nH="<<nHadr<<",lPDG="<<theNew->GetPDGCode()<<G4endl;
#endif
      theQHadrons.pop_back();                  // the last QHadron is excluded from OUTPUT
      delete theLast;//*!!When kill,DON'T forget to delete theLastQHadron as an instance!!*
      nHadr--;                                // TheLastHadron is only virtually exists now
      G4LorentzVector exRes4M=theNew->Get4Momentum()+sum;//Icrease 4Mom of theLast by "sum"
      G4QNucleus exResidN(exRes4M,theNew->GetPDGCode());
      //G4double mGamEva=2700.;               // Threshold for the evaporation
      G4double mGamEva=1700.;                 // Threshold for the evaporation
	  if(exResidN.SplitBaryon())
	  //if(2>3)                               //CloseTheFirstPriorityResN+gamSumEvaporation
	  {
        theNew->Set4Momentum(exRes4M);        // Icrease 4Mom of theLast by "sum" to Evapor
#ifdef fdebug
        G4cout<<"G4QE::FSI:Before E(1),nH="<<nHadr<<",nPDG="<<theNew->GetPDGCode()<<G4endl;
#endif
        EvaporateResidual(theNew);            // Try to evaporate the Nucl.(@@DecDib)(d.e.)
      }
      else if(theNew->GetPDGCode()==90002002 && exRes4M.m()>mHe3+mNeut
                                             && G4UniformRand()>.5)
	  {
        theNew->Set4Momentum(exRes4M);        // Icrease 4Mom of theLast by "sum" to Evapor
        G4LorentzVector n4M(0.,0.,0.,mNeut);
        G4LorentzVector h4M(0.,0.,0.,mHe3);
        if(!theNew->DecayIn2(n4M,h4M))
        {
          G4cerr<<"***G4QE::FSI:GamSup, tM="<<exRes4M.m()<<"<n+He3="<<mNeut+mHe3<<G4endl;
          throw G4QException("***G4QEnvironment::FSInter:GamSUPPRES DecIn2(n+He3) error");
        }
#ifdef fdebug
        G4cout<<"G4QE::FSI:Gamma Suppression succided, n="<<n4M<<", He3="<<h4M<<G4endl;
#endif
        theNew->Set4Momentum(n4M);
        theNew->SetQPDG(nQPDG);               // convert the alpha to the neutron
        theQHadrons.push_back(theNew);        // (delete equivalent for theHad=neutron)
        G4QHadron* theHe3 = new G4QHadron(90002001,h4M);// Make a New Hadr for the He3
        theQHadrons.push_back(theHe3);        // (delete equivalent for the proton)
      }
	  else if(nHadr)                           // Get LastHadrBefResNuc, absorb gam & decay
	  //else if(2>3)                           // Close the pair absorbtion of gamma
	  {
        if(nHadr>1)for(unsigned sh=0; sh<theQHadrons.size()-1; sh++)//Ord:theMinE isTheLast
        {
          G4QHadron* curHadr = theQHadrons[sh];// Get a pointer to the current Hadron
          G4QHadron* thePrev = theQHadrons[theQHadrons.size()-1];//Get Ptr to theLastHadron
          G4LorentzVector h4M= curHadr->Get4Momentum();
          G4LorentzVector l4M= thePrev->Get4Momentum();
#ifdef fdebug
          G4cout<<"G4QE::FSI:SO,h="<<sh<<"<"<<nHadr<<",PDG/LV="<<curHadr->GetPDGCode()<<h4M
                <<G4endl;
#endif
          G4double hM=h4M.m();
          G4double hT=h4M.e()-hM;
          G4double lT=l4M.e()-l4M.m();
#ifdef fdebug
          G4cout<<"G4QE::FSI:hT="<<hT<<"<lT="<<lT<<".lPDG="<<thePrev->GetPDGCode()<<G4endl;
#endif
          if(hM>mGamEva&&lT>hT)                // Must be swapped as the current is smaller
	      {
            G4QPDGCode   hQPDG = curHadr->GetQPDG();
            curHadr->Set4Momentum(l4M);
            curHadr->SetQPDG(thePrev->GetQPDG());
            thePrev->Set4Momentum(h4M);
            thePrev->SetQPDG(hQPDG);
	      }
        }
        nHadr=theQHadrons.size();
        G4QHadron* thePrev = theQHadrons[nHadr-1]; // Get Ptr to the BeforeResidNuclHadron
        if(thePrev->Get4Momentum().m()>mGamEva)
        {
          G4QHadron* theHad  = new G4QHadron(thePrev);// Make NewHadr of theBeforeResNuclH
#ifdef fdebug
          G4cout<<"G4QE::FSI:BeforeResidNucHadr nH="<<nHadr<<",hPDG="<<theHad->GetPDGCode()
                <<G4endl;
#endif
          theQHadrons.pop_back();                    // theLastQHadron excluded from OUTPUT
          delete thePrev;//*!When kill,DON'T forget to delete theLastQHadron asAnInstance!*
          G4LorentzVector n4M=theNew->Get4Momentum();// 4Mom of the Last (biggest nucleus)
          G4LorentzVector h4M=theHad->Get4Momentum();// 4Mom of the previous Hadron in HV
          G4LorentzVector dh4M=exRes4M+h4M;          // 4Mom of LH+PH+sum(gam) for theDecay
          G4double dhM=dh4M.m();                     // M of LH+PH+sum(gammas) for theDecay
          if(theHad->GetPDGCode()==90001001&&dhM>n4M.m()+mProt+mNeut && G4UniformRand()>.5)
		  //if(2>3)                                  // Close Possibility toSplitDeuteron
          {
            G4double nuM=n4M.m();
            h4M=G4LorentzVector(0.,0.,0.,mNeut);
            G4LorentzVector p4M(0.,0.,0.,mProt);
            G4double sum=nuM+mNeut+mProt;
            if(fabs(dhM-sum)<eps)
            {
              n4M=dh4M*(nuM/sum);
              h4M=dh4M*(mNeut/sum);
              p4M=dh4M*(mProt/sum);
            }
            else if(dhM<sum || !G4QHadron(dh4M).DecayIn3(n4M,h4M,p4M))
            {
              G4cerr<<"***G4QE::FSI:GamSupByD,M="<<dhM<<"<A+p+n="<<sum<<G4endl;
              throw G4QException("G4QEnviron::FSInt:Gamma SUPPRESSSION by D DecIn3 error");
            }
#ifdef fdebug
            G4cout<<"G4QE::FSI:GamSuppression by d succided, h="<<h4M<<", A="<<n4M<<G4endl;
#endif
            theHad->Set4Momentum(h4M);
            theHad->SetQPDG(nQPDG);               // convert the deuteron to the neutron
            theQHadrons.push_back(theHad);        // (delete equivalent for theHad=neutron)
            G4QHadron* theProt = new G4QHadron(90001000,p4M);// Make NewHadr for the proton
            theQHadrons.push_back(theProt);       // (delete equivalent for the proton)
            theNew->Set4Momentum(n4M);
            EvaporateResidual(theNew);    // Try to evaporate theResNuc once more (del.eq.)
          }
          else
          {
            if(!G4QHadron(dh4M).DecayIn2(n4M,h4M))
            {
              G4cerr<<"***G4QE::FSI:GamSup,M="<<dh4M.m()<<"<A+h="<<n4M.m()+h4M.m()<<G4endl;
              throw G4QException("G4QEnviron::FSInt:GamSUPPRES (3) DecIn2 didn't succeed");
            }
#ifdef fdebug
            G4cout<<"G4QE::FSI:Gamma Suppression succided, h="<<h4M<<", A="<<n4M<<G4endl;
#endif
            theHad->Set4Momentum(h4M);
            theQHadrons.push_back(theHad);          // (delete equivalent for theHad)
            theNew->Set4Momentum(n4M);
            EvaporateResidual(theNew);              // Try to evaporate theResNuc (del.eq.)
          }
	    }
        else
	    {
          theNew->Set4Momentum(exRes4M);// Icrease 4Mom of the Last by the "sum" for Evap
#ifdef fdebug
          G4cout<<"G4QE::FSI:BeforeE(2),nH="<<nHadr<<",PDG="<<theNew->GetPDGCode()<<G4endl;
#endif
          EvaporateResidual(theNew);    // Try to evaporate the Nucl.(@@DecDib)(delete eq.)
        }
      }
	  else                              // Absorb gammas to theResidNucleus and evaporateIt
	  {
        theNew->Set4Momentum(exRes4M);  // Icrease 4Mom of the Last by the "sum" for Evap
#ifdef fdebug
        G4cout<<"G4QE::FSI: BeforeE(3),nH="<<nHadr<<",nPDG="<<theNew->GetPDGCode()<<G4endl;
#endif
        EvaporateResidual(theNew);      // Try to evaporate the Nucl.(@@DecDib)(delete eq.)
      }
      //G4int onH=nHadr;
      nHadr=theQHadrons.size();
      //if(nHadr>onH) bfAct=true;
    } // End of "the last is the nucleus" case
  } // End of the While-LOOOP for the Back Fusion
  // Final attempt to alpha-decay the residual nucleus, suppressing the gamma ===========
  G4int gamcnt=0; // Counter of the residual gammas at this level
  unsigned maxB=nHadr-1;
  lHadr=theQHadrons[maxB]->GetBaryonNumber();
  if(nHadr>1)for(unsigned ipo=0; ipo<theQHadrons.size(); ipo++) // Find theBiggestNuclFragm
  {
    G4int hPDG = theQHadrons[ipo]->GetPDGCode();
    if(hPDG==22) gamcnt++;
    else
    {
      G4int hBN  = theQHadrons[ipo]->GetBaryonNumber();
#ifdef fdebug
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
#ifdef fdebug
  G4cout<<"G4QE::FSI:maxB#"<<maxB<<", gamcnt="<<gamcnt<<G4endl;
#endif
  if(gamcnt)                                    // Only if there are gammas one should act
  {
    if(maxB+1<nHadr)                            // If maxB<Last, swap theCurH and theLastH
    {
      G4QHadron* theCurr = theQHadrons[maxB];   // Pointer to the Current Hadron
      G4QHadron* theLast = theQHadrons[nHadr-1];// Pointer to the Last Hadron
      G4QHadron* curHadr = new G4QHadron(theCurr);//Remember theCurrentHadron to put on top
      theCurr->SetQPDG(theLast->GetQPDG());     // theCurHadr is substituted by theLastHadr
      theCurr->Set4Momentum(theLast->Get4Momentum()); // ... continue substitution
      theQHadrons.pop_back();                   // Rnt to theLastHadron is excluded from HV
      delete theLast;//*!!When kill,DON'T forget to delete the last QHadron as an inst. !!*
      theQHadrons.push_back(curHadr);           // The current Hadron, which is the Biggest
    }
    nHadr=theQHadrons.size();                   // Must be the same
    // Now it is necessary to absorb the photon (photons) and try to radiate alpha or evap.
    G4LorentzVector gamSum(0.,0.,0.,0.);
    if(nHadr>1)for(unsigned gp=0; gp<nHadr-1; gp++)// Find Gumma, remember and kill
    {
      G4QHadron* theCurr = theQHadrons[gp];       // Pointer to the Current Hadron
      G4int hPDG=theCurr->GetPDGCode();
#ifdef fdebug
      G4cout<<"G4QE::FSI:gp#"<<gp<<", PDG="<<hPDG<<", is found"<<G4endl;
#endif
      if(hPDG==22)
	  {
        gamSum=gamSum+theCurr->Get4Momentum();    // Accumulate the 4Momenta of photons
        G4QHadron* theLast = theQHadrons[nHadr-1];// Pointer to the Last Hadron
        G4QPDGCode theLQPDG=theLast->GetQPDG();
        theCurr->SetQPDG(theLQPDG);               // theCurHad is substituted by theLastHad
        theCurr->Set4Momentum(theLast->Get4Momentum()); // ... continue substitution
        theQHadrons.pop_back();                   // Pnt to theLastHadr.is excluded from HV
        delete theLast;//*!!When kill,DON'T forget to delete theLastQHadron as an inst.!!*
        nHadr=theQHadrons.size();
#ifdef fdebug
        G4cout<<"G4QE::FSI:Rep. by lPDG="<<theLQPDG<<",nH="<<nHadr<<",gS="<<gamSum<<G4endl;
#endif
	  }
    }
    // @@ Now it is necessary to try to emit alpha or evaporate the residual nucleus
    G4QHadron* theLast = theQHadrons[nHadr-1];   // Pointer to the Last Hadron
    G4QHadron* curHadr = new G4QHadron(theLast); // Pnt to theCurrentHadron is theLastCopy
    theQHadrons.pop_back();                 // Pnt to theLastHadron is excluded from OUTPUT
    delete theLast;//*!! When kill,DON'T forget to delete theLas QHadron as an instance !!*
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
#ifdef fdebug
          G4cout<<"G4QE::FSI: Final A+alpha, A="<<APDG<<lr4M<<", a="<<al4M<<G4endl;
#endif
          G4QHadron* lrH = new G4QHadron(APDG,lr4M);
          theQHadrons.push_back(lrH);      // (delete equivalent for newH)
          G4QHadron* alH = new G4QHadron(90002002,al4M);
          theQHadrons.push_back(alH);      // (delete equivalent for newH)
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
      G4cout<<"G4QE::FSI: After Evap (3) nH="<<theQHadrons.size()<<G4endl;
#endif
	}
  }
  //Now just fill the output theFravment vector (User is responsible to ClearAndDestroy it)
  nHadr=theQHadrons.size();
  if(nHadr) for(unsigned hd=0; hd<theQHadrons.size(); hd++)
  {
    G4QHadron* curHadr = new G4QHadron(theQHadrons[hd]);
    G4int hPDG=curHadr->GetPDGCode();
    if(hPDG==22 && fabs(curHadr->Get4Momentum().e())<.00001) // E=0, gamma in the OUTPUT
    {
      unsigned lin=theQHadrons.size()-1;
      if(lin>hd)
      {
        G4QHadron* theLast = theQHadrons[lin];// Pointer to theLastHadron in theQHadrVector
        curHadr->SetQPDG(theLast->GetQPDG());//theCurHadron is substituted by theLastHadron
        curHadr->Set4Momentum(theLast->Get4Momentum()); // ... continue substitution (4Mom)
        theQHadrons.pop_back();              //pointer to theLastHadron is excluded from HV
        delete theLast;//*!!When kill, DON'T forget to delete theLastQHadron asAnInstance!!
        hPDG=curHadr->GetPDGCode();          // Redefine thePDGCode of theCurrentHadron
      }
      else break; // It's not necessary to delete(~): do not copy to theFragments is enough
    }
#ifdef fdebug
    G4cout<<"G4QE::FSI: LOOP starts nH="<<nHadr<<", h#"<<hd<<", PDG="<<hPDG<<G4endl;
#endif
    if(hPDG==89999003||hPDG==90002999) G4cout<<"G4QEnv::FSI:nD-/pD++(0)="<<hPDG<<G4endl;
    if(hPDG==89999004||hPDG==90003999) G4cout<<"G4QEnv::FSI:nnD-/ppD++(0)="<<hPDG<<G4endl;
#ifdef pdebug
    G4cout<<"G4QE::FSI:Copy is made with PDG="<<hPDG<<G4endl;//To compare with Original
#endif
    G4int hBN=curHadr->GetBaryonNumber();
    G4int hCG=curHadr->GetCharge();
    G4int hST=curHadr->GetStrangeness();
    if(hPDG==91000000) curHadr->SetQPDG(G4QPDGCode(3122));// Move to theNextLoop @@NucToHad
    else if(hBN>1 && (hBN==hCG&&!hST || !hCG&&!hST || !hCG&&hST==hBN))
	{
      if     (!hCG&&!hST)     hPDG=90000001;
      else if(hBN==hCG&&!hST) hPDG=90001000;
      else if(!hCG&&hST==hBN) hPDG=91000000;
      else throw G4QException("***G4QEnvironment::FSInteract: MultyDibaryon cant be here");
      G4LorentzVector newLV=(curHadr->Get4Momentum())/hBN;
      curHadr->Set4Momentum(newLV);
      curHadr->SetQPDG(G4QPDGCode(hPDG));
      for(G4int bd=1; bd<hBN; bd++)
      {
        G4QHadron* secHadr = new G4QHadron(curHadr);
        theFragments->push_back(secHadr);// (del.equivalent - user is responsible for that)
      }
    }
    else if(curHadr->GetStrangeness()<0&&curHadr->GetBaryonNumber()>0)// AntistrangeNucleus
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
      G4bool fOK=true;
      G4LorentzVector n4M(0.,0.,0.,mR);
      G4LorentzVector k4M(0.,0.,0.,mKaon);
      G4double sum=mR+mKaon;
      if(fabs(reM-sum)<eps)
	  {
        n4M=r4M*(mR/sum);
        k4M=r4M*(mKaon/sum);
      }
      else if(reM<sum || !G4QHadron(r4M).DecayIn2(n4M,k4M))
      {
        G4cerr<<"---Worning---G4QE::FSI: Try to recover ASN="<<hPDG<<","<<reM<<"<"<<mR<<"+"
              <<mKaon<<"="<<sum<<G4endl;
        if(!theEnvironment.GetA())
        {
          G4QHadron* theLast = curHadr;          // Prototype of Pointer to theLastHadron
          G4QHadron* qH = new G4QHadron(curHadr);// Copy of the Current Hadron
          if(hd+1<theQHadrons.size())           // If ipo<Last, swap CurHadr & LastHadr
          {
            theLast = theQHadrons[theQHadrons.size()-1]; // Pointer to LastHadr(ipo<Last)
            curHadr->SetQPDG(theLast->GetQPDG());// CurHadrPDG's substituted by LastHaPDG
            curHadr->Set4Momentum(theLast->Get4Momentum()); // ... 4Momentum substitution
          }
          theQHadrons.pop_back();        // exclude theLastHadron pointer from the OUTPUT
          delete theLast; //*!! When killing, DON'T forget to delete the last QHadron !!*
          G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());//Fake Quasmon
          if(!CheckGroundState(quasH,true))      // Try to correct by other hadrons
          {
            G4cerr<<"---Worning---G4QE::FSI:AntiStraN Failed LeaveAsItIs 4m="<<r4M<<G4endl;
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
          G4cerr<<"***G4QEnv::FSI: No Final AntiSN recovery, E="<<theEnvironment<<G4endl;
          throw G4QException("***G4QEnvironment::FSInter:AntistrangeNucleus decayError");
        }
        if(fOK)
        {
#ifdef fdebug
          G4cout<<"G4QE::FSI:AntiSN==>"<<r4M<<"->A="<<RPDG<<n4M<<"+K="<<kPDG<<k4M<<G4endl;
#endif
          curHadr->Set4Momentum(n4M);
          curHadr->SetQPDG(G4QPDGCode(RPDG));  // converted antiStrangeNucleus
          G4QHadron* theKaon = new G4QHadron(kPDG,k4M);// Make a New Hadr for the kaon
          theFragments->push_back(theKaon);    // (delete equivalent for the kaon)
        }
        hPDG = curHadr->GetPDGCode();          // Change PDG Code of theFirstFragment ?
      }
    }
    else if(hPDG>90500000)            // Hypernucleus (@@ Only for the GEANT4 applications)
	{
#ifdef fdebug
      G4cerr<<"**G4QEnvironment::FSI: Hypernucleus PDG="<<hPDG<<" must decay"<<G4endl;
#endif
      G4int nL=curHadr->GetStrangeness();      // A number of Lambdas in the Hypernucleus
      G4LorentzVector r4M=curHadr->Get4Momentum(); // Real 4-momentum of the hypernucleus
      G4double reM=r4M.m();                    // Real mass of the hypernucleus
      G4int rnPDG = hPDG-nL*999999;            // Convert Lambdas to neutrons
      G4double rnM=G4QNucleus(rnPDG).GetMZNS();// Mass of the new nonstrange nucleus
      // @@ Take into account that there can be Iso-Hypernucleus (Add PI+,R & Pi-,R decays)
      if(reM<rnM+mPi0)
	  {
        G4double d=rnM+mPi0-reM;
		G4cerr<<"***G4QE::FSI:Small Hypernuclear M="<<reM<<">"<<rnM+mPi0<<",d="<<d<<G4endl;
        throw G4QException("G4QEnvironment::FSInteract: Hypernuclear conversion");
      }
      else
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
          G4cerr<<"***G4QE::FSI:Hypern, M="<<reM<<"<A+n*Pi0="<<sum<<",d="<<sum-reM<<G4endl;
          throw G4QException("***G4QEnvironment::FSInter: Hypernuclear decay error");
        }
#ifdef fdebug
        G4cout<<"*G4QE::FSI:HypN "<<r4M<<"->A="<<rnPDG<<n4M<<",n*Pi0="<<nPi<<h4M<<G4endl;
#endif
        curHadr->Set4Momentum(n4M);
        curHadr->SetQPDG(G4QPDGCode(rnPDG));  // converted hypernucleus
        if(nPi>1) h4M=h4M/nPi;                // split the 4-mom if necessary
        for(G4int ihn=0; ihn<nPi; ihn++)      // loop over additional pions
        {
          G4QHadron* thePion = new G4QHadron(111,h4M);// Make a New Hadr for the pion
          theFragments->push_back(thePion);   // (delete equivalent for the pion)
        }
        if(rnPDG==90000002)                   // Additional action with curHadr change
	    {
          G4LorentzVector newLV=n4M/2.;
          curHadr->Set4Momentum(newLV);
          curHadr->SetQPDG(G4QPDGCode(90000001));
          G4QHadron* secHadr = new G4QHadron(curHadr);
          theFragments->push_back(secHadr);   // (del.equiv.- user is responsible for del)
        }
        else if(rnPDG==90002000)              // Additional action with curHadr change
	    {
          G4LorentzVector newLV=n4M/2.;
          curHadr->Set4Momentum(newLV);
          curHadr->SetQPDG(G4QPDGCode(90001000));
          G4QHadron* secHadr = new G4QHadron(curHadr);
          theFragments->push_back(secHadr);   // (del.equiv.- user is responsible for del)
        }
      }
    }
    theFragments->push_back(curHadr);  // (delete equivalent - user is responsible for del)
  } // 
#ifdef fdebug
  G4cout<<"G4QE::FSI:==>OUT,nH="<<theQHadrons.size()<<",nF="<<theFragments->size()<<G4endl;
#endif
  return theFragments;
} // End of "FSInteraction"

//The public Quasmons duplication with delete responsibility of User (!)
G4QuasmonVector* G4QEnvironment::GetQuasmons()
{//              =============================
  G4int nQ=theQuasmons.size();
#ifdef pdebug
  G4cout<<"G4QEnvironment::GetQuasmons is called nQ="<<nQ<<G4endl;
#endif
  G4QuasmonVector* quasmons = new G4QuasmonVector; // Intermediate
  if(nQ) for(G4int iq=0; iq<nQ; iq++)
  {
#ifdef pdebug
    G4cout<<"G4QEnv::GetQuasm:Q#"<<iq<<",QQPDG="<<theQuasmons[iq]->GetQPDG()<<",QQC="
          <<theQuasmons[iq]->GetQC()<<",M="<<theQuasmons[iq]->Get4Momentum().m()<<G4endl;
#endif
    G4Quasmon* curQ = new G4Quasmon(theQuasmons[iq]);
    quasmons->push_back(curQ);                 // (delete equivalent - user is responsible)
  }
#ifdef pdebug
  G4cout<<"G4QEnvironment::GetQuasmons ===OUT==="<<G4endl;
#endif
  return quasmons;
} // End of GetQuasmons

//The public Quasmons duplication with delete responsibility of User (!)
G4QHadronVector* G4QEnvironment::GetQHadrons()
{//              =============================
  G4int nH=theQHadrons.size();
#ifdef pdebug
  G4cout<<"G4QEnvironment::GetQHadrons is called nH="<<nH<<G4endl;
#endif
  G4QHadronVector* hadrons = new G4QHadronVector; // Intermediate
  if(nH) for(G4int ih=0; ih<nH; ih++)
  {
#ifdef pdebug
    G4cout<<"G4QEnv::GetQHadrons:H#"<<ih<<",HQPDG="<<theQHadrons[ih]->GetQPDG()<<",HQC="
          <<theQHadrons[ih]->GetQC()<<",HM="<<theQHadrons[ih]->GetMass()<<G4endl;
#endif
    G4QHadron* curH = new G4QHadron(theQHadrons[ih]);
    hadrons->push_back(curH);                       // (del. equiv. - user is responsibile)
  }
#ifdef pdebug
  G4cout<<"G4QEnvironment::GetQHadrons ===OUT==="<<G4endl;
#endif
  return hadrons;
} // End of GetQHadrons

//Unavoidable decay of the Isonucleus in nP+(Pi+) or nN+(Pi-)
void G4QEnvironment::DecayIsonucleus(G4QHadron* qH)
{//  ============================================
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double eps  = 0.003;
  G4LorentzVector q4M = qH->Get4Momentum();      // Get 4-momentum of the Isonucleus
  G4double qM=q4M.m();                           // Real mass of the Isonucleus
  G4QContent qQC = qH->GetQC();                  // Get QuarcContent of the Isonucleus
  G4int qBN=qQC.GetBaryonNumber();               // Baryon number of the IsoNucleus
  G4int qC=qQC.GetCharge();                      // Charge of the IsoNucleus
  G4int qPN=qC-qBN;                              // Number of pions in the Isonucleus
  G4int          fPDG = 2212;                    // Prototype for nP+(Pi+) case
  G4int          sPDG = 211;
  G4double       fMass= mProt;
  G4double       sMass= mPi;
  if(qC<0)
  {
    qPN  = -qC;
    fPDG = 2112;                                 // nN+(Pi-) case
    sPDG = -211;
    fMass= mNeut;
  }
  G4double tfM=qBN*fMass;
  G4double tsM=qPN*sMass;
  G4LorentzVector f4Mom(0.,0.,0.,tfM);
  G4LorentzVector s4Mom(0.,0.,0.,tsM);
  G4double sum=tfM+tsM;
  if(fabs(qM-sum)<eps)
  {
    f4Mom=q4M*(tfM/sum);
    s4Mom=q4M*(tsM/sum);
  }
  else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
  {
#ifdef pdebug
    G4cerr<<"***G4QEnv::DecIsonuc:fPDG="<<fPDG<<"*"<<qBN<<"(fM="<<fMass<<") + sPDG="<<sPDG
          <<"*"<<qPN<<"(sM="<<sMass<<")"<<"="<<sum<<" > TotM="<<qM<<q4M<<qQC<<G4endl;
#endif
    if(!theEnvironment.GetA())
    {
      G4Quasmon* quasH = new G4Quasmon(qQC,q4M);
	  if(!CheckGroundState(quasH,true)) theQHadrons.push_back(qH);// Cor or fill as it is
      else delete qH;  
      delete quasH;
      return;
    }
    else
    {
      G4cerr<<"***G4QE::DecayIsonucleus:Cann't Correct*EmptyEnv*="<<theEnvironment<<G4endl;
      throw G4QException("G4QEnvir::DecayIsonucleus: Isonucleus DecIn2 didn't succeed");
    }
  }
#ifdef pdebug
  G4cout<<"G4QEnv::DecayIsonucleus: *DONE* sum4M="<<f4Mom<<",nPDG="<<fPDG
        <<", m4M="<<s4Mom<<",mPDG="<<sPDG<<G4endl;
#endif
  delete qH;
  f4Mom/=qBN;
  for(G4int ih=0; ih<qBN; ih++)
  {
    G4QHadron* Hi = new G4QHadron(fPDG,f4Mom);   // Create a Hadron for the baryon
    theQHadrons.push_back(Hi);                   // Fill "Hi" (delete equivalent)
  }
  s4Mom/=qPN;
  for(G4int ip=0; ip<qPN; ip++)
  {
    G4QHadron* Hj = new G4QHadron(sPDG,s4Mom);     // Create a Hadron for the meson
    theQHadrons.push_back(Hj);                     // Fill "Hj" (delete equivalent)
  }
} // End of DecayIsonucleus

//Decay of the excited dibayon in two baryons
void G4QEnvironment::DecayDibaryon(G4QHadron* qH)
{//  ============================================
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mSigM= G4QPDGCode(3112).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mSigP= G4QPDGCode(3222).GetMass();
  static const G4double mKsiM= G4QPDGCode(3312).GetMass();
  static const G4double mKsiZ= G4QPDGCode(3322).GetMass();
  static const G4double mDeut= G4QPDGCode(2112).GetNuclMass(1,1,0);
  static const G4double mPiN = mPi+mNeut;
  static const G4double mPiP = mPi+mProt;
  static const G4double dmPiN= mPiN+mPiN;
  static const G4double dmPiP= mPiP+mPiP;
  static const G4double nnPi = mNeut+mPiN;
  static const G4double ppPi = mProt+mPiP;
  static const G4double lnPi = mLamb+mPiN;
  static const G4double lpPi = mLamb+mPiP;
  static const G4double dNeut= mNeut+mNeut;
  static const G4double dProt= mProt+mProt;
  static const G4double dLamb= mLamb+mLamb;
  static const G4double dLaNe= mLamb+mNeut;
  static const G4double dLaPr= mLamb+mProt;
  static const G4double dSiPr= mSigP+mProt;
  static const G4double dSiNe= mSigM+mNeut;
  static const G4double dKsPr= mKsiZ+mProt;
  static const G4double dKsNe= mKsiM+mNeut;
  static const G4double eps  = 0.003;
  static const G4QNucleus vacuum(90000000);
#ifdef ppdebug
  G4cerr<<"***G4QEnvironment::DecayDibaryon: *** TEST OF EXCEPTION ***"<<G4endl;
  throw G4QException("G4QEnv::DecayDibaryon:Unknown DBary PDG code or small Mass of DB");
#endif
  G4bool four=false;                           // defFALSE for 4-particle decay of diDelta
  G4LorentzVector q4M = qH->Get4Momentum();    // Get 4-momentum of the Dibaryon
  G4int          qPDG = qH->GetPDGCode();      // PDG Code of the decayin dybaryon
  G4double         qM = q4M.m();
  G4double         rM = qM+eps;                // Just to avoid the computer accuracy
#ifdef pdebug
  G4cout<<"G4QEnv::DecayDibaryon: *Called* PDG="<<qPDG<<",4Mom="<<q4M<<", M="<<qM<<G4endl;
#endif
  // Select a chanel of the dibaryon decay (including Delta+Delta-> 4 particle decay
  G4int          fPDG = 2212;                  // Prototype for pp case
  G4int          sPDG = 2212;
  G4int          tPDG = 0;                     // Zero prototype to separate 3 from 2 
  G4double       fMass= mProt;
  G4double       sMass= mProt;
  G4double       tMass= mPi;
  if     (qPDG==90003998 && rM>=dmPiP)         // "diDelta++" case
  {
    sPDG = 211;
    sMass= mPi;
    four = true;
  }
  else if(qPDG==89998004 && rM>=dmPiN)         // "diDelta--" case
  {
    sPDG = -211;
    fPDG = 2112;
    sMass= mPi;
    fMass= mNeut;
    four = true;
  }
  else if(qPDG==90000002 && rM>=dNeut)         // "dineutron" case
  {
    fPDG = 2112;
    sPDG = 2112;
    fMass= mNeut;
    sMass= mNeut;    
  }
  else if(qPDG==90001001 && rM>=mDeut)         // "exited deutron" case
  {
    if(abs(qM-mDeut)<eps)
	{
      theQHadrons.push_back(qH);               // Fill as it is (delete equivalent)
      return;
	}
    else if(mProt+mNeut<rM)
    {
      fPDG = 2112;
      fMass= mNeut;    
    }
    else
    {
      fPDG = 22;
      sPDG = 90001001;
      fMass= 0.;
      sMass= mDeut;    
    }
  }
  else if(qPDG==91000001 && rM>=dLaNe)         // "Lambda-neutron" case
  {
    fPDG = 2112;
    sPDG = 3122;
    fMass= mNeut;
    sMass= mLamb;    
  }
  else if(qPDG==91001000 && rM>=dLaPr)         // "Lambda-proton" case
  {
    sPDG = 3122;
    sMass= mLamb;    
  }
  else if(qPDG==89999003 && rM>=nnPi)         // "neutron/Delta-" case
  {
    fPDG = 2112;
    sPDG = 2112;
    tPDG = -211;
    fMass= mNeut;
    sMass= mNeut;    
  }
  else if(qPDG==90002999 && rM>=ppPi)         // "proton/Delta++" case
  {
    tPDG = 211;
  }
  else if(qPDG==90999002 && rM>=lnPi)         // "lambda/Delta-" case
  {
    fPDG = 2112;
    sPDG = 3122;
    tPDG = -211;
    fMass= mNeut;
    sMass= mLamb;    
  }
  else if(qPDG==91001999 && rM>=lpPi)         // "lambda/Delta+" case
  {
    sPDG = 3122;
    tPDG = 211;
    sMass= mLamb;    
  }
  else if(qPDG==90999002 && rM>=dSiNe)         // "Sigma-/neutron" case
  {
    fPDG = 2112;
    sPDG = 3112;
    fMass= mNeut;
    sMass= mSigM;    
  }
  else if(qPDG==91001999 && rM>=dSiPr)         // "Sigma+/proton" case
  {
    sPDG = 3222;
    sMass= mSigP;    
  }
  else if(qPDG==92000000 && rM>=dLamb)         // "diLambda" case
  {
    fPDG = 3122;
    sPDG = 3122;
    fMass= mLamb;
    sMass= mLamb;    
  }
  else if(qPDG==91999001 && rM>=dKsNe)         // "Ksi-/neutron" case
  {
    fPDG = 2112;
    sPDG = 3312;
    fMass= mNeut;
    sMass= mKsiM;    
  }
  else if(qPDG==92000999 && rM>=dKsPr)         // "Ksi0/proton" case
  {
    sPDG = 3322;
    sMass= mKsiZ;    
  }
  else if(qPDG!=90002000|| rM<dProt)           // Other possibilities (if not a default)
  {
    G4int qS = qH->GetStrangeness();
    G4int qB = qH->GetBaryonNumber();
    if(qB>0&&qS<0)                             // Antistrange diBarion
    {
      DecayAntiStrange(qH);
      return;
    }
    else
    {
#ifdef pdebug
      G4cerr<<"***G4QE::DecayDibaryon:PDG="<<qPDG<<",QC="<<qH->GetQC()<<",M="<<qM<<",Env="
            <<theEnvironment<<"nQ="<<theQuasmons.size()<<G4endl;
#endif
      if(!theEnvironment.GetA())
      {
        G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());
	    if(!CheckGroundState(quasH,true))
        {
#ifdef pdebug
          G4cerr<<"***G4QE::DecDiB:GSCorFailed. FillAsItIs,n="<<theQHadrons.size()<<G4endl;
#endif
          theQHadrons.push_back(qH);        // Correction failed: fill as it is
        }
        else delete qH;  
        delete quasH;
        return;
      }
      else
      {
        G4cerr<<"***G4QE::DecayDiBaryon:Cann't Correct*EmptyEnv*="<<theEnvironment<<G4endl;
        throw G4QException("G4QEnvir::DecayDibar:Unknown PDG code or small Mass of DB");
      }
    }
  }
  G4LorentzVector f4Mom(0.,0.,0.,fMass);
  G4LorentzVector s4Mom(0.,0.,0.,sMass);
  G4LorentzVector t4Mom(0.,0.,0.,tMass);
  if(!tPDG&&!four)
  {
    G4double sum=fMass+sMass;
    if(fabs(qM-sum)<eps)
	{
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cerr<<"---Worning---G4QE::DecDib:fPDG="<<fPDG<<"(M="<<fMass<<")+sPDG="<<sPDG<<"(M="
            <<sMass<<")="<<sum<<" >? TotM="<<q4M.m()<<q4M<<",Env="<<theEnvironment<<G4endl;
      if(!theEnvironment.GetA())
      {
        G4QContent qQC=qH->GetQC();
        G4Quasmon* quasH = new G4Quasmon(qQC,q4M);
	    if(!CheckGroundState(quasH,true))
        {
          G4cerr<<"---Worning---G4QE::DecDiBaryon: FillAsItIs diBaryon="<<qQC<<q4M<<G4endl;
          theQHadrons.push_back(qH);      // Fill as it is
        }
        else delete qH;  
        delete quasH;
        return;
      }
      else
      {
        G4cerr<<"***G4QE::DecayDiBar:Cann't Correct*EmptyEnv*="<<theEnvironment<<G4endl;
        throw G4QException("***G4QEnv::DecayDibaryon: DiBar DecayIn2 error");
      }
    }
#ifdef pdebug
    G4cout<<"G4QEnv::DecayDibaryon:(2) *DONE* f4M="<<f4Mom<<",fPDG="<<fPDG
          <<", s4M="<<s4Mom<<",sPDG="<<sPDG<<G4endl;
#endif
    delete qH;
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
    theQHadrons.push_back(H1);                 // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the 2-nd baryon
    theQHadrons.push_back(H2);                 // Fill "H2" (delete equivalent)
  }
  else if(four)
  {
    q4M=q4M/2.;                                // Divided in 2 !!!
    qM/=2.;                                    // Divide the mass in 2 !
    G4double sum=fMass+sMass;
    if(fabs(qM-sum)<eps)
	{
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cerr<<"---Worning---G4QE::DecDib:fPDG="<<fPDG<<"(M="<<fMass<<")+sPDG="<<sPDG<<"(M="
            <<sMass<<")"<<"="<<sum<<">tM="<<q4M.m()<<q4M<<",Env="<<theEnvironment<<G4endl;
      if(!theEnvironment.GetA())
      {
        // do not chande to q4M as it is devided by two !!
        G4QContent qQC=qH->GetQC();
        G4Quasmon* quasH = new G4Quasmon(qQC,qH->Get4Momentum());
	    if(!CheckGroundState(quasH,true))
        {
          G4cerr<<"---Worning---G4QE::DecDibar:FillAsItIs DelDel QC/4M="<<qQC<<q4M<<G4endl;
          theQHadrons.push_back(qH); // Fill as it is
        }
        else delete qH;  
        delete quasH;
        return;
      }
      else
      {
        G4cerr<<"***G4QE::DecayDibaryon:Cann't Correct*EmptyEnv*="<<theEnvironment<<G4endl;
        throw G4QException("***G4QEnv::DecDibaryon: Dibaryon DecayIn2 error");
      }
    }
#ifdef pdebug
    G4cout<<"G4QEnv::DecayDibaryon:(3) *DONE* f4M="<<f4Mom<<",fPDG="<<fPDG
          <<", s4M="<<s4Mom<<",sPDG="<<sPDG<<G4endl;
#endif
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
    theQHadrons.push_back(H1);                 // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the 2-nd baryon
    theQHadrons.push_back(H2);                 // Fill "H2" (delete equivalent)
    // Now the second pair mus be decayed
    if(fabs(qM-sum)<eps)
	{
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
    }
    else if(!G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      // Should not be here as sum was already compared with qM above for the first delta
      G4cerr<<"***G4QE::DecDibar:fPDG="<<fPDG<<"(fM="<<fMass<<") + sPDG="<<sPDG<<"(sM="
            <<sMass<<")="<<sum<<" >? (DD2,Can't be here) TotM="<<q4M.m()<<q4M<<G4endl;
      throw G4QException("G4QEnv::DecayDibaryon: General DecayIn2 error");
    }
#ifdef pdebug
    G4cout<<"G4QEnv::DecayDibaryon:(4) *DONE* f4M="<<f4Mom<<",fPDG="<<fPDG
          <<", s4M="<<s4Mom<<",sPDG="<<sPDG<<G4endl;
#endif
    G4QHadron* H3 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
    theQHadrons.push_back(H3);                 // Fill "H1" (delete equivalent)
    G4QHadron* H4 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the 2-nd baryon
    theQHadrons.push_back(H4);                 // Fill "H2" (delete equivalent)
    delete qH;
  }
  else
  {
    G4double sum=fMass+sMass+tMass;
    if(fabs(qM-sum)<eps)
	{
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
      t4Mom=q4M*(tMass/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn3(f4Mom, s4Mom, t4Mom))
    {
      G4cerr<<"---Worning---G4QE::DecDib:fPDG="<<fPDG<<"(M="<<fMass<<")+sPDG="<<sPDG<<"(M="
            <<sMass<<")+tPDG="<<tPDG<<"(tM="<<tMass<<")="<<sum<<">TotM="<<q4M.m()<<G4endl;
      //if(theEnvironment==vacuum)
      if(!theEnvironment.GetA())
      {
        G4Quasmon* quasH = new G4Quasmon(qH->GetQC(),qH->Get4Momentum());
	    if(!CheckGroundState(quasH,true)) theQHadrons.push_back(qH); // Cor or fill asItIs
        else delete qH;  
        delete quasH;
        return;
      }
      else
      {
        G4cerr<<"***G4QEnv::DecayDibaryon:Cann't Corr.*EmptyEnv*="<<theEnvironment<<G4endl;
        throw G4QException("G4QEnv::DecayDibaryon: diBar DecayIn3 error");
      }
    }
#ifdef pdebug
    G4cout<<"G4QEnv::DecayDibaryon:(5) *DONE* f4M="<<f4Mom<<",fPDG="<<fPDG<<", s4M="<<s4Mom
          <<",sPDG="<<sPDG<<", t4M="<<t4Mom<<",tPDG="<<tPDG<<G4endl;
#endif
    //qH->SetNFragments(2);                    // Fill a#of fragments to decaying Dibaryon
    //theQHadrons.push_back(qH);               // Fill hadron with nf=2 (delete equivalent)
    // Instead
    delete qH;
    //
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the 1-st baryon
    theQHadrons.push_back(H1);                 // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(sPDG,s4Mom); // Create a Hadron for the 2-nd baryon
    theQHadrons.push_back(H2);                 // Fill "H2" (delete equivalent)
    G4QHadron* H3 = new G4QHadron(tPDG,t4Mom); // Create a Hadron for the meson
    theQHadrons.push_back(H3);                 // Fill "H3" (delete equivalent)
  }
} // End of DecayDibaryon

//Decay of the nuclear states with antistrangeness (K:/K0)
void G4QEnvironment::DecayAntiStrange(G4QHadron* qH)
{//  ===============================================
  static const G4double mK    = G4QPDGCode(321).GetMass();
  static const G4double mK0   = G4QPDGCode(311).GetMass();
  static const G4double eps   = 0.003;
  G4LorentzVector q4M = qH->Get4Momentum();    // Get 4-mom of the AntiStrangeNuclearState
  G4double         qM = q4M.m();               // Real mass of the AntiStrangeNuclearState
  G4QContent       qQC= qH->GetQC();           // PDGCode of theDecayingAntiStrangeNuclSt.
  G4int            qS = qH->GetStrangeness();  // Strangness of the AntiStrangeNuclearState
  G4int            qB = qH->GetBaryonNumber(); // BaryonNumber of the AntiStrangeNuclearSt.
  G4int            qP = qH->GetCharge();       // Charge of the AntiStranNuclState (a#of p)
#ifdef pdebug
  G4cout<<"G4QEnv::DecAntS:QC="<<qQC<<",S="<<qS<<",B="<<qB<<",C="<<qP<<",4M="<<q4M<<G4endl;
#endif
  G4int            qN = qB-qP-qS;              // a#of neuterons
  if(qS>=0 || qB<1)
  {
    G4cerr<<"G4QEnv::DecayAntiStrange:QC="<<qQC<<",S="<<qS<<",B="<<qB<<",4M="<<q4M<<G4endl;
    throw G4QException("G4QEnvironment::DecayAntiStrange: not an Anti Strange Nucleus");
  }
  G4int n1=1;         // prototype of a#of K0's
  G4double k1M=mK0;
  G4int k1PDG=311;    // K0 (as a prototype)
  G4int n2=0;         // prototype of a#of K+'s
  G4double k2M=mK;
  G4int k2PDG=321;    // K+
  G4int aS=-qS;       // -Strangness = antistrangeness
  G4int sH=aS/2;     // a small half of the antistrangeness
  G4int bH=aS-sH;    // a big half to take into account all the antistrangeness
  if(qP>0 && qP>qN)  // a#of protons > a#of neutrons
  {
    if(qP>=bH)                       // => "Enough protons in nucleus" case
    {
      if(qN>=sH)
      {
        n1=sH;
        n2=bH;
      }
      else
      {
        G4int dN=qP-qN;
        if(dN>=aS)
        {
          n1=0;
          n2=aS;
        }
        else
        {
          G4int sS=(aS-dN)/2;
          G4int bS=aS-dN-sS;
          sS+=dN;
          if(qP>=sS&&qN>=bS)
          {
            n1=bS;
            n2=sS;
          }
          else if(qP<sS)
          {
            G4int dS=aS-qP;
            n1=dS;
            n2=qP;
          }
          else
          {
            G4int dS=aS-qN;
            n1=qN;
            n2=dS;
          }
        }
      }
    }
  }
  else if(qN>=bH)
  {
    if(qP>=sH)
    {
      n2=sH;
      n1=bH;
    }
    else
    {
      G4int dN=qN-qP;
      if(dN>=aS)
      {
        n1=aS;
        n2=0;
      }
      else
      {
        G4int sS=(aS-dN)/2;
        G4int bS=aS-sS;
        if(qN>=bS&&qP>=sS)
        {
          n1=bS;
          n2=sS;
        }
        else if(qN<bS)
        {
          G4int dS=aS-qN;
          n1=qN;
          n2=dS;
        }
        else
        {
          G4int dS=aS-qP;
          n1=dS;
          n2=qP;
        }
      }
    }
  }
  G4int qPDG=90000000+(qP-n2)*1000+(qN-n1);     // PDG of the Residual Nucleus
  G4double nucM = G4QNucleus(qPDG).GetGSMass(); // Mass of the Residual Nucleus
#ifdef pdebug
  G4cout<<"G4QEnv::DecayAnStran:nK0="<<n1<<",nK+="<<n2<<", nucM="<<nucM<<G4endl;
#endif
  G4int m1=0;                        // prototype of a#of K0's
  G4int m2=qP;                       // prototype of a#of K+'s
  if(qP>=-qS)   m2=-qS;              // Enough charge for K+'s
  else if(qP>0) m1=-qS-qP;           // Anti-Lambdas are partially compensated by neutrons
  G4int sPDG=90000000+(qP-m2)*1000+(qN-m1);     // PDG of the Residual Nucleus
  G4double mucM = G4QNucleus(sPDG).GetGSMass(); // Mass of the Residual Nucleus
  if(mucM+m1*mK+m2*mK0<nucM+n1*mK+n2*mK0)       // New is smaller
  {
    qPDG=sPDG;
    nucM=mucM;
    n1=m1;
    n2=m2;
  }
#ifdef pdebug
  G4cout<<"G4QEnv::DecayAnStran:n1="<<n1<<",n2="<<n2<<", nM="<<nucM<<G4endl;
#endif
  if(!n1||!n2)                            // AntiKaons of only one sort are found
  {
    if(!n1)                               // No K0's only K+'s
	{
      if(n2==1 && mK+nucM>qM+.0001)  // Mass limit: switch to K0
	  {
        k1M=mK0;
        n1=1;
        qPDG=90000000+qP*1000+qN-1;       // PDG of the Residual Nucleus
        nucM = G4QNucleus(qPDG).GetGSMass(); // Mass of the Residual Nucleus
      }
      else
      {
        k1M=mK;
        k1PDG=321;                        // Only K+'s (default K0's)
        n1=n2;                            // only n1 is used
	  }
	}
    else                                  // No K+'s only K0's
    {
      if(n1==1 && mK0+nucM>qM+.0001) // Mass limit: switch to K+
	  {
        k1M=mK;
        k1PDG=321;                        // K+ instead of K0
        qPDG=90000000+(qP-1)*1000+qN;     // PDG of the Residual Nucleus
        nucM = G4QNucleus(qPDG).GetGSMass(); // Mass of the Residual Nucleus
      }
      else k1M=mK0;                      // Only anti-K0's (default k1PDG)
    }
#ifdef pdebug
    G4int naPDG=90000000+(qP-1)*1000+qN; // Prot PDG of the Alternative Residual Nucleus
    G4double naM=G4QNucleus(naPDG).GetGSMass(); // Prot Mass of the Alt. Residual Nucleus
    G4double kaM=mK;                     // Prot Mass of the Alternative kaon (K+)
    if(k1PDG==321)                       // Calculate alternative to K+
    {
      naPDG=90000000+qP*1000+qN-1;       // PDG of the Alternative Residual Nucleus
      naM=G4QNucleus(naPDG).GetGSMass(); // Mass of the Alt. Residual Nucleus
      kaM=mK0;                           // Prot Mass of the Alternative kaon (K0)
    }
    G4cout<<"G4QEnv::DecayAnStran:M="<<qM<<",kM="<<k1M<<"+nM="<<nucM<<"="<<k1M+nucM
          <<",m="<<kaM<<"+n="<<naM<<"="<<kaM+naM<<G4endl;
#endif
    G4double n1M=n1*k1M;
    G4LorentzVector f4Mom(0.,0.,0.,n1M);
    G4LorentzVector s4Mom(0.,0.,0.,nucM);
    G4double sum=nucM+n1M;
    if(sum>qM+eps && n1==1)              // Try to use another K if this is the only kaon
    {
      G4int naPDG=90000000+(qP-1)*1000+qN; // Prot PDG of the Alternative Residual Nucleus
      G4double naM=G4QNucleus(naPDG).GetGSMass(); // Prot Mass of the Alt. Residual Nucleus
      G4int akPDG=321;                   // Prototype PDGCode of the AlternativeKaon (K+)
      G4double kaM=mK;                   // Prototype Mass of the AlternativeKaon (K+)
      if(k1PDG==321)                     // Calculate alternative to the K+ meson
      {
        naPDG=90000000+qP*1000+qN-1;     // PDG of the Alternative Residual Nucleus
        naM=G4QNucleus(naPDG).GetGSMass(); // Mass of the Alt. Residual Nucleus
        akPDG=311;                       // PDG Code of the Alternative kaon (K0)
        kaM=mK0;                         // Mass of the Alternative kaon (K0)
      }
      G4double asum=naM+kaM;
      if(asum<sum)                       // Make a KSwap correction
      {
        nucM=naM;
        n1M=kaM;
        k1M=kaM;
        k1PDG=akPDG;
        qPDG=naPDG;
        f4Mom=G4LorentzVector(0.,0.,0.,n1M);
        s4Mom=G4LorentzVector(0.,0.,0.,nucM);
      }
    }
    if(fabs(qM-sum)<eps)
	{
      f4Mom=q4M*(n1M/sum);
      s4Mom=q4M*(nucM/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
#ifdef pdebug
	  G4cerr<<"***G4QE::DecayAntiS:fPDG="<<n1<<"*"<<k1PDG<<"(M="<<k1M<<") + sPDG= "<<qPDG
            <<"(M="<<nucM<<" = "<<sum<<" > TotM="<<qM<<q4M<<G4endl;
#endif
      if(theEnvironment.GetA()==0)
	  {
        G4Quasmon* quasH = new G4Quasmon(qQC,q4M);
		if(!CheckGroundState(quasH,true))
        {
          G4cerr<<"---Worning---G4QE::DecAntiStran:Failed FillAsItIs h="<<qQC<<q4M<<G4endl;
          theQHadrons.push_back(qH);// @@ Can cause problems with ParticleConversion in G4
        }
        else delete qH;  
        delete quasH;
        return;
      }
	  else
      {
        G4cerr<<"---Worning---G4QE::DAS:AsItIsE="<<theEnvironment<<",h="<<qQC<<q4M<<G4endl;
        theQHadrons.push_back(qH);  // @@ Can cause problems with particle conversion in G4
        return;
        //throw G4QException("G4QE::DecayAntiStrange:AntStrangeNuc DecIn2 didn't succeed");
      }
    }
#ifdef pdebug
    G4cout<<"G4QEnv::DecAntiS: nK+N "<<n1<<"*K="<<k1PDG<<f4Mom<<",N="<<qPDG<<s4Mom<<G4endl;
#endif
    delete qH;
    //
    f4Mom/=n1;
    for(G4int i1=0; i1<n1; i1++)
    {
      G4QHadron* H1 = new G4QHadron(k1PDG,f4Mom); // Create a Hadron for the Kaon
      theQHadrons.push_back(H1);                  // Fill "H1" (delete equivalent)
	}
    G4QHadron* H2 = new G4QHadron(qPDG,s4Mom);   // Create a Hadron for the Nucleus
    //theQHadrons.push_back(H2);                 // Fill "H2" (delete equivalent)
    EvaporateResidual(H2);                       // Fill "H2" (delete equivalent)
#ifdef debug
    G4cout<<"G4QEnv::DecAntiS:*** After EvaporateResidual nH="<<theQHadrons.size()<<G4endl;
#endif
  }
  else
  {
    G4double n1M=n1*k1M;
    G4double n2M=n2*k2M;
    G4LorentzVector f4Mom(0.,0.,0.,n1M);
    G4LorentzVector s4Mom(0.,0.,0.,n2M);
    G4LorentzVector t4Mom(0.,0.,0.,nucM);
    G4double sum=nucM+n1M+n2M;
    if(fabs(qM-sum)<eps)
	{
      f4Mom=q4M*(n1M/sum);
      s4Mom=q4M*(n2M/sum);
      t4Mom=q4M*(nucM/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn3(f4Mom, s4Mom, t4Mom))
    {
      G4cerr<<"---Worning---G4QE::DeAS:nPDG="<<qPDG<<"(M="<<nucM<<")+1="<<k1PDG<<"(M="<<k1M
            <<")+2="<<k2PDG<<"(M="<<k2M<<")="<<nucM+n1*k1M+n2*k2M<<">tM="<<qM<<q4M<<G4endl;
      if(theEnvironment.GetA()==0)
	  {
        G4Quasmon* quasH = new G4Quasmon(qQC,q4M);
		if(!CheckGroundState(quasH,true))
        {
          G4cerr<<"---Worning---G4QE::DecAntiS:NoCorDone,FillAsItIs h="<<q4M<<qQC<<G4endl;
          theQHadrons.push_back(qH); //@@ Can cause problems with particle conversion in G4
        }
        else delete qH;  
        delete quasH;
        return;
      }
	  else
      {
        G4cerr<<"---Worning---G4QE::DAS: NoCorrDone,FillAsItIs E="<<theEnvironment<<G4endl;
        theQHadrons.push_back(qH); // @@ Can cause problems with particle conversion in G4
        return;
        //throw G4QException("G4QEnv::DecayAntiStrange:AntiStrangeNucleus DecIn3 error");
      }
    }
#ifdef pdebug
    G4cout<<"G4QEnv::DecAntiS:*DONE* nPDG="<<qPDG<<",1="<<f4Mom<<",2="<<s4Mom<<",n="<<t4Mom
          <<G4endl;
#endif
    delete qH;
    //
    f4Mom/=n1;
    for(G4int i1=0; i1<n1; i1++)
    {
      G4QHadron* H1 = new G4QHadron(k1PDG,f4Mom); // Create a Hadron for the K0
      theQHadrons.push_back(H1);               // Fill "H1" (delete equivalent)
	}
    s4Mom/=n2;
    for(G4int i2=0; i2<n2; i2++)
    {
      G4QHadron* H2 = new G4QHadron(k2PDG,s4Mom); // Create a Hadron for the K+
      theQHadrons.push_back(H2);                  // Fill "H2" (delete equivalent)
	}
    G4QHadron* H3 = new G4QHadron(qPDG,t4Mom);    // Create a Hadron for the nucleus
    //theQHadrons.push_back(H3);                  // Fill "H3" (delete equivalent)
    EvaporateResidual(H3);                        // Fill "H3" (delete equivalent)
  }
#ifdef debug
    G4cout<<"G4QEnv::DecAntiS: ===> End of DecayAntiStrangness"<<G4endl;
#endif
} // End of DecayAntiStrange

//Decay of the excited 3p or 3n systems in three baryons
void G4QEnvironment::DecayMultyBaryon(G4QHadron* qH)
{//  ===============================================
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double eps=0.003;
  G4LorentzVector q4M = qH->Get4Momentum();       // Get 4-momentum of the MultyBaryon
  G4double         qM = q4M.m();                  // Mass of the Multybaryon
  G4int          qPDG = qH->GetPDGCode();         // PDG Code of the decaying multybar
  G4QContent      qQC = qH->GetQC();              // PDG Code of the decaying multibar
#ifdef pdebug
  G4cout<<"G4QEnv::DecayMultyBaryon: *Called* PDG="<<qPDG<<",4M="<<q4M<<qQC<<G4endl;
#endif
  G4int totS=qQC.GetStrangeness();  //  Total Strangeness       (L)                ^
  G4int totC=qQC.GetCharge();       //  Total Charge            (p)                ^
  G4int totBN=qQC.GetBaryonNumber();// Total Baryon Number      (A)                ^
  G4int totN=totBN-totS-totC;       // Total Number of Neutrons (n)                ^
  G4int          fPDG = 3122;       // Prototype for A lambdas case
  G4double       fMass= mLamb;
  if     (totN==totBN)              // "A-neutron" case
  {
    fPDG = 2112;
    fMass= mNeut;
  }
  else if(totC==totBN)              // "A-protons" case
  {
    fPDG = 2212;
    fMass= mProt;
  }
  else if(totS!=totBN)            // "Bad call" case
  {
    G4cerr<<"***G4QEnv::DecayMultyBaryon: PDG="<<qPDG<<G4endl;
    throw G4QException("***G4QEnv::DecayMultyBaryon: Can not decay this PDG Code");
  }
#ifdef pdebug
  else
  {
    G4cerr<<"**G4QEnv::DecayMultyBaryon: PDG="<<qPDG<<G4endl;
    throw G4QException("***G4QEnv::DecayMultyBaryon: Unknown PDG code of the MultiBaryon");
  }
#endif
  if(totBN==1)theQHadrons.push_back(qH);
  else if(totBN==2)
  {
    G4LorentzVector f4Mom(0.,0.,0.,fMass);
    G4LorentzVector s4Mom(0.,0.,0.,fMass);
    G4double sum=fMass+fMass;
    if(fabs(qM-sum)<eps)
	{
      f4Mom=q4M/2.;
      s4Mom=f4Mom;
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cerr<<"---Worning---G4QEnv::DecayMultyBar:fPDG="<<fPDG<<"(fM="<<fMass<<")*2="<<sum
            <<" > TotM="<<q4M.m()<<q4M<<G4endl;
      if(!theEnvironment.GetA())
      {
        G4QContent qQC=qH->GetQC();
        G4Quasmon* quasH = new G4Quasmon(qQC,q4M);
	    if(!CheckGroundState(quasH,true))
        {
          G4cerr<<"---Worning---G4QE::DecMultyBar:FillAsItIsDiBar(I=1)="<<qQC<<q4M<<G4endl;
          theQHadrons.push_back(qH);               // Fill as it is
        }
        else delete qH;  
        delete quasH;
        return;
      }
      else
      {
        G4cerr<<"***G4QEnv::DecayMultyBaryon:Cann'tCor*EmptyEnv*="<<theEnvironment<<G4endl;
        throw G4QException("G4QEnv::DecayMultyBaryon:diBaryon DecayIn2 didn't succeed");
      }
    }
#ifdef pdebug
    G4cout<<"G4QEnv::DecMBar:*DONE* fPDG="<<fPDG<<",f="<<f4Mom<<",s="<<s4Mom<<G4endl;
#endif
    delete qH;
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom);   // Create a Hadron for the 1-st baryon
    theQHadrons.push_back(H1);                   // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(fPDG,s4Mom);   // Create a Hadron for the 2-nd baryon
    theQHadrons.push_back(H2);                   // Fill "H2" (delete equivalent)
  }
  else if(totBN==3)
  {
    G4LorentzVector f4Mom(0.,0.,0.,fMass);
    G4LorentzVector s4Mom(0.,0.,0.,fMass);
    G4LorentzVector t4Mom(0.,0.,0.,fMass);
    G4double sum=fMass+fMass+fMass;
    if(fabs(qM-sum)<eps)
	{
      f4Mom=q4M/3.;
      s4Mom=f4Mom;
      t4Mom=f4Mom;
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn3(f4Mom, s4Mom, t4Mom))
    {
      G4cerr<<"---Worning---G4QEnv::DecayMultyBaryon: fPDG="<<fPDG<<"(fM="<<fMass<<")*3 = "
            <<3*fMass<<" >? TotM="<<q4M.m()<<q4M<<G4endl;
      if(!theEnvironment.GetA())
      {
        G4QContent qQC=qH->GetQC();
        G4Quasmon* quasH = new G4Quasmon(qQC,q4M);
	    if(!CheckGroundState(quasH,true))
        {
          G4cerr<<"---Worning---G4QEnv::DecMultyBar:FillAsItIs 3Baryon="<<qQC<<q4M<<G4endl;
          theQHadrons.push_back(qH);             // Fill as it is
        }
        else delete qH;  
        delete quasH;
        return;
      }
      else
      {
        G4cerr<<"***G4QE::DecayMultyBar:Cann't Correct*EmptyEnv*="<<theEnvironment<<G4endl;
        throw G4QException("G4QEnv::DecayMultyBar:ThreeBaryon DecayIn3 didn't succeed");
      }
    }
#ifdef pdebug
    G4cout<<"G4QEnv::DecMBar:*DONE*, fPDG="<<fPDG<<",f="<<f4Mom<<",s="<<s4Mom<<",t="
          <<t4Mom<<G4endl;
#endif
    delete qH;
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom);   // Create a Hadron for the 1-st baryon
    theQHadrons.push_back(H1);                   // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(fPDG,s4Mom);   // Create a Hadron for the 2-nd baryon
    theQHadrons.push_back(H2);                   // Fill "H2" (delete equivalent)
    G4QHadron* H3 = new G4QHadron(fPDG,t4Mom);   // Create a Hadron for the 3-d baryon
    theQHadrons.push_back(H3);                   // Fill "H3" (delete equivalent)
  }
  else
  {
    // @@It must be checked, that they are not under the mass shell
    // !! OK !! Checked by the worning print that they are mostly in the Ground State !!
    G4LorentzVector f4Mom=q4M/totBN; // @@ Too simple solution (split in two parts!)
#ifdef pdebug
    // Worning for the future development
    G4cout<<"**G4QE::DecMulBar:SplitMultiBar inEqParts M="<<totBN<<"*"<<f4Mom.m()<<G4endl;
    G4cout<<"G4QEnv::DecMBar:*DONE* fPDG="<<fPDG<<",f="<<f4Mom<<G4endl;
#endif
    delete qH;
    for(G4int h=0; h<totBN; h++)
    {
      G4QHadron* H1 = new G4QHadron(fPDG,f4Mom); // Create a Hadron for the baryon
      theQHadrons.push_back(H1);                 // Fill "H1" (delete equivalent)
	}
  }
} // End of DecayMultyBaryon

//Decay of the excited alpha+2p or alpha+2n systems
void G4QEnvironment::DecayAlphaDiN(G4QHadron* qH)
{//  ============================================
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double mHel6= G4QPDGCode(2112).GetNuclMass(2,4,0);
  static const G4double eps=0.003;
  G4LorentzVector q4M = qH->Get4Momentum();      // Get 4-momentum of the AlphaDibaryon
  G4double         qM = q4M.m();                 // Real mass of the AlphaDibaryon
  G4int          qPDG = qH->GetPDGCode();        // PDG Code of the decayin AlphaDybaryon
#ifdef pdebug
  G4cout<<"G4QEnv::DecayAlphaDiN: *Called* PDG="<<qPDG<<",4M="<<q4M<<G4endl;
#endif
  G4int          fPDG = 2212;                    // Prototype for alpha+pp case
  G4double       fMass= mProt;
  G4int          sPDG = 90002002;
  G4double       sMass= mAlph;
  if     (qPDG==90002004)                        // "alpha+2neutrons" case
  {
    if(abs(qM-mHel6)<eps)
	{
      theQHadrons.push_back(qH);                 // Fill as it is (delete equivalent)
      return;
	}
    else if(mNeut+mNeut+mAlph<qM)
    {
      fPDG = 2112;
      fMass= mNeut;
	}
    else
    {
      G4cerr<<"***G4QEn::DecAlDiN:M(He6="<<mHel6<<")="<<qM<<"<"<<mNeut+mNeut+mAlph<<G4endl;
      throw G4QException("G4QEnv::DecayAlphaDiN: Cannot decay excited He6 with this mass");
    }
  }
  else if(qPDG!=90004002)                         // "Bad call" case
  {
    G4cerr<<"***G4QEnv::DecayAlphaDiN: PDG="<<qPDG<<G4endl;
    throw G4QException("G4QEnv::DecayAlphaDiN: Can not decay this PDG Code");
  }
  G4LorentzVector f4Mom(0.,0.,0.,fMass);
  G4LorentzVector s4Mom(0.,0.,0.,fMass);
  G4LorentzVector t4Mom(0.,0.,0.,sMass);
  G4double sum=fMass+fMass+sMass;
  if(fabs(qM-sum)<eps)
  {
    f4Mom=q4M*(fMass/sum);
    s4Mom=f4Mom;
    t4Mom=q4M*(sMass/sum);
  }
  else if(qM<sum || !G4QHadron(q4M).DecayIn3(f4Mom, s4Mom, t4Mom))
  {
    G4int eA=theEnvironment.GetA();
    G4cerr<<"---Worning---G4QEnv::DecayAlphaDiN:fPDG="<<fPDG<<"(M="<<fMass<<")*2+mAlpha = "
          <<sum<<" >? TotM="<<qM<<q4M<<", d="<<sum-qM<<", envA="<<eA<<G4endl;
    if(!eA)
    {
      G4QContent qQC=qH->GetQC();
      G4Quasmon* quasH = new G4Quasmon(qQC,q4M);
	  if(!CheckGroundState(quasH,true))
      {
        G4cerr<<"---Worning---G4QEnv::DecayAlphaDiN:FillAsItIs AlphaNN="<<qQC<<q4M<<G4endl;
        theQHadrons.push_back(qH);      // Fill as it is
      }
      else delete qH;  
      delete quasH;
      return;
    }
    else
    {
      G4cerr<<"***G4QEnv::DecayAlphaDiN:Cann't Correct*EmptyEnv*="<<theEnvironment<<G4endl;
      throw G4QException("G4QEnv::DecayAlphaDiN: Alpha+N+N DecayIn3 error");
    }
  }
#ifdef pdebug
  G4cout<<"G4QE::DecAl2N: fPDG="<<fPDG<<",f="<<f4Mom<<",s="<<s4Mom<<",t="<<t4Mom<<G4endl;
#endif
  delete qH;
  G4QHadron* H1 = new G4QHadron(fPDG,f4Mom);    // Create a Hadron for the 1-st baryon
  theQHadrons.push_back(H1);                    // Fill "H1" (delete equivalent)
  G4QHadron* H2 = new G4QHadron(fPDG,s4Mom);    // Create a Hadron for the 2-nd baryon
  theQHadrons.push_back(H2);                    // Fill "H2" (delete equivalent)
  G4QHadron* H3 = new G4QHadron(sPDG,t4Mom);    // Create a Hadron for the alpha
  theQHadrons.push_back(H3);                    // Fill "H3" (delete equivalent)
} // End of DecayAlphaDiN

//Decay of the excited alpha+bayon state in alpha and baryons
void G4QEnvironment::DecayAlphaBar(G4QHadron* qH)
{//  ============================================
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double mTrit= G4QPDGCode(2112).GetNuclMass(1,2,0);
  static const G4double mHe3 = G4QPDGCode(2112).GetNuclMass(2,1,0);
  static const G4double eps=0.003;
  G4LorentzVector q4M = qH->Get4Momentum();     // Get 4-momentum of the Alpha-Baryon
  G4double         qM = q4M.m();                // Mass of Alpha-Baryon
  G4int          qPDG = qH->GetPDGCode();       // PDG Code of the decayin Alpha-Baryon
  G4QContent      qQC = qH->GetQC();            // PDG Code of the decaying Alpha-Bar
#ifdef pdebug
  G4cout<<"G4QEnv::DecayAlphaBar: *Called* PDG="<<qPDG<<",4M="<<q4M<<qQC<<G4endl;
#endif
  G4int totS=qQC.GetStrangeness();              //  Total Strangeness       (L)
  G4int totC=qQC.GetCharge();                   //  Total Charge            (p)
  G4int totBN=qQC.GetBaryonNumber();            // Total Baryon Number      (A)
  G4int totN=totBN-totS-totC;                   // Total Number of Neutrons (n)
  if((totN==totBN||totC==totBN||totS==totBN)&&totBN>1) DecayMultyBaryon(qH);
  else if(qPDG==92001002||qPDG==92002001) theQHadrons.push_back(qH);
  else if(qPDG==91003001||qPDG==91001003||qPDG==93001001)theQHadrons.push_back(qH);
  else if(qPDG==92000003||qPDG==92003000||qPDG==93000002||qPDG==93002000)
  {
    G4int          fPDG = 3122;                 // 1st Prototype for 2L+3n case
    G4double       fMass= mLamb;
    G4int          sPDG = 2112;
    G4double       sMass= mNeut;
    if     (qPDG==92003000)                     // "2L+3p" case
    {
      sPDG = 2212;
      sMass= mProt;
    }
    else if(qPDG==93000002)                     // "2n+3L" case
    {
      fPDG = 2112;
      fMass= mNeut;
      sPDG = 3122;
      sMass= mLamb;
    }
    else if(qPDG==93002000)                     // "2p+3L" case
    {
      fPDG = 2212;
      fMass= mProt;
      sPDG = 3122;
      sMass= mLamb;
    }
    G4double tfM=fMass+fMass;
    G4double tsM=sMass+sMass+sMass;
    G4LorentzVector f4Mom(0.,0.,0.,tfM);
    G4LorentzVector s4Mom(0.,0.,0.,tsM);
    G4double sum=tfM+tsM;
    if(fabs(qM-sum)<eps)
	{
      f4Mom=q4M*(tfM/sum);
      s4Mom=q4M*(tsM/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cerr<<"---Worning---G4QE::DecAlB:fPDG="<<fPDG<<"(M="<<fMass<<")*2="<<2*fMass<<",s="
            <<sPDG<<"(sM="<<sMass<<")*3="<<3*sMass<<"="<<sum<<">M="<<q4M.m()<<q4M<<G4endl;
      if(!theEnvironment.GetA())
      {
        G4QContent qQC=qH->GetQC();
        G4Quasmon* quasH = new G4Quasmon(qQC,q4M);
	    if(!CheckGroundState(quasH,true))
        {
          G4cerr<<"---Worning---G4QE::DecAlphBar:FillAsItIsAlphaBaryon="<<qQC<<q4M<<G4endl;
          theQHadrons.push_back(qH);            // Fill as it is
        }
        else delete qH;  
        delete quasH;
        return;
      }
      else
      {
        G4cerr<<"***G4QE::DecayAlphaBar:Cann't Correct*EmptyEnv*="<<theEnvironment<<G4endl;
        throw G4QException("G4QEnv::DecayAlphaBar: DecayIn2 didn't succeed for 3/2");
      }
    }
#ifdef pdebug
    G4cout<<"G4QEnv::DecAB:*DONE*, fPDG="<<fPDG<<f4Mom<<",sPDG="<<sPDG<<s4Mom<<G4endl;
#endif
    delete qH;
    G4LorentzVector rf4Mom=f4Mom/2;
    G4QHadron* H1 = new G4QHadron(fPDG,rf4Mom); // Create a Hadron for the 1-st baryon
    theQHadrons.push_back(H1);                  // Fill "H1" (delete equivalent)
    theQHadrons.push_back(H1);                  // Fill "H1" (delete equivalent)
    G4LorentzVector rs4Mom=s4Mom/3;
    G4QHadron* H2 = new G4QHadron(sPDG,rs4Mom); // Create a Hadron for the 2-nd baryon
    theQHadrons.push_back(H2);                  // Fill "H2" (delete equivalent)
    theQHadrons.push_back(H2);                  // Fill "H2" (delete equivalent)
    theQHadrons.push_back(H2);                  // Fill "H2" (delete equivalent)
  }
  else if(qPDG==90004001||qPDG==90001004)
  {
    G4int          fPDG = 90002001;             // Prototype for "He3+2p" case
    G4double       fMass= mHe3;
    G4int          sPDG = 2212;
    G4double       sMass= mProt;
    if     (qPDG==90001004)                     // "t+2n" case
    {
      fPDG = 90001002;
      fMass= mTrit;
      sPDG = 2112;
      sMass= mNeut;
    }
    G4LorentzVector f4Mom(0.,0.,0.,fMass);
    G4LorentzVector s4Mom(0.,0.,0.,sMass);
    G4LorentzVector t4Mom(0.,0.,0.,sMass);
    G4double sum=fMass+sMass+sMass;
    if(fabs(qM-sum)<eps)
	{
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
      t4Mom=s4Mom;
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn3(f4Mom, s4Mom, t4Mom))
    {
      G4cerr<<"---Worning---G4QE::DecAlB:fPDG="<<fPDG<<",M="<<fMass<<",sPDG="<<sPDG<<",sM="
            <<sMass<<",2sM+fM="<<2*sMass+fMass<<" > TotM="<<q4M.m()<<q4M<<G4endl;
      if(!theEnvironment.GetA())
      {
        G4QContent qQC=qH->GetQC();
        G4Quasmon* quasH = new G4Quasmon(qQC,q4M);
	    if(!CheckGroundState(quasH,true))
        {
          G4cerr<<"---Worning---G4QE::DecAlphBar:FillAsItIsAlphaBaryon="<<qQC<<q4M<<G4endl;
          theQHadrons.push_back(qH);      // Fill as it is
        }
        else delete qH;  
        delete quasH;
        return;
      }
      else
      {
        G4cerr<<"***G4QE::DecayAlphaBar:Cann't Correct*EmptyEnv*="<<theEnvironment<<G4endl;
        throw G4QException("G4QEnv::DecayAlphaBar: t/nn,He3/pp DecayIn3 didn't");
      }
    }
#ifdef pdebug
    G4cout<<"G4QE::DecAlB:*DONE*, f="<<fPDG<<f4Mom<<", s="<<sPDG<<s4Mom<<t4Mom<<G4endl;
#endif
    delete qH;
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom);   // Create a Hadron for the 1-st baryon
    theQHadrons.push_back(H1);                   // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(sPDG,s4Mom);   // Create a Hadron for the 2-nd baryon
    theQHadrons.push_back(H2);                   // Fill "H2" (delete equivalent)
    G4QHadron* H3 = new G4QHadron(sPDG,t4Mom);   // Create a Hadron for the 3-d baryon
    theQHadrons.push_back(H3);                   // Fill "H3" (delete equivalent)
  }
  else if(qPDG==94000001||qPDG==94001000||qPDG==91000004||qPDG==91004000)
  {
    G4int          fPDG = 3122;                 // Prototype for "4L+n" case
    G4double       fMass= mLamb+mLamb;
    G4int          sPDG = 2112;
    G4double       sMass= mNeut;
    if     (qPDG==94001000)                     // "4L+p" case
    {
      sPDG = 2212;
      sMass= mProt;
    }
    else if(qPDG==91000004)                     // "4n+L" case
    {
      fPDG = 2112;
      fMass= mNeut+mNeut;
      sPDG = 3122;
      sMass= mLamb;
    }
    else if(qPDG==91004000)                     // "4p+L" case
    {
      fPDG = 2212;
      fMass= mProt+mProt;
      sPDG = 3122;
      sMass= mLamb;
    }
    G4LorentzVector f4Mom(0.,0.,0.,fMass+fMass);
    G4LorentzVector s4Mom(0.,0.,0.,sMass);
    G4double sum=fMass+fMass+sMass;
    if(fabs(qM-sum)<eps)
	{
      f4Mom=q4M*((fMass+fMass)/sum);
      s4Mom=q4M*(sMass/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cerr<<"---Worning---G4QE::DecAlphBar:fPDG="<<fPDG<<"(2*fM="<<fMass<<")*2="<<2*fMass
            <<",sPDG="<<sPDG<<"(sM="<<sMass<<" > TotM="<<q4M.m()<<q4M<<G4endl;
      if(!theEnvironment.GetA())
      {
        G4QContent qQC=qH->GetQC();
        G4Quasmon* quasH = new G4Quasmon(qQC,q4M);
	    if(!CheckGroundState(quasH,true))
        {
          G4cerr<<"---Worning---G4QE::DecAlphBar:FillAsItIsAlphaBaryon="<<qQC<<q4M<<G4endl;
          theQHadrons.push_back(qH);            // Fill as it is
        }
        else delete qH;  
        delete quasH;
        return;
      }
      else
      {
        G4cerr<<"***G4QE::DecayAlphaBar:Cann't Correct*EmptyEnv*="<<theEnvironment<<G4endl;
        throw G4QException("G4QEnv::DecayAlphaBar:QuintBaryon DecayIn2 didn't succeed");
      }
    }
#ifdef pdebug
    G4cout<<"G4QEnv::Dec5B:*DONE*, fPDG="<<fPDG<<f4Mom<<",sPDG="<<sPDG<<s4Mom<<G4endl;
#endif
    delete qH;
    G4LorentzVector rf4Mom=f4Mom/4;
    G4QHadron* H1 = new G4QHadron(fPDG,rf4Mom); // Create a Hadron for the 1-st baryon
    theQHadrons.push_back(H1);                  // Fill "H1" (delete equivalent)
    theQHadrons.push_back(H1);                  // Fill "H1" (delete equivalent)
    theQHadrons.push_back(H1);                  // Fill "H1" (delete equivalent)
    theQHadrons.push_back(H1);                  // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(sPDG,s4Mom);  // Create a Hadron for the 2-nd baryon
    theQHadrons.push_back(H2);                  // Fill "H2" (delete equivalent)
  }
  else if(qPDG==90003002||qPDG==90002003||qPDG==91002002)
  {
    G4int          fPDG = 90002002;             // Prototype for "alpha+n" case
    G4int          sPDG = 2112;
    G4double       fMass= mAlph;
    G4double       sMass= mNeut;
    if(qPDG==90003002)                          // "alpha+p" case
    {
      sPDG = 2212;
      sMass= mProt;    
    }
    else if(qPDG==9100202)                      // "alpha+l" case
    {
      sPDG = 3122;
      sMass= mLamb;    
    }
    else if(qPDG!=90002003)
    {
      //theQHadrons.push_back(qH);              // Fill hadron as it is (delete equivalent)
      EvaporateResidual(qH);                    // Evaporate ResNuc (delete equivivalent)
      return;
    }
    G4double dM=fMass+sMass-qM;
    if(dM>0.&&dM<1.)
    {
#ifdef pdebug
      G4cerr<<"***Corrected***G4QEnv::DecayAlphaBar:fPDG="<<fPDG<<"(fM="<<fMass<<")+ sPDG="
            <<sPDG<<"(sM="<<sMass<<")="<<fMass+sMass<<" > TotM="<<qM<<q4M<<G4endl;
#endif
      G4double hdM=dM/2;
      fMass-=hdM;
      sMass-=hdM;
    }
    G4LorentzVector f4Mom(0.,0.,0.,fMass);
    G4LorentzVector s4Mom(0.,0.,0.,sMass);      // Mass is random since probab. time
    G4double sum=fMass+sMass;
    if(fabs(qM-sum)<eps)
	{
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cerr<<"---Worning---G4QE::DecAlphaBar:fPDG="<<fPDG<<"(fM="<<fMass<<")+sPDG="<<sPDG
            <<"(sM="<<sMass<<")="<<fMass+sMass<<"="<<sum<<" > TotM="<<q4M.m()<<q4M<<G4endl;
      if(!theEnvironment.GetA())
      {
        G4QContent qQC=qH->GetQC();
        G4Quasmon* quasH = new G4Quasmon(qQC,q4M);
	    if(!CheckGroundState(quasH,true))
        {
          G4cerr<<"---Worning---G4QE::DecAlphBar:FillAsItIsAlphaBaryon="<<qQC<<q4M<<G4endl;
          theQHadrons.push_back(qH);      // Fill as it is
        }
        else delete qH;  
        delete quasH;
        return;
      }
      else
      {
        G4cerr<<"***G4QE::DecayAlphaBar:Cann't Correct*EmptyEnv*="<<theEnvironment<<G4endl;
        throw G4QException("***G4QE::DecayAlphaBar:Alpha+Baryon DecIn2 didn't succeed");
      }
    }
#ifdef pdebug
    G4cout<<"G4QEnv::DecAlBar:*DONE*a4M="<<f4Mom<<",s4M="<<s4Mom<<",sPDG="<<sPDG<<G4endl;
#endif
    delete qH;
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom);      // Create a Hadron for the alpha
    theQHadrons.push_back(H1);                      // Fill "H1" (delete equivalent)
    G4QHadron* H2 = new G4QHadron(sPDG,s4Mom);      // Create a Hadron for the baryon
    theQHadrons.push_back(H2);                      // Fill "H2" (delete equivalent)
  }
  else G4cerr<<"---Worning---G4QEnv::DecayAlphaBar: Unknown PDG="<<qPDG<<G4endl;
} // End of DecayAlphaBar

//Decay of the excited alpha+alpha state in 2 alphas
void G4QEnvironment::DecayAlphaAlpha(G4QHadron* qH)
{//  ==============================================
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);
  static const G4double aaGSM= G4QPDGCode(2112).GetNuclMass(4,4,0);
  static const G4double eps=0.003;
  G4int          qPDG = qH->GetPDGCode();         // PDG Code of the decayin dialpha
  if(qPDG!=90004004)
  {
    G4cerr<<"***G4QEnv::DecayAlphaAlpha: qPDG="<<qPDG<<G4endl;
    throw G4QException("***G4QEnv::DecayAlphaAlpha: Not Be8 state decais in 2 alphas");
  }
  G4LorentzVector q4M = qH->Get4Momentum();       // Get 4-momentum of the Dibaryon
  G4double qM=q4M.m();
#ifdef pdebug
  G4cout<<"G4QEnv::DecayAlAl: *Called* PDG="<<qPDG<<",M="<<qM<<q4M<<">"<<aaGSM<<G4endl;
#endif
  //if(qM>aaGSM+.01)  // @@ Be8*->gamma+Be8 (as in evaporation)
  if(2>3)
  {
    G4int          fPDG = 22;
    G4int          sPDG = 90004004;
    G4double       fMass= 0.;
    G4double       sMass= aaGSM;
    G4LorentzVector f4Mom(0.,0.,0.,fMass);
    G4LorentzVector s4Mom(0.,0.,0.,sMass);          // Mass is random since probab. time
    G4double sum=fMass+sMass;
    if(fabs(qM-sum)<eps)
	{
      f4Mom=q4M*(fMass/sum);
      s4Mom=q4M*(sMass/sum);
    }
    else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
    {
      G4cerr<<"---Worning---G4QEnv::DecayAlphaAlpha:gPDG="<<fPDG<<"(gM="<<fMass<<")+PDG="
            <<sPDG<<"(sM="<<sMass<<")="<<sum<<" > TotM="<<q4M.m()<<q4M<<G4endl;
      if(!theEnvironment.GetA())
      {
        G4QContent qQC=qH->GetQC();
        G4Quasmon* quasH = new G4Quasmon(qQC,q4M);
	    if(!CheckGroundState(quasH,true))
        {
          G4cerr<<"---Worning---G4QE::DecAlphAlph:FillAsItIsAlphaAlpha="<<qQC<<q4M<<G4endl;
          theQHadrons.push_back(qH);      // Fill as it is
        }
        else delete qH;  
        delete quasH;
        return;
      }
      else
      {
        G4cerr<<"***G4QEnv::DecayAlphAlph:Cann't Corr.*EmptyEnv*="<<theEnvironment<<G4endl;
        throw G4QException("G4QEnv::DecayAlphaAlpha: g+diAlph DecayIn2 didn't succeed");
      }
    }
#ifdef pdebug
    G4cout<<"G4QEnv::DecayAlphaAlpha: *DONE* gam4M="<<f4Mom<<", aa4M="<<s4Mom<<G4endl;
#endif
    G4QHadron* H1 = new G4QHadron(fPDG,f4Mom);      // Create a Hadron for the 1-st alpha
    theQHadrons.push_back(H1);                      // Fill "H1" (delete equivalent)
    qH->Set4Momentum(s4Mom);
    q4M=s4Mom;
  }
  G4int          fPDG = 90002002;
  G4int          sPDG = 90002002;
  G4double       fMass= mAlph;
  G4LorentzVector f4Mom(0.,0.,0.,fMass);
  G4LorentzVector s4Mom(0.,0.,0.,fMass);
  G4double sum=fMass+fMass;
  if(fabs(qM-sum)<eps)
  {
    f4Mom=q4M*(fMass/sum);
    s4Mom=f4Mom;
  }
  else if(qM<sum || !G4QHadron(q4M).DecayIn2(f4Mom, s4Mom))
  {
    G4cerr<<"---Worning---G4QEnv::DecayAlphaAlpha:fPDG="<<fPDG<<"(fM="<<fMass<<")*2="<<sum
          <<" > TotM="<<q4M.m()<<q4M<<G4endl;
    if(!theEnvironment.GetA())
    {
      G4QContent qQC=qH->GetQC();
      G4Quasmon* quasH = new G4Quasmon(qQC,q4M);
	  if(!CheckGroundState(quasH,true))
      {
        G4cerr<<"---Worning---G4QE::DecAlphaAlpha:FillAsItIsAlphaAlpha="<<qQC<<q4M<<G4endl;
        theQHadrons.push_back(qH);      // Fill as it is
      }
      else delete qH;  
      delete quasH;
      return;
    }
    else
    {
      G4cerr<<"***G4QEnv::DecayAlphAlph:Cann't Correct*EmptyEnv*="<<theEnvironment<<G4endl;
      throw G4QException("G4QEnv::DecayAlphaAlpha: diAlpha DecayIn2 didn't succeed");
    }
  }
#ifdef pdebug
  G4cout<<"G4QEnv::DecayAlphaAlpha: *DONE* fal4M="<<f4Mom<<", sal4M="<<s4Mom<<G4endl;
#endif
  delete qH;
  G4QHadron* H1 = new G4QHadron(fPDG,f4Mom);      // Create a Hadron for the 1-st alpha
  theQHadrons.push_back(H1);                      // Fill "H1" (delete equivalent)
  G4QHadron* H2 = new G4QHadron(sPDG,s4Mom);      // Create a Hadron for the 2-nd alpha
  theQHadrons.push_back(H2);                      // Fill "H2" (delete equivalent)
} // End of DecayAlphaAlpha

// Check that it's possible to decay the TotalResidualNucleus in Quasmon+Environ & correct
G4bool G4QEnvironment::CheckGroundState(G4Quasmon* quasm, G4bool corFlag)
{ //   ==================================================================
  static const G4QContent neutQC(2,1,0,0,0,0);
  static const G4QContent protQC(1,2,0,0,0,0);
  static const G4QContent lambQC(1,1,1,0,0,0);
  static const G4QContent pimQC(1,0,0,0,1,0);
  static const G4QContent pipQC(0,1,0,1,0,0);
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mLamb= G4QPDGCode(3122).GetMass();
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  ///@@@///
  ///////////////corFlag=true;
  ///@@@///
  G4QContent valQ=quasm->GetQC();               // Quark content of the Quasmon
  G4int    resQPDG=valQ.GetSPDGCode();          // Reachable in a member function
  if(resQPDG==10&&valQ.GetBaryonNumber()>0) resQPDG=valQ.GetZNSPDGCode();
  G4double resQMa=G4QPDGCode(resQPDG).GetMass();// GS Mass of the Residual Quasmon
  G4double resEMa=0.;                           // GS Mass of the EmptyResidualEnvironment
  G4bool   bsCond=false;                        // FragSeparatCondition for QuasmonInVacuum
  G4LorentzVector enva4M=G4LorentzVector(0.,0.,0.,0.);
  G4QContent reTQC=valQ;                        // Prototype QuarkContent of the ResidNucl
  G4LorentzVector reTLV=quasm->Get4Momentum();  // Prototyoe 4-Mom of the Residual Nucleus
#ifdef pdebug
  G4cout<<"G4QE::CGS: "<<valQ<<",PDG="<<resQPDG<<",GSM="<<resQMa<<",M="<<reTLV.m()<<G4endl;
#endif
  G4double resSMa=resQMa;                       // Prototype MinimalSplitMass of ResidNucl
  G4int envPDG=theEnvironment.GetPDG();
  if(envPDG!=90000000)                          // => "Environment is not vacuum" case
  { // @@@@@@@@@@@@@@@@@@@ CALL SUBROUTINE @@@@@@@@@
    resEMa=theEnvironment.GetMZNS();            // GSMass of the Residual Environment
    enva4M=theEnvironment.Get4Momentum();       // 4-Mom of the Residual Environment
#ifdef pdebug
	G4cout<<"G4QE::CGS: Environ exists gsM="<<resEMa<<",4M="<<enva4M<<enva4M.m()<<G4endl;
#endif
    reTQC+=theEnvironment.GetQCZNS();           // Quark content of the Residual Nucleus
    reTLV+=enva4M;                              // 4-Mom of Residual Nucleus
    //resSMa+=resEMa;                           // Minimal Split Mass of Residual Nucleus
    resSMa=G4QPDGCode(reTQC).GetMass();         // GS Mass of the Residual Quasmon+Environ
  }
  else                                          // Calculate BaryonSeparatCond for vacQuasm
  {
    G4double resQM=reTLV.m();                   // CM Mass of the Residual vacQuasmon
    G4int   baryn=valQ.GetBaryonNumber();       // Baryon Number of the Residual vacQuasmon
    if(baryn>1)                                 // => "Can split baryon" case
	{
      if(valQ.GetN())                           // ===> "Can split neutron" case
	  {
        G4QContent resQC=valQ-neutQC;           // QC of Residual for the Neutron
        G4int    resPDG=resQC.GetSPDGCode();    // PDG of Residual for the Neutron
        if(resPDG==10&&resQC.GetBaryonNumber()>0) resPDG=resQC.GetZNSPDGCode();
        G4double resMas=G4QPDGCode(resPDG).GetMass(); // GS Mass of the Residual
        if(resMas+mNeut<resQM) bsCond=true;
	  }
      else if(valQ.GetP())                      // ===> "Can split proton" case
	  {
        G4QContent resQC=valQ-protQC;           // QC of Residual for the Proton
        G4int    resPDG=resQC.GetSPDGCode();    // PDG of Residual for the Proton
        if(resPDG==10&&resQC.GetBaryonNumber()>0) resPDG=resQC.GetZNSPDGCode();
        G4double resMas=G4QPDGCode(resPDG).GetMass(); // GS Mass of the Residual
        if(resMas+mProt<resQM) bsCond=true;
	  }
      else if(valQ.GetL())                      // ===> "Can split lambda" case
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
#ifdef pdebug
  G4int    reTPDG=reTQC.GetSPDGCode();
  if(reTPDG==10&&reTQC.GetBaryonNumber()>0) reTPDG=reTQC.GetZNSPDGCode();
  G4cout<<"G4QEnv::CheckGS:(tM="<<resTMa<<"<rQM+rEM="<<resSMa<<",d="<<resSMa-resTMa
        <<" || rEM="<<resEMa<<"=0 & "<<!bsCond<<"=1) & n="<<nOfOUT<<">0 & F="<<corFlag
        <<" then the correction must be done for PDG="<<reTPDG<<G4endl;
#endif
  if((resTMa<resSMa || !resEMa&&!bsCond) && nOfOUT>0 && corFlag) // *** CORRECTION ***
  {
    G4QHadron*  theLast = theQHadrons[nOfOUT-1];
    G4int cNf=theLast->GetNFragments();
    G4int crPDG=theLast->GetPDGCode();
#ifdef pdebug
	G4cout<<"G4QE::CGS:***Correction***, lastNF="<<cNf<<",lastPDG="<<crPDG<<G4endl;
#endif
    G4LorentzVector hadr4M = theLast->Get4Momentum();
    G4double  hadrMa=hadr4M.m();
    G4LorentzVector tmpTLV=reTLV+hadr4M;        // Tot (ResidNucl+LastHadron) 4-Mom
#ifdef pdebug
	G4cout<<"G4QE::CGS:YES,4M/M="<<tmpTLV<<tmpTLV.m()<<" > rSM+hM="<<resSMa+hadrMa<<G4endl;
#endif
    G4double tmpTM=tmpTLV.m();
    if(tmpTM>resSMa+hadrMa &&!cNf && crPDG!=22) // Q(E)L contain QM(+EM)+lM ***Last CORR***
    {
      if(resEMa)                                // => "NonVacuumEnvironment exists" case
      {
        G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,resQMa); // GS Mass of Quasmon
        G4QHadron* quasH = new G4QHadron(valQ, quas4M);
        G4QHadron* envaH = new G4QHadron(theEnvironment.GetQCZNS(),enva4M);
        if(tmpTM>=resQMa+resEMa+hadrMa && G4QHadron(tmpTLV).DecayIn3(hadr4M,quas4M,enva4M))
        {
          //@@CHECK CoulBar (only for ResQuasmon in respect to ResEnv) and may be evaporate
          theLast->Set4Momentum(hadr4M);
          quasH->Set4Momentum(quas4M);
          if(resQPDG==92000000||resQPDG==90002000||resQPDG==90000002)DecayDibaryon(quasH);
          else if(resQPDG==93000000||resQPDG==90003000||resQPDG==90000003)
                                                                   DecayMultyBaryon(quasH);
          else if(resQPDG==90004002) DecayAlphaDiN(quasH);//Decay alph+2p(alph+2n isStable)
          else if(resQPDG==90002003||resQPDG==90003002) DecayAlphaBar(quasH); //DelEqu
          else if(resQPDG==90004004) DecayAlphaAlpha(quasH); //DelEqu
          else theQHadrons.push_back(quasH);    // Fill ResidQ as QHadron (delete equiv.)
          envaH->Set4Momentum(enva4M);
          if(envPDG==92000000||envPDG==90002000||envPDG==90000002) DecayDibaryon(envaH);
          else if(envPDG==93000000||envPDG==90003000||envPDG==90000003)
                                                                   DecayMultyBaryon(envaH);
          else if(envPDG==90002003||envPDG==90003002) DecayAlphaBar(envaH); //DelEqu
          else if(envPDG==90004002) DecayAlphaDiN(envaH);//Decay alph+2p(alph+2n IsStable)
          else if(envPDG==90004004) DecayAlphaAlpha(envaH); //DelEqu
          else theQHadrons.push_back(envaH);    // Fill 2nd Hadron (delete equivalent)
		}
        else
        {
          delete quasH;                         // Delete "Quasmon Hadron"
          delete envaH;                         // Delete "Environ Hadron"
#ifdef pdebug
          G4cout<<"***G4QEnv::CheckGS: Decay in Frag+ResQ+ResE did not succeeded"<<G4endl;
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
          quasH->Set4Momentum(quas4M);
          if(resQPDG==92000000||resQPDG==90002000||resQPDG==90000002)DecayDibaryon(quasH);
          else if(resQPDG==93000000||resQPDG==90003000||resQPDG==90000003)
                                                                   DecayMultyBaryon(quasH);
          else if(resQPDG==90002003||resQPDG==90003002) DecayAlphaBar(quasH); //DelEqu
          else if(resQPDG==90004002) DecayAlphaDiN(quasH);//Decay alph+2p(alph+2n IsStable)
          else if(resQPDG==90004004) DecayAlphaAlpha(quasH); // DelEqu
          else theQHadrons.push_back(quasH);    // Fill ResidQuasmHadron (del. equiv.)
		}
        else
        {
          delete quasH;                         // Delete "Quasmon Hadron"
#ifdef pdebug
          G4cerr<<"***G4QEnv::CheckGS: Decay in Fragm+ResQ did not succeeded"<<G4endl;
#endif
          return false;
        }
	  }
    }
    else                                        // *** Try Last+Previous CORRECTION ***
    {
#ifdef pdebug
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
#ifdef pdebug
	    G4cout<<"G4QE::CGS:M="<<tmpTLV<<tmpTLV.m()<<">t+h+pM="<<tQMa+hadrMa+prevMa<<G4endl;
#endif
        if(tmpTLV.m()>tQMa+hadrMa+prevMa)
        {
          G4LorentzVector nuc4M = G4LorentzVector(0.,0.,0.,tQMa);// 4-Mom of ResidNucAtRest
          G4QHadron* nucH = new G4QHadron(reTQC, nuc4M);
          if(!G4QHadron(tmpTLV).DecayIn3(hadr4M,prev4M,nuc4M))
          {
            delete nucH;                        // Delete "Residual Nucleus Hadron"
            G4cerr<<"---Worning---G4QE::CGS:DecayIn ResNuc+LastH+PrevH Error"<<G4endl;
            return false;
          }
          else
          {
            theLast->Set4Momentum(hadr4M);
            thePrev->Set4Momentum(prev4M);
            nucH->Set4Momentum(nuc4M);
#ifdef pdebug
	        G4cout<<"G4QE::CGS:**SUCCESS**>CHECK,D4M="<<tmpTLV-hadr4M-prev4M-nuc4M<<G4endl;
#endif
            if(totPDG==92000000||totPDG==90002000||totPDG==90000002) DecayDibaryon(nucH);
            else if(totPDG==93000000||totPDG==90003000||totPDG==90000003)
                                                                  DecayMultyBaryon(nucH);
            else if(totPDG==90004002) DecayAlphaDiN(nucH);// Dec alph+2p (alph+2n isStable)
            else if(totPDG==90002003||totPDG==90003002) DecayAlphaBar(nucH); //DelEqu
            else if(totPDG==90004004) DecayAlphaAlpha(nucH); //DelEqu
            else theQHadrons.push_back(nucH);             // Fill ResidNuclHadron (del.eq.)
	      }
		}
        else
        {
#ifdef pdebug
		  G4cout<<"G4QE::CGS: P&L didn't help, nH="<<nOfOUT<<">2"<<",MRQ="<<resQMa<<G4endl;
#endif
		  if(nOfOUT>2)                   // Try to find the appropriate partner
          {
            G4int nphot=-1;
            for(G4int id=nOfOUT-3; id>=0; id--) // Search for a photon
  			{
               G4QHadron* curHadr = theQHadrons[id];
               G4int hPDG=curHadr->GetPDGCode();
               if(hPDG==22) nphot=id;
            }
            if(nphot>=0)                 // Photon is found, try to use it to resolve PANIC
			{
              G4QHadron* curHadr = theQHadrons[nphot];        // Pointer to the photon
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
#ifdef pdebug
                  G4cerr<<"***G4QEnv::CheckGS: Decay in Photon+ResQ tM="<<ttM<<G4endl;
#endif
                }
                else
                {
                  curHadr->Set4Momentum(ch4M);// Change 4M of Photon (reduced by decay)
                  rqH->Set4Momentum(quas4M);  // Fill 4M of the GS Residual Quasmon
                  G4double hB=rqH->GetBaryonNumber();
                  G4double hS=rqH->GetStrangeness();
                  G4double hC=rqH->GetCharge();
#ifdef pdebug
                  G4cout<<"G4QE::CGS: n="<<nphot<<",ph4M="<<ch4M<<"+RQ4M="<<quas4M<<G4endl;
#endif
                  if(totPDG==92000000||totPDG==90002000||totPDG==90000002)
                                                                  DecayDibaryon(rqH);
                  else if(hB>2 && (hB==hS || hB==hC || !hC&&!hS)) DecayMultyBaryon(rqH);
                  else if(totPDG==90004002) DecayAlphaDiN(rqH);//Dec al+2p (al+2n isStable)
                  else if(totPDG==90002003||totPDG==90003002) DecayAlphaBar(rqH); //DelEqu
                  else if(totPDG==90004004) DecayAlphaAlpha(rqH); //DelEqu
                  else if(hB>0 && (hC<0 || hC>hB)) DecayIsonucleus(rqH); //DelEqu
                  else theQHadrons.push_back(rqH); // Fill ResidQuasmHadron (del eq)
                  if(envPDG!=90000000) theEnvironment=G4QNucleus(G4QContent(0,0,0,0,0,0),
                                                            G4LorentzVector(0.,0.,0.,0.));
                  return true;
		        }
              } // End of the KINEMATIC CHECK FOR THE PHOTON if
            }
            // > Photons did not help, try to find an appropriate partner to join and decay
            G4int    reTBN=reTQC.GetBaryonNumber();    // Baryon number of theHadronicState
            G4int    reTCH=reTQC.GetCharge();          // Charge of theHadronicState
            G4bool isoN = reTCH-reTBN==1 || reTCH==-1; // UnavoidableIsonucleus (Delta)cond
            G4bool norN = reTCH<=reTBN || reTCH>=0;    // "Regular nucleus" condition
            G4double nnM=resSMa;               // Fake prototype of the NormalNucleusMass
            G4QContent ipiQC=pipQC;            // Prototype of QCont for the Residual Pion+
            G4QContent nnQC=reTQC-ipiQC;       // Prototype for the NormalNucleus (Delta++)
            G4int nnPDG=nnQC.GetSPDGCode();    // Prot. PDGCode of the ResidNormalNucleus
            if(nnPDG==10&&nnQC.GetBaryonNumber()>0) nnPDG=nnQC.GetZNSPDGCode();
            if(isoN)                           // Calculations for the Isonuclear Residual
			{
              if(reTCH<0)                      // "Delta-" isostate
              {
                ipiQC=pimQC;                   // Change QCont for the Residual Pion-
                nnQC=reTQC-ipiQC;              // Change QCont for theNormalNucleus(Delta-)
                nnPDG=nnQC.GetSPDGCode();      // Change PDGCode of theResidNormalNucleus
                if(nnPDG==10&&nnQC.GetBaryonNumber()>0) nnPDG=nnQC.GetZNSPDGCode();
			  }
              G4QPDGCode nnQPDG(nnPDG);
              if(nnPDG<80000000) nnM=nnQPDG.GetMass(); // Mass for the Fundamental Hadron
              else               nnM=nnQPDG.GetNuclMass(nnPDG); // Mass for the Nucleus
            }
            for(G4int hd=nOfOUT-3; hd>=0; hd--)// Try to use any hadron to resolve PANIC
  			{
              G4QHadron* curHadr = theQHadrons[hd];
              G4int chNF=curHadr->GetNFragments();
              G4int chCH=curHadr->GetCharge();
              G4int chBN=curHadr->GetBaryonNumber();
              if(!chNF)
              {
                G4LorentzVector ch4M=curHadr->Get4Momentum(); // 4Mom of the Current Hadron
                G4LorentzVector tt4M=ch4M+reTLV;// resSMa=GSMass of the ResidQuasmon(+Env)
                G4double chM=ch4M.m();          // Mass of the CurrentHadron from theOUTPUT
                G4double ttM=tt4M.m();          // TotalMass of CurHadr+ResidQuasm compaund
                if(isoN)                        // "1 Delta Isonucleus" case
		        {
                  if(nnM+mPi+chM<ttM)           // PANIC can be resolved with thisCurHadron
			      {
                    G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,nnM); // GSMass of RQ
                    G4QHadron* rqH = new G4QHadron(nnQC,quas4M); // Prototype OutputResidQ
                    G4LorentzVector ipi4M = G4LorentzVector(0.,0.,0.,mPi);// GSMass of Pion
                    G4QHadron* rpH = new G4QHadron(ipiQC,ipi4M); // Prototype ResidualPion
                    if(!G4QHadron(tt4M).DecayIn3(ch4M,ipi4M,quas4M))
                    {
                      delete rqH;               // Delete tmp "Residual Quasmon Hadron"
                      delete rpH;               // Delete tmp "Residual Pion"
#ifdef pdebug
                      G4cerr<<"***G4QE::CGS: DecayIn3 CurH+ResQ+Pion dM="<<ttM-chM<<G4endl;
#endif
                    }
                    else
                    {
                      rpH->Set4Momentum(ipi4M);   // Change 4M of the Residual Pion(@@Why?)
                      curHadr->Set4Momentum(ch4M);// Change 4M of the Current Hadron
                      rqH->Set4Momentum(quas4M);  // Fill 4M of the GS Residual Quasmon
                      theQHadrons.push_back(rpH); // Fill Resid Pion (delete equivalent)
                      G4double hB=rqH->GetBaryonNumber();
                      G4double hS=rqH->GetStrangeness();
                      G4double hC=rqH->GetCharge();
#ifdef pdebug
                      G4cout<<"G4QE::CGS:#"<<hd<<",h="<<ch4M<<"+rq="<<quas4M<<"+pi="<<ipi4M
                            <<G4endl;
#endif
                      if(nnPDG==92000000||nnPDG==90002000||nnPDG==90000002)
                                                                     DecayDibaryon(rqH);
                      else if(hB>2 && (hB==hS || hB==hC || !hC&&!hS))DecayMultyBaryon(rqH);
                      else if(nnPDG==90004002) DecayAlphaDiN(rqH);//Dec al+2p(al+2n isStab)
                      else if(nnPDG==90002003 || nnPDG==90003002) DecayAlphaBar(rqH);//DelE
                      else if(nnPDG==90004004) DecayAlphaAlpha(rqH); // Delete Equivalent
                      else if(hB>0 && (hC<0 || hC>hB)) DecayIsonucleus(rqH); //DelEqu
                      else theQHadrons.push_back(rqH); // Fill ResidQuasmHadron (del eq.)
                      if(envPDG!=90000000)theEnvironment=G4QNucleus(G4QContent(0,0,0,0,0,0)
                                                        ,G4LorentzVector(0.,0.,0.,0.));
                      return true;
		            }
                  }
                  if(reTCH<0&&chCH>0 || reTCH+chCH<=reTBN+chBN) // Isoexchange can help
			      {
                    G4QContent chQC=curHadr->GetQC(); // QuarkCont of the CurrentHadron
                    if(reTCH<0) chQC+=pimQC;          // Add the negativPion QC to CurHadr
                    else        chQC+=pipQC;          // Add the positivePion QC to CurHadr
                    G4QPDGCode chQPDG=G4QPDGCode(chQC.GetSPDGCode());// New QPDG of CurHadr
                    chM=chQPDG.GetMass();             // New Mass of the CurHadron
                    if(nnM+chM<ttM)
                    {
                      G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,nnM);//GSMass of RQ
                      G4QHadron* rqH = new G4QHadron(nnQC,quas4M);// Prototype OutputResidQ
                      if(!G4QHadron(tt4M).DecayIn2(ch4M,quas4M))
                      {
                        delete rqH;               // Delete tmp "Residual Quasmon Hadron"
#ifdef pdebug
                        G4cerr<<"***G4QE::CGS: DecayIn2 CurH+ResQ d="<<ttM-chM-nnM<<G4endl;
#endif
                      }
                      else
                      {
                        curHadr->Set4Momentum(ch4M);// Change 4M of the Current Hadron
                        curHadr->SetQPDG(chQPDG);   // Change QPDG of the Current Hadron
                        rqH->Set4Momentum(quas4M);  // Fill 4M of the GS Residual Quasmon
                        G4double hB=rqH->GetBaryonNumber();
                        G4double hS=rqH->GetStrangeness();
                        G4double hC=rqH->GetCharge();
#ifdef pdebug
                        G4cout<<"G4QE::CGS:#"<<hd<<",h="<<ch4M<<"+rq="<<quas4M<<G4endl;
#endif
                        if(nnPDG==92000000||nnPDG==90002000||nnPDG==90000002)
                                                                       DecayDibaryon(rqH);
                        else if(hB>2 && (hB==hS||hB==hC || !hC&&!hS))DecayMultyBaryon(rqH);
                        else if(nnPDG==90004002) DecayAlphaDiN(rqH);//Dec a+2p(a+2n isStab)
                        else if(nnPDG==90002003 || nnPDG==90003002)DecayAlphaBar(rqH); //DE
                        else if(nnPDG==90004004) DecayAlphaAlpha(rqH); // Delete Equivalent
                        else if(hB>0 && (hC<0 || hC>hB)) DecayIsonucleus(rqH); //DelEqu
                        else theQHadrons.push_back(rqH); // Fill ResidQuasmHadron (del eq.)
                        if(envPDG!=90000000) theEnvironment=
                          G4QNucleus(G4QContent(0,0,0,0,0,0),G4LorentzVector(0.,0.,0.,0.));
                        return true;
		              }
                    }
                  }
                } // End of the KINEMATIC CHECK FOR THE CURRENT HADRON if (Isonuclear Case)
                else if(norN)
		        {
                  if(resSMa+chM<ttM)          // PANIC can be resolved with this CurHadron
			      {
                    G4LorentzVector quas4M = G4LorentzVector(0.,0.,0.,resSMa);// GSMass ofQ
                    G4QHadron* rqH = new G4QHadron(reTQC,quas4M);// Prototype OutputResidQ
                    if(!G4QHadron(tt4M).DecayIn2(ch4M,quas4M))
                    {
                      delete rqH;                   // Delete tmp "Residual Quasmon Hadron"
#ifdef pdebug
                      G4cerr<<"***G4QE::CheckGS: Decay in CurH+ResQ dM="<<ttM-chM<<G4endl;
#endif
                    }
                    else
                    {
                      curHadr->Set4Momentum(ch4M);  // Change 4M of the Current Hadron
                      rqH->Set4Momentum(quas4M);    // Fill 4M of the GS Residual Quasmon
                      G4double hB=rqH->GetBaryonNumber();
                      G4double hS=rqH->GetStrangeness();
                      G4double hC=rqH->GetCharge();
#ifdef pdebug
                      G4cout<<"G4QEnv::CheckGS:#"<<hd<<",ch4M="<<curHadr->GetPDGCode()
                            <<ch4M<<" + ResQ4M="<<totPDG<<quas4M<<G4endl;
#endif
                      if(totPDG==92000000||totPDG==90002000||totPDG==90000002)
                                                                     DecayDibaryon(rqH);
                      else if(hB>2 && (hB==hS || hB==hC || !hC&&!hS))DecayMultyBaryon(rqH);
                      else if(totPDG==90004002) DecayAlphaDiN(rqH);//Dec a+2p(a+2n isStab.)
                      else if(totPDG==90002003||totPDG==90003002) DecayAlphaBar(rqH);//DelE
                      else if(totPDG==90004004) DecayAlphaAlpha(rqH); //DelEqu
                      else if(hB>0 && (hC<0 || hC>hB)) DecayIsonucleus(rqH); //DelEqu
                      else theQHadrons.push_back(rqH); // Fill ResidQuasmHadron (del. eq.)
                      if(envPDG!=90000000)theEnvironment=G4QNucleus(G4QContent(0,0,0,0,0,0)
                                                        ,G4LorentzVector(0.,0.,0.,0.));
                      return true;
		            }
                  } // End of the KINEMATIC CHECK FOR THE CURRENT HADRON if (NormNuclCase)
                }
              } // End of the NumberOfFragments=0 (NOT DECAYED PARTICLE) if
			} // End of the LOOP over hadrons and all attempts to resolve PANIC
#ifdef pdebug
            G4cout<<"G4QEnv::CheckGS: *** Any hadron from the OUTPUT did not help"<<G4endl;
#endif
            return false;
		  } // End of the POSSIBILITY OF MORE THAN L&P CORRECTION if
		  else return false; // If L&P is not possible 
        } // End of the KINEMATIC LIMIT FOR THE L&P CORRECTION if/els
	  } // End of the POSSIBILITY OF PREV+LAST (OR MORE) CORRECTION if
      else return false;
    } // End of the CORRECTION WITH THE LAST if/else
  } // End of the CORRECTION IS POSSIBLE if
  else return false;                     // Correction can not be done
  return true;                           // If correction was done successfully
} // End of "CheckGroundState"

// Try to decay the Total Residual Nucleus in Environ+Quasmon
void G4QEnvironment::CopyAndDeleteHadronVector(G4QHadronVector* HV)
{ // ==============================================================
  G4int nHadrons = HV->size();
  if(nHadrons)
  {
    for(G4int ih=0; ih<nHadrons; ih++) // LOOP over output QHadrons
    {
      G4QHadron* inH = HV->operator[](ih);       // Pointer to the i-th QHadron
      G4int hNF  = inH->GetNFragments();            // A#of secondary fragments
      if(!hNF)                                      // Fill only final hadrons
      {
#ifdef pdebug
	     G4cout<<"G4QEnv::Copy&DeleteHV:#"<<ih<<", hPDG="<<inH->GetPDGCode()<<G4endl;
#endif
        G4QHadron* curH = new G4QHadron(inH);      // will be deleted with allQHadronVector
        theQHadrons.push_back(curH);               // Fill hadron-copy (delete equivalent)
      }
    }
    std::for_each(HV->begin(), HV->end(), DeleteQHadron()); // Delete instances
    HV->clear();                                            // Delete pointers
  }
#ifdef pdebug
  else G4cerr<<"***G4QEnv::Kopy&DelHV: No hadrons in the QHadronVector"<<G4endl;
#endif
  delete HV;                                       // Delete the inputQHadronVector
} // End of "CopyAndDeleteHadronVector"

// Try to decay the Total Residual Nucleus in Environ+Quasmon
G4bool G4QEnvironment::DecayInEnvQ(G4Quasmon* quasm)
{ //   =============================================
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
	//#ifdef pdebug
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
        G4cerr<<"---Worning---G4Q::DecInEnvQ:Decay in Environment+ResidualQuasmon"<<G4endl;
        return false;
      }
      else
      {
        quasH->Set4Momentum(quas4M);
        if(resQPDG==92000000||resQPDG==90002000||resQPDG==90000002) DecayDibaryon(quasH);
        else if(resQPDG==93000000||resQPDG==90003000||resQPDG==90000003)
                                                                   DecayMultyBaryon(quasH);
        else if(resQPDG==90004002) DecayAlphaDiN(quasH);//Decay alph+2p (alph+2n is stable)
        else if(resQPDG==90002003||resQPDG==90003002) DecayAlphaBar(quasH); //DelEqu
        else if(resQPDG==90004004) DecayAlphaAlpha(quasH); //DelEqu
        else theQHadrons.push_back(quasH);        // Fill ResidQuasmHadron (del.equivalent)
        envaH->Set4Momentum(enva4M);
        if(envPDG==92000000||envPDG==90002000||envPDG==90000002) DecayDibaryon(envaH);
        else if(envPDG==93000000||envPDG==90003000||envPDG==90000003)
                                                                   DecayMultyBaryon(envaH);
        else if(envPDG==90004002) DecayAlphaDiN(envaH);// Decay al+2p (a+2n is stable)
        else if(envPDG==90002003||envPDG==90003002) DecayAlphaBar(envaH); //DelEqu
        else if(envPDG==90004004) DecayAlphaAlpha(envaH); //DelEqu
        else theQHadrons.push_back(envaH);        // Fill Environment Hadron (del.Equiv.)
	  }
    }
    else return false;
  }
  else return false;                              // => "Environment is vacuum" case
  return true;
} // End of "DecayInEnvQ"
