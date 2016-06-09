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
// $Id: G4QCaptureAtRest.cc,v 1.2 2010-06-25 09:46:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QCaptureAtRest class -----------------
//                 by Mikhail Kossov, December 2003.
// G4QCaptureAtRest class of the CHIPS Simulation Branch in GEANT4
// ---------------------------------------------------------------
// ****************************************************************************************
// This Header is a part of the CHIPS physics package (author: M. Kosov)
// ****************************************************************************************
// Short Description: This is a universal process for nuclear capture
// (including annihilation) of all negative particles (negative hadrons,
// negative leptons: mu- & tau-). It can be used for the cold neutron
// capture, but somebody should decide what is the probability (defined
// by the capture cross-section and atomic material properties) to switch
// the cold neutron to the at-rest neutron. - M.K.2009.
// ----------------------------------------------------------------------


//#define debug
//#define pdebug
//#define tdebug

#include "G4QCaptureAtRest.hh"

G4QCaptureAtRest::G4QCaptureAtRest(const G4String& processName)
  : G4VRestProcess(processName, fHadronic), Time(0.), EnergyDeposition(0.)
{
#ifdef debug
  G4cout<<"G4QCaptureAtRest::Constructor is called"<<G4endl;
#endif
  if (verboseLevel>0) G4cout << GetProcessName() << " is created "<< G4endl;

  //G4QCHIPSWorld::Get()->GetParticles(nPartCWorld); // Create CHIPS World (234 part. max)
  G4QNucleus::SetParameters(freeNuc,freeDib,clustProb,mediRatio); // Clusterization param's
  G4Quasmon::SetParameters(Temperature,SSin2Gluons,EtaEtaprime);  // Hadronic parameters
  G4QEnvironment::SetParameters(SolidAngle); // SolAngle of pbar-A secondary mesons capture
}

G4bool   G4QCaptureAtRest::manualFlag=false; // If false then standard parameters are used
G4double G4QCaptureAtRest::Temperature=180.; // Critical Temperature (sensitive at High En)
G4double G4QCaptureAtRest::SSin2Gluons=0.3;  // Supression of s-quarks (in respect to u&d)
G4double G4QCaptureAtRest::EtaEtaprime=0.3;  // Supression of eta mesons (gg->qq/3g->qq)
G4double G4QCaptureAtRest::freeNuc=0.5;      // Percentage of free nucleons on the surface
G4double G4QCaptureAtRest::freeDib=0.05;     // Percentage of free diBaryons on the surface
G4double G4QCaptureAtRest::clustProb=5.;     // Nuclear clusterization parameter
G4double G4QCaptureAtRest::mediRatio=10.;    // medium/vacuum hadronization ratio
G4int    G4QCaptureAtRest::nPartCWorld=152;  // The#of particles initialized in CHIPS World
G4double G4QCaptureAtRest::SolidAngle=0.5;   // Part of Solid Angle to capture (@@A-dep.)
G4bool   G4QCaptureAtRest::EnergyFlux=false; // Flag for Energy Flux use (not MultyQuasmon)
G4double G4QCaptureAtRest::PiPrThresh=141.4; // Pion Production Threshold for gammas
G4double G4QCaptureAtRest::M2ShiftVir=20000.;// Shift for M2=-Q2=m_pi^2 of the virtualGamma
G4double G4QCaptureAtRest::DiNuclMass=1880.; // DoubleNucleon Mass for VirtualNormalization

void G4QCaptureAtRest::SetManual()   {manualFlag=true;}
void G4QCaptureAtRest::SetStandard() {manualFlag=false;}

// Fill the private parameters
void G4QCaptureAtRest::SetParameters(G4double temper, G4double ssin2g, G4double etaetap,
                                     G4double fN, G4double fD, G4double cP, G4double mR,
                                     G4int nParCW, G4double solAn, G4bool efFlag,
                                     G4double piThresh, G4double mpisq, G4double dinum)
{//  =============================================================================
  Temperature=temper;
  SSin2Gluons=ssin2g;
  EtaEtaprime=etaetap;
  freeNuc=fN;
  freeDib=fD;
  clustProb=cP;
  mediRatio=mR;
  nPartCWorld = nParCW;
  EnergyFlux=efFlag;
  SolidAngle=solAn;
  PiPrThresh=piThresh;
  M2ShiftVir=mpisq;
  DiNuclMass=dinum;
  G4QCHIPSWorld::Get()->GetParticles(nPartCWorld); // Create CHIPS World with 234 particles
  G4QNucleus::SetParameters(freeNuc,freeDib,clustProb,mediRatio); // Clusterization param's
  G4Quasmon::SetParameters(Temperature,SSin2Gluons,EtaEtaprime);  // Hadronic parameters
  G4QEnvironment::SetParameters(SolidAngle); // SolAngle of pbar-A secondary mesons capture
}

// Destructor

G4QCaptureAtRest::~G4QCaptureAtRest()
{}

G4LorentzVector G4QCaptureAtRest::GetEnegryMomentumConservation()
{
  return EnMomConservation;
}

G4int G4QCaptureAtRest::GetNumberOfNeutronsInTarget()
{
  return nOfNeutrons;
}

G4bool G4QCaptureAtRest::IsApplicable(const G4ParticleDefinition& particle) 
{
  if      (particle == *(    G4PionMinus::PionMinus()    )) return true;
  else if (particle == *(    G4KaonMinus::KaonMinus()    )) return true;
  else if (particle == *(   G4AntiProton::AntiProton()   )) return true;
  else if (particle == *(    G4MuonMinus::MuonMinus()    )) return true;
  else if (particle == *(     G4TauMinus::TauMinus()     )) return true;
  else if (particle == *(   G4SigmaMinus::SigmaMinus()   )) return true;
  else if (particle == *(      G4XiMinus::XiMinus()      )) return true;
  else if (particle == *(   G4OmegaMinus::OmegaMinus()   )) return true;
  else if (particle == *(      G4Neutron::Neutron()      )) return true;
  else if (particle == *(  G4AntiNeutron::AntiNeutron()  )) return true;
  else if (particle == *(G4AntiSigmaPlus::AntiSigmaPlus())) return true;
#ifdef debug
  G4cout<<"***G4QCaptureAtRest::IsApplicable: PDG="<<particle.GetPDGEncoding()<<G4endl;
#endif
  return false;
}

G4VParticleChange* G4QCaptureAtRest::AtRestDoIt(const G4Track& track, const G4Step& step)
{
  static const G4double mNeut = G4QPDGCode(2112).GetMass();
  static const G4double mNeut2= mNeut*mNeut;
  static const G4double dmNeut= mNeut+mNeut;
  static const G4double mProt = G4QPDGCode(2212).GetMass();
  static const G4double dmProt= mProt+mProt;
  static const G4double mPi0  = G4QPDGCode(111).GetMass();
  static const G4double mDeut = G4QPDGCode(2112).GetNuclMass(1,1,0);
  static const G4double mAlph = G4QPDGCode(2112).GetNuclMass(2,2,0);
  //static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mMu  = G4QPDGCode(13).GetMass();
  //static const G4double mMu2 = mMu*mMu;
  //static const G4double dmMu = mMu+mMu;
  // The best parameters for mu->e+nu+nu decay
  static const G4double b1=2.784;
  static const G4double ga=.000015;
  static const G4double rb=1./b1;
  //static const G4double mTau = G4QPDGCode(15).GetMass();
  static const G4double mEl  = G4QPDGCode(11).GetMass();
  static const G4double mEl2 = mEl*mEl;
  //-------------------------------------------------------------------------------------
  static G4bool CWinit = true;                       // CHIPS Warld needs to be initted
  if(CWinit)
  {
    CWinit=false;
    G4QCHIPSWorld::Get()->GetParticles(nPartCWorld); // Create CHIPS World (234 part.max)
  }
  //-------------------------------------------------------------------------------------
  const G4DynamicParticle* stoppedHadron = track.GetDynamicParticle();
  const G4ParticleDefinition* particle=stoppedHadron->GetDefinition();
  Time=0.;
  EnergyDeposition=0.;
#ifdef debug
  G4cout<<"G4QCaptureAtRest::AtRestDoIt is called,EnDeposition="<<EnergyDeposition<<G4endl;
#endif
  if (!IsApplicable(*particle))  // Check applicability
  {
    G4cerr<<"G4QCaptureAtRest::AtRestDoIt: Only mu-,pi-,K-,S-,X-,O-,aP,aN,aS+."<< G4endl;
    return 0;
  }
  const G4Material* material = track.GetMaterial();      // Get the current material
  G4int Z=0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int i=0;
  G4double sum=0.;
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QCaptureAtRest::AtRestDoIt: "<<nE<<" elements in the material."<<G4endl;
#endif

  //VI===== protection against super high energy - start of a loop
  G4double primEnergy = 0.0; 
  G4double secEnergy = 0.0;
  G4bool productOK = true;
  const G4double limEnergy = 2*GeV;
  G4int counter = 0;
  do {
  secEnergy = 0.0;
  ++counter;
  G4int countA = particle->GetBaryonNumber();
  //VI===== 


  G4int projPDG=0;                           // PDG Code prototype for the captured hadron
  if      (particle ==     G4MuonMinus::MuonMinus()    ) projPDG=   13;
  else if (particle ==      G4TauMinus::TauMinus()     ) projPDG=   15; // @@AtomicRad?
  else if (particle ==     G4PionMinus::PionMinus()    ) projPDG= -211; // @@AtomicRad?
  else if (particle ==     G4KaonMinus::KaonMinus()    ) projPDG= -321;
  else if (particle ==    G4AntiProton::AntiProton()   ) projPDG=-2212;
  else if (particle ==    G4SigmaMinus::SigmaMinus()   ) projPDG= 3112;
  else if (particle ==       G4XiMinus::XiMinus()      ) projPDG= 3312;
  else if (particle ==    G4OmegaMinus::OmegaMinus()   ) projPDG= 3334;
  else if (particle ==       G4Neutron::Neutron()      ) projPDG= 2112;
  else if (particle ==   G4AntiNeutron::AntiNeutron()  ) projPDG=-2112;
  else if (particle == G4AntiSigmaPlus::AntiSigmaPlus()) projPDG=-3222;
#ifdef debug
  G4cout<<"G4QCaptureAtRest::AtRestDoIt: projPDG="<<projPDG<<G4endl;
#endif
  if(!projPDG)
  {
    G4cerr<<"---Worning---G4QCaptureAtRest::AtRestDoIt: Undefined captured hadron"<<G4endl;
    return 0;
  }
  std::vector<G4double> sumfra;
  for(i=0; i<nE; ++i)
  {
    G4double frac=material->GetFractionVector()[i];
    if(projPDG==13||projPDG==15) frac*=(*theElementVector)[i]->GetZ();
    sum+=frac;
    sumfra.push_back(sum);             // remember the summation steps
  }
  G4double rnd = sum*G4UniformRand();
  for(i=0; i<nE; ++i) if (rnd<sumfra[i]) break;
  G4Element* pElement=(*theElementVector)[i];
  Z=static_cast<G4int>(pElement->GetZ());
  if(Z<=0)
  {
    G4cerr<<"---Worning---G4QCaptureAtRest::AtRestDoIt:Element with Z="<<Z<< G4endl;
    if(Z<0) return 0;
  }
  G4QIsotope* Isotopes = G4QIsotope::Get(); // Pointer to the G4QIsotopes singleton
  G4int N = Z;
  G4int isoSize=0;                         // The default for the isoVectorLength is 0
  G4int indEl=0;                           // Index of non-natural element or 0 (default)
  G4IsotopeVector* isoVector=pElement->GetIsotopeVector();
  if(isoVector) isoSize=isoVector->size(); // Get real size of the isotopeVector if exists
#ifdef debug
  G4cout<<"G4QCaptureAtRest::AtRestDoIt: isovectorLength="<<isoSize<<G4endl;
#endif
  if(isoSize)                              // The Element has not trivial abumdance set
  {
    indEl=pElement->GetIndex()+1;          // Index of the non-trivial element is an order
#ifdef debug
    G4cout<<"G4QCapAR::GetMFP: iE="<<indEl<<", def="<<Isotopes->IsDefined(Z,indEl)<<G4endl;
#endif
    if(!Isotopes->IsDefined(Z,indEl))      // This index is not defined for this Z: define
    {
      std::vector<std::pair<G4int,G4double>*>* newAbund =
                                               new std::vector<std::pair<G4int,G4double>*>;
      G4double* abuVector=pElement->GetRelativeAbundanceVector();
      for(G4int j=0; j<isoSize; j++)
      {
        N=pElement->GetIsotope(j)->GetN()-Z;
        if(pElement->GetIsotope(j)->GetZ()!=Z) G4cerr<<"*G4QCaptureAtRest::AtRestDoIt: Z="
                                        <<pElement->GetIsotope(j)->GetZ()<<"#"<<Z<<G4endl;
        G4double abund=abuVector[j];
        std::pair<G4int,G4double>* pr= new std::pair<G4int,G4double>(N,abund);
#ifdef debug
        G4cout<<"G4QCaptureAtRest::AtRestDoIt:pair#="<<j<<", N="<<N<<",ab="<<abund<<G4endl;
#endif
        newAbund->push_back(pr);
      }
#ifdef debug
      G4cout<<"G4QCaptureAtRest::AtRestDoIt: pairVectorLength="<<newAbund->size()<<G4endl;
#endif
      indEl=Isotopes->InitElement(Z,indEl,newAbund); // redefinie newInd (if exists)
      for(G4int k=0; k<isoSize; k++) delete (*newAbund)[k];
      delete newAbund;
    }
    // @@ ^^^^^^^^^^ End of the temporary solution ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    N = Isotopes->GetNeutrons(Z,indEl);
  }
  else  N = Isotopes->GetNeutrons(Z);
  nOfNeutrons=N;                                       // Remember it for energy-mom. check
  G4double dd=0.025;
  G4double am=Z+N;
  G4double sr=std::sqrt(am);
  G4double dsr=0.01*(sr+sr);
  if(dsr<dd)dsr=dd;
  if(manualFlag) G4QNucleus::SetParameters(freeNuc,freeDib,clustProb,mediRatio); // ManualP
  else if(projPDG==-2212) G4QNucleus::SetParameters(1.-dsr-dsr,dd+dd,5.,10.);//aP ClustPars
  else if(projPDG==-211)  G4QNucleus::SetParameters(.67-dsr,.32-dsr,5.,9.);//Pi- ClustPars
#ifdef debug
  G4cout<<"G4QCaptureAtRest::AtRestDoIt: N="<<N<<" for element with Z="<<Z<<G4endl;
#endif
  if(N<0)
  {
    G4cerr<<"---Worning---G4QCaptureAtRest::AtRestDoIt:Element with N="<<N<< G4endl;
    return 0;
  }
  G4bool lepChan=true;
  if(projPDG==13)
  {
     CalculateEnergyDepositionOfMuCapture(Z);// Fills the "EnergyDeposition" value
     lepChan=RandomizeMuDecayOrCapture(Z, N);// Fills the "Time" value
  }
  else if(projPDG==15)
  {
     CalculateEnergyDepositionOfTauCapture(Z);// Fills the "EnergyDeposition" value
     lepChan=RandomizeTauDecayOrCapture(Z,N);// Fills the "Time" value
  }
  G4double mp=G4QPDGCode(projPDG).GetMass(); // Mass of the captured hadron
  G4int targPDG=90000000+Z*1000+N;           // PDG Code of the target nucleus

  //VI===== save primary energy
  primEnergy = mp+G4QPDGCode(targPDG).GetMass();
  
  G4QHadronVector* output=new G4QHadronVector; // Prototype of the output G4QHadronVector
#ifdef debug
  G4cout<<"G4QCaptureAtRest::AtRestDoIt: projPDG="<<projPDG<<", targPDG="<<targPDG<<G4endl;
#endif
  G4double      weight    = track.GetWeight();
  G4double      localtime = track.GetGlobalTime();
  G4ThreeVector position  = track.GetPosition();
#ifdef debug
  G4cout<<"G4QCaptureAtRest::AtRestDoIt: t="<<localtime<<", p="<<position<<G4endl;
#endif
  G4TouchableHandle trTouchable = track.GetTouchableHandle();
#ifdef debug
  G4cout<<"G4QCaptureAtRest::AtRestDoIt: touch="<<trTouchable<<G4endl;
#endif
  localtime += Time;
  G4bool neutronElastic = false;             // Flag of elastic neutro-nucleeus Scattering
  G4bool chargExElastic = false;             // Flag of charge exchange Quasi-Elastic Scat.
  G4LorentzVector proj4M;
  G4double mAT=0.;                           // Prototype of mass of the GS target nucleus
  G4double totNE=0.;                         // Prototype of total neutron energy
  G4double mAP=0.;                           // Prototype of M_GSCompoundNucleus-proton
  if(projPDG==2112)
  {
    proj4M=stoppedHadron->Get4Momentum();
    totNE=proj4M.e();
    // For low energy neutrons the threshold of the capture must be checked
    G4QPDGCode QPDGbase;
    mAT=QPDGbase.GetNuclMass(Z,N,0);          // mass of the GS target nucleus
    G4double mAR=std::sqrt(mAT*(mAT+totNE+totNE)+mNeut2); // Real compound mass
    G4double mAN=QPDGbase.GetNuclMass(Z,N+1,0); // mass of the GS compound nucleus
    mAP=QPDGbase.GetNuclMass(Z-1,N+1,0);      // M_GSCompoundNucleus-proton
    G4double mAA=1000000.; // Default (light nuclei) mass of the GSCompoundNucleus-alpha
    if(Z>=2 && N>=1) mAA=QPDGbase.GetNuclMass(Z-2,N-1,0); // mass of GSCompNucleus-alpha
    G4double eProt=mAR-mAP-mProt;             // Possible kin Enrgy of residual proton
    if(mAR<mAN && eProt>0.)                   // Compound is impossible but ChEx's possible
    {
#ifdef debug
      G4cout<<"G4QCaptureAtRest::AtRestDoIt: n-Capture isn't possible mC="<<mAR<<" < mGS="
            <<mAN<<", Ep="<<eProt<<G4endl;
#endif
      G4double eNeut=totNE-mNeut;             // Kinetic energy of the projectile neutron
      if(eNeut<0.) eNeut=0.;                  // This is just an accuracy correction
      if(eNeut<.0001) chargExElastic=true;    // neutron is too soft -> charge Exchange
      else
      {
        G4double probP=std::sqrt(eProt*(dmProt+eProt));
        G4double probN=std::sqrt(eNeut*(dmNeut+eNeut));
        if((probN+probP)*G4UniformRand()<probN) neutronElastic=true; // nElastic Scattering
        else chargExElastic=true; // proton's phase space is bigger -> chargeExchange
      }
    }
    else if(mAR<=mAN||(mAR<=mAP+mProt&&mAR<=mAA+mAlph)) // Impossible to radiate n or Alpha
    {
#ifdef debug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt: n-Capture only elastic is possible"<<G4endl;
#endif
      neutronElastic=true; // nElaScattering
    }
#ifdef debug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt: n-Capture El="<<neutronElastic<<", Ex="
          <<chargExElastic<<G4endl;
#endif
  }
  G4int           nuPDG=14;                  // Prototype for weak decay
  if(projPDG==15) nuPDG=16;
#ifdef debug
  G4int CV=0;
  G4cout<<"G4QCaptureAtRest::AtRestDoIt:DecayIf is reached CV="<<CV<<G4endl;
#endif
  if(projPDG==2112 && neutronElastic)        // Elastic scattering of low energy neutron
  {
#ifdef debug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt:nElast, 4M="<<proj4M<<",Z="<<Z<<",N="<<N<<G4endl;
#endif
    G4LorentzVector a4Mom(0.,0.,0.,mAT);     // 4-momentum of the target nucleus
    G4LorentzVector totLV=proj4M+a4Mom;      // 4-momentum of the compound system
    G4LorentzVector n4Mom(0.,0.,0.,mNeut);   // mass of the secondary neutron
    if(!G4QHadron(totLV).DecayIn2(a4Mom,n4Mom))
    {
      G4cerr<<"---Worning---G4QCaptureAtRest::AtRestDoIt:n+A=>n+A,Z="<<Z<<",N="<<N<<G4endl;
      return 0;
    }
    G4QHadron* secnuc = new G4QHadron(targPDG,a4Mom); // Create recoil nucleus
    output->push_back(secnuc);               // Fill recoil nucleus to the output
    G4QHadron* neutron = new G4QHadron(2112,n4Mom);    // Create Hadron for the Neutron
    output->push_back(neutron);              // Fill the neutron to the output
#ifdef debug
    CV=27;
    G4cout<<"G4QCaptureAtRest::AtRestDoIt:ElasN="<<n4Mom<<",A="<<a4Mom<<",CV="<<CV<<G4endl;
#endif
  }
  else if(projPDG==2112 && chargExElastic)   // ChargeEx from neutron to proton: (n,p) reac
  {
#ifdef debug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt:npChEx, 4M="<<proj4M<<",Z="<<Z<<",N="<<N<<G4endl;
#endif
    G4LorentzVector totLV(0.,0.,0.,mAT);     // 4-momentum of the target nucleus
    totLV+=proj4M;                           // 4-momentum of the compound system
    G4LorentzVector a4Mom(0.,0.,0.,mAP);     // 4-momentum of the residual to (n,p)
    G4LorentzVector p4Mom(0.,0.,0.,mProt);   // mass of the secondary proton
    if(!G4QHadron(totLV).DecayIn2(a4Mom,p4Mom))
    {
      G4cerr<<"---Worning---G4QCaptureAtRest::AtRestDoIt:n+A=p+A',Z="<<Z<<",N="<<N<<G4endl;
      return 0;
    }
    G4QHadron* secnuc = new G4QHadron(targPDG-999,a4Mom); // Create recharged nucleus
    output->push_back(secnuc);               // Fill recharged nucleus to the output
    G4QHadron* proton = new G4QHadron(2212,p4Mom); // Create Hadron for the Proton
    output->push_back(proton) ;              // Fill the proton to the output
#ifdef debug
    CV=21;
    G4cout<<"G4QCaptureAtRest::AtRestDoIt:ChExP="<<p4Mom<<",A="<<a4Mom<<",CV="<<CV<<G4endl;
#endif
  }
  else if(projPDG==-211 && targPDG==90001000)// Use Panofsky Ratio for (p+pi-) system decay
  {                                          // (p+pi-=>n+pi0)/p+pi-=>n+gamma) = 3/2
#ifdef debug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt: Panofsky targPDG="<<targPDG<<G4endl;
#endif
    G4LorentzVector totLV(0.,0.,0.,mp+mProt);// 4-momentum of the compound system
    G4int pigamPDG=111;                      // Prototype is for pi0
    G4double pigamM=mPi0;
    if(G4UniformRand()>0.6)
    {
      pigamPDG=22;
      pigamM=0.;
    }
    G4LorentzVector g4Mom(0.,0.,0.,pigamM);  // mass of the photon/Pi0
    G4LorentzVector n4Mom(0.,0.,0.,mNeut);   // mass of the secondary neutron
    if(!G4QHadron(totLV).DecayIn2(g4Mom,n4Mom))
    {
      G4cerr<<"---Worning---G4QCaptureAtRest::AtRestDoIt: H1+(pi-)=>n+"<<pigamPDG<<G4endl;
      return 0;
    }
    G4QHadron* pigamma = new G4QHadron(pigamPDG,g4Mom); // Creation Hadron for Pi0/Gamma
    output->push_back(pigamma);              // Fill pi0 or Gamma in the output
    G4QHadron* neutron = new G4QHadron(2112,n4Mom);    // Create Hadron for the Neutron
    output->push_back(neutron);              // Fill the neutron to the output
  }
  // @@ For pi-,d reactions one can use just nn dedcay (see above) ?
  // @@ For K- Capture the quasifree n+Lambda+(A-n-p) reaction can be applyed as well...
  else if(projPDG==-211 && G4UniformRand()>1.&& Z>0&&N>0) // @@Quasi-Free PiCapture => tune
  {
    G4double mt=G4QPDGCode(targPDG).GetMass();// Mass of the target Nucleus
    G4LorentzVector totLV(0.,0.,0.,mp+mt);   // 4-momentum of the (A+pi-) compound system
    if(Z==1 && N==1)                         // Quasi-Free process on Deuteron
    {
      G4LorentzVector f4Mom(0.,0.,0.,mNeut); // First neutron 
      G4LorentzVector s4Mom(0.,0.,0.,mNeut); // Second neutron
      if(!G4QHadron(totLV).DecayIn2(f4Mom,s4Mom))
      {
        G4cerr<<"---Worning---G4QCaptureAtRest::AtRestDoIt: H2+(pi-)=>n+n"<<G4endl;
        return 0;
      }
      G4QHadron* neutr1 = new G4QHadron(2112,f4Mom); // Create Hadron for the 1st Neutron
      output->push_back(neutr1);             // Fill pi0 or Gamma in the output
      G4QHadron* neutr2 = new G4QHadron(2112,s4Mom); // Create Hadron for the 2nd Neutron
      output->push_back(neutr2);             // Fill the neutron to the output
    }
    else
    {
      G4int rPDG=targPDG-1001;
      G4double mr=G4QPDGCode(rPDG).GetMass();// Mass of the residual Nucleus
      G4LorentzVector f4Mom(0.,0.,0.,mNeut); // First neutron 
      G4LorentzVector s4Mom(0.,0.,0.,mNeut); // Second neutron
      G4LorentzVector r4Mom(0.,0.,0.,mr);    // Residual nucleus
      if(!G4QHadron(totLV).DecayIn3(f4Mom,s4Mom,r4Mom))
      {
        G4cerr<<"---Worning---G4QCaptureAtRest::AtRestDoIt: A+(pi-)=>n+n+(A-n-p)"<<G4endl;
        return 0;
      }
      G4QHadron* neutr1 = new G4QHadron(2112,f4Mom); // Create Hadron for the 1st Neutron
      output->push_back(neutr1);             // Fill the first neutron in the output
      G4QHadron* neutr2 = new G4QHadron(2112,s4Mom); // Create Hadron for the 2nd Neutron
      output->push_back(neutr2);             // Fill the second neutron to the output
      G4QHadron* resnuc = new G4QHadron(rPDG,r4Mom); // Create Hadron for the ResidualNucl
      output->push_back(resnuc);             // Fill the Residual Nucleus to the output
    }
  }
  else if((projPDG==13||projPDG==15) && !lepChan)//Normal BoundLepton->e+nu+anti_nu_e decay
  {
    G4double mbm=mp-EnergyDeposition;        // Mass of the bounded muon
    G4LorentzVector totLV(0.,0.,0.,mbm);     // 4-momentum of bounded muon/tau
    if(projPDG==13) // now for tau it is only energy deposition, for mu this EMCascade
    {
#ifdef debug
      G4cout<<"G4QCaptureAtRest::AtRestDoIt: e+nu+nu decay 4M="<<totLV<<totLV.m()<<G4endl;
#endif
      G4double mt=G4QPDGCode(targPDG).GetMass();// Mass of the target Nucleus
      G4double ee=mEl+RandomizeDecayElectron(Z);// Randomized electron total energy
      totLV+=G4LorentzVector(0.,0.,0.,mt);      // 4-momentum of the compound system
      G4double mmt=totLV.m();                   // Total energy of the compound system
      if(ee>=mbm*(mmt-mbm/2)/mmt)
      {
        G4cout<<"-W-G4QCaptureAtRest::AtRestDoIt: Unrealistic E="<<ee<<", m="<<mMu<<G4endl;
        G4LorentzVector f4Mom(0.,0.,0.,mEl);    // Electron
        G4LorentzVector s4Mom(0.,0.,0.,mt);     // Quark-A
        if(!G4QHadron(totLV).DecayIn2(f4Mom,s4Mom))
        {
          G4cerr<<"---Worning---G4QCaptureAtRest::AtRestDoIt: (mubound+A)->mu+A"<<G4endl;
          return 0;
        }
        G4QHadron* electron = new G4QHadron(11,f4Mom); // Create Electron
        output->push_back(electron);            // Fill the electron to the output
        G4QHadron* resnuc = new G4QHadron(targPDG,s4Mom); // Create Residual Nucleus
        output->push_back(resnuc);              // Fill the Residual Nucleus to the output
      }
      else
      {
        G4double mr=std::sqrt(mmt*(mmt-ee-ee)+mEl2); // Mass of the M+nu+nu system
        if(mr<mt) mr=mt+.000001;                  // To be shure with the next decay
        G4LorentzVector f4Mom(0.,0.,0.,mEl);      // Electron 
        G4LorentzVector s4Mom(0.,0.,0.,mr);       // Nucleus+2neutrinos
        if(!G4QHadron(totLV).DecayIn2(f4Mom,s4Mom))
        {
          G4cerr<<"---Warning---G4QCaptureAtRest::AtRestDoIt: (mu+A)=>e+(A+2nu)"<<G4endl;
          return 0;
        }
#ifdef debug
        G4double fe=f4Mom.e();
        G4cout<<"G4QCaptureAtRest::AtRestDoIt: Ei="<<ee<<",Ef="<<fe<<",de="<<fe-ee<<G4endl;
#endif
        G4QHadron* electron = new G4QHadron(11,f4Mom); // Create Electron
        output->push_back(electron);             // ==> Fill the electron to the output
        G4LorentzVector r4Mom(0.,0.,0.,mt);      // residual nucleus
        G4LorentzVector n4Mom(0.,0.,0.,0.);      // muon neutrino
        G4LorentzVector a4Mom(0.,0.,0.,0.);      // electron anti-nutrino
        if(!G4QHadron(s4Mom).DecayIn3(r4Mom,n4Mom,a4Mom))
        {
          G4cerr<<"-Warning-G4QCaptureAtRest::AtRestDoIt: (A+2nu)=>A+NuMu+antiNuE"<<G4endl;
          return 0;
        }
#ifdef debug
        G4cout<<"G4QCaptureAtRest::AtRestDoIt: (A+2nu) Decay is successful - 2"<<G4endl;
#endif
        G4QHadron* resnuc = new G4QHadron(targPDG,r4Mom); // Creation Hadron for ResidNucl
#ifdef debug
        G4cout<<"G4QCaptureAtRest::AtRestDoIt: ResNuc 4M="<<mt<<r4Mom<<r4Mom.m()<<G4endl;
#endif
        output->push_back(resnuc);               // Fill the Residual Nucleus to the output
#ifdef debug
        G4cout<<"G4QCaptureAtRest::AtRestDoIt: ResNuc is filled nu="<<n4Mom<<nuPDG<<G4endl;
#endif
        G4QHadron* numu = new G4QHadron(nuPDG,n4Mom); // Create Hadron for LeptonicNeutrino
        output->push_back(numu);                      // Fill Muonic Neutrino to the output
#ifdef debug
        G4cout<<"G4QCaptureAtRest::AtRestDoIt:Nu is filled,anu="<<a4Mom<<a4Mom.m()<<G4endl;
#endif
        G4QHadron* anue = new G4QHadron(-12,a4Mom); // Create Hadron for the AntiE Neutrino
        output->push_back(anue);                    // Fill the AntiENeutrino to the output
#ifdef debug
        G4cout<<"G4QCaptureAtRest::AtRestDoIt: anu is filled == Success of Mu-cap"<<G4endl;
#endif
      }
    }
    else                                            // @@Should be developed for tau-lepton
    {
      G4int deL=11;                                 // Prototype of decay lepton
      G4int deN=-12;                                // Prototype of decay neutrino
      G4double mdl=mEl;                             // Prototype of the decay mass
      if(G4UniformRand()>.55)                       // Use mu decay instead of e-decay
      {
        deL=13;
        deN=-14;
        mdl=mMu;
      }
      G4LorentzVector e4Mom(0.,0.,0.,mdl);          // mass of the electron
      G4LorentzVector n4Mom(0.,0.,0.,0.);           // muon neutrino
      G4LorentzVector a4Mom(0.,0.,0.,0.);           // electron anti-nutrino
      if(!G4QHadron(totLV).DecayIn3(e4Mom,n4Mom,a4Mom))
      {
        G4cerr<<"--Worning--G4QCaptureAtRest::AtRestDoIt:Tau_b=>L+Nu_tau+anti_NuL"<<G4endl;
        return 0;
      }
#ifdef debug
      G4cout<<"G4QCaptureAtRest::AtRestDoIt: Tau Decay is successful"<<G4endl;
#endif
      G4QHadron* lept = new G4QHadron(deL,e4Mom);   // Creation Hadron for the Electron
#ifdef debug
      G4cout<<"G4QCaptureAtRest::AtRestDoIt: lepton 4M="<<e4Mom<<e4Mom.m()<<G4endl;
#endif
      output->push_back(lept);                      // Fill the Electron in the output
#ifdef debug
      G4cout<<"G4QCaptureAtRest::AtRestDoIt: lepton is filled nu="<<n4Mom<<nuPDG<<G4endl;
#endif
      G4QHadron* nuta = new G4QHadron(nuPDG,n4Mom); // Create Hadron for LeptonicNeutrino
#ifdef debug
      G4cout<<"G4QCaptureAtRest::AtRestDoIt: nu 4M="<<n4Mom<<n4Mom.m()<<G4endl;
#endif
      output->push_back(nuta);                      // Fill Muonic Neutrino to the output
      G4QHadron* anul = new G4QHadron(deN,a4Mom);   // Create Hadron for the AntiE Neutrino
#ifdef debug
      G4cout<<"G4QCaptureAtRest::AtRestDoIt: antiNu 4M="<<a4Mom<<a4Mom.m()<<G4endl;
#endif
      output->push_back(anul);                      // Fill the AntiENeutrino to the output
    }
  }
  else if((projPDG==13||projPDG==15)&&lepChan&&targPDG==90001000)// LeptonCapture on Proton
  {
    G4LorentzVector totLV(0.,0.,0.,mp+mProt-EnergyDeposition);// 4-mom of theCompoundSystem
#ifdef debug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt:CapOnProton decay 4M="<<totLV<<totLV.m()<<G4endl;
#endif
    G4LorentzVector g4Mom(0.,0.,0.,0.);      // mass of the muon neutrino
    G4LorentzVector n4Mom(0.,0.,0.,mNeut);   // mass of the secondary neutron
    if(!G4QHadron(totLV).DecayIn2(g4Mom,n4Mom))
    {
      G4cerr<<"---Worning---G4QCaptureAtRest::AtRestDoIt: H1+(mu-)=>n+nu_mu"<<G4endl;
      return 0;
    }
    G4QHadron* neutrino = new G4QHadron(nuPDG,g4Mom); // Creation Hadron for neutrino
    output->push_back(neutrino);             // Fill pi0 or Gamma in the output
    G4QHadron* neutron = new G4QHadron(2112,n4Mom);    // Create Hadron for the Neutron
    output->push_back(neutron);              // Fill the neutron to the output
  }
  else if((projPDG==13||projPDG==15)&&lepChan&&targPDG==90001001)//LeptonCapture onDeuteron
  {
    G4LorentzVector totLV(0.,0.,0.,mp+mDeut-EnergyDeposition);// 4-mom of theCompoundSystem
#ifdef debug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt: CapOnDeutr decay 4M="<<totLV<<totLV.m()<<G4endl;
#endif
    G4LorentzVector g4Mom(0.,0.,0.,0.);      // mass of the muon neutrino
    G4LorentzVector n4Mom(0.,0.,0.,mNeut);   // mass of the first neutron
    G4LorentzVector s4Mom(0.,0.,0.,mNeut);   // mass of the second neutron
    if(!G4QHadron(totLV).DecayIn3(g4Mom,n4Mom,s4Mom))
    {
      G4cerr<<"---Worning---G4QCaptureAtRest::AtRestDoIt: D+(mu-)=>n+n+nu_mu"<<G4endl;
      return 0;
    }
    G4QHadron* neutrino = new G4QHadron(nuPDG,g4Mom); // Creation Hadron for the Neutrino
    output->push_back(neutrino);             // Fill pi0 or Gamma in the output
    G4QHadron* neut1 = new G4QHadron(2112,n4Mom);    // Create Hadron for the FirstNeutron
    output->push_back(neut1);                // Fill the neutron to the output
    G4QHadron* neut2 = new G4QHadron(2112,s4Mom);    // Create Hadron for the SecondNeutron
    output->push_back(neut2);                // Fill the neutron to the output
  }
  //                                           *** Closed ***
  else if((projPDG==13||projPDG==15)&&lepChan&&G4UniformRand()>1&&Z>0&&N>0)//@@QuasiFreeCap
  {
    G4double mt=G4QPDGCode(targPDG).GetMass();// Mass of the target Nucleus
    G4LorentzVector totLV(0.,0.,0.,mp+mt-EnergyDeposition);// 4-mom of the(A+mu-) compound
#ifdef debug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt: Quasi-Free decay 4M="<<totLV<<totLV.m()<<G4endl;
#endif
    G4int rPDG=targPDG-1000;                  // Subtract one proton from the nucleus
    G4double mr=G4QPDGCode(rPDG).GetMass();   // Mass of the residual Nucleus
    G4LorentzVector f4Mom(0.,0.,0.,0.);       // Muon neutrino 
    G4LorentzVector s4Mom(0.,0.,0.,mNeut);    // Second neutron
    G4LorentzVector r4Mom(0.,0.,0.,mr);       // Residual nucleus
    if(!G4QHadron(totLV).DecayIn3(f4Mom,s4Mom,r4Mom))
    {
      G4cerr<<"---Worning---G4QCaptureAtRest::AtRestDoIt: A+(mu-)=>nu_mu+n+(A-p)"<<G4endl;
      return 0;
    }
    G4QHadron* neutrino = new G4QHadron(nuPDG,f4Mom); // Create Hadron for the 1st Neutron
    output->push_back(neutrino);              // Fill nutrino_mu in the output
    G4QHadron* neutron = new G4QHadron(2112,s4Mom);// Create Hadron for the 2nd Neutron
    output->push_back(neutron);               // Fill the neutron to the output
    G4QHadron* resnuc = new G4QHadron(rPDG,r4Mom); // Create Hadron for the ResidualNucl
    output->push_back(resnuc);                // Fill the Residual Nucleus to the output
  }
  else                                        // *** THIS works for all particles ! ***
  //else if(1>2)// !! Immediately change back - Very temporary to avoid nuclear capture !!
  {
    G4QHadron* neutr = 0; // Create Neutrino
    if(projPDG==13||projPDG==15) mp-=EnergyDeposition;//TheEnergyDeposit is only for LepCap
#ifdef debug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt: CHIPS M="<<mp<<",dE="<<EnergyDeposition<<G4endl;
#endif
    G4LorentzVector projLV(0.,0.,0.,mp);
    if(projPDG==13) // now for tau it is only energy deposition @@ add similar qqnu decay
    {
      if(G4UniformRand()<.04) // .04 is a parameter !
      {
        projPDG=-211;
        // Phase space decay of mu->nu+q+aq with matrix element
        G4double mmu2=projLV.m2();
        G4double mmu=std::sqrt(mmu2);
        G4double hmm=mmu/2.;
        G4double dmm=mmu+mmu;
        G4double qp=std::pow((std::pow(1.+ga*std::pow(hmm,b1),G4UniformRand())-1.)/ga,rb);
        G4double mqq=0.;
        if(qp<hmm) mqq=std::sqrt(mmu2-dmm*qp);
        G4LorentzVector f4Mom(0.,0.,0.,0.);      // Muon neutrino 
        G4LorentzVector s4Mom(0.,0.,0.,mqq);     // Quark-Antiquark
        if(!G4QHadron(projLV).DecayIn2(f4Mom,s4Mom))
        {
          G4cerr<<"---Worning---G4QCaptureAtRest::AtRestDoIt: (mu-)=>q+aq+nu, M#1"<<G4endl;
          return 0;
        }
        neutr = new G4QHadron(nuPDG,f4Mom);      // Create Neutrino
        projLV-=f4Mom;
        // end of --- ??? ---
      }
    }
    G4QPDGCode targQPDG(targPDG);
    G4double tM=mp+targQPDG.GetMass();
    EnMomConservation=G4LorentzVector(0.,0.,0.,tM);         // Total 4-mom of the reaction

#ifdef tdebug
    G4cout<<"====>G4QCapAR:E/MCons, p="<<mp<<","<<projPDG<<",t="<<tM<<","<<targPDG<<",t4M="
          <<EnMomConservation<<G4endl;
#endif
    G4QHadron* pH = new G4QHadron(projPDG,projLV);          // ---> DELETED---->---->----+
    G4QHadronVector projHV;                                 //                           |
    projHV.push_back(pH);                                   // DESTROYED over 2 lines -+ |
    G4QEnvironment* pan= new G4QEnvironment(projHV,targPDG);// ---> DELETED --->-----+ | |
    std::for_each(projHV.begin(), projHV.end(), DeleteQHadron()); // <---<------<----+-+-+
    projHV.clear(); // <------------<---------------<-------------------<------------+-+
#ifdef debug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt: pPDG="<<projPDG<<", m="<<mp<<G4endl; //   |
#endif
    try                                                           //                 |
    {                                                             //                 |
      delete output;                                              //                 |
      output = pan->Fragment();// DESTROYED in the end of the LOOP work space        |
    }                                                             //                 |
    catch (G4QException& error)//                                                    |
    {                                                             //                 |
      //#ifdef pdebug
      G4cerr<<"***G4QCaptureAtRest::AtRestDoIt: Exception is catched"<<G4endl; //    |
      //#endif
      G4Exception("G4QCaptureAtRest::AtRestDoIt:","27",FatalException,"Gen.CHIPS Except.");
    }                                                             //                 |
    delete pan;                              // Delete the Nuclear Environment <--<--+
#ifdef debug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt: CHIPS fragmentation is done, CV="<<CV<<G4endl;
#endif
    if(neutr) output->push_back(neutr);    // Fill nutrino_mu in the output
  }
  aParticleChange.Initialize(track);
  G4int tNH = output->size();              // A#of hadrons in the output without EM Cascade
  aParticleChange.SetNumberOfSecondaries(tNH); 
  // Now add nuclear fragments
#ifdef debug
  G4cout<<"G4QCaptureAtRest::AtRestDoIt: "<<tNH<<" particles are generated"<<G4endl;
#endif
  // Deal with ParticleChange final state interface to GEANT4 output of the process

  for(i=0; i<tNH; i++)
  {
    // Note that one still has to take care of Hypernuclei (with Lambda or Sigma inside)
    // Hypernucleus mass calculation and ion-table interface upgrade => work for Hisaya @@
    // The decau process for hypernuclei must be developed in GEANT4 (change CHIPS body)
    G4QHadron* hadr=output->operator[](i);   // Pointer to the output hadron    
    if(hadr->GetNFragments())                // Intermediate hadron
    {
#ifdef debug
      G4cout<<"G4QCaptureAtRest::AtRestDoIt: Intermediate particle is found i="<<i<<G4endl;
#endif
      delete hadr;
      continue;
    }
    //VI count Z and A
    countA -= hadr->GetBaryonNumber();
    //VI

    G4DynamicParticle* theSec = new G4DynamicParticle;  
    G4int PDGCode = hadr->GetPDGCode();
#ifdef pdebug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt:#"<<i<<",PDG="<<PDGCode<<G4endl;
#endif
    G4ParticleDefinition* theDefinition=0;
    if     (PDGCode==90000001) theDefinition = G4Neutron::Neutron();
    else if(PDGCode==90001000) theDefinition = G4Proton::Proton();//While it can be in ions
    else if(PDGCode==91000000) theDefinition = G4Lambda::Lambda();
    else if(PDGCode==311 || PDGCode==-311)
    {
      if(G4UniformRand()>.5) theDefinition = G4KaonZeroLong::KaonZeroLong();   // K_L
      else                   theDefinition = G4KaonZeroShort::KaonZeroShort(); // K_S
    }
    else if(PDGCode==91000999) theDefinition = G4SigmaPlus::SigmaPlus();
    else if(PDGCode==90999001) theDefinition = G4SigmaMinus::SigmaMinus();
    else if(PDGCode==91999000) theDefinition = G4XiMinus::XiMinus();
    else if(PDGCode==91999999) theDefinition = G4XiZero::XiZero();
    else if(PDGCode==92998999) theDefinition = G4OmegaMinus::OmegaMinus();
    else if(PDGCode >80000000) // Defines hypernuclei as normal nuclei (N=N+S Correction!)
    {
      G4int aZ = hadr->GetCharge();
      G4int aA = hadr->GetBaryonNumber();
#ifdef pdebug
      G4cout<<"G4QCaptureAtRest::AtRestDoIt:Ion Z="<<aZ<<", A="<<aA<<G4endl;
#endif
      //if      (PDGCode==90001001) theDefinition = G4Deuteron::Deuteron();
      //else if (PDGCode==90001002) theDefinition = G4Triton::Triton();
      //else if (PDGCode==90002001) theDefinition = G4He3::He3();
      //else if (PDGCode==90002002) theDefinition = G4Alpha::Alpha();
      //else
      theDefinition = G4ParticleTable::GetParticleTable()->FindIon(aZ,aA,0,aZ);
    }
    else
    {
#ifdef pdebug
      G4cout<<"G4QCaptureAtRest::AtRestDoIt:Define particle with PDG="<<PDGCode<<G4endl;
#endif
      theDefinition = G4QPDGToG4Particle::Get()->GetParticleDefinition(PDGCode);
#ifdef pdebug
      G4cout<<"G4QCaptureAtRest::AtRestDoIt:AfterParticleDefinition PDG="<<PDGCode<<G4endl;
#endif
    }
    if(!theDefinition)
    {
      G4cout<<"---Worning---G4QCaptureAtRest::AtRestDoIt: drop PDG="<<PDGCode<<G4endl;
      delete hadr;
      continue;
    }
#ifdef pdebug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt:Name="<<theDefinition->GetParticleName()<<G4endl;
#endif
    theSec->SetDefinition(theDefinition);
    G4LorentzVector h4M=hadr->Get4Momentum();
    EnMomConservation-=h4M;

    //VI===== check on final energy
    secEnergy += h4M.e();
    if(h4M.e() > 200*GeV) {
      G4cout<<"G4QCaptureAtRest: wrong product"<<i<<" imax= "<<tNH
	    <<" " << theDefinition->GetParticleName() << "  4-mom= " 
	    << h4M <<G4endl;
    }
    //VI===== 

#ifdef tdebug
    G4cout<<"G4QCapAR::ARDoIt:"<<i<<","<<PDGCode<<h4M<<h4M.m()<<EnMomConservation<<G4endl;
#endif
#ifdef debug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt:#"<<i<<",PDG="<<PDGCode<<",4M="<<h4M<<G4endl;
#endif
    theSec->Set4Momentum(h4M);
    delete hadr;
#ifdef debug
    G4ThreeVector curD=theSec->GetMomentumDirection();
    G4double curM=theSec->GetMass();
    G4double curE=theSec->GetKineticEnergy()+curM;
    G4cout<<"G4QCapAtRest::AtRDoIt:p="<<curD<<curD.mag()<<",e="<<curE<<",m="<<curM<<G4endl;
#endif
    G4Track* aNewTrack = new G4Track(theSec, localtime, position );
    aNewTrack->SetWeight(weight);                                   //    weighted
    aNewTrack->SetTouchableHandle(trTouchable);
    aParticleChange.AddSecondary( aNewTrack );
#ifdef debug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt:#"<<i<<" is done."<<G4endl;
#endif
  }
  delete output;
  //VI=====
  productOK = true;
  countA += Z + N;
  secEnergy += neutron_mass_c2*countA;
  if(EnergyDeposition > limEnergy || std::fabs(secEnergy - primEnergy) > limEnergy) {
    productOK = false;
    G4cout<<"G4QCaptureAtRest::AtRestDoIt: Big energy non-concervation we need redo sampling"<<G4endl;
    G4cout << " Z= " << Z << "  N= " << N << G4endl;
    G4cout<<"G4QCaptureAtRest::AtRestDoIt: the EnergyDeposition(GeV)="<<EnergyDeposition/GeV<<G4endl;
    G4cout<<"  primEnergy(GeV)= " <<primEnergy/GeV << " secEnergy(GeV)="<< secEnergy/GeV <<G4endl;
    for(i=0; i<tNH; i++)
      {
	delete aParticleChange.GetSecondary(i);
      }
    aParticleChange.Clear();
    if(counter >= 100) {
      G4cout<<"G4QCaptureAtRest::AtRestDoIt: Cannot sample final state after " 
	    << counter << "  iterations" <<G4endl;
      G4Exception("G4QCaptureAtRest::AtRestDoIt:","VI",FatalException,"Cannot sample final state");
    }
  }
  } while(!productOK);
#ifdef debug
  G4cout<<"G4QCaptureAtRest::AtRestDoIt: the EnergyDeposition="<<EnergyDeposition<<G4endl;
#endif
  aParticleChange.ProposeLocalEnergyDeposit(EnergyDeposition);// Fill EnergyDeposition
  aParticleChange.ProposeTrackStatus(fStopAndKill);           // Kill the absorbed particle
  //return &aParticleChange;                               // This is not enough (ClearILL)
  return G4VRestProcess::AtRestDoIt(track, step);
}

// The MeanLifeTime (before NucCapture) exists only for MuonCapture, which is a WeakProcess
G4double G4QCaptureAtRest::GetMeanLifeTime(const G4Track& aTrack, G4ForceCondition*)
{
  const G4DynamicParticle* stoppedHadron = aTrack.GetDynamicParticle();
#ifdef debug
  G4cout<<"G4QCaptureAtRest::GetMeanLifeTime is called"<<G4endl;
#endif
  if (*(stoppedHadron->GetDefinition())==*(G4MuonMinus::MuonMinus()) ||
      *(stoppedHadron->GetDefinition())==*(G4TauMinus::TauMinus())      ) return Time;
  else return 0.;
}

// Muon can decay or to be captured by the nucleus (Z,N): true=MuCapture, false=MuDecay
G4bool G4QCaptureAtRest::RandomizeMuDecayOrCapture(G4int Z, G4int N)
{
  static G4int mZ=0;                          // static memory about the last Z
  static G4int mN=0;                          // static memory about the last N
  static G4double mH=0.;                      // static memory about the last Huff(Z,N)
  static G4double mR=0.;                      // static memory about the last rate(Z,N)
  static const G4int nAZ=17;                  // total number of tabulated isotopes
  static const G4int nZm=10;                  // (maximumZ)+1 for which rates are tabulated
  static const G4int rin[nZm]={0,0,1,1,1,2,3,4,5,6}; // i=rin[Z]+N for the tabulated rate
  static const G4double rate[nAZ]={.00000045, .00000047, .00000215, .000000356, .00000468,
                                   .00000226, .00000610, .00002750, .000023500, .00003790,
                                   .00003500, .00006600, .00006200, .000102500, .00009500,
                                   .00008800, .00022900};
#ifdef debug
  G4cout<<"G4QCaptureAtRest::RandomizeMuDecayOrCapture is called"<<G4endl;
#endif
  G4double z=Z;
  G4double Huff=1.;
  if(Z==mZ && N==mN)    Huff=mH;   // Use already calculated value
  else if(Z==8 || Z==9) Huff=.998; // Use nontrivial values for O and F
  // @@ Cash this values ! (M.K.)
  else if(Z>9) Huff=1.-.000394*std::pow(z,2.19)/(1.+12.18*std::exp(z*.01373));
  G4double pD=.00045516*Huff;    // 1/MeanLifeTime of muon in atoms (in ns^-1)?
  G4double pC=1.e99;             // Default 1/MeanLifeTime of muon NuclCapture(in ns^-1)
  // @@ Use Primakov correction (1-1.5625*N/(Z+N)) for isotopes. (M.K.)
  if(Z==mZ && N==mN)    pC=mR;   // Use already calculated value
  else if(Z>9) pC=.0000001256*std::pow(z,3.5)/(1.+.00429*std::exp(1.67*std::pow(z,.39)));
  else if(Z>0) pC=rate[rin[Z]+N];                // Tabulated light isotopes
  else G4cout<<"--Warning--G4QCaptureAtRest::RandomizeMuDecayOrCapture: NEG Z="<<Z<<G4endl;
  mZ=Z; mN=N; mH=Huff; mR=pC;    // Remember the last calculated parameters
  //G4double DLifeT=-std::log(G4UniformRand())/pD;// Time of the muon decay inside the atom
  //G4double CLifeT=-std::log(G4UniformRand())/pC;// Time of the muon capture by nucleus
  //if(DLifeT<CLifeT)
  //{
  //  Time=DLifeT;
#ifdef debug
  //  G4cout<<"G4QCaptureAtRest::RandomizeMuDecayOrCapture: DecayLifeTime="<<Time<<G4endl;
#endif
  //  return false;
  //}
  //else
  //{
  //  Time=CLifeT;
#ifdef debug
  //  G4cout<<"G4QCaptureAtRest::RandomizeMuDecayOrCapture:CaptureLifeTime="<<Time<<G4endl;
#endif
  //  return true;
  //}
  if((pD+pC)*G4UniformRand()>pD) // CAPTURE @@ Cash pD+pC
  {
     Time=-std::log(G4UniformRand())/pC;
     return true;
  }
  else
  {
     Time=-std::log(G4UniformRand())/pD;
     return false;
  }
}

// Calculate the TotalEnergyDeposition for the AtomicCascadeDecay of MuMesoAtom to K-shell
void G4QCaptureAtRest::CalculateEnergyDepositionOfMuCapture(G4int Z)
{
  EnergyDeposition = Z*Z/(3556.+403.4/Z/(2.15+.0039*Z)); // MeV
#ifdef debug
  G4cout<<"G4QCaptureAtR::CalculateEnergyDepositionOfMuCapture="<<EnergyDeposition<<G4endl;
#endif
}

// Muon can decay or to be captured by the nucleus (Z,N): true=TauCapture, false=TauDecay
G4bool G4QCaptureAtRest::RandomizeTauDecayOrCapture(G4int Z, G4int N)
{
#ifdef debug
  G4cout<<"G4QCaptureAtRest::RandomizeMuDecayOrCapture is called"<<G4endl;
#endif
  G4double Z27 =0.002727*Z;
  G4double Z227=Z27*Z27;
  G4double Z427=Z227*Z227;
  G4double Zeff=(Z-0.13782)*(1.2162-(0.09118-Z427)*std::sqrt((G4double)Z));// EffNuclCharge
  G4double Ze2=Zeff*Zeff;      // Squared effective charge of the Nucleus
  G4double pD=3436.*(1.-Ze2*.00014658);     //@@ 1./MeanLifeTime of Tau in atoms (in ns^-1)
  G4double pC=227.*Ze2*Ze2/(33.563+N);      //@@1./MeanLifeTime of TauNuclCapture(in ns^-1)
  if(Z==1&&N==0) pC=10.;                    // @@
  if(Z==1&&N==1) pC=.2;                     // @@
  G4double DLifeT=-std::log(G4UniformRand())/pD; // Time of the muon decay inside the atom
  G4double CLifeT=-std::log(G4UniformRand())/pC; // Time of the muon capture by nucleus
  if(DLifeT<CLifeT)
  {
    Time=DLifeT;
#ifdef debug
    G4cout<<"G4QCaptureAtRest::RandomizeTauDecayOrCapture: DecayLifeTime="<<Time<<G4endl;
#endif
    return false;
  }
  else
  {
    Time=CLifeT;
#ifdef debug
    G4cout<<"G4QCaptureAtRest::RandomizeTauDecayOrCapture: CaptureLifeTime="<<Time<<G4endl;
#endif
    return true;
  }
}

// Calculate the TotalEnergyDeposition for the AtomicCascadeDecay of TauMesoAtom to K-shell
void G4QCaptureAtRest::CalculateEnergyDepositionOfTauCapture(G4int Z)
{
  EnergyDeposition = Z*Z/(21.144+40.38/Z/(2.15+.0039*Z)); // MeV
#ifdef debug
  G4cout<<"G4QCapAtRest::CalculateEnergyDepositionOfTauCapture="<<EnergyDeposition<<G4endl;
#endif
}

// Calculate the TotalEnergyDeposition for the AtomicCascadeDecay of TauMesoAtom to K-shell
G4double G4QCaptureAtRest::RandomizeDecayElectron(G4int z) // E_kin in MeV
{
  static const G4int nZ=12;        // A#of tabulated nuclei
  static const G4int nE=200;       // A#of tabulated energies for each nucleus
  static const G4int nEl=nE-1;     // The last tabulated energy for nuclei
  static const G4int nEb=nE-2;     // The before last tabulated energy for nuclei
  static G4double tZ[nZ]={1.,2.,4.,6.,8.,11.,15.,20.,30.,45.,65.,92.};
  //H1(Z=1)
  static G4double P0[nE]={
          0.00000,7.28739,9.21005,10.5820,11.6892,12.6350,13.4701,14.2238,14.9144,15.5546,
          16.1533,16.7171,17.2511,17.7593,18.2450,18.7106,19.1584,19.5901,20.0073,20.4112,
          20.8031,21.1839,21.5544,21.9153,22.2675,22.6114,22.9476,23.2766,23.5988,23.9146,
          24.2243,24.5284,24.8270,25.1205,25.4091,25.6931,25.9726,26.2479,26.5192,26.7865,
          27.0502,27.3103,27.5670,27.8204,28.0706,28.3178,28.5621,28.8035,29.0422,29.2783,
          29.5118,29.7428,29.9715,30.1978,30.4219,30.6439,30.8637,31.0815,31.2972,31.5111,
          31.7231,31.9332,32.1416,32.3482,32.5531,32.7564,32.9581,33.1583,33.3569,33.5540,
          33.7497,33.9439,34.1368,34.3283,34.5185,34.7074,34.8951,35.0815,35.2667,35.4507,
          35.6335,35.8152,35.9958,36.1753,36.3537,36.5311,36.7075,36.8828,37.0572,37.2306,
          37.4030,37.5745,37.7451,37.9147,38.0835,38.2515,38.4185,38.5847,38.7501,38.9147,
          39.0785,39.2415,39.4038,39.5652,39.7260,39.8859,40.0452,40.2038,40.3616,40.5188,
          40.6753,40.8311,40.9862,41.1407,41.2946,41.4478,41.6004,41.7524,41.9038,42.0546,
          42.2048,42.3545,42.5035,42.6520,42.8000,42.9474,43.0942,43.2405,43.3863,43.5316,
          43.6764,43.8206,43.9644,44.1076,44.2504,44.3927,44.5345,44.6759,44.8168,44.9572,
          45.0972,45.2367,45.3758,45.5145,45.6527,45.7905,45.9279,46.0648,46.2014,46.3375,
          46.4733,46.6086,46.7436,46.8781,47.0123,47.1461,47.2796,47.4126,47.5453,47.6776,
          47.8096,47.9412,48.0724,48.2034,48.3339,48.4641,48.5940,48.7236,48.8528,48.9817,
          49.1103,49.2385,49.3665,49.4941,49.6214,49.7484,49.8751,50.0015,50.1276,50.2534,
          50.3789,50.5041,50.6290,50.7536,50.8780,51.0021,51.1259,51.2494,51.3727,51.4957,
          51.6184,51.7409,51.8633,51.9856,52.1080,52.2311,52.3561,52.4857,52.6273,52.8065};
  //He4(Z=2)
  static G4double P1[nE]={
          0.00000,7.27695,9.19653,10.5662,11.6715,12.6156,13.4493,14.2016,14.8910,15.5300,
          16.1275,16.6903,17.2233,17.7306,18.2153,18.6800,19.1269,19.5578,19.9742,20.3774,
          20.7685,21.1485,21.5182,21.8785,22.2299,22.5732,22.9087,23.2370,23.5585,23.8737,
          24.1828,24.4862,24.7843,25.0772,25.3652,25.6486,25.9275,26.2023,26.4729,26.7398,
          27.0029,27.2625,27.5186,27.7715,28.0212,28.2679,28.5116,28.7525,28.9907,29.2263,
          29.4593,29.6899,29.9180,30.1439,30.3675,30.5890,30.8083,31.0256,31.2410,31.4544,
          31.6659,31.8756,32.0835,32.2897,32.4942,32.6970,32.8983,33.0980,33.2962,33.4929,
          33.6881,33.8820,34.0744,34.2655,34.4553,34.6438,34.8310,35.0170,35.2018,35.3854,
          35.5678,35.7492,35.9293,36.1085,36.2865,36.4635,36.6395,36.8144,36.9884,37.1614,
          37.3334,37.5046,37.6748,37.8441,38.0125,38.1800,38.3467,38.5126,38.6776,38.8418,
          39.0053,39.1679,39.3298,39.4909,39.6512,39.8109,39.9698,40.1280,40.2855,40.4423,
          40.5984,40.7539,40.9087,41.0629,41.2164,41.3693,41.5215,41.6732,41.8242,41.9747,
          42.1246,42.2739,42.4226,42.5707,42.7184,42.8654,43.0119,43.1579,43.3034,43.4483,
          43.5928,43.7367,43.8801,44.0231,44.1655,44.3075,44.4490,44.5900,44.7306,44.8707,
          45.0104,45.1496,45.2884,45.4267,45.5646,45.7021,45.8392,45.9758,46.1121,46.2479,
          46.3833,46.5184,46.6530,46.7873,46.9211,47.0546,47.1878,47.3205,47.4529,47.5849,
          47.7166,47.8479,47.9788,48.1094,48.2397,48.3696,48.4992,48.6285,48.7574,48.8860,
          49.0143,49.1422,49.2699,49.3972,49.5242,49.6509,49.7773,49.9034,50.0292,50.1548,
          50.2800,50.4050,50.5297,50.6541,50.7783,50.9023,51.0261,51.1498,51.2735,51.3973,
          51.5215,51.6462,51.7721,51.8999,52.0308,52.1668,52.3115,52.4720,52.6642,52.9396};
  //Be9(Z=4)
  static G4double P2[nE]={
          0.00000,7.23479,9.14334,10.5051,11.6039,12.5425,13.3712,14.1191,14.8044,15.4396,
          16.0336,16.5930,17.1228,17.6270,18.1088,18.5707,19.0150,19.4432,19.8571,20.2579,
          20.6466,21.0243,21.3918,21.7499,22.0992,22.4404,22.7739,23.1002,23.4198,23.7330,
          24.0403,24.3419,24.6381,24.9292,25.2155,25.4971,25.7744,26.0474,26.3165,26.5817,
          26.8432,27.1012,27.3558,27.6071,27.8553,28.1005,28.3428,28.5823,28.8190,29.0531,
          29.2847,29.5139,29.7407,29.9652,30.1874,30.4076,30.6256,30.8416,31.0556,31.2677,
          31.4779,31.6863,31.8930,32.0979,32.3012,32.5028,32.7029,32.9014,33.0984,33.2939,
          33.4880,33.6806,33.8719,34.0619,34.2505,34.4379,34.6240,34.8088,34.9925,35.1750,
          35.3563,35.5366,35.7157,35.8937,36.0707,36.2466,36.4215,36.5954,36.7683,36.9403,
          37.1113,37.2814,37.4506,37.6189,37.7863,37.9528,38.1185,38.2834,38.4475,38.6107,
          38.7731,38.9348,39.0957,39.2559,39.4153,39.5740,39.7319,39.8892,40.0457,40.2016,
          40.3568,40.5114,40.6652,40.8185,40.9711,41.1231,41.2744,41.4252,41.5754,41.7249,
          41.8739,42.0223,42.1702,42.3174,42.4642,42.6104,42.7560,42.9012,43.0458,43.1899,
          43.3334,43.4765,43.6191,43.7612,43.9028,44.0440,44.1846,44.3248,44.4646,44.6039,
          44.7427,44.8811,45.0191,45.1566,45.2937,45.4304,45.5666,45.7025,45.8379,45.9730,
          46.1076,46.2419,46.3757,46.5092,46.6423,46.7750,46.9074,47.0394,47.1710,47.3023,
          47.4332,47.5637,47.6939,47.8238,47.9533,48.0826,48.2114,48.3400,48.4683,48.5962,
          48.7239,48.8512,48.9783,49.1052,49.2318,49.3581,49.4843,49.6103,49.7362,49.8620,
          49.9878,50.1137,50.2397,50.3660,50.4927,50.6200,50.7481,50.8774,51.0083,51.1412,
          51.2769,51.4163,51.5607,51.7120,51.8729,52.0476,52.2435,52.4742,52.7705,53.2313};
  //C12(Z=6)
  static G4double P3[nE]={
          0.00000,7.18201,9.07814,10.4312,11.5230,12.4557,13.2793,14.0225,14.7036,15.3349,
          15.9253,16.4814,17.0080,17.5092,17.9882,18.4474,18.8890,19.3148,19.7263,20.1247,
          20.5112,20.8867,21.2521,21.6082,21.9555,22.2948,22.6264,22.9509,23.2687,23.5802,
          23.8858,24.1857,24.4803,24.7699,25.0546,25.3347,25.6105,25.8821,26.1497,26.4135,
          26.6737,26.9303,27.1836,27.4336,27.6805,27.9244,28.1654,28.4037,28.6392,28.8721,
          29.1026,29.3306,29.5562,29.7796,30.0007,30.2198,30.4367,30.6516,30.8646,31.0756,
          31.2848,31.4923,31.6979,31.9019,32.1041,32.3048,32.5039,32.7015,32.8975,33.0921,
          33.2853,33.4770,33.6674,33.8565,34.0443,34.2308,34.4160,34.6000,34.7829,34.9645,
          35.1451,35.3245,35.5028,35.6800,35.8562,36.0313,36.2055,36.3786,36.5508,36.7220,
          36.8923,37.0616,37.2301,37.3976,37.5643,37.7302,37.8952,38.0593,38.2227,38.3852,
          38.5470,38.7080,38.8683,39.0277,39.1865,39.3445,39.5018,39.6585,39.8144,39.9696,
          40.1242,40.2781,40.4314,40.5840,40.7360,40.8874,41.0382,41.1884,41.3379,41.4869,
          41.6353,41.7832,41.9305,42.0772,42.2234,42.3690,42.5141,42.6587,42.8028,42.9463,
          43.0894,43.2319,43.3740,43.5156,43.6567,43.7973,43.9375,44.0772,44.2164,44.3552,
          44.4936,44.6315,44.7690,44.9060,45.0426,45.1789,45.3147,45.4501,45.5851,45.7197,
          45.8539,45.9877,46.1211,46.2542,46.3869,46.5193,46.6512,46.7829,46.9142,47.0451,
          47.1758,47.3061,47.4361,47.5658,47.6952,47.8243,47.9532,48.0819,48.2103,48.3385,
          48.4666,48.5945,48.7223,48.8501,48.9778,49.1056,49.2335,49.3616,49.4899,49.6187,
          49.7479,49.8779,50.0087,50.1406,50.2738,50.4087,50.5457,50.6854,50.8282,50.9752,
          51.1272,51.2858,51.4528,51.6310,51.8244,52.0390,52.2854,52.5832,52.9766,53.6078};
  //O16(Z=8)
  static G4double P4[nE]={
          0.00000,7.11944,9.00194,10.3456,11.4301,12.3566,13.1748,13.9133,14.5901,15.2176,
          15.8044,16.3570,16.8806,17.3789,17.8550,18.3116,18.7507,19.1741,19.5833,19.9796,
          20.3640,20.7375,21.1010,21.4552,21.8007,22.1382,22.4681,22.7910,23.1072,23.4172,
          23.7213,24.0198,24.3130,24.6011,24.8845,25.1633,25.4378,25.7082,25.9746,26.2372,
          26.4962,26.7517,27.0039,27.2528,27.4987,27.7415,27.9816,28.2188,28.4534,28.6854,
          28.9149,29.1420,29.3667,29.5892,29.8095,30.0277,30.2438,30.4579,30.6701,30.8803,
          31.0888,31.2954,31.5004,31.7036,31.9052,32.1051,32.3035,32.5004,32.6958,32.8898,
          33.0823,33.2734,33.4632,33.6517,33.8388,34.0247,34.2094,34.3929,34.5751,34.7563,
          34.9363,35.1151,35.2929,35.4696,35.6453,35.8200,35.9936,36.1663,36.3380,36.5088,
          36.6786,36.8475,37.0155,37.1827,37.3489,37.5144,37.6789,37.8427,38.0057,38.1678,
          38.3292,38.4899,38.6497,38.8089,38.9673,39.1250,39.2819,39.4382,39.5938,39.7488,
          39.9030,40.0566,40.2096,40.3619,40.5136,40.6647,40.8152,40.9651,41.1144,41.2632,
          41.4113,41.5589,41.7059,41.8524,41.9984,42.1438,42.2887,42.4330,42.5769,42.7202,
          42.8631,43.0054,43.1473,43.2887,43.4296,43.5701,43.7101,43.8496,43.9887,44.1274,
          44.2656,44.4034,44.5407,44.6777,44.8142,44.9503,45.0861,45.2214,45.3564,45.4910,
          45.6252,45.7590,45.8925,46.0256,46.1584,46.2909,46.4231,46.5549,46.6865,46.8178,
          46.9488,47.0796,47.2101,47.3404,47.4706,47.6005,47.7304,47.8601,47.9898,48.1194,
          48.2491,48.3788,48.5087,48.6388,48.7692,48.8999,49.0311,49.1629,49.2955,49.4289,
          49.5635,49.6993,49.8367,49.9760,50.1175,50.2617,50.4091,50.5605,50.7166,50.8785,
          51.0475,51.2256,51.4151,51.6195,51.8439,52.0962,52.3897,52.7496,53.2325,54.0205};
  //Na23(Z=11)
  static G4double P5[nE]={
          0.00000,7.00695,8.86633,10.1943,11.2667,12.1832,12.9928,13.7237,14.3938,15.0151,
          15.5963,16.1438,16.6625,17.1563,17.6283,18.0810,18.5164,18.9362,19.3421,19.7352,
          20.1165,20.4872,20.8479,21.1994,21.5424,21.8775,22.2051,22.5257,22.8397,23.1476,
          23.4497,23.7462,24.0375,24.3239,24.6055,24.8826,25.1555,25.4243,25.6892,25.9503,
          26.2078,26.4619,26.7127,26.9604,27.2050,27.4466,27.6854,27.9215,28.1550,28.3859,
          28.6143,28.8404,29.0641,29.2857,29.5050,29.7223,29.9375,30.1508,30.3621,30.5716,
          30.7793,30.9852,31.1894,31.3919,31.5928,31.7921,31.9899,32.1862,32.3810,32.5743,
          32.7663,32.9569,33.1462,33.3341,33.5208,33.7063,33.8905,34.0735,34.2554,34.4361,
          34.6157,34.7942,34.9717,35.1481,35.3235,35.4978,35.6712,35.8436,36.0150,36.1856,
          36.3552,36.5239,36.6917,36.8586,37.0247,37.1900,37.3544,37.5181,37.6809,37.8430,
          38.0043,38.1648,38.3246,38.4837,38.6420,38.7997,38.9566,39.1129,39.2685,39.4234,
          39.5777,39.7314,39.8844,40.0368,40.1885,40.3397,40.4903,40.6403,40.7897,40.9385,
          41.0868,41.2345,41.3817,41.5284,41.6745,41.8201,41.9652,42.1097,42.2538,42.3974,
          42.5405,42.6831,42.8253,42.9670,43.1082,43.2490,43.3893,43.5292,43.6687,43.8078,
          43.9464,44.0847,44.2225,44.3600,44.4971,44.6338,44.7702,44.9062,45.0418,45.1772,
          45.3122,45.4469,45.5814,45.7155,45.8494,45.9831,46.1165,46.2497,46.3828,46.5156,
          46.6484,46.7810,46.9136,47.0462,47.1787,47.3113,47.4439,47.5768,47.7098,47.8430,
          47.9767,48.1107,48.2452,48.3804,48.5163,48.6530,48.7908,48.9297,49.0700,49.2119,
          49.3556,49.5015,49.6499,49.8012,49.9559,50.1146,50.2780,50.4470,50.6226,50.8063,
          50.9997,51.2054,51.4263,51.6671,51.9343,52.2382,52.5959,53.0400,53.6438,54.6445};
  //P31(Z=15)
  static G4double P6[nE]={
          0.00000,6.82422,8.64751,9.95140,11.0052,11.9065,12.7032,13.4228,14.0829,14.6952,
          15.2682,15.8082,16.3200,16.8074,17.2734,17.7205,18.1506,18.5656,18.9668,19.3555,
          19.7327,20.0995,20.4564,20.8044,21.1441,21.4759,21.8005,22.1182,22.4294,22.7347,
          23.0342,23.3283,23.6173,23.9014,24.1809,24.4560,24.7270,24.9939,25.2570,25.5164,
          25.7723,26.0249,26.2742,26.5204,26.7637,27.0040,27.2416,27.4765,27.7088,27.9387,
          28.1661,28.3912,28.6141,28.8347,29.0533,29.2698,29.4843,29.6969,29.9076,30.1165,
          30.3236,30.5290,30.7327,30.9348,31.1353,31.3342,31.5317,31.7277,31.9222,32.1153,
          32.3071,32.4975,32.6867,32.8745,33.0611,33.2465,33.4308,33.6138,33.7957,33.9765,
          34.1562,34.3349,34.5125,34.6890,34.8646,35.0392,35.2128,35.3855,35.5573,35.7281,
          35.8980,36.0671,36.2353,36.4027,36.5692,36.7350,36.8999,37.0640,37.2274,37.3900,
          37.5519,37.7130,37.8734,38.0331,38.1921,38.3505,38.5081,38.6651,38.8214,38.9771,
          39.1322,39.2866,39.4404,39.5937,39.7463,39.8983,40.0498,40.2007,40.3511,40.5009,
          40.6502,40.7989,40.9471,41.0948,41.2420,41.3887,41.5349,41.6806,41.8259,41.9707,
          42.1150,42.2589,42.4024,42.5454,42.6880,42.8302,42.9720,43.1134,43.2545,43.3951,
          43.5354,43.6754,43.8150,43.9542,44.0932,44.2319,44.3702,44.5084,44.6462,44.7838,
          44.9212,45.0584,45.1954,45.3323,45.4690,45.6057,45.7422,45.8787,46.0152,46.1517,
          46.2882,46.4249,46.5617,46.6987,46.8360,46.9735,47.1115,47.2499,47.3888,47.5284,
          47.6687,47.8098,47.9518,48.0950,48.2394,48.3852,48.5326,48.6819,48.8332,48.9869,
          49.1432,49.3027,49.4657,49.6327,49.8044,49.9815,50.1649,50.3557,50.5552,50.7653,
          50.9882,51.2267,51.4850,51.7686,52.0858,52.4496,52.8817,53.4232,54.1672,55.4144};
  //Ca40(Z=20)
  static G4double P7[nE]={
          0.00000,6.55304,8.32352,9.59220,10.6190,11.4982,12.2762,12.9795,13.6251,14.2244,
          14.7856,15.3149,15.8168,16.2951,16.7526,17.1918,17.6146,18.0227,18.4174,18.8000,
          19.1714,19.5327,19.8846,20.2277,20.5627,20.8902,21.2106,21.5244,21.8320,22.1337,
          22.4298,22.7208,23.0067,23.2880,23.5648,23.8373,24.1058,24.3704,24.6313,24.8886,
          25.1426,25.3933,25.6409,25.8855,26.1272,26.3661,26.6023,26.8360,27.0671,27.2959,
          27.5223,27.7465,27.9685,28.1884,28.4062,28.6221,28.8360,29.0481,29.2584,29.4669,
          29.6737,29.8789,30.0824,30.2844,30.4848,30.6838,30.8813,31.0774,31.2721,31.4654,
          31.6575,31.8482,32.0377,32.2260,32.4131,32.5990,32.7838,32.9675,33.1501,33.3316,
          33.5120,33.6915,33.8699,34.0474,34.2238,34.3994,34.5740,34.7478,34.9206,35.0926,
          35.2637,35.4339,35.6034,35.7720,35.9399,36.1070,36.2733,36.4389,36.6037,36.7678,
          36.9312,37.0939,37.2560,37.4173,37.5780,37.7380,37.8974,38.0562,38.2144,38.3719,
          38.5289,38.6853,38.8411,38.9963,39.1510,39.3051,39.4587,39.6117,39.7643,39.9163,
          40.0679,40.2189,40.3695,40.5196,40.6692,40.8184,40.9671,41.1154,41.2633,41.4108,
          41.5579,41.7045,41.8508,41.9967,42.1423,42.2875,42.4324,42.5769,42.7212,42.8651,
          43.0087,43.1521,43.2953,43.4381,43.5808,43.7233,43.8655,44.0077,44.1496,44.2915,
          44.4332,44.5749,44.7166,44.8582,44.9999,45.1416,45.2834,45.4253,45.5675,45.7098,
          45.8524,45.9953,46.1386,46.2824,46.4267,46.5716,46.7172,46.8635,47.0107,47.1589,
          47.3083,47.4588,47.6108,47.7643,47.9195,48.0768,48.2362,48.3981,48.5628,48.7306,
          48.9019,49.0772,49.2570,49.4420,49.6329,49.8306,50.0362,50.2510,50.4768,50.7155,
          50.9699,51.2436,51.5414,51.8702,52.2401,52.6667,53.1765,53.8195,54.7092,56.213};
  //Zn64(Z=30)
  static G4double P8[nE]={
          0.00000,5.93631,7.58588,8.77313,9.73699,10.5643,11.2979,11.9624,12.5734,13.1415,
          13.6742,14.1774,14.6552,15.1111,15.5477,15.9674,16.3719,16.7627,17.1412,17.5084,
          17.8654,18.2129,18.5517,18.8825,19.2058,19.5221,19.8319,20.1356,20.4335,20.7260,
          21.0135,21.2961,21.5742,21.8479,22.1175,22.3833,22.6453,22.9037,23.1588,23.4106,
          23.6593,23.9050,24.1479,24.3880,24.6255,24.8605,25.0931,25.3232,25.5512,25.7769,
          26.0005,26.2221,26.4417,26.6594,26.8753,27.0893,27.3017,27.5123,27.7214,27.9288,
          28.1348,28.3392,28.5422,28.7437,28.9439,29.1428,29.3404,29.5367,29.7318,29.9257,
          30.1184,30.3100,30.5005,30.6899,30.8782,31.0655,31.2518,31.4371,31.6215,31.8049,
          31.9874,32.1690,32.3497,32.5296,32.7086,32.8869,33.0643,33.2409,33.4168,33.5919,
          33.7663,33.9399,34.1129,34.2851,34.4567,34.6277,34.7979,34.9676,35.1366,35.3050,
          35.4728,35.6400,35.8067,35.9727,36.1383,36.3033,36.4677,36.6317,36.7951,36.9580,
          37.1205,37.2825,37.4440,37.6050,37.7656,37.9258,38.0855,38.2449,38.4038,38.5623,
          38.7205,38.8782,39.0356,39.1927,39.3494,39.5058,39.6618,39.8176,39.9730,40.1281,
          40.2830,40.4376,40.5920,40.7461,40.9000,41.0537,41.2072,41.3606,41.5137,41.6667,
          41.8196,41.9724,42.1251,42.2777,42.4303,42.5829,42.7354,42.8880,43.0407,43.1934,
          43.3462,43.4992,43.6524,43.8058,43.9595,44.1135,44.2679,44.4226,44.5779,44.7336,
          44.8900,45.0470,45.2047,45.3632,45.5227,45.6831,45.8447,46.0075,46.1716,46.3372,
          46.5045,46.6735,46.8446,47.0179,47.1936,47.3720,47.5535,47.7382,47.9267,48.1193,
          48.3166,48.5191,48.7275,48.9426,49.1653,49.3968,49.6383,49.8917,50.1588,50.4424,
          50.7458,51.0735,51.4316,51.8286,52.2771,52.7968,53.4209,54.2121,55.3131,57.1875};
  //Rh103(Z=45)
  static G4double P9[nE]={
          0.00000,5.00701,6.48187,7.55018,8.42075,9.17008,9.83596,10.4402,10.9967,11.5149,
          12.0017,12.4619,12.8995,13.3176,13.7185,14.1042,14.4765,14.8365,15.1856,15.5248,
          15.8548,16.1764,16.4903,16.7972,17.0974,17.3914,17.6797,17.9627,18.2406,18.5137,
          18.7825,19.0470,19.3075,19.5643,19.8176,20.0674,20.3141,20.5576,20.7983,21.0362,
          21.2714,21.5041,21.7344,21.9623,22.1880,22.4116,22.6331,22.8527,23.0703,23.2862,
          23.5003,23.7127,23.9234,24.1327,24.3404,24.5466,24.7514,24.9549,25.1571,25.3580,
          25.5577,25.7562,25.9535,26.1497,26.3449,26.5390,26.7322,26.9243,27.1155,27.3058,
          27.4952,27.6837,27.8714,28.0583,28.2445,28.4298,28.6144,28.7984,28.9816,29.1641,
          29.3460,29.5273,29.7079,29.8880,30.0675,30.2464,30.4248,30.6026,30.7799,30.9568,
          31.1332,31.3091,31.4845,31.6595,31.8341,32.0083,32.1821,32.3555,32.5286,32.7013,
          32.8736,33.0456,33.2173,33.3887,33.5598,33.7306,33.9011,34.0713,34.2414,34.4111,
          34.5806,34.7500,34.9191,35.0880,35.2567,35.4252,35.5936,35.7618,35.9299,36.0978,
          36.2656,36.4333,36.6008,36.7683,36.9357,37.1031,37.2704,37.4376,37.6048,37.7720,
          37.9391,38.1063,38.2735,38.4407,38.6080,38.7754,38.9428,39.1103,39.2779,39.4457,
          39.6136,39.7817,39.9500,40.1184,40.2872,40.4562,40.6255,40.7951,40.9650,41.1354,
          41.3061,41.4773,41.6490,41.8213,41.9941,42.1675,42.3417,42.5165,42.6922,42.8687,
          43.0462,43.2246,43.4042,43.5849,43.7669,43.9504,44.1353,44.3218,44.5102,44.7004,
          44.8928,45.0875,45.2847,45.4846,45.6876,45.8939,46.1038,46.3178,46.5362,46.7596,
          46.9885,47.2236,47.4656,47.7155,47.9743,48.2434,48.5242,48.8187,49.1292,49.4588,
          49.8113,50.1919,50.6075,51.0680,51.5878,52.1895,52.9113,53.8254,55.0955,57.2547};
  //Th159(Z=65)
  static G4double PA[nE]={
          0.00000,3.72463,4.98461,5.91463,6.68054,7.34437,7.93719,8.47714,8.97585,9.44132,
          9.87929,10.2941,10.6890,11.0666,11.4290,11.7780,12.1150,12.4412,12.7576,13.0651,
          13.3644,13.6563,13.9413,14.2199,14.4925,14.7597,15.0217,15.2789,15.5315,15.7800,
          16.0244,16.2651,16.5023,16.7361,16.9667,17.1943,17.4190,17.6410,17.8604,18.0774,
          18.2920,18.5044,18.7147,18.9229,19.1292,19.3336,19.5362,19.7371,19.9364,20.1341,
          20.3303,20.5251,20.7185,20.9106,21.1014,21.2910,21.4794,21.6667,21.8529,22.0381,
          22.2223,22.4055,22.5879,22.7693,22.9500,23.1298,23.3088,23.4871,23.6647,23.8416,
          24.0179,24.1935,24.3685,24.5430,24.7170,24.8904,25.0633,25.2358,25.4078,25.5794,
          25.7506,25.9214,26.0919,26.2620,26.4319,26.6014,26.7706,26.9396,27.1084,27.2770,
          27.4453,27.6135,27.7815,27.9494,28.1171,28.2847,28.4522,28.6197,28.7871,28.9545,
          29.1218,29.2891,29.4564,29.6237,29.7911,29.9586,30.1260,30.2936,30.4613,30.6291,
          30.7970,30.9651,31.1333,31.3018,31.4704,31.6392,31.8082,31.9775,32.1471,32.3169,
          32.4870,32.6574,32.8281,32.9992,33.1706,33.3425,33.5147,33.6873,33.8603,34.0338,
          34.2078,34.3822,34.5572,34.7326,34.9087,35.0852,35.2624,35.4402,35.6187,35.7978,
          35.9775,36.1580,36.3393,36.5213,36.7041,36.8877,37.0722,37.2576,37.4439,37.6312,
          37.8196,38.0089,38.1994,38.3910,38.5838,38.7778,38.9731,39.1698,39.3680,39.5676,
          39.7688,39.9716,40.1762,40.3827,40.5910,40.8014,41.0140,41.2289,41.4463,41.6662,
          41.8890,42.1147,42.3437,42.5761,42.8122,43.0524,43.2969,43.5462,43.8008,44.0610,
          44.3276,44.6012,44.8826,45.1727,45.4727,45.7838,46.1077,46.4465,46.8024,47.1788,
          47.5795,48.0100,48.4774,48.9920,49.5687,50.2309,51.0179,52.0041,53.3576,55.6238};

  //U238(Z=92)
  static G4double PB[nE]={
          0.00000,1.93933,2.80877,3.48857,4.06872,4.58456,5.05440,5.48917,5.89607,6.28010,
          6.64494,6.99335,7.32751,7.64913,7.95964,8.26019,8.55175,8.83516,9.11111,9.38023,
          9.64304,9.90002,10.1516,10.3981,10.6399,10.8772,11.1105,11.3398,11.5654,11.7875,
          12.0063,12.2220,12.4347,12.6446,12.8518,13.0564,13.2586,13.4585,13.6561,13.8516,
          14.0450,14.2364,14.4260,14.6137,14.7998,14.9841,15.1668,15.3480,15.5278,15.7060,
          15.8829,16.0585,16.2328,16.4059,16.5778,16.7485,16.9181,17.0867,17.2543,17.4208,
          17.5864,17.7511,17.9149,18.0779,18.2400,18.4013,18.5619,18.7218,18.8809,19.0394,
          19.1972,19.3544,19.5110,19.6670,19.8225,19.9774,20.1318,20.2858,20.4393,20.5924,
          20.7450,20.8973,21.0492,21.2007,21.3519,21.5028,21.6534,21.8038,21.9539,22.1037,
          22.2534,22.4029,22.5522,22.7013,22.8503,22.9992,23.1481,23.2968,23.4455,23.5942,
          23.7428,23.8915,24.0402,24.1889,24.3377,24.4866,24.6355,24.7846,24.9339,25.0833,
          25.2329,25.3828,25.5328,25.6832,25.8337,25.9846,26.1359,26.2874,26.4394,26.5917,
          26.7444,26.8976,27.0513,27.2055,27.3602,27.5155,27.6713,27.8278,27.9849,28.1427,
          28.3012,28.4604,28.6204,28.7812,28.9429,29.1055,29.2690,29.4334,29.5989,29.7655,
          29.9331,30.1019,30.2719,30.4431,30.6156,30.7895,30.9648,31.1416,31.3200,31.4999,
          31.6816,31.8650,32.0502,32.2374,32.4266,32.6178,32.8113,33.0071,33.2053,33.4060,
          33.6094,33.8156,34.0248,34.2370,34.4525,34.6714,34.8939,35.1203,35.3508,35.5855,
          35.8249,36.0691,36.3185,36.5734,36.8343,37.1015,37.3756,37.6571,37.9466,38.2448,
          38.5525,38.8705,39.2000,39.5422,39.8984,40.2703,40.6601,41.0701,41.5034,41.9638,
          42.4562,42.9870,43.5648,44.2017,44.9151,45.7322,46.6981,47.8972,49.5188,52.1676};
  static G4double* PP[nZ]= {P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,PA,PB};
  static const G4int mZ1=93;      // MaxPossibleZ+1
  static const G4int mZ=mZ1-1;    // MaxPossibleZ
  static G4double* P[mZ1]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  //---------------------------------------------------------------------------------------
  if(z<1 || z>mZ)
  {
    G4cout<<"-W-G4QCapAtRest::RandomizeDecayElectron: <=0 or big(>"<<mZ<<") Z="<<z<<G4endl;
    return 0.;
  }
  G4double  Z  = z;
  G4double* nP = 0;
  if(!P[z])          // The table for this element must be created
  {
    for(G4int i=0; i<nZ; i++)
    {
      G4double fZ=tZ[i];
      if(Z<=fZ)      // The extrapolation bin is found
      {
#ifdef debug
        G4cout<<"G4QCapAtRest::RandomizeDecayElectron: Z="<<Z<<", fZ="<<fZ<<G4endl;
#endif
        nP = new G4double[nE];
        G4double* fP=PP[i];
        if(Z==fZ) for(G4int j=0; j<nE; j++) nP[j]=fP[j];
        else
        {
          G4int i1=i-1;  // i>2, asthe first tabilated Z are 1,2,4,... and min_i=3
          G4double  iZ=tZ[i1];
          G4double* iP=PP[i1];
          G4double rZ=(Z-iZ)/(fZ-iZ);
          for(G4int j=0; j<nE; j++) nP[j]=iP[j]+(fP[j]-iP[j])*rZ;
        }
#ifdef debug
        for(G4int k=0; k<nE; k++)G4cout<<"G4QCAR::RandDecEle:P["<<k<<"]="<<nP[k]<<G4endl;
#endif
        P[z]=nP;
        break;
      }
    }
  }
  else nP=P[z];
  // At this point the randomization table for the element is in nP
  G4double R=G4UniformRand()*nE;
  G4int iR=static_cast<int>(R);
  if(iR > nEl) iR=nEl;
  else if(iR<0) iR=0;
  G4double nPi=nP[iR];
  G4double nPf=0.;
  if(iR<nEl) nPf=nP[iR+1];
  else       nPf=nP[nEl]+(nP[nEl]-nP[nEb])/2; // An artificial tail
#ifdef debug
  G4cout<<"G4QCapAtR::RaDEl:R="<<R<<",Ei="<<nPi<<",E="<<MeV*(nPi+(R-iR)*(nPf-nPi))<<G4endl;
#endif
  return MeV*(nPi+(R-iR)*(nPf-nPi));
}
