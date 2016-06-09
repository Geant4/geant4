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
// $Id: G4QCaptureAtRest.cc,v 1.20 2005/06/27 15:30:29 gunter Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//      ---------------- G4QCaptureAtRest class -----------------
//                 by Mikhail Kossov, December 2003.
// G4QCaptureAtRest class of the CHIPS Simulation Branch in GEANT4
// ---------------------------------------------------------------
// ****************************************************************************************
// ********** This CLASS is temporary moved from the photolepton_hadron directory *********
// ******* DO NOT MAKE ANY CHANGE! With time it'll move back to photolepton...(M.K.) ******
// ****************************************************************************************

//#define debug
//#define pdebug
//#define tdebug

#include "G4QCaptureAtRest.hh"

G4QCaptureAtRest::G4QCaptureAtRest(const G4String& processName)
  : G4VRestProcess(processName), Time(0.), EnergyDeposition(0.)
{
#ifdef debug
  G4cout<<"G4QCaptureAtRest::Constructor is called"<<G4endl;
#endif
  if (verboseLevel>0) G4cout << GetProcessName() << " is created "<< G4endl;

  G4QCHIPSWorld::Get()->GetParticles(nPartCWorld); // Create CHIPS World with 234 particles
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
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mPi0 = G4QPDGCode(111).GetMass();
  static const G4double mDeut= G4QPDGCode(2112).GetNuclMass(1,1,0);
  //static const G4double mPi  = G4QPDGCode(211).GetMass();
  //static const G4double mMu  = G4QPDGCode(13).GetMass();
  //static const G4double mTau = G4QPDGCode(15).GetMass();
  static const G4double mEl  = G4QPDGCode(11).GetMass();
  //
  const G4DynamicParticle* stoppedHadron = track.GetDynamicParticle();
  const G4ParticleDefinition* particle=stoppedHadron->GetDefinition();
  Time=0.;
  EnergyDeposition=0.;
#ifdef debug
  G4cout<<"G4QCaptureAtRest::AtRestDoIt is called,EnDeposition="<<EnergyDeposition<<G4endl;
#endif
  if (! IsApplicable(*particle))  // Check applicability
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
    if(projPDG==13||projPDG==15)
    {
      G4int cZ=static_cast<G4int>((*theElementVector)[i]->GetZ());
      frac*=cZ;
      if(cZ==9||cZ==35||cZ==53||cZ==85) frac*=.66;
      else if                  (cZ== 3) frac*=.50;
      else if          (cZ==24||cZ==28) frac*=.90;
      else if          (cZ== 5||cZ==17) frac*=.70;
      else if                  (cZ== 8) frac*=.56;
    }
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
  G4int N = Z;
  G4int isoSize=0;                         // The default for the isoVectorLength is 0
  G4IsotopeVector* isoVector=pElement->GetIsotopeVector();
  if(isoVector) isoSize=isoVector->size(); // Get real size of the isotopeVector if exists
#ifdef debug
  G4cout<<"G4QCaptureAtRest::AtRestDoIt: isovectorLength="<<isoSize<<G4endl;
#endif
  if(isoSize)                         // The Element has not trivial abumdance set
  {
    // @@ the following solution is temporary till G4Element can contain the QIsotopIndex
    G4int curInd=G4QIsotope::Get()->GetLastIndex(Z);
    if(!curInd)                       // The new artificial element must be defined 
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
      curInd=G4QIsotope::Get()->InitElement(Z,1,newAbund);
      for(G4int k=0; k<isoSize; k++) delete (*newAbund)[k];
      delete newAbund;
    }
    // @@ ^^^^^^^^^^ End of the temporary solution ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    N = G4QIsotope::Get()->GetNeutrons(Z,curInd);
  }
  else  N = G4QIsotope::Get()->GetNeutrons(Z);
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
  G4QHadronVector* output=new G4QHadronVector; // Prototype of the output G4QHadronVector
#ifdef debug
  G4cout<<"G4QCaptureAtRest::AtRestDoIt: projPDG="<<projPDG<<", targPDG="<<targPDG<<G4endl;
#endif
  G4double      localtime = track.GetGlobalTime();
  G4ThreeVector position  = track.GetPosition();
  localtime += Time;
	 std::vector<G4double>* cascE = new std::vector<G4double>;
	 std::vector<G4Track*>* cascT = new std::vector<G4Track*>;
  G4int           nuPDG=14;                  // Prototype for weak decay
  if(projPDG==15) nuPDG=16;
  if(projPDG==-211 && targPDG==90001000)     // Use Panofsky ratio for (p+pi-) system decay
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
    G4LorentzVector totLV(0.,0.,0.,mp-EnergyDeposition);// 4-momentum of the bounded muon
#ifdef debug
	   G4cout<<"G4QCaptureAtRest::AtRestDoIt: e+nu+nu decay 4M="<<totLV<<totLV.m()<<G4endl;
#endif
    // @@ Should be developed for tau-lepton
    G4LorentzVector e4Mom(0.,0.,0.,mEl);     // mass of the electron
    G4LorentzVector n4Mom(0.,0.,0.,0.);      // muon neutrino
    G4LorentzVector a4Mom(0.,0.,0.,0.);      // electron anti-nutrino
    if(!G4QHadron(totLV).DecayIn3(e4Mom,n4Mom,a4Mom))
    {
      G4cerr<<"---Worning---G4QCaptureAtRest::AtRestDoIt: Mu_b=>E+Nu_mu+anti_Nu_e"<<G4endl;
      return 0;
    }
#ifdef debug
	   G4cout<<"G4QCaptureAtRest::AtRestDoIt: Decay is successful"<<G4endl;
#endif
    G4QHadron* elect = new G4QHadron(11,e4Mom);   // Creation Hadron for the Electron
#ifdef debug
	   G4cout<<"G4QCaptureAtRest::AtRestDoIt: electron 4M="<<e4Mom<<e4Mom.m()<<G4endl;
#endif
    output->push_back(elect);                     // Fill the Electron in the output
#ifdef debug
	   G4cout<<"G4QCaptureAtRest::AtRestDoIt: electron is filled nu4M="<<n4Mom<<nuPDG<<G4endl;
#endif
    G4QHadron* numu = new G4QHadron(nuPDG,n4Mom); // Create Hadron for the LeptonicNeutrino
#ifdef debug
	   G4cout<<"G4QCaptureAtRest::AtRestDoIt: nu 4M="<<n4Mom<<n4Mom.m()<<G4endl;
#endif
    output->push_back(numu);                      // Fill the Muonic Neutrino to the output
    G4QHadron* anue = new G4QHadron(-12,a4Mom);   // Create Hadron for the AntiE Neutrino
#ifdef debug
	   G4cout<<"G4QCaptureAtRest::AtRestDoIt: antiNu 4M="<<a4Mom<<a4Mom.m()<<G4endl;
#endif
    output->push_back(anue);                      // Fill the AntiE Neutrino to the output
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
  else
  {
    if(projPDG==13||projPDG==15) mp-=EnergyDeposition;//TheEnergyDeposit is only for LepCap
#ifdef debug
	   G4cout<<"G4QCaptureAtRest::AtRestDoIt: CHIPS M="<<mp<<",dE="<<EnergyDeposition<<G4endl;
#endif
    G4LorentzVector projLV(0.,0.,0.,mp);
    // @@ it can be a flag @@ now for tau it is only energy deposition, for mu +EMCascade
    if(projPDG==13)
    {
      MuCaptureEMCascade(Z, N, cascE);
      G4int nsec=cascE->size();
      G4DynamicParticle* theSec = 0; // Prototype to fill particle in the G4ParticleChange
	     for(G4int is=0; is<nsec; is++)
	     {
        G4double ener=cascE->operator[](is);
        if(ener>0) theSec = new G4DynamicParticle(G4Electron::Electron(),RndmDir(),ener);
        else       theSec = new G4DynamicParticle(G4Gamma::Gamma(),RndmDir(),-ener);
        projLV-=theSec->Get4Momentum();
        G4Track* aNewTrack = new G4Track(theSec, localtime, position );
        cascT->push_back(aNewTrack);
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
  }
  aParticleChange.Initialize(track);
  G4int tNH = output->size(); // A#of hadrons in the output without EM Cascade
  G4int nsec=cascE->size();
  aParticleChange.SetNumberOfSecondaries(tNH+nsec); 
	 for(G4int is=0; is<nsec; is++) aParticleChange.AddSecondary((*cascT)[is]);
  cascE->clear();
  delete cascE;
  cascT->clear();
  delete cascT;
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
    aParticleChange.AddSecondary( aNewTrack );
#ifdef debug
    G4cout<<"G4QCaptureAtRest::AtRestDoIt:#"<<i<<" is done."<<G4endl;
#endif
  }
  delete output;
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
#ifdef debug
  G4cout<<"G4QCaptureAtRest::RandomizeMuDecayOrCapture is called"<<G4endl;
#endif
  G4double Z27 =0.002727*Z;
  G4double Z227=Z27*Z27;
  G4double Z427=Z227*Z227;
  G4double Zeff=(Z-0.13782)*(1.2162-(0.09118-Z427)*std::sqrt((G4double)Z)); // EffNuclear Charge
  G4double Ze2=Zeff*Zeff;      // Squared effective charge of the Nucleus
  G4double pD=.00045516*(1.-Ze2*.00014658);// 1./MeanLifeTime of muon in atoms (in ns^-1)
  G4double pC=.00001637*Ze2*Ze2/(33.563+N);// 1./MeanLifeTime of muon NuclCapture(in ns^-1)
  if(Z==1&&N==0) pC=.0000007;
  if(Z==1&&N==1) pC=.000000012;
  G4double DLifeT=-std::log(G4UniformRand())/pD; // Time of the muon decay inside the atom
  G4double CLifeT=-std::log(G4UniformRand())/pC; // Time of the muon capture by nucleus
  if(DLifeT<CLifeT)
  {
    Time=DLifeT;
#ifdef debug
    G4cout<<"G4QCaptureAtRest::RandomizeMuDecayOrCapture: DecayLifeTime="<<Time<<G4endl;
#endif
    return false;
  }
  else
  {
    Time=CLifeT;
#ifdef debug
    G4cout<<"G4QCaptureAtRest::RandomizeMuDecayOrCapture: CaptureLifeTime="<<Time<<G4endl;
#endif
    return true;
  }
}

// Calculate the TotalEnergyDeposition for the AtomicCascadeDecay of MuMesoAtom to K-shell
void G4QCaptureAtRest::CalculateEnergyDepositionOfMuCapture(G4int Z) // (2p->1s) in MeV
{
  EnergyDeposition = .0029035*Z*Z*(1.-.0056817*Z)-.0006343; // MeV
#ifdef debug
  G4cout<<"G4QCaptureAtR::CalculateEnergyDepositionOfMuCapture="<<EnergyDeposition<<G4endl;
#endif
}

// Calculate gamma cascade from high (14th level) to the K(1s)-shell (follows V.Ivanchenko)
void G4QCaptureAtRest::MuCaptureEMCascade(G4int Z, G4int N, std::vector<G4double>* dV)
{ 
  static const G4double mEl  = G4Electron::Electron()->GetPDGMass(); // GEANT4 style
  static const G4double mMu  = G4MuonMinus::MuonMinus()->GetPDGMass();
  //static const G4double mEl  = G4QPDGCode(11).GetMass(); // CHIPS style
  //static const G4double mMu  = G4QPDGCode(13).GetMass();
  static const G4double vEl  = .0000136/mEl;
  //static const G4double dElM = mEl+mEl;
  // Inicialization - cascade start from 14th level (N.C.Mukhopadhyay Phy.Rep. 30 (1977) 1)
  G4double EnergyLevel[14];
  G4double dZ=Z;
  G4double nucM=G4NucleiProperties::GetNuclearMass(dZ+N,dZ);
  if(nucM<900.) nucM=G4QPDGCode(2112).GetNuclMass(Z,N,0); // CHIPS style

  G4double mass = mMu*nucM/(mMu+nucM);        //equivalemtMassOfMuon in C.M. muA-system

  G4double Z2=Z*Z;
  G4double KEnergy = vEl*Z2*mass;             // Finaite nuclear size (?)

  EnergyLevel[0] = EnergyDeposition;
#ifdef debug
  G4cout<<"G4QCapAtR::MuCapEMCascade:E="<<EnergyDeposition<<",e="<<mEl<<",m="<<mMu<<G4endl;
#endif
  for(G4int i=2; i<15; i++) EnergyLevel[i-1]=KEnergy/i/i; // To simple to be right (? M.K.)

  G4int nAuger = 1;
  G4int nGamma = 0;
  G4int nLevel = 13;
  G4double DeltaE=0.;
  G4double pGamma = Z2*Z2;

  // Capture on 14-th level
  G4double energy=EnergyLevel[13];
  //G4double ptot = sqrt(energy*(energy + dElM));
  //G4ThreeVector moment = ptot * RndmDir();
#ifdef debug
  G4cout<<"G4QCaptureAtR::MuCaptureEMCascade: first electron E="<<energy<<G4endl;
#endif
  dV->push_back(energy);
#ifdef debug
  G4cout<<"G4QCaptureAtR::MuCaptureEMCascade:before while nl="<<nLevel<<G4endl;
#endif
  // Algorithm of Vladimir Ivanchenko
  while(nLevel>0)                     // Radiative transitions and  Auger electron emission
  {
#ifdef debug
    G4cout<<"G4QCaptureAtR::MuCaptureEMCascade: in while nLevel="<<nLevel<<G4endl;
#endif
    // case of Auger electrons
    if((nAuger < Z) && ((pGamma + 10000.0) * G4UniformRand() < 10000.0) ) // 10000 (? M.K.)
    {
   	  nAuger++;                       // Radiate one more Auger electron
      DeltaE = EnergyLevel[nLevel-1] - EnergyLevel[nLevel];
      nLevel--;
#ifdef debug
	  G4cout<<"G4QCaptureAtR::MuCaptureEMCascade: Auger_e E="<<DeltaE<<G4endl;
#endif
      dV->push_back(DeltaE);
    }
    else // Rad transitions from C.S.Wu and L.Wilets, Ann. Rev. Nuclear Sci. 19 (1969) 527.
    {
      G4int iLevel = nLevel - 1 ;
      G4double var = 10.0 + iLevel * G4UniformRand(); // 10.0 (? M.K.)
      if(var > 10.0) iLevel -= G4int(var-10.0) + 1;
      if( iLevel < 0 ) iLevel = 0;
      DeltaE = EnergyLevel[iLevel] - EnergyLevel[nLevel];
      nLevel = iLevel;
#ifdef debug
	     G4cout<<"G4QCaptureAtR::MuCaptureEMCascade: photon E="<<DeltaE<<G4endl;
#endif
      dV->push_back(-DeltaE);
      nGamma++;
    }
  }
#ifdef debug
  G4cout<<"G4QCaptureAtR::MuCaptureEMCascade: nElect="<<nAuger<<", nGamm="<<nGamma<<G4endl;
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
  G4double Zeff=(Z-0.13782)*(1.2162-(0.09118-Z427)*std::sqrt((G4double)Z)); // EffNuclear Charge
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
void G4QCaptureAtRest::CalculateEnergyDepositionOfTauCapture(G4int Z) // (2p->1s) in MeV
{
  EnergyDeposition = .05*Z*Z*(1.-.0056817*Z)-.01; // MeV (@@ Must be improved)
#ifdef debug
  G4cout<<"G4QCapAtRest::CalculateEnergyDepositionOfTauCapture="<<EnergyDeposition<<G4endl;
#endif
}
