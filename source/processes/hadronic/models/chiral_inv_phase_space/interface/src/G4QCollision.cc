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
// $Id: G4QCollision.cc,v 1.11 2006/06/29 20:08:30 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
//      ---------------- G4QCollision class -----------------
//                 by Mikhail Kossov, December 2003.
// G4QCollision class of the CHIPS Simulation Branch in GEANT4
// ---------------------------------------------------------------
// ****************************************************************************************
// ********** This CLASS is temporary moved from the photolepton_hadron directory *********
// ****************************************************************************************

//#define debug
//#define pdebug

#include "G4QCollision.hh"

// Initialization of static vectors
std::vector<G4int> G4QCollision::ElementZ;            // Z of the element(i) in theLastCalc
std::vector<G4double> G4QCollision::ElProbInMat;      // SumProbabilityElements in Material
std::vector<std::vector<G4int>*> G4QCollision::ElIsoN;    // N of isotope(j) of Element(i)
std::vector<std::vector<G4double>*>G4QCollision::IsoProbInEl;//SumProbabIsotopes inElementI

G4QCollision::G4QCollision(const G4String& processName) : G4VDiscreteProcess(processName)
{
#ifdef debug
  G4cout<<"G4QCollision::Constructor is called"<<G4endl;
#endif
  if (verboseLevel>0) G4cout << GetProcessName() << " process is created "<< G4endl;

  G4QCHIPSWorld::Get()->GetParticles(nPartCWorld); // Create CHIPS World with 234 particles
  G4QNucleus::SetParameters(freeNuc,freeDib,clustProb,mediRatio); // Clusterization param's
  G4Quasmon::SetParameters(Temperature,SSin2Gluons,EtaEtaprime);  // Hadronic parameters
  G4QEnvironment::SetParameters(SolidAngle); // SolAngle of pbar-A secondary mesons capture
  //@@ Initialize here the G4QuasmonString parameters
}

G4bool   G4QCollision::manualFlag=false; // If false then standard parameters are used
G4double G4QCollision::Temperature=180.; // Critical Temperature (sensitive at High En)
G4double G4QCollision::SSin2Gluons=0.3;  // Supression of s-quarks (in respect to u&d)
G4double G4QCollision::EtaEtaprime=0.3;  // Supression of eta mesons (gg->qq/3g->qq)
G4double G4QCollision::freeNuc=0.5;      // Percentage of free nucleons on the surface
G4double G4QCollision::freeDib=0.05;     // Percentage of free diBaryons on the surface
G4double G4QCollision::clustProb=5.;     // Nuclear clusterization parameter
G4double G4QCollision::mediRatio=10.;    // medium/vacuum hadronization ratio
G4int    G4QCollision::nPartCWorld=152;  // The#of particles initialized in CHIPS World
G4double G4QCollision::SolidAngle=0.5;   // Part of Solid Angle to capture (@@A-dep.)
G4bool   G4QCollision::EnergyFlux=false; // Flag for Energy Flux use (not MultyQuasmon)
G4double G4QCollision::PiPrThresh=141.4; // Pion Production Threshold for gammas
G4double G4QCollision::M2ShiftVir=20000.;// Shift for M2=-Q2=m_pi^2 of the virtualGamma
G4double G4QCollision::DiNuclMass=1880.; // DoubleNucleon Mass for VirtualNormalization
G4double G4QCollision::photNucBias=1.;   // BiasingParameter for photo(e,mu,tau)Nuclear
G4double G4QCollision::weakNucBias=1.;   // BiasingParameter for ChargedCurrents(nu,mu) 

void G4QCollision::SetManual()   {manualFlag=true;}
void G4QCollision::SetStandard() {manualFlag=false;}

// Fill the private parameters
void G4QCollision::SetParameters(G4double temper, G4double ssin2g, G4double etaetap,
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

void G4QCollision::SetPhotNucBias(G4double phnB) {photNucBias=phnB;}
void G4QCollision::SetWeakNucBias(G4double ccnB) {weakNucBias=ccnB;}

// Destructor

G4QCollision::~G4QCollision() {}


G4LorentzVector G4QCollision::GetEnegryMomentumConservation()
{
  return EnMomConservation;
}

G4int G4QCollision::GetNumberOfNeutronsInTarget()
{
  return nOfNeutrons;
}

// output of the function must be in units of length! L=1/sig_V,sig_V=SUM(n(j,i)*sig(j,i)),
// where n(i,j) is a number of nuclei of the isotop j of the element i in V=1(lengtUnit^3)
// ********** All CHIPS cross sections are calculated in the surface units ************
G4double G4QCollision::GetMeanFreePath(const G4Track& aTrack,G4double,G4ForceCondition* Fc)
{
  *Fc = NotForced;
  const G4DynamicParticle* incidentParticle = aTrack.GetDynamicParticle();
  G4ParticleDefinition* incidentParticleDefinition=incidentParticle->GetDefinition();
  if( !IsApplicable(*incidentParticleDefinition))
    G4cout<<"-W-G4QCollision::GetMeanFreePath called for not implemented particle"<<G4endl;
  // Calculate the mean Cross Section for the set of Elements(*Isotopes) in the Material
  G4double Momentum = incidentParticle->GetTotalMomentum(); // 3-momentum of the Particle
  const G4Material* material = aTrack.GetMaterial();        // Get the current material
  const G4double* NOfNucPerVolume = material->GetVecNbOfAtomsPerVolume();
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QCollision::GetMeanFreePath:"<<nE<<" Elem's in theMaterial"<<G4endl;
#endif
  G4bool leptoNuc=false;       // By default the reaction is not lepto-nuclear
  G4VQCrossSection* CSmanager=0;
  G4int pPDG=0;
  if(incidentParticleDefinition == G4Proton::Proton())
  {
    CSmanager=G4QProtonNuclearCrossSection::GetPointer();
    pPDG=2212;
  }
  else if(incidentParticleDefinition == G4Gamma::Gamma())
  {
    CSmanager=G4QPhotonNuclearCrossSection::GetPointer();
    pPDG=22;
  }
  else if(incidentParticleDefinition == G4Electron::Electron() ||
          incidentParticleDefinition == G4Positron::Positron())
  {
    CSmanager=G4QElectronNuclearCrossSection::GetPointer();
    leptoNuc=true;
    pPDG=11;
  }
  else if(incidentParticleDefinition == G4MuonPlus::MuonPlus() ||
          incidentParticleDefinition == G4MuonMinus::MuonMinus())
  {
    CSmanager=G4QMuonNuclearCrossSection::GetPointer();
    leptoNuc=true;
    pPDG=13;
  }
  else if(incidentParticleDefinition == G4TauPlus::TauPlus() ||
          incidentParticleDefinition == G4TauMinus::TauMinus())
  {
    CSmanager=G4QTauNuclearCrossSection::GetPointer();
    leptoNuc=true;
    pPDG=15;
  }
  else if(incidentParticleDefinition == G4NeutrinoMu::NeutrinoMu() )
  {
    CSmanager=G4QNuMuNuclearCrossSection::GetPointer();
    leptoNuc=true;
    pPDG=14;
  }
  else if(incidentParticleDefinition == G4AntiNeutrinoMu::AntiNeutrinoMu() )
  {
    CSmanager=G4QANuMuNuclearCrossSection::GetPointer();
    leptoNuc=true;
    pPDG=-14;
  }
  else G4cout<<"G4QCollision::GetMeanFreePath:Particle isn't implemented in CHIPS"<<G4endl;
  
  G4QIsotope* Isotopes = G4QIsotope::Get(); // Pointer to the G4QIsotopes singleton
  G4double sigma=0.;                        // Sums over elements for the material
  G4int IPIE=IsoProbInEl.size();            // How many old elements?
  if(IPIE) for(G4int ip=0; ip<IPIE; ++ip)   // Clean up the SumProb's of Isotopes (SPI)
  {
    std::vector<G4double>* SPI=IsoProbInEl[ip]; // Pointer to the SPI vector
    SPI->clear();
    delete SPI;
    std::vector<G4int>* IsN=ElIsoN[ip];     // Pointer to the N vector
    IsN->clear();
    delete IsN;
  }
  ElProbInMat.clear();                      // Clean up the SumProb's of Elements (SPE)
  ElementZ.clear();                         // Clear the body vector for Z of Elements
  IsoProbInEl.clear();                      // Clear the body vector for SPI
  ElIsoN.clear();                           // Clear the body vector for N of Isotopes
  for(G4int i=0; i<nE; ++i)
  {
    G4Element* pElement=(*theElementVector)[i]; // Pointer to the current element
    G4int Z = static_cast<G4int>(pElement->GetZ()); // Z of the Element
    ElementZ.push_back(Z);                  // Remember Z of the Element
    G4int isoSize=0;                        // The default for the isoVectorLength is 0
    G4int indEl=0;                          // Index of non-trivial element or 0(default)
    G4IsotopeVector* isoVector=pElement->GetIsotopeVector(); // Get the predefined IsoVect
    if(isoVector) isoSize=isoVector->size();// Get size of the existing isotopeVector
#ifdef debug
    G4cout<<"G4QCollision::GetMeanFreePath: isovectorLength="<<isoSize<<G4endl; // Result
#endif
    if(isoSize)                             // The Element has non-trivial abumdance set
    {
      indEl=pElement->GetIndex();           // Index of the non-trivial element
      if(!Isotopes->IsDefined(Z,indEl))     // This index is not defined for this Z: define
      {
        std::vector<std::pair<G4int,G4double>*>* newAbund =
                                               new std::vector<std::pair<G4int,G4double>*>;
        G4double* abuVector=pElement->GetRelativeAbundanceVector();
        for(G4int j=0; j<isoSize; j++)      // Calculation of abundance vector for isotopes
        {
          G4int N=pElement->GetIsotope(j)->GetN()-Z; // N means A=N+Z !
          if(pElement->GetIsotope(j)->GetZ()!=Z)G4cerr<<"G4QCaptureAtRest::GetMeanFreePath"
																																	<<": Z="<<pElement->GetIsotope(j)->GetZ()<<"#"<<Z<<G4endl;
          G4double abund=abuVector[j];
								  std::pair<G4int,G4double>* pr= new std::pair<G4int,G4double>(N,abund);
#ifdef debug
          G4cout<<"G4QCollision::PostStepDoIt:pair#="<<j<<", N="<<N<<",ab="<<abund<<G4endl;
#endif
          newAbund->push_back(pr);
						  }
#ifdef debug
        G4cout<<"G4QCollision::PostStepDoIt: pairVectorLength="<<newAbund->size()<<G4endl;
#endif
        indEl=G4QIsotope::Get()->InitElement(Z,indEl,newAbund); // definition of the newInd
        for(G4int k=0; k<isoSize; k++) delete (*newAbund)[k];   // Cleaning temporary
        delete newAbund; // Was "new" in the beginning of the name space
      }
    }
    std::vector<std::pair<G4int,G4double>*>* cs= Isotopes->GetCSVector(Z,indEl);//CSPointer
    std::vector<G4double>* SPI = new std::vector<G4double>; // Pointer to the SPI vector
    IsoProbInEl.push_back(SPI);
    std::vector<G4int>* IsN = new std::vector<G4int>; // Pointer to the N vector
    ElIsoN.push_back(IsN);
    G4int nIs=cs->size();                   // A#Of Isotopes in the Element
    G4double susi=0.;                       // sum of CS over isotopes
    if(nIs) for(G4int j=0; j<nIs; j++)      // Calculate CS for eachIsotope of El
    {
      std::pair<G4int,G4double>* curIs=(*cs)[j]; // A pointer, which is used twice
      G4int N=curIs->first;                 // #of Neuterons in the isotope j of El i
      IsN->push_back(N);                    // Remember Min N for the Element
      G4double CSI=CSmanager->GetCrossSection(true,Momentum,Z,N,pPDG);//CS(j,i) for isotope
      curIs->second = CSI;
      susi+=CSI;                            // Make a sum per isotopes
      SPI->push_back(susi);                 // Remember summed cross-section
    } // End of temporary initialization of the cross sections in the G4QIsotope singeltone
    sigma+=Isotopes->GetMeanCrossSection(Z,indEl)*NOfNucPerVolume[i];//SUM(MeanCS*NOfNperV)
    ElProbInMat.push_back(sigma);
  } // End of LOOP over Elements

  // Check that cross section is not zero and return the mean free path
  if(photNucBias!=1.) if(incidentParticleDefinition == G4MuonPlus::MuonPlus()   ||
                         incidentParticleDefinition == G4MuonMinus::MuonMinus() ||
                         incidentParticleDefinition == G4Electron::Electron()   ||
                         incidentParticleDefinition == G4Positron::Positron()   ||  
                         incidentParticleDefinition == G4TauMinus::TauMinus()   ||
                         incidentParticleDefinition == G4TauPlus::TauPlus()       )
                                                                        sigma*=photNucBias;
  if(photNucBias!=1.) if(incidentParticleDefinition==G4NeutrinoE::NeutrinoE()            ||
                         incidentParticleDefinition==G4AntiNeutrinoE::AntiNeutrinoE()    ||
                         incidentParticleDefinition==G4NeutrinoTau::NeutrinoTau()        ||
                         incidentParticleDefinition==G4AntiNeutrinoTau::AntiNeutrinoTau()||
                         incidentParticleDefinition==G4NeutrinoMu::NeutrinoMu()          ||
                         incidentParticleDefinition==G4AntiNeutrinoMu::AntiNeutrinoMu()   )
                                                                        sigma*=weakNucBias;
  if(sigma > 0.) return 1./sigma;                 // Mean path [distance] 
  return DBL_MAX;
}


G4bool G4QCollision::IsApplicable(const G4ParticleDefinition& particle) 
{
  if      (particle == *(      G4MuonPlus::MuonPlus()      )) return true;
  else if (particle == *(     G4MuonMinus::MuonMinus()     )) return true; 
  else if (particle == *(       G4TauPlus::TauPlus()       )) return true;
  else if (particle == *(      G4TauMinus::TauMinus()      )) return true;
  else if (particle == *(      G4Electron::Electron()      )) return true;
  else if (particle == *(      G4Positron::Positron()      )) return true;
  else if (particle == *(         G4Gamma::Gamma()         )) return true;
  //else if (particle == *(        G4Proton::Proton()        )) return true;
  else if (particle == *(G4AntiNeutrinoMu::AntiNeutrinoMu())) return true;
  else if (particle == *(   G4NeutrinoMu::NeutrinoMu()   )) return true;
  //else if (particle == *(       G4Neutron::Neutron()       )) return true;
  //else if (particle == *(     G4PionMinus::PionMinus()     )) return true;
  //else if (particle == *(      G4PionPlus::PionPlus()      )) return true;
  //else if (particle == *(      G4KaonPlus::KaonPlus()      )) return true;
  //else if (particle == *(     G4KaonMinus::KaonMinus()     )) return true;
  //else if (particle == *(  G4KaonZeroLong::KaonZeroLong()  )) return true;
  //else if (particle == *( G4KaonZeroShort::KaonZeroShort() )) return true;
  //else if (particle == *(        G4Lambda::Lambda()        )) return true;
  //else if (particle == *(     G4SigmaPlus::SigmaPlus()     )) return true;
  //else if (particle == *(    G4SigmaMinus::SigmaMinus()    )) return true;
  //else if (particle == *(     G4SigmaZero::SigmaZero()     )) return true;
  //else if (particle == *(       G4XiMinus::XiMinus()       )) return true;
  //else if (particle == *(        G4XiZero::XiZero()        )) return true;
  //else if (particle == *(    G4OmegaMinus::OmegaMinus()    )) return true;
  //else if (particle == *(   G4AntiNeutron::AntiNeutron()   )) return true;
  //else if (particle == *(    G4AntiProton::AntiProton()    )) return true;
#ifdef debug
  G4cout<<"***G4QCollision::IsApplicable: PDG="<<particle.GetPDGEncoding()<<G4endl;
#endif
  return false;
}

G4VParticleChange* G4QCollision::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  static const G4double mu=G4MuonMinus::MuonMinus()->GetPDGMass(); // muon mass
  static const G4double mu2=mu*mu;                                 // squared muon mass
  //static const G4double dpi=M_PI+M_PI;   // 2*pi (for Phi distr.) ***changed to twopi***
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mNeut2= mNeut*mNeut;
  static const G4double muN= mNeut+mu;
  static const G4double muN2= muN*muN;
  static const G4double fmuN= 4*mNeut2*mu2;
  static const G4double musN= mNeut2+mu2;
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mProt2= mProt*mProt;
  static const G4double muP= mProt+mu;
  static const G4double muP2= muP*muP;
  static const G4double fmuP= 4*mProt2*mu2;
  static const G4double musP= mProt2+mu2;
  static const G4double dM=mProt+mNeut;                            // doubled nucleon mass
  static const G4double mudM=mu2/dM;                               // for x limit
  static const G4double hdM=dM/2.;                                 // M of the "nucleon"
  //static const G4double hdM2=hdM*hdM;                              // M2 of the "nucleon"
  //static const G4double mPi0 = G4QPDGCode(111).GetMass();
  //static const G4double mDeut= G4QPDGCode(2112).GetNuclMass(1,1,0);
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double mPPi = mPi+mProt;   // Delta threshold
  //static const G4double mPPi2= mPPi*mPPi; // Delta low threshold for W2
  //static const G4double mDel2= 1400*1400;   // Delta up threshold for W2 (in MeV^2)
  static const G4double muD  = mPPi+mu;     // Multiperipheral threshold
  static const G4double muD2 = muD*muD;
  //static const G4double mMu  = G4QPDGCode(13).GetMass();
  //static const G4double mTau = G4QPDGCode(15).GetMass();
  //static const G4double mEl  = G4QPDGCode(11).GetMass();
  //
  const G4DynamicParticle* projHadron = track.GetDynamicParticle();
  const G4ParticleDefinition* particle=projHadron->GetDefinition();
#ifdef debug
  G4cout<<"G4QCollision::PostStepDoIt: Before the GetMeanFreePath is called"<<G4endl;
#endif
  G4ForceCondition cond=NotForced;
  GetMeanFreePath(track, 1., &cond);
#ifdef debug
  G4cout<<"G4QCollision::PostStepDoIt: After the GetMeanFreePath is called"<<G4endl;
#endif
  G4bool scat=false;
  G4int  scatPDG=0;                                   // Must be filled if true
  G4LorentzVector proj4M=projHadron->Get4Momentum();
  G4LorentzVector scat4M=proj4M;                      // Must be filled if true
  G4double momentum = projHadron->GetTotalMomentum(); // 3-momentum of the Particle
  G4double Momentum=proj4M.rho();
  if(std::fabs(Momentum-momentum)>.001)
    G4cerr<<"G4QCollision::PostStepDoIt: P="<<Momentum<<"="<<momentum<<G4endl;
#ifdef debug
  G4double mp=proj4M.m();
  G4cout<<"G4QCollision::PostStepDoIt is called, P="<<Momentum<<"="<<momentum<<G4endl;
#endif
  if (!IsApplicable(*particle))  // Check applicability
  {
    G4cerr<<"G4QCollision::PostStepDoIt:Only gam,e+,e-,mu+,mu-,t+,t-,p are implemented."
          <<G4endl;
    return 0;
  }
  const G4Material* material = track.GetMaterial();      // Get the current material
  G4int Z=0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QCollision::PostStepDoIt: "<<nE<<" elements in the material."<<G4endl;
#endif
  G4int projPDG=0;                           // PDG Code prototype for the captured hadron
  // Not all these particles are implemented yet (see Is Applicable)
  if      (particle ==        G4MuonPlus::MuonPlus()       ) projPDG=  -13;
  else if (particle ==       G4MuonMinus::MuonMinus()      ) projPDG=   13;
  else if (particle ==      G4NeutrinoMu::NeutrinoMu()     ) projPDG=   14;
  else if (particle ==  G4AntiNeutrinoMu::AntiNeutrinoMu() ) projPDG=  -14;
  else if (particle ==        G4Electron::Electron()       ) projPDG=   11;
  else if (particle ==        G4Positron::Positron()       ) projPDG=  -11;
  else if (particle ==       G4NeutrinoE::NeutrinoE()      ) projPDG=   12;
  else if (particle ==   G4AntiNeutrinoE::AntiNeutrinoE()  ) projPDG=  -12;
  else if (particle ==           G4Gamma::Gamma()          ) projPDG=   22;
  else if (particle ==          G4Proton::Proton()         ) projPDG= 2212;
  else if (particle ==         G4Neutron::Neutron()        ) projPDG= 2112;
  else if (particle ==       G4PionMinus::PionMinus()      ) projPDG= -211;
  else if (particle ==        G4PionPlus::PionPlus()       ) projPDG=  211;
  else if (particle ==        G4KaonPlus::KaonPlus()       ) projPDG= 2112;
  else if (particle ==       G4KaonMinus::KaonMinus()      ) projPDG= -321;
  else if (particle ==    G4KaonZeroLong::KaonZeroLong()   ) projPDG=  130;
  else if (particle ==   G4KaonZeroShort::KaonZeroShort()  ) projPDG=  310;
  else if (particle ==         G4TauPlus::TauPlus()        ) projPDG=  -15;
  else if (particle ==        G4TauMinus::TauMinus()       ) projPDG=   15;
  else if (particle ==     G4NeutrinoTau::NeutrinoTau()    ) projPDG=   16;
  else if (particle == G4AntiNeutrinoTau::AntiNeutrinoTau()) projPDG=  -16;
  else if (particle ==          G4Lambda::Lambda()         ) projPDG= 3122;
  else if (particle ==       G4SigmaPlus::SigmaPlus()      ) projPDG= 3222;
  else if (particle ==      G4SigmaMinus::SigmaMinus()     ) projPDG= 3112;
  else if (particle ==       G4SigmaZero::SigmaZero()      ) projPDG= 3212;
  else if (particle ==         G4XiMinus::XiMinus()        ) projPDG= 3312;
  else if (particle ==          G4XiZero::XiZero()         ) projPDG= 3322;
  else if (particle ==      G4OmegaMinus::OmegaMinus()     ) projPDG= 3334;
  else if (particle ==     G4AntiNeutron::AntiNeutron()    ) projPDG=-2112;
  else if (particle ==      G4AntiProton::AntiProton()     ) projPDG=-2212;
  G4int aProjPDG=std::abs(projPDG);
#ifdef debug
  G4int prPDG=particle->GetPDGEncoding();
		G4cout<<"G4QCollision::PostStepDoIt: projPDG="<<projPDG<<", stPDG="<<prPDG<<G4endl;
#endif
  if(!projPDG)
  {
    G4cerr<<"---Warning---G4QCollision::PostStepDoIt:Undefined interacting hadron"<<G4endl;
    return 0;
  }
  G4int EPIM=ElProbInMat.size();
#ifdef debug
		G4cout<<"G4QCollis::PostStDoIt: m="<<EPIM<<",n="<<nE<<",T="<<ElProbInMat[EPIM-1]<<G4endl;
#endif
  G4int i=0;
  if(EPIM>1)
  {
    G4double rnd = ElProbInMat[EPIM-1]*G4UniformRand();
    for(i=0; i<nE; ++i)
		  {
#ifdef debug
				  G4cout<<"G4QCollision::PostStepDoIt:E["<<i<<"]="<<ElProbInMat[i]<<",r="<<rnd<<G4endl;
#endif
      if (rnd<ElProbInMat[i]) break;
    }
    if(i>=nE) i=nE-1;                        // Top limit for the Element
  }
  G4Element* pElement=(*theElementVector)[i];
  Z=static_cast<G4int>(pElement->GetZ());
#ifdef debug
				G4cout<<"G4QCollision::PostStepDoIt: i="<<i<<", Z(element)="<<Z<<G4endl;
#endif
  if(Z<=0)
  {
    G4cerr<<"---Warning---G4QCollision::PostStepDoIt: Element with Z="<<Z<<G4endl;
    if(Z<0) return 0;
  }
  std::vector<G4double>* SPI = IsoProbInEl[i];// Vector of summedProbabilities for isotopes
  std::vector<G4int>* IsN = ElIsoN[i];     // Vector of "#of neutrons" in the isotope El[i]
  G4int nofIsot=SPI->size();               // #of isotopes in the element i
#ifdef debug
		G4cout<<"G4QCollis::PosStDoIt:n="<<nofIsot<<",T="<<(*SPI)[nofIsot-1]<<",r="<<rnd<<G4endl;
#endif
  G4int j=0;
  if(nofIsot>1)
  {
    G4double rndI=(*SPI)[nofIsot-1]*G4UniformRand(); // Randomize the isotop of the Element
    for(j=0; j<nofIsot; ++j)
    {
#ifdef debug
				  G4cout<<"G4QCollision::PostStepDoIt: SP["<<j<<"]="<<(*SPI)[j]<<", r="<<rndI<<G4endl;
#endif
      if(rndI < (*SPI)[j]) break;
    }
    if(j>=nofIsot) j=nofIsot-1;            // Top limit for the isotope
  }
  G4int N =(*IsN)[j]; ;                    // Randomized number of neutrons
#ifdef debug
		G4cout<<"G4QCollision::PostStepDoIt: j="<<i<<", N(isotope)="<<N<<G4endl;
#endif
  if(N<0)
  {
    G4cerr<<"-Warning-G4QCollision::PostStepDoIt: Isotope with Z="<<Z<<", 0>N="<<N<<G4endl;
    return 0;
  }
  nOfNeutrons=N;                           // Remember it for the energy-momentum check
  G4double dd=0.025;
  G4double am=Z+N;
  G4double sr=std::sqrt(am);
  G4double dsr=0.01*(sr+sr);
  if(dsr<dd)dsr=dd;
  if(manualFlag) G4QNucleus::SetParameters(freeNuc,freeDib,clustProb,mediRatio); // ManualP
		else if(projPDG==-2212) G4QNucleus::SetParameters(1.-dsr-dsr,dd+dd,5.,10.);//aP ClustPars
  else if(projPDG==-211)  G4QNucleus::SetParameters(.67-dsr,.32-dsr,5.,9.);//Pi- ClustPars
#ifdef debug
  G4cout<<"G4QCollision::PostStepDoIt: N="<<N<<" for element with Z="<<Z<<G4endl;
#endif
  if(N<0)
  {
    G4cerr<<"---Warning---G4QCollision::PostStepDoIt:Element with N="<<N<< G4endl;
    return 0;
  }
  G4int targPDG=90000000+Z*1000+N;       // PDG Code of the target nucleus
  G4QPDGCode targQPDG(targPDG);
  G4double tM=targQPDG.GetMass();
		if(aProjPDG==11 || aProjPDG==13 || aProjPDG==15) // leptons with photonuclear
		{ // Lepto-nuclear case with the equivalent photon algorithm. @@InFuture + neutrino & QE
    G4double kinEnergy= projHadron->GetKineticEnergy();
    G4ParticleMomentum dir = projHadron->GetMomentumDirection();
    G4VQCrossSection* CSmanager=G4QElectronNuclearCrossSection::GetPointer();
    if(aProjPDG== 13) CSmanager=G4QMuonNuclearCrossSection::GetPointer();
    if(aProjPDG== 15) CSmanager=G4QTauNuclearCrossSection::GetPointer();
    // @@ Probably this is not necessary any more
    G4double xSec=CSmanager->GetCrossSection(false,Momentum,Z,N);//Recalculate CrossSection
    // @@ check a possibility to separate p, n, or alpha (!)
    G4double photonEnergy = CSmanager->GetExchangeEnergy(); // Energy of EqivExchangePart
    if(xSec <= 0.) // The cross-section iz 0 -> Do Nothing
    {
      G4cerr<<"-Warning-G4QCollision::PSDoIt: IsStillCalled photE="<<photonEnergy<<G4endl;
      //Do Nothing Action insead of the reaction
      aParticleChange.ProposeEnergy(kinEnergy);
      aParticleChange.ProposeLocalEnergyDeposit(0.);
      aParticleChange.ProposeMomentumDirection(dir) ;
      return G4VDiscreteProcess::PostStepDoIt(track,step);
    }
    if( kinEnergy < photonEnergy )
    {
      //Do Nothing Action insead of the reaction
      G4cerr<<"G4QCollision::PSDoIt: photE="<<photonEnergy<<">leptE="<<kinEnergy<<G4endl;
      aParticleChange.ProposeEnergy(kinEnergy);
      aParticleChange.ProposeLocalEnergyDeposit(0.);
      aParticleChange.ProposeMomentumDirection(dir) ;
      return G4VDiscreteProcess::PostStepDoIt(track,step);
    }
    G4double photonQ2 = CSmanager->GetExchangeQ2(photonEnergy);// Q2(t) of EqivExchangePart
    G4double W=photonEnergy-photonQ2/dM;// HadronicEnergyFlow (W-energy) for virtual photon
    if(W<0.) 
    {
      //Do Nothing Action insead of the reaction
      G4cout << "G4QCollision::PostStepDoIt:(lN) negative equivalent energy W="<<W<<G4endl;
      aParticleChange.ProposeEnergy(kinEnergy);
      aParticleChange.ProposeLocalEnergyDeposit(0.);
      aParticleChange.ProposeMomentumDirection(dir) ;
      return G4VDiscreteProcess::PostStepDoIt(track,step);
    }
    // Update G4VParticleChange for the scattered muon
    G4VQCrossSection* thePhotonData=G4QPhotonNuclearCrossSection::GetPointer();
    G4double sigNu=thePhotonData->GetCrossSection(true,photonEnergy,Z,N);// IntegratedCrSec
    G4double sigK =thePhotonData->GetCrossSection(true, W, Z, N);        // Real CrossSect.
    G4double rndFraction = CSmanager->GetVirtualFactor(photonEnergy, photonQ2);
    if(sigNu*G4UniformRand()>sigK*rndFraction) 
    {
      //Do NothingToDo Action insead of the reaction
      G4cout << "G4QCollision::PostStepDoIt: probability correction - DoNothing"<<G4endl;
      aParticleChange.ProposeEnergy(kinEnergy);
      aParticleChange.ProposeLocalEnergyDeposit(0.);
      aParticleChange.ProposeMomentumDirection(dir) ;
      return G4VDiscreteProcess::PostStepDoIt(track,step);
    }
    G4double iniE=kinEnergy+mu;          // Initial total energy of the muon
    G4double finE=iniE-photonEnergy;     // Final total energy of the muon
    if(finE>0) aParticleChange.ProposeEnergy(finE) ;
    else
    {
      aParticleChange.ProposeEnergy(0.) ;
      aParticleChange.ProposeTrackStatus(fStopAndKill);
    }
    // Scatter the muon
    G4double EEm=iniE*finE-mu2;          // Just an intermediate value to avoid "2*"
    G4double iniP=std::sqrt(iniE*iniE-mu2);   // Initial momentum of the electron
    G4double finP=std::sqrt(finE*finE-mu2);   // Final momentum of the electron
    G4double cost=(EEm+EEm-photonQ2)/iniP/finP; // cos(theta) for the electron scattering
    if(cost>1.) cost=1.;                 // To avoid the accuracy of calculation problem
    //else if(cost>1.001)                // @@ error report can be done, but not necessary
    if(cost<-1.) cost=-1.;               // To avoid the accuracy of calculation problem
    //else if(cost<-1.001)               // @@ error report can be done, but not necessary
    // --- Example from electromagnetic physics --
    //G4ThreeVector newMuonDirection(dirx,diry,dirz);
    //newMuonDirection.rotateUz(dir);
    //aParticleChange.ProposeMomentumDirection(newMuonDirection1) ;
    // The scattering in respect to the derection of the incident muon is made impicitly:
    G4ThreeVector ort=dir.orthogonal();  // Not normed orthogonal vector (!) (to dir)
    G4ThreeVector ortx = ort.unit();     // First unit vector orthogonal to the direction
    G4ThreeVector orty = dir.cross(ortx);// Second unit vector orthoganal to the direction
    G4double sint=std::sqrt(1.-cost*cost);    // Perpendicular component
    G4double phi=twopi*G4UniformRand();  // phi of scattered electron
    G4double sinx=sint*std::sin(phi);         // x-component
    G4double siny=sint*std::cos(phi);         // y-component
    G4ThreeVector findir=cost*dir+sinx*ortx+siny*orty;
    aParticleChange.ProposeMomentumDirection(findir); // new direction for the muon
    const G4ThreeVector photon3M=iniP*dir-finP*findir;
    projPDG=22;
    proj4M=G4LorentzVector(photon3M,photon3M.mag()); //@@ photon is real?
  }
		else if(aProjPDG==14)// ** neutrino nuclear interactions (only nu_mu/anu_mu & only CC) **
		{
    G4double kinEnergy= projHadron->GetKineticEnergy();// For neutrino this is total energy
    G4double dKinE=kinEnergy+kinEnergy;  // doubled energy for s calculation
    G4ParticleMomentum dir = projHadron->GetMomentumDirection(); // unit vector
    G4VQCrossSection* CSmanager=G4QNuMuNuclearCrossSection::GetPointer();
    proj4M=G4LorentzVector(dir*kinEnergy,kinEnergy);   // temporary
    G4bool nuanu=true;
    scatPDG=13;                          // Prototype = secondary scattered mu-
    if(projPDG==-14)
    {
      nuanu=false;
      CSmanager=G4QANuMuNuclearCrossSection::GetPointer(); // @@ open
      scatPDG=-13;                       // secondary scattered mu+
    }
    // @@ Probably this is not necessary any more
    G4double xSec=CSmanager->GetCrossSection(false,Momentum,Z,N);//Recalculate CrossSection
    // @@ check a possibility to separate p, n, or alpha (!)
    if(xSec <= 0.) // The cross-section = 0 -> Do Nothing
    {
      G4cerr<<"-Warning-G4QCollision::PSDoIt: IsStillCalled nuE="<<kinEnergy<<G4endl;
      //Do Nothing Action insead of the reaction
      aParticleChange.ProposeEnergy(kinEnergy);
      aParticleChange.ProposeLocalEnergyDeposit(0.);
      aParticleChange.ProposeMomentumDirection(dir) ;
      return G4VDiscreteProcess::PostStepDoIt(track,step);
    }
    scat=true;                                  // event with changed scattered projectile
    G4double totCS = CSmanager->GetLastTOTCS(); // the last total cross section (isotope?)
    if(std::fabs(xSec-totCS)/xSec>.0001)
          G4cout<<"-Warning-G4QCollision::PostStepDoIt: xS="<<xSec<<"# CS="<<totCS<<G4endl;
    G4double qelCS = CSmanager->GetLastQELCS(); // the last total cross section
    if(totCS - qelCS < 0.) totCS = qelCS;       // only at low energies
    // make different definitions for neutrino and antineutrino (inFuture for nue/anue too)
    G4double mIN=mProt;                         // Just a prototype (for anu, Z=1, N=0)
    G4double mOT=mNeut;
    G4double OT=muN2;
    G4double mOT2=mNeut2;
    G4double muOT=fmuN;
    G4double musOT=musN;
    if(nuanu)
    {
      targPDG-=1;
      G4QPDGCode targQPDG(targPDG);
      G4double rM=targQPDG.GetMass();
      mIN=tM-rM;                                // bounded in-mass of the neutron
      tM=rM;
      mOT=mProt;
      OT=muP2;
      mOT2=mProt2;
      muOT=fmuP;
      musOT=musP;
      projPDG=2212;                             // proton is going out
    }
    else
    {
      if(Z>1||N>0)                              // Calculate the splitted mass
						{
        targPDG-=1000;
        G4QPDGCode targQPDG(targPDG);
        G4double rM=targQPDG.GetMass();
        mIN=tM-rM;                              // bounded in-mass of the proton
        tM=rM;
      }
      else targPDG=0;
      projPDG=2112;                             // neutron is going out
    }
    G4double s=mIN*(mIN+dKinE);                 // s=(M_cm)^2
    if(s<=OT)                                   // *** Do nothing solution ***
    {
      //Do NothingToDo Action insead of the reaction (@@ Can we make it common?)
      G4cout << "G4QCollision::PostStepDoIt: probability correction - DoNothing"<<G4endl;
      aParticleChange.ProposeEnergy(kinEnergy);
      aParticleChange.ProposeLocalEnergyDeposit(0.);
      aParticleChange.ProposeMomentumDirection(dir) ;
      return G4VDiscreteProcess::PostStepDoIt(track,step);
    }
    if((!nuanu||N)&&totCS*G4UniformRand()<qelCS||s<muD2)// ****** Quasi-Elastic interaction
    {
      G4double Q2=CSmanager->GetQEL_ExchangeQ2(); // OK, im MeV^2
      G4double ds=s+s;                          // dpubled s
      G4double sqs=std::sqrt(s);                // M_cm
      G4double pi=(s-mIN*mIN)/(sqs+sqs);        // initial momentum in CMS
      G4double dpi=pi+pi;                       // doubled initial momentum in CMS
      G4double sd=s-musOT;                      // s-mu2-mOT2
      G4double qo2=(sd*sd-muOT)/(ds+ds);        // squared momentum of secondaries in CMS
      G4double qo=std::sqrt(qo2);               // momentum of secondaries in CMS
      G4double cost=(dpi*std::sqrt(qo2+mu2)-Q2-mu2)/dpi/qo; // cos(theta) in CMS
      G4LorentzVector t4M(0.,0.,0.,mIN);        // 4mom of the effective target
      G4LorentzVector c4M=t4M+proj4M;           // 4mom of the compound system
      t4M.setT(mOT);                            // now it is 4mom of the outgoing nucleon
      scat4M=G4LorentzVector(0.,0.,0.,mu);      // 4mom of the scattered muon
      if(!G4QHadron(c4M).RelDecayIn2(scat4M, t4M, proj4M, cost, cost))
      {
        G4cerr<<"G4QCol::PSD:c4M="<<c4M<<sqs<<",mM="<<mu<<",tM="<<mOT<<",c="<<cost<<G4endl;
        throw G4QException("G4QCollision::HadronizeQuasm: Can't dec QE nu,mu Compound");
      }
      proj4M=t4M;                               // 4mom of the new projectile nucleon
    }
    else                                        // ***** Non Quasi Elastic interaction
    {
      G4double Q2=CSmanager->GetNQE_ExchangeQ2();
      projPDG=CSmanager->GetExchangePDGCode();
      //@@ Temporary made only for direct interaction and for N=3 (good for small Q2)
      //@@ inFuture use N=GetNPartons and directFraction=GetDirectPart, @@ W2...
      G4double r=G4UniformRand();
      G4double r1=0.5;                                  // (1-x)
      if(r<0.5)      r1=std::sqrt(r+r)*(.5+.1579*(r-.5));
      else if(r>0.5) r1=1.-std::sqrt(2.-r-r)*(.5+.1579*(.5-r));
      G4double xn=1.-mudM/Momentum;             // Normalization of (1-x) [x>mudM/Mom]
      G4double x1=xn*r1;                        // (1-x)
      G4double x=1.-x1;                         // x=2k/M
      //G4double W2=(hdM2+Q2/x)*x1;               // W2 candidate
      G4double mx=hdM*x;
      G4double we=(Q2/mx-mx)/2;                 // transfered energy
      G4double muQ2=(mu2+Q2)/2;
      G4double cost=(kinEnergy*we-muQ2)/kinEnergy/std::sqrt(we*we+Q2);
      if(std::fabs(cost)>1)
      {
        if(cost>1.) cost=1.;
        else        cost=-1.;
        we=(muQ2*muQ2-kinEnergy*kinEnergy*Q2)/kinEnergy/(mu2+Q2); // minimum we
      }
      scat4M=G4LorentzVector(0.,0.,0.,mu);      // 4mom of the scattered muon
      G4LorentzVector t4M(0.,0.,0.,-Q2);        // 4mom of the virtual W
      G4LorentzVector dir4M=proj4M-G4LorentzVector(0.,0.,0.,proj4M.e()*.1);// projDirection
      if(!G4QHadron(proj4M).RelDecayIn2(scat4M, t4M, dir4M, cost, cost))
      {
        G4cerr<<"G4QCol::PSD:4M="<<proj4M<<",mM="<<mu<<",Q2="<<Q2<<",c="<<cost<<G4endl;
        throw G4QException("G4Quasmon::HadronizeQuasm: Can't dec nu->mu+W");
      }
      proj4M=t4M;                               // 4m of the pion
    }
    aParticleChange.ProposeEnergy(0.) ;
    aParticleChange.ProposeTrackStatus(fStopAndKill); // the initial neutrino is killed
  }
  EnMomConservation=proj4M+G4LorentzVector(0.,0.,0.,tM);    // Total 4-mom of the reaction
  G4QHadronVector* output=new G4QHadronVector; // Prototype of the output G4QHadronVector
  // @@@@@@@@@@@@@@ Temporary for the testing purposes --- Begin
  //G4bool elF=false; // Flag of the ellastic scattering is "false" by default
  //G4double eWei=1.;
  // @@@@@@@@@@@@@@ Temporary for the testing purposes --- End  
#ifdef debug
  G4cout<<"G4QCollision::PostStepDoIt: projPDG="<<projPDG<<", targPDG="<<targPDG<<G4endl;
#endif
  G4QHadron* pH = new G4QHadron(projPDG,proj4M);                // ---> DELETED -->--  -+
  //if(momentum<1000.) // Condition for using G4QEnvironment (not G4QuasmonString)      |
		{ //                                                                                  |
    G4QHadronVector projHV;                                 //                          |
    projHV.push_back(pH);                                   // DESTROYED over 2 lines-+ |
    G4QEnvironment* pan= new G4QEnvironment(projHV,targPDG);// ---> DELETED --->----+ | |
    std::for_each(projHV.begin(), projHV.end(), DeleteQHadron()); // <---<------<---+-+-+
    projHV.clear(); // <------------<---------------<-------------------<-----------+-+ .
#ifdef debug
    G4cout<<"G4QCollision::PostStepDoIt: pPDG="<<projPDG<<", mp="<<mp<<G4endl; //   |   .
#endif
    try                                                           //                |   .
	   {                                                             //                |   .
	     delete output;                                              //                |   .
      output = pan->Fragment();// DESTROYED in the end of the LOOP work space       |   .
    }                                                             //                |   .
    catch (G4QException& error)//                                                   |   .
	   {                                                             //                |   .
	     //#ifdef pdebug
      G4cerr<<"***G4QCollision::PostStepDoIt: G4QE Exception is catched"<<G4endl;// |   .
	     //#endif
      G4Exception("G4QCollision::PostStepDoIt:","27",FatalException,"CHIPSCrash");//|   .
    }                                                             //                |   .
    delete pan;                              // Delete the Nuclear Environment <-<--+   .
  } //                                                                                  .
  //else               // Use G4QuasmonString                                             .
		//{ //                                                                                  ^
  //  G4QuasmonString* pan= new G4QuasmonString(pH,false,targPDG,false);//-> DELETED --+  |
  //  delete pH;                                                    // --------<-------+--+
#ifdef debug
  //  G4double mp=G4QPDGCode(projPDG).GetMass();   // Mass of the projectile particle  |
  //  G4cout<<"G4QCollision::PostStepDoIt: pPDG="<<projPDG<<", pM="<<mp<<G4endl; //    |
#endif
  //  //G4int tNH=0;                    // Prototype of the number of secondaries inOut|
  //  try                                                           //                 |
	 //  {                                                             //                 |
		//		  delete output;                                            //                   |
  //    output = pan->Fragment();// DESTROYED in the end of the LOOP work space        |
  //    // @@@@@@@@@@@@@@ Temporary for the testing purposes --- Begin                 |
  //    //tNH=pan->GetNOfHadrons();     // For the test purposes of the String         |
  //    //if(tNH==2)                    // At least 2 hadrons are in the Constr.Output |
		//		  //{//                                                                          |
  //    //  elF=true;                   // Just put a flag for the ellastic Scattering |
	 //    //  delete output;              // Delete a prototype of dummy G4QHadronVector |
  //    //  output = pan->GetHadrons(); // DESTROYED in the end of the LOOP work space |
  //    //}//                                                                          |
  //    //eWei=pan->GetWeight();        // Just an example for the weight of the event |
#ifdef debug
  //    //G4cout<<"=====>>G4QCollision::PostStepDoIt: elF="<<elF<<",n="<<tNH<<G4endl;//|
#endif
  //    // @@@@@@@@@@@@@@ Temporary for the testing purposes --- End                   |
  //  }                                                             //                 |
  //  catch (G4QException& error)//                                                    |
	 //  {                                                             //                 |
	 //    //#ifdef pdebug
  //    G4cerr<<"***G4QCollision::PostStepDoIt: GEN Exception is catched"<<G4endl; //  |
	 //    //#endif
  //    G4Exception("G4QCollision::AtRestDoIt:","27",FatalException,"QString Excep");//|
  //  }                                                             //                 |
  //  delete pan;                              // Delete the Nuclear Environment ---<--+
  //}
  aParticleChange.Initialize(track);
  G4double localtime = track.GetGlobalTime();
  G4ThreeVector position = track.GetPosition();
  // --- the scattered hadron with changed nature can be added here ---
  if(scat)
  {
    G4QHadron* scatHadron = new G4QHadron(scatPDG,scat4M);
    output->push_back(scatHadron);
  }
  // ------------- From here the secondaries are filled -------------------------
  G4int tNH = output->size();       // A#of hadrons in the output
  aParticleChange.SetNumberOfSecondaries(tNH); 
  // Now add nuclear fragments
#ifdef debug
  G4cout<<"G4QCollision::PostStepDoIt: "<<tNH<<" particles are generated"<<G4endl;
#endif
  G4int nOut=output->size();               // Real length of the output @@ Temporary
  if(tNH==1 && !scat)                      // @@ Temporary. Find out why it happened!
  {
    G4cout<<"-Warning-G4QCollision::PostStepDoIt: only one secondary! Make 0."<<G4endl;
    tNH=0;
    delete output->operator[](0);          // delete the creazy hadron
    output->pop_back();                     // clean up the output vector
  }
  if(tNH==2&&2!=nOut) G4cout<<"--Warning--G4QCollision::PostStepDoIt: 2 # "<<nOut<<G4endl;
  // Deal with ParticleChange final state interface to GEANT4 output of the process
  //if(tNH==2) for(i=0; i<tNH; i++) // @@ Temporary tNH==2 instead of just tNH
  if(tNH) for(i=0; i<tNH; i++) // @@ Temporary tNH==2 instead of just tNH
  {
    // Note that one still has to take care of Hypernuclei (with Lambda or Sigma inside)
    // Hypernucleus mass calculation and ion-table interface upgrade => work for Hisaya @@
    // The decau process for hypernuclei must be developed in GEANT4 (change CHIPS body)
    G4QHadron* hadr=output->operator[](i);   // Pointer to the output hadron    
    G4int PDGCode = hadr->GetPDGCode();
    G4int nFrag   = hadr->GetNFragments();
#ifdef pdebug
    G4cout<<"G4QCollision::AtRestDoIt: H#"<<i<<",PDG="<<PDGCode<<",nF="<<nFrag<<G4endl;
#endif
    if(nFrag)                // Skip intermediate (decayed) hadrons
    {
#ifdef debug
	     G4cout<<"G4QCollision::PostStepDoIt: Intermediate particle is found i="<<i<<G4endl;
#endif
      delete hadr;
      continue;
    }
    G4DynamicParticle* theSec = new G4DynamicParticle;  
    G4ParticleDefinition* theDefinition;
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
						G4cout<<"G4QCollision::AtRestDoIt:Ion Z="<<aZ<<", A="<<aA<<G4endl;
#endif
      theDefinition = G4ParticleTable::GetParticleTable()->FindIon(aZ,aA,0,aZ);
    }
    //else theDefinition = G4ParticleTable::GetParticleTable()->FindParticle(PDGCode);
    else
    {
#ifdef pdebug
						G4cout<<"G4QCollision::PostStepDoIt:Define particle with PDG="<<PDGCode<<G4endl;
#endif
      theDefinition = G4QPDGToG4Particle::Get()->GetParticleDefinition(PDGCode);
#ifdef pdebug
						G4cout<<"G4QCollision::PostStepDoIt:AfterParticleDefinition PDG="<<PDGCode<<G4endl;
#endif
    }
    if(!theDefinition)
    {
      G4cout<<"---Warning---G4QCollision::PostStepDoIt: drop PDG="<<PDGCode<<G4endl;
      delete hadr;
      continue;
    }
#ifdef pdebug
    G4cout<<"G4QCollision::PostStepDoIt:Name="<<theDefinition->GetParticleName()<<G4endl;
#endif
    theSec->SetDefinition(theDefinition);
    G4LorentzVector h4M=hadr->Get4Momentum();
    EnMomConservation-=h4M;
#ifdef tdebug
    G4cout<<"G4QCollis::PSDI:"<<i<<","<<PDGCode<<h4M<<h4M.m()<<EnMomConservation<<G4endl;
#endif
#ifdef debug
    G4cout<<"G4QCollision::PostStepDoIt:#"<<i<<",PDG="<<PDGCode<<",4M="<<h4M<<G4endl;
#endif
    theSec->Set4Momentum(h4M); //                                                         ^
    delete hadr; // <-----<-----------<-------------<---------------------<---------<-----+
#ifdef debug
    G4ThreeVector curD=theSec->GetMomentumDirection();              //                    ^
    G4double curM=theSec->GetMass();                                //                    |
    G4double curE=theSec->GetKineticEnergy()+curM;                  //                    ^
    G4cout<<"G4QCollis::PSDoIt:p="<<curD<<curD.mag()<<",e="<<curE<<",m="<<curM<<G4endl;// |
#endif
    G4Track* aNewTrack = new G4Track(theSec, localtime, position ); //                    ^
    aParticleChange.AddSecondary( aNewTrack ); //                                         |
#ifdef debug
    G4cout<<"G4QCollision::PostStepDoIt:#"<<i<<" is done."<<G4endl; //                    |
#endif
  } //                                                                                    |
  delete output; // instances of the G4QHadrons from the output are already deleted above +
  aParticleChange.ProposeTrackStatus(fStopAndKill);        // Kill the absorbed particle
  //return &aParticleChange;                               // This is not enough (ClearILL)
#ifdef debug
    G4cout<<"G4QCollision::PostStepDoIt: **** PostStepDoIt is done ****"<<G4endl;
#endif
  return G4VDiscreteProcess::PostStepDoIt(track, step);
}
