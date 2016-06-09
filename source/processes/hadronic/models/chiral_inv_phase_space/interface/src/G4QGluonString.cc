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
// $Id: G4QGluonString.cc,v 1.1 2006/10/30 10:33:38 mkossov Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//      ---------------- G4QGluonString class -----------------
//                 by Mikhail Kossov, December 2003.
// G4QGluonString class of the CHIPS Simulation Branch in GEANT4
// ---------------------------------------------------------------
// ****************************************************************************************
// ********** This CLASS is temporary moved from the photolepton_hadron directory *********
// ****************************************************************************************

//#define debug
//#define pdebug
//#define ppdebug

#include "G4QGluonString.hh"

// Initialization of static vectors
std::vector<G4int> G4QGluonString::ElementZ;       // Z of the element(i) in theLastCalc
std::vector<G4double> G4QGluonString::ElProbInMat; // SumProbabilityElements in Material
std::vector<std::vector<G4int>*> G4QGluonString::ElIsoN; // N of isotope(j) of Element(i)
std::vector<std::vector<G4double>*>G4QGluonString::IsoProbInEl;//SumProbabIsotopeInElementI

G4QGluonString::G4QGluonString(const G4String& processName):G4VDiscreteProcess(processName)
{
#ifdef debug
  G4cout<<"G4QGluonString::Constructor is called"<<G4endl;
#endif
  if (verboseLevel>0) G4cout<<GetProcessName()<<" process is created by CHIPS"<<G4endl;

  G4QCHIPSWorld::Get()->GetParticles(nPartCWorld); // Create CHIPS World with 234 particles
  G4QNucleus::SetParameters(freeNuc,freeDib,clustProb,mediRatio); // Clusterization param's
  G4Quasmon::SetParameters(Temperature,SSin2Gluons,EtaEtaprime);  // Hadronic parameters
  G4QEnvironment::SetParameters(SolidAngle); // SolAngle of pbar-A secondary mesons capture
  //@@ Initialize here other parameters
}

G4bool   G4QGluonString::manualFlag=false; // If false then standard parameters are used
G4double G4QGluonString::Temperature=180.; // Critical Temperature (sensitive at High En)
G4double G4QGluonString::SSin2Gluons=0.3;  // Supression of s-quarks (in respect to u&d)
G4double G4QGluonString::EtaEtaprime=0.3;  // Supression of eta mesons (gg->qq/3g->qq)
G4double G4QGluonString::freeNuc=0.5;      // Percentage of free nucleons on the surface
G4double G4QGluonString::freeDib=0.05;     // Percentage of free diBaryons on the surface
G4double G4QGluonString::clustProb=5.;     // Nuclear clusterization parameter
G4double G4QGluonString::mediRatio=10.;    // medium/vacuum hadronization ratio
G4int    G4QGluonString::nPartCWorld=152;  // The#of particles initialized in CHIPS World
G4double G4QGluonString::SolidAngle=0.5;   // Part of Solid Angle to capture (@@A-dep.)
G4bool   G4QGluonString::EnergyFlux=false; // Flag for Energy Flux use (not MultyQuasmon)
G4double G4QGluonString::PiPrThresh=141.4; // Pion Production Threshold for gammas
G4double G4QGluonString::M2ShiftVir=20000.;// Shift for M2=-Q2=m_pi^2 of the virtualGamma
G4double G4QGluonString::DiNuclMass=1880.; // DoubleNucleon Mass for VirtualNormalization

void G4QGluonString::SetManual()   {manualFlag=true;}
void G4QGluonString::SetStandard() {manualFlag=false;}

// Fill the private parameters
void G4QGluonString::SetParameters(G4double temper, G4double ssin2g, G4double etaetap,
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

G4QGluonString::~G4QGluonString() {}

// Internal E/P conservation 4-Mom & the selected number of neutrons in the Element

G4LorentzVector G4QGluonString::GetEnegryMomentumConservation()
{
  return EnMomConservation;
}

G4int G4QGluonString::GetNumberOfNeutronsInTarget()
{
  return nOfNeutrons;
}

// output of the function must be in units of length! L=1/sig_V,sig_V=SUM(n(j,i)*sig(j,i)),
// where n(i,j) is a number of nuclei of the isotop j of the element i in V=1(lengtUnit^3)
// ********** All CHIPS cross sections are calculated in the surface units ************
G4double G4QGluonString::GetMeanFreePath(const G4Track& aTrack,
                                         G4double, G4ForceCondition* Fc)
{
#ifdef debug
  G4cout<<"G4QGluonString::GetMeanFreePath: Called Fc="<<*Fc<<G4endl;
#endif
  *Fc = NotForced;
#ifdef debug
  G4cout<<"G4QGluonString::GetMeanFreePath: Before GetDynPart"<<G4endl;
#endif
  const G4DynamicParticle* incidentParticle = aTrack.GetDynamicParticle();
#ifdef debug
  G4cout<<"G4QGluonString::GetMeanFreePath: Before GetDef"<<G4endl;
#endif
  G4ParticleDefinition* incidentParticleDefinition=incidentParticle->GetDefinition();
  if( !IsApplicable(*incidentParticleDefinition))
    G4cout<<"-W-G4QGluonString::GetMeanFreePath called for NotImplementedParticle"<<G4endl;
  // Calculate the mean Cross Section for the set of Elements(*Isotopes) in the Material
  G4double Momentum = incidentParticle->GetTotalMomentum(); // 3-momentum of the Particle
#ifdef debug
  G4cout<<"G4QCollis::GetMeanFreePath: BeforeGetMaterial"<<G4endl;
#endif
  const G4Material* material = aTrack.GetMaterial();        // Get the current material
  const G4double* NOfNucPerVolume = material->GetVecNbOfAtomsPerVolume();
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QGluonString::GetMeanFreePath:"<<nE<<" Elem's in theMaterial"<<G4endl;
#endif
  G4VQCrossSection* CSmanager=0;
  G4int pPDG=0;
  if(incidentParticleDefinition == G4Proton::Proton())
  {
    CSmanager=G4QProtonNuclearCrossSection::GetPointer();
    pPDG=2212;
  }                            //@@ Make cross-section mahnagers for other mesons & baryons
  else
  {
    G4cerr<<"***G4QGluonString::GetMeanFreePath: Particle isn't implemented"<<G4endl;
    G4Exception("G4QGluonString::PostStepDoIt:","72",FatalException,"BadProjectile");
  }
  
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
    G4cout<<"G4QGluonString::GetMeanFreePath: isovectorLength="<<isoSize<<G4endl; // Result
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
          G4cout<<"G4QGluonString::GetMeanFreePath:p#"<<j<<",N="<<N<<",ab="<<abund<<G4endl;
#endif
          newAbund->push_back(pr);
						  }
#ifdef debug
        G4cout<<"G4QGluonString::PostStepDoIt:pairVectorLength="<<newAbund->size()<<G4endl;
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
#ifdef debug
      G4cout<<"GQC::GMF:X="<<CSI<<",M="<<Momentum<<",Z="<<Z<<",N="<<N<<",P="<<pPDG<<G4endl;
#endif
      curIs->second = CSI;                  // Remenber the calculated cross-section
      susi+=CSI;                            // Make a sum per isotopes
      SPI->push_back(susi);                 // Remember summed cross-section
    } // End of temporary initialization of the cross sections in the G4QIsotope singeltone
    sigma+=Isotopes->GetMeanCrossSection(Z,indEl)*NOfNucPerVolume[i];//SUM(MeanCS*NOfNperV)
    ElProbInMat.push_back(sigma);
  } // End of LOOP over Elements
#ifdef debug
  G4cout<<"G4QCol::GetMeanFrPa: S="<<sigma<<",e="<<photNucBias<<",w="<<weakNucBias<<G4endl;
#endif
  // Check that cross section is not zero and return the mean free path
  if(sigma > 0.) return 1./sigma;                 // Mean path [distance] 
  return DBL_MAX;                                 // If Sigma=0, return max value for PATH
}

// Check applicability of the process
G4bool G4QGluonString::IsApplicable(const G4ParticleDefinition& particle) 
{
  if      (particle == *(        G4Proton::Proton()        )) return true;
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
  //else if (particle == *(    G4AntiLambda::AntiLambda()    )) return true;
  //else if (particle == *( G4AntiSigmaPlus::AntiSigmaPlus() )) return true;
  //else if (particle == *(G4AntiSigmaMinus::AntiSigmaMinus())) return true;
  //else if (particle == *( G4AntiSigmaZero::AntiSigmaZero() )) return true;
  //else if (particle == *(   G4AntiXiMinus::AntiXiMinus()   )) return true;
  //else if (particle == *(    G4AntiXiZero::AntiXiZero()    )) return true;
  //else if (particle == *(G4AntiOmegaMinus::AntiOmegaMinus())) return true;
  else G4cerr<<"***G4QGluonString::IsApplicable: PDG?="<<particle.GetPDGEncoding()<<G4endl;
#ifdef debug
  G4cout<<"***G4QGluonString::IsApplicable: PDG="<<particle.GetPDGEncoding()<<G4endl;
#endif
  return false;
}

G4VParticleChange* G4QGluonString::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  //static const G4double dpi=M_PI+M_PI;   // 2*pi (for Phi distr.) ***changed to twopi***
  //static const G4double mNeut= G4QPDGCode(2112).GetMass();
  //static const G4double mNeut2= mNeut*mNeut;          // Squared neutron mass
  //static const G4double mProt= G4QPDGCode(2212).GetMass();
  //static const G4double mProt2= mProt*mProt;          // Squared Proton mass
  //static const G4double dM=mProt+mNeut;               // doubled nucleon mass
  //static const G4double hdM=dM/2.;                    // M of the "nucleon"
  //static const G4double hdM2=hdM*hdM;                 // M2 of the "nucleon"
  //static const G4double mPi0 = G4QPDGCode(111).GetMass();
  //static const G4double mDeut= G4QPDGCode(2112).GetNuclMass(1,1,0);
  //static const G4double mPi  = G4QPDGCode(211).GetMass();
  //static const G4double tmPi = mPi+mPi;               // DoubledMass of the charged pion
  //static const G4double stmPi= tmPi*tmPi;             // SquareDoubledMass of ChargedPion
  //static const G4double mPPi = mPi+mProt;             // Delta threshold
  //static const G4double mPPi2= mPPi*mPPi;             // Delta low threshold for W2
  //-------------------------------------------------------------------------------------
  const G4DynamicParticle* projHadron = track.GetDynamicParticle();
  const G4ParticleDefinition* particle=projHadron->GetDefinition();
#ifdef debug
  G4cout<<"G4QGluonString::PostStepDoIt: Before the GetMeanFreePath is called"<<G4endl;
#endif
  G4ForceCondition cond=NotForced;
  GetMeanFreePath(track, 1., &cond);                  // Just to check that still sig>0
#ifdef debug
  G4cout<<"G4QGluonString::PostStepDoIt: After the GetMeanFreePath is called"<<G4endl;
#endif
  G4LorentzVector proj4M=projHadron->Get4Momentum();
  G4double momentum = projHadron->GetTotalMomentum(); // 3-momentum of the Particle
  G4double Momentum=proj4M.rho();
  if(std::fabs(Momentum-momentum)>.001)
    G4cerr<<"G4QGluonString::PostStepDoIt: P="<<Momentum<<"="<<momentum<<G4endl;
#ifdef debug
  G4double mp=proj4M.m();
  G4cout<<"G4QGluonString::PostStepDoIt is called, P="<<Momentum<<"="<<momentum<<G4endl;
#endif
  if (!IsApplicable(*particle))  // Check applicability
  {
    G4cerr<<"G4QGluonString::PostStepDoIt:Only gam,e+,e-,mu+,mu-,t+,t-,p are implemented."
          <<G4endl;
    return 0;
  }
  const G4Material* material = track.GetMaterial();      // Get the current material
  G4int Z=0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QGluonString::PostStepDoIt: "<<nE<<" elements in the material."<<G4endl;
#endif
  G4int projPDG=0;                           // PDG Code prototype for the captured hadron
  // Not all these particles are implemented yet (see Is Applicable)
  if      (particle ==          G4Proton::Proton()         ) projPDG= 2212;
  //else if (particle ==         G4Neutron::Neutron()        ) projPDG= 2112;
  //else if (particle ==       G4PionMinus::PionMinus()      ) projPDG= -211;
  //else if (particle ==        G4PionPlus::PionPlus()       ) projPDG=  211;
  //else if (particle ==        G4KaonPlus::KaonPlus()       ) projPDG= 2112;
  //else if (particle ==       G4KaonMinus::KaonMinus()      ) projPDG= -321;
  //else if (particle ==    G4KaonZeroLong::KaonZeroLong()   ) projPDG=  130;
  //else if (particle ==   G4KaonZeroShort::KaonZeroShort()  ) projPDG=  310;
  //else if (particle ==         G4TauPlus::TauPlus()        ) projPDG=  -15;
  //else if (particle ==        G4TauMinus::TauMinus()       ) projPDG=   15;
  //else if (particle ==     G4NeutrinoTau::NeutrinoTau()    ) projPDG=   16;
  //else if (particle == G4AntiNeutrinoTau::AntiNeutrinoTau()) projPDG=  -16;
  //else if (particle ==          G4Lambda::Lambda()         ) projPDG= 3122;
  //else if (particle ==       G4SigmaPlus::SigmaPlus()      ) projPDG= 3222;
  //else if (particle ==      G4SigmaMinus::SigmaMinus()     ) projPDG= 3112;
  //else if (particle ==       G4SigmaZero::SigmaZero()      ) projPDG= 3212;
  //else if (particle ==         G4XiMinus::XiMinus()        ) projPDG= 3312;
  //else if (particle ==          G4XiZero::XiZero()         ) projPDG= 3322;
  //else if (particle ==      G4OmegaMinus::OmegaMinus()     ) projPDG= 3334;
  //else if (particle ==     G4AntiNeutron::AntiNeutron()    ) projPDG=-2112;
  //else if (particle ==      G4AntiProton::AntiProton()     ) projPDG=-2212;
  //else if (particle ==      G4AntiLambda::AntiLambda()     ) projPDG=-3122;
  //else if (particle ==   G4AntiSigmaPlus::AntiSigmaPlus()  ) projPDG=-3222;
  //else if (particle ==  G4AntiSigmaMinus::AntiSigmaMinus() ) projPDG=-3112;
  //else if (particle ==   G4AntiSigmaZero::AntiSigmaZero()  ) projPDG=-3212;
  //else if (particle ==     G4AntiXiMinus::AntiXiMinus()    ) projPDG=-3312;
  //else if (particle ==      G4AntiXiZero::AntiXiZero()     ) projPDG=-3322;
  //else if (particle ==  G4AntiOmegaMinus::AntiOmegaMinus() ) projPDG=-3334;
  //G4int aProjPDG=std::abs(projPDG);
#ifdef debug
  G4int prPDG=particle->GetPDGEncoding();
		G4cout<<"G4QGluonString::PostStepDoIt: projPDG="<<projPDG<<", stPDG="<<prPDG<<G4endl;
#endif
  if(!projPDG)
  {
    G4cerr<<"--Warning--G4QGluonString::PostStepDoIt:Undefined interacting hadron"<<G4endl;
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
				  G4cout<<"G4QGluonString::PostStepDoIt:E["<<i<<"]="<<ElProbInMat[i]<<",r="<<rnd<<G4endl;
#endif
      if (rnd<ElProbInMat[i]) break;
    }
    if(i>=nE) i=nE-1;                        // Top limit for the Element
  }
  G4Element* pElement=(*theElementVector)[i];
  Z=static_cast<G4int>(pElement->GetZ());
#ifdef debug
				G4cout<<"G4QGluonString::PostStepDoIt: i="<<i<<", Z(element)="<<Z<<G4endl;
#endif
  if(Z<=0)
  {
    G4cerr<<"---Warning---G4QGluonString::PostStepDoIt: Element with Z="<<Z<<G4endl;
    if(Z<0) return 0;
  }
  std::vector<G4double>* SPI = IsoProbInEl[i];// Vector of summedProbabilities for isotopes
  std::vector<G4int>* IsN = ElIsoN[i];     // Vector of "#of neutrons" in the isotope El[i]
  G4int nofIsot=SPI->size();               // #of isotopes in the element i
#ifdef debug
		G4cout<<"G4QCollis::PosStDoIt:n="<<nofIsot<<",T="<<(*SPI)[nofIsot-1]<<G4endl;
#endif
  G4int j=0;
  if(nofIsot>1)
  {
    G4double rndI=(*SPI)[nofIsot-1]*G4UniformRand(); // Randomize the isotop of the Element
    for(j=0; j<nofIsot; ++j)
    {
#ifdef debug
				  G4cout<<"G4QGluonString::PostStepDoIt: SP["<<j<<"]="<<(*SPI)[j]<<", r="<<rndI<<G4endl;
#endif
      if(rndI < (*SPI)[j]) break;
    }
    if(j>=nofIsot) j=nofIsot-1;            // Top limit for the isotope
  }
  G4int N =(*IsN)[j]; ;                    // Randomized number of neutrons
#ifdef debug
		G4cout<<"G4QGluonString::PostStepDoIt: j="<<i<<", N(isotope)="<<N<<G4endl;
#endif
  if(N<0)
  {
    G4cerr<<"-Warning-G4QGluonString::PostStepDoIt:Isotope with N="<<Z<<"<0,Z="<<Z<<G4endl;
    return 0;
  }
  nOfNeutrons=N;                           // Remember it for the energy-momentum check
  G4double dd=0.025;
  G4double am=Z+N;
  G4double sr=std::sqrt(am);
  G4double dsr=0.01*(sr+sr);
  if(dsr<dd)dsr=dd;
  if(manualFlag) G4QNucleus::SetParameters(freeNuc,freeDib,clustProb,mediRatio);// ManualPa
		//else if(projPDG==-2212) G4QNucleus::SetParameters(1.-dsr-dsr,dd+dd,5.,10.);//aP CluPars
  //else if(projPDG==-211)  G4QNucleus::SetParameters(.67-dsr,.32-dsr,5.,9.); //Pi- CluPars
#ifdef debug
  G4cout<<"G4QGluonString::PostStepDoIt: N="<<N<<" for element with Z="<<Z<<G4endl;
#endif
  if(N<0)
  {
    G4cerr<<"---Warning---G4QGluonString::PostStepDoIt:Element with N="<<N<< G4endl;
    return 0;
  }
  aParticleChange.Initialize(track);
  G4double localtime = track.GetGlobalTime();
  G4ThreeVector position = track.GetPosition();
  G4TouchableHandle trTouchable = track.GetTouchableHandle();
  //
  G4int targPDG=90000000+Z*1000+N;            // PDG Code of the target nucleus
  //    =========================
  G4QPDGCode targQPDG(targPDG);
  G4double tM=targQPDG.GetMass();
  G4QHadronVector* output=new G4QHadronVector;// Prototype of EnvironOutput G4QHadronVector
  G4double absMom = 0.;                       // Prototype of absorbed by nucleus Moment
  G4QHadronVector* leadhs=new G4QHadronVector;// Prototype of QuasmOutput G4QHadronVectorum
  G4LorentzVector lead4M(0.,0.,0.,0.);        // Prototype of LeadingQ 4-momentum
  EnMomConservation=proj4M+G4LorentzVector(0.,0.,0.,tM);    // Total 4-mom of the reaction
  if(absMom) EnMomConservation+=lead4M;       // Add E/M of leading System
#ifdef debug
  G4cout<<"G4QGluonString::PostStepDoIt: projPDG="<<projPDG<<", targPDG="<<targPDG<<G4endl;
#endif
  //G4QHadron* pH = new G4QHadron(projPDG,proj4M);
  // @@@@@@@@@@@@@@@@@@@@@@@@@@ Hrere the CHIPS_QGS must be implemented @@@@@@@@@@@@@@@
  G4int qNH=leadhs->size();
  if(absMom)
  {
    if(qNH) for(G4int iq=0; iq<qNH; iq++)
    {
      G4QHadron* loh=(*leadhs)[iq];   // Pointer to the output hadron
      output->push_back(loh);
    }
    delete leadhs;
  }
  // ------------- From here the secondaries are filled -------------------------
  G4int tNH = output->size();       // A#of hadrons in the output
  aParticleChange.SetNumberOfSecondaries(tNH); 
  // Now add nuclear fragments
#ifdef debug
  G4cout<<"G4QGluonString::PostStepDoIt: "<<tNH<<" particles are generated"<<G4endl;
#endif
#ifdef ppdebug
  if(absMom)G4cout<<"G4QGluonString::PostStepDoIt: t="<<tNH<<", q="<<qNH<<G4endl;
#endif
  G4int nOut=output->size();               // Real length of the output @@ Temporary
  if(tNH==1)                               // @@ Temporary. Find out why it happened!
  {
    G4cout<<"-Warning-G4QGluonString::PostStepDoIt: 1 secondary! absMom="<<absMom;
    if(absMom) G4cout<<", qNH="<<qNH;
    G4cout<<", PDG0="<<(*output)[0]->GetPDGCode();
    G4cout<<G4endl;
    tNH=0;
    delete output->operator[](0);          // delete the creazy hadron
    output->pop_back();                    // clean up the output vector
  }
  if(tNH==2&&2!=nOut) G4cout<<"--Warning--G4QGluonString::PostStepDoIt:2 # "<<nOut<<G4endl;
  // Deal with ParticleChange final state interface to GEANT4 output of the process
  //if(tNH==2) for(i=0; i<tNH; i++) // @@ Temporary tNH==2 instead of just tNH
  if(tNH) for(i=0; i<tNH; i++) // @@ Temporary tNH==2 instead of just tNH
  {
    // Note that one still has to take care of Hypernuclei (with Lambda or Sigma inside)
    // Hypernucleus mass calculation and ion-table interface upgrade => work for Hisaya @@
    // The decau process for hypernuclei must be developed in GEANT4 (change CHIPS body)
    G4QHadron* hadr=(*output)[i];          // Pointer to the output hadron    
    G4int PDGCode = hadr->GetPDGCode();
    G4int nFrag   = hadr->GetNFragments();
#ifdef pdebug
    G4cout<<"G4QGluonString::AtRestDoIt: H#"<<i<<",PDG="<<PDGCode<<",nF="<<nFrag<<G4endl;
#endif
    if(nFrag)                // Skip intermediate (decayed) hadrons
    {
#ifdef debug
	     G4cout<<"G4QGluonString::PostStepDoIt: Intermediate particle is found i="<<i<<G4endl;
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
						G4cout<<"G4QGluonString::AtRestDoIt:Ion Z="<<aZ<<", A="<<aA<<G4endl;
#endif
      theDefinition = G4ParticleTable::GetParticleTable()->FindIon(aZ,aA,0,aZ);
    }
    //else theDefinition = G4ParticleTable::GetParticleTable()->FindParticle(PDGCode);
    else
    {
#ifdef pdebug
						G4cout<<"G4QGluonString::PostStepDoIt:Define particle with PDG="<<PDGCode<<G4endl;
#endif
      theDefinition = G4QPDGToG4Particle::Get()->GetParticleDefinition(PDGCode);
#ifdef pdebug
						G4cout<<"G4QGluonString::PostStepDoIt:AfterParticleDefinition PDG="<<PDGCode<<G4endl;
#endif
    }
    if(!theDefinition)
    {
#ifdef debug
      G4cout<<"---Warning---G4QGluonString::PostStepDoIt: drop PDG="<<PDGCode<<G4endl;
#endif
      delete hadr;
      continue;
    }
#ifdef pdebug
    G4cout<<"G4QGluonString::PostStepDoIt:Name="<<theDefinition->GetParticleName()<<G4endl;
#endif
    theSec->SetDefinition(theDefinition);
    G4LorentzVector h4M=hadr->Get4Momentum();
    EnMomConservation-=h4M;
#ifdef tdebug
    G4cout<<"G4QCollis::PSDI:"<<i<<","<<PDGCode<<h4M<<h4M.m()<<EnMomConservation<<G4endl;
#endif
#ifdef debug
    G4cout<<"G4QGluonString::PostStepDoIt:#"<<i<<",PDG="<<PDGCode<<",4M="<<h4M<<G4endl;
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
    aNewTrack->SetTouchableHandle(trTouchable);                     //                    |
    aParticleChange.AddSecondary( aNewTrack );                      //                    |
#ifdef debug
    G4cout<<"G4QGluonString::PostStepDoIt:#"<<i<<" is done"<<G4endl;//                    |
#endif
  } //                                                                                    |
  delete output; // instances of the G4QHadrons from the output are already deleted above +
#ifdef debug
		G4cout<<"G4QGluonString::PostStDoIt: afterSt="<<aParticleChange.GetTrackStatus()<<G4endl;
#endif
  aParticleChange.ProposeTrackStatus(fStopAndKill);        // Kill the absorbed particle
#ifdef debug
		G4cout<<"G4QGluonString::PostStepDoIt:*** PostStepDoIt is done ***, P="<<aProjPDG
        <<", St="<<aParticleChange.GetTrackStatus()<<G4endl;
#endif
  return G4VDiscreteProcess::PostStepDoIt(track, step);
}
