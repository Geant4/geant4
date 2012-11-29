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
// $Id$
//
//      ---------------- G4QInelastic class -----------------
//                 by Mikhail Kossov, December 2003.
// G4QInelastic class of the CHIPS Simulation Branch in GEANT4
// ---------------------------------------------------------------
// ****************************************************************************************
// This Header is a part of the CHIPS physics package (author: M. Kosov)
// ****************************************************************************************
// Short description: This is a universal class for the incoherent (inelastic)
// nuclear interactions in the CHIPS model.
// ---------------------------------------------------------------------------
//#define debug
//#define pdebug
//#define pickupd
//#define ldebug
//#define ppdebug
//#define qedebug

#include "G4QInelastic.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicDeprecate.hh"


// Initialization of static vectors
std::vector<G4int> G4QInelastic::ElementZ;            // Z of the element(i) in theLastCalc
std::vector<G4double> G4QInelastic::ElProbInMat;      // SumProbabilityElements in Material
std::vector<std::vector<G4int>*> G4QInelastic::ElIsoN;// N of isotope(j) of Element(i)
std::vector<std::vector<G4double>*>G4QInelastic::IsoProbInEl;//SumProbabIsotopes inElementI

G4QInelastic::G4QInelastic(const G4String& processName): 
 G4VDiscreteProcess(processName, fHadronic)
{
  G4HadronicDeprecate("G4QInelastic");

  EnMomConservation = G4LorentzVector(0.,0.,0.,0.);
  nOfNeutrons       = 0;
#ifdef debug
  G4cout<<"G4QInelastic::Constructor is called"<<G4endl;
#endif
  if (verboseLevel>0) G4cout << GetProcessName() << " process is created "<< G4endl;
  //G4QCHIPSWorld::Get()->GetParticles(nPartCWorld); // Create CHIPSWorld (234 part.max)
  G4QNucleus::SetParameters(freeNuc,freeDib,clustProb,mediRatio); // Clusterization param's
  G4Quasmon::SetParameters(Temperature,SSin2Gluons,EtaEtaprime);  // Hadronic parameters
  G4QEnvironment::SetParameters(SolidAngle); // SolAngle of pbar-A secondary mesons capture
  //@@ Initialize here the G4QuasmonString parameters
}

G4bool   G4QInelastic::manualFlag=false; // If false then standard parameters are used
G4double G4QInelastic::Temperature=140.; // Critical Temperature (sensitive at High En)
G4double G4QInelastic::SSin2Gluons=0.3;  // Supression of s-quarks (in respect to u&d)
G4double G4QInelastic::EtaEtaprime=0.3;  // Supression of eta mesons (gg->qq/3g->qq)
G4double G4QInelastic::freeNuc=.04;       // Percentage of free nucleons on the surface
G4double G4QInelastic::freeDib=.08;      // Percentage of free diBaryons on the surface
G4double G4QInelastic::clustProb=4.;     // Nuclear clusterization parameter
G4double G4QInelastic::mediRatio=1.;     // medium/vacuum hadronization ratio
//G4int    G4QInelastic::nPartCWorld=152;// The#of particles initialized in CHIPS World
//G4int    G4QInelastic::nPartCWorld=122;// The#of particles initialized in CHIPS World
G4int    G4QInelastic::nPartCWorld=85;   // The#of particles initialized in CHIPS World Red
G4double G4QInelastic::SolidAngle=0.5;   // Part of Solid Angle to capture (@@A-dep.)
G4bool   G4QInelastic::EnergyFlux=false; // Flag for Energy Flux use (not MultyQuasmon)
G4double G4QInelastic::PiPrThresh=141.4; // Pion Production Threshold for gammas
G4double G4QInelastic::M2ShiftVir=20000.;// Shift for M2=-Q2=m_pi^2 of the virtualGamma
G4double G4QInelastic::DiNuclMass=1870.; // DoubleNucleon Mass for VirtualNormalization
G4double G4QInelastic::photNucBias=1.;   // BiasingParameter for photo(e,mu,tau)Nuclear
G4double G4QInelastic::weakNucBias=1.;   // BiasingParameter for ChargedCurrents(nu,mu) 

void G4QInelastic::SetManual()   {manualFlag=true;}
void G4QInelastic::SetStandard() {manualFlag=false;}

// Fill the private parameters
void G4QInelastic::SetParameters(G4double temper, G4double ssin2g, G4double etaetap,
                                     G4double fN, G4double fD, G4double cP, G4double mR,
                                     G4int nParCW, G4double solAn, G4bool efFlag,
                                     G4double piThresh, G4double mpisq, G4double dinum)
{
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

void G4QInelastic::SetPhotNucBias(G4double phnB) {photNucBias=phnB;}
void G4QInelastic::SetWeakNucBias(G4double ccnB) {weakNucBias=ccnB;}

// Destructor

G4QInelastic::~G4QInelastic()
{
  // The following is just a copy of what is done in PostStepDoIt every interaction !
  // The correction is if(IPIE), so just for(...;ip<IPIE;...) does not work ! @@
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
}


G4LorentzVector G4QInelastic::GetEnegryMomentumConservation()
{
  return EnMomConservation;
}

G4int G4QInelastic::GetNumberOfNeutronsInTarget()
{
  return nOfNeutrons;
}

// output of the function must be in units of length! L=1/sig_V,sig_V=SUM(n(j,i)*sig(j,i)),
// where n(i,j) is a number of nuclei of the isotop j of the element i in V=1(lengtUnit^3)
// ********** All CHIPS cross sections are calculated in the surface units ************
G4double G4QInelastic::GetMeanFreePath(const G4Track& aTrack,G4double,G4ForceCondition* Fc)
{
#ifdef debug
  G4cout<<"G4QInelastic::GetMeanFreePath: Called Fc="<<*Fc<<G4endl;
#endif
  *Fc = NotForced;
#ifdef debug
  G4cout<<"G4QInelastic::GetMeanFreePath: Before GetDynPart"<<G4endl;
#endif
  const G4DynamicParticle* incidentParticle = aTrack.GetDynamicParticle();
#ifdef debug
  G4cout<<"G4QInelastic::GetMeanFreePath: Before GetDef"<<G4endl;
#endif
  G4ParticleDefinition* incidentParticleDefinition=incidentParticle->GetDefinition();
  if( !IsApplicable(*incidentParticleDefinition))
  {
    G4cout<<"-W-G4QInelastic::GetMeanFreePath called for not implemented particle"<<G4endl;
    return DBL_MAX;
  }
  // Calculate the mean Cross Section for the set of Elements(*Isotopes) in the Material
  G4double Momentum = incidentParticle->GetTotalMomentum(); // 3-momentum of the Particle
#ifdef debug
  G4cout<<"G4QInelastic::GetMeanFreePath: BeforeGetMaterial"<<G4endl;
#endif
  const G4Material* material = aTrack.GetMaterial();        // Get the current material
  const G4double* NOfNucPerVolume = material->GetVecNbOfAtomsPerVolume();
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QInelastic::GetMeanFreePath:"<<nE<<" Elem's in theMaterial"<<G4endl;
#endif
  //G4bool leptoNuc=false; // By default the reaction is not lepto-nuclear *Growing point*
  G4VQCrossSection* CSmanager=0;
  G4VQCrossSection* CSmanager2=0;
  G4QNeutronCaptureRatio* capMan=0;
  G4int pPDG=0;
  G4int pZ = incidentParticleDefinition->GetAtomicNumber();
  G4int pA = incidentParticleDefinition->GetAtomicMass();
  if(incidentParticleDefinition == G4Neutron::Neutron())
  {
    CSmanager=G4QNeutronNuclearCrossSection::GetPointer();
    capMan=G4QNeutronCaptureRatio::GetPointer();
#ifdef debug
    G4cout<<"G4QInelastic::GetMeanFreePath: CSmanager is defined for the neutron"<<G4endl;
#endif
    pPDG=2112;
  }
  else if(incidentParticleDefinition == G4Proton::Proton())
  {
    CSmanager=G4QProtonNuclearCrossSection::GetPointer();
    pPDG=2212;
  }
  else if(incidentParticleDefinition == G4PionMinus::PionMinus())
  {
    CSmanager=G4QPionMinusNuclearCrossSection::GetPointer();
    pPDG=-211;
  }
  else if(incidentParticleDefinition == G4PionPlus::PionPlus())
  {
    CSmanager=G4QPionPlusNuclearCrossSection::GetPointer();
    pPDG=211;
  }
  else if(incidentParticleDefinition == G4KaonMinus::KaonMinus())
  {
    CSmanager=G4QKaonMinusNuclearCrossSection::GetPointer();
    pPDG=-321;
  }
  else if(incidentParticleDefinition == G4KaonPlus::KaonPlus())
  {
    CSmanager=G4QKaonPlusNuclearCrossSection::GetPointer();
    pPDG=321;
  }
  else if(incidentParticleDefinition == G4KaonZeroLong::KaonZeroLong()   ||
          incidentParticleDefinition == G4KaonZeroShort::KaonZeroShort() ||
          incidentParticleDefinition == G4KaonZero::KaonZero()           ||
          incidentParticleDefinition == G4AntiKaonZero::AntiKaonZero()   )
  {
    CSmanager=G4QKaonZeroNuclearCrossSection::GetPointer();
    if(G4UniformRand() > 0.5) pPDG= 311;
    else                      pPDG=-311;
  }
  else if(incidentParticleDefinition == G4Lambda::Lambda())
  {
    CSmanager=G4QHyperonNuclearCrossSection::GetPointer();
    pPDG=3122;
  }
  else if(pZ > 0 && pA > 1) // Ions (not implemented yet, should not be used)
  {
    G4cout<<"-Warning-G4QInelastic::GetMeanFreePath: G4QInelastic called for ions"<<G4endl;
    CSmanager=G4QProtonNuclearCrossSection::GetPointer();
    pPDG=90000000+999*pZ+pA;
  }
  else if(incidentParticleDefinition == G4SigmaPlus::SigmaPlus())
  {
    CSmanager=G4QHyperonPlusNuclearCrossSection::GetPointer();
    pPDG=3222;
  }
  else if(incidentParticleDefinition == G4SigmaMinus::SigmaMinus())
  {
    CSmanager=G4QHyperonNuclearCrossSection::GetPointer();
    pPDG=3112;
  }
  else if(incidentParticleDefinition == G4SigmaZero::SigmaZero())
  {
    CSmanager=G4QHyperonNuclearCrossSection::GetPointer();
    pPDG=3212;
  }
  else if(incidentParticleDefinition == G4XiMinus::XiMinus())
  {
    CSmanager=G4QHyperonNuclearCrossSection::GetPointer();
    pPDG=3312;
  }
  else if(incidentParticleDefinition == G4XiZero::XiZero())
  {
    CSmanager=G4QHyperonNuclearCrossSection::GetPointer();
    pPDG=3322;
  }
  else if(incidentParticleDefinition == G4OmegaMinus::OmegaMinus())
  {
    CSmanager=G4QHyperonNuclearCrossSection::GetPointer();
    pPDG=3334;
  }
  else if(incidentParticleDefinition == G4MuonPlus::MuonPlus() ||
          incidentParticleDefinition == G4MuonMinus::MuonMinus())
  {
    CSmanager=G4QMuonNuclearCrossSection::GetPointer();
    //leptoNuc=true;
    pPDG=13;
  }
  else if(incidentParticleDefinition == G4AntiNeutron::AntiNeutron())
  {
    CSmanager=G4QAntiBaryonNuclearCrossSection::GetPointer();
    pPDG=-2112;
  }
  else if(incidentParticleDefinition == G4AntiProton::AntiProton())
  {
    CSmanager=G4QAntiBaryonNuclearCrossSection::GetPointer();
    pPDG=-2212;
  }
  else if(incidentParticleDefinition == G4AntiLambda::AntiLambda())
  {
    CSmanager=G4QAntiBaryonNuclearCrossSection::GetPointer();
    pPDG=-3122;
  }
  else if(incidentParticleDefinition == G4AntiSigmaPlus::AntiSigmaPlus())
  {
    CSmanager=G4QAntiBaryonNuclearCrossSection::GetPointer();
    pPDG=-3222;
  }
  else if(incidentParticleDefinition == G4AntiSigmaMinus::AntiSigmaMinus())
  {
    CSmanager=G4QAntiBaryonPlusNuclearCrossSection::GetPointer();
    pPDG=-3112;
  }
  else if(incidentParticleDefinition == G4AntiSigmaZero::AntiSigmaZero())
  {
    CSmanager=G4QAntiBaryonNuclearCrossSection::GetPointer();
    pPDG=-3212;
  }
  else if(incidentParticleDefinition == G4AntiXiMinus::AntiXiMinus())
  {
    CSmanager=G4QAntiBaryonPlusNuclearCrossSection::GetPointer();
    pPDG=-3312;
  }
  else if(incidentParticleDefinition == G4AntiXiZero::AntiXiZero())
  {
    CSmanager=G4QAntiBaryonNuclearCrossSection::GetPointer();
    pPDG=-3322;
  }
  else if(incidentParticleDefinition == G4AntiOmegaMinus::AntiOmegaMinus())
  {
    CSmanager=G4QAntiBaryonPlusNuclearCrossSection::GetPointer();
    pPDG=-3334;
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
    //leptoNuc=true;
    pPDG=11;
  }
  else if(incidentParticleDefinition == G4TauPlus::TauPlus() ||
          incidentParticleDefinition == G4TauMinus::TauMinus())
  {
    CSmanager=G4QTauNuclearCrossSection::GetPointer();
    //leptoNuc=true;
    pPDG=15;
  }
  else if(incidentParticleDefinition == G4NeutrinoMu::NeutrinoMu() )
  {
    CSmanager=G4QNuMuNuclearCrossSection::GetPointer();
    CSmanager2=G4QNuNuNuclearCrossSection::GetPointer();
    //leptoNuc=true;
    pPDG=14;
  }
  else if(incidentParticleDefinition == G4AntiNeutrinoMu::AntiNeutrinoMu() )
  {
    CSmanager=G4QANuMuNuclearCrossSection::GetPointer();
    CSmanager2=G4QANuANuNuclearCrossSection::GetPointer();
    //leptoNuc=true;
    pPDG=-14;
  }
  else if(incidentParticleDefinition == G4NeutrinoE::NeutrinoE() )
  {
    CSmanager=G4QNuENuclearCrossSection::GetPointer();
    CSmanager2=G4QNuNuNuclearCrossSection::GetPointer();
    //leptoNuc=true;
    pPDG=12;
  }
  else if(incidentParticleDefinition == G4AntiNeutrinoE::AntiNeutrinoE() )
  {
    CSmanager=G4QANuENuclearCrossSection::GetPointer();
    CSmanager2=G4QANuANuNuclearCrossSection::GetPointer();
    //leptoNuc=true;
    pPDG=-12;
  }
  else
  {
    G4cout<<"-Warning-G4QInelastic::GetMeanFreePath:Particle "
          <<incidentParticleDefinition->GetPDGEncoding()<<" isn't supported"<<G4endl;
    return DBL_MAX;
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
    G4cout<<"G4QInelastic::GetMeanFreePath: isovectorLength="<<isoSize<<G4endl; // Result
#endif
    if(isoSize)                             // The Element has non-trivial abundance set
    {
      indEl=pElement->GetIndex()+1;         // Index of the non-trivial element
      if(!Isotopes->IsDefined(Z,indEl))     // This index is not defined for this Z: define
      {
        std::vector<std::pair<G4int,G4double>*>* newAbund =
                                               new std::vector<std::pair<G4int,G4double>*>;
        G4double* abuVector=pElement->GetRelativeAbundanceVector();
        for(G4int j=0; j<isoSize; j++)      // Calculation of abundance vector for isotopes
        {
          G4int N=pElement->GetIsotope(j)->GetN()-Z; // N means A=N+Z !
          if(pElement->GetIsotope(j)->GetZ()!=Z)G4cerr<<"G4QInelastic::GetMeanFreePath"
                                 <<": Z="<<pElement->GetIsotope(j)->GetZ()<<"#"<<Z<<G4endl;
          G4double abund=abuVector[j];
          std::pair<G4int,G4double>* pr= new std::pair<G4int,G4double>(N,abund);
#ifdef debug
          G4cout<<"G4QInelastic::GetMeanFreePath: p#="<<j<<",N="<<N<<",ab="<<abund<<G4endl;
#endif
          newAbund->push_back(pr);
        }
#ifdef debug
        G4cout<<"G4QInelastic::GetMeanFreePath: pairVectLength="<<newAbund->size()<<G4endl;
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
#ifdef debug
    G4cout<<"G4QInelastic::GetMeanFreePath: Before Loop nIs="<<nIs<<G4endl;
#endif
    if(nIs) for(G4int j=0; j<nIs; j++)      // Calculate CS for eachIsotope of El
    {
      std::pair<G4int,G4double>* curIs=(*cs)[j]; // A pointer, which is used twice
      G4int N=curIs->first;                 // #of Neuterons in the isotope j of El i
      IsN->push_back(N);                    // Remember Min N for the Element
#ifdef debug
      G4cout<<"G4QInel::GetMeanFrP: Before CS, P="<<Momentum<<",Z="<<Z<<",N="<<N<<G4endl;
#endif
      if(!pPDG) G4cout<<"-Warning-G4QInelastic::GetMeanFrP: (1) projectile PDG=0"<<G4endl;
      G4double CSI=CSmanager->GetCrossSection(true,Momentum,Z,N,pPDG);//CS(j,i) for isotope
      if(CSmanager2)CSI+=CSmanager2->GetCrossSection(true,Momentum,Z,N,pPDG);//CS(j,i)nu,nu
      if(capMan) CSI*=(1.-capMan->GetRatio(Momentum, Z, N));
#ifdef debug
      G4cout<<"GQC::GMF:X="<<CSI<<",M="<<Momentum<<",Z="<<Z<<",N="<<N<<",P="<<pPDG<<G4endl;
#endif
      curIs->second = CSI;
      susi+=CSI;                            // Make a sum per isotopes
      SPI->push_back(susi);                 // Remember summed cross-section
    } // End of temporary initialization of the cross sections in the G4QIsotope singeltone
    sigma+=Isotopes->GetMeanCrossSection(Z,indEl)*NOfNucPerVolume[i];//SUM(MeanCS*NOfNperV)
    ElProbInMat.push_back(sigma);
  } // End of LOOP over Elements
#ifdef debug
  G4cout<<"G4QInel::GetMeanFrPa:S="<<sigma<<",e="<<photNucBias<<",w="<<weakNucBias<<G4endl;
#endif
  // Check that cross section is not zero and return the mean free path
  if(photNucBias!=1.) if(incidentParticleDefinition == G4Gamma::Gamma()         ||
                         incidentParticleDefinition == G4MuonPlus::MuonPlus()   ||
                         incidentParticleDefinition == G4MuonMinus::MuonMinus() ||
                         incidentParticleDefinition == G4Electron::Electron()   ||
                         incidentParticleDefinition == G4Positron::Positron()   ||  
                         incidentParticleDefinition == G4TauMinus::TauMinus()   ||
                         incidentParticleDefinition == G4TauPlus::TauPlus()       )
                                                                        sigma*=photNucBias;
  if(weakNucBias!=1.) if(incidentParticleDefinition==G4NeutrinoE::NeutrinoE()            ||
                         incidentParticleDefinition==G4AntiNeutrinoE::AntiNeutrinoE()    ||
                         incidentParticleDefinition==G4NeutrinoTau::NeutrinoTau()        ||
                         incidentParticleDefinition==G4AntiNeutrinoTau::AntiNeutrinoTau()||
                         incidentParticleDefinition==G4NeutrinoMu::NeutrinoMu()          ||
                         incidentParticleDefinition==G4AntiNeutrinoMu::AntiNeutrinoMu()   )
                                                                        sigma*=weakNucBias;
  if(sigma > 0.) return 1./sigma;                 // Mean path [distance] 
  return DBL_MAX;
}


G4bool G4QInelastic::IsApplicable(const G4ParticleDefinition& particle) 
{
  G4int Z=particle.GetAtomicNumber();
  G4int A=particle.GetAtomicMass();
  if      (particle == *(       G4Neutron::Neutron()       )) return true; 
  else if (particle == *(        G4Proton::Proton()        )) return true;
  else if (particle == *(      G4MuonPlus::MuonPlus()      )) return true;
  else if (particle == *(     G4MuonMinus::MuonMinus()     )) return true;
  else if (particle == *(         G4Gamma::Gamma()         )) return true;
  else if (particle == *(      G4Electron::Electron()      )) return true;
  else if (particle == *(      G4Positron::Positron()      )) return true;
  else if (particle == *(     G4PionMinus::PionMinus()     )) return true;
  else if (particle == *(      G4PionPlus::PionPlus()      )) return true;
  else if (particle == *(      G4KaonPlus::KaonPlus()      )) return true;
  else if (particle == *(     G4KaonMinus::KaonMinus()     )) return true;
  else if (particle == *(  G4KaonZeroLong::KaonZeroLong()  )) return true;
  else if (particle == *( G4KaonZeroShort::KaonZeroShort() )) return true;
  else if (particle == *(        G4Lambda::Lambda()        )) return true;
  else if (particle == *(     G4SigmaPlus::SigmaPlus()     )) return true;
  else if (particle == *(    G4SigmaMinus::SigmaMinus()    )) return true;
  else if (particle == *(     G4SigmaZero::SigmaZero()     )) return true;
  else if (particle == *(       G4XiMinus::XiMinus()       )) return true;
  else if (particle == *(        G4XiZero::XiZero()        )) return true;
  else if (particle == *(    G4OmegaMinus::OmegaMinus()    )) return true;
  else if (particle == *(G4GenericIon::GenericIon()) || (Z > 0 && A > 1)) return true;
  else if (particle == *(   G4AntiNeutron::AntiNeutron()   )) return true;
  else if (particle == *(    G4AntiProton::AntiProton()    )) return true;
  else if (particle == *(    G4AntiLambda::AntiLambda()    )) return true;
  else if (particle == *( G4AntiSigmaPlus::AntiSigmaPlus() )) return true;
  else if (particle == *(G4AntiSigmaMinus::AntiSigmaMinus())) return true;
  else if (particle == *( G4AntiSigmaZero::AntiSigmaZero() )) return true;
  else if (particle == *(   G4AntiXiMinus::AntiXiMinus()   )) return true;
  else if (particle == *(    G4AntiXiZero::AntiXiZero()    )) return true;
  else if (particle == *(G4AntiOmegaMinus::AntiOmegaMinus())) return true;
  else if (particle == *(       G4TauPlus::TauPlus()       )) return true;
  else if (particle == *(      G4TauMinus::TauMinus()      )) return true;
  else if (particle == *( G4AntiNeutrinoE::AntiNeutrinoE() )) return true;
  else if (particle == *(     G4NeutrinoE::NeutrinoE()     )) return true;
  else if (particle == *(G4AntiNeutrinoMu::AntiNeutrinoMu())) return true;
  else if (particle == *(    G4NeutrinoMu::NeutrinoMu()    )) return true;
  //else if (particle == *(G4AntiNeutrinoTau::AntiNeutrinoTau())) return true;
  //else if (particle == *(    G4NeutrinoTau::NeutrinoTau()    )) return true;
#ifdef debug
  G4cout<<"***G4QInelastic::IsApplicable: PDG="<<particle.GetPDGEncoding()<<G4endl;
#endif
  return false;
}

G4VParticleChange* G4QInelastic::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  static const G4double third = 1./3.;
  static const G4double me=G4Electron::Electron()->GetPDGMass();   // electron mass
  static const G4double me2=me*me;                                 // squared electron mass
  static const G4double mu=G4MuonMinus::MuonMinus()->GetPDGMass(); // muon mass
  static const G4double mu2=mu*mu;                                 // squared muon mass
  static const G4double mt=G4TauMinus::TauMinus()->GetPDGMass();   // tau mass
  static const G4double mt2=mt*mt;                                 // squared tau mass
  //static const G4double dpi=M_PI+M_PI;   // 2*pi (for Phi distr.) ***changed to twopi***
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
  static const G4double mNeut2= mNeut*mNeut;
  static const G4double mProt= G4QPDGCode(2212).GetMass();
  static const G4double mProt2= mProt*mProt;
  static const G4double dM=mProt+mNeut;                            // doubled nucleon mass
  static const G4double hdM=dM/2.;                                 // M of the "nucleon"
  static const G4double hdM2=hdM*hdM;                              // M2 of the "nucleon"
  static const G4double mLamb= G4QPDGCode(3122).GetMass();         // Mass of Lambda/antiL
  static const G4double mPi0 = G4QPDGCode(111).GetMass();
  static const G4double mPi0s= mPi0*mPi0;
  static const G4double mDeut= G4QPDGCode(2112).GetNuclMass(1,1,0);// Mass of deuteron
  //static const G4double mDeut2=mDeut*mDeut;                      // SquaredMassOfDeuteron
  static const G4double mTrit= G4QPDGCode(2112).GetNuclMass(1,2,0);// Mass of tritium
  static const G4double mHel3= G4QPDGCode(2112).GetNuclMass(2,1,0);// Mass of Helium3
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0);// Mass of alpha
  static const G4double mPi  = G4QPDGCode(211).GetMass();
  static const G4double tmPi = mPi+mPi;     // Doubled mass of the charged pion
  static const G4double stmPi= tmPi*tmPi;   // Squared Doubled mass of the charged pion
  static const G4double mPPi = mPi+mProt;   // Delta threshold
  static const G4double mPPi2= mPPi*mPPi; // Delta low threshold for W2
  //static const G4double mDel2= 1400*1400; // Delta up threshold for W2 (in MeV^2)
  // Static definitions for electrons (nu,e) -----------------------------------------
  static const G4double meN = mNeut+me;
  static const G4double meN2= meN*meN;
  static const G4double fmeN= 4*mNeut2*me2;
  static const G4double mesN= mNeut2+me2;
  static const G4double meP = mProt+me;
  static const G4double meP2= meP*meP;
  static const G4double fmeP= 4*mProt2*me2;
  static const G4double mesP= mProt2+me2;
  static const G4double medM= me2/dM;       // for x limit
  static const G4double meD = mPPi+me;      // Multiperipheral threshold
  static const G4double meD2= meD*meD;
  // Static definitions for muons (nu,mu) -----------------------------------------
  static const G4double muN = mNeut+mu;
  static const G4double muN2= muN*muN;
  static const G4double fmuN= 4*mNeut2*mu2;
  static const G4double musN= mNeut2+mu2;
  static const G4double muP = mProt+mu;
  static const G4double muP2= muP*muP;      // +
  static const G4double fmuP= 4*mProt2*mu2; // +
  static const G4double musP= mProt2+mu2;
  static const G4double mudM= mu2/dM;       // for x limit
  static const G4double muD = mPPi+mu;      // Multiperipheral threshold
  static const G4double muD2= muD*muD;
  // Static definitions for muons (nu,nu) -----------------------------------------
  //static const G4double nuN = mNeut;
  //static const G4double nuN2= mNeut2;
  //static const G4double fnuN= 0.;
  //static const G4double nusN= mNeut2;
  //static const G4double nuP = mProt;
  //static const G4double nuP2= mProt2;
  //static const G4double fnuP= 0.;
  //static const G4double nusP= mProt2;
  //static const G4double nudM= 0.;           // for x limit
  //static const G4double nuD = mPPi;         // Multiperipheral threshold
  //static const G4double nuD2= mPPi2;
  static const G4LorentzVector vacuum4M(0.,0.,0.,0.);
  //-------------------------------------------------------------------------------------
  static G4bool CWinit = true;              // CHIPS Warld needs to be initted
  if(CWinit)
  {
    CWinit=false;
    G4QCHIPSWorld::Get()->GetParticles(nPartCWorld); // Create CHIPS World (234 part.max)
  }
  //-------------------------------------------------------------------------------------
  const G4DynamicParticle* projHadron = track.GetDynamicParticle();
  const G4ParticleDefinition* particle=projHadron->GetDefinition();
#ifdef debug
  G4cout<<"G4QInelastic::PostStepDoIt: Before the GetMeanFreePath is called"<<G4endl;
#endif
  G4ForceCondition cond=NotForced;
  GetMeanFreePath(track, 1., &cond);
#ifdef debug
  G4cout<<"G4QInelastic::PostStepDoIt: After the GetMeanFreePath is called"<<G4endl;
#endif
  G4bool scat=false;                                  // No CHEX in proj scattering
  G4int  scatPDG=0;                                   // Must be filled if true (CHEX)
  G4LorentzVector proj4M=projHadron->Get4Momentum();  // 4-momentum of the projectile (IU?)
  G4LorentzVector scat4M=proj4M;                      // Must be filled if true
  G4double momentum = projHadron->GetTotalMomentum(); // 3-momentum of the Particle
  G4double Momentum=proj4M.rho();
  if(std::fabs(Momentum-momentum)>.001)
    G4cerr<<"*G4QInelastic::PostStepDoIt: P="<<Momentum<<"#"<<momentum<<G4endl;
#ifdef debug
  G4double mp=proj4M.m();
  G4cout<<"-->G4QInel::PostStDoIt:*called*, 4M="<<proj4M<<", P="<<Momentum<<"="<<momentum
        <<",m="<<mp<<G4endl;
#endif
  if (!IsApplicable(*particle))  // Check applicability
  {
    G4cerr<<"G4QInelastic::PostStepDoIt:Only gam,e+,e-,mu+,mu-,t+,t-,p are implemented."
          <<G4endl;
    return 0;
  }
  const G4Material* material = track.GetMaterial();      // Get the current material
  G4int Z=0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QInelastic::PostStepDoIt: "<<nE<<" elements in the material."<<G4endl;
#endif
  G4int projPDG=0;                           // PDG Code prototype for the captured hadron
  G4int pZ=particle->GetAtomicNumber();
  G4int pA=particle->GetAtomicMass();
  if      (particle ==         G4Neutron::Neutron()        ) projPDG= 2112;
  else if (particle ==          G4Proton::Proton()         ) projPDG= 2212;
  else if (particle ==       G4PionMinus::PionMinus()      ) projPDG= -211;
  else if (particle ==        G4PionPlus::PionPlus()       ) projPDG=  211;
  else if (particle ==        G4KaonPlus::KaonPlus()       ) projPDG=  321;
  else if (particle ==       G4KaonMinus::KaonMinus()      ) projPDG= -321;
  else if (particle ==    G4KaonZeroLong::KaonZeroLong()  ||
           particle ==   G4KaonZeroShort::KaonZeroShort() ||
           particle ==        G4KaonZero::KaonZero()      ||
           particle ==    G4AntiKaonZero::AntiKaonZero()   )
  {
    if(G4UniformRand() > 0.5)                                projPDG=  311;
    else                                                     projPDG= -311;
  }
  else if (                      pZ > 0 && pA > 1          ) projPDG = 90000000+999*pZ+pA;
  else if (particle ==          G4Lambda::Lambda()         ) projPDG= 3122;
  else if (particle ==       G4SigmaPlus::SigmaPlus()      ) projPDG= 3222;
  else if (particle ==      G4SigmaMinus::SigmaMinus()     ) projPDG= 3112;
  else if (particle ==       G4SigmaZero::SigmaZero()      ) projPDG= 3212;
  else if (particle ==         G4XiMinus::XiMinus()        ) projPDG= 3312;
  else if (particle ==          G4XiZero::XiZero()         ) projPDG= 3322;
  else if (particle ==      G4OmegaMinus::OmegaMinus()     ) projPDG= 3334;
  else if (particle ==     G4AntiNeutron::AntiNeutron()    ) projPDG=-2112;
  else if (particle ==      G4AntiProton::AntiProton()     ) projPDG=-2212;
  else if (particle ==        G4MuonPlus::MuonPlus()       ) projPDG=  -13;
  else if (particle ==       G4MuonMinus::MuonMinus()      ) projPDG=   13;
  else if (particle ==           G4Gamma::Gamma()          ) projPDG=   22;
  else if (particle ==        G4Electron::Electron()       ) projPDG=   11;
  else if (particle ==        G4Positron::Positron()       ) projPDG=  -11;
  else if (particle ==      G4NeutrinoMu::NeutrinoMu()     ) projPDG=   14;
  else if (particle ==  G4AntiNeutrinoMu::AntiNeutrinoMu() ) projPDG=  -14;
  else if (particle ==      G4AntiLambda::AntiLambda()     ) projPDG=-3122;
  else if (particle ==   G4AntiSigmaPlus::AntiSigmaPlus()  ) projPDG=-3222;
  else if (particle ==  G4AntiSigmaMinus::AntiSigmaMinus() ) projPDG=-3112;
  else if (particle ==   G4AntiSigmaZero::AntiSigmaZero()  ) projPDG=-3212;
  else if (particle ==     G4AntiXiMinus::AntiXiMinus()    ) projPDG=-3312;
  else if (particle ==      G4AntiXiZero::AntiXiZero()     ) projPDG=-3322;
  else if (particle ==  G4AntiOmegaMinus::AntiOmegaMinus() ) projPDG=-3334;
  else if (particle ==       G4NeutrinoE::NeutrinoE()      ) projPDG=   12;
  else if (particle ==   G4AntiNeutrinoE::AntiNeutrinoE()  ) projPDG=  -12;
  else if (particle ==         G4TauPlus::TauPlus()        ) projPDG=  -15;
  else if (particle ==        G4TauMinus::TauMinus()       ) projPDG=   15;
  else if (particle ==     G4NeutrinoTau::NeutrinoTau()    ) projPDG=   16;
  else if (particle == G4AntiNeutrinoTau::AntiNeutrinoTau()) projPDG=  -16;
  G4int aProjPDG=std::abs(projPDG);
#ifdef debug
  G4int prPDG=particle->GetPDGEncoding();
  G4cout<<"G4QInelastic::PostStepDoIt: projPDG="<<projPDG<<", stPDG="<<prPDG<<G4endl;
#endif
  if(!projPDG)
  {
    G4cerr<<"---Warning---G4QInelastic::PostStepDoIt:Undefined interacting hadron"<<G4endl;
    return 0;
  }
  G4int EPIM=ElProbInMat.size();
#ifdef debug
  G4cout<<"G4QInel::PostStDoIt: m="<<EPIM<<",n="<<nE<<",T="<<ElProbInMat[EPIM-1]<<G4endl;
#endif
  G4int i=0;
  if(EPIM>1)
  {
    G4double rnd = ElProbInMat[EPIM-1]*G4UniformRand();
    for(i=0; i<nE; ++i)
    {
#ifdef debug
      G4cout<<"G4QInelastic::PostStepDoIt:E["<<i<<"]="<<ElProbInMat[i]<<",r="<<rnd<<G4endl;
#endif
      if (rnd<ElProbInMat[i]) break;
    }
    if(i>=nE) i=nE-1;                        // Top limit for the Element
  }
  G4Element* pElement=(*theElementVector)[i];
  Z=static_cast<G4int>(pElement->GetZ());
#ifdef debug
    G4cout<<"G4QInelastic::PostStepDoIt: i="<<i<<", Z(element)="<<Z<<G4endl;
#endif
  if(Z<=0)
  {
    G4cerr<<"---Warning---G4QInelastic::PostStepDoIt: Element with Z="<<Z<<G4endl;
    if(Z<0) return 0;
  }
  std::vector<G4double>* SPI = IsoProbInEl[i];// Vector of summedProbabilities for isotopes
  std::vector<G4int>* IsN = ElIsoN[i];     // Vector of "#of neutrons" in the isotope El[i]
  G4int nofIsot=SPI->size();               // #of isotopes in the element i
#ifdef debug
  G4cout<<"G4QInel::PosStDoIt:n="<<nofIsot<<",T="<<(*SPI)[nofIsot-1]<<G4endl;
#endif
  G4int j=0;
  if(nofIsot>1)
  {
    G4double rndI=(*SPI)[nofIsot-1]*G4UniformRand(); // Randomize the isotop of the Element
    for(j=0; j<nofIsot; ++j)
    {
#ifdef debug
      G4cout<<"G4QInelastic::PostStepDoIt: SP["<<j<<"]="<<(*SPI)[j]<<", r="<<rndI<<G4endl;
#endif
      if(rndI < (*SPI)[j]) break;
    }
    if(j>=nofIsot) j=nofIsot-1;            // Top limit for the isotope
  }
  G4int N =(*IsN)[j]; ;                    // Randomized number of neutrons
#ifdef debug
  G4cout<<"G4QInelastic::PostStepDoIt: Z="<<Z<<", j="<<i<<", N(isotope)="<<N<<G4endl;
#endif
  G4double kinEnergy= projHadron->GetKineticEnergy();
  G4ParticleMomentum dir = projHadron->GetMomentumDirection();
  if(projPDG==22 && Z==1 && !N && Momentum<150.)
  {
#ifdef debug
    G4cerr<<"---LowPHOTON---G4QInelastic::PSDoIt: Called for zero Cross-section"<<G4endl;
#endif
    //Do Nothing Action insead of the reaction
    aParticleChange.ProposeEnergy(kinEnergy);
    aParticleChange.ProposeLocalEnergyDeposit(0.);
    aParticleChange.ProposeMomentumDirection(dir);
    aParticleChange.ProposeTrackStatus(fAlive);
    return G4VDiscreteProcess::PostStepDoIt(track,step);
  }
  if(N<0)
  {
    G4cerr<<"-Warning-G4QInelastic::PostStepDoIt: Isotope with Z="<<Z<<", 0>N="<<N<<G4endl;
    return 0;
  }
  nOfNeutrons=N;                           // Remember it for the energy-momentum check
  //G4double am=Z+N;
  G4int am=Z+N;
  //G4double dd=0.025;
  //G4double sr=std::sqrt(am);
  //G4double dsr=0.01*(sr+sr);
  //if(dsr<dd)dsr=dd;
  //G4double medRA=mediRatio*G4QThd(am);
  G4double medRA=mediRatio;
  if(manualFlag) G4QNucleus::SetParameters(freeNuc,freeDib,clustProb,medRA); // ManualP
  //else if(projPDG==-2212)G4QNucleus::SetParameters(1.-dsr-dsr,dd+dd,5.,9.);//aP ClustPars
  //else if(projPDG==-211) G4QNucleus::SetParameters(.67-dsr,.32-dsr,5.,9.); //Pi-ClustPars
  //else if(projPDG ==-2212 || projPDG ==-2112)
  //                       G4QNucleus::SetParameters(1.-dsr-dsr,dd+dd,5.,1.);//aP ClustPars
  //else if(projPDG ==-211 ||projPDG == 211 )
  //                       G4QNucleus::SetParameters(.67-dsr,.32-dsr,5.,1.); //Pi-ClustPars
#ifdef debug
  G4cout<<"G4QInelastic::PostStepDoIt: N="<<N<<" for element with Z="<<Z<<G4endl;
#endif
  aParticleChange.Initialize(track);
  G4double weight = track.GetWeight();
#ifdef debug
  G4cout<<"G4QInelastic::PostStepDoIt: weight="<<weight<<G4endl;
#endif
  if(photNucBias!=1.)      weight/=photNucBias;
  else if(weakNucBias!=1.) weight/=weakNucBias;
  G4double localtime = track.GetGlobalTime();
#ifdef debug
  G4cout<<"G4QInelastic::PostStepDoIt: localtime="<<localtime<<G4endl;
#endif
  G4ThreeVector position = track.GetPosition();
  G4TouchableHandle trTouchable = track.GetTouchableHandle();
#ifdef debug
  G4cout<<"G4QInelastic::PostStepDoIt: position="<<position<<G4endl;
#endif
  //
  G4int targPDG=90000000+Z*1000+N;            // PDG Code of the target nucleus
  G4QPDGCode targQPDG(targPDG);
  G4double tgM=targQPDG.GetMass();            // Target mass
  G4double tM=tgM;                            // Target mass (copy to be changed)
  G4QHadronVector* output=new G4QHadronVector;// Prototype of Fragm Output G4QHadronVector
  G4double absMom = 0.;                       // Prototype of absorbed by nucleus Moment
  G4QHadronVector* leadhs=new G4QHadronVector;// Prototype of QuasmOutput G4QHadronVector
  G4LorentzVector lead4M(0.,0.,0.,0.);        // Prototype of LeadingQ 4-momentum
#ifdef debug
  G4cout<<"G4QInelastic::PostStepDoIt: projPDG="<<aProjPDG<<", targPDG="<<targPDG<<G4endl;
#endif
  //  
  // Leptons with photonuclear
  // Lepto-nuclear case with the equivalent photon algorithm. @@InFuture + NC (?)
  //
  if (aProjPDG == 11 || aProjPDG == 13 || aProjPDG == 15)
  {
#ifdef debug
    G4cout<<"G4QInelastic::PostStDoIt:startSt="<<aParticleChange.GetTrackStatus()<<G4endl;
#endif
    G4VQCrossSection* CSmanager=G4QElectronNuclearCrossSection::GetPointer();
    G4double ml=me;
    G4double ml2=me2;
    if(aProjPDG== 13)
    {
      CSmanager=G4QMuonNuclearCrossSection::GetPointer();
      ml=mu;
      ml2=mu2;
    }
    if(aProjPDG== 15)
    {
      CSmanager=G4QTauNuclearCrossSection::GetPointer();
      ml=mt;
      ml2=mt2;
    }
    // @@ Probably this is not necessary any more (?)
    if(!aProjPDG) G4cout<<"-Warning-G4QInelast::PostStepDoIt:(2) projectile PDG=0"<<G4endl;
    G4double xSec=CSmanager->GetCrossSection(false,Momentum,Z,N,aProjPDG);// Recalculate XS
    // @@ check a possibility to separate p, n, or alpha (!)
    if(Z==1 && !N && Momentum<150.) xSec=0.;
#ifdef debug
    G4cout<<"-Forse-G4QInel::PStDoIt:P="<<Momentum<<",PDG="<<projPDG<<",xS="<<xSec<<G4endl;
#endif
    if(xSec <= 0.) // The cross-section is 0 -> Do Nothing
    {
#ifdef debug
      G4cerr<<"---OUT---G4QInelastic::PSDoIt: Called for zero Cross-section"<<G4endl;
#endif
      //Do Nothing Action insead of the reaction
      aParticleChange.ProposeEnergy(kinEnergy);
      aParticleChange.ProposeLocalEnergyDeposit(0.);
      aParticleChange.ProposeMomentumDirection(dir);
      aParticleChange.ProposeTrackStatus(fAlive);
      delete output;
      return G4VDiscreteProcess::PostStepDoIt(track,step);
    }
    G4double photonEnergy = CSmanager->GetExchangeEnergy(); // Energy of EqivExchangePart
#ifdef debug
    G4cout<<"G4QInel::PStDoIt:kE="<<kinEnergy<<",dir="<<dir<<",phE="<<photonEnergy<<G4endl;
#endif
    if( kinEnergy < photonEnergy || photonEnergy < 0.)
    {
      //Do Nothing Action insead of the reaction
      G4cerr<<"--G4QInelastic::PSDoIt: photE="<<photonEnergy<<">leptE="<<kinEnergy<<G4endl;
      aParticleChange.ProposeEnergy(kinEnergy);
      aParticleChange.ProposeLocalEnergyDeposit(0.);
      aParticleChange.ProposeMomentumDirection(dir);
      aParticleChange.ProposeTrackStatus(fAlive);
      delete output;
      return G4VDiscreteProcess::PostStepDoIt(track,step);
    }
    G4double photonQ2 = CSmanager->GetExchangeQ2(photonEnergy);// Q2(t) of EqivExchangePart
    G4double W=photonEnergy-photonQ2/dM;// HadronicEnergyFlow (W-energy) for virtual photon
    if(tM<999.) W-=mPi0+mPi0s/dM;       // Pion production threshold for a nucleon target
    if(W<0.) 
    {
      //Do Nothing Action insead of the reaction
#ifdef debug
      G4cout<<"--G4QInelastic::PostStepDoIt:(lN) negative equivalent energy W="<<W<<G4endl;
#endif
      aParticleChange.ProposeEnergy(kinEnergy);
      aParticleChange.ProposeLocalEnergyDeposit(0.);
      aParticleChange.ProposeMomentumDirection(dir);
      aParticleChange.ProposeTrackStatus(fAlive);
      delete output;
      return G4VDiscreteProcess::PostStepDoIt(track,step);
    }
    // Update G4VParticleChange for the scattered muon
    G4VQCrossSection* thePhotonData=G4QPhotonNuclearCrossSection::GetPointer();
    G4double sigNu=thePhotonData->GetCrossSection(true,photonEnergy,Z,N,22);//Integrated XS
    G4double sigK =thePhotonData->GetCrossSection(true, W, Z, N, 22);       // Real XS
    G4double rndFraction = CSmanager->GetVirtualFactor(photonEnergy, photonQ2);
    if(sigNu*G4UniformRand()>sigK*rndFraction) 
    {
      //NothingToDo Action insead of the reaction
#ifdef debug
      G4cout<<"-DoNoth-G4QInelastic::PostStepDoIt: probab. correction - DoNothing"<<G4endl;
#endif
      aParticleChange.ProposeEnergy(kinEnergy);
      aParticleChange.ProposeLocalEnergyDeposit(0.);
      aParticleChange.ProposeMomentumDirection(dir);
      aParticleChange.ProposeTrackStatus(fAlive);
      delete output;
      return G4VDiscreteProcess::PostStepDoIt(track,step);
    }
    G4double iniE=kinEnergy+ml;          // Initial total energy of the lepton
    G4double finE=iniE-photonEnergy;     // Final total energy of the lepton
#ifdef pdebug
    G4cout<<"G4QInelastic::PoStDoIt:E="<<iniE<<",lE="<<finE<<"-"<<ml<<"="<<finE-ml<<G4endl;
#endif
    aParticleChange.ProposeEnergy(finE-ml);
    if(finE<=ml)                         // Secondary lepton (e/mu/tau) at rest disappears
    {
      aParticleChange.ProposeEnergy(0.);
      if(aProjPDG== 11) aParticleChange.ProposeTrackStatus(fStopAndKill);
      else aParticleChange.ProposeTrackStatus(fStopButAlive);
      aParticleChange.ProposeMomentumDirection(dir);
    }
    else aParticleChange.ProposeTrackStatus(fAlive);
    G4double iniP=std::sqrt(iniE*iniE-ml2); // Initial momentum of the electron
    G4double finP=std::sqrt(finE*finE-ml2); // Final momentum of the electron
    G4double cost=(iniE*finE-ml2-photonQ2/2)/iniP/finP; // cos(scat_ang_of_lepton)
#ifdef pdebug
    G4cout<<"G4QI::PSDoIt:Q2="<<photonQ2<<",ct="<<cost<<",Pi="<<iniP<<",Pf="<<finP<<G4endl;
#endif
    if(cost>1.) cost=1.;                 // To avoid the accuracy of calculation problem
    if(cost<-1.) cost=-1.;               // To avoid the accuracy of calculation problem
    //
    // Scatter the lepton ( @@ make the same thing for real photons)
    // At this point we have photonEnergy and photonQ2 (with notDefinedPhi)->SelectProjPart
    G4double absEn = G4QThd(am)*GeV;     // @@(b) Mean Energy Absorbed by a Nucleus
    //G4double absEn = std::pow(am,third)*GeV;  // @@(b) Mean Energy Absorbed by a Nucleus
    if(am>1 && absEn < photonEnergy)     // --> the absorption of energy can happen
    //if(absEn < photonEnergy)             // --> the absorption of energy can happen
    {
      G4double abtEn = absEn+hdM;        // @@(b) MeanEnergyAbsorbed by a nucleus (+M_N)
      G4double abEn2 = abtEn*abtEn;      // Squared absorbed Energy + MN
      G4double abMo2 = abEn2-hdM2;       // Squared absorbed Momentum of compound system
      G4double phEn2 = photonEnergy*photonEnergy;
      G4double phMo2 = phEn2+photonQ2;   // Squared momentum of primary virtual photon
      G4double phMo  = std::sqrt(phMo2); // Momentum of the primary virtual photon
      absMom         = std::sqrt(abMo2); // Absorbed Momentum
      if(absMom < phMo)                  // --> the absorption of momentum can happen
      {
        G4double dEn = photonEnergy - absEn; // Leading energy
        G4double dMo = phMo - absMom;    // Leading momentum
        G4double sF  = dEn*dEn - dMo*dMo;// s of leading particle
#ifdef ppdebug
        G4cout<<"-PhotoAbsorption-G4QIn::PStDoIt: sF="<<sF<<",phEn="<<photonEnergy<<G4endl;
#endif
        if(sF > stmPi)                   // --> Leading fragmentation is possible
        {
          photonEnergy = absEn;          // New value of the photon energy
          photonQ2=abMo2-absEn*absEn;    // New value of the photon Q2
          absEn = dEn;                   // Put energy of leading particle to absEn (!)
        }
        else absMom=0.;                  // Flag that nothing has happened
      }
      else absMom=0.;                    // Flag that nothing has happened
    }
    // ------------- End of ProjPart selection
    //
    // Scattering in respect to the derection of the incident muon is made impicitly:
    G4ThreeVector ort=dir.orthogonal();  // Not normed orthogonal vector (!) (to dir)
    G4ThreeVector ortx = ort.unit();     // First unit vector orthogonal to the direction
    G4ThreeVector orty = dir.cross(ortx);// Second unit vector orthoganal to the direction
    G4double sint=std::sqrt(1.-cost*cost); // Perpendicular component
    G4double phi=twopi*G4UniformRand();  // phi of scattered electron
    G4double sinx=sint*std::sin(phi);    // x perpendicular component
    G4double siny=sint*std::cos(phi);    // y perpendicular component
    G4ThreeVector findir=cost*dir+sinx*ortx+siny*orty;
    aParticleChange.ProposeMomentumDirection(findir); // new direction for the lepton
#ifdef pdebug
    G4cout<<"G4QInelastic::PostStepDoIt: E="<<aParticleChange.GetEnergy()<<"="<<finE<<"-"
          <<ml<<", d="<<*aParticleChange.GetMomentumDirection()<<","<<findir<<G4endl;
#endif
    G4ThreeVector photon3M=iniP*dir-finP*findir;// 3D total momentum of photon
    if(absMom)                           // Photon must be reduced & LeadingSyst fragmented
    {
      G4double ptm=photon3M.mag();                    // 3M of the virtual photon
#ifdef ppdebug
      G4cout<<"-Absorption-G4QInelastic::PostStepDoIt: ph3M="<<photon3M<<", eIn3M="
            <<iniP*dir<<", eFin3M="<<finP*findir<<", abs3M="<<absMom<<"<ptm="<<ptm<<G4endl;
#endif
      G4ThreeVector lead3M=photon3M*(ptm-absMom)/ptm; // Keep the direction for leading Q
      photon3M-=lead3M; // Reduced photon Momentum (photEn already = absEn)
      proj4M=G4LorentzVector(lead3M,absEn); // 4-momentum of leading System
#ifdef ppdebug
      G4cout<<"-->G4QI::PoStDoIt: new sF="<<proj4M.m2()<<", lead4M="<<proj4M<<G4endl;
#endif
      lead4M=proj4M;                        // Remember 4-mom for the total 4-momentum
      G4Quasmon* pan= new G4Quasmon(G4QContent(1,1,0,1,1,0),proj4M);// ---> DELETED -->---+
      try                                                           //                    |
      {                                                             //                    |
        if(leadhs) delete leadhs;                                   //                    |
        G4QNucleus vac(90000000);                                   //                    |
        leadhs=pan->Fragment(vac,1);  // DELETED after it is copied to output vector      |
      }                                                             //                    |
      catch (G4QException& error)                                   //                    |
      {                                                             //                    |
        G4cerr<<"***G4QInelastic::PostStepDoIt: G4Quasmon Exception is catched"<<G4endl;//|
        // G4Exception("G4QInelastic::PostStepDoIt:","72",FatalException,"QuasmonCrash");  //|
        G4Exception("G4QInelastic::PostStepDoIt()","HAD_CHPS_0072",
                    FatalException, "QuasmonCrash");
      }                                                             //                    |
      delete pan;                              // Delete the Nuclear Environment <----<---+
#ifdef ppdebug
      G4cout<<"G4QInel::PStDoIt: l4M="<<proj4M<<proj4M.m2()<<", N="<<leadhs->size()<<",pt="
            <<ptm<<",pa="<<absMom<<",El="<<absEn<<",Pl="<<ptm-absMom<<G4endl;
#endif
    }
    else
    {
      G4int qNH=0;
      if(leadhs) qNH=leadhs->size();
      if(qNH) for(G4int iq=0; iq<qNH; iq++) delete (*leadhs)[iq];
      if(leadhs) delete leadhs;
      leadhs=0;
    }
    projPDG=22;
    proj4M=G4LorentzVector(photon3M,photonEnergy);
#ifdef debug
    G4cout<<"G4QInelastic::PostStDoIt: St="<<aParticleChange.GetTrackStatus()<<", g4m="
          <<proj4M<<", lE="<<finE<<", lP="<<finP*findir<<", d="<<findir.mag2()<<G4endl;
#endif

  //
  // neutrinoNuclear interactions (nu_e, nu_mu only)
  //
  }
  else if (aProjPDG == 12 || aProjPDG == 14)
  {
    kinEnergy= projHadron->GetKineticEnergy()/MeV; // Total energy of the neutrino
    G4double dKinE=kinEnergy+kinEnergy;  // doubled energy for s calculation
#ifdef debug
    G4cout<<"G4QInelastic::PostStDoIt: 2*nuEnergy="<<dKinE<<"(MeV), PDG="<<projPDG<<G4endl;
#endif
    dir = projHadron->GetMomentumDirection(); // unit vector
    G4double ml  = mu;
    G4double ml2 = mu2;
    //G4double mlN = muN;
    G4double mlN2= muN2;
    G4double fmlN= fmuN;
    G4double mlsN= musN;
    //G4double mlP = muP;
    G4double mlP2= muP2;
    G4double fmlP= fmuP;
    G4double mlsP= musP;
    G4double mldM= mudM;
    //G4double mlD = muD;
    G4double mlD2= muD2;
    if(aProjPDG==12)
    {
      ml  = me;
      ml2 = me2;
      //mlN = meN;
      mlN2= meN2;
      fmlN= fmeN;
      mlsN= mesN;
      //mlP = meP;
      mlP2= meP2;
      fmlP= fmeP;
      mlsP= mesP;
      mldM= medM;
      //mlD = meD;
      mlD2= meD2;
    }
    G4VQCrossSection* CSmanager =G4QNuMuNuclearCrossSection::GetPointer(); // (nu,l)
    G4VQCrossSection* CSmanager2=G4QNuNuNuclearCrossSection::GetPointer(); // (nu,nu)
    proj4M=G4LorentzVector(dir*kinEnergy,kinEnergy);   // temporary
    G4bool nuanu=true;
    scatPDG=13;                          // Prototype = secondary scattered mu-
    if(projPDG==-14)
    {
      nuanu=false;                       // Anti-neutrino
      CSmanager=G4QANuMuNuclearCrossSection::GetPointer();  // (anu,mu+) CC @@ open
      CSmanager=G4QANuANuNuclearCrossSection::GetPointer(); // (anu,anu) NC @@ open
      scatPDG=-13;                       // secondary scattered mu+
    }
    else if(projPDG==12)
    {
      CSmanager=G4QNuENuclearCrossSection::GetPointer(); // @@ open (only CC is changed)
      scatPDG=11;                        // secondary scattered e-
    }
    else if(projPDG==-12)
    {
      nuanu=false;                       // anti-neutrino
      CSmanager=G4QANuENuclearCrossSection::GetPointer();   // (anu,e+) CC @@ open
      CSmanager=G4QANuANuNuclearCrossSection::GetPointer(); // (anu,anu) NC @@ open
      scatPDG=-11;                       // secondary scattered e+
    }
    // @@ Probably this is not necessary any more
    if(!projPDG) G4cout<<"-Warning-G4QInelastic::PostStepDoIt:(3)projectile PDG=0"<<G4endl;
    G4double xSec1=CSmanager->GetCrossSection(false,Momentum,Z,N,projPDG); //Recalculate XS
    G4double xSec2=CSmanager2->GetCrossSection(false,Momentum,Z,N,projPDG);//Recalculate XS
    G4double xSec=xSec1+xSec2;
    // @@ check a possibility to separate p, n, or alpha (!)
    if(xSec <= 0.) // The cross-section = 0 -> Do Nothing
    {
      G4cerr<<"G4QInelastic::PSDoIt:nuE="<<kinEnergy<<",X1="<<xSec1<<",X2="<<xSec2<<G4endl;
      //Do Nothing Action insead of the reaction
      aParticleChange.ProposeEnergy(kinEnergy);
      aParticleChange.ProposeLocalEnergyDeposit(0.);
      aParticleChange.ProposeMomentumDirection(dir);
      aParticleChange.ProposeTrackStatus(fAlive);
      delete output;
      return G4VDiscreteProcess::PostStepDoIt(track,step);
    }
    G4bool secnu=false;
    if(xSec*G4UniformRand()>xSec1)               // recover neutrino/antineutrino
    {
      if(scatPDG>0) scatPDG++;
      else          scatPDG--;
      secnu=true;
    }
    scat=true;                                   // event with changed scattered projectile
    G4double totCS1 = CSmanager->GetLastTOTCS(); // the last total cross section1(isotope?)
    G4double totCS2 = CSmanager2->GetLastTOTCS();// the last total cross section2(isotope?)
    G4double totCS  = totCS1+totCS2;             // the last total cross section (isotope?)
    if(std::fabs(xSec-totCS*millibarn)/xSec>.0001)
          G4cout<<"-Warning-G4QInelastic::PostStepDoIt: xS="<<xSec<<"# CS="<<totCS<<G4endl;
    G4double qelCS1 = CSmanager->GetLastQELCS(); // the last quasi-elastic cross section1
    G4double qelCS2 = CSmanager2->GetLastQELCS();// the last quasi-elastic cross section2
    G4double qelCS  = qelCS1+qelCS2;             // the last quasi-elastic cross section
    if(totCS - qelCS < 0.)                       // only at low energies
    {
      totCS  = qelCS;
      totCS1 = qelCS1;
      totCS2 = qelCS2;
    }
    // make different definitions for neutrino and antineutrino
    G4double mIN=mProt;                          // Just a prototype (for anu, Z=1, N=0)
    G4double mOT=mNeut;
    G4double OT=mlN2;
    //G4double mOT2=mNeut2;                        // Other formula: sd=s-mlsOT -Not used?-
    G4double mlOT=fmlN;
    G4double mlsOT=mlsN;
    if(secnu)
    {
      if(am*G4UniformRand()>Z)                   // Neutron target
      {
        targPDG-=1;                              // subtract neutron
        projPDG=2112;                            // neutron is going out
        mIN =mNeut;                     
        OT  =mNeut2;
        //mOT2=mNeut2;
        mlOT=0.;
        mlsOT=mNeut2;
      }
      else
      {
        targPDG-=1000;                           // subtract neutron
        projPDG=2212;                            // neutron is going out
        mOT  =mProt;
        OT   =mProt2;
        //mOT2 =mProt2;
        mlOT =0.;
        mlsOT=mProt2;
      }
      ml=0.;
      ml2=0.;
      mldM=0.;
      mlD2=mPPi2;
      G4QPDGCode temporary_targQPDG(targPDG); 
      targQPDG = temporary_targQPDG;
      G4double rM=targQPDG.GetMass();
      mIN=tM-rM;                                 // bounded in-mass of the neutron
      tM=rM;
    }
    else if(nuanu)
    {
      targPDG-=1;                                // Neutrino -> subtract neutron
      G4QPDGCode temporary_targQPDG(targPDG); 
      targQPDG = temporary_targQPDG;
      G4double rM=targQPDG.GetMass();
      mIN=tM-rM;                                 // bounded in-mass of the neutron
      tM=rM;
      mOT=mProt;
      OT=mlP2;
      //mOT2=mProt2;
      mlOT=fmlP;
      mlsOT=mlsP;
      projPDG=2212;                              // proton is going out
    }
    else
    {
      if(Z>1||N>0)                               // Calculate the splitted mass
      {
        targPDG-=1000;                           // Anti-Neutrino -> subtract proton
        G4QPDGCode temporary_targQPDG(targPDG); 
        targQPDG = temporary_targQPDG;
        G4double rM=targQPDG.GetMass();
        mIN=tM-rM;                               // bounded in-mass of the proton
        tM=rM;
      }
      else
      {
        targPDG=0;
        mIN=tM;
        tM=0.;
      }
      projPDG=2112;                              // neutron is going out
    }
    G4double s_value=mIN*(mIN+dKinE);            // s_value=(M_cm)^2=m2+2mE (m=targetMass,E=E_nu)
#ifdef debug
    G4cout<<"G4QInelastic::PostStDoIt: s="<<s_value<<" >? OT="<<OT<<", mlD2="<<mlD2<<G4endl;
#endif
    if(s_value<=OT)                                    // *** Do nothing solution ***
    {
      //Do NothingToDo Action insead of the reaction (@@ Can we make it common?)
      G4cout<<"G4QInelastic::PostStepDoIt: tooSmallFinalMassOfCompound: DoNothing"<<G4endl;
      aParticleChange.ProposeEnergy(kinEnergy);
      aParticleChange.ProposeLocalEnergyDeposit(0.);
      aParticleChange.ProposeMomentumDirection(dir);
      aParticleChange.ProposeTrackStatus(fAlive);
      delete output;
      return G4VDiscreteProcess::PostStepDoIt(track,step);
    }
#ifdef debug
    G4cout<<"G4QInelastic::PostStDoIt: Stop and kill the projectile neutrino"<<G4endl;
#endif
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill); // the initial neutrino is killed
    // There is no way back from here !
    if ( ((secnu || !nuanu || N) && totCS*G4UniformRand() < qelCS) || s_value < mlD2 ) 
    {   // Quasi-Elastic interaction
      G4double Q2=0.;                           // Simulate transferred momentum, in MeV^2
      if(secnu) Q2=CSmanager2->GetQEL_ExchangeQ2();
      else      Q2=CSmanager->GetQEL_ExchangeQ2();
#ifdef debug
      G4cout<<"G4QInelastic::PostStDoIt:QuasiEl(nu="<<secnu<<"),s="<<s_value<<",Q2="<<Q2<<G4endl;
#endif
      //G4double ds=s_value+s_value;              // doubled s_value
      G4double sqs=std::sqrt(s_value);          // M_cm
      G4double dsqs=sqs+sqs;                    // 2*M_cm
      G4double p_init=(s_value-mIN*mIN)/dsqs;   // initial momentum in CMS (checked MK)
      G4double dpi=p_init+p_init;               // doubled initial momentum in CMS
      G4double sd=s_value-mlsOT;                // s_value-ml2-mOT2 (mlsOT=m^2_neut+m^2_lept)
      G4double qo2=(sd*sd-mlOT)/dsqs;           // squared momentum of secondaries in CMS
      G4double qo=std::sqrt(qo2);               // momentum of secondaries in CMS
      G4double cost=(dpi*std::sqrt(qo2+ml2)-Q2-ml2)/dpi/qo; // cos(theta) in CMS (chck MK)
      G4LorentzVector t4M(0.,0.,0.,mIN);        // 4mom of the effective target
      G4LorentzVector c4M=t4M+proj4M;           // 4mom of the compound system
      t4M.setT(mOT);                            // now it is 4mom of the outgoing nucleon
      scat4M=G4LorentzVector(0.,0.,0.,ml);      // 4mom of the scattered muon
      if(!G4QHadron(c4M).RelDecayIn2(scat4M, t4M, proj4M, cost, cost))
      {
        G4cerr<<"G4QIn::PStD:c4M="<<c4M<<sqs<<",mM="<<ml<<",tM="<<mOT<<",c="<<cost<<G4endl;
        // throw G4QException("G4QInelastic::HadronizeQuasm: Can't dec QE nu,lept Compound");
        G4Exception("G4QInelastic::PostStepDoIt()", "HAD_CHPS_0000",
                    FatalException, "Hadronize quasmon: Can't dec QE nu,lept Compound");
      }
      proj4M=t4M;                               // 4mom of the new projectile nucleon
    }
    else                                        // ***** Non Quasi Elastic interaction
    {
      // Recover the target PDG Code if the quasi-elastic scattering did not happened
      if      ( (secnu && projPDG == 2212) || (!secnu && projPDG == 2112) ) targPDG+=1;
      else if ( (secnu && projPDG == 2112) || (!secnu && projPDG == 2212) ) targPDG+=1000;
      G4double Q2=0;                            // Simulate transferred momentum, in MeV^2
      if(secnu) Q2=CSmanager->GetNQE_ExchangeQ2();
      else      Q2=CSmanager2->GetNQE_ExchangeQ2();
#ifdef debug
      G4cout<<"G4QInel::PStDoIt: MultiPeriferal s="<<s_value<<",Q2="<<Q2<<",T="<<targPDG<<G4endl;
#endif
      if(secnu) projPDG=CSmanager2->GetExchangePDGCode();// PDG Code of the effective gamma
      else      projPDG=CSmanager->GetExchangePDGCode(); // PDG Code of the effective pion
      //@@ Temporary made only for direct interaction and for N=3 (good for small Q2)
      //@@ inFuture use N=GetNPartons and directFraction=GetDirectPart, @@ W2...
      G4double r=G4UniformRand();
      G4double r1=0.5;                          // (1-x)
      if(r<0.5)      r1=std::sqrt(r+r)*(.5+.1579*(r-.5));
      else if(r>0.5) r1=1.-std::sqrt(2.-r-r)*(.5+.1579*(.5-r));
      G4double xn=1.-mldM/Momentum;             // Normalization of (1-x) [x>mldM/Mom]
      G4double x1=xn*r1;                        // (1-x)
      G4double x=1.-x1;                         // x=2k/M
      //G4double W2=(hdM2+Q2/x)*x1;               // W2 candidate
      G4double mx=hdM*x;                        // Part of the target to interact with
      G4double we=Q2/(mx+mx);                   // transfered energy
      if(we>=kinEnergy-ml-.001) we=kinEnergy-ml-.0001; // safety to avoid nan=sqrt(neg)
      G4double pot=kinEnergy-we;                // energy of the secondary lepton
      G4double mlQ2=ml2+Q2;
      G4double cost=(pot-mlQ2/dKinE)/std::sqrt(pot*pot-ml2); // LS cos(theta)
      if(std::fabs(cost)>1)
      {
#ifdef debug
        G4cout<<"*G4QInelastic::PostStDoIt: cost="<<cost<<", Q2="<<Q2<<", nu="<<we<<", mx="
              <<mx<<", pot="<<pot<<", 2KE="<<dKinE<<G4endl;
#endif
        if(cost>1.) cost=1.;
        else        cost=-1.;
        pot=mlQ2/dKinE+dKinE*ml2/mlQ2;          // extreme output momentum
      }
      G4double lEn=std::sqrt(pot*pot+ml2);      // Lepton energy
      G4double lPl=pot*cost;                    // Lepton longitudinal momentum
      G4double lPt=pot*std::sqrt(1.-cost*cost); // Lepton transverse momentum
      std::pair<G4double,G4double> d2d=Random2DDirection(); // Randomize phi
      G4double lPx=lPt*d2d.first;
      G4double lPy=lPt*d2d.second;
      G4ThreeVector vdir=proj4M.vect();         // 3D momentum of the projectile
      G4ThreeVector vz= vdir.unit();            // Ort in the direction of the projectile
      G4ThreeVector vv= vz.orthogonal();        // Not normed orthogonal vector (!)
      G4ThreeVector vx= vv.unit();              // First ort orthogonal to the direction
      G4ThreeVector vy= vz.cross(vx);           // Second ort orthoganal to the direction
      G4ThreeVector lP= lPl*vz+lPx*vx+lPy*vy;   // 3D momentum of the scattered lepton
      scat4M=G4LorentzVector(lP,lEn);           // 4mom of the scattered lepton
      proj4M-=scat4M;                           // 4mom of the W/Z (effective pion/gamma)
#ifdef debug
      G4cout<<"G4QInelastic::PostStDoIt: proj4M="<<proj4M<<", ml="<<ml<<G4endl;
#endif
      // Check that the en/mom transfer is possible, if not -> elastic
      G4int fintPDG=targPDG;                    // Prototype for the compound nucleus
      if(!secnu)
      {
        if(projPDG<0) fintPDG-= 999;
        else          fintPDG+= 999;
      }
      G4double fM=G4QPDGCode(fintPDG).GetMass();// compound nucleus Mass (MeV)
      G4double fM2=fM*fM;
      G4LorentzVector tg4M=G4LorentzVector(0.,0.,0.,tgM);
      G4LorentzVector c4M=tg4M+proj4M;
#ifdef debug
      G4cout<<"G4QInelastic::PSDI:fM2="<<fM2<<"<? mc4M="<<c4M.m2()<<",dM="<<fM-tgM<<G4endl;
#endif
      if(fM2>=c4M.m2())                         // Elastic scattering should be done
      {
        G4LorentzVector tot4M=tg4M+proj4M+scat4M; // recover the total 4-momentum
        s_value=tot4M.m2();
        G4double fs=s_value-fM2-ml2;
        G4double fMl=fM2*ml2;
        G4double hQ2max=(fs*fs/2-fMl-fMl)/s_value; // Maximum possible Q2/2
        cost=1.-Q2/hQ2max;                // cos(theta) in CMS (use MultProd Q2)
#ifdef debug
        G4cout<<"G4QI::PSDI:ct="<<cost<<",Q2="<<Q2<<",hQ2="<<hQ2max<<",4M="<<tot4M<<G4endl;
#endif
        G4double acost=std::fabs(cost);
        if(acost>1.)
        {
          if(acost>1.001) G4cout<<"-Warning-G4QInelastic::PostStDoIt: cost="<<cost<<G4endl;
          if     (cost> 1.) cost= 1.;
          else if(cost<-1.) cost=-1.;
        }
        G4LorentzVector reco4M=G4LorentzVector(0.,0.,0.,fM); // 4mom of the recoilNucleus
        scat4M=G4LorentzVector(0.,0.,0.,ml); // 4mom of the scatteredLepton
        G4LorentzVector dir4M=tot4M-G4LorentzVector(0.,0.,0.,(tot4M.e()-ml)*.01);
        if(!G4QHadron(tot4M).RelDecayIn2(scat4M, reco4M, dir4M, cost, cost))
        {
          G4cerr<<"G4QI::PSDI:t4M="<<tot4M<<",lM="<<ml<<",rM="<<fM<<",cost="<<cost<<G4endl;
          //G4Exception("G4QInelastic::PostStepDoIt:","027",FatalException,"ElasticDecay");
        }
#ifdef debug
        G4cout<<"G4QIn::PStDoI:l4M="<<scat4M<<"+r4M="<<reco4M<<"="<<scat4M+reco4M<<G4endl;
#endif
        // ----------------------------------------------------
        G4ParticleDefinition* theDefinition=0; // Prototype of a particle for E-Secondaries
        // Fill scattered lepton
        if     (scatPDG==-11) theDefinition = G4Positron::Positron();
        else if(scatPDG== 11) theDefinition = G4Electron::Electron();
        else if(scatPDG== 13) theDefinition = G4MuonMinus::MuonMinus();
        else if(scatPDG==-13) theDefinition = G4MuonPlus::MuonPlus();
        //else if(scatPDG== 15) theDefinition = G4TauMinus::TauMinus();
        //else if(scatPDG==-15) theDefinition = G4TauPlus::TauPlus();
        if     (scatPDG==-12) theDefinition = G4AntiNeutrinoE::AntiNeutrinoE();
        else if(scatPDG== 12) theDefinition = G4NeutrinoE::NeutrinoE();
        else if(scatPDG== 14) theDefinition = G4NeutrinoMu::NeutrinoMu();
        else if(scatPDG==-14) theDefinition = G4AntiNeutrinoMu::AntiNeutrinoMu();
        //else if(scatPDG== 16) theDefinition = G4NeutrinoTau::NeutrinoTau();
        //else if(scatPDG==-16) theDefinition = G4AntiNeutrinoTau::AntiNeutrinoTau();
        else  G4cout<<"-Warning-G4QInelastic::PostStDoIt: UnknownLepton="<<scatPDG<<G4endl;
        G4DynamicParticle* theScL = new G4DynamicParticle(theDefinition,scat4M);
        G4Track* scatLep = new G4Track(theScL, localtime, position ); //    scattered
        scatLep->SetWeight(weight);                                   //    weighted
        scatLep->SetTouchableHandle(trTouchable);                     //    residual
        aParticleChange.AddSecondary(scatLep);                        //    lepton
        // Fill residual nucleus
        if     (fintPDG==90000001) theDefinition = G4Neutron::Neutron(); // neutron
        else if(fintPDG==90001000) theDefinition = G4Proton::Proton();   // proton
        else                                                             // ion
        {
          G4int fm=static_cast<G4int>(fintPDG/1000000);               // Strange part
          G4int ZN=fintPDG-1000000*fm;
          G4int rZ=static_cast<G4int>(ZN/1000);
          G4int rA=ZN-999*rZ;
          theDefinition = G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,0,rZ);
        }
        G4DynamicParticle* theReN = new G4DynamicParticle(theDefinition,reco4M);
        G4Track* scatReN = new G4Track(theReN, localtime, position ); //    scattered
        scatReN->SetWeight(weight);                                   //    weighted
        scatReN->SetTouchableHandle(trTouchable);                     //    residual
        aParticleChange.AddSecondary(scatReN);                        //    nucleus
        delete output;
        return G4VDiscreteProcess::PostStepDoIt(track, step);
      }
    } // End of non-quasi-elastic interaction
  //
  // quasi-elastic (& pickup process) for p+A(Z,N)
  //
  }
  //else if ((projPDG == 2212 || projPDG == 2112) && Z > 0 && N > 0) // Never (in Fragm!)
  else if(2>3) // Possibility is closed as quasi-free scattering is made in fragmentation
  {
    if(momentum<500. && projPDG == 2112) // @@ It's reasonable to add proton capture too !
    {
      G4LorentzVector tot4M=G4LorentzVector(0.,0.,0.,tgM)+proj4M;
      G4double totM2=tot4M.m2();
      G4int tZ=Z;
      G4int tN=N;
      if(projPDG == 2212) tZ++;
      else                tN++; // @@ Now a neutron is the only alternative
      if(momentum<100.)         // Try the production threshold (@@ Why 5 MeV?)
      {
        G4QPDGCode FakeP(2112);
        //G4int m1p=tZ-1;
        G4int m2n=tN-2;
        // @@ For the secondary proton the Coulomb barrier should be taken into account
        //G4double mRP= FakeP.GetNuclMass(m1p,tN,0);    // Residual to the secondary proton
        G4double mR2N= FakeP.GetNuclMass(tZ,m2n,0);   // Residual to two secondary neutrons
        //G4double mRAl= FakeP.GetNuclMass(tZ-2,m2n,0); // Residual to the secondary alpha
        //G4double mR2N= FakeP.GetNuclMass(m1p,tN-1,0); // Residual to the p+n secondaries
        G4double minM=mR2N+mNeut+mNeut;
        // Compare with other possibilities (sciped for acceleration)
        if(totM2 < minM*minM) //  DoNothing (under reasonable threshold)
        {
          aParticleChange.ProposeEnergy(kinEnergy);
          aParticleChange.ProposeLocalEnergyDeposit(0.);
          aParticleChange.ProposeMomentumDirection(dir);
          aParticleChange.ProposeTrackStatus(fAlive);
          return G4VDiscreteProcess::PostStepDoIt(track,step);
        }
      }
    }
    // Now the quasi-free and PickUp processes should be done:
    G4QuasiFreeRatios* qfMan=G4QuasiFreeRatios::GetPointer();
    std::pair<G4double,G4double> fief=qfMan->GetRatios(momentum, projPDG, Z, N);
    G4double qepart=fief.first*fief.second; // @@ charge exchange is not included @@
    // Keep QE part for the quasi-free nucleons
    const G4int lCl=3; // The last clProb[lCl]==1. by definition, MUST be increasing
    G4double clProb[lCl]={.65,.85,.95};// N/P,D,t/He3,He4, integroProbab for .65,.2,.1,.05
    if(qepart>0.45) qepart=.45; // Quasielastic part is too large - shrink
    //else qepart/=clProb[0];     // Add corresponding number of 2N, 3N, & 4N clusters
    qepart=qepart/clProb[0]-qepart;// Add QE for 2N, 3N, & 4N clusters (N is made in G4QF)
    G4double pickup=1.-qepart;  // Estimate the rest of the cross-section
    G4double thresh=100.;
    //if(momentum >thresh) pickup*=50./momentum/std::pow(G4double(Z+N),third);// 50 is Par
    if(momentum > thresh) pickup*=50./momentum/G4QThd(Z+N);// 50 is Par
    // pickup = 0.;               // To exclude the pickup process
    if (N) pickup+=qepart;
    else   pickup =qepart;
    G4double rnd=G4UniformRand();
#ifdef ldebug
    G4cout<<"-->G4QInelastic::PSD:QE[p("<<proj4M<<")+(Z="<<Z<<",N="<<N<<")]="<<qepart
          <<", pickup="<<pickup<<G4endl;
#endif
    if(rnd<pickup) // Make a quasi free scattering (out:A-1,h,N) @@ KinLim,CoulBar,PauliBl
    {
      G4double dmom=91.;  // Fermi momentum (proto default for a deuteron)
      G4double fmm=287.;  // @@ Can be A-dependent
      if(Z>1||N>1) dmom=fmm*std::pow(-std::log(G4UniformRand()),third);
      // Values to be redefined for quasi-elastic
      G4LorentzVector r4M(0.,0.,0.,0.);      // Prototype of 4mom of the residual nucleus
      G4LorentzVector n4M(0.,0.,0.,0.);      // Prototype of 4mom of the quasi-cluster
      G4int nPDG=90000001;                   // Prototype for quasiCluster mass calculation
      G4int restPDG=targPDG;                 // Prototype should be reduced by quasiCluster
      G4int rA=Z+N-1;                        // Prototype for the residualNucl definition
      G4int rZ=Z;                            // residZ: OK for the quasi-free neutron
      G4int nA=1;                            // Prototype for the quasi-cluster definition
      G4int nZ=0;                            // nA=1,nZ=0: OK for the quasi-free neutron
      G4double qM=mNeut;                     // Free mass of the quasi-free cluster
      G4int A = Z + N;                       // Baryon number of the nucleus
      if(rnd<qepart)
      {
       // ===> First decay a nucleus in a nucleon and a residual (A-1) nucleus
       // Calculate clusterProbabilities (n,p,d,t,he3,he4 now only, can use UpdateClusters)
       G4double base=1.;  // Base for randomization (can be reduced by totZ & totN)
       G4int max=lCl;     // Number of boundaries (can be reduced by totZ & totN)
       // Take into account that at least one nucleon must be left !
       if(Z<2 || N<2 || A < 6) base = clProb[--max]; // Alpha cluster is impossible
       if((Z > 1 && N < 2) || (Z < 2 && N > 1)) base=(clProb[max]+clProb[max-1])/2;// t/He3
       if ( (Z < 2 && N < 2) || A < 5) base=clProb[--max];// He3&t clusters are impossible
       if(A<3)           base=clProb[--max]; // Deuteron cluster is impossible
       //G4int cln=0;                        // Cluster#0 (Default for the selectedNucleon)
       G4int cln=1;                          // Cluster#1 (Default for the selectedDeutron)
       //if(max)                             // Not only nucleons are possible
       if(max>1)                             // Not only deuterons are possible
       {
         base-=clProb[0];                   // Exclude scattering on QF Nucleon
         G4double ran=+clProb[0]+base*G4UniformRand(); // Base can be reduced
         G4int ic=1;                        // Start from the smallest cluster boundary
         if(max>1) while(ic<max) if(ran>clProb[ic++]) cln=ic;
       }
       G4ParticleDefinition* theDefinition=0;// Prototype for qfNucleon
       G4bool cp1 = cln+2==A;               // A=ClusterBN+1 condition
       if(!cln || cp1)                      // Split in nucleon + (A-1) with Fermi momentum
       { 
        G4int nln=0;
        if(cln==2) nln=1;                   // @@ only for cp1: t/He3 choice from A=4
        // mass(A)=tM. Calculate masses of A-1 (rM) and mN (mNeut or mProt bounded mass)
        if ( ((!cln || cln == 2) && G4UniformRand()*(A-cln) > (N-nln)) || 
             ((cln == 3 || cln == 1) && Z > N) )
        {
          nPDG=90001000;                          // Update quasi-free nucleon PDGCode to P
          nZ=1;                                   // Change charge of the quasiFree nucleon
          qM=mProt;                               // Update quasi-free nucleon mass
          rZ--;                                   // Reduce the residual Z
          restPDG-=1000;                          // Reduce the residual PDGCode
        }
        else restPDG--;
        G4LorentzVector t4M(0.,0.,0.,tM);         // 4m of the target nucleus to be decayed
        G4double rM=G4QPDGCode(restPDG).GetMass();// Mass of the residual nucleus
        r4M=G4LorentzVector(0.,0.,0.,rM);         // 4mom of the residual nucleus
        G4double rM2=rM*rM;
        G4double nM=std::sqrt(rM2+tM*tM-(tM+tM)*std::sqrt(rM2+dmom*dmom));// M of q-nucleon
        n4M=G4LorentzVector(0.,0.,0.,nM);         // 4mom of the quasi-nucleon
#ifdef qedebug
        G4cout<<"G4QInel::PStDoIt:QE,p="<<dmom<<",tM="<<tM<<",R="<<rM<<",N="<<nM<<G4endl;
#endif
        if(!G4QHadron(t4M).DecayIn2(r4M, n4M))
        {
          //G4cerr<<"G4QInel::PostStDoIt:M="<<tM<<"<r="<<rM<<"+n="<<nM<<"="<<rM+nM<<G4endl;
          //G4Exception("G4QInelastic::PostStepDoIt()", "HAD_CHPS_0002",
          //           FatalException,"Hadronize quasmon: Can't Dec totNuc->QENuc+ResNuc");
          G4cerr<<"-W-G4QInel::PoStDoI:M="<<tM<<"<rM="<<rM<<"+nM="<<nM<<"="<<rM+nM<<G4endl;
          aParticleChange.ProposeEnergy(kinEnergy);
          aParticleChange.ProposeLocalEnergyDeposit(0.);
          aParticleChange.ProposeMomentumDirection(dir);
          aParticleChange.ProposeTrackStatus(fAlive);
          return G4VDiscreteProcess::PostStepDoIt(track,step);
        }
#ifdef qedebug
        G4cout<<"G4QIn::PStDoIt:QE-N, RA="<<r4M.rho()<<r4M<<",QN="<<n4M.rho()<<n4M<<G4endl;
#endif
        if(cp1 && cln)                           // Quasi-cluster case: swap the output
        {
          qM=rM;                                 // Scattering will be made on a cluster
          nln=nPDG;
          nPDG=restPDG;
          restPDG=nln;
          t4M=n4M;
          n4M=r4M;
          r4M=t4M;
          nln=nZ;
          nZ=rZ;
          rZ=nln;
          nln=nA;
          nA=rA;
          rA=nln;
        }
       }
       else // Split the cluster (with or w/o "Fermi motion" and "Fermi decay")
       {
        if(cln==1)
        {
          nPDG=90001001;                  // Deuteron
          qM=mDeut;
          nA=2;
          nZ=1;
          restPDG-=1001;
          // Recalculate dmom
          G4ThreeVector sumv=G4ThreeVector(0.,0.,dmom) +
                        fmm*std::pow(-std::log(G4UniformRand()),third)*G4RandomDirection();
          dmom=sumv.mag();
        }
        else if(cln==2)
        {
          nA=3;
          if(G4UniformRand()*(A-2)>(N-1)) // He3
          {
            nPDG=90002001;
            qM=mHel3;
            nZ=2;
            restPDG-=2001;
          }
          else                            // tritium
          {
            nPDG=90001002;
            qM=mTrit;
            nZ=1;
            restPDG-=1002;
          }
          // Recalculate dmom
          G4ThreeVector sumv=G4ThreeVector(0.,0.,dmom) +
                        fmm*std::pow(-std::log(G4UniformRand()),third)*G4RandomDirection()+
                        fmm*std::pow(-std::log(G4UniformRand()),third)*G4RandomDirection();
          dmom=sumv.mag();
        }
        else
        {
          nPDG=90002002;                  // Alpha
          qM=mAlph;
          nA=4;
          nZ=2;
          restPDG-=2002;
          G4ThreeVector sumv=G4ThreeVector(0.,0.,dmom) +
                        fmm*std::pow(-std::log(G4UniformRand()),third)*G4RandomDirection()+
                        fmm*std::pow(-std::log(G4UniformRand()),third)*G4RandomDirection()+
                        fmm*std::pow(-std::log(G4UniformRand()),third)*G4RandomDirection();
          dmom=sumv.mag();
        }
        rA=A-nA;
        rZ=Z-nZ;
        // This is a simple case of cluster at rest
        //G4double rM=G4QPDGCode(restPDG).GetMass();// Mass of the residual nucleus
        //r4M=G4LorentzVector(0.,0.,0.,rM);         // 4mom of the residual nucleus
        //n4M=G4LorentzVector(0.,0.,0.,tM-rM);      // 4mom of the quasi-free cluster
        // --- End of the "simple case of cluster at rest"
        // Make a fake quasi-Fermi distribution for clusters (clusters are not at rest) 
        G4LorentzVector t4M(0.,0.,0.,tM);         // 4m of the target nucleus to be decayed
        G4double rM=G4QPDGCode(restPDG).GetMass();// Mass of the residual nucleus
        r4M=G4LorentzVector(0.,0.,0.,rM);         // 4mom of the residual nucleus
        G4double rM2=rM*rM;
        G4double nM=std::sqrt(rM2+tM*tM-(tM+tM)*std::sqrt(rM2+dmom*dmom));// M of q-cluster
        n4M=G4LorentzVector(0.,0.,0.,nM);         // 4mom of the quasi-nucleon
#ifdef qedebug
        G4cout<<"G4QInel::PStDoIt:QEC, p="<<dmom<<",T="<<tM<<",R="<<rM<<",N="<<nM<<G4endl;
#endif
        if(!G4QHadron(t4M).DecayIn2(r4M, n4M))
        {
          //G4cerr<<"G4QInel::PostStDoIt:M="<<tM<<"<r="<<rM<<"+c="<<nM<<"="<<rM+nM<<G4endl;
          //G4Exception("G4QInelastic::PostStepDoIt()", "HAD_CHPS_0003",
          //           FatalException,"Hadronize quasmon: Can't Dec totNuc->QEClu+ResNuc");
          G4cerr<<"-W-G4QInel::PoStDoI:M="<<tM<<"<rM="<<rM<<"+cM="<<nM<<"="<<rM+nM<<G4endl;
          aParticleChange.ProposeEnergy(kinEnergy);
          aParticleChange.ProposeLocalEnergyDeposit(0.);
          aParticleChange.ProposeMomentumDirection(dir);
          aParticleChange.ProposeTrackStatus(fAlive);
          return G4VDiscreteProcess::PostStepDoIt(track,step);
        }
        // --- End of the moving cluster implementation ---
#ifdef qedebug
        G4cout<<"G4QIn::PStDoIt:QEC, RN="<<r4M.rho()<<r4M<<",QCl="<<n4M.rho()<<n4M<<G4endl;
#endif
       }
       G4LorentzVector s4M=n4M+proj4M;             // Tot 4-momentum for scattering
       G4double prjM2 = proj4M.m2();               // Squared mass of the projectile
       G4double prjM = std::sqrt(prjM2);           // @@ Get from pPDG (?)
       G4double minM = prjM+qM;                    // Min mass sum for the final products
       G4double cmM2 =s4M.m2();
       if(cmM2>minM*minM)
       {
#ifdef qedebug
        G4cout<<"G4QInel::PStDoIt:***Enter***,cmM2="<<cmM2<<" > minM2="<<minM*minM<<G4endl;
#endif
        // Estimate and randomize charge-exchange with a quasi-free cluster
        G4bool chex=false;                        // Flag of the charge exchange scattering
        G4ParticleDefinition* projpt=G4Proton::Proton(); // Prototype, only for chex=true
        // This is reserved for the charge-exchange scattering (to help neutron spectras)
        //if(cln&&!cp1 &&(projPDG==2212&&rA>rZ || projPDG==2112&&rZ>1))// @@ Use proj chex
        if(2>3) // @@ charge exchange is not implemented yet @@
        {
#ifdef qedebug
          G4cout<<"G4QIn::PStDoIt:-Enter, P="<<projPDG<<",cln="<<cln<<",cp1="<<cp1<<G4endl;
#endif
          G4double tprM=mProt;
          G4double tprM2=mProt2;
          G4int tprPDG=2212;
          G4int tresPDG=restPDG+999;
          if(projPDG==2212)
          {
            projpt=G4Neutron::Neutron();
            tprM=mNeut;
            tprM2=mNeut2;
            tprPDG=2112;
            tresPDG=restPDG-999;
          }
          minM=tprM+qM;
          G4double efE=(cmM2-tprM2-qM*qM)/(qM+qM);
          G4double efP=std::sqrt(efE*efE-tprM2);
          G4double chl=qfMan->ChExElCoef(efP*MeV, nZ, nA-nZ, projPDG); // ChEx/Elast(pPDG!)
#ifdef qedebug
          G4cout<<"G4QInel::PStDoIt:chl="<<chl<<",P="<<efP<<",nZ="<<nZ<<",nA="<<nA<<G4endl;
#endif
          if(chl>0.&&cmM2>minM*minM&&G4UniformRand()<chl/(1.+chl))     // minM is redefined
          {
            projPDG=tprPDG;
            prjM=tprM;
            G4double rM=G4QPDGCode(tresPDG).GetMass();// Mass of the residual nucleus
            r4M=G4LorentzVector(0.,0.,0.,rM);         // 4mom of the residual nucleus
            n4M=G4LorentzVector(0.,0.,0.,tM-rM);      // 4mom of the quasi-free cluster
            chex=true;                                // Confirm charge exchange scattering
          }
        }
        //
        std::pair<G4LorentzVector,G4LorentzVector> sctout=qfMan->Scatter(nPDG, n4M,
                                                                         projPDG, proj4M);
#ifdef qedebug
        G4cout<<"G4QInelastic::PStDoIt:QElS,proj="<<prjM<<sctout.second<<",qfCl="<<qM
              <<sctout.first<<",chex="<<chex<<",nA="<<nA<<",nZ="<<nZ<<G4endl;
#endif
        aParticleChange.ProposeLocalEnergyDeposit(0.); // Everything is in particles
        // @@ @@ @@ Coulomb barriers must be checked !! @@ @@ @@ Skip if not
        if(chex)               // ==> Projectile is changed: fill everything to secondaries
        {
          aParticleChange.ProposeEnergy(0.);                // @@ ??
          aParticleChange.ProposeTrackStatus(fStopAndKill); // projectile nucleon is killed
          aParticleChange.SetNumberOfSecondaries(3); 
          G4DynamicParticle* thePrH = new G4DynamicParticle(projpt,sctout.second);
          G4Track* scatPrH = new G4Track(thePrH, localtime, position ); // scattered & chex
          scatPrH->SetWeight(weight);                                   //    weighted
          scatPrH->SetTouchableHandle(trTouchable);                     //   projectile
          aParticleChange.AddSecondary(scatPrH);                        //     hadron
        }
        else               // ==> The leading particle is filled to the updated projectilee
        {
          aParticleChange.SetNumberOfSecondaries(2);        // @@ if proj=leading 
          G4double ldT=(sctout.second).e()-prjM;            // kin Energy of scat project.
          aParticleChange.ProposeEnergy(ldT);               // Change the kin Energy
          G4ThreeVector ldV=(sctout.second).vect();         // Change momentum direction
          aParticleChange.ProposeMomentumDirection(ldV/ldV.mag());
          aParticleChange.ProposeTrackStatus(fAlive);
        }
        // ---------------------------------------------------------
        // Fill scattered quasi-free nucleon/fragment
        if     (nPDG==90000001) theDefinition = G4Neutron::Neutron();
        else if(nPDG==90001000) theDefinition = G4Proton::Proton();
        else if(nZ>0 && nA>1)
                  theDefinition = G4ParticleTable::GetParticleTable()->FindIon(nZ,nA,0,nZ);
#ifdef debug
        else G4cout<<"-Warning_G4QIn::PSD:scatqfPDG="<<nPDG<<", Z="<<nZ<<",A="<<nA<<G4endl;
#endif
        if(nZ>0 && nA>0)
        {
          G4DynamicParticle* theQFN = new G4DynamicParticle(theDefinition,sctout.first);
          G4Track* scatQFN = new G4Track(theQFN, localtime, position ); //   scattered
          scatQFN->SetWeight(weight);                                   //    weighted
          scatQFN->SetTouchableHandle(trTouchable);                     //   quasi-free
          aParticleChange.AddSecondary(scatQFN);                        //  nucleon/cluster
        }
        // ----------------------------------------------------
        // Fill residual nucleus
        if     (restPDG==90000001) theDefinition = G4Neutron::Neutron();
        else if(restPDG==90001000) theDefinition = G4Proton::Proton();
        else if(rZ>0 && rA>1)
                  theDefinition = G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,0,rZ);
#ifdef debug
        else G4cout<<"-Warning_G4QIn::PSD: resPDG="<<restPDG<<",Z="<<rZ<<",A="<<rA<<G4endl;
#endif
        if(rZ>0 && rA>0)
        {
          G4DynamicParticle* theReN = new G4DynamicParticle(theDefinition,r4M);
          G4Track* scatReN = new G4Track(theReN, localtime, position ); //    scattered
          scatReN->SetWeight(weight);                                   //    weighted
          scatReN->SetTouchableHandle(trTouchable);                     //    residual
          aParticleChange.AddSecondary(scatReN);                        //    nucleus
        }
        delete output;
        return G4VDiscreteProcess::PostStepDoIt(track, step);
       }
#ifdef qedebug
       else G4cout<<"G4QInel::PSD:OUT,M2="<<s4M.m2()<<"<"<<minM*minM<<", N="<<nPDG<<G4endl;
#endif
      } // end of proton quasi-elastic (including QE on NF)
      else // @@ make cost-condition @@ Pickup process (pickup 1 or 2 n and make d or t)
      {
       if(projPDG==2212) restPDG--; // Res. nucleus PDG is one neutron less than targ. PDG
       else
       {
         restPDG-=1000;          // Res. nucleus PDG is one proton less than targ. PDG
         rZ--;                   // Reduce the charge of the residual nucleus
       }
       G4double mPUF=mDeut;      // Prototype of the mass of the created pickup fragment
       G4ParticleDefinition* theDefinition = G4Deuteron::Deuteron(); // Default: make d
       //theDefinition = G4ParticleTable::GetParticleTable()->FindIon(nZ,nA,0,nZ);//ion
       G4bool din=false;         // Flag of picking up 2 neutrons to create t (tritium)(p)
       G4bool dip=false;         // Flag of picking up 2 protons to create t (tritium) (n)
       G4bool pin=false;         // Flag of picking up 1 n + 1 p to create He3/t      (p/n)
       G4bool hin=false;         // Flag of PickUp creation of alpha (He4)            (p/n)
       G4double frn=G4UniformRand();
       if(N>2 && frn > 0.95)     // .95 is a parameter to pickup 2N (to cteate t or He3)
       {
         if(projPDG==2212)
         {
           if(N>1 && G4UniformRand()*(Z+.5*(N-1)) > Z)
           {
             theDefinition = G4Triton::Triton();
             mPUF=mTrit;         // The mass of the created pickup fragment is triton
             restPDG--;          // Res. nucleus PDG is two neutrons less than target PDG
             din=true;
           }
           else
           {
             theDefinition = G4He3::He3();
             mPUF=mHel3;         // The mass of the created pickup fragment is Helium3
             restPDG-=1000;      // Res. nucleus PDG is two neutrons less than target PDG
             rZ--;               // Reduce the charge of the residual nucleus
             pin=true;
           }
         }
         else // Neutron projectile (2112)
         {
           if(Z>1 && G4UniformRand()*(N+.5*(Z-1)) > N)
           {
             theDefinition = G4He3::He3();
             mPUF=mHel3;         // The mass of the created pickup fragment is Helium3
             restPDG-=1000;      // Res. nucleus PDG is two neutrons less than target PDG
             rZ--;               // Reduce the charge of the residual nucleus
             dip=true;
           }
           else
           {
             theDefinition = G4Triton::Triton();
             mPUF=mTrit;         // The mass of the created pickup fragment is triton
             restPDG--;          // Res. nucleus PDG is two neutrons less than target PDG
             pin=true;
           }
         }
         rA--;                   // Reduce the baryon number of the residual nucleus
         // Recalculate dmom
         G4ThreeVector sumv=G4ThreeVector(0.,0.,dmom) +
                       fmm*std::pow(-std::log(G4UniformRand()),third)*G4RandomDirection();
         dmom=sumv.mag();
         if(Z>1 && frn > 0.99)  // .99 is a parameter to pickup 2n1p || 1n2p & create alpha
         {
           theDefinition = G4Alpha::Alpha();
           if((din && projPDG==2212) || (pin && projPDG==2112))
           {
             restPDG-=1000;        // Residual nucleus PDG is reduced to create alpha
             rZ--;                 // Reduce the charge of the residual nucleus
           }
           else if((pin && projPDG==2212) || (dip && projPDG==2112)) restPDG--;
           else G4cout<<"-Warning-G4QIn::PSD: PickUp logic error, proj="<<projPDG<<G4endl;
           hin=true;
           mPUF=mAlph;           // The mass of the created pickup fragment (alpha)
           rA--;                 // Reduce the baryon number of the residual nucleus
           // Recalculate dmom
           sumv += fmm*std::pow(-std::log(G4UniformRand()),third)*G4RandomDirection();
           dmom=sumv.mag();
         }
       }
       G4double mPUF2=mPUF*mPUF;
       G4LorentzVector t4M(0.,0.,0.,tM);     // 4m of the target nucleus to be decayed
       G4double rM=G4QPDGCode(restPDG).GetMass();// Mass of the residual nucleus
       G4double rM2=rM*rM;                   // Squared mass of the residual nucleus
       G4double nM2=rM2+tM*tM-(tM+tM)*std::sqrt(rM2+dmom*dmom);// M2 of boundQF-nucleon(2N)
       if(nM2 < 0) nM2=0.;
       G4double den2=(dmom*dmom+nM2);        // squared energy of bQF-nucleon
       G4double den=std::sqrt(den2);         // energy of bQF-nucleon
#ifdef qedebug
       G4double nM=std::sqrt(nM2);           // M of bQF-nucleon
       G4cout<<"G4QInel::PStDoIt:PiU, p="<<dmom<<",tM="<<tM<<", R="<<rM<<",N="<<nM<<G4endl;
#endif
       G4double qp=momentum*dmom;
       G4double srm=proj4M.e()*den;
       G4double mNucl2=mProt2;
       if(projPDG == 2112) mNucl2=mNeut2;
       G4double cost=(nM2+mNucl2+srm+srm-mPUF2)/(qp+qp);
       G4int ict=0;
       while(std::fabs(cost)>1. && ict<3)
       {
         dmom=91.;                           // Fermi momentum (proDefault for a deuteron)
         if(Z>1||N>1) dmom=fmm*std::pow(-std::log(G4UniformRand()),third);
         if(din || pin || dip)               // Recalculate dmom for the final t/He3
         {
          G4ThreeVector sumv=G4ThreeVector(0.,0.,dmom) +
                        fmm*std::pow(-std::log(G4UniformRand()),third)*G4RandomDirection();
          if(hin)
                  sumv+=fmm*std::pow(-std::log(G4UniformRand()),third)*G4RandomDirection();
          dmom=sumv.mag();
         }
         nM2=rM2+tM*tM-(tM+tM)*std::sqrt(rM2+dmom*dmom); // M2 of bQF-nucleon
         if(nM2 < 0) nM2=0.;
         //nM=std::sqrt(nM2);                  // M of bQF-nucleon --Not used?--
         den2=(dmom*dmom+nM2);               // squared energy of bQF-nucleon
         den=std::sqrt(den2);                // energy of bQF-nucleon
         qp=momentum*dmom;
         srm=proj4M.e()*den;
         cost=(nM2+mNucl2+srm+srm-mPUF2)/(qp+qp);
         ict++;
#ifdef ldebug
         if(ict>2)G4cout<<"G4QInelast::PStDoIt:i="<<ict<<",d="<<dmom<<",ct="<<cost<<G4endl;
#endif
       }
       if(std::fabs(cost)<=1.)
       {// ---- @@ From here can be a MF for QF nucleon extraction (if used by others)
        G4ThreeVector vdir = proj4M.vect();  // 3-Vector in the projectile direction
        G4ThreeVector vx(0.,0.,1.);          // ProtoOrt in theDirection of the projectile
        G4ThreeVector vy(0.,1.,0.);          // First ProtoOrt orthogonal to the direction
        G4ThreeVector vz(1.,0.,0.);          // SecondProtoOrt orthoganal to the direction
        if(vdir.mag2() > 0.)                 // the projectile isn't at rest
        {
          vx = vdir.unit();                  // Ort in the direction of the projectile
          G4ThreeVector vv= vx.orthogonal(); // Not normed orthogonal vector (!)
          vy = vv.unit();                    // First ort orthogonal to the proj direction
          vz = vx.cross(vy);                 // Second ort orthoganal to the projDirection
        }
        // ---- @@ End of possible MF (Similar is in G4QEnvironment)
        G4double tdom=dmom*std::sqrt(1-cost*cost);// Transversal component of the Fermi-Mom
        G4double phi=twopi*G4UniformRand();  // Phi of the Fermi-Mom
        G4ThreeVector vedm=vx*dmom*cost+vy*tdom*std::sin(phi)+vz*tdom*std::cos(phi);// bQFN
        G4LorentzVector bqf(vedm,den);      // 4-mom of the bQF nucleon
        r4M=t4M-bqf;                         // 4mom of the residual nucleus
#ifdef debug
        if(std::fabs(r4M.m()-rM)>.001) G4cout<<"G4QIn::PSD: rM="<<rM<<"#"<<r4M.m()<<G4endl;
#endif
        G4LorentzVector f4M=proj4M+bqf;      // Prototype of 4-mom of Deuteron
#ifdef debug
        if(std::fabs(f4M.m()-mPUF)>.001)G4cout<<"G4QI::PSD:m="<<mPUF<<"#"<<f4M.m()<<G4endl;
#endif
        aParticleChange.ProposeEnergy(0.);   // @@ ??
        aParticleChange.ProposeTrackStatus(fStopAndKill);// projectile nucleon is killed
        aParticleChange.SetNumberOfSecondaries(2); 
        // Fill pickup hadron
        G4DynamicParticle* theQFN = new G4DynamicParticle(theDefinition,f4M);
        G4Track* scatQFN = new G4Track(theQFN, localtime, position ); //    pickup
        scatQFN->SetWeight(weight);                                   //    weighted
        scatQFN->SetTouchableHandle(trTouchable);                     //    nuclear
        aParticleChange.AddSecondary(scatQFN);                        //    cluster
#ifdef pickupd
        G4cout<<"G4QInelastic::PStDoIt: f="<<theDefinition<<",f4M="<<f4M.m()<<f4M<<G4endl;
#endif
        // ----------------------------------------------------
        // Fill residual nucleus
        if     (restPDG==90000001) theDefinition = G4Neutron::Neutron();
        else if(restPDG==90001000) theDefinition = G4Proton::Proton();
        else theDefinition = G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,0,rZ);//ion
        G4DynamicParticle* theReN = new G4DynamicParticle(theDefinition,r4M);
        G4Track* scatReN = new G4Track(theReN, localtime, position ); //    scattered
        scatReN->SetWeight(weight);                                   //    weighted
        scatReN->SetTouchableHandle(trTouchable);                     //    residual
        aParticleChange.AddSecondary(scatReN);                        //    nucleus
#ifdef pickupd
        G4cout<<"G4QInelastic::PStDoIt:rZ="<<rZ<<",rA="<<rA<<",r4M="<<r4M.m()<<r4M<<G4endl;
#endif
        delete output;
        return G4VDiscreteProcess::PostStepDoIt(track, step);
	 }
      }
    } // end of the quasi-elastic decoupled process for the neucleon projectile
  }  // end lepto-nuclear, neutrino-nuclear, proton/neutron decoupled processes
  EnMomConservation=proj4M+G4LorentzVector(0.,0.,0.,tM);    // Total 4-mom of the reaction
  if(absMom) EnMomConservation+=lead4M;         // Add E/M of leading System
#ifdef debug
  G4cout<<"G4QInelastic::PostStDoIt:before St="<<aParticleChange.GetTrackStatus()<<G4endl;
#endif

  // @@@@@@@@@@@@@@ Temporary for the testing purposes --- Begin
  //G4bool elF=false; // Flag of the ellastic scattering is "false" by default
  //G4double eWei=1.;
  // @@@@@@@@@@@@@@ Temporary for the testing purposes --- End  
#ifdef debug
  G4cout<<"^^G4QInelastic::PostStepDoIt: projPDG="<<projPDG<<", targPDG="<<targPDG<<G4endl;
#endif
  const G4QNucleus targNuc(Z,N);                          // Define the target nucleus
  if(projPDG>9000000)
  {
    delete output;                                        // Before the output creation
    G4QNucleus projNuc(proj4M,projPDG);                   // Define the projectile nucleus
    G4QIonIonCollision IonIon(projNuc, targNuc);  // Define deep-inelastic ion-ion reaction
#ifdef debug
    G4cout<<"G4QInel::PStDoIt: ProjNuc="<<projPDG<<proj4M<<", TargNuc="<<targPDG<<G4endl;
#endif
    try {output = IonIon.Fragment();}                     // DINR AA interaction products
    catch (G4QException& error)
    {
      G4cerr<<"***G4QInelastic::PostStepDoIt: G4QE Exception is catched in hA"<<G4endl;
      // G4Exception("G4QInelastic::PostStepDoIt:","27",FatalException,"CHIPS hA crash");
      G4Exception("G4QInelastic::PostStepDoIt()", "HAD_CHPS_0027",
                  FatalException, "CHIPS hA crash");
    }
  }
  else                                                    // --> A projectile hadron
  {
    // *** Let to produce elastic final states for acceleration ***
#ifdef debug
    G4int maxCn=7;
#endif
    G4int atCn=0;                                         // Attempts counter
    //G4bool inel=false;
    //while (inel==false && ++atCn <= maxCn)                // To evoid elastic final state
    //{
#ifdef debug
      G4cout<<"G4QInelastic::PStDoIt: atCn="<<atCn<<" <= maxCn="<<maxCn<<G4endl;
#endif
      G4int outN=output->size();
      if(outN)                                            // The output is not empty
      {
        std::for_each(output->begin(), output->end(), DeleteQHadron());
        output->clear();
      }
      delete output;                                      // Before the new output creation
      const G4QHadron projH(projPDG,proj4M);              // Define the projectile particle
#ifdef debug
      G4cout<<"G4QInelastic::PStDoIt: Before Fragmentation"<<G4endl;
#endif
      //if(aProjPDG==22 && projPDG==22 && proj4M.m2()<.01 && proj4M.e()<150.)
      //{
      //  G4QHadron* tgH= new G4QHadron(targNuc.GetPDGCode(),targNuc.Get4Momentum());
      //  G4QHadron* prH= new G4QHadron(projPDG,proj4M);
      //  output = new G4QHadronVector();
      //  output->push_back(prH);
      //  output->push_back(tgH);
      //}
      //else
      //{
        G4QFragmentation DINR(targNuc, projH);            // Define deep-inel. h-A reaction
#ifdef debug
        G4cout<<"G4QInelastic::PStDoIt:Proj="<<projPDG<<proj4M<<",TgNuc="<<targPDG<<G4endl;
#endif
        try {output = DINR.Fragment();}                   // DINR hA fragmentation
        catch (G4QException& error)
        {
          G4cerr<<"***G4QInelastic::PostStepDoIt: G4QE Exception is catched in hA"<<G4endl;
          // G4Exception("G4QInelastic::PostStepDoIt:","27",FatalException,"CHIPS hA crash");
          G4Exception("G4QInelastic::PostStepDoIt()", "HAD_CHPS_0027",
                      FatalException, "CHIPS hA crash");
        }
	//}
      outN=output->size();
#ifdef debug
      G4cout<<"G4QInelastic::PostStepDoIt:  At# "<<atCn<<", nSec="<<outN<<", fPDG="
            <<(*output)[0]->GetPDGCode()<<", pPDG="<< projPDG<<G4endl;
#endif
      //inel=true;                                        // For while
      if(outN < 2)
      {
        G4cout<<"-Warning-G4QInelastic::PostStepDoIt: nSec="<<outN<<", At# "<<atCn<<G4endl;
        //inel=false;                                     // For while
      }
      //else if(outN==2 && (*output)[0]->GetPDGCode() == projPDG) inel=false; // For while
#ifdef debug
      if(atCn==maxCn)G4cout<<"-Warning-G4QI::PostStDoIt:mAt="<<atCn<<" is reached"<<G4endl;
#endif
    //}
  }
  //
  // --- the scattered hadron with changed nature can be added here ---
  if(scat)
  {
    G4QHadron* scatHadron = new G4QHadron(scatPDG,scat4M);
    output->push_back(scatHadron);
  }
  G4int qNH=0;
  if(leadhs) qNH=leadhs->size();
  if(absMom)
  {
    if(qNH) for(G4int iq=0; iq<qNH; iq++)
    {
      G4QHadron* loh=(*leadhs)[iq];   // Pointer to the output hadron
      output->push_back(loh);
    }
    if(leadhs) delete leadhs;
    leadhs=0;
  }
  else
  {
    if(qNH) for(G4int iq=0; iq<qNH; iq++) delete (*leadhs)[iq];
    if(leadhs) delete leadhs;
    leadhs=0;
  }
  // ------------- From here the secondaries are filled -------------------------
  G4int tNH = output->size();       // A#of hadrons in the output
  aParticleChange.SetNumberOfSecondaries(tNH); // @@ tNH can be changed (drop/decay!)
  // Now add nuclear fragments
#ifdef debug
  G4cout<<"G4QInelastic::PostStepDoIt: "<<tNH<<" particles are generated"<<G4endl;
#endif
#ifdef ppdebug
  if(absMom)G4cout<<"G4QInelastic::PostStepDoIt: t="<<tNH<<", q="<<qNH<<G4endl;
#endif
  G4int nOut=output->size();               // Real length of the output @@ Temporary
  if(tNH==1 && !scat)                      // @@ Temporary. Find out why it happened!
  {
    G4cout<<"-Warning-G4QInelastic::PostStepDoIt: 1 secondary! absMom="<<absMom;
    if(absMom) G4cout<<", qNH="<<qNH;
    G4cout<<", PDG0="<<(*output)[0]->GetPDGCode();
    G4cout<<G4endl;
    tNH=0;
    delete output->operator[](0);          // delete the creazy hadron
    output->pop_back();                    // clean up the output vector
  }
  if(tNH==2&&2!=nOut) G4cout<<"--Warning--G4QInelastic::PostStepDoIt: 2 # "<<nOut<<G4endl;
  // Deal with ParticleChange final state interface to GEANT4 output of the process
  //if(tNH==2) for(i=0; i<tNH; i++)          // @@ Temporary tNH==2 instead of just tNH
  if(tNH) for(i=0; i<tNH; i++)             // @@ Temporary tNH==2 instead of just tNH
  {
    // Note that one still has to take care of Hypernuclei (with Lambda or Sigma inside)
    // Hypernucleus mass calculation and ion-table interface upgrade => work for Hisaya @@
    // The decau process for hypernuclei must be developed in GEANT4 (change CHIPS body)
    G4QHadron* hadr=(*output)[i];          // Pointer to the output hadron    
    G4int PDGCode = hadr->GetPDGCode();
    G4int nFrag   = hadr->GetNFragments();
#ifdef pdebug
    G4cout<<"G4QInelastic::PostStepDoIt: H#"<<i<<",PDG="<<PDGCode<<",nF="<<nFrag
          <<", 4Mom="<<hadr->Get4Momentum()<<G4endl;
#endif
    if(nFrag)                              // Skip intermediate (decayed) hadrons
    {
#ifdef debug
      G4cout<<"G4QInelastic::PostStepDoIt: Intermediate particle is found i="<<i<<G4endl;
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
    else if(PDGCode==90000000)
    {
#ifdef pdebug
      G4cout<<"-Warning-G4QInelastic::PostStepDoIt:Vacuum PDG="<<PDGCode
            <<", 4Mom="<<hadr->Get4Momentum()<<G4endl;
#endif
      theDefinition = G4Gamma::Gamma();
      G4double hM2=(hadr->Get4Momentum()).m2();
      if(std::fabs(hM2)>.01) G4cout<<"-Warning-G4QInel::PoStDoIt:90000000M2="<<hM2<<G4endl;
      else if(hadr->Get4Momentum()==vacuum4M)
      {
        delete hadr;
        continue;
      }
    }
    else if(PDGCode >80000000) // Defines hypernuclei as normal nuclei (N=N+S Correction!)
    {
      G4int aZ = hadr->GetCharge();
      G4int aA = hadr->GetBaryonNumber();
#ifdef pdebug
      G4cout<<"G4QInelastic::PostStepDoIt:Ion Z="<<aZ<<", A="<<aA<<G4endl;
#endif
      theDefinition = G4ParticleTable::GetParticleTable()->FindIon(aZ,aA,0,aZ);
    }
    //else theDefinition = G4ParticleTable::GetParticleTable()->FindParticle(PDGCode);
    else if(PDGCode==89999998 || PDGCode==89998000 || PDGCode==88000000) // anti-di-baryon
    {
      G4double rM=mNeut;                                // Prototype of the baryon mass
      G4int  rPDG=-2112;                                // Prototype of the baryon PDG
      G4ParticleDefinition* BarDef=G4Neutron::Neutron();// Prototype of theBaryonDefinition
      if     (PDGCode==89998000)
      {
        rM=mProt;
        rPDG=-2212;
        BarDef=G4Proton::Proton();
      }
      else if(PDGCode==88000000)
      {
        rM=mLamb;
        rPDG=-3122;
        BarDef=G4Lambda::Lambda();
      }
      G4LorentzVector t4M=hadr->Get4Momentum();         // 4m of the di-baryon
      G4LorentzVector f4M=G4LorentzVector(0.,0.,0.,rM); // 4mom of the 1st anti-baryon
      G4LorentzVector s4M=G4LorentzVector(0.,0.,0.,rM); // 4mom of the 2nd anti-baryon
#ifdef qedebug
      G4cout<<"G4QInel::PStDoIt:AntiDiBar,t4M="<<tM<<",m="<<rM<<",PDG="<<PDGCode<<G4endl;
#endif
      if(!G4QHadron(t4M).DecayIn2(f4M, s4M))
      {
        G4cerr<<"G4QIn::PostStDoIt: ADB, M="<<t4M.m()<<" < 2*rM="<<rM<<" = "<<2*rM<<G4endl;
        // throw G4QException("G4QInelastic::HadronizeQuasm:Can't decay anti-dibaryon");
        G4Exception("G4QInelastic::PostStepDoIt()", "HAD_CHPS_0004",
                    FatalException, "Hadronize quasmon: Can't decay anti-dibaryon");
      }
      // --- End of the moving cluster implementation ---
#ifdef qedebug
      G4cout<<"G4QInelastic::PStDoIt: ADB, H1="<<rM<<f4M<<", H2="<<rM<<s4M<<G4endl;
#endif
      G4QHadron fH(rPDG,f4M);
      hadr->Set4Momentum(f4M);
      hadr->SetQPDG(fH.GetQPDG());
      theDefinition=BarDef;
#ifdef debug
      G4cout<<"G4QInel::PostStDoIt:Anti-DiBar, DecayIn2, h1="<<rPDG<<f4M<<G4endl;
#endif
      G4QHadron* sH = new G4QHadron(rPDG,s4M);
      output->push_back(sH);
      ++tNH;
#ifdef debug
      G4cout<<"G4QInel::PostStDoIt:Anti-DiBar, DecayIn2, h2="<<rPDG<<s4M<<G4endl;
#endif
    }
    else if(PDGCode==90000997 || PDGCode==89997001) // anti-NDelta
    {
      G4double rM=mNeut;                                // Prototype of the baryon mass
      G4int  rPDG=-2112;                                // Prototype of the baryon PDG
      G4double iM=mPi;                                  // Prototype of the pion mass
      G4int  iPDG= 211;                                 // Prototype of the pion PDG
      G4ParticleDefinition* BarDef=G4Neutron::Neutron();// Prototype of theBaryonDefinition
      if(PDGCode==90000997)
      {
        rM=mProt;
        rPDG=-2212;
        iPDG=-211;
        BarDef=G4Proton::Proton();
      }
      G4LorentzVector t4M=hadr->Get4Momentum();         // 4m of the di-baryon
      G4LorentzVector f4M=G4LorentzVector(0.,0.,0.,rM); // 4mom of the 1st anti-baryon
      G4LorentzVector s4M=G4LorentzVector(0.,0.,0.,rM); // 4mom of the 2nd anti-baryon
      G4LorentzVector u4M=G4LorentzVector(0.,0.,0.,iM); // 4mom of the pion
#ifdef qedebug
      G4cout<<"G4QInel::PStDoIt:AntiNDelta, t4M="<<tM<<",m="<<rM<<",PDG="<<PDGCode<<G4endl;
#endif
      if(!G4QHadron(t4M).DecayIn3(f4M, s4M, u4M))
      {
        G4cerr<<"G4QIn::PostStDoIt: AND, tM="<<t4M.m()<<" < 2*mB+mPi="<<2*rM+iM<<G4endl;
        // throw G4QException("G4QInelastic::HadronizeQuasm:Can't decay anti-NDelta");
        G4Exception("G4QInelastic::PostStepDoIt()", "HAD_CHPS_0005",
                    FatalException, "Hadronize quasmon: Can't decay anti-NDelta");
      }
      // --- End of the moving cluster implementation ---
#ifdef qedebug
      G4cout<<"G4QInel::PStDoI:AND,B1="<<rM<<f4M<<",B2="<<rM<<s4M<<",Pi="<<iM<<u4M<<G4endl;
#endif
      G4QHadron fH(rPDG,f4M);
      hadr->Set4Momentum(f4M);
      hadr->SetQPDG(fH.GetQPDG());
      theDefinition=BarDef;
#ifdef debug
      G4cout<<"G4QInel::PostStDoIt:Anti-NDelta, DecayIn2, h1="<<rPDG<<f4M<<G4endl;
#endif
      G4QHadron* sH = new G4QHadron(rPDG,s4M);
      output->push_back(sH);
      ++tNH;
#ifdef debug
      G4cout<<"G4QInel::PostStDoIt:Anti-NDelta, DecayIn2, h2="<<rPDG<<s4M<<G4endl;
#endif
      G4QHadron* uH = new G4QHadron(iPDG,u4M);
      output->push_back(uH);
      ++tNH;
#ifdef debug
      G4cout<<"G4QInel::PostStDoIt:Anti-NDelta, DecayIn2, h2="<<rPDG<<s4M<<G4endl;
#endif
    }
    else
    {
#ifdef pdebug
      G4cout<<"G4QInelastic::PostStepDoIt:Define particle with PDG="<<PDGCode<<G4endl;
#endif
      theDefinition = G4QPDGToG4Particle::Get()->GetParticleDefinition(PDGCode);
#ifdef pdebug
      G4cout<<"G4QInelastic::PostStepDoIt:AfterParticleDefinition PDG="<<PDGCode<<G4endl;
#endif
    }
    if(!theDefinition)
    { // !! Commenting of this print costs days of searching for E/mom nonconservation !!
      if(PDGCode!=90000000 || hadr->Get4Momentum()!=vacuum4M)
        G4cout<<"---Warning---G4QInelastic::PostStepDoIt: drop PDG="<<PDGCode<<", 4Mom="
              <<hadr->Get4Momentum()<<G4endl;
      delete hadr;
      continue;
    }
#ifdef pdebug
    G4cout<<"G4QInelastic::PostStepDoIt:Name="<<theDefinition->GetParticleName()<<G4endl;
#endif
    theSec->SetDefinition(theDefinition);
    G4LorentzVector h4M=hadr->Get4Momentum();
    EnMomConservation-=h4M;
#ifdef tdebug
    G4cout<<"G4QInelast::PSDI: "<<i<<","<<PDGCode<<h4M<<h4M.m()<<EnMomConservation<<G4endl;
#endif
#ifdef debug
    G4cout<<"G4QInelastic::PostStepDoIt:#"<<i<<",PDG="<<PDGCode<<",4M="<<h4M<<G4endl;
#endif
    theSec->Set4Momentum(h4M); //                                                         ^
    delete hadr; // <-----<-----------<-------------<---------------------<---------<-----+
#ifdef debug
    G4ThreeVector curD=theSec->GetMomentumDirection();              //                    ^
    G4double curM=theSec->GetMass();                                //                    |
    G4double curE=theSec->GetKineticEnergy()+curM;                  //                    ^
    G4cout<<"G4QInelast::PSDoIt:p="<<curD<<curD.mag()<<",e="<<curE<<",m="<<curM<<G4endl;//|
#endif
    G4Track* aNewTrack = new G4Track(theSec, localtime, position ); //                    ^
    aNewTrack->SetWeight(weight);                                   //    weighted        |
    aNewTrack->SetTouchableHandle(trTouchable);                     //                    |
    aParticleChange.AddSecondary( aNewTrack );                      //                    |
#ifdef debug
    G4cout<<"G4QInelastic::PostStepDoIt:#"<<i<<" is done."<<G4endl; //                    |
#endif
  } //                                                                                    |
  delete output; // instances of the G4QHadrons from the output are already deleted above +
  if(leadhs)     // To satisfy Valgrind ( How can that be?)
  {
    qNH=leadhs->size();
    if(qNH) for(G4int iq=0; iq<qNH; iq++) delete (*leadhs)[iq];
    delete leadhs;
    leadhs=0;
  }
#ifdef debug
  G4cout<<"G4QInelastic::PostStDoIt: after St="<<aParticleChange.GetTrackStatus()<<G4endl;
#endif
  if(aProjPDG!=11 && aProjPDG!=13 && aProjPDG!=15)
    aParticleChange.ProposeTrackStatus(fStopAndKill);        // Kill the absorbed particle
#ifdef pdebug
    G4cout<<"G4QInelastic::PostStepDoIt: E="<<aParticleChange.GetEnergy()
          <<", d="<<*aParticleChange.GetMomentumDirection()<<G4endl;
#endif
#ifdef ldebug
  G4cout<<"G4QInelastic::PostStepDoIt:*** PostStepDoIt is done ***, P="<<aProjPDG<<", St="
        <<aParticleChange.GetTrackStatus()<<G4endl;
#endif
  return G4VDiscreteProcess::PostStepDoIt(track, step);
}

std::pair<G4double,G4double> G4QInelastic::Random2DDirection()
{
  G4double sp=0;              // sin(phi)
  G4double cp=1.;             // cos(phi)
  G4double r2=2.;             // to enter the loop
  while(r2>1. || r2<.0001)    // pi/4 efficiency
  {
    G4double sine=G4UniformRand();
    G4double cosine=G4UniformRand();
    sp=1.-sine-sine;
    cp=1.-cosine-cosine;
    r2=sp*sp+cp*cp;
  }
  G4double norm=std::sqrt(r2);
  return std::make_pair(sp/norm,cp/norm);
}
