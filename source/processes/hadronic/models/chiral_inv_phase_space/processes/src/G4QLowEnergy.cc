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
// $Id: G4QLowEnergy.cc,v 1.5 2010/06/14 16:11:27 mkossov Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
//      ---------------- G4QLowEnergy class -----------------
//                 by Mikhail Kossov, Aug 2007.
// G4QLowEnergy class of the CHIPS Simulation Branch in GEANT4
// ---------------------------------------------------------------
// Short description: This is a fast low energy algorithm for the
// inelastic interactions of nucleons and nuclei (ions) with nuclei.
// This is a fase-space algorithm, but not quark level. Provides
// nuclear fragments upto alpha only. Never was tumed (but can be).
// ---------------------------------------------------------------
// Andrea Dotti (andrea.dotti@cern.ch) 19 January 2010:
//		-	Fix for crash seen on p@450GeV/c + Be interaction
//			due to a fragment of the target composed of two neutrons.
//			The fixt consists of copying the (missing) handling of
//			special cases (Z=0,A>=2) when doing the "splitting" of the
//			target. Code copied from the section dealing with the
//			"splitting" of the projectile. See comments at line 1197
//			Code added is between comments //ANDREA.
//		-	raise a G4QException when theDefinition is a null-pointer
//			(before was only issuing a WARNING and continuing)
//#define debug
//#define pdebug
//#define edebug
//#define tdebug
//#define nandebug
//#define ppdebug

#include "G4QLowEnergy.hh"

// Initialization of static vectors
G4int    G4QLowEnergy::nPartCWorld=152;  // #of particles initialized in CHIPS
std::vector<G4int> G4QLowEnergy::ElementZ; // Z of element(i) in theLastCalc
std::vector<G4double> G4QLowEnergy::ElProbInMat; // SumProbOfElem in Material
std::vector<std::vector<G4int>*> G4QLowEnergy::ElIsoN;// N of isotope(j), E(i)
std::vector<std::vector<G4double>*>G4QLowEnergy::IsoProbInEl;//SumProbIsotE(i)

// Constructor
G4QLowEnergy::G4QLowEnergy(const G4String& processName):
  G4VDiscreteProcess(processName, fHadronic), evaporate(true)
{
#ifdef debug
  G4cout<<"G4QLowEnergy::Constructor is called processName="<<processName<<G4endl;
#endif
  if (verboseLevel>0) G4cout<<GetProcessName()<<" process is created "<<G4endl;
  G4QCHIPSWorld::Get()->GetParticles(nPartCWorld); // Create CHIPS World (234 part. max)
}

// Destructor
G4QLowEnergy::~G4QLowEnergy() {}


G4LorentzVector G4QLowEnergy::GetEnegryMomentumConservation(){return EnMomConservation;}

G4int G4QLowEnergy::GetNumberOfNeutronsInTarget() {return nOfNeutrons;}

// output of the function must be in units of length! L=1/sig_V,sig_V=SUM(n(j,i)*sig(j,i)),
// where n(i,j) is a number of nuclei of the isotop j of the element i in V=1(lengtUnit^3)
// ********** All CHIPS cross sections are calculated in the surface units ************
G4double G4QLowEnergy::GetMeanFreePath(const G4Track&Track, G4double, G4ForceCondition*F)
{
  *F = NotForced;
  const G4DynamicParticle* incidentParticle = Track.GetDynamicParticle();
  G4ParticleDefinition* incidentParticleDefinition=incidentParticle->GetDefinition();
  if( !IsApplicable(*incidentParticleDefinition))
    G4cout<<"-Warning-G4QLowEnergy::GetMeanFreePath for notImplemented Particle"<<G4endl;
  // Calculate the mean Cross Section for the set of Elements(*Isotopes) in the Material
  G4double Momentum = incidentParticle->GetTotalMomentum(); // 3-momentum of the Particle
#ifdef debug
  G4double KinEn = incidentParticle->GetKineticEnergy();
  G4cout<<"G4QLowEnergy::GetMeanFreePath:Prpj, kinE="<<KinEn<<", Mom="<<Momentum<<G4endl;
#endif
  const G4Material* material = Track.GetMaterial();        // Get the current material
  const G4double* NOfNucPerVolume = material->GetVecNbOfAtomsPerVolume();
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QLowEnergy::GetMeanFreePath:"<<nE<<" Elems"<<G4endl;
#endif
  G4int pPDG=0;
  G4int Z=static_cast<G4int>(incidentParticleDefinition->GetPDGCharge());
  G4int A=incidentParticleDefinition->GetBaryonNumber();
  if      ( incidentParticleDefinition == G4Proton::Proton()     ) pPDG = 2212;
  else if ( incidentParticleDefinition == G4Deuteron::Deuteron() ) pPDG = 1000010020;
  else if ( incidentParticleDefinition == G4Alpha::Alpha()       ) pPDG = 1000020040;
  else if ( incidentParticleDefinition == G4Triton::Triton()     ) pPDG = 1000010030;
  else if ( incidentParticleDefinition == G4He3::He3()           ) pPDG = 1000020030;
  else if ( incidentParticleDefinition == G4GenericIon::GenericIon() || (Z > 0 && A > 0))
  {
    pPDG=incidentParticleDefinition->GetPDGEncoding();
#ifdef debug
    G4int prPDG=1000000000+10000*A+10*Z;
    G4cout<<"G4QIonIonElastic::GetMeanFreePath: PDG="<<prPDG<<"="<<pPDG<<G4endl;
#endif
  }
  else G4cout<<"-Warning-G4QLowEnergy::GetMeanFreePath: only AA & pA implemented"<<G4endl;
  G4VQCrossSection* CSmanager=G4QIonIonCrossSection::GetPointer();
  if(pPDG == 2212) CSmanager=G4QProtonNuclearCrossSection::GetPointer(); // just to test
  Momentum/=incidentParticleDefinition->GetBaryonNumber(); // Divide Mom by projectile A
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
    G4int indEl=0;                          // Index of non-natural element or 0(default)
    G4IsotopeVector* isoVector=pElement->GetIsotopeVector(); // Get the predefined IsoVect
    if(isoVector) isoSize=isoVector->size();// Get size of the existing isotopeVector
#ifdef debug
    G4cout<<"G4QLowEnergy::GetMeanFreePath: isovector Length="<<isoSize<<G4endl;
#endif
    if(isoSize)                             // The Element has non-trivial abundance set
    {
      indEl=pElement->GetIndex()+1;         // Index of the non-trivial element is an order
#ifdef debug
      G4cout<<"G4QLowEn::GetMFP:iE="<<indEl<<",def="<<Isotopes->IsDefined(Z,indEl)<<G4endl;
#endif
      if(!Isotopes->IsDefined(Z,indEl))     // This index is not defined for this Z: define
      {
        std::vector<std::pair<G4int,G4double>*>* newAbund =
                                               new std::vector<std::pair<G4int,G4double>*>;
        G4double* abuVector=pElement->GetRelativeAbundanceVector();
        for(G4int j=0; j<isoSize; j++)      // Calculation of abundance vector for isotopes
        {
          G4int N=pElement->GetIsotope(j)->GetN()-Z; // N means A=N+Z !
          if(pElement->GetIsotope(j)->GetZ()!=Z)G4cerr<<"G4QDiffract::GetMeanFreePath: Z="
                                         <<pElement->GetIsotope(j)->GetZ()<<"#"<<Z<<G4endl;
          G4double abund=abuVector[j];
          std::pair<G4int,G4double>* pr= new std::pair<G4int,G4double>(N,abund);
#ifdef debug
          G4cout<<"G4QLowEnergy::GetMeanFreePath:pair#"<<j<<",N="<<N<<",a="<<abund<<G4endl;
#endif
          newAbund->push_back(pr);
        }
#ifdef debug
        G4cout<<"G4QLowEnergy::GetMeanFreePath: pairVectLength="<<newAbund->size()<<G4endl;
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
#ifdef debug
    G4cout<<"G4QLowEnergy::GetMFP:***=>,#isot="<<nIs<<", Z="<<Z<<", indEl="<<indEl<<G4endl;
#endif
    G4double susi=0.;                       // sum of CS over isotopes
    if(nIs) for(G4int j=0; j<nIs; j++)      // Calculate CS for eachIsotope of El
    {
      std::pair<G4int,G4double>* curIs=(*cs)[j]; // A pointer, which is used twice
      G4int N=curIs->first;                 // #of Neuterons in the isotope j of El i
      IsN->push_back(N);                    // Remember Min N for the Element
#ifdef debug
      G4cout<<"G4QLoE::GMFP:true,P="<<Momentum<<",Z="<<Z<<",N="<<N<<",pPDG="<<pPDG<<G4endl;
#endif
      G4bool ccsf=true;                    // Extract inelastic Ion-Ion cross-section
#ifdef debug
      G4cout<<"G4QLowEnergy::GMFP: GetCS #1 j="<<j<<G4endl;
#endif
      G4double CSI=CSmanager->GetCrossSection(ccsf,Momentum,Z,N,pPDG);//CS(j,i) for isotope
#ifdef debug
      G4cout<<"G4QLowEnergy::GetMeanFreePath: jI="<<j<<", Zt="<<Z<<", Nt="<<N<<", Mom="
            <<Momentum<<", XSec="<<CSI/millibarn<<G4endl;
#endif
      curIs->second = CSI;
      susi+=CSI;                            // Make a sum per isotopes
      SPI->push_back(susi);                 // Remember summed cross-section
    } // End of temporary initialization of the cross sections in the G4QIsotope singeltone
    sigma+=Isotopes->GetMeanCrossSection(Z,indEl)*NOfNucPerVolume[i];//SUM(MeanCS*NOfNperV)
#ifdef debug
    G4cout<<"G4QLowEnergy::GetMeanFreePath:<XS>="<<Isotopes->GetMeanCrossSection(Z,indEl)
          <<",AddSigm="<<Isotopes->GetMeanCrossSection(Z,indEl)*NOfNucPerVolume[i]<<G4endl;
#endif
    ElProbInMat.push_back(sigma);
  } // End of LOOP over Elements
  // Check that cross section is not zero and return the mean free path
#ifdef debug
  G4cout<<"G4QLowEnergy::GetMeanFreePath: MeanFreePath="<<1./sigma<<G4endl;
#endif
  if(sigma > 0.000000001) return 1./sigma;                 // Mean path [distance] 
  return DBL_MAX;
}

G4bool G4QLowEnergy::IsApplicable(const G4ParticleDefinition& particle) 
{
  G4int Z=static_cast<G4int>(particle.GetPDGCharge());
  G4int A=particle.GetBaryonNumber();
  if      (particle == *(     G4Proton::Proton()     )) return true;
  else if (particle == *(    G4Neutron::Neutron()    )) return true;
  else if (particle == *(   G4Deuteron::Deuteron()   )) return true;
  else if (particle == *(      G4Alpha::Alpha()      )) return true;
  else if (particle == *(     G4Triton::Triton()     )) return true;
  else if (particle == *(        G4He3::He3()        )) return true;
  else if (particle == *( G4GenericIon::GenericIon() )) return true;
  else if (Z > 0 && A > 0)                              return true;
#ifdef debug
  G4cout<<"***>>G4QLowEnergy::IsApplicable: projPDG="<<particle.GetPDGEncoding()<<", A="
        <<A<<", Z="<<Z<<G4endl;
#endif
  return false;
}

G4VParticleChange* G4QLowEnergy::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  static const G4double mProt= G4QPDGCode(2212).GetMass()/MeV; // CHIPS proton Mass in MeV
  static const G4double mPro2= mProt*mProt;                    // CHIPS sq proton Mass
  static const G4double mNeut= G4QPDGCode(2112).GetMass()/MeV; // CHIPS neutron Mass in MeV
  static const G4double mNeu2= mNeut*mNeut;                    // CHIPS sq neutron Mass
  static const G4double mLamb= G4QPDGCode(3122).GetMass()/MeV; // CHIPS Lambda Mass in MeV
  static const G4double mDeut= G4QPDGCode(2112).GetNuclMass(1,1,0)/MeV;
  static const G4double mTrit= G4QPDGCode(2112).GetNuclMass(1,2,0)/MeV;
  static const G4double mHel3= G4QPDGCode(2112).GetNuclMass(2,1,0)/MeV;
  static const G4double mAlph= G4QPDGCode(2112).GetNuclMass(2,2,0)/MeV;
  static const G4double mFm= 250*MeV;
  static const G4double third= 1./3.;
  static const G4ThreeVector zeroMom(0.,0.,0.);
  static G4ParticleDefinition* aGamma    = G4Gamma::Gamma();
  static G4ParticleDefinition* aPiZero   = G4PionZero::PionZero();
  static G4ParticleDefinition* aPiPlus   = G4PionPlus::PionPlus();
  static G4ParticleDefinition* aPiMinus  = G4PionMinus::PionMinus();
  static G4ParticleDefinition* aProton   = G4Proton::Proton();
  static G4ParticleDefinition* aNeutron  = G4Neutron::Neutron();
  static G4ParticleDefinition* aLambda   = G4Lambda::Lambda();
  static G4ParticleDefinition* aDeuteron = G4Deuteron::Deuteron();
  static G4ParticleDefinition* aTriton   = G4Triton::Triton();
  static G4ParticleDefinition* aHe3      = G4He3::He3();
  static G4ParticleDefinition* anAlpha   = G4Alpha::Alpha();
  static const G4int nCh=26;                         // #of combinations
  static G4QNucleus Nuc;                             // A fake nucleus to call Evaporation
  //
  //-------------------------------------------------------------------------------------
  static G4bool CWinit = true;                       // CHIPS Warld needs to be initted
  if(CWinit)
  {
    CWinit=false;
    G4QCHIPSWorld::Get()->GetParticles(nPartCWorld); // Create CHIPS World (234 part.max)
  }
  //-------------------------------------------------------------------------------------
  const G4DynamicParticle* projHadron = track.GetDynamicParticle();
  const G4ParticleDefinition* particle=projHadron->GetDefinition();
#ifdef pdebug
  G4cout<<"G4QLowEnergy::PostStepDoIt: *** Called *** In4M="
        <<projHadron->Get4Momentum()<<" of PDG="<<particle->GetPDGEncoding()<<", Type="
        <<particle->GetParticleType()<<",SubType="<<particle->GetParticleSubType()<<G4endl;
#endif
  //G4ForceCondition cond=NotForced;
  //GetMeanFreePath(track, -27., &cond);              // @@ ?? jus to update parameters?
#ifdef debug
  G4cout<<"G4QLowEnergy::PostStepDoIt: After GetMeanFreePath is called"<<G4endl;
#endif
  std::vector<G4Track*> result;
  G4LorentzVector proj4M=(projHadron->Get4Momentum())/MeV; // Convert to MeV!
  G4double momentum = projHadron->GetTotalMomentum()/MeV; // 3-momentum of the Proj in MeV
  G4double Momentum = proj4M.rho();                   // @@ Just for the test purposes
  if(std::fabs(Momentum-momentum)>.000001)
    G4cerr<<"-Warning-G4QLowEnergy::PostStepDoIt:P_IU="<<Momentum<<"#"<<momentum<<G4endl;
#ifdef debug
  G4cout<<"G4QLowEnergy::PostStepDoIt: pP(IU)="<<Momentum<<"="<<momentum<<",proj4M,m="
        <<proj4M<<proj4M.m()<<G4endl;
#endif
  if (!IsApplicable(*particle))  // Check applicability
  {
    G4cerr<<"G4QLowEnergy::PostStepDoIt: Only NA is implemented."<<G4endl;
    return 0;
  }
  const G4Material* material = track.GetMaterial();      // Get the current material
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QLowEnergy::PostStepDoIt: "<<nE<<" elements in the material."<<G4endl;
#endif
  G4int projPDG=0;                           // PDG Code prototype for the captured hadron
  // Not all these particles are implemented yet (see Is Applicable)
  G4int Z=static_cast<G4int>(particle->GetPDGCharge());
  G4int A=particle->GetBaryonNumber();
  if      (particle ==      G4Proton::Proton()     ) projPDG= 2212;
  else if (particle ==     G4Neutron::Neutron()    ) projPDG= 2112;
  else if (particle ==    G4Deuteron::Deuteron()   ) projPDG= 1000010020;
  else if (particle ==       G4Alpha::Alpha()      ) projPDG= 1000020040;
  else if (particle ==      G4Triton::Triton()     ) projPDG= 1000010030;
  else if (particle ==         G4He3::He3()        ) projPDG= 1000020030;
  else if (particle ==  G4GenericIon::GenericIon() || (Z > 0 && A > 0))
  {
    projPDG=particle->GetPDGEncoding();
#ifdef debug
    G4int prPDG=1000000000+10000*Z+10*A;
    G4cout<<"G4QLowEnergy::PostStepDoIt: PDG="<<prPDG<<"="<<projPDG<<G4endl;
#endif
  }
  else G4cout<<"-Warning-G4QLowEnergy::PostStepDoIt:Unknown projectile Ion"<<G4endl;
#ifdef pdebug
  G4int prPDG=particle->GetPDGEncoding();
  G4cout<<"G4QLowEnergy::PostStepDoIt: projPDG="<<projPDG<<", stPDG="<<prPDG<<G4endl;
#endif
  if(!projPDG)
  {
    G4cerr<<"-Warning-G4QLowEnergy::PostStepDoIt:UndefProjHadron(PDG=0) ->ret 0"<<G4endl;
    return 0;
  }
  // Element treatment
  G4int EPIM=ElProbInMat.size();
#ifdef debug
  G4cout<<"G4QLowEn::PostStDoIt: m="<<EPIM<<", n="<<nE<<",T="<<ElProbInMat[EPIM-1]<<G4endl;
#endif
  G4int i=0;
  if(EPIM>1)
  {
    G4double rnd = ElProbInMat[EPIM-1]*G4UniformRand();
    for(i=0; i<nE; ++i)
    {
#ifdef debug
      G4cout<<"G4QLowEn::PostStepDoIt: EPM["<<i<<"]="<<ElProbInMat[i]<<", r="<<rnd<<G4endl;
#endif
      if (rnd<ElProbInMat[i]) break;
    }
    if(i>=nE) i=nE-1;                        // Top limit for the Element
  }
  G4Element* pElement=(*theElementVector)[i];
  G4int tZ=static_cast<G4int>(pElement->GetZ());
#ifdef debug
    G4cout<<"G4QLowEnergy::PostStepDoIt: i="<<i<<", Z(element)="<<tZ<<G4endl;
#endif
  if(tZ<=0)
  {
    G4cerr<<"-Warning-G4QLowEnergy::PostStepDoIt: Element with Z="<<tZ<<G4endl;
    if(tZ<0) return 0;
  }
  std::vector<G4double>* SPI = IsoProbInEl[i];// Vector of summedProbabilities for isotopes
  std::vector<G4int>* IsN = ElIsoN[i];     // Vector of "#of neutrons" in the isotope El[i]
  G4int nofIsot=SPI->size();               // #of isotopes in the element i
#ifdef debug
  G4cout<<"G4QLowEnergy::PostStepDoIt: nI="<<nofIsot<<", T="<<(*SPI)[nofIsot-1]<<G4endl;
#endif
  G4int j=0;
  if(nofIsot>1)
  {
    G4double rndI=(*SPI)[nofIsot-1]*G4UniformRand(); // Randomize the isotop of the Element
    for(j=0; j<nofIsot; ++j)
    {
#ifdef debug
      G4cout<<"G4QLowEnergy::PostStepDoIt: SP["<<j<<"]="<<(*SPI)[j]<<",r="<<rndI<<G4endl;
#endif
      if(rndI < (*SPI)[j]) break;
    }
    if(j>=nofIsot) j=nofIsot-1;            // Top limit for the isotope
  }
  G4int tN =(*IsN)[j]; ;                    // Randomized number of neutrons
#ifdef debug
  G4cout<<"G4QLowEnergy::PostStepDoIt: j="<<i<<", N(isotope)="<<tN<<", MeV="<<MeV<<G4endl;
#endif
  if(tN<0)
  {
    G4cerr<<"-Warning-G4QLowEnergy::PostStepDoIt: Isotope Z="<<tZ<<" has 0>N="<<tN<<G4endl;
    return 0;
  }
  nOfNeutrons=tN;                           // Remember it for the energy-momentum check
#ifdef debug
  G4cout<<"G4QLowEnergy::PostStepDoIt: N="<<tN<<" for element with Z="<<tZ<<G4endl;
#endif
  if(tN<0)
  {
    G4cerr<<"***Warning***G4QLowEnergy::PostStepDoIt: Element with N="<<tN<< G4endl;
    return 0;
  }
  aParticleChange.Initialize(track);
#ifdef debug
  G4cout<<"G4QLowEnergy::PostStepDoIt: track is initialized"<<G4endl;
#endif
  G4double weight        = track.GetWeight();
  G4double localtime     = track.GetGlobalTime();
  G4ThreeVector position = track.GetPosition();
#ifdef debug
  G4cout<<"G4QLowEnergy::PostStepDoIt: before Touchable extraction"<<G4endl;
#endif
  G4TouchableHandle trTouchable = track.GetTouchableHandle();
#ifdef debug
  G4cout<<"G4QLowEnergy::PostStepDoIt: Touchable is extracted"<<G4endl;
#endif
  G4QPDGCode targQPDG(90000000+tZ*1000+tN);  // @@ use G4Ion and get rid of CHIPS World
  G4double tM=targQPDG.GetMass();            // CHIPS target nucleus mass in MeV
  G4int pL=particle->GetQuarkContent(3)-particle->GetAntiQuarkContent(3); // Strangeness
  G4int pZ=static_cast<G4int>(particle->GetPDGCharge());  // Charge of the projectile
  G4int pN=particle->GetBaryonNumber()-pZ-pL;// #of neutrons in projectile
  G4double pM=targQPDG.GetNuclMass(pZ,pN,0); // CHIPS projectile nucleus mass in MeV
  G4double cosp=-14*Momentum*(tM-pM)/tM/pM;  // Asymmetry power for angular distribution
#ifdef debug
  G4cout<<"G4QLowEnergy::PoStDoIt: Proj("<<pZ<<","<<pN<<","<<pL<<")p="<<pM<<",Targ=("<<tZ
        <<","<<tN<<"), cosp="<<cosp<<G4endl;
#endif
  G4double kinEnergy= projHadron->GetKineticEnergy()*MeV; // Kin energy in MeV (Is *MeV n?)
  G4ParticleMomentum dir = projHadron->GetMomentumDirection();// It is a unit three-vector
  G4LorentzVector targ4M=G4LorentzVector(0.,0.,0.,tM); // Target's 4-mom
  G4LorentzVector tot4M =proj4M+targ4M;      // Total 4-mom of the reaction
#ifdef pdebug
  G4cout<<"G4QLowEnergy::PostStepDoIt: tM="<<tM<<",p4M="<<proj4M<<",t4M="<<tot4M<<G4endl;
#endif
  EnMomConservation=tot4M;                 // Total 4-mom of reaction for E/M conservation
  // @@ Probably this is not necessary any more
#ifdef debug
  G4cout<<"G4QLE::PSDI:false,P="<<Momentum<<",Z="<<pZ<<",N="<<pN<<",PDG="<<projPDG<<G4endl;
#endif
  G4double xSec=CalculateXS(Momentum, tZ, tN, projPDG); // Recalculate CrossSection
#ifdef pdebug
  G4cout<<"G4QLowEn::PSDI:PDG="<<projPDG<<",P="<<Momentum<<",tZ="<<tZ<<",N="<<tN<<",XS="
        <<xSec/millibarn<<G4endl;
#endif
#ifdef nandebug
  if(xSec>0. || xSec<0. || xSec==0);
  else  G4cout<<"-Warning-G4QLowEnergy::PostSDI: *NAN* xSec="<<xSec/millibarn<<G4endl;
#endif
  // @@ check a possibility to separate p, n, or alpha (!)
  if(xSec <= 0.) // The cross-section iz 0 -> Do Nothing
  {
#ifdef debug
    G4cerr<<"-Warning-G4QLowEnergy::PSDoIt:*Zero cross-section* PDG="<<projPDG
          <<",Z="<<tZ<<",tN="<<tN<<",P="<<Momentum<<G4endl;
#endif
    //Do Nothing Action insead of the reaction
    aParticleChange.ProposeEnergy(kinEnergy);
    aParticleChange.ProposeLocalEnergyDeposit(0.);
    aParticleChange.ProposeMomentumDirection(dir) ;
    return G4VDiscreteProcess::PostStepDoIt(track,step);
  }
  // Kill interacting hadron
  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  G4int tB=tZ+tN;
  G4int pB=pZ+pN;
#ifdef pdebug
  G4cout<<"G4QLowEn::PSDI: Projectile track is killed"<<", tA="<<tB<<", pA="<<pB<<G4endl;
#endif
  // algorithm implementation STARTS HERE --- All calculations are in IU --------
  G4double tA=tB;
  G4double pA=pB;
  G4double tR=1.1;                    // target nucleus R in fm
  if(tB > 1) tR*=std::pow(tA,third);  // in fm
  G4double pR=1.1;                    // projectile nucleus R in fm
  if(pB > 1) pR*=std::pow(pA,third);  // in fm
  G4double R=tR+pR;                   // total radius
  G4double R2=R*R;
  G4int tD=0;
  G4int pD=0;
  G4int nAt=0;
  G4int nAtM=27;
  G4int nSec = 0;
  G4double tcM=0.;
  G4double tnM=1.;
#ifdef edebug
  G4int totChg=0;
  G4int totBaN=0;
  G4LorentzVector tch4M =tot4M;       // Total 4-mom of the reaction
#endif
  while((!tD && !pD && ++nAt<nAtM) || tcM<tnM)
  {
#ifdef edebug
    totChg=tZ+pZ;
    totBaN=tB+pB;
    tch4M =tot4M;
    G4cout<<">G4QLEn::PSDI:#"<<nAt<<",tC="<<totChg<<",tA="<<totBaN<<",t4M="<<tch4M<<G4endl;
#endif
    G4LorentzVector tt4M=tot4M;
    G4int resN=result.size();
    if(resN)
    {
      for(G4int i=0; i<resN; ++i) delete result[i];
      result.clear();
    }
    G4double D=std::sqrt(R2*G4UniformRand());
#ifdef pdebug
    G4cout<<"G4QLowEn::PSDI: D="<<D<<", tR="<<tR<<", pR="<<pR<<G4endl;
#endif
    if(D > std::fabs(tR-pR))          // leading parts can be separated
    {
      nSec = 0;
      G4double DtR=D-tR;
      G4double DpR=D-pR;
      G4double DtR2=DtR*DtR;
      G4double DpR2=DpR*DpR;
      //G4double DtR3=DtR2*DtR;
      G4double DpR3=DpR2*DpR;
      //G4double DtR4=DtR3*DtR;
      G4double DpR4=DpR3*DpR;
      G4double tR2=tR*tR;
      G4double pR2=pR*pR;
      G4double tR3=tR2*tR;
      //G4double pR3=pR2*pR;
      G4double tR4=tR3*tR;
      //G4double pR4=pR3*pR;
      G4double HD=16.*D;
      G4double tF=tA*(6*tR2*(pR2-DtR2)-4*D*(tR3-DpR3)+3*(tR4-DpR4))/HD/tR3;
      G4double pF=tF;
      tD=static_cast<G4int>(tF);
      pD=static_cast<G4int>(pF);
      //if(G4UniformRand() < tF-tD) ++tD; // Simple solution
      //if(G4UniformRand() < pF-pD) ++pD;
      // Enhance alphas solution
      if(std::fabs(tF-4.)<1.) tD=4;       // @@ 1. is the enhancement parameter
      else if(G4UniformRand() < tF-tD) ++tD;
      if(std::fabs(pF-4.)<1.) pD=4;
      else if(G4UniformRand() < pF-pD) ++pD;
      if(tD > tB) tD=tB;
      if(pD > pB) pD=tB;
      // @@ Quasi Free is not debugged @@ The following close it
      if(tD < 1) tD=1;
      if(pD < 1) pD=1;
#ifdef pdebug
      G4cout<<"G4QLE::PSDI:#"<<nAt<<",pF="<<pF<<",tF="<<tF<<",pD="<<pD<<",tD="<<tD<<G4endl;
#endif
      G4int pC=0;                     // charge of the projectile fraction
      G4int tC=0;                     // charge of the target fraction
      if((tD || pD) && tD<tB && pD<pB)// Periferal interaction
      {
        if(!tD || !pD)                // Quasi-Elastic: B+A->(B-1)+N+A || ->B+N+(A-1)
        {
          G4VQCrossSection* PCSmanager=G4QProtonElasticCrossSection::GetPointer();
          G4VQCrossSection* NCSmanager=G4QNeutronElasticCrossSection::GetPointer();
          G4int pPDG=2112;            // Proto of the nucleon PDG (proton)
          G4double prM =mNeut;        // Proto of the nucleon mass
          G4double prM2=mNeu2;        // Proto of the nucleon sq mass
          if     (!tD)                // Quasi-elastic scattering of the proj QF nucleon
          {
            if(pD > 1) pD=1;
            if(!pN || (pZ && pA*G4UniformRand() < pZ) ) // proj QF proton
            {
              pPDG=2212;
              prM=mProt;
              prM2=mPro2;
              pC=1;
            }
            G4double Mom=0.;
            G4LorentzVector com4M = targ4M; // Proto of 4mom of compound
            G4double tgM=targ4M.e();
            G4double rNM=0.;
            G4LorentzVector rNuc4M(0.,0.,0.,0.);
            if(pD==pB)
            {
              Mom=proj4M.rho();
              com4M += proj4M;        // Total 4-mom for scattering
              rNM=prM;
            }
            else                      // It is necessary to split the nucleon (with fermiM)
            {
              G4ThreeVector fm=mFm*std::pow(G4UniformRand(),third)*G4RandomDirection();
              rNM=G4QPDGCode(2112).GetNuclMass(pZ-pC, pN, 0);
              G4double rNE=std::sqrt(fm*fm+rNM*rNM);
              rNuc4M=G4LorentzVector(fm,rNE);
              G4ThreeVector boostV=proj4M.vect()/proj4M.e();
              rNuc4M.boost(boostV);
              G4LorentzVector qfN4M=proj4M - rNuc4M;// 4Mom of the quasi-free nucleon in LS
              com4M += qfN4M;         // Calculate Total 4Mom for NA scattering
              G4double pNE = qfN4M.e(); // Energy of the QF nucleon in LS
              if(rNE >= prM) Mom = std::sqrt(pNE*pNE-prM2); // Mom(s) fake value
              else break;             // Break the while loop
            }
            G4double xSec=0.;
            if(pPDG==2212) xSec=PCSmanager->GetCrossSection(false, Mom, tZ, tN, pPDG);
            else           xSec=NCSmanager->GetCrossSection(false, Mom, tZ, tN, pPDG);
            if( xSec <= 0. ) break;   // Break the while loop
            G4double mint=0.;                 // Prototype of functional randomized -t
            if(pPDG==2212) mint=PCSmanager->GetExchangeT(tZ,tN,pPDG); // randomized -t
            else           mint=NCSmanager->GetExchangeT(tZ,tN,pPDG); // randomized -t
            G4double maxt=0.;                           // Prototype of maximum -t
            if(pPDG==2212) maxt=PCSmanager->GetHMaxT(); // maximum -t
            else           maxt=NCSmanager->GetHMaxT(); // maximum -t
            G4double cost=1.-mint/maxt;       // cos(theta) in CMS
            if(cost>1. || cost<-1.) break; // Break the while loop
            G4LorentzVector reco4M=G4LorentzVector(0.,0.,0.,tgM); // 4mom of recoil target
            G4LorentzVector scat4M=G4LorentzVector(0.,0.,0.,rNM); // 4mom of scattered N
            G4LorentzVector dir4M=tt4M-G4LorentzVector(0.,0.,0.,(com4M.e()-rNM-prM)*.01);
            if(!G4QHadron(com4M).RelDecayIn2(scat4M, reco4M, dir4M, cost, cost))
            {
              G4cout<<"G4QLE::Pt="<<com4M.m()<<",p="<<prM<<",r="<<rNM<<",c="<<cost<<G4endl;
              break;                  // Break the while loop
            }
            G4Track* projSpect = 0;
            G4Track* aNucTrack = 0;
            if(pB > pD)               // Fill the proj residual nucleus
            {
              G4int rZ=pZ-pC;
              G4int rA=pB-1;
              G4ParticleDefinition* theDefinition;// Prototype of residualNucleusDefinition
              if(rA==1)
              {
                if(!rZ) theDefinition = aNeutron;
                else    theDefinition = aProton;
              }
              else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,0,rZ);
              G4DynamicParticle* resN = new G4DynamicParticle(theDefinition, rNuc4M);
              projSpect = new G4Track(resN, localtime, position );
              projSpect->SetWeight(weight);                    //    weighted
              projSpect->SetTouchableHandle(trTouchable);
#ifdef pdebug
             G4cout<<"G4QLowEn::PSDI:-->ProjQFSA="<<rA<<",rZ="<<rZ<<",4M="<<rNuc4M<<G4endl;
#endif
              ++nSec;
            }
            G4ParticleDefinition* theDefinition = G4Neutron::Neutron();
            if(pPDG==2212) theDefinition = G4Proton::Proton();
            G4DynamicParticle* scatN = new G4DynamicParticle(theDefinition, scat4M);
            aNucTrack = new G4Track(scatN, localtime, position );
            aNucTrack->SetWeight(weight);                    //    weighted
            aNucTrack->SetTouchableHandle(trTouchable);
#ifdef pdebug
             G4cout<<"G4QLowEn::PSDI:-->TgQFNPDG="<<pPDG<<", 4M="<<scat4M<<G4endl;
#endif
            ++nSec;
            G4Track* aFraTrack=0;
            if(tA==1)
            {
              if(!tZ) theDefinition = aNeutron;
              else    theDefinition = aProton;
            }
            else if(tA==8 && tC==4)             // Be8 decay
            {
              theDefinition = anAlpha;
              G4LorentzVector f4M=G4LorentzVector(0.,0.,0.,mAlph); // 4mom of 1st alpha
              G4LorentzVector s4M=G4LorentzVector(0.,0.,0.,mAlph); // 4mom of 2nd alpha
              if(!G4QHadron(reco4M).DecayIn2(f4M, s4M))
              {
                G4cout<<"*G4QLE::TS->A+A,t="<<reco4M.m()<<" >? 2*MAlpha="<<2*mAlph<<G4endl;
              }
              G4DynamicParticle* pAl = new G4DynamicParticle(theDefinition, f4M);
              aFraTrack = new G4Track(pAl, localtime, position );
              aFraTrack->SetWeight(weight);                    //    weighted
              aFraTrack->SetTouchableHandle(trTouchable);
#ifdef pdebug
              G4cout<<"G4QLowEn::PSDI:-->TgRQFA4M="<<f4M<<G4endl;
#endif
              ++nSec;
              reco4M=s4M;
            }
            else if(tA==5 && (tC==2 || tC==3))   // He5/Li5 decay
            {
              theDefinition = aProton;
              G4double mNuc = mProt;
              if(tC==2)
              {
                theDefinition = aNeutron;
                mNuc = mNeut;
              }
              G4LorentzVector f4M=G4LorentzVector(0.,0.,0.,mNuc);  // 4mom of the nucleon
              G4LorentzVector s4M=G4LorentzVector(0.,0.,0.,mAlph); // 4mom of the alpha
              if(!G4QHadron(reco4M).DecayIn2(f4M, s4M))
              {
                G4cout<<"*G4QLE::TS->N+A,t="<<reco4M.m()<<" >? MA+MN="<<mAlph+mNuc<<G4endl;
              }
              G4DynamicParticle* pAl = new G4DynamicParticle(theDefinition, f4M);
              aFraTrack = new G4Track(pAl, localtime, position );
              aFraTrack->SetWeight(weight);                    //    weighted
              aFraTrack->SetTouchableHandle(trTouchable);
#ifdef pdebug
              G4cout<<"G4QLowEn::PSDI:-->TgQFRN4M="<<f4M<<G4endl;
#endif
              ++nSec;
              theDefinition = anAlpha;
              reco4M=s4M;
            }
            else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(tZ,tB,0,tZ);
            ++nSec;
#ifdef pdebug
            G4cout<<"G4QLowEn::PSDI:-->PQF_nSec="<<nSec<<G4endl;
#endif
            aParticleChange.SetNumberOfSecondaries(nSec); 
            if(projSpect) aParticleChange.AddSecondary(projSpect);
            if(aNucTrack) aParticleChange.AddSecondary(aNucTrack);
            if(aFraTrack) aParticleChange.AddSecondary(aFraTrack);
            G4DynamicParticle* resA = new G4DynamicParticle(theDefinition, reco4M);
            G4Track* aResTrack = new G4Track(resA, localtime, position );
            aResTrack->SetWeight(weight);                    //    weighted
            aResTrack->SetTouchableHandle(trTouchable);
#ifdef pdebug
            G4cout<<"G4QLowEn::PSDI:-->TgR4M="<<reco4M<<", checkNSec="<<nSec<<G4endl;
#endif
            aParticleChange.AddSecondary(aResTrack);
          }
          else                         // !pD : QF target Nucleon on the whole Projectile
          {
            if(tD > 1) tD=1;
            if(!tN || (tZ && tA*G4UniformRand() < tZ) ) // target QF proton
            {
              pPDG=2212;
              prM=mProt;
              prM2=mPro2;
              tC=1;
            }
            G4double Mom=0.;
            G4LorentzVector com4M=proj4M; // Proto of 4mom of compound
            G4double prM=proj4M.m();
            G4double rNM=0.;
            G4LorentzVector rNuc4M(0.,0.,0.,0.);
            if(tD==tB)
            {
              Mom=prM*proj4M.rho()/proj4M.m();
              com4M += targ4M;        // Total 4-mom for scattering
              rNM=prM;
            }
            else                      // It is necessary to split the nucleon (with fermiM)
            {
              G4ThreeVector fm=250.*std::pow(G4UniformRand(),third)*G4RandomDirection();
              rNM=G4QPDGCode(2112).GetNuclMass(tZ-tC, tN, 0)/MeV;
              G4double rNE=std::sqrt(fm*fm+rNM*rNM);
              rNuc4M=G4LorentzVector(fm,rNE);
              G4LorentzVector qfN4M=targ4M - rNuc4M;// 4Mom of the quasi-free nucleon in LS
              com4M += qfN4M;                 // Calculate Total 4Mom for NA scattering
              G4ThreeVector boostV=proj4M.vect()/proj4M.e();
              qfN4M.boost(-boostV);
              G4double tNE = qfN4M.e();       // Energy of the QF nucleon in LS
              if(rNE >= prM) Mom = std::sqrt(tNE*tNE-prM2); // Mom(s) fake value
              else break;                     // Break the while loop
            }
            G4double xSec=0.;
            if(pPDG==2212) xSec=PCSmanager->GetCrossSection(false, Mom, tZ, tN, pPDG);
            else           xSec=NCSmanager->GetCrossSection(false, Mom, tZ, tN, pPDG);
            if( xSec <= 0. ) break;           // Break the while loop
            G4double mint=0.;                 // Prototype of functional randomized -t
            if(pPDG==2212) mint=PCSmanager->GetExchangeT(tZ,tN,pPDG); // randomized -t
            else           mint=NCSmanager->GetExchangeT(tZ,tN,pPDG); // randomized -t
            G4double maxt=0.;                           // Prototype of maximum -t
            if(pPDG==2212) maxt=PCSmanager->GetHMaxT(); // maximum -t
            else           maxt=NCSmanager->GetHMaxT(); // maximum -t
            G4double cost=1.-mint/maxt;                 // cos(theta) in CMS
            if(cost>1. || cost<-1.) break;    // Break the while loop
            G4LorentzVector reco4M=G4LorentzVector(0.,0.,0.,prM); // 4mom of recoil target
            G4LorentzVector scat4M=G4LorentzVector(0.,0.,0.,rNM); // 4mom of scattered N
            G4LorentzVector dir4M=tt4M-G4LorentzVector(0.,0.,0.,(com4M.e()-rNM-prM)*.01);
            if(!G4QHadron(com4M).RelDecayIn2(scat4M, reco4M, dir4M, cost, cost))
            {
              G4cout<<"G4QLE::Tt="<<com4M.m()<<",p="<<prM<<",r="<<rNM<<",c="<<cost<<G4endl;
              break;                          // Break the while loop
            }
            G4Track* targSpect = 0;
            G4Track* aNucTrack = 0;
            if(tB > tD)                       // Fill the residual nucleus
            {
              G4int rZ=tZ-tC;
              G4int rA=tB-1;
              G4ParticleDefinition* theDefinition;// Prototype of residualNucleusDefinition
              if(rA==1)
              {
                if(!rZ) theDefinition = aNeutron;
                else    theDefinition = aProton;
              }
              else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,0,rZ);
              G4DynamicParticle* resN = new G4DynamicParticle(theDefinition, rNuc4M);
              targSpect = new G4Track(resN, localtime, position );
              targSpect->SetWeight(weight);                    //    weighted
              targSpect->SetTouchableHandle(trTouchable);
#ifdef pdebug
             G4cout<<"G4QLowEn::PSDI:-->TargQFSA="<<rA<<",rZ="<<rZ<<",4M="<<rNuc4M<<G4endl;
#endif
              ++nSec;
            }
            G4ParticleDefinition* theDefinition = G4Neutron::Neutron();
            if(pPDG==2212) theDefinition = G4Proton::Proton();
            G4DynamicParticle* scatN = new G4DynamicParticle(theDefinition, scat4M);
            aNucTrack = new G4Track(scatN, localtime, position );
            aNucTrack->SetWeight(weight);                    //    weighted
            aNucTrack->SetTouchableHandle(trTouchable);
#ifdef pdebug
             G4cout<<"G4QLowEn::PSDI:-->PrQFNPDG="<<pPDG<<", 4M="<<scat4M<<G4endl;
#endif
            ++nSec;
            G4Track* aFraTrack=0;
            if(pA==1)
            {
              if(!pZ) theDefinition = aNeutron;
              else    theDefinition = aProton;
            }
            else if(pA==8 && pC==4)             // Be8 decay
            {
              theDefinition = anAlpha;
              G4LorentzVector f4M=G4LorentzVector(0.,0.,0.,mAlph); // 4mom of 1st alpha
              G4LorentzVector s4M=G4LorentzVector(0.,0.,0.,mAlph); // 4mom of 2nd alpha
              if(!G4QHadron(reco4M).DecayIn2(f4M, s4M))
              {
                G4cout<<"*G4QLE::PS->A+A,t="<<reco4M.m()<<" >? 2*MAlpha="<<2*mAlph<<G4endl;
              }
              G4DynamicParticle* pAl = new G4DynamicParticle(theDefinition, f4M);
              aFraTrack = new G4Track(pAl, localtime, position );
              aFraTrack->SetWeight(weight);                    //    weighted
              aFraTrack->SetTouchableHandle(trTouchable);
#ifdef pdebug
              G4cout<<"G4QLowEn::PSDI:-->PrRQFA4M="<<f4M<<G4endl;
#endif
              ++nSec;
              reco4M=s4M;
            }
            else if(pA==5 && (pC==2 || pC==3))   // He5/Li5 decay
            {
              theDefinition = aProton;
              G4double mNuc = mProt;
              if(pC==2)
              {
                theDefinition = aNeutron;
                mNuc = mNeut;
              }
              G4LorentzVector f4M=G4LorentzVector(0.,0.,0.,mNuc);  // 4mom of the nucleon
              G4LorentzVector s4M=G4LorentzVector(0.,0.,0.,mAlph); // 4mom of the alpha
              if(!G4QHadron(reco4M).DecayIn2(f4M, s4M))
              {
                G4cout<<"*G4QLE::PS->N+A,t="<<reco4M.m()<<" >? MA+MN="<<mAlph+mNuc<<G4endl;
              }
              G4DynamicParticle* pAl = new G4DynamicParticle(theDefinition, f4M);
              aFraTrack = new G4Track(pAl, localtime, position );
              aFraTrack->SetWeight(weight);                    //    weighted
              aFraTrack->SetTouchableHandle(trTouchable);
#ifdef pdebug
              G4cout<<"G4QLowEn::PSDI:-->PrQFRN4M="<<f4M<<G4endl;
#endif
              ++nSec;
              theDefinition = anAlpha;
              reco4M=s4M;
            }
            else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(pZ,pB,0,pZ);
            ++nSec;
#ifdef pdebug
            G4cout<<"G4QLowEn::PSDI:-->TQF_nSec="<<nSec<<G4endl;
#endif
            aParticleChange.SetNumberOfSecondaries(nSec); 
            if(targSpect) aParticleChange.AddSecondary(targSpect);
            if(aNucTrack) aParticleChange.AddSecondary(aNucTrack);
            if(aFraTrack) aParticleChange.AddSecondary(aFraTrack);
            G4DynamicParticle* resA = new G4DynamicParticle(theDefinition, reco4M);
            G4Track* aResTrack = new G4Track(resA, localtime, position );
            aResTrack->SetWeight(weight);                    //    weighted
            aResTrack->SetTouchableHandle(trTouchable);
#ifdef pdebug
            G4cout<<"G4QLowEn::PSDI:-->TgR4M="<<reco4M<<", checkNSec="<<nSec<<G4endl;
#endif
            aParticleChange.AddSecondary( aResTrack );
          }
#ifdef debug
          G4cout<<"G4QLowEnergy::PostStepDoIt:***PostStepDoIt is done:Quasi-El***"<<G4endl;
#endif
          return G4VDiscreteProcess::PostStepDoIt(track, step);
        }
        else                          // The cental region compound can be created
        {
          // First calculate the isotopic state of the parts of the compound
          if(!pZ) pC=0;
          else if(!pN) pC=pD;
          else
          {
#ifdef pdebug
            G4cout<<"G4QLowEn::PSDI: pD="<<pD<<", pZ="<<pZ<<", pA="<<pA<<G4endl;
#endif
            G4double C=pD*pZ/pA;
            pC=static_cast<G4int>(C); 
            if(G4UniformRand() < C-pC) ++pC;
          }
          if(!tZ) tC=0;
          else if(!tN) tC=tD;
          else
          {
#ifdef pdebug
            G4cout<<"G4QLowEn::PSDI: tD="<<tD<<", tZ="<<tZ<<", tA="<<tA<<G4endl;
#endif
            G4double C=tD*tZ/tA;
            tC=static_cast<G4int>(C); 
            if(G4UniformRand() < C-tC) ++tC;
          }
          // calculate the transferred momentum
          G4ThreeVector pFM(0.,0.,0.);
          if(pD<pB)                    // The projectile nucleus must be splitted
          {
            G4int nc=pD;
            if(pD+pD>pB) nc=pB-pD;
            pFM = mFm*std::pow(G4UniformRand(),third)*G4RandomDirection();
            for(G4int i=1; i < nc; ++i)
                             pFM+= mFm*std::pow(G4UniformRand(),third)*G4RandomDirection();
          }
          G4ThreeVector tFM(0.,0.,0.);
          if(tD<tB)                    // The projectile nucleus must be splitted
          {
            G4int nc=pD;
            if(tD+tD>tB) nc=tB-tD;
            tFM = mFm*std::pow(G4UniformRand(),third)*G4RandomDirection();
            for(G4int i=1; i < nc; ++i)
                             tFM+= mFm*std::pow(G4UniformRand(),third)*G4RandomDirection();
          }
#ifdef pdebug
          G4cout<<"G4QLE::PSDI:pC="<<pC<<", tC="<<tC<<", pFM="<<pFM<<", tFM="<<tFM<<G4endl;
#endif
          // Split the projectile spectator
          G4int rpZ=pZ-pC;            // number of protons in the projectile spectator
          G4int pF=pD-pC;             // number of neutrons in the projectile part of comp
          G4int rpN=pN-pF;            // number of neutrons in the projectile spectator
          G4double rpNM=G4QPDGCode(2112).GetNuclMass(rpZ, rpN, 0); // Mass of the spectator
          G4ThreeVector boostV=proj4M.vect()/proj4M.e(); // Antilab Boost Vector
          G4double rpE=std::sqrt(rpNM*rpNM+pFM.mag2());
          G4LorentzVector rp4M(pFM,rpE);
#ifdef pdebug
	    G4cout<<"G4QLE::PSDI: boostV="<<boostV<<",rp4M="<<rp4M<<",pr4M="<<proj4M<<G4endl;
#endif
          rp4M.boost(boostV);
#ifdef pdebug
	    G4cout<<"G4QLE::PSDI: After boost, rp4M="<<rp4M<<G4endl;
#endif
          G4ParticleDefinition* theDefinition; // Prototype of projSpectatorNuclDefinition
          G4int rpA=rpZ+rpN;
          G4Track* aFraPTrack = 0;
          theDefinition = 0;
          if(rpA==1)
          {
            if(!rpZ) theDefinition = G4Neutron::Neutron();
            else     theDefinition = G4Proton::Proton();
#ifdef pdebug
            G4cout<<"G4QLE::PSDI: rpA=1, rpZ"<<rpZ<<G4endl;
#endif
          }
          else if(rpA==2 && rpZ==0)            // nn decay
          {
            theDefinition = aNeutron;
            G4LorentzVector f4M=G4LorentzVector(0.,0.,0.,mNeut); // 4mom of 1st neutron
            G4LorentzVector s4M=G4LorentzVector(0.,0.,0.,mNeut); // 4mom of 2nd neutron
#ifdef pdebug
            G4cout<<"G4QLE::CPS->n+n,nn="<<rp4M.m()<<" >? 2*MNeutron="<<2*mNeutron<<G4endl;
#endif
            if(!G4QHadron(rp4M).DecayIn2(f4M, s4M))
            {
              G4cout<<"*W*G4QLE::CPS->n+n,t="<<rp4M.m()<<" >? 2*Neutron="<<2*mAlph<<G4endl;
            }
            G4DynamicParticle* pNeu = new G4DynamicParticle(theDefinition, f4M);
            aFraPTrack = new G4Track(pNeu, localtime, position );
            aFraPTrack->SetWeight(weight);                    //    weighted
            aFraPTrack->SetTouchableHandle(trTouchable);
            tt4M-=f4M;
#ifdef edebug
            totBaN-=2;
            tch4M -=f4M;
            G4cout<<">>G4QLEn::PSDI:n,tZ="<<totChg<<",tB="<<totBaN<<",t4M="<<tch4M<<G4endl;
#endif
#ifdef pdebug
            G4cout<<"G4QLowEn::PSDI:-->ProjSpectA4M="<<f4M<<G4endl;
#endif
            ++nSec;
            rp4M=s4M;
          }
          else if(rpA>2 && rpZ==0)            // Z=0 decay
          {
            theDefinition = aNeutron;
            G4LorentzVector f4M=rp4M/rpA;     // 4mom of 1st neutron
#ifdef pdebug
            G4cout<<"G4QLE::CPS->Nn,M="<<rp4M.m()<<" >? N*MNeutron="<<rpA*mNeutron<<G4endl;
#endif
            for(G4int it=1; it<rpA; ++it)     // Fill (N-1) neutrons to output
            {
              G4DynamicParticle* pNeu = new G4DynamicParticle(theDefinition, f4M);
              G4Track* aNTrack = new G4Track(pNeu, localtime, position );
              aNTrack->SetWeight(weight);                    //    weighted
              aNTrack->SetTouchableHandle(trTouchable);
              result.push_back(aNTrack);
            }
            G4int nesc = rpA-1;
            tt4M-=f4M*nesc;
#ifdef edebug
            totBaN-=nesc;
            tch4M -=f4M*nesc;
            G4cout<<">G4QLEn::PSDI:Nn,tZ="<<totChg<<",tB="<<totBaN<<",t4M="<<tch4M<<G4endl;
#endif
#ifdef pdebug
            G4cout<<"G4QLowEn::PSDI:-->ProjSpectA4M="<<f4M<<G4endl;
#endif
            nSec+=nesc;
            rp4M=f4M;
          }
          else if(rpA==8 && rpZ==4)            // Be8 decay
          {
            theDefinition = anAlpha;
            G4LorentzVector f4M=G4LorentzVector(0.,0.,0.,mAlph); // 4mom of 1st alpha
            G4LorentzVector s4M=G4LorentzVector(0.,0.,0.,mAlph); // 4mom of 2nd alpha
#ifdef pdebug
            G4cout<<"G4QLE::CPS->A+A,mBe8="<<rp4M.m()<<" >? 2*MAlpha="<<2*mAlph<<G4endl;
#endif
            if(!G4QHadron(rp4M).DecayIn2(f4M, s4M))
            {
              G4cout<<"*W*G4QLE::CPS->A+A,t="<<rp4M.m()<<" >? 2*MAlpha="<<2*mAlph<<G4endl;
            }
            G4DynamicParticle* pAl = new G4DynamicParticle(theDefinition, f4M);
            aFraPTrack = new G4Track(pAl, localtime, position );
            aFraPTrack->SetWeight(weight);                    //    weighted
            aFraPTrack->SetTouchableHandle(trTouchable);
            tt4M-=f4M;
#ifdef edebug
            totChg-=2;
            totBaN-=4;
            tch4M -=f4M;
            G4cout<<">>G4QLEn::PSDI:1,tZ="<<totChg<<",tB="<<totBaN<<",t4M="<<tch4M<<G4endl;
#endif
#ifdef pdebug
            G4cout<<"G4QLowEn::PSDI:-->ProjSpectA4M="<<f4M<<G4endl;
#endif
            ++nSec;
            rp4M=s4M;
          }
          else if(rpA==5 && (rpZ==2 || rpZ==3)) // He5/Li5 decay
          {
            theDefinition = aProton;
            G4double mNuc = mProt;
            if(rpZ==2)
            {
              theDefinition = aNeutron;
              mNuc = mNeut;
            }
            G4LorentzVector f4M=G4LorentzVector(0.,0.,0.,mNuc);  // 4mom of the nucleon
            G4LorentzVector s4M=G4LorentzVector(0.,0.,0.,mAlph); // 4mom of the alpha
#ifdef pdebug
            G4cout<<"G4QLowE::CPS->N+A, tM5="<<rp4M.m()<<" >? MA+MN="<<mAlph+mNuc<<G4endl;
#endif
            if(!G4QHadron(rp4M).DecayIn2(f4M, s4M))
            {
              G4cout<<"*W*G4QLE::CPS->N+A,t="<<rp4M.m()<<" >? MA+MN="<<mAlph+mNuc<<G4endl;
            }
            G4DynamicParticle* pAl = new G4DynamicParticle(theDefinition, f4M);
            aFraPTrack = new G4Track(pAl, localtime, position );
            aFraPTrack->SetWeight(weight);                    //    weighted
            aFraPTrack->SetTouchableHandle(trTouchable);
            tt4M-=f4M;
#ifdef edebug
            if(theDefinition == aProton) totChg-=1;
            totBaN-=1;
            tch4M -=f4M;
            G4cout<<">>G4QLEn::PSDI:2,tZ="<<totChg<<",tB="<<totBaN<<",t4M="<<tch4M<<G4endl;
#endif
#ifdef pdebug
            G4cout<<"G4QLowEn::PSDI:-->ProjSpectN4M="<<f4M<<G4endl;
#endif
            ++nSec;
            theDefinition = anAlpha;
            rp4M=s4M;
          }
          else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rpZ,rpA,0,rpZ);
          if(!theDefinition) //ANDREA: Add raise of exception
	    {
	      G4cout<<"-Warning-G4QLowEn::PSDI: pDef=0, Z="<<rpZ<<", A="<<rpA<<G4endl;
	      throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
	    } 
#ifdef edebug
          if(theDefinition == anAlpha)
          {
            totChg-=2;
            totBaN-=4;
          }
          else
          {
            totChg-=rpZ;
            totBaN-=rpA;
          }
          tch4M -=rp4M;
          G4cout<<">>G4QLEn::PSDI:3, tZ="<<totChg<<",tB="<<totBaN<<", t4M="<<tch4M<<G4endl;
#endif
          G4DynamicParticle* rpD = new G4DynamicParticle(theDefinition, rp4M);
          G4Track* aNewPTrack = new G4Track(rpD, localtime, position);
          aNewPTrack->SetWeight(weight);//    weighted
          aNewPTrack->SetTouchableHandle(trTouchable);
          tt4M-=rp4M;
#ifdef pdebug
          G4cout<<"G4QLowEn::PSDI:-->ProjSpectR4M="<<rp4M<<",Z="<<rpZ<<",A="<<rpA<<G4endl;
#endif
          ++nSec;
          //
          // Split the target spectator
          G4int rtZ=tZ-tC;            // number of protons in the target spectator
          G4int tF=tD-tC;             // number of neutrons in the target part of comp
          G4int rtN=tN-tF;            // number of neutrons in the target spectator
          G4double rtNM=G4QPDGCode(2112).GetNuclMass(rtZ, rtN, 0); // Mass of the spectator
          G4double rtE=std::sqrt(rtNM*rtNM+tFM.mag2());
          G4LorentzVector rt4M(tFM,rtE);
          G4int rtA=rtZ+rtN;
          G4Track* aFraTTrack = 0;
          theDefinition = 0;
          if(rtA==1)
          {
            if(!rtZ) theDefinition = G4Neutron::Neutron();
            else     theDefinition = G4Proton::Proton();
#ifdef pdebug
            G4cout<<"G4QLE::PSDI: rtA=1, rtZ"<<rtZ<<G4endl;
#endif
          }
          //AND 19JAN2010->: handling of Z=0,A>=2 case (see lines 949 and following
          else if(rtA==2 && rtZ==0)            // nn decay
          {
            theDefinition = aNeutron;
            G4LorentzVector f4M=G4LorentzVector(0.,0.,0.,mNeut); // 4mom of 1st neutron
            G4LorentzVector s4M=G4LorentzVector(0.,0.,0.,mNeut); // 4mom of 2nd neutron
#ifdef pdebug
            G4cout<<"G4QLE::CPS->n+n,nn="<<rptM.m()<<" >? 2*MNeutron="<<2*mNeutron<<G4endl;
#endif
            if(!G4QHadron(rt4M).DecayIn2(f4M, s4M))
            {
              G4cout<<"*W*G4QLE::CPS->n+n,t="<<rt4M.m()<<" >? 2*Neutron="<<2*mAlph<<G4endl;
            }
            G4DynamicParticle* pNeu = new G4DynamicParticle(theDefinition, f4M);
            aFraPTrack = new G4Track(pNeu, localtime, position );
            aFraPTrack->SetWeight(weight);                    //    weighted
            aFraPTrack->SetTouchableHandle(trTouchable);
            tt4M-=f4M;
#ifdef edebug
            totBaN-=2;
            tch4M -=f4M;
            G4cout<<">>G4QLEn::PSDI:n,tZ="<<totChg<<",tB="<<totBaN<<",t4M="<<tch4M<<G4endl;
#endif
#ifdef pdebug
            G4cout<<"G4QLowEn::PSDI:-->ProjSpectA4M="<<f4M<<G4endl;
#endif
            ++nSec;
            rt4M=s4M;
          }
          else if(rtA>2 && rtZ==0)            // Z=0 decay
          {
            theDefinition = aNeutron;
            G4LorentzVector f4M=rt4M/rtA;     // 4mom of 1st neutron
#ifdef pdebug
            G4cout<<"G4QLE::CPS->Nn,M="<<rt4M.m()<<" >? N*MNeutron="<<rtA*mNeutron<<G4endl;
#endif
            for(G4int it=1; it<rtA; ++it)     // Fill (N-1) neutrons to output
            {
              G4DynamicParticle* pNeu = new G4DynamicParticle(theDefinition, f4M);
              G4Track* aNTrack = new G4Track(pNeu, localtime, position );
              aNTrack->SetWeight(weight);                    //    weighted
              aNTrack->SetTouchableHandle(trTouchable);
              result.push_back(aNTrack);
            }
            G4int nesc = rtA-1;
            tt4M-=f4M*nesc;
#ifdef edebug
            totBaN-=nesc;
            tch4M -=f4M*nesc;
            G4cout<<">G4QLEn::PSDI:Nn,tZ="<<totChg<<",tB="<<totBaN<<",t4M="<<tch4M<<G4endl;
#endif
#ifdef pdebug
            G4cout<<"G4QLowEn::PSDI:-->ProjSpectA4M="<<f4M<<G4endl;
#endif
            nSec+=nesc;
            rt4M=f4M;
          }
          //AND 19JAN2010<-
          else if(rtA==8 && rtZ==4)            // Be8 decay
          {
            theDefinition = anAlpha;
            G4LorentzVector f4M=G4LorentzVector(0.,0.,0.,mAlph); // 4mom of 1st alpha
            G4LorentzVector s4M=G4LorentzVector(0.,0.,0.,mAlph); // 4mom of 2nd alpha
#ifdef pdebug
            G4cout<<"G4QLE::CPS->A+A,mBe8="<<rt4M.m()<<" >? 2*MAlpha="<<2*mAlph<<G4endl;
#endif
            if(!G4QHadron(rt4M).DecayIn2(f4M, s4M))
            {
              G4cout<<"*W*G4QLE::CTS->A+A,t="<<rt4M.m()<<" >? 2*MAlpha="<<2*mAlph<<G4endl;
            }
            G4DynamicParticle* pAl = new G4DynamicParticle(theDefinition, f4M);
            aFraTTrack = new G4Track(pAl, localtime, position );
            aFraTTrack->SetWeight(weight);                        // weighted
            aFraTTrack->SetTouchableHandle(trTouchable);
            tt4M-=f4M;
#ifdef edebug
            totChg-=2;
            totBaN-=4;
            tch4M -=f4M;
            G4cout<<">>G4QLEn::PSDI:4,tZ="<<totChg<<",tB="<<totBaN<<",t4M="<<tch4M<<G4endl;
#endif
#ifdef pdebug
            G4cout<<"G4QLowEn::PSDI:-->TargSpectA4M="<<f4M<<G4endl;
#endif
            ++nSec;
            rt4M=s4M;
          }
          else if(rtA==5 && (rtZ==2 || rtZ==3)) // He5/Li5 decay
          {
            theDefinition = aProton;
            G4double mNuc = mProt;
            if(rpZ==2)
            {
              theDefinition = aNeutron;
              mNuc = mNeut;
            }
            G4LorentzVector f4M=G4LorentzVector(0.,0.,0.,mNuc);  // 4mom of the nucleon
            G4LorentzVector s4M=G4LorentzVector(0.,0.,0.,mAlph); // 4mom of the alpha
#ifdef pdebug
            G4cout<<"G4QLowE::CPS->N+A, tM5="<<rt4M.m()<<" >? MA+MN="<<mAlph+mNuc<<G4endl;
#endif
            if(!G4QHadron(rt4M).DecayIn2(f4M, s4M))
            {
              G4cout<<"*W*G4QLE::CTS->N+A,t="<<rt4M.m()<<" >? MA+MN="<<mAlph+mNuc<<G4endl;
            }
            G4DynamicParticle* pAl = new G4DynamicParticle(theDefinition, f4M);
            aFraTTrack = new G4Track(pAl, localtime, position );
            aFraTTrack->SetWeight(weight);                        // weighted
            aFraTTrack->SetTouchableHandle(trTouchable);
            tt4M-=f4M;
#ifdef edebug
            if(theDefinition == aProton) totChg-=1;
            totBaN-=1;
            tch4M -=f4M;
            G4cout<<">>G4QLEn::PSDI:5,tZ="<<totChg<<",tB="<<totBaN<<",t4M="<<tch4M<<G4endl;
#endif
#ifdef pdebug
            G4cout<<"G4QLowEn::PSDI:-->TargSpectN4M="<<f4M<<G4endl;
#endif
            ++nSec;
            theDefinition = anAlpha;
            rt4M=s4M;
          }
          else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rtZ,rtA,0,rtZ);
          if(!theDefinition) //ANDREA Adding throw of exception
	  {
            G4cout<<"-Warning-G4QLowEn::PSDI: tDef=0, Z="<<rtZ<<", A="<<rtA<<G4endl;
	    throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
	  }
#ifdef edebug
          if(theDefinition == anAlpha)
          {
            totChg-=2;
            totBaN-=4;
          }
          else
          {
            totChg-=rtZ;
            totBaN-=rtA;
          }
          tch4M -=rt4M;
          G4cout<<">>G4QLEn::PSDI:6, tZ="<<totChg<<",tB="<<totBaN<<", t4M="<<tch4M<<G4endl;
#endif
          G4DynamicParticle* rtD = new G4DynamicParticle(theDefinition, rt4M);
          G4Track* aNewTTrack = new G4Track(rtD, localtime, position);
          aNewTTrack->SetWeight(weight);                         // weighted
          aNewTTrack->SetTouchableHandle(trTouchable);
          tt4M-=rt4M;
#ifdef pdebug
          G4cout<<"G4QLowEn::PSDI:-->TargSpectR4M="<<rt4M<<",Z="<<rtZ<<",A="<<rtA<<G4endl;
#endif
          ++nSec;
          nSec+=3;
          aParticleChange.SetNumberOfSecondaries(nSec); 
          if(aFraPTrack) result.push_back(aFraPTrack);
          if(aNewPTrack) result.push_back(aNewPTrack);
          if(aFraTTrack) result.push_back(aFraTTrack);
          if(aNewTTrack) result.push_back(aNewTTrack);
          tcM = tt4M.m();            // total CMS mass of the compound
          G4int sN=tF+pF;
          G4int sZ=tC+pC;
          tnM = targQPDG.GetNuclMass(sZ,sN,0); // GSM
#ifdef pdebug
          G4cout<<"G4QLEn::PSDI:At#"<<nAt<<",totM="<<tcM<<",gsM="<<tnM<<",Z="<<sZ<<",N="
                <<sN<<G4endl;
#endif
          if(tcM>tnM)
          {
            pZ=pC;
            pN=pF;
            tZ=tC;
            tN=tF;
            tot4M=tt4M;
          }
        }
      }                               // At least one is splitted
      else if(tD==tB || pD==pB)       // Total absorption
      {
        tD=tZ+tN;
        pD=pZ+pN;
        tcM=tnM+1.;
      }
    }
    else                              // Total compound (define tD to go out of the while
    {
      tD=tZ+tN;
      pD=pZ+pN;
     tcM=tnM+1.;
    }
  } // End of the interaction WHILE
  G4double totM=tot4M.m();            // total CMS mass of the reaction
  G4int totN=tN+pN;
  G4int totZ=tZ+pZ;
#ifdef edebug
  G4cout<<">>>G4QLEn::PSDI: dZ="<<totChg-totZ<<", dB="<<totBaN-totN-totZ<<", d4M="
        <<tch4M-tot4M<<",N="<<totN<<",Z="<<totZ<<G4endl;
#endif
  // @@ Here mass[i] can be calculated if mass=0
  G4double mass[nCh]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                      0.,0.,0.,0.,0.,0.};
  mass[0] = tM = targQPDG.GetNuclMass(totZ,totN,0); // gamma+gamma
#ifdef pdebug
  G4cout<<"G4QLEn::PSDI:TotM="<<totM<<",NucM="<<tM<<",totZ="<<totZ<<",totN="<<totN<<G4endl;
#endif
  if (totN>0 && totZ>0)
  {
    mass[1] = targQPDG.GetNuclMass(totZ,totN-1,0); // gamma+neutron
    mass[2] = targQPDG.GetNuclMass(totZ-1,totN,0); // gamma+proton
  }
  if ( totZ > 1 && totN > 1 ) mass[3] = targQPDG.GetNuclMass(totZ-1,totN-1,0); //g+d
  if ( totZ > 1 && totN > 2 ) mass[4] = targQPDG.GetNuclMass(totZ-1,totN-2,0); //g+t
  if ( totZ > 2 && totN > 1 ) mass[5] = targQPDG.GetNuclMass(totZ-2,totN-1,0); //g+3
  if ( totZ > 2 && totN > 2 ) mass[6] = targQPDG.GetNuclMass(totZ-2,totN-2,0); //g+a
  if ( totZ > 0 && totN > 2 ) mass[7] = targQPDG.GetNuclMass(totZ  ,totN-2,0); //n+n
    mass[ 8] = mass[3]; // neutron+proton (the same as a deuteron)
  if ( totZ > 2 )             mass[9] = targQPDG.GetNuclMass(totZ-2,totN  ,0); //p+p
    mass[10] = mass[5]; // proton+deuteron (the same as He3)
    mass[11] = mass[4]; // neutron+deuteron (the same as t)
    mass[12] = mass[6]; // deuteron+deuteron (the same as alpha)
    mass[13] = mass[6]; // proton+tritium (the same as alpha)
  if ( totZ > 1 && totN > 3 ) mass[14] = targQPDG.GetNuclMass(totZ-1,totN-3,0);//n+t
  if ( totZ > 3 && totN > 1 ) mass[15] = targQPDG.GetNuclMass(totZ-3,totN-1,0);//He3+p
    mass[16] = mass[6]; // neutron+He3 (the same as alpha)
  if ( totZ > 3 && totN > 2 ) mass[17] = targQPDG.GetNuclMass(totZ-3,totN-2,0);//pa
  if ( totZ > 2 && totN > 3 ) mass[18] = targQPDG.GetNuclMass(totZ-2,totN-3,0);//na
  if(pL>0) // @@ Not debugged @@
  {
    G4int pL1=pL-1;
    if(totN>0||totZ>0) mass[19] = targQPDG.GetNuclMass(totZ  ,totN  ,pL1);// Lambda+gamma
    if( (totN > 0 && totZ > 0) ||        totZ > 1 ) 
      mass[20]=targQPDG.GetNuclMass(totZ-1,totN  ,pL1);//Lp
    if( (totN > 0 && totZ > 0) || totN > 0        ) 
      mass[21]=targQPDG.GetNuclMass(totZ  ,totN-1,pL1);//Ln
    if( (totN > 1 && totZ > 0) || (totN > 0 && totZ > 1) ) 
      mass[22]=targQPDG.GetNuclMass(totZ-1,totN-1,pL1);//Ld
    if( (totN > 2 && totZ > 0) || (totN > 1 && totZ > 1) ) 
      mass[23]=targQPDG.GetNuclMass(totZ-1,totN-2,pL1);//Lt
    if( (totN > 0 && totZ > 2) || (totN > 1 && totZ > 1) ) 
      mass[24]=targQPDG.GetNuclMass(totZ-2,totN-1,pL1);//L3
    if( (totN > 1 && totZ > 2) || (totN > 2 && totZ > 1) ) 
      mass[25]=targQPDG.GetNuclMass(totZ-2,totN-2,pL1);//La
  }
#ifdef debug
  G4cout<<"G4QLowEn::PSDI: Residual masses are calculated"<<G4endl;
#endif
  tA=tZ+tN; 
  pA=pZ+pN; 
  G4double prZ=pZ/pA+tZ/tA;
  G4double prN=pN/pA+tN/tA;
  G4double prD=prN*prZ;
  G4double prA=prD*prD;
  G4double prH=prD*prZ;
  G4double prT=prD*prN;
  G4double fhe3=6.*std::sqrt(tA);
  G4double prL=0.;
  if(pL>0) prL=pL/pA;
  G4double qval[nCh]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                      0.,0.,0.,0.,0.,0.};
  qval[ 0] = (totM - mass[ 0])/137./137.;
  qval[ 1] = (totM - mass[ 1] - mNeut)*prN/137.;
  qval[ 2] = (totM - mass[ 2] - mProt)*prZ/137.;
  qval[ 3] = (totM - mass[ 3] - mDeut)*prD/3./137.;
  qval[ 4] = (totM - mass[ 4] - mTrit)*prT/6./137.;
  qval[ 5] = (totM - mass[ 5] - mHel3)*prH/fhe3/137.;
  qval[ 6] = (totM - mass[ 6] - mAlph)*prA/9./137.;
  qval[ 7] = (totM - mass[ 7] - mNeut - mNeut)*prN*prN;
  qval[ 8] = (totM - mass[ 8] - mNeut - mProt)*prD;
  qval[ 9] = (totM - mass[ 9] - mProt - mProt)*prZ*prZ;
  qval[10] = (totM - mass[10] - mProt - mDeut)*prH/3.;
  qval[11] = (totM - mass[11] - mNeut - mDeut)*prT/3.;
  qval[12] = (totM - mass[12] - mDeut - mDeut)*prA/3./3.;
  qval[13] = (totM - mass[13] - mProt - mTrit)*prA/6.;
  qval[14] = (totM - mass[14] - mNeut - mTrit)*prT*prN/6.;
  qval[15] = (totM - mass[15] - mProt - mHel3)*prH*prZ/fhe3;
  qval[16] = (totM - mass[16] - mNeut - mHel3)*prA/fhe3;
  qval[17] = (totM - mass[17] - mProt - mAlph)*prZ*prA/9.;
  qval[18] = (totM - mass[18] - mNeut - mAlph)*prN*prA/9.;
  if(pZ>0)
  {
    qval[19] = (totM - mass[19] - mLamb)*prL;
    qval[20] = (totM - mass[20] - mProt - mLamb)*prL*prZ;
    qval[21] = (totM - mass[21] - mNeut - mLamb)*prL*prN;
    qval[22] = (totM - mass[22] - mDeut - mLamb)*prL*prD/2.;
    qval[23] = (totM - mass[23] - mTrit - mLamb)*prL*prT/3.;
    qval[24] = (totM - mass[24] - mHel3 - mLamb)*prL*prH/fhe3;
    qval[25] = (totM - mass[25] - mAlph - mLamb)*prL*prA/4;
  }
#ifdef debug
  G4cout<<"G4QLowEn::PSDI: Q-values are calculated, tgA="<<tA<<", prA="<<pA<<G4endl;
#endif
  
  G4double qv = 0.0;                        // Total sum of probabilities (q-values)
  for(G4int i=0; i<nCh; ++i )
  {
#ifdef sdebug
    G4cout<<"G4QLowEn::PSDI: i="<<i<<", q="<<qval[i]<<",m="<<mass[i]<<G4endl;
#endif
    if( mass[i] < 500.*MeV ) qval[i] = 0.0; // Close A/Z impossible channels
    if( qval[i] < 0.0 )      qval[i] = 0.0; // Close the splitting impossible channels
    qv += qval[i];
  }
  // Select the channel
  G4double qv1 = 0.0;
  G4double ran = G4UniformRand();
  G4int index  = 0;
  for( index=0; index<nCh; ++index )
  {
    if( qval[index] > 0.0 )
    {
      qv1 += qval[index]/qv;
      if( ran <= qv1 ) break;
    }
  }
#ifdef debug
  G4cout<<"G4QLowEn::PSDI: index="<<index<<" < "<<nCh<<G4endl;
#endif
  if(index == nCh)
  {
    G4cout<<"***G4QLowEnergy::PoStDI:Decay is impossible,totM="<<totM<<",GSM="<<tM<<G4endl;
    throw G4QException("G4QLowEnergy::PostStepDoIt: Can't decay the Compound");
  }
  // @@ Convert it to G4QHadrons    
  G4DynamicParticle* ResSec = new G4DynamicParticle;
  G4DynamicParticle* FstSec = new G4DynamicParticle;
  G4DynamicParticle* SecSec = new G4DynamicParticle;
#ifdef debug
  G4cout<<"G4QLowEn::PSDI: Dynamic particles are created pL="<<pL<<G4endl;
#endif

  G4LorentzVector res4Mom(zeroMom,mass[index]*MeV); // The recoil nucleus prototype
  G4double mF=0.;
  G4double mS=0.;
  G4int    rA=totZ+totN;
  G4int    rZ=totZ;
  G4int    rL=pL;
  G4int complete=3;
  G4ParticleDefinition* theDefinition;  // Prototype for qfNucleon
  switch( index )
  {
   case 0:
     if(!evaporate || rA<2)
     {
       if(!rZ) theDefinition=aNeutron;
       else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
       if(!theDefinition) //ANDREA: Adding throw of exception
       {
         G4cerr<<"-Warning-G4LE::PSDI: notDef(1), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
	 throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
       }
       ResSec->SetDefinition( theDefinition );
       FstSec->SetDefinition( aGamma );
       SecSec->SetDefinition( aGamma );
     }
     else
     {
       delete ResSec;
       delete FstSec;
       delete SecSec;
       complete=0;
     }
     break;
   case 1:
     rA-=1;                                              // gamma+n
     if(!evaporate || rA<2)
     {
       if(!rZ) theDefinition=aNeutron;
       else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
       if(!theDefinition) //ANDREA: Adding throw of an excpetion
       {
         G4cerr<<"-Warning-G4LE::PSDI: notDef(2), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
	 throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
       }
       ResSec->SetDefinition( theDefinition );
       SecSec->SetDefinition( aGamma );
     }
     else
     {
       delete ResSec;
       delete SecSec;
       complete=1;
     }
     FstSec->SetDefinition( aNeutron );
     mF=mNeut; // First hadron 4-momentum
     break;
   case 2:
     rA-=1;
     rZ-=1;                                             // gamma+p
     if(!evaporate || rA<2)
     {
       if(!rZ) theDefinition=aNeutron;
       else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
       if(!theDefinition) //ANDREA: Adding throw of an exception
       {
         G4cerr<<"-Warning-G4LE::PSDI: notDef(3), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
	 throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
       }
       ResSec->SetDefinition( theDefinition );
       SecSec->SetDefinition( aGamma );
     }
     else
     {
       delete ResSec;
       delete SecSec;
       complete=1;
     }
     FstSec->SetDefinition( aProton );
     mF=mProt; // First hadron 4-momentum
     break;
   case 3:
     rA-=2;
     rZ-=1;                                             // gamma+d
     if(!evaporate || rA<2)
     {
       if(!rZ) theDefinition=aNeutron;
       else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
       if(!theDefinition) //ANDREA: Adding throw of an exception
       {
         G4cerr<<"-Warning-G4LE::PSDI: notDef(4), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
	 throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
       }
       ResSec->SetDefinition( theDefinition );
       SecSec->SetDefinition( aGamma );
     }
     else
     {
       delete ResSec;
       delete SecSec;
       complete=1;
     }
     FstSec->SetDefinition( aDeuteron );
     mF=mDeut; // First hadron 4-momentum
     break;
   case 4:
     rA-=3;                                             // gamma+t
     rZ-=1;
     if(!evaporate || rA<2)
     {
       if(!rZ) theDefinition=aNeutron;
       else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
       if(!theDefinition) //ANDREA: Adding throw of an exception
       {
         G4cerr<<"-Warning-G4LE::PSDI: notDef(5), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
	 throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
       }
       ResSec->SetDefinition( theDefinition );
       SecSec->SetDefinition( aGamma );
     }
     else
     {
       delete ResSec;
       delete SecSec;
       complete=1;
     }
     FstSec->SetDefinition( aTriton );
     mF=mTrit; // First hadron 4-momentum
     break;
  case 5:                                            // gamma+He3
     rA-=3;
     rZ-=2;
     if(!evaporate || rA<2)
     {
       if(!rZ) theDefinition=aNeutron;
       else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
       if(!theDefinition) //ANDREA: Adding throw of an excetpion
       {
         G4cerr<<"-Warning-G4LE::PSDI: notDef(6), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
	 throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
       }
       ResSec->SetDefinition( theDefinition );
       SecSec->SetDefinition( aGamma );
     }
     else
     {
       delete ResSec;
       delete SecSec;
       complete=1;
     }
     FstSec->SetDefinition( aHe3);
     mF=mHel3; // First hadron 4-momentum
     break;
   case 6:
     rA-=4;
     rZ-=2;                                         // gamma+He4
     if(!evaporate || rA<2)
     {
       if(!rZ) theDefinition=aNeutron;
       else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
       if(!theDefinition) //ANDREA: Adding throw of an excetpion
       {
         G4cerr<<"-Warning-G4LE::PSDI: notDef(7), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
	 throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
       }
       ResSec->SetDefinition( theDefinition );
       SecSec->SetDefinition( aGamma );
     }
     else
     {
       delete ResSec;
       delete SecSec;
       complete=1;
     }
     FstSec->SetDefinition( anAlpha );
     mF=mAlph; // First hadron 4-momentum
     break;
   case 7:
     rA-=2;                                          // n+n
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition) //ANDREA: Adding thorw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(8), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
       throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aNeutron );
     SecSec->SetDefinition( aNeutron );
     mF=mNeut; // First hadron 4-momentum
     mS=mNeut; // Second hadron 4-momentum
     break;
   case 8:
     rZ-=1;
     rA-=2;                                           // n+p
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else if(rA==2 && !rZ)
     {
       index=7;
       ResSec->SetDefinition( aDeuteron);
       FstSec->SetDefinition( aNeutron );
       SecSec->SetDefinition( aNeutron );
       mF=mNeut; // First hadron 4-momentum
       mS=mNeut; // Second hadron 4-momentum
       break;
     }
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition) //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(9), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
       throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aNeutron );
     SecSec->SetDefinition( aProton );
     mF=mNeut; // First hadron 4-momentum
     mS=mProt; // Second hadron 4-momentum
     break;
   case 9:
     rZ-=2;
     rA-=2;                                           // p+p
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition) //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(10), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
       throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aProton );
     SecSec->SetDefinition( aProton );
     mF=mProt; // First hadron 4-momentum
     mS=mProt; // Second hadron 4-momentum
     break;
   case 10:
     rZ-=2;
     rA-=3;                                            // p+d
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition)//ANDREA: Adding throw of an excpetion
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(11), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
       throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aProton );
     SecSec->SetDefinition( aDeuteron );
     mF=mProt; // First hadron 4-momentum
     mS=mDeut; // Second hadron 4-momentum
     break;
   case 11:
     rZ-=1;
     rA-=3;                                            // n+d
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition) //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(12), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
       throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aNeutron );
     SecSec->SetDefinition( aDeuteron );
     mF=mNeut; // First hadron 4-momentum
     mS=mDeut; // Second hadron 4-momentum
     break;
   case 12:
     rZ-=2;
     rA-=4;                                            // d+d
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition) //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(13), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
       throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aDeuteron );
     SecSec->SetDefinition( aDeuteron );
     mF=mDeut; // First hadron 4-momentum
     mS=mDeut; // Second hadron 4-momentum
     break;
   case 13:
     rZ-=2;
     rA-=4;                                            // p+t
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition) //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(14), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
       throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aProton );
     SecSec->SetDefinition( aTriton );
     mF=mProt; // First hadron 4-momentum
     mS=mTrit; // Second hadron 4-momentum
     break;
   case 14:
     rZ-=1;
     rA-=4;                                             // n+t
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition) //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(15), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
       throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aNeutron );
     SecSec->SetDefinition( aTriton );
     mF=mNeut; // First hadron 4-momentum
     mS=mTrit; // Second hadron 4-momentum
     break;
   case 15:
     rZ-=3;
     rA-=4;                                             // p+He3
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition)  //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(16), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
       throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aProton);
     SecSec->SetDefinition( aHe3 );
     mF=mProt; // First hadron 4-momentum
     mS=mHel3; // Second hadron 4-momentum
     break;
   case 16:
     rZ-=2;
     rA-=4;                                             // n+He3
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition)   //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(17), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
       throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aNeutron );
     SecSec->SetDefinition( aHe3 );
     mF=mNeut; // First hadron 4-momentum
     mS=mHel3; // Second hadron 4-momentum
     break;
   case 17:
     rZ-=3;
     rA-=5;                                            // p+alph
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition)   //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(18), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
       throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aProton );
     SecSec->SetDefinition( anAlpha );
     mF=mProt; // First hadron 4-momentum
     mS=mAlph; // Second hadron 4-momentum
     break;
   case 18:
     rZ-=2;
     rA-=5;                                              // n+alph
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition)   //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(19), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
      throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aNeutron );
     SecSec->SetDefinition( anAlpha );
     mF=mNeut; // First hadron 4-momentum
     mS=mAlph; // Second hadron 4-momentum
     break;
   case 19:
     rL-=1;                                              // L+gamma (@@ rA inludes rL?)
     rA-=1;
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition)   //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(20), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
      throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aLambda );
     SecSec->SetDefinition( aGamma );
     mF=mLamb; // First hadron 4-momentum
     break;
   case 20:
     rL-=1;                                              // L+p (@@ rA inludes rL?)
     rZ-=1;
     rA-=2;
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition)   //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(21), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
      throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aProton );
     SecSec->SetDefinition( aLambda );
     mF=mProt; // First hadron 4-momentum
     mS=mLamb; // Second hadron 4-momentum
     break;
   case 21:
     rL-=1;                                              // L+n (@@ rA inludes rL?)
     rA-=2;
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition)   //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(22), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
      throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aNeutron );
     SecSec->SetDefinition( aLambda );
     mF=mNeut; // First hadron 4-momentum
     mS=mLamb; // Second hadron 4-momentum
     break;
   case 22:
     rL-=1;                                              // L+d (@@ rA inludes rL?)
     rZ-=1;
     rA-=3;
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition)   //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(23), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
       throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
    ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aDeuteron );
     SecSec->SetDefinition( aLambda );
     mF=mDeut; // First hadron 4-momentum
     mS=mLamb; // Second hadron 4-momentum
     break;
   case 23:
     rL-=1;                                              // L+t (@@ rA inludes rL?)
     rZ-=1;
     rA-=4;
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition)   //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(24), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
       throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
    ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aTriton );
     SecSec->SetDefinition( aLambda );
     mF=mTrit; // First hadron 4-momentum
     mS=mLamb; // Second hadron 4-momentum
     break;
   case 24:
     rL-=1;                                              // L+He3 (@@ rA inludes rL?)
     rZ-=2;
     rA-=4;
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition)   //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(25), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
      throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( aHe3 );
     SecSec->SetDefinition( aLambda );
     mF=mHel3; // First hadron 4-momentum
     mS=mLamb; // Second hadron 4-momentum
     break;
   case 25:
     rL-=1;                                              // L+alph (@@ rA inludes rL?)
     rZ-=2;
     rA-=5;
     if(rA==1 && !rZ) theDefinition=aNeutron;
     else theDefinition=G4ParticleTable::GetParticleTable()->FindIon(rZ,rA,rL,rZ);
     if(!theDefinition)   //ANDREA: Adding throw of an exception
     {
       G4cerr<<"-Warning-G4LE::PSDI: notDef(26), Z="<<rZ<<", A="<<rA<<", L="<<rL<<G4endl;
       throw G4QException("G4QLowEnergy::PostStepDoIt particle definition is a null pointer");
     }
     ResSec->SetDefinition( theDefinition );
     FstSec->SetDefinition( anAlpha );
     SecSec->SetDefinition( aLambda );
     mF=mAlph; // First hadron 4-momentum
     mS=mLamb; // Second hadron 4-momentum
     break;
  }
#ifdef debug
  G4cout<<"G4QLowEn::PSDI:F="<<mF<<",S="<<mS<<",com="<<complete<<",ev="<<evaporate<<G4endl;
#endif
  G4LorentzVector fst4Mom(zeroMom,mF); // Prototype of the first hadron 4-momentum
  G4LorentzVector snd4Mom(zeroMom,mS); // Prototype of the second hadron 4-momentum
  G4LorentzVector dir4Mom=tot4M;       // Prototype of the resN decay direction 4-momentum
  dir4Mom.setE(tot4M.e()/2.);          // Get half energy and total 3-momentum
  // @@ Can be repeated to take into account the Coulomb Barrier
  if(!G4QHadron(tot4M).CopDecayIn3(fst4Mom,snd4Mom,res4Mom,dir4Mom,cosp))
  {                   //                                          
    G4cerr<<"**G4LowEnergy::PoStDoIt:i="<<index<<",tM="<<totM<<"->M1="<<res4Mom.m()<<"+M2="
     <<fst4Mom.m()<<"+M3="<<snd4Mom.m()<<"=="<<res4Mom.m()+fst4Mom.m()+snd4Mom.m()<<G4endl;
    throw G4QException("G4QLowEnergy::PostStepDoIt: Can't decay the Compound");
  }                   //                                          
#ifdef debug
  G4cout<<"G4QLowEn::PSDI:r4M="<<res4Mom<<",f4M="<<fst4Mom<<",s4M="<<snd4Mom<<G4endl;
#endif
  G4Track* aNewTrack = 0;
  if(complete)
  {
    FstSec->Set4Momentum(fst4Mom);
    aNewTrack = new G4Track(FstSec, localtime, position );
    aNewTrack->SetWeight(weight);                                   //    weighted
    aNewTrack->SetTouchableHandle(trTouchable);
    result.push_back( aNewTrack );
    EnMomConservation-=fst4Mom;
#ifdef debug
    G4cout<<"G4QLowEn::PSDI: ***Filled*** 1stH4M="<<fst4Mom
          <<", PDG="<<FstSec->GetDefinition()->GetPDGEncoding()<<G4endl;
#endif
    if(complete>2)                     // Final solution
    {
      ResSec->Set4Momentum(res4Mom);
      aNewTrack = new G4Track(ResSec, localtime, position );
      aNewTrack->SetWeight(weight);                                   //    weighted
      aNewTrack->SetTouchableHandle(trTouchable);
      result.push_back( aNewTrack );
      EnMomConservation-=res4Mom;
#ifdef debug
      G4cout<<"G4QLowEn::PSDI: ***Filled*** rA4M="<<res4Mom<<",rZ="<<rZ<<",rA="<<rA<<",rL="
            <<rL<<G4endl;
#endif
      SecSec->Set4Momentum(snd4Mom);
      aNewTrack = new G4Track(SecSec, localtime, position );
      aNewTrack->SetWeight(weight);                                   //    weighted
      aNewTrack->SetTouchableHandle(trTouchable);
      result.push_back( aNewTrack );
      EnMomConservation-=snd4Mom;
#ifdef debug
      G4cout<<"G4QLowEn::PSDI: ***Filled*** 2ndH4M="<<snd4Mom
          <<", PDG="<<SecSec->GetDefinition()->GetPDGEncoding()<<G4endl;
#endif
    }
    else res4Mom+=snd4Mom;
  }
  else res4Mom=tot4M;
  if(complete<3)                        // Evaporation of the residual must be done
  {
    G4QHadron* rHadron = new G4QHadron(90000000+999*rZ+rA,res4Mom); // Input hadron-nucleus
    G4QHadronVector* evaHV = new G4QHadronVector; // Output vector of hadrons (delete!)
    Nuc.EvaporateNucleus(rHadron, evaHV); // here a pion can appear !
    G4int nOut=evaHV->size();
    for(G4int i=0; i<nOut; i++)
    {
      G4QHadron* curH = (*evaHV)[i];
      G4int hPDG=curH->GetPDGCode();
      G4LorentzVector h4Mom=curH->Get4Momentum();
      EnMomConservation-=h4Mom;
#ifdef debug
      G4cout<<"G4QLowEn::PSDI: ***FillingCand#"<<i<<"*** evaH="<<hPDG<<h4Mom<<G4endl;
#endif
      if     (hPDG==90000001 || hPDG==2112) theDefinition = aNeutron;
      else if(hPDG==90001000 || hPDG==2212) theDefinition = aProton;
      else if(hPDG==91000000 || hPDG==3122) theDefinition = aLambda;
      else if(hPDG==     22 )               theDefinition = aGamma;
      else if(hPDG==     111)               theDefinition = aPiZero;
      else if(hPDG==     211)               theDefinition = aPiPlus;
      else if(hPDG==    -211)               theDefinition = aPiMinus;
      else
      {
        G4int hZ=curH->GetCharge();
        G4int hA=curH->GetBaryonNumber();
        G4int hS=curH->GetStrangeness();
        theDefinition = G4ParticleTable::GetParticleTable()->FindIon(hZ,hA,hS,0); // ion
      }
      if(theDefinition)
      {
        G4DynamicParticle* theEQH = new G4DynamicParticle(theDefinition,h4Mom);
        G4Track* evaQH = new G4Track(theEQH, localtime, position );
        evaQH->SetWeight(weight);                                   //    weighted
        evaQH->SetTouchableHandle(trTouchable);
        result.push_back( evaQH );         
      }
      else G4cerr<<"-Warning-G4QLowEnergy::PostStepDoIt: Bad secondary PDG="<<hPDG<<G4endl;
    }
  }
  // algorithm implementation --- STOPS HERE
  G4int nres=result.size();
  aParticleChange.SetNumberOfSecondaries(nres);
  for(G4int i=0; i<nres; ++i) aParticleChange.AddSecondary(result[i]);
#ifdef debug
  G4cout<<"G4QLowEnergy::PostStepDoIt:*** PostStepDoIt is done ***"<<G4endl;
#endif
  return G4VDiscreteProcess::PostStepDoIt(track, step);
}

G4double G4QLowEnergy::CalculateXS(G4double p, G4int Z, G4int N, G4int PDG) 
{
  static G4bool first=true;
  static G4VQCrossSection* CSmanager;
  if(first)                              // Connection with a singletone
  {
    CSmanager=G4QIonIonCrossSection::GetPointer();
    if(PDG == 2212) CSmanager=G4QProtonNuclearCrossSection::GetPointer();
    first=false;
  }
#ifdef debug
  G4cout<<"G4QLowE::CXS: *DONE* p="<<p<<",Z="<<Z<<",N="<<N<<",PDG="<<PDG<<G4endl;
#endif
  return CSmanager->GetCrossSection(true, p, Z, N, PDG);
}
