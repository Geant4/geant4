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
// $Id: G4QDiffraction.cc,v 1.1 2009-11-17 10:36:55 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//      ---------------- G4QDiffraction class -----------------
//                 by Mikhail Kossov, Aug 2007.
// G4QDiffraction class of the CHIPS Simulation Branch in GEANT4
// ---------------------------------------------------------------
// Short description: This is a process, which describes the diffraction
// excitation of the projectile and the nucleus. On nuclei in addition there
// can be a coherent diffraction process for the projectile, but it is
// comparably small. The most important part of the diffraction is the
// progectile diffraction excitation, as in this interaction proton can lose
// only a small part of its energy and make the shower longer. This is because
// only 1-2 (n) pions are produce in the diffraction escitation, and the mean
// kept energy of the nucleon is (1-n/7)=80%. For kaons the kept energy is much
// smaller (1-n/3.5)=60%, and for pions it is less important (about 40%).
// ----------------------------------------------------------------------------

//#define debug
//#define pdebug
//#define tdebug
//#define nandebug
//#define ppdebug

#include "G4QDiffraction.hh"

// Initialization of static vectors
G4int    G4QDiffraction::nPartCWorld=152;  // #of particles initialized in CHIPS
std::vector<G4int> G4QDiffraction::ElementZ; // Z of element(i) in theLastCalc
std::vector<G4double> G4QDiffraction::ElProbInMat; // SumProbOfElem in Material
std::vector<std::vector<G4int>*> G4QDiffraction::ElIsoN;// N of isotope(j), E(i)
std::vector<std::vector<G4double>*>G4QDiffraction::IsoProbInEl;//SumProbIsotE(i)

// Constructor
G4QDiffraction::G4QDiffraction(const G4String& processName):
 G4VDiscreteProcess(processName, fHadronic)
{
#ifdef debug
  G4cout<<"G4QDiffraction::Constructor is called processName="<<processName<<G4endl;
#endif
  if (verboseLevel>0) G4cout << GetProcessName() << " process is created "<< G4endl;
  G4QCHIPSWorld::Get()->GetParticles(nPartCWorld); // Create CHIPS World (234 part. max)
}

// Destructor
G4QDiffraction::~G4QDiffraction() {}


G4LorentzVector G4QDiffraction::GetEnegryMomentumConservation(){return EnMomConservation;}

G4int G4QDiffraction::GetNumberOfNeutronsInTarget() {return nOfNeutrons;}

// output of the function must be in units of length! L=1/sig_V,sig_V=SUM(n(j,i)*sig(j,i)),
// where n(i,j) is a number of nuclei of the isotop j of the element i in V=1(lengtUnit^3)
// ********** All CHIPS cross sections are calculated in the surface units ************
G4double G4QDiffraction::GetMeanFreePath(const G4Track&Track,G4double Q,G4ForceCondition*F)
{
  *F = NotForced;
  const G4DynamicParticle* incidentParticle = Track.GetDynamicParticle();
  G4ParticleDefinition* incidentParticleDefinition=incidentParticle->GetDefinition();
  if( !IsApplicable(*incidentParticleDefinition))
    G4cout<<"-Warning-G4QDiffraction::GetMeanFreePath for notImplemented Particle"<<G4endl;
  // Calculate the mean Cross Section for the set of Elements(*Isotopes) in the Material
  G4double Momentum = incidentParticle->GetTotalMomentum(); // 3-momentum of the Particle
#ifdef debug
  G4double KinEn = incidentParticle->GetKineticEnergy();
  G4cout<<"G4QDiffraction::GetMeanFreePath:Prpj, kinE="<<KinEn<<", Mom="<<Momentum<<G4endl;
#endif
  const G4Material* material = Track.GetMaterial();        // Get the current material
  const G4double* NOfNucPerVolume = material->GetVecNbOfAtomsPerVolume();
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QDiffraction::GetMeanFreePath:"<<nE<<" Elems in Material="<<*material<<G4endl;
#endif
  G4int pPDG=0;
  // @@ At present it is made only for n & p, but can be extended if inXS are available
  if     (incidentParticleDefinition == G4Proton::Proton()  ) pPDG=2212;
  else if(incidentParticleDefinition == G4Neutron::Neutron()) pPDG=2112;
  else G4cout<<"G4QDiffraction::GetMeanFreePath: only nA & pA are implemented"<<G4endl;
  
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
    G4cout<<"G4QDiffraction::GetMeanFreePath: isovector Length="<<isoSize<<G4endl;
#endif
    if(isoSize)                             // The Element has non-trivial abundance set
    {
      indEl=pElement->GetIndex()+1;         // Index of the non-trivial element is an order
#ifdef debug
      G4cout<<"G4QDiffr::GetMFP:iE="<<indEl<<",def="<<Isotopes->IsDefined(Z,indEl)<<G4endl;
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
          G4cout<<"G4QDiffract::GetMeanFreePath:pair#"<<j<<",N="<<N<<",ab="<<abund<<G4endl;
#endif
          newAbund->push_back(pr);
        }
#ifdef debug
        G4cout<<"G4QDiffract::GetMeanFreePath:pairVectorLength="<<newAbund->size()<<G4endl;
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
    G4cout<<"G4QDiffract::GetMFP:=***=>,#isot="<<nIs<<", Z="<<Z<<", indEl="<<indEl<<G4endl;
#endif
    G4double susi=0.;                       // sum of CS over isotopes
    if(nIs) for(G4int j=0; j<nIs; j++)      // Calculate CS for eachIsotope of El
    {
      std::pair<G4int,G4double>* curIs=(*cs)[j]; // A pointer, which is used twice
      G4int N=curIs->first;                 // #of Neuterons in the isotope j of El i
      IsN->push_back(N);                    // Remember Min N for the Element
#ifdef debug
      G4cout<<"G4QDiff::GMFP:true,P="<<Momentum<<",Z="<<Z<<",N="<<N<<",PDG="<<pPDG<<G4endl;
#endif
      G4bool ccsf=true;
      if(Q==-27.) ccsf=false;
#ifdef debug
      G4cout<<"G4QDiffraction::GMFP: GetCS #1 j="<<j<<G4endl;
#endif
      G4double CSI=CalculateXS(Momentum, Z, N, pPDG); // XS(j,i) for theIsotope

#ifdef debug
      G4cout<<"G4QDiffraction::GetMeanFreePath: jI="<<j<<", Zt="<<Z<<", Nt="<<N<<", Mom="
            <<Momentu<<", XSec="<<CSI/millibarn<<G4endl;
#endif
      curIs->second = CSI;
      susi+=CSI;                            // Make a sum per isotopes
      SPI->push_back(susi);                 // Remember summed cross-section
    } // End of temporary initialization of the cross sections in the G4QIsotope singeltone
    sigma+=Isotopes->GetMeanCrossSection(Z,indEl)*NOfNucPerVolume[i];//SUM(MeanCS*NOfNperV)
#ifdef debug
    G4cout<<"G4QDiffraction::GetMeanFreePath:<XS>="<<Isotopes->GetMeanCrossSection(Z,indEl)
          <<",AddSigm="<<Isotopes->GetMeanCrossSection(Z,indEl)*NOfNucPerVolume[i]<<G4endl;
#endif
    ElProbInMat.push_back(sigma);
  } // End of LOOP over Elements
  // Check that cross section is not zero and return the mean free path
#ifdef debug
  G4cout<<"G4QDiffraction::GetMeanFreePath: MeanFreePath="<<1./sigma<<G4endl;
#endif
  if(sigma > 0.) return 1./sigma;                 // Mean path [distance] 
  return DBL_MAX;
}

G4bool G4QDiffraction::IsApplicable(const G4ParticleDefinition& particle) 
{
  if      (particle == *(        G4Proton::Proton()        )) return true;
  else if (particle == *(       G4Neutron::Neutron()       )) return true;
  //else if (particle == *(     G4MuonMinus::MuonMinus()     )) return true; 
  //else if (particle == *(       G4TauPlus::TauPlus()       )) return true;
  //else if (particle == *(      G4TauMinus::TauMinus()      )) return true;
  //else if (particle == *(      G4Electron::Electron()      )) return true;
  //else if (particle == *(      G4Positron::Positron()      )) return true;
  //else if (particle == *(         G4Gamma::Gamma()         )) return true;
  //else if (particle == *(      G4MuonPlus::MuonPlus()      )) return true;
  //else if (particle == *(G4AntiNeutrinoMu::AntiNeutrinoMu())) return true;
  //else if (particle == *(    G4NeutrinoMu::NeutrinoMu()    )) return true;
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
  G4cout<<"***>>G4QDiffraction::IsApplicable: projPDG="<<particle.GetPDGEncoding()<<G4endl;
#endif
  return false;
}

G4VParticleChange* G4QDiffraction::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  static const G4double mProt= G4QPDGCode(2212).GetMass();     // CHIPS proton Mass in MeV
  static const G4double mNeut= G4QPDGCode(2112).GetMass();     // CHIPS neutron Mass in MeV
  static const G4double mPion= G4QPDGCode(111).GetMass();      // CHIPS Pi0 Mass in MeV
  static G4QDiffractionRatio* diffRatio;
  //
  //-------------------------------------------------------------------------------------
  static G4bool CWinit = true;                       // CHIPS Warld needs to be initted
  if(CWinit)
  {
    CWinit=false;
    G4QCHIPSWorld::Get()->GetParticles(nPartCWorld); // Create CHIPS World (234 part.max)
    diffRatio=G4QDiffractionRatio::GetPointer();
  }
  //-------------------------------------------------------------------------------------
  const G4DynamicParticle* projHadron = track.GetDynamicParticle();
  const G4ParticleDefinition* particle=projHadron->GetDefinition();
#ifdef debug
  G4cout<<"G4QDiffraction::PostStepDoIt: Before the GetMeanFreePath is called In4M="
        <<projHadron->Get4Momentum()<<" of PDG="<<particle->GetPDGEncoding()<<", Type="
        <<particle->GetParticleType()<<",SubType="<<particle->GetParticleSubType()<<G4endl;
#endif
  G4ForceCondition cond=NotForced;
  GetMeanFreePath(track, -27., &cond);                  // @@ ?? jus to update parameters?
#ifdef debug
  G4cout<<"G4QDiffraction::PostStepDoIt: After GetMeanFreePath is called"<<G4endl;
#endif
  G4LorentzVector proj4M=(projHadron->Get4Momentum())/MeV; // Convert to MeV!
  G4double momentum = projHadron->GetTotalMomentum()/MeV; // 3-momentum of the Proj in MeV
  G4double Momentum = proj4M.rho();                   // @@ Just for the test purposes
  if(std::fabs(Momentum-momentum)>.000001)
    G4cerr<<"-Warning-G4QDiffraction::PostStepDoIt:P_IU="<<Momentum<<"#"<<momentum<<G4endl;
#ifdef pdebug
  G4cout<<"G4QDiffraction::PostStepDoIt: pP(IU)="<<Momentum<<"="<<momentum
        <<",proj4M="<<proj4M<<", projM="<<proj4M.m()<<G4endl;
#endif
  if (!IsApplicable(*particle))  // Check applicability
  {
    G4cerr<<"G4QDiffraction::PostStepDoIt: Only NA is implemented."<<G4endl;
    return 0;
  }
  const G4Material* material = track.GetMaterial();      // Get the current material
  G4int Z=0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QDiffraction::PostStepDoIt: "<<nE<<" elements in the material."<<G4endl;
#endif
  G4int projPDG=0;                           // PDG Code prototype for the captured hadron
  // Not all these particles are implemented yet (see Is Applicable)
  if      (particle ==          G4Proton::Proton()         ) projPDG= 2212;
  else if (particle ==         G4Neutron::Neutron()        ) projPDG= 2112;
  //else if (particle ==       G4PionMinus::PionMinus()      ) projPDG= -211;
  //else if (particle ==        G4PionPlus::PionPlus()       ) projPDG=  211;
  //else if (particle ==        G4KaonPlus::KaonPlus()       ) projPDG=  321;
  //else if (particle ==       G4KaonMinus::KaonMinus()      ) projPDG= -321;
  //else if (particle ==    G4KaonZeroLong::KaonZeroLong()   ) projPDG=  130;
  //else if (particle ==   G4KaonZeroShort::KaonZeroShort()  ) projPDG=  310;
  //else if (particle ==        G4MuonPlus::MuonPlus()       ) projPDG=  -13;
  //else if (particle ==       G4MuonMinus::MuonMinus()      ) projPDG=   13;
  //else if (particle ==      G4NeutrinoMu::NeutrinoMu()     ) projPDG=   14;
  //else if (particle ==  G4AntiNeutrinoMu::AntiNeutrinoMu() ) projPDG=  -14;
  //else if (particle ==        G4Electron::Electron()       ) projPDG=   11;
  //else if (particle ==        G4Positron::Positron()       ) projPDG=  -11;
  //else if (particle ==       G4NeutrinoE::NeutrinoE()      ) projPDG=   12;
  //else if (particle ==   G4AntiNeutrinoE::AntiNeutrinoE()  ) projPDG=  -12;
  //else if (particle ==           G4Gamma::Gamma()          ) projPDG=   22;
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
#ifdef debug
  G4int prPDG=particle->GetPDGEncoding();
  G4cout<<"G4QDiffraction::PostStepDoIt: projPDG="<<projPDG<<", stPDG="<<prPDG<<G4endl;
#endif
  if(!projPDG)
  {
    G4cerr<<"-Warning-G4QDiffraction::PostStepDoIt:UndefProjHadron(PDG=0) ->ret 0"<<G4endl;
    return 0;
  }
  //G4double pM2=proj4M.m2();        // in MeV^2
  //G4double pM=std::sqrt(pM2);      // in MeV
  G4double pM=mNeut;
  G4int    fPDG=2112;
  if(projPDG==2112)
  {
    pM=mProt;
    fPDG=2212;
  }
  // Element treatment
  G4int EPIM=ElProbInMat.size();
#ifdef debug
  G4cout<<"G4QDiffra::PostStDoIt: m="<<EPIM<<",n="<<nE<<",T="<<ElProbInMat[EPIM-1]<<G4endl;
#endif
  G4int i=0;
  if(EPIM>1)
  {
    G4double rnd = ElProbInMat[EPIM-1]*G4UniformRand();
    for(i=0; i<nE; ++i)
    {
#ifdef debug
      G4cout<<"G4QDiffra::PostStepDoIt: EPM["<<i<<"]="<<ElProbInMat[i]<<",r="<<rnd<<G4endl;
#endif
      if (rnd<ElProbInMat[i]) break;
    }
    if(i>=nE) i=nE-1;                        // Top limit for the Element
  }
  G4Element* pElement=(*theElementVector)[i];
  Z=static_cast<G4int>(pElement->GetZ());
#ifdef debug
    G4cout<<"G4QDiffraction::PostStepDoIt: i="<<i<<", Z(element)="<<Z<<G4endl;
#endif
  if(Z<=0)
  {
    G4cerr<<"-Warning-G4QDiffraction::PostStepDoIt: Element with Z="<<Z<<G4endl;
    if(Z<0) return 0;
  }
  std::vector<G4double>* SPI = IsoProbInEl[i];// Vector of summedProbabilities for isotopes
  std::vector<G4int>* IsN = ElIsoN[i];     // Vector of "#of neutrons" in the isotope El[i]
  G4int nofIsot=SPI->size();               // #of isotopes in the element i
#ifdef debug
  G4cout<<"G4QDiffraction::PostStepDoIt: nI="<<nofIsot<<", T="<<(*SPI)[nofIsot-1]<<G4endl;
#endif
  G4int j=0;
  if(nofIsot>1)
  {
    G4double rndI=(*SPI)[nofIsot-1]*G4UniformRand(); // Randomize the isotop of the Element
    for(j=0; j<nofIsot; ++j)
    {
#ifdef debug
      G4cout<<"G4QDiffraction::PostStepDoIt: SP["<<j<<"]="<<(*SPI)[j]<<",r="<<rndI<<G4endl;
#endif
      if(rndI < (*SPI)[j]) break;
    }
    if(j>=nofIsot) j=nofIsot-1;            // Top limit for the isotope
  }
  G4int N =(*IsN)[j]; ;                    // Randomized number of neutrons
#ifdef debug
  G4cout<<"G4QDiffraction::PostStepDoIt: j="<<i<<", N(isotope)="<<N<<", MeV="<<MeV<<G4endl;
#endif
  if(N<0)
  {
    G4cerr<<"-Warning-G4QDiffraction::PostStepDoIt: Isotope Z="<<Z<<" has 0>N="<<N<<G4endl;
    return 0;
  }
  nOfNeutrons=N;                           // Remember it for the energy-momentum check
#ifdef debug
  G4cout<<"G4QDiffraction::PostStepDoIt: N="<<N<<" for element with Z="<<Z<<G4endl;
#endif
  if(N<0)
  {
    G4cerr<<"*Warning*G4QDiffraction::PostStepDoIt:Element with N="<<N<< G4endl;
    return 0;
  }
  aParticleChange.Initialize(track);
#ifdef debug
  G4cout<<"G4QDiffraction::PostStepDoIt: track is initialized"<<G4endl;
#endif
  G4double      weight    = track.GetWeight();
  G4double      localtime = track.GetGlobalTime();
  G4ThreeVector position  = track.GetPosition();
#ifdef debug
  G4cout<<"G4QDiffraction::PostStepDoIt: before Touchable extraction"<<G4endl;
#endif
  G4TouchableHandle trTouchable = track.GetTouchableHandle();
#ifdef debug
  G4cout<<"G4QDiffraction::PostStepDoIt: Touchable is extracted"<<G4endl;
#endif
  G4int targPDG=90000000+Z*1000+N;         // CHIPS PDG Code of the target nucleus
  G4QPDGCode targQPDG(targPDG);            // @@ use G4Ion and get rid of CHIPS World
  G4double tM=targQPDG.GetMass();          // CHIPS final nucleus mass in MeV
  G4double kinEnergy= projHadron->GetKineticEnergy()*MeV; // Kin energy in MeV (Is *MeV n?)
  G4ParticleMomentum dir = projHadron->GetMomentumDirection();// It is a unit three-vector
  G4LorentzVector tot4M=proj4M+G4LorentzVector(0.,0.,0.,tM); // Total 4-mom of the reaction
#ifdef debug
  G4cout<<"G4QDiffraction::PostStepDoIt: tM="<<tM<<",p4M="<<proj4M<<",t4M="<<tot4M<<G4endl;
#endif
  EnMomConservation=tot4M;                 // Total 4-mom of reaction for E/M conservation
  // @@ Probably this is not necessary any more
#ifdef debug
  G4cout<<"G4QDiff::PSDI:false,P="<<Momentum<<",Z="<<Z<<",N="<<N<<",PDG="<<projPDG<<G4endl;
#endif
  G4double xSec=CalculateXS(Momentum, Z, N, projPDG); // Recalculate CrossSection
#ifdef debug
  G4cout<<"G4QDiffra::PSDI:PDG="<<projPDG<<",P="<<Momentum<<",XS="<<xSec/millibarn<<G4endl;
#endif
#ifdef nandebug
  if(xSec>0. || xSec<0. || xSec==0);
  else  G4cout<<"-Warning-G4QDiffraction::PostSDI: *NAN* xSec="<<xSec/millibarn<<G4endl;
#endif
  // @@ check a possibility to separate p, n, or alpha (!)
  if(xSec <= 0.) // The cross-section iz 0 -> Do Nothing
  {
#ifdef pdebug
    G4cerr<<"*Warning*G4QDiffraction::PSDoIt:*Zero cross-section* PDG="<<projPDG
          <<",tPDG="<<targPDG<<",P="<<Momentum<<G4endl;
#endif
    //Do Nothing Action insead of the reaction
    aParticleChange.ProposeEnergy(kinEnergy);
    aParticleChange.ProposeLocalEnergyDeposit(0.);
    aParticleChange.ProposeMomentumDirection(dir) ;
    return G4VDiscreteProcess::PostStepDoIt(track,step);
  }
  G4double totCMMass=tot4M.m(); // Total CM mass, pM=projectileMass, tM=targetMass
  if(totCMMass < mPion+pM+tM) // The diffraction reaction is impossible -> Do Nothing
  {
#ifdef pdebug
    G4cerr<<"*Warning*G4QDiffraction::PSDoIt:*Below Diffraction Threshold* cmM="<<totCMMass
          <<">pM="<<pM<<"+tM="<<tM<<"+pi0="<<mPion<<"=="<<pM+tM+mPion<<G4endl;
#endif
    //Do Nothing Action insead of the reaction
    aParticleChange.ProposeEnergy(kinEnergy);
    aParticleChange.ProposeLocalEnergyDeposit(0.);
    aParticleChange.ProposeMomentumDirection(dir) ;
    return G4VDiscreteProcess::PostStepDoIt(track,step);
  }
  // Kill interacting hadron
  aParticleChange.ProposeTrackStatus(fStopAndKill);
  G4QHadronVector* out=diffRatio->TargFragment(projPDG, proj4M, Z, N);
  G4int nSec=out->size();             // #of secondaries in the diffraction reaction
  G4DynamicParticle* theSec=0;        // A prototype for secondary for the secondary
  G4LorentzVector dif4M(0.,0.,0.,0.); // Prototype for the secondary 4-momentum
  G4int difPDG=0;                     // PDG code of the secondary
  G4QHadron* difQH=0;                 // Prototype for a Q-secondary
#ifdef pdebug
  G4cout<<"G4QDiffraction::PostStepDoIt: =====found===== nSecondaries="<<nSec<<G4endl;
#endif
  for(G4int i=0; i<nSec; i++)
  {
    difQH = (*out)[i];
    difPDG= difQH->GetPDGCode();
    G4ParticleDefinition*  theDefinition=0;
    if     (difPDG==2212 || difPDG==90001000) theDefinition=G4Proton::Proton();
    else if(difPDG==2112 || difPDG==90000001) theDefinition=G4Neutron::Neutron();
    else if(difPDG==  22) theDefinition=G4Gamma::Gamma();
    else if(difPDG== 111) theDefinition=G4PionZero::PionZero();
    else if(difPDG==-211 || difPDG==89999001) theDefinition=G4PionMinus::PionMinus();
    else if(difPDG== 211 || difPDG==90000999) theDefinition=G4PionPlus::PionPlus();
    else if(difPDG== 321 || difPDG==89001000) theDefinition=G4KaonPlus::KaonPlus();
    else if(difPDG==-321 || difPDG==90999000) theDefinition=G4KaonMinus::KaonMinus();
    else if(difPDG== 130 || difPDG==-311 || difPDG==89000001)
                                              theDefinition=G4KaonZeroLong::KaonZeroLong();
    else if(difPDG== 310 || difPDG== 311 || difPDG==90999999)
                                            theDefinition=G4KaonZeroShort::KaonZeroShort();
    else if(difPDG==3122 || difPDG==91000000) theDefinition=G4Lambda::Lambda();
    else if(difPDG== 3222) theDefinition=G4SigmaPlus::SigmaPlus();
    else if(difPDG== 3112) theDefinition=G4SigmaMinus::SigmaMinus();
    else if(difPDG== 3212) theDefinition=G4SigmaZero::SigmaZero();
    else if(difPDG== 3312) theDefinition=G4XiMinus::XiMinus();
    else if(difPDG== 3322) theDefinition=G4XiZero::XiZero();
    else if(difPDG== 3334) theDefinition=G4OmegaMinus::OmegaMinus();
    else if(difPDG==-2112) theDefinition=G4AntiNeutron::AntiNeutron();
    else if(difPDG==-2212) theDefinition=G4AntiProton::AntiProton();
    else if(difPDG==-3122) theDefinition=G4AntiLambda::AntiLambda();
    else if(difPDG==-3222) theDefinition=G4AntiSigmaPlus::AntiSigmaPlus();
    else if(difPDG==-3112) theDefinition=G4AntiSigmaMinus::AntiSigmaMinus();
    else if(difPDG==-3212) theDefinition=G4AntiSigmaZero::AntiSigmaZero();
    else if(difPDG==-3312) theDefinition=G4AntiXiMinus::AntiXiMinus();
    else if(difPDG==-3322) theDefinition=G4AntiXiZero::AntiXiZero();
    else if(difPDG==-3334) theDefinition=G4AntiOmegaMinus::AntiOmegaMinus();
    else if(difPDG==  -11) theDefinition=G4Electron::Electron();
    else if(difPDG==  -13) theDefinition=G4MuonMinus::MuonMinus();
    else if(difPDG==   11) theDefinition=G4Positron::Positron();
    else if(difPDG==   13) theDefinition=G4MuonPlus::MuonPlus();
    else
    {
      G4int Z = difQH->GetCharge();
      G4int B = difQH->GetBaryonNumber();
      G4int S = difQH->GetStrangeness();
      if(S||Z>B||Z<0)G4cout<<"-Warning-G4QDif::PoStDoIt:Z="<<Z<<",A="<<B<<",S="<<S<<G4endl;
      theDefinition = G4ParticleTable::GetParticleTable()->FindIon(Z,B,0,0);
#ifdef pdebug
      G4cout<<"G4QDiffraction::PoStDoIt:Ion,Z="<<Z<<",A="<<B<<",D="<<theDefinition<<G4endl;
#endif
    }
    if(theDefinition)
    {
      theSec = new G4DynamicParticle;       // A secondary for the recoil hadron 
      theSec->SetDefinition(theDefinition);
      dif4M  = difQH->Get4Momentum();
      EnMomConservation-=dif4M;
      theSec->Set4Momentum(dif4M);
      G4Track* aNewTrack = new G4Track(theSec, localtime, position );
      aNewTrack->SetWeight(weight);                                   //    weighted
      aNewTrack->SetTouchableHandle(trTouchable);
      aParticleChange.AddSecondary( aNewTrack );
#ifdef pdebug
      G4cout<<"G4QDiffraction::PostStepDoIt: Filled 4M="<<dif4M<<", PDG="<<difPDG<<G4endl;
#endif
    }
    else G4cout<<"-Warning-G4QDif::PSDI: Lost PDG="<<difPDG<<", Z="<<difQH->GetCharge()
               <<", A="<<difQH->GetBaryonNumber()<<",S ="<<difQH->GetStrangeness()<<G4endl;
    delete difQH; // Clean up the output QHadrons
  }
  delete out;     // Delete the output QHadron-vector
#ifdef debug
  G4cout<<"G4QDiffraction::PostStepDoIt:*** PostStepDoIt is done ***"<<G4endl;
#endif
  return G4VDiscreteProcess::PostStepDoIt(track, step);
}

G4double G4QDiffraction::CalculateXS(G4double p, G4int Z, G4int N, G4int PDG) 
{
  static G4bool first=true;
  static G4VQCrossSection* CSmanager;
  static G4QDiffractionRatio* diffRatio;
  if(first)                              // Connection with a singletone
  {
    CSmanager=G4QProtonNuclearCrossSection::GetPointer();
    diffRatio=G4QDiffractionRatio::GetPointer();
    first=false;
  }
  //G4double x=CSmanager->GetCrossSection(true, p, Z, N, PDG); // inelastic XS
  //G4double pIU=p*GeV;                                        // IndependentUnistMomentum
  //G4double r=diffRatio->GetRatio(pIU, PDG, Z, N);            // Proj. Diffraction Part
  //G4double s=x*r;                                            // XS for proj. diffraction
  G4double s=diffRatio->GetTargSingDiffXS(p, PDG, Z, N);     // XS for target diffraction
#ifdef debug
  G4cout<<"G4QDiff::CXS:p="<<p<<",Z="<<Z<<",N="<<N<<",C="<<PDG<<",XS="<<s<<G4endl;
#endif
  return s;
}
