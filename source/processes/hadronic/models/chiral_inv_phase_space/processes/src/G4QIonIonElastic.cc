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
//      ---------------- G4QIonIonElastic class -----------------
//                 by Mikhail Kossov, December 2006.
// G4QIonIonElastic class of the CHIPS Simulation Branch in GEANT4
// ---------------------------------------------------------------
// ****************************************************************************************
// ********** This CLASS is temporary moved from the photolepton_hadron directory *********
// ****************************************************************************************
// Short description: a simple process for the Ion-Ion elastic scattering.
// For heavy by heavy ions it can reach 50% of the total cross-section.
// -----------------------------------------------------------------------

//#define debug
//#define pdebug
//#define tdebug
//#define nandebug
//#define ppdebug

#include "G4QIonIonElastic.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronicDeprecate.hh"


// Initialization of static vectors
//G4int G4QIonIonElastic::nPartCWorld=152;   // The#of particles initialized in CHIPS World
//G4int G4QIonIonElastic::nPartCWorld=122;   // The#of particles initialized in CHIPS World
G4int G4QIonIonElastic::nPartCWorld=85; // The#of particles initialized in CHIPS World Red.
std::vector<G4int> G4QIonIonElastic::ElementZ;        // Z of the element(i) in theLastCalc
std::vector<G4double> G4QIonIonElastic::ElProbInMat;  // SumProbabilityElements in Material
std::vector<std::vector<G4int>*> G4QIonIonElastic::ElIsoN; // N of isotope(j) of Element(i)
std::vector<std::vector<G4double>*>G4QIonIonElastic::IsoProbInEl;//SumProbabIsotopes in ElI

// Constructor
G4QIonIonElastic::G4QIonIonElastic(const G4String& processName):
  G4VDiscreteProcess(processName, fHadronic)
{
  G4HadronicDeprecate("G4QIonIonElastic");

#ifdef debug
  G4cout<<"G4QIonIonElastic::Constructor is called processName="<<processName<<G4endl;
#endif
  if (verboseLevel>0) G4cout << GetProcessName() << " process is created "<< G4endl;
  //G4QCHIPSWorld::Get()->GetParticles(nPartCWorld); // Create CHIPS World (234 part. max)
}

// Destructor
G4QIonIonElastic::~G4QIonIonElastic() {}

// output of the function must be in units of length! L=1/sig_V,sig_V=SUM(n(j,i)*sig(j,i)),
// where n(i,j) is a number of nuclei of the isotop j of the element i in V=1(lengtUnit^3)
// ********** All CHIPS cross sections are calculated in the surface units ************
G4double G4QIonIonElastic::GetMeanFreePath(const G4Track& aTrack, G4double,
                                           G4ForceCondition* Fc)
{
  static const G4double mProt = G4QPDGCode(2212).GetMass()/MeV;// CHIPS proton Mass in MeV
  *Fc = NotForced;
  const G4DynamicParticle* incidentParticle = aTrack.GetDynamicParticle();
  G4ParticleDefinition* incidentParticleDefinition=incidentParticle->GetDefinition();
  if( !IsApplicable(*incidentParticleDefinition))
    G4cout<<"-Warning-G4QIonIonElastic::GetMeanFreePath: notImplementedParticle"<<G4endl;
  // Calculate the mean Cross Section for the set of Elements(*Isotopes) in the Material
  G4double Momentum = incidentParticle->GetTotalMomentum(); // 3-momentum of the Particle
#ifdef debug
  G4double KinEn = incidentParticle->GetKineticEnergy();
  G4cout<<"G4QIonIonElastic::GetMeanFreePath: kinE="<<KinEn<<",Mom="<<Momentum<<G4endl;
#endif
  const G4Material* material = aTrack.GetMaterial();        // Get the current material
  const G4double* NOfNucPerVolume = material->GetVecNbOfAtomsPerVolume();
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QIonIonElastic::GetMeanFreePath:"<<nE<<" Elem's in theMaterial"<<G4endl;
#endif
  G4VQCrossSection* CSmanager=G4QIonIonCrossSection::GetPointer();
  G4VQCrossSection* PCSmanager=G4QProtonElasticCrossSection::GetPointer();
  G4int pPDG=0;
  // Probably enough: pPDG=incidentParticleDefinition->GetPDGEncoding();
  G4int pZ=static_cast<G4int>(incidentParticleDefinition->GetPDGCharge());
  G4int pA=incidentParticleDefinition->GetBaryonNumber();
  if      (incidentParticleDefinition ==  G4Deuteron::Deuteron()     ) pPDG = 1000010020;
  else if (incidentParticleDefinition ==  G4Alpha::Alpha()           ) pPDG = 1000020040;
  else if (incidentParticleDefinition ==  G4Triton::Triton()         ) pPDG = 1000010030;
  else if (incidentParticleDefinition ==  G4He3::He3()               ) pPDG = 1000020030;
  else if (incidentParticleDefinition ==  G4GenericIon::GenericIon() || (pZ > 0 && pA > 0))
  {
    pPDG=incidentParticleDefinition->GetPDGEncoding();
#ifdef debug
    G4int prPDG=1000000000+10000*pZ+10*pA;
    G4cout<<"G4QIonIonElastic::GetMeanFreePath: PDG="<<prPDG<<"="<<pPDG<<G4endl;
#endif
  }
  else G4cout<<"-Warning-G4QIonIonElastic::GetMeanFreePath:Unknown projectile Ion"<<G4endl;
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
    G4cout<<"G4QIonIonElastic::GetMeanFreePath: isovectorLength="<<isoSize<<G4endl;
#endif
    if(isoSize)                             // The Element has non-trivial abundance set
    {
      indEl=pElement->GetIndex()+1;         // Index of the non-trivial element is an order
#ifdef debug
      G4cout<<"G4QIIEl::GetMFP:iE="<<indEl<<", def="<<Isotopes->IsDefined(Z,indEl)<<G4endl;
#endif
      if(!Isotopes->IsDefined(Z,indEl))     // This index is not defined for this Z: define
      {
        std::vector<std::pair<G4int,G4double>*>* newAbund =
                                               new std::vector<std::pair<G4int,G4double>*>;
        G4double* abuVector=pElement->GetRelativeAbundanceVector();
        for(G4int j=0; j<isoSize; j++)      // Calculation of abundance vector for isotopes
        {
          G4int N=pElement->GetIsotope(j)->GetN()-Z; // N means A=N+Z !
          if(pElement->GetIsotope(j)->GetZ()!=Z) G4cerr<<"G4QIonIonEl::GetMeanFreePath Z="
                                         <<pElement->GetIsotope(j)->GetZ()<<"#"<<Z<<G4endl;
          G4double abund=abuVector[j];
          std::pair<G4int,G4double>* pr= new std::pair<G4int,G4double>(N,abund);
#ifdef debug
          G4cout<<"G4QIonIonElastic::GetMeanFP:pair#="<<j<<",N="<<N<<",ab="<<abund<<G4endl;
#endif
          newAbund->push_back(pr);
        }
#ifdef debug
        G4cout<<"G4QIonIonElastic::GetMeanFP: pairVectorLength="<<newAbund->size()<<G4endl;
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
    G4cout<<"G4QIonIonEl::GetMFP:=***=> #isot="<<nIs<<", Z="<<Z<<", indEl="<<indEl<<G4endl;
#endif
    G4double susi=0.;                       // sum of CS over isotopes
    if(nIs) for(G4int j=0; j<nIs; j++)      // Calculate CS for eachIsotope of El
    {
      std::pair<G4int,G4double>* curIs=(*cs)[j]; // A pointer, which is used twice
      G4int N=curIs->first;                 // #of Neuterons in the isotope j of El i
      IsN->push_back(N);                    // Remember Min N for the Element
#ifdef debug
      G4cout<<"G4QIIE::GMFP:true,P="<<Momentum<<",Z="<<Z<<",N="<<N<<",pPDG="<<pPDG<<G4endl;
#endif
      G4bool ccsf=false;                    // Extract elastic Ion-Ion cross-section
#ifdef debug
      G4cout<<"G4QIonIonElastic::GMFP: GetCS #1 j="<<j<<G4endl;
#endif
      G4double CSI=0.;
      if(Z==1 && !N)                        // Proton target. Reversed kinematics.
      {
        G4double projM=G4QPDGCode(2212).GetNuclMass(pZ,pA-pZ,0); // Mass of the projNucleus
        CSI=PCSmanager->GetCrossSection(true, Momentum*mProt/projM, Z, N, 2212); // Ap CS
      }
      else CSI=CSmanager->GetCrossSection(ccsf,Momentum,Z,N,pPDG); // CS(j,i) for isotope
#ifdef debug
      G4cout<<"G4QIonIonElastic::GMFP: jI="<<j<<", Zt="<<Z<<", Nt="<<N<<", Mom="<<Momentum
            <<", XSec="<<CSI/millibarn<<G4endl;
#endif
      curIs->second = CSI;
      susi+=CSI;                            // Make a sum per isotopes
      SPI->push_back(susi);                 // Remember summed cross-section
    } // End of temporary initialization of the cross sections in the G4QIsotope singeltone
    sigma+=Isotopes->GetMeanCrossSection(Z,indEl)*NOfNucPerVolume[i];//SUM(MeanCS*NOfNperV)
#ifdef debug
    G4cout<<"G4QIonIonEl::GMFP:<S>="<<Isotopes->GetMeanCrossSection(Z,indEl)<<",AddToSig="
          <<Isotopes->GetMeanCrossSection(Z,indEl)*NOfNucPerVolume[i]<<G4endl;
#endif
    ElProbInMat.push_back(sigma);
  } // End of LOOP over Elements
  // Check that cross section is not zero and return the mean free path
#ifdef debug
  G4cout<<"G4QIonIonElastic::GetMeanFreePath: MeanFreePath="<<1./sigma<<G4endl;
#endif
  if(sigma > 0.00000000001) return 1./sigma;                 // Mean path [distance] 
  return DBL_MAX;
}


G4bool G4QIonIonElastic::IsApplicable(const G4ParticleDefinition& particle) 
{
  G4int Z=static_cast<G4int>(particle.GetPDGCharge());
  G4int A=particle.GetBaryonNumber();
  if      (particle == *( G4Deuteron::Deuteron()     )) return true;
  else if (particle == *( G4Alpha::Alpha()           )) return true;
  else if (particle == *( G4Triton::Triton()         )) return true;
  else if (particle == *( G4He3::He3()               )) return true;
  else if (particle == *( G4GenericIon::GenericIon() )) return true;
  else if (Z > 0 && A > 0)                              return true;
#ifdef debug
  G4cout<<"***>>G4QIonIonElastic::IsApplicable: PDG="<<particle.GetPDGEncoding()<<", A="
        <<A<<", Z="<<Z<<G4endl;
#endif
  return false;
}

G4VParticleChange* G4QIonIonElastic::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  static const G4double mProt= G4QPDGCode(2212).GetMass(); // CHIPS proton mass in MeV
  static const G4double fm2MeV2 = 3*38938./1.09; // (3/1.09)*(hc)^2 in fm^2*MeV^2
  static G4bool CWinit = true;                   // CHIPS Warld needs to be initted
  if(CWinit)
  {
    CWinit=false;
    G4QCHIPSWorld::Get()->GetParticles(nPartCWorld); // Create CHIPS World (234 part.max)
  }
  //-------------------------------------------------------------------------------------
  const G4DynamicParticle* projHadron = track.GetDynamicParticle();
  const G4ParticleDefinition* particle=projHadron->GetDefinition();
#ifdef debug
  G4cout<<"G4QIonIonElastic::PostStepDoIt: Before the GetMeanFreePath is called In4M="
        <<projHadron->Get4Momentum()<<" of PDG="<<particle->GetPDGEncoding()<<", Type="
        <<particle->GetParticleType()<<", Subtp="<<particle->GetParticleSubType()<<G4endl;
#endif
  G4LorentzVector proj4M=(projHadron->Get4Momentum())/MeV; // Convert to MeV!
  G4double momentum = projHadron->GetTotalMomentum()/MeV; // 3-momentum of the Proj in MeV
  G4double Momentum = proj4M.rho();                   // @@ Just for the test purposes
  if(std::fabs(Momentum-momentum)>.000001)
    G4cerr<<"-Warn-G4QIonIonElastic::PostStepDoIt:P(IU)="<<Momentum<<"#"<<momentum<<G4endl;
  G4double pM2=proj4M.m2();        // in MeV^2
  G4double pM=std::sqrt(pM2);      // in MeV
#ifdef pdebug
  G4cout<<"G4QIonIonEl::PostStepDoIt: pP(IU)="<<Momentum<<"="<<momentum<<",pM="<<pM<<G4endl;
#endif
  if (!IsApplicable(*particle))  // Check applicability
  {
    G4cerr<<"G4QIonIonElastic::PostStepDoIt: Only NA elastic is implemented."<<G4endl;
    return 0;
  }
  const G4Material* material = track.GetMaterial();      // Get the current material
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QIonIonElastic::PostStepDoIt: "<<nE<<" elements in the material."<<G4endl;
#endif
  // Probably enough: projPDG=particle->GetPDGEncoding();
  G4int projPDG=0;                           // CHIPS PDG Code for the captured hadron
  G4int pZ=static_cast<G4int>(particle->GetPDGCharge());
  G4int pA=particle->GetBaryonNumber();
  if      (particle ==  G4Deuteron::Deuteron() ) projPDG= 1000010020;
  else if (particle ==  G4Alpha::Alpha()       ) projPDG= 1000020040;
  else if (particle ==  G4Triton::Triton()     ) projPDG= 1000010030;
  else if (particle ==  G4He3::He3()           ) projPDG= 1000020030;
  else if (particle ==  G4GenericIon::GenericIon() || (pZ > 0 && pA > 0))
  {
    projPDG=particle->GetPDGEncoding();
#ifdef debug
    G4int prPDG=1000000000+10000*pZ+10*pA;
    G4cout<<"G4QIonIonElastic::PostStepDoIt: PDG="<<prPDG<<"="<<projPDG<<G4endl;
#endif
  }
  else G4cout<<"-Warning-G4QIonIonElastic::PostStepDoIt:Unknown projectile Ion"<<G4endl;
#ifdef debug
  G4int prPDG=particle->GetPDGEncoding();
  G4cout<<"G4QIonIonElastic::PostStepDoIt: projPDG="<<projPDG<<", stPDG="<<prPDG<<G4endl;
#endif
  if(!projPDG)
  {
    G4cerr<<"-Warning-G4QIonIonElastic::PostStepDoIt:Undefined interactingNucleus"<<G4endl;
    return 0;
  }
  G4int pN=pA-pZ;                           // Projectile N
  G4int EPIM=ElProbInMat.size();
#ifdef debug
  G4cout<<"G4QIonIonElastic::PSDI:m="<<EPIM<<",n="<<nE<<",T="<<ElProbInMat[EPIM-1]<<G4endl;
#endif
  G4int i=0;
  if(EPIM>1)
  {
    G4double rnd = ElProbInMat[EPIM-1]*G4UniformRand();
    for(i=0; i<nE; ++i)
    {
#ifdef debug
      G4cout<<"G4QIonIonElastic::PSDI: EPM["<<i<<"]="<<ElProbInMat[i]<<", r="<<rnd<<G4endl;
#endif
      if (rnd<ElProbInMat[i]) break;
    }
    if(i>=nE) i=nE-1;                        // Top limit for the Element
  }
  G4Element* pElement=(*theElementVector)[i];
  G4int tZ=static_cast<G4int>(pElement->GetZ());
#ifdef debug
    G4cout<<"G4QIonIonElastic::PostStepDoIt: i="<<i<<", Z(element)="<<tZ<<G4endl;
#endif
  if(tZ<=0)
  {
    G4cerr<<"---Warning---G4QIonIonElastic::PostStepDoIt: Element with Z="<<tZ<<G4endl;
    if(tZ<0) return 0;
  }
  std::vector<G4double>* SPI = IsoProbInEl[i];// Vector of summedProbabilities for isotopes
  std::vector<G4int>* IsN = ElIsoN[i];     // Vector of "#of neutrons" in the isotope El[i]
  G4int nofIsot=SPI->size();               // #of isotopes in the element i
#ifdef debug
  G4cout<<"G4QIonIonElastic::PosStDoIt: nI="<<nofIsot<<",T="<<(*SPI)[nofIsot-1]<<G4endl;
#endif
  G4int j=0;
  if(nofIsot>1)
  {
    G4double rndI=(*SPI)[nofIsot-1]*G4UniformRand(); // Randomize the isotop of the Element
    for(j=0; j<nofIsot; ++j)
    {
#ifdef debug
      G4cout<<"G4QIonIonElastic::PostStDI: SP["<<j<<"]="<<(*SPI)[j]<<", r="<<rndI<<G4endl;
#endif
      if(rndI < (*SPI)[j]) break;
    }
    if(j>=nofIsot) j=nofIsot-1;            // Top limit for the isotope
  }
  G4int tN =(*IsN)[j]; ;                   // Randomized number of neutrons
#ifdef debug
  G4cout<<"G4QIonIonElastic::PostStepDoIt:j="<<i<<",N(isotope)="<<tN<<",MeV="<<MeV<<G4endl;
#endif
  if(tN<0)
  {
    G4cerr<<"-Warning-G4QIonIonElastic::PostStepDoIt:IsotopeZ="<<tZ<<" & 0>N="<<tN<<G4endl;
    return 0;
  }
  nOfNeutrons=tN;                           // Remember it for the energy-momentum check
#ifdef debug
  G4cout<<"G4QIonIonElastic::PostStepDoIt: N="<<tN<<" for element with Z="<<tZ<<G4endl;
#endif
  aParticleChange.Initialize(track);
#ifdef debug
  G4cout<<"G4QIonIonElastic::PostStepDoIt: track is initialized"<<G4endl;
#endif
  G4double weight        = track.GetWeight();
  G4double localtime     = track.GetGlobalTime();
  G4ThreeVector position = track.GetPosition();
#ifdef debug
  G4cout<<"G4QIonIonElastic::PostStepDoIt: before Touchable extraction"<<G4endl;
#endif
  G4TouchableHandle trTouchable = track.GetTouchableHandle();
#ifdef debug
  G4cout<<"G4QIonIonElastic::PostStepDoIt: Touchable is extracted"<<G4endl;
#endif
  // Definitions for do nothing
  G4double kinEnergy= projHadron->GetKineticEnergy()*MeV; // Kin energy in MeV (Is *MeV n?)
  G4ParticleMomentum dir = projHadron->GetMomentumDirection();// It is a unit three-vector
  // !! Exception for the proton target !!
  G4bool revkin=false;
  G4ThreeVector bvel(0.,0.,0.);
  G4int tA=tZ+tN;
  G4int targPDG=90000000+tZ*1000+tN;       // CHIPS PDG Code of the target nucleus
  G4QPDGCode targQPDG(targPDG);            // @@ one can use G4Ion & get rid of CHIPS World
  G4double tM=targQPDG.GetMass();          // CHIPS target mass in MeV
  G4LorentzVector targ4M(0.,0.,0.,tM);
  if(tZ==1 && !tN)                         // Update parameters in DB and trans to Antilab
  {
    G4ForceCondition cond=NotForced;
    GetMeanFreePath(track, -27., &cond);   // @@ ?? jus to update parameters?
#ifdef debug
    G4cout<<"G4QIonIonElastic::PostStepDoIt: After the GetMeanFreePath is called"<<G4endl;
#endif
    revkin=true;
    tZ=pZ;
    tN=pN;
    tA=tZ+tN;
    tM=pM;
    pZ=1;                                  // @@ Is that necessary ??
    pN=0;                                  // @@ Is that necessary ??
    pA=1;                                  // @@ Is that necessary ??
    pM=mProt;
    bvel=proj4M.vect()/proj4M.e();         // Lab->Antilab transition boost velocity
    proj4M=targ4M.boost(-bvel);            // Proton 4-mom in Antilab
    targ4M=G4LorentzVector(0.,0.,0.,tM);   // Projectile nucleus 4-mom in Antilab
    Momentum = proj4M.rho();               // Recalculate Momentum in Antilab
  }
  G4LorentzVector tot4M=proj4M+targ4M;     // Total 4-mom of the reaction
#ifdef debug
  G4cout<<"G4QIonIonElastic::PostStDI: tM="<<tM<<", p4M="<<proj4M<<", t4M="<<tot4M<<G4endl;
#endif
  EnMomConservation=tot4M;                 // Total 4-mom of reaction for E/M conservation
  G4VQCrossSection* PELmanager=G4QProtonElasticCrossSection::GetPointer();
  G4VQCrossSection* NELmanager=G4QNeutronElasticCrossSection::GetPointer();
  G4VQCrossSection* CSmanager=G4QIonIonCrossSection::GetPointer();
  // @@ Probably this is not necessary any more
#ifdef debug
  G4cout<<"G4QIIEl::PSDI:f, P="<<Momentum<<",Z="<<tZ<<",N="<<tN<<",tPDG="<<projPDG<<G4endl;
#endif
  // false means elastic cross-section
  G4double xSec=0.;                           // Proto of Recalculated Cross Section
  if(revkin) xSec=PELmanager->GetCrossSection(false, Momentum, tZ, tN, 2212);
  else       xSec=CSmanager->GetCrossSection(false, Momentum, tZ, tN, projPDG);
#ifdef debug
  G4cout<<"G4QIIEl::PSDI: pPDG="<<projPDG<<",P="<<Momentum<<",CS="<<xSec/millibarn<<G4endl;
#endif
#ifdef nandebug
  if(xSec>0. || xSec<0. || xSec==0);
  else  G4cout<<"-NaN-Warning-G4QIonIonElastic::PostStDoIt: xSec="<<xSec/millibarn<<G4endl;
#endif
  // @@ check a possibility to separate p, n, or alpha (!)
  if(xSec <= 0.) // The cross-section iz 0 -> Do Nothing
  {
#ifdef pdebug
    G4cerr<<"-Warning-G4QIonIonElastic::PostStDoIt: *Zero cross-section* PDG="<<projPDG
          <<",tPDG="<<targPDG<<",P="<<Momentum<<G4endl;
#endif
    //Do Nothing Action insead of the reaction
    aParticleChange.ProposeEnergy(kinEnergy);
    aParticleChange.ProposeLocalEnergyDeposit(0.);
    aParticleChange.ProposeMomentumDirection(dir) ;
    return G4VDiscreteProcess::PostStepDoIt(track,step);
  }
  G4double mint=0;
  G4double maxt=0;
  G4double dtM=tM+tM;
  if(revkin)
  {
    mint=PELmanager->GetExchangeT(tZ,tN,2212); // functional randomized -t in MeV^2
    maxt=PELmanager->GetHMaxT();
  }
  else
  {
    G4double PA=Momentum*pA;
    G4double PA2=PA*PA;
    maxt=dtM*PA2/(std::sqrt(PA2+pM2)+tM/2+pM2/dtM);
#ifdef pdebug
    G4cout<<"G4QIonIonElastic::PostStDoIt:pPDG="<<projPDG<<",tPDG="<<targPDG<<",P="
          <<Momentum<<",CS="<<xSec<<",maxt="<<maxt<<G4endl;
#endif
    xSec=PELmanager->GetCrossSection(false, Momentum, 1, 0, 2212);// pp=nn
    G4double B1=PELmanager->GetSlope(1,0,2212); // slope for pp=nn
    xSec=NELmanager->GetCrossSection(false, Momentum, 1, 0, 2112);// np=pn
    G4double B2 =NELmanager->GetSlope(1,0,2112); // slope for np=pn
    G4double mB =((pZ*tZ+pN*tN)*B1+(pZ*tN+pN*tZ)*B2)/(pA+tA);
    G4double pR2=std::pow(pA+4.,.305)/fm2MeV2;
    G4double tR2=std::pow(tA+4.,.305)/fm2MeV2;
    G4double eB =mB+pR2+tR2;
    mint=-std::log(1.-G4UniformRand()*(1.-std::exp(-eB*maxt)))/eB;
    mint+=mint;
#ifdef pdebug
    G4cout<<"G4QIonIonElastic::PostStDoIt:B1="<<B1<<",B2="<<B2<<",mB="<<mB
          <<",pR2="<<pR2<<",tR2="<<tR2<<",eB="<<eB<<",mint="<<mint<<G4endl;
#endif
  }
#ifdef nandebug
  if(mint>-.0000001);
  else  G4cout<<"-Warning-G4QIonIonElastic::PostStDoIt:-t="<<mint<<G4endl;
#endif
  G4double cost=1.-mint/maxt; // cos(theta) in CMS
  // 
#ifdef ppdebug
  G4cout<<"G4QIonIonElastic::PoStDoI:t="<<mint<<",dpcm2="<<CSmanager->GetHMaxT()<<",Ek="
        <<kinEnergy<<",tM="<<tM<<",pM="<<pM<<",cost="<<cost<<G4endl;
#endif
  if(cost>1. || cost<-1. || !(cost>-1. || cost<1.))
  {
    if(cost>1.000001 || cost<-1.000001 || !(cost>-1. || cost<1.))
    {
      G4double tM2=tM*tM;                         // Squared target mass
      G4double pEn=pM+kinEnergy;                  // tot projectile Energy in MeV
      G4double sM=dtM*pEn+tM2+pM2;                // Mondelstam s
      G4double twop2cm=(tM2+tM2)*(pEn*pEn-pM2)/sM;// Max_t/2 (2*p^2_cm)
      G4cout<<"-Warning-G4QIonIonElastic::PoStDI:cos="<<cost<<",t="<<mint<<",T="<<kinEnergy
            <<",tM="<<tM<<",tmax="<<2*kinEnergy*tM<<",p="<<projPDG<<",t="<<targPDG<<G4endl;
      G4cout<<"G4QIonIonElastic::PSDI:dpcm2="<<twop2cm<<"="<<CSmanager->GetHMaxT()<<G4endl;
    }
    if     (cost>1.)  cost=1.;
    else if(cost<-1.) cost=-1.;
  }
  G4LorentzVector scat4M(0.,0.,0.,pM);            // 4-mom of the scattered projectile
  G4LorentzVector reco4M(0.,0.,0.,tM);            // 4-mom of the recoil target
  G4LorentzVector dir4M=tot4M-G4LorentzVector(0.,0.,0.,(tot4M.e()-tM-pM)*.01);
  if(!G4QHadron(tot4M).RelDecayIn2(scat4M, reco4M, dir4M, cost, cost))
  {
    G4cout<<"-Warning-G4QIonIonE::PSDI:t4M="<<tot4M<<",pM="<<pM<<",tM="<<tM<<",cost="
          <<cost<<G4endl;
  }
  if(revkin)
  {
    G4LorentzVector tmpLV=scat4M.boost(bvel);     // Recoil target Back to LS
    scat4M=reco4M.boost(bvel);                    // Scattered Progectile Back to LS
    reco4M=tmpLV;
    pM=tM;
    tM=mProt;
  }
#ifdef debug
  G4cout<<"G4QIonIonElast::PSDI:s4M="<<scat4M<<"+r4M="<<reco4M<<"="<<scat4M+reco4M<<G4endl;
  G4cout<<"G4QIonIonElastic::PSDI:scatE="<<scat4M.e()-pM<<",recoE="<<reco4M.e()-tM<<",d4M="
        <<tot4M-scat4M-reco4M<<G4endl;
#endif
  // Update G4VParticleChange for the scattered projectile
  G4double finE=scat4M.e()-pM;             // Final kinetic energy of the scattered proton
  if(finE>0.0) aParticleChange.ProposeEnergy(finE);
  else
  {
    if(finE<-1.e-8 || !(finE>-1.||finE<1.)) // NAN or negative
      G4cerr<<"*Warning*G4QIonIonElastic::PostStDoIt: Zero or negative scattered E="<<finE
            <<", s4M="<<scat4M<<", r4M="<<reco4M<<", d4M="<<tot4M-scat4M-reco4M<<G4endl;
    aParticleChange.ProposeEnergy(0.) ;
    aParticleChange.ProposeTrackStatus(fStopButAlive);
  }
  G4ThreeVector findir=scat4M.vect()/scat4M.rho();  // Unit vector in new direction
  aParticleChange.ProposeMomentumDirection(findir); // new direction for the scattered part
  EnMomConservation-=scat4M;                        // It must be initialized by (pE+tM,pP)
  // This is how in general the secondary should be identified
  G4DynamicParticle* theSec = new G4DynamicParticle; // A secondary for the recoil hadron 
#ifdef pdebug
  G4cout<<"G4QIonIonElastic::PostStepDoIt: Ion tZ="<<tZ<<", tA="<<tA<<G4endl;
#endif
  G4ParticleDefinition* theDefinition=0;
  if(revkin) theDefinition = G4Proton::Proton();
  else       theDefinition = G4ParticleTable::GetParticleTable()->FindIon(tZ,tA,0,tZ);
  if(!theDefinition)G4cout<<"-Warning-G4QIonIonElastic::PoStDI:drop PDG="<<targPDG<<G4endl;
#ifdef pdebug
  G4cout<<"G4QIonIonElastic::PoStDI:RecoilName="<<theDefinition->GetParticleName()<<G4endl;
#endif
  theSec->SetDefinition(theDefinition);
  EnMomConservation-=reco4M;
#ifdef tdebug
  G4cout<<"G4QIonIonElastic::PSD:"<<targPDG<<reco4M<<reco4M.m()<<EnMomConservation<<G4endl;
#endif
  theSec->Set4Momentum(reco4M);
#ifdef debug
  G4ThreeVector curD=theSec->GetMomentumDirection();
  G4double curM=theSec->GetMass();
  G4double curE=theSec->GetKineticEnergy()+curM;
  G4cout<<"G4QIonIonElastic::PSDI: p="<<curD<<curD.mag()<<",e="<<curE<<",m="<<curM<<G4endl;
#endif
  // Make a recoil nucleus
  G4Track* aNewTrack = new G4Track(theSec, localtime, position );
  aNewTrack->SetWeight(weight);                                   //    weighted
  aNewTrack->SetTouchableHandle(trTouchable);
  aParticleChange.AddSecondary( aNewTrack );
#ifdef debug
    G4cout<<"G4QIonIonElastic::PostStepDoIt: **** PostStepDoIt is done ****"<<G4endl;
#endif
  return G4VDiscreteProcess::PostStepDoIt(track, step);
}
