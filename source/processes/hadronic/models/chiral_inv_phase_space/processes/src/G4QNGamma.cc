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
//      ---------------- G4QNGamma class -----------------
//                 by Mikhail Kossov, December 2003.
// G4QNGamma class of the CHIPS Simulation Branch in GEANT4
// ---------------------------------------------------------------
// **************************************************************************
// This Header is a part of the CHIPS Physics Package (author: M. Kosov)
// **************************************************************************
// Short description: This is a universal class for the incoherent (inelastic)
// nuclear (n,gamma) interactions (neutron capture) in the CHIPS model.
// @@ At present the gamma cascade is not simulated (one final photon)
// ---------------------------------------------------------------------------
//#define debug
//#define pdebug

#include "G4QNGamma.hh"
#include "G4HadronicDeprecate.hh"


// Initialization of static Material/Element/Isotope vectors
std::vector<G4int> G4QNGamma::ElementZ;            // Z of the element(i) in the Last Calc
std::vector<G4double> G4QNGamma::ElProbInMat;      // ProbabilitySum ofElements inMaterial
std::vector<std::vector<G4int>*> G4QNGamma::ElIsoN;// # of isotope(j), # of Element(i)
std::vector<std::vector<G4double>*>G4QNGamma::IsoProbInEl;//SumProbabIsotopes in Element I

// Constructor
G4QNGamma::G4QNGamma(const G4String& processName)
 : G4VDiscreteProcess(processName, fHadronic)
{
  G4HadronicDeprecate("G4QNGamma");

  EnMomConservation = G4LorentzVector(0.,0.,0.,0.);
  nOfNeutrons       = 0;
#ifdef debug
  G4cout<<"G4QNGamma::Constructor is called"<<G4endl;
#endif
  if (verboseLevel>0) G4cout << GetProcessName() << " process is created "<< G4endl;
}

// Destructor (standard procedure @@ to be moved to G4VQProcess)
G4QNGamma::~G4QNGamma()
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


G4LorentzVector G4QNGamma::GetEnegryMomentumConservation() // @@ move to G4VQProcess
{
  return EnMomConservation;
}

G4int G4QNGamma::GetNumberOfNeutronsInTarget() // @@ move to G4VQProcess
{
  return nOfNeutrons;
}

// output of the function must be in units of length! L=1/sig_V,sig_V=SUM(n(j,i)*sig(j,i)),
// where n(i,j) is a number of nuclei of the isotop j of the element i in V=1(lengtUnit^3)
// ********** All CHIPS cross sections are calculated in the surface units ************
// @@ Can demand 3 internal functions when G4VQProcess is used @@ future plans
G4double G4QNGamma::GetMeanFreePath(const G4Track& aTrack,G4double,G4ForceCondition* Fc)
{
#ifdef debug
  G4cout<<"G4QNGamma::GetMeanFreePath: Called Fc="<<*Fc<<G4endl;
#endif
  *Fc = NotForced;
#ifdef debug
  G4cout<<"G4QNGamma::GetMeanFreePath: Before GetDynPart"<<G4endl;
#endif
  const G4DynamicParticle* incidentParticle = aTrack.GetDynamicParticle();
#ifdef debug
  G4cout<<"G4QNGamma::GetMeanFreePath: Before GetDef"<<G4endl;
#endif
  G4ParticleDefinition* incidentParticleDefinition=incidentParticle->GetDefinition();
  G4double Momentum = incidentParticle->GetTotalMomentum(); // 3-momentum of the Particle
  if( !IsApplicable(*incidentParticleDefinition)) // @@ Unique for all QProcesses
  {
    G4cout<<"-W-G4QNGamma::GetMeanFreePath called for not implemented particle"<<G4endl;
    return DBL_MAX;
  }
#ifdef debug
  G4cout<<"G4QNGamma::GetMeanFreePath: BeforeGetMaterial P="<<Momentum<<G4endl;
#endif
  // @@ Can be additional condition internal function of G4VQProcess
  if(Momentum > 500.) return DBL_MAX; // @@ Temporary cut (QInternal=MeV -> IU!)
  // @@ This is a standard procedure, which can be moved to G4VQProcess (above is a funct)
  const G4Material* material = aTrack.GetMaterial();        // Get the current material
  const G4double* NOfNucPerVolume = material->GetVecNbOfAtomsPerVolume();
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QNGamma::GetMeanFreePath:"<<nE<<" Elem's in theMaterial"<<G4endl;
#endif
  // @@ Can be internal function called by GetMeanFreePath (above Isotope LOOP)
  G4VQCrossSection* CSmanager    = 0;       // @@ Reference modified in the function
  G4QNeutronCaptureRatio* capMan = 0;       // @@ Reference modified in the function
  G4int pPDG     =0;                        // @@ Reference modified in the function
  G4double sigma =0.; // CS mean over isotopes @@ Reference modified in the function
  if(incidentParticleDefinition == G4Neutron::Neutron())
  {
    CSmanager=G4QNeutronNuclearCrossSection::GetPointer();
    capMan=G4QNeutronCaptureRatio::GetPointer(); // @@ can be CSmanager2
#ifdef debug
    G4cout<<"G4QNGamma::GetMeanFreePath: CSmanager is defined for neutrons"<<G4endl;
#endif
    pPDG=2112;
  }
  else
  {
    G4cout<<"-Warning-G4QNGamma::GetMeanFreePath:Particle "
          <<incidentParticleDefinition->GetPDGEncoding()<<" isn't a neutron"<<G4endl;
    return DBL_MAX;                         // can be returned in sigma
  }
  // @@ End of possible internal function
  G4QIsotope* Isotopes = G4QIsotope::Get(); // Pointer to the G4QIsotopes singleton
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
    G4cout<<"G4QNGamma::GetMeanFreePath: isovectorLength="<<isoSize<<G4endl; // Result
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
          if(pElement->GetIsotope(j)->GetZ()!=Z)G4cerr<<"G4QNGamma::GetMeanFreePath"
                                 <<": Z="<<pElement->GetIsotope(j)->GetZ()<<"#"<<Z<<G4endl;
          G4double abund=abuVector[j];
          std::pair<G4int,G4double>* pr= new std::pair<G4int,G4double>(N,abund);
#ifdef debug
          G4cout<<"G4QNGamma::GetMeanFreePath: p#="<<j<<",N="<<N<<",ab="<<abund<<G4endl;
#endif
          newAbund->push_back(pr);
        }
#ifdef debug
        G4cout<<"G4QNGamma::GetMeanFreePath: pairVectLength="<<newAbund->size()<<G4endl;
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
    G4cout<<"G4QNGamma::GetMeanFreePath: Before Loop nIs="<<nIs<<G4endl;
#endif
    if(nIs) for(G4int j=0; j<nIs; j++)      // Calculate CS for eachIsotope of El
    {
      std::pair<G4int,G4double>* curIs=(*cs)[j]; // A pointer, which is used twice
      G4int N=curIs->first;                 // #of Neuterons in the isotope j of El i
      IsN->push_back(N);                    // Remember Min N for the Element
#ifdef debug
      G4cout<<"G4QNGam::GetMeanFrP: Before CS, P="<<Momentum<<",Z="<<Z<<",N="<<N<<G4endl;
#endif
      // @@ Can be a function, depanding on CSm1, CSm2, Momentum, Z, N, pPDG
      G4double CSI=CSmanager->GetCrossSection(true, Momentum, Z, N, pPDG) *
	             capMan->GetRatio(Momentum, Z, N); // CS(j,i) for the isotope
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
  G4cout<<"G4QNGam::GetMeanFrPa: Sigma="<<sigma<<G4endl;
#endif
  if(sigma > 0.) return 1./sigma;                 // Mean path [distance] 
  return DBL_MAX;
}

// Original in any G4QProcess inheriting G4Process
G4bool G4QNGamma::IsApplicable(const G4ParticleDefinition& particle) 
{
  if ( particle == *( G4Neutron::Neutron() ) ) return true; 
#ifdef debug
  G4cout<<"***G4QNGamma::IsApplicable: PDG="<<particle.GetPDGEncoding()<<G4endl;
#endif
  return false;
}

G4VParticleChange* G4QNGamma::PostStepDoIt(const G4Track& track, const G4Step& step)
{
#ifdef debug
  static const G4double mNeut= G4QPDGCode(2112).GetMass();
#endif
  static const G4LorentzVector vacuum4M(0.,0.,0.,0.);
  //-------------------------------------------------------------------------------------
  const G4DynamicParticle* projHadron = track.GetDynamicParticle();
  const G4ParticleDefinition* particle=projHadron->GetDefinition();
#ifdef debug
  G4cout<<"G4QNGamma::PostStepDoIt: Before the GetMeanFreePath is called"<<G4endl;
#endif
  G4ForceCondition cond=NotForced;
  GetMeanFreePath(track, 1., &cond);
#ifdef debug
  G4cout<<"G4QNGamma::PostStepDoIt: After the GetMeanFreePath is called"<<G4endl;
#endif
  G4LorentzVector proj4M=projHadron->Get4Momentum();  // 4-momentum of the projectile (IU?)
  G4double momentum = projHadron->GetTotalMomentum(); // 3-momentum of the Particle
  G4double Momentum=proj4M.rho();
  if(std::fabs(Momentum-momentum)>.001)
                   G4cerr<<"*G4QNGamma::PostStepDoIt: P="<<Momentum<<"#"<<momentum<<G4endl;
#ifdef debug
  G4double mp=proj4M.m(); // @@ must be just the neutron mass
  if(std::fabs(mp-mNeut)>.001)G4cerr<<"*G4QNGamma::PostStDoIt: M="<<mp<<"#"<<mNeut<<G4endl;
  G4cout<<"->G4QNGam::PostStDoIt:*called*,4M="<<proj4M<<",P="<<Momentum<<",m="<<mp<<G4endl;
#endif
  // The same cut function can be used as in MeanFreePath (500)
  if (!IsApplicable(*particle) || Momentum > 500.)  // Check applicability (@@ IU?)
  {
    G4cerr<<"G4QNGamma::PostStepDoIt: Only neutrons with P="<<Momentum<<" < 500"<<G4endl;
    return 0;
  }
  const G4Material* material = track.GetMaterial();      // Get the current material
  G4int Z=0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QNGamma::PostStepDoIt: "<<nE<<" elements in the material."<<G4endl;
#endif
  G4int EPIM=ElProbInMat.size();
#ifdef debug
  G4cout<<"G4QNGam::PostStDoIt: m="<<EPIM<<",n="<<nE<<",T="<<ElProbInMat[EPIM-1]<<G4endl;
#endif
  G4int i=0;
  if(EPIM>1)
  {
    G4double rnd = ElProbInMat[EPIM-1]*G4UniformRand();
    for(i=0; i<nE; ++i)
    {
#ifdef debug
      G4cout<<"G4QNGamma::PostStepDoIt:E["<<i<<"]="<<ElProbInMat[i]<<",r="<<rnd<<G4endl;
#endif
      if (rnd<ElProbInMat[i]) break;
    }
    if(i>=nE) i=nE-1;                        // Top limit for the Element
  }
  G4Element* pElement=(*theElementVector)[i];
  Z=static_cast<G4int>(pElement->GetZ());
#ifdef debug
    G4cout<<"G4QNGamma::PostStepDoIt: i="<<i<<", Z(element)="<<Z<<G4endl;
#endif
  if(Z <= 0)
  {
    G4cerr<<"---Warning---G4QNGamma::PostStepDoIt: Element with Z="<<Z<<G4endl;
    if(Z<0) return 0;
  }
  std::vector<G4double>* SPI = IsoProbInEl[i];// Vector of summedProbabilities for isotopes
  std::vector<G4int>* IsN = ElIsoN[i];     // Vector of "#of neutrons" in the isotope El[i]
  G4int nofIsot=SPI->size();               // #of isotopes in the element i
#ifdef debug
  G4cout<<"G4QNGam::PosStDoIt:n="<<nofIsot<<",T="<<(*SPI)[nofIsot-1]<<G4endl;
#endif
  G4int j=0;
  if(nofIsot>1)
  {
    G4double rndI=(*SPI)[nofIsot-1]*G4UniformRand(); // Randomize the isotop of the Element
    for(j=0; j<nofIsot; ++j)
    {
#ifdef debug
      G4cout<<"G4QNGamma::PostStepDoIt: SP["<<j<<"]="<<(*SPI)[j]<<", r="<<rndI<<G4endl;
#endif
      if(rndI < (*SPI)[j]) break;
    }
    if(j>=nofIsot) j=nofIsot-1;            // Top limit for the isotope
  }
  G4int N =(*IsN)[j];                      // Randomized number of neutrons
#ifdef debug
  G4cout<<"G4QNGamma::PostStepDoIt: Z="<<Z<<", j="<<i<<", N(isotope)="<<N<<G4endl;
#endif
  G4double kinEnergy= projHadron->GetKineticEnergy();
  G4ParticleMomentum dir = projHadron->GetMomentumDirection();
  //if() //DoNothing Action insead of the reaction
  //{
  //  aParticleChange.ProposeEnergy(kinEnergy);
  //  aParticleChange.ProposeLocalEnergyDeposit(0.);
  //  aParticleChange.ProposeMomentumDirection(dir);
  //  aParticleChange.ProposeTrackStatus(fAlive);
  //  return G4VDiscreteProcess::PostStepDoIt(track,step);
  //}
  if(N<0)
  {
    G4cerr<<"-Warning-G4QNGamma::PostStepDoIt: Isotope with Z="<<Z<<", 0>N="<<N<<G4endl;
    return 0;
  }
  nOfNeutrons=N;                           // Remember it for the energy-momentum check
#ifdef debug
  G4cout<<"G4QNGamma::PostStepDoIt: N="<<N<<" for element with Z="<<Z<<G4endl;
#endif
  aParticleChange.Initialize(track);
  G4double weight = track.GetWeight();
#ifdef debug
  G4cout<<"G4QNGamma::PostStepDoIt: weight="<<weight<<G4endl;
#endif
  G4double localtime = track.GetGlobalTime();
#ifdef debug
  G4cout<<"G4QNGamma::PostStepDoIt: localtime="<<localtime<<G4endl;
#endif
  G4ThreeVector position = track.GetPosition();
  G4TouchableHandle trTouchable = track.GetTouchableHandle();
#ifdef debug
  G4cout<<"G4QNGamma::PostStepDoIt: position="<<position<<G4endl;
#endif
  G4int targPDG = 90000000 + Z*1000 + N;                  // PDG Code of the target nucleus
  G4QPDGCode targQPDG(targPDG);
  G4double tM = targQPDG.GetMass();                       // Target mass
#ifdef debug
  G4cout<<"G4QNGamma::PostStepDoIt: n + targPDG="<<targPDG<<G4endl;
#endif
  // @@ All above is universal for all processes except for the additional condition (500)
  G4LorentzVector tot4M=G4LorentzVector(0.,0.,0.,tM)+proj4M;
  G4double totM2=tot4M.m2();
  G4int tZ=Z;
  G4int tN=N+1;
  G4int resPDG = targPDG + 1;                             // Final ++N nucleus PDG
  G4double rM=G4QPDGCode(resPDG).GetMass();               // Mass of the final nucleus
  G4LorentzVector r4M=G4LorentzVector(0.,0.,0.,rM);       // 4mom of the final nucleus
  G4LorentzVector g4M=G4LorentzVector(0.,0.,0.,0.);       // 4mom of the gamma
#ifdef debug
  G4cout<<"G4QNGamma::PostStepDoIt: tM="<<tM << ", rM="<<rM << ", Q="<<tM+mNeut-rM<<G4endl;
#endif
  if(!G4QHadron(tot4M).DecayIn2(r4M, g4M)) // The compoun decay din't succeed
  {
    //G4cerr<<"G4QNGamma::PostStDoIt: tM="<<std::sqrt(totM2)<<" < rM="<<rM<<G4endl;
    //G4Exception("G4QNGamma::PostStepDoIt()", "HAD_CHPS_0001",
    //            FatalException, "Hadronize quasmon: Can't decay TotNuc->ResNuc+gam");
    G4cerr<<"-Warning-G4QNGamma::PostStDoIt: tM="<<std::sqrt(totM2)<<" < rM="<<rM<<G4endl;
    aParticleChange.ProposeEnergy(kinEnergy);
    aParticleChange.ProposeLocalEnergyDeposit(0.);
    aParticleChange.ProposeMomentumDirection(dir);
    aParticleChange.ProposeTrackStatus(fAlive);
    return G4VDiscreteProcess::PostStepDoIt(track,step);
  }
#ifdef debug
  G4cout<<"G4QNGam::PStDoIt: RA="<<r4M.rho()<<r4M<<", Gamma="<<g4M.rho()<<g4M<<G4endl;
#endif
  EnMomConservation = tot4M - r4M - g4M;           // EM conservation check 4mom
  aParticleChange.ProposeEnergy(0.);               // A standard procedure of killing proj.
  aParticleChange.ProposeTrackStatus(fStopAndKill);// projectile neutron is killed
  aParticleChange.SetNumberOfSecondaries(2);       // Fix a#of secondaries
  // Fill the gamma
  G4ParticleDefinition* theDefinition = G4Gamma::Gamma();
  G4DynamicParticle* theGam = new G4DynamicParticle(theDefinition, g4M);
  G4Track* capGamma = new G4Track(theGam, localtime, position );
  capGamma->SetWeight(weight);
  capGamma->SetTouchableHandle(trTouchable);
  aParticleChange.AddSecondary(capGamma);
  // ----------------------------------------------------
  // Fill the final nucleus
  G4int tA=tZ+tN;
  if     (resPDG==90000001) theDefinition = G4Neutron::Neutron();
  else if(resPDG==90001000) theDefinition = G4Proton::Proton();
  else theDefinition = G4ParticleTable::GetParticleTable()->FindIon(tZ, tA, 0, tZ);
  G4DynamicParticle* theReN = new G4DynamicParticle(theDefinition, r4M);
  G4Track* scatReN = new G4Track(theReN, localtime, position );
  scatReN->SetWeight(weight);
  scatReN->SetTouchableHandle(trTouchable);
  aParticleChange.AddSecondary(scatReN);

  return G4VDiscreteProcess::PostStepDoIt(track, step);
}
