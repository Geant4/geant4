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
// $Id: G4QInelastic.cc,v 1.1 2004-03-05 13:28:34 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: M. Kossov (Mikhail.Kossov@cern.ch)
//
// History:
// ----------------------------------------------------------------------------------------
// 20 Dec 2003   M. Kossov  "Hadronic package independent" photo-nuclear process is created
// ----------------------------------------------------------------------------------------

#include "G4QInelastic.hh"

G4QInelastic::G4QInelastic(const G4String& processName)
  : G4VDiscreteProcess(processName), lowEnergyLimit(2.2*MeV) // Preinitialization of CHIPS
{
  if(verboseLevel>0)
    G4cout<<GetProcessName()<<" is created Thresh="<<lowEnergyLimit/MeV<<" MeV"<<G4endl;
  G4QCHIPSWorld::Get()->GetParticles(234);           // Create CHIPS World of 234 particles
  G4QNucleus::SetParameters(0.,0.,1.,1.);            // Nuclear clusterization parameters
  G4Quasmon::SetParameters(180.,.09,.3);             // Temperature, s-antis, eta suppress
  G4QEnvironment::SetParameters(.5);                 // SolAngle (pbar-A secondary capture)
  thePhotonCS = new G4QPhotoNuclearCrossSection();   // Define the CrossSection Meneger
}

G4QInelastic::~G4QInelastic()
{
  delete thePhotonCS;
}

// Here the cross-section of lepto-nuclear interaction is used ============================
G4double G4QInelastic::GetMeanFreePath(const G4Track& aTrack,
                                                  G4double, G4ForceCondition*) //Not used
{
  const G4DynamicParticle* projectile = aTrack.GetDynamicParticle();
  const G4ParticleDefinition* particle= projectile->GetDefinition();
  if     (particle==G4Gamma::Gamma()) theInelCS=thePhotonCS;
  if(!IsApplicable(*particle))
    G4cout<<"-Wor-G4QInel::GetMeanFreePath calledFor "<<particle->GetPDGEncoding()<<G4endl;
  // Calculate the mean Cross Section for the set of Elements(*Isotopes) in the Material
  G4double kineticEnergy = projectile->GetKineticEnergy(); // Kinetic energy of the lepton
  const G4Material* material = aTrack.GetMaterial();         // Get the current material
  const G4double* NOfNucPerVolume = material->GetVecNbOfAtomsPerVolume();
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QInelastic::GetMeanFreePath:"<<nE<<" Elems in theMaterial"<<G4endl;
#endif
  G4QIsotope* Isotopes = G4QIsotope::Get(); // Pointer to the G4QIsotopes singelton
  G4double sigma=0.;
  for(G4int i=0; i<nE; ++i)
  {
    G4int Z = static_cast<G4int>((*theElementVector)[i]->GetZ()); // Z of the Element
    std::vector<std::pair<G4int,G4double>*>* cs= Isotopes->GetCSVector(Z); // Pointer to CS
    G4int nIs=cs->size();                         // A#Of Isotopes in the Element
    if(nIs) for(G4int j=0; j<nIs; j++)            // Calculate CS for eachIsotope of El
    {
      G4int N=cs->at(j)->first;                   // #ofNeuterons in the isotope
      cs->at(j)->second = theInelCS->GetCrossSection(kineticEnergy,Z,N);// CS calculation
    } // End of temporary initialization of the cross sections in the G4QIsotope singeltone
    sigma+=Isotopes->GetMeanCrossSection(Z)*NOfNucPerVolume[i]; // SUM(MeanCS*NOFNperV)
  } // End of LOOP over Elements

  // Check that cross section is not zero and return the mean free path
  if(sigma > 0.) return 1./sigma;                 // Mean path [distance] 
  return DBL_MAX;
}

G4VParticleChange* G4QInelastic::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);

  // Dynamic particle quantities
  const G4DynamicParticle* projectile = aTrack.GetDynamicParticle();
  const G4ParticleDefinition* particle= projectile->GetDefinition();
  G4double mass=0., m2=0.; // For the particle with m#0 mass & m2 must be redefined in "if"
  if     (particle==G4Gamma::Gamma()) theInelCS=thePhotonCS; 
  G4double kineticEnergy = projectile->GetKineticEnergy();
  G4ParticleMomentum dir = projectile->GetMomentumDirection();
  if(kineticEnergy <= lowEnergyLimit) // The reaction is below the hadronic threshold
  {
    //Do Nothing Action insead of the reaction
    aParticleChange.SetEnergyChange(kineticEnergy);
    aParticleChange.SetLocalEnergyDeposit(0.);
    aParticleChange.SetMomentumChange(dir) ;
    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }

  // Select randomly one element in the current material and isotope in the element
  const G4Material* material = aTrack.GetMaterial();      // Get the current material
  G4int Z=0;
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4int i=0;
  G4double sum=0.;
  G4int nE=material->GetNumberOfElements();
#ifdef debug
  G4cout<<"G4QInelastic::PostStepDoIt:"<<nE<<" Elem's in the Material."<<G4endl;
#endif
  G4QIsotope* Isotopes = G4QIsotope::Get();         // Pointer to the G4QIsotopes singelton
  std::vector<G4double> sumfra;
  for(i=0; i<nE; ++i)
  {
	G4double frac=material->GetFractionVector()[i];   // fraction of the Element
    G4int Z = static_cast<G4int>((*theElementVector)[i]->GetZ()); // Z of the Element
    std::vector<std::pair<G4int,G4double>*>* cs= Isotopes->GetCSVector(Z); // Pointer to CS
    G4int nIs=cs->size();                             // A#Of Isotopes in the Element
    if(nIs) for(G4int j=0; j<nIs; j++)                // Calculate CS for eachIsotope of El
    {
      G4int N=cs->at(j)->first;                       // #ofNeuterons in the isotope
      cs->at(j)->second = theInelCS->GetCrossSection(kineticEnergy,Z,N);// CS calculation
    } // End of temporary initialization of the cross sections in the G4QIsotope singeltone
	frac*=Isotopes->GetMeanCrossSection(Z);           // (MeanCS for theElement)*(fraction)
    //G4int cZ=static_cast<G4int>((*theElementVector)[i]->GetZ());
    sum+=frac;
    sumfra.push_back(sum);                            // remember the summation steps
  }
  G4double rnd = sum*G4UniformRand();
  for(i=0; i<nE; ++i)
  {
    G4int cZ=static_cast<G4int>((*theElementVector)[i]->GetZ());
    sum=sumfra[i];
    if (rnd<sum)  
    { 
	  Z = cZ;
      break;
	}
  }
  if(Z<=0)
  {
    G4cerr<<"--Worning--G4QInelastic::PostStepDoIt:Element with Z="<<Z<< G4endl;
    if(Z<0) return 0;
  }
  // Z is known (CS's are initialized), now N value must be found for the selected Element
  G4int N = Isotopes->GetCSNeutrons(Z);                // Randomize N of the Isotope withCS
  if(Z+N>20) G4QNucleus::SetParameters(.18,.06,6.,1.); // HeavyNuclei NuclearClusterization
  else       G4QNucleus::SetParameters(0.0,0.0,1.,1.); // LightNuclei NuclearClusterization
#ifdef debug
  G4cout<<"G4QInelastic::PostStepDoIt:N="<<N<<" for Element with Z="<<Z<<G4endl;
#endif
  if(N<0)
  {
    G4cerr<<"---Worning--G4QInelastic::PostStepDoIt:Element with N="<<N<<G4endl;
    return 0;
  }
  G4double xSec=theInelCS->GetCrossSection(kineticEnergy,Z,N);// Recalculate CrossSection
  // @@ check a possibility to separate p, n, or alpha (!)
  if(xSec <= 0.) // The cross-section iz 0 -> Do Nothing
  {
    //Do Nothing Action insead of the reaction
    aParticleChange.SetEnergyChange(kineticEnergy);
    aParticleChange.SetLocalEnergyDeposit(0.);
    aParticleChange.SetMomentumChange(dir) ;
    return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
  }
  else // Reaction is possible -> Kill the absorbed particle
  {
    aParticleChange.SetEnergyChange(0.) ;
    aParticleChange.SetStatusChange(fStopAndKill);
  }
  // Scatter the lepton
  G4double energy=kineticEnergy;
  G4double momentum=kineticEnergy;
  if(mass){energy+=mass;momentum=sqrt(energy*energy-m2);}
  G4ThreeVector threeMom=momentum*dir;
  // --- CHIPS interaction of the virtual photon
  G4int targPDG=90000000+Z*1000+N;         // PDG Code of the target nucleus
  G4QHadronVector* output=0;               // Prototype of the output G4QHadronVector
  G4QHadron* pH = new G4QHadron(22,G4LorentzVector(threeMom, energy)); // -------------+
  G4QHadronVector projHV;                  //                                          |
  projHV.push_back(pH);                    // DESTROYED over 1 line -----------------+ |
  G4QEnvironment* pan= new G4QEnvironment(projHV,targPDG);// ---> DELETED ---------+ | |
  std::for_each(projHV.begin(), projHV.end(), DeleteQHadron()); // ----------------+-+-+
  projHV.clear(); // --------------------------------------------------------------+-+
#ifdef debug
  G4cout<<"G4QInelastic::PostStepDoIt:W="<<W<<G4endl; //                 |
#endif
  try                                                           //                 |
  {                                                             //                 |
    output = pan->Fragment();// DESTROYED in the end of the LOOP work space        |
  }                                                             //                 |
  catch (G4QException& error)//                                                    |
  {                                                             //                 |
    G4cerr<<"G4QInelastic::PostStepDoIt:Exception from CHIPS"<<G4endl;//           |
    G4Exception("G4QInelastic::PostStepDoIt:","27",FatalException,"CHIPSError");// /
  }                                                             //                 |
  delete pan;                              // Delete the Nuclear Environment ------+
  // Fill the generated hadrons to the G4ParticleChange's "Secondaries"
  G4int tNH = output->size();                // A#of hadrons in the output
#ifdef debug
  G4cout<<"G4QInelastic::PostStepDoIt:"<<tNH<<" hadronss are generated"<<G4endl;
#endif
  // Deal with ParticleChange final state interface to GEANT4 output of the process
  aParticleChange.Initialize(aTrack);          // @@ (?)
  aParticleChange.SetNumberOfSecondaries(tNH); 
  for(i=0; i<tNH; i++)
  {
    // Note that one still has to take care of Hypernuclei (with Lambda or Sigma inside)
    // Hypernucleus mass calculation and ion-table interface upgrade => work for Hisaya @@
    // The decau process for hypernuclei must be developed in GEANT4 (change CHIPS body)
    G4double localtime = aTrack.GetGlobalTime();
    G4ThreeVector   position = aTrack.GetPosition();
    G4QHadron* hadr=output->operator[](i);   // Pointer to the output hadron    
    if(hadr->GetNFragments())                // Intermediate hadron
    {
#ifdef debug
	  G4cout<<"G4QInelastic::PostStepDoIt:Intermediate particle i="<<i<<G4endl;
#endif
      delete hadr;
      continue;
    }
    G4DynamicParticle* theSec = new G4DynamicParticle;  
    G4int PDGCode = hadr->GetPDGCode();
#ifdef pdebug
    G4cout<<"G4QInelastic::PostStepDoIt:#"<<i<<",PDG="<<PDGCode<<G4endl;
#endif
    G4ParticleDefinition * theDefinition;
    if     (PDGCode==90000001) theDefinition = G4Neutron::Neutron();
    else if(PDGCode==91000000) theDefinition = G4Lambda::Lambda();
    else if(PDGCode==91000999) theDefinition = G4SigmaPlus::SigmaPlus();
    else if(PDGCode==90999001) theDefinition = G4SigmaMinus::SigmaMinus();
    else if(PDGCode==91999000) theDefinition = G4XiMinus::XiMinus();
    else if(PDGCode==91999999) theDefinition = G4XiZero::XiZero();
    else if(PDGCode==92998999) theDefinition = G4OmegaMinus::OmegaMinus();
	else if(PDGCode >80000000) // Defines hypernuclei as normal nuclei (N=N+S Correction!)
    {
      G4int aZ = hadr->GetCharge();
      G4int aA = hadr->GetBaryonNumber();
      theDefinition = G4ParticleTable::GetParticleTable()->FindIon(aZ,aA,0,aZ);
    }
    else theDefinition = G4ParticleTable::GetParticleTable()->FindParticle(PDGCode);
    if(!theDefinition)
    {
      G4cout<<"--Worning--G4QInelastic::PostStepDoIt:dropPDG="<<PDGCode<<G4endl;
      delete hadr;
      continue;
    }
    theSec->SetDefinition(theDefinition);
    G4LorentzVector h4M=hadr->Get4Momentum();
#ifdef debug
    G4cout<<"G4QInelastic::PostStepDoIt:"<<i<<",PDG/4M="<<PDGCode<<<h4M<<G4endl;
#endif
    theSec->SetMomentum(h4M.vect());
    delete hadr;
    G4Track* aNewTrack = new G4Track(theSec, localtime, position );
    aParticleChange.AddSecondary( aNewTrack );
  }
  delete output;
  // Return the standard G4ParticleChange with standard actions (reset etc.)
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// Just check that this is a lepton for fhich the cross-sections are defined
G4bool G4QInelastic::IsApplicable(const G4ParticleDefinition& particle)
//---- ==========================-------------------------------------- 
{
  if(&particle==G4Gamma::Gamma()) return true;
  return false;
}
