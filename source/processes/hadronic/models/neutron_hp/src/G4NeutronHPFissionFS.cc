// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPFissionFS.hh"
#include "G4Nucleus.hh"
#include "G4DynamicParticleVector.hh"
#include "G4NeutronHPFissionERelease.hh"

  G4NeutronHPFissionFS::G4NeutronHPFissionFS(){ hasXsec = false; }
  G4NeutronHPFissionFS::~G4NeutronHPFissionFS(){}

  void G4NeutronHPFissionFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
  {
     theFS.Init(A, Z, dirName, aFSType);
     theFC.Init(A, Z, dirName, aFSType);
     theSC.Init(A, Z, dirName, aFSType);
     theTC.Init(A, Z, dirName, aFSType);
     theLC.Init(A, Z, dirName, aFSType);
  }
  G4ParticleChange * G4NeutronHPFissionFS::ApplyYourself(const G4Track & theTrack)
  {  
    theResult.Initialize(theTrack); 

// prepare neutron
    G4double eKinetic = theTrack.GetKineticEnergy();
    const G4DynamicParticle *incidentParticle = theTrack.GetDynamicParticle();
    G4ReactionProduct theNeutron( incidentParticle->GetDefinition() );
    theNeutron.SetMomentum( incidentParticle->GetMomentum() );
    theNeutron.SetKineticEnergy( eKinetic );

// prepare target
    G4Nucleus aNucleus;
    G4ReactionProduct theTarget; 
    G4double targetMass = theFS.GetMass();
    theTarget = aNucleus.GetThermalNucleus(targetMass);
    
// set neutron and target in the FS classes 
   theFS.SetNeutron(theNeutron);
   theFS.SetTarget(theTarget);
   theFC.SetNeutron(theNeutron);
   theFC.SetTarget(theTarget);
   theSC.SetNeutron(theNeutron);
   theSC.SetTarget(theTarget);
   theTC.SetNeutron(theNeutron);
   theTC.SetTarget(theTarget);
   theLC.SetNeutron(theNeutron);
   theLC.SetTarget(theTarget);
    
// boost to target rest system and decide on channel.
    theNeutron.Lorentz(theNeutron, -1*theTarget);

// dice the photons

    G4DynamicParticleVector * thePhotons;    
    thePhotons = theFS.GetPhotons();

// select the FS in charge

    eKinetic = theNeutron.GetKineticEnergy();    
    G4double xSec[4];
    xSec[0] = theFC.GetXsec(eKinetic);
    xSec[1] = xSec[0]+theSC.GetXsec(eKinetic);
    xSec[2] = xSec[1]+theTC.GetXsec(eKinetic);
    xSec[3] = xSec[2]+theLC.GetXsec(eKinetic);
    G4int i, it;
    G4double random = G4UniformRand();
    for(i=0; i<4; i++)
    {
      it =i;
      if(random<xSec[i]/xSec[3]) break;
    }
    if(xSec[3]==0) it=-1;
    
// dice neutron multiplicities, energies and momenta in Lab. @@
// no energy conservation on an event-to-event basis. we rely on the data to be ok. @@
// also for mean, we rely on the consistancy of the data. @@

    G4int Prompt=0, delayed=0, all=0;
    G4DynamicParticleVector * theNeutrons = NULL;
    switch(it) // check logic, and ask, if partials can be assumed to correspond to individual particles @@@
    {
      case 0:
        theFS.SampleNeutronMult(all, Prompt, delayed, eKinetic, 0);
        if(Prompt==0&&delayed==0) Prompt=all;
        theNeutrons = theFC.ApplyYourself(Prompt); // delayed always in FS 
	// take 'U' into account explicitely (see 5.4) in the sampling of energy @@@@
        break;
      case 1:
        theFS.SampleNeutronMult(all, Prompt, delayed, eKinetic, 1);
        if(Prompt==0&&delayed==0) Prompt=all;
        theNeutrons = theSC.ApplyYourself(Prompt); // delayed always in FS, off done in FSFissionFS
        break;
      case 2:
        theFS.SampleNeutronMult(all, Prompt, delayed, eKinetic, 2);
        if(Prompt==0&&delayed==0) Prompt=all;
        theNeutrons = theTC.ApplyYourself(Prompt); // delayed always in FS
        break;
      case 3:
        theFS.SampleNeutronMult(all, Prompt, delayed, eKinetic, 3);
        if(Prompt==0&&delayed==0) Prompt=all;
        theNeutrons = theLC.ApplyYourself(Prompt); // delayed always in FS
        break;
      default:
        break;
    }
// dice delayed neutrons and photons, and fallback 
// for Prompt in case channel had no FS data; add all paricles to FS.

    G4double * theDecayConstants;
    
    if(theNeutrons != NULL)
    {
      theDecayConstants = new G4double[delayed];
      G4int nPhotons = 0;
      if(thePhotons!=NULL) nPhotons = thePhotons->length();
      theResult.SetNumberOfSecondaries(nPhotons+Prompt+delayed);
      for(i=0; i<theNeutrons->length(); i++)
      {
        theResult.AddSecondary(theNeutrons->at(i));
      }
      delete theNeutrons;  
      
      G4DynamicParticleVector * theDelayed = NULL;
      theDelayed = theFS.ApplyYourself(0, delayed, theDecayConstants);
      for(i=0; i<theDelayed->length(); i++)
      {
        G4double time = -log(G4UniformRand())/theDecayConstants[i];
        time += theResult.GetTimeChange();
        theResult.AddSecondary(theDelayed->at(i), time);
      }
      delete theDelayed;                  
    }
    else
    {
//    cout << " all = "<<all<<endl;
      theFS.SampleNeutronMult(all, Prompt, delayed, eKinetic, 0);
      theDecayConstants = new G4double[delayed];
      if(Prompt==0&&delayed==0) Prompt=all;
      theNeutrons = theFS.ApplyYourself(Prompt, delayed, theDecayConstants);
      G4int nPhotons = 0;
      if(thePhotons!=NULL) nPhotons = thePhotons->length();
      theResult.SetNumberOfSecondaries(nPhotons+Prompt+delayed);
      for(i=0; i<Prompt; i++)
      {
        theResult.AddSecondary(theNeutrons->at(i));
      }
      for(i=Prompt; i<Prompt+delayed; i++)
      {
        G4double time = -log(G4UniformRand())/theDecayConstants[i-Prompt];
        time += theResult.GetTimeChange();        
        theResult.AddSecondary(theNeutrons->at(i), time);
      }
      delete theNeutrons;   
    }
    delete [] theDecayConstants;
//    cout << "all delayed "<<delayed<<endl; 
    G4int nPhotons = 0;
    if(thePhotons!=NULL)
    {
      nPhotons = thePhotons->length();
      for(i=0; i<thePhotons->length(); i++)
      {
        theResult.AddSecondary(thePhotons->at(i));
      }
      delete thePhotons; 
    }

// do some rotating, if that helps to conserve momentum @@@@
         
// finally deal with local energy depositions.
//    G4cout <<"Number of secondaries = "<<theResult.GetNumberOfSecondaries()<< endl;
//    G4cout <<"Number of Prompt = "<<Prompt<<endl;
//    G4cout <<"Number of delayed = "<<delayed<<endl;
//    G4cout <<"Number of photons = "<<nPhotons<<endl;
    G4NeutronHPFissionERelease * theERelease;
    theERelease = theFS.GetEnergyRelease();
    G4double eDepByFragments = theERelease->GetFragmentKinetic();
    theResult.SetLocalEnergyDeposit(eDepByFragments);
//    cout << "local energy deposit" << eDepByFragments<<endl;
// clean up the primary neutron
    theResult.SetStatusChange(fStopAndKill);
    return &theResult;
  }
