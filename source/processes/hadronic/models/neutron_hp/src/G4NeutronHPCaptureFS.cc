// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "G4NeutronHPCaptureFS.hh"
#include "G4Gamma.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
#include "G4PhotonEvaporation.hh"
#include "G4Fragment.hh"
#include "G4ParticleTable.hh" 
#include "G4NeutronHPDataUsed.hh"

  G4ParticleChange * G4NeutronHPCaptureFS::ApplyYourself(const G4Track & theTrack)
  {
    G4int i;
    theResult.Initialize(theTrack);   
        
// prepare neutron
    G4double eKinetic = theTrack.GetKineticEnergy();
    const G4DynamicParticle *incidentParticle = theTrack.GetDynamicParticle();
    G4ReactionProduct theNeutron( incidentParticle->GetDefinition() );
    theNeutron.SetMomentum( incidentParticle->GetMomentum() );
    theNeutron.SetKineticEnergy( eKinetic );

// prepare target
    G4ReactionProduct theTarget; 
    G4Nucleus aNucleus;
    G4double eps = 0.0001;
    if(targetMass<500*MeV)
      targetMass = ( G4NucleiPropertiesTable::GetAtomicMass(theBaseZ+eps, theBaseA+eps)-
                       theBaseZ* G4Electron::Electron()->GetPDGMass() ) /
                     G4Neutron::Neutron()->GetPDGMass();
    theTarget = aNucleus.GetThermalNucleus(targetMass);

// go to nucleus rest system
    theNeutron.Lorentz(theNeutron, -1*theTarget);
    eKinetic = theNeutron.GetKineticEnergy();

// dice the photons

    G4ReactionProductVector * thePhotons = NULL;    
    if (HasFSData()) 
    { 
      thePhotons = theFinalStatePhotons.GetPhotons(eKinetic);
    }
    else
    {
      G4ThreeVector aCMSMomentum = theNeutron.GetMomentum()+theTarget.GetMomentum();
      G4LorentzVector p4(aCMSMomentum, theTarget.GetTotalEnergy() + theNeutron.GetTotalEnergy());
      G4Fragment nucleus(theBaseA+1, theBaseZ ,p4);
      G4PhotonEvaporation photonEvaporation;
      G4FragmentVector* products = photonEvaporation.BreakItUp(nucleus);
      G4int i;
      thePhotons = new G4ReactionProductVector;
      for(i=0; i<products->entries(); i++)
      {
        G4ReactionProduct * theOne = new G4ReactionProduct;
        theOne->SetDefinition( G4Gamma::Gamma() );
        G4ParticleTable* theTable = G4ParticleTable::GetParticleTable();
        if(products->at(i)->GetMomentum().mag() > 10*MeV) 
                 theOne->SetDefinition( 
                 theTable->FindIon(theBaseZ, theBaseA+1, 0, theBaseZ) );
        theOne->SetMomentum( products->at(i)->GetMomentum().vect() ) ;
        theOne->SetTotalEnergy( products->at(i)->GetMomentum().t() );
        thePhotons->insert(theOne);
        delete products->at(i);
      } 
      delete products;
    }

// add them to the final state

    G4int nPhotons = 0;
    if(thePhotons!=NULL) nPhotons=thePhotons->length();
    G4int nParticles = nPhotons;
    if(1==nPhotons) nParticles = 2;
    theResult.SetNumberOfSecondaries(nParticles);

    // back to lab system
    for(i=0; i<nPhotons; i++)
    {
      thePhotons->at(i)->Lorentz(*(thePhotons->at(i)), theTarget);
    }
    
    // Recoil, if only one gamma
    if (1==nPhotons)
    {
       G4DynamicParticle * theOne = new G4DynamicParticle;
       G4ParticleDefinition * aRecoil = G4ParticleTable::GetParticleTable()
                                        ->FindIon(theBaseZ, theBaseA+1, 0, theBaseZ);
       theOne->SetDefinition(aRecoil);
       // Now energy; assume that momentum gets taken by 'surrounding'.
       G4ThreeVector aMomentum =  theTrack.GetMomentum()
                                 +theTarget.GetMomentum()
				 -thePhotons->at(0)->GetMomentum();
       G4ThreeVector theMomUnit = aMomentum.unit();
       G4double aKinEnergy =  theTrack.GetKineticEnergy()
                             +theTarget.GetKineticEnergy()
			     -thePhotons->at(0)->GetKineticEnergy();
       G4double theResMass = aRecoil->GetPDGMass();
       G4double theResE = aRecoil->GetPDGMass()+aKinEnergy;
       G4double theAbsMom = sqrt(theResE*theResE - theResMass*theResMass);
       G4ThreeVector theMomentum = theAbsMom*theMomUnit;
       theOne->SetMomentum(theMomentum);
       theResult.AddSecondary(theOne);
    }

    // Now fill in the gammas.
    for(i=0; i<nPhotons; i++)
    {
      // back to lab system
      G4DynamicParticle * theOne = new G4DynamicParticle;
      theOne->SetDefinition(thePhotons->at(i)->GetDefinition());
      theOne->SetMomentum(thePhotons->at(i)->GetMomentum());
      theResult.AddSecondary(theOne);
      delete thePhotons->at(i);
    }
    delete thePhotons; 
// clean up the primary neutron
    theResult.SetStatusChange(fStopAndKill);
    return &theResult;
  }

  void G4NeutronHPCaptureFS::Init (G4double A, G4double Z, G4String & dirName, G4String & aFSType)
  {
    G4String tString = "/FS/";
    G4bool dbool;
    G4NeutronHPDataUsed aFile = theNames.GetName(A, Z, dirName, tString, dbool);
    G4String filename = aFile.GetName();
    theBaseA = A;
    theBaseZ = G4int(Z+.5);
    if(!dbool)
    {
      hasAnyData = false;
      hasFSData = false; 
      hasXsec = false;
      return;
    }
    G4std::ifstream theData(filename, G4std::ios::in);
    
    hasFSData = theFinalStatePhotons.InitMean(theData); 
    if(hasFSData)
    {
      targetMass = theFinalStatePhotons.GetTargetMass();
      theFinalStatePhotons.InitAngular(theData); 
      theFinalStatePhotons.InitEnergies(theData); 
    }
  }
