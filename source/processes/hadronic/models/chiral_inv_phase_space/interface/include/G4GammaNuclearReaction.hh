#ifndef G4GammaNuclearReaction_h
#define G4GammaNuclearReaction_h

#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4ParticleTable.hh"
#include "G4Quasmon.hh"
#include "G4QHadronVector.hh"
#include "G4ParticleChange.hh"
#include "G4LorentzVector.hh"
#include "G4DynamicParticle.hh"
#include "G4IonTable.hh"
#include "G4Neutron.hh"

class G4GammaNuclearReaction : public G4HadronicInteraction
{
  public: 
    G4VParticleChange * ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus);

  private:
    G4ParticleChange theResult;
};

inline
G4VParticleChange * G4GammaNuclearReaction::
ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus)
{
  //projectile properties needed in constructor of quasmon
  G4LorentzVector proj4Mom;
  proj4Mom = aTrack.GetDynamicParticle()->Get4Momentum();
  G4int projectilePDGCode = aTrack.GetDynamicParticle()
                                  ->GetDefinition()
				  ->GetPDGEncoding();
  
  
  //target properties needed in constructor of quasmon
  G4int targetZ = G4int(aTargetNucleus.GetZ()+0.5);
  G4int targetA = G4int(aTargetNucleus.GetN()+0.5);
  G4int targetPDGCode = 90000000 + 1000*targetZ + (targetA-targetZ);
  G4double targetMass = G4ParticleTable::GetParticleTable()->GetIonTable()
                                                           ->GetIonMass(targetZ, targetA);
  G4LorentzVector targ4Mom(0.,0.,0.,targetMass);  
  
  G4int nop = 223; // ??????
  
  // construct and fragment the quasmon
  G4Quasmon* pan= new G4Quasmon(projectilePDGCode, targetPDGCode, 1./MeV*proj4Mom, 1./MeV*targ4Mom, nop);
  G4QHadronVector output = pan->HadronizeQuasmon();
  
  // Fill the particle change.
  theResult.Initialize(aTrack);
  theResult.SetStatusChange(fStopAndKill);
  theResult.SetNumberOfSecondaries(output.length());
  G4DynamicParticle * theSec;
  G4cout << "NEW EVENT"<<endl;
  for(G4int particle = 0; particle < output.length(); particle++)
  {
    if(output[particle]->GetNFragments() != 0) 
    {
      delete output[particle];
      continue;
    }
    theSec = new G4DynamicParticle;  
    G4int pdgCode = output[particle]->GetPDGCode();
    G4ParticleDefinition * theDefinition;
    if(pdgCode>90000000) 
    {
      G4int aZ = (pdgCode-90000000)/1000;
      G4int anN = pdgCode-90000000-1000*aZ;
      theDefinition = G4ParticleTable::GetParticleTable()->FindIon(aZ,anN+aZ,0,aZ);
      if(aZ == 0 && anN == 1) theDefinition = G4Neutron::Neutron();
    }
    else
    {
      theDefinition = G4ParticleTable::GetParticleTable()->FindParticle(output[particle]->GetPDGCode());
    }
    G4cout << "Particle code produced = "<< pdgCode <<endl;
    theSec->SetDefinition(theDefinition);
    theSec->SetMomentum(output[particle]->Get4Momentum().vect());
    theResult.AddSecondary(theSec); 
    delete output[particle];
  }
  return & theResult;
}

#endif
