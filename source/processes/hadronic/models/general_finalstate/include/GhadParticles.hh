#ifndef GhadParticles_h
#define GhadParticles_h

#include <vector>
#include "GhadTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4LorentzVector.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4KineticTrack.hh"
#include "GhadParticles.hh"

class GhadParticles : public std::vector<GhadTrack>
{
  public:
    virtual ~GhadParticles() {}
    void Initialize (const G4Track & it, G4ThreeVector aPos)
    {
      clear();
      G4DynamicParticle * it1 = const_cast<G4DynamicParticle *>(it.GetDynamicParticle());
      G4LorentzVector theM = it1->Get4Momentum();
      G4KineticTrack * aPro = new G4KineticTrack(it1->GetDefinition(), 0, aPos, theM );
      push_back(GhadTrack(aPro));
    }
    
    void Initialize (G4KineticTrackVector * it)
    {
      G4int i=0;
      clear();
      for(i=0; i<it->size(); i++)
      {
        push_back(GhadTrack( it->operator[](i) ));
      }
    }
    
    G4int GetNCharged()
    {
      G4int result = 0;
      GhadParticles::iterator i=0;
      for(i=begin(); i!= end(); i++)
      {
        if(i->GetDefinition()->GetPDGCharge() != 0)
	{
	  result ++;
	}
      }
      return result;
    }
    
    G4LorentzVector GetTotal4Momentum()
    {
      G4LorentzVector result(0,0,0,0);
      GhadParticles::iterator i=0;
      for(i=begin(); i!= end(); i++)
      {
	result +=i->GetMom();
      }
      return result;
    }

    G4double GetTotalKinetic()
    {
      G4double result=0;
      GhadParticles::iterator i=0;
      for(i=begin(); i!= end(); i++)
      {
	result +=i->GetMom().t() - i->GetDefinition()->GetPDGMass();
      }
      return result;
    }
};

#endif
