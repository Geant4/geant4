#ifndef GhadNucleus_h
#define GhadNucleus_h

#include "G4V3DNucleus.hh"
#include "GhadParticles.hh"
#include "GhadFS.hh"
#include "G4Nucleus.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4Scatterer.hh"
#include "GhadNuclearGeometry.hh"


class GhadNucleus 
{
  public:
    virtual ~GhadNucleus() {}
    
    void Initialize(G4Nucleus & aNuc) 
    {
      theTarget = new G4Fancy3DNucleus();
      theTarget->Init(aNuc.GetN(), aNuc.GetZ());
      theShellModel.Init(theTarget);
      nHits = 0;
      initialized = -1111;
    }
    
    void Initialize(G4V3DNucleus * aNuc) 
    {
      theTarget = aNuc;
      theShellModel.Init(theTarget);
      const std::vector<G4Nucleon *> & nucs =  aNuc->GetNucleons();
      nHits = 0;
      for(G4int i=0; i<nucs.size(); i++) if(nucs[i]->AreYouHit()) nHits++;
      initialized = -1112;
    }
    
    GhadFS * YouAreTraversedBy(GhadParticles & aSetOfPro)
    {
      vector<G4Nucleon *> theNuc = theTarget->GetNucleons();
      GhadProjectileSearch aSearch(theShellModel, theNuc);
      
      for_each(aSetOfPro.begin(), aSetOfPro.end(), aSearch);
      
      GhadParticles::iterator thePro = aSearch.GetProjectile();
      G4Nucleon * theTarg = aSearch.GetTarget();
      GhadAction limit(aSearch.GetLimitation());
      G4double aStep = aSearch.GetStep();
      
      GhadFS * result = 0;
      
      // scatter thePro and theTarg
      G4KineticTrackVector* secondaries = 0;
      G4bool block = false;
      G4double etot = 0;
      G4double eaft = 0;
      
      if(HIT==limit)
      {
        G4LorentzVector mom = theTarg->Get4Momentum();
        G4KineticTrack Target(theTarg->GetDefinition(), 0, theTarg->GetPosition(), mom);
        secondaries = theScattering.Scatter(*thePro->GetOrg(), Target);
	etot = (*thePro).GetMom().t()+mom.t();
        G4int u=0;
	G4int buffer = nHits;
	if(secondaries)
	{
	  for(u=0; u<secondaries->size(); u++)
	  {
	    if(secondaries->operator[](u)->Get4Momentum().vect().mag() < theTarget->GetFermiMomentum(theTarg->GetPosition()))
	    {
	      if(buffer!=0)
	      {
		buffer--;
	      }
	      else
	      {
		if(buffer<nHits) buffer++;
		block=true;
		limit = FLIGHT;
		//G4cout << "BLOCKED "
	        //       <<secondaries->operator[](u)->Get4Momentum().vect().mag()<<" "
		//       <<theTarget->GetFermiMomentum(theTarg->GetPosition())<<" "
		//       <<G4endl;
		break;
	      }
	    }
	  }
	}
	else
	{
	  limit = FLIGHT;
	}
	nHits = buffer;
      }
      result = new GhadFS(limit, aStep, thePro, theTarg);
      if(HIT==limit)
      {
	nHits++;
	for(G4int u=0; u<secondaries->size(); u++)
	{
          //G4cout << "We have  REAL hit ! "<<G4endl;
	  result->push_back( GhadTrack(secondaries->operator[](u) ) );
    	  theTarg->Hit();
	  eaft +=secondaries->operator[](u)->Get4Momentum().t();
	  G4double ekAft = secondaries->operator[](u)->Get4Momentum().t()
	                  -secondaries->operator[](u)->GetDefinition()->GetPDGMass();
	  G4double est = secondaries->operator[](u)->Get4Momentum().t()
	                 -secondaries->operator[](u)->GetActualMass();
	  if(ekAft<0) 
	  {
	    G4cerr << "ekin = "<<ekAft<<" "<<est<<G4endl;
	    G4cerr << "Scattering unphysical !!!"<<G4endl;
	  }
	}
//        G4cout << "EBALANCE PER SCATTER"<<etot << " "<<eaft<<" "<<eaft-etot<<G4endl;
      }
      if(GEO==limit&& G4UniformRand()>0.001) // @@@@ some protection needed
      {
	G4LorentzVector mom4 = thePro->GetMom();
	G4ThreeVector mom = mom4.vect();
	G4ThreeVector pos = thePro->GetPosition();
	G4double r = pos.mag();
	G4double pr = mom*pos/r; 

	thePro->Go(aStep);
	G4int dummy;
	dummy++;
	G4double potentialStep = theShellModel.GetPotentialStep(thePro);

	G4double newRadial = 0;;
	G4double qv = potentialStep*potentialStep 
                      - 2.*potentialStep*mom4.t() 
		      + pr*pr;

	if(qv <= 0.0) 
	{ 
	  newRadial = -pr;
	} else 
	{ 
	  newRadial = sqrt(qv)*pr/abs(pr);
	};

	G4double prr = (newRadial - pr) / r;  
	G4LorentzVector aNew4Mom(mom4.t(), mom4.vect()+prr*pos);
//	G4LorentzVector aNew4Mom(mom4);

	G4KineticTrack * aNew = 
	    new G4KineticTrack(thePro->GetDefinition(), 0, pos, aNew4Mom);
	GhadTrack theNew(aNew);
	theNew.Go(aStep);
	result->push_back(theNew);
	//G4cout << "BOUNDARY CROSSED"<<G4endl;
      }
      return result;
    }
    
    G4double GetOuterRadius() {return theTarget->GetOuterRadius();}
    G4double GetFermiMomentum(G4ThreeVector it)
    {
      return theTarget->GetFermiMomentum(it);
    }
    const std::vector<G4Nucleon *> & GetNucleons()
    {
      return theTarget->GetNucleons();
    }
    
    G4Pair<G4double, G4double> ChooseImpactXandY(G4double maxImpact)
    {
      return theTarget->ChooseImpactXandY(maxImpact);
    }
    
    G4int GetNHoles() {return nHits;}

  private:
    G4V3DNucleus * theTarget;
    GhadNuclearGeometry theShellModel;
    
    G4Scatterer theScattering;
    G4int nHits;
    
    G4int initialized;
};

#endif
