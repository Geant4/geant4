#include "Ghad_hp.hh"
#include "Randomize.hh"

void Ghad_hp::Stepping()
{
	G4double maxStep = 100.*fermi/c_light;
	int time = 0;
	G4double tCut = 30*MeV;
G4double maxE = 0;
	do
	{
    cout << "time-stepping "<<time<<endl;
		G4double aStep = maxStep;
		// scatter, decay, boundary, blocking, etc..
		GhadFS * theAction = theTarget.YouAreTraversedBy(theProjectiles, aStep); 

		aStep = theAction->GetTime();		
    cout << "stepping "<<aStep<<endl;
		GhadParticles::iterator p;
		for(p=theProjectiles.begin(); p!=theProjectiles.end(); p++)
		{
		  p->Go(aStep);
		}

  	if(!theAction->ImVoid())
  	{
    	theProjectiles.erase(find(theProjectiles.begin(), theProjectiles.end(),
			                          theAction->GetProjectile()));
    	GhadFS::iterator i;
    	for(i=theAction->begin(); i!=theAction->end(); i++)
    	{
      	theProjectiles.push_back(*i);
    	}
  	}
	  maxE = 0;
		for(p=theProjectiles.begin(); p!=theProjectiles.end(); p++)
		{
		  if(p->GetPosition().mag()>4.5*theTarget.GetOuterRadius())
			{
			  
				theEscaped.push_back(*p);
				p=theProjectiles.erase(find(theProjectiles.begin(), theProjectiles.end(), p));
				p--;
			}
			else
			{
			  G4double tmp = p->GetMom().t() - p->GetActualMass();
				if(tmp>maxE) maxE = tmp;
			}
		} 
	}
	while(maxE>0.1*MeV);
	// it is all in the angular distribution; make it a paremetrization,
	// with a steepness parameter....and tune.
	
	for(G4int it=0; it<theEscaped.size(); it++)
  {
    G4LorentzVector v = theEscaped.operator[](it).GetMom();
		G4ParticleDefinition * p = theEscaped.operator[](it).GetDefinition();
		G4double pfermi = theTarget.GetFermiMomentum(theEscaped.operator[](it).GetPosition());
//		if(v.vect().mag()>pfermi)
		if(theProjectiles.size()+theEscaped.size()>1)
		  G4cout << "SEC "
             << v.vect().x()<<" "
             << v.vect().y()<<" "
             << v.vect().z()<<" "
				  	 << v.t()-p->GetPDGMass()<<" "
             << p->GetPDGEncoding()<<" "
	           <<G4endl;
  }

}
