#include "PCTCompositeNucleus.hh"

#include "G4Proton.hh"
#include "G4Fragment.hh"



//#define test

#ifdef test
#include "G4ParticleTable.hh"
#endif

const G4Fragment * PCTCompositeNucleus::GetNewCNucleus()
{
  if ( theTarget->IsNatural() || !lastCNucleus )
    {
      G4int Ap = theProjectile->GetA();
      G4int Zp = theProjectile->GetZ();
      G4int At = theTarget->GetA();
      G4int Zt = theTarget->GetZ();
      
      G4LorentzVector compoundMomentum = theTarget->GetMomentum() + theProjectile->GetMomentum();
      
      if (lastCNucleus) delete lastCNucleus;
      lastCNucleus = new G4Fragment(At+Ap,Zt+Zp,compoundMomentum);
      if (randExcitons)
	{
	  // This is for projectile=nucleon
	  G4int extraExciton  = Ap * G4int(0.5+G4UniformRand()); 
	  lastCNucleus->SetNumberOfParticles(1+extraExciton);
	  lastCNucleus->SetNumberOfHoles(extraExciton);
	  if (G4UniformRand() < G4double(Zt)/G4double(At))
	    {
	      lastCNucleus->SetNumberOfCharged(extraExciton);
	    }
	  else
	    {
	      lastCNucleus->SetNumberOfCharged(0);
	    }
	}
      else
	{
	  lastCNucleus->SetNumberOfParticles(particles);
	  lastCNucleus->SetNumberOfHoles(holes);
	  lastCNucleus->SetNumberOfCharged(charged);
	}
      //
      lastCNucleus->SetParticleDefinition(G4Proton::ProtonDefinition());
#ifdef test
      G4double compoundMass = G4ParticleTable::GetParticleTable()->
	GetIonTable()->GetIonMass(Zt+Zp,At+Ap);
      G4LorentzVector tmp(compoundMomentum.vect(),
			  sqrt(compoundMomentum.vect().mag2()+compoundMass*compoundMass));
      tmp -= compoundMomentum;
      G4double U = -tmp.e();
      
      G4cout << "Excitation Energy = " << U*MeV << '\n';
      G4cout << *lastCNucleus << '\n'
	     << "U = " << lastCNucleus->GetExcitationEnergy() << '\n';

#endif
    }
  return lastCNucleus;
}


