#include "precoreaction.h"


ClassImp(precoreaction)

precoreaction::precoreaction() : 
   ReactionNumber(-1), ProjectileA(-1), ProjectileZ(-1),
   ProjectileP(0.0,0.0,0.0,0.0), ProjectileM(0.0),
   TargetA(-1), TargetZ(-1), TargetM(0.0),
   CompoundU(0.0),CompoundP(0.0,0.0,0.0), 
   CompoundM(-1.0), NFragments(0)
{
    Fragments = new TClonesArray("precofragment",15);
}

precoreaction::precoreaction(const int n) : 
   ReactionNumber(-1), ProjectileA(-1), ProjectileZ(-1), 
   ProjectileP(0.0,0.0,0.0,0.0), ProjectileM(0.0),
   TargetA(-1), TargetZ(-1), TargetM(0.0),
   CompoundU(0.0),CompoundP(0.0,0.0,0.0),  
   CompoundM(-1.0), NFragments(0)
{
    Fragments = new TClonesArray("precofragment",n);
}


void 
precoreaction::AddFragment(const std::string& frag_name,
			   const int a,
			   const int z,
			   const double m,
			   const TLorentzVector& p,
			   const std::string& process)
{
    TClonesArray & frag = *Fragments;
    new (frag[NFragments++])
	precofragment(frag_name,a,z,m,p,process);

    return;
}


int 
precoreaction::CheckConservationA() const 
{
   int A = this->GetCompoundA();
   for (int i = 0; i < NFragments; i++)
   {
      precofragment * frag = (precofragment*)Fragments->At(i);
      A -= frag->GetA();
   }
   return A;
}

int 
precoreaction::CheckConservationZ() const 
{
   int Z = this->GetCompoundZ();
   for (int i = 0; i < NFragments; i++)
   {
      precofragment * frag = (precofragment*)Fragments->At(i);
      Z -= frag->GetZ();
   }
   return Z;
}


double 
precoreaction::CheckConservationE() const 
{
  double E = this->GetCompoundE();
  for (int i = 0; i < NFragments; i++)
    {
      precofragment * frag = (precofragment*)Fragments->At(i);
      E -= frag->GetTotalE();
    }
  return E;
}

TVector3
precoreaction::CheckConservationP() const
{
  TVector3 P0 = this->GetProjectileMomentum().Vect();
  TVector3 P1(0.0,0.0,0.0);
  for (int i = 0; i < NFragments; i++)
  {
     precofragment * frag = (precofragment*)Fragments->At(i);
     P1 +=  frag->GetMomentum().Vect();
  }
  return P0-P1;
}

double
precoreaction::CheckConservationPMag() const
{
   double P = this->GetProjectileMomentum().Vect().Mag();
   TVector3 P1(0.0,0.0,0.0);
   for (int i = 0; i < NFragments; i++)
   {
      precofragment * frag = (precofragment*)Fragments->At(i);
      P1 +=  frag->GetMomentum().Vect();
   }
   return P - P1.Mag();
}








