#ifndef precoreaction_h
#define precoreaction_h

#include "precofragment.h"

#include <string>

#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector3.h"

class precoreaction
{
private:
  int      ReactionNumber;

  double   CrossSection;

  int      ProjectileA;
  int      ProjectileZ;
  TLorentzVector ProjectileP;
  double   ProjectileM;

  int      TargetA;
  int      TargetZ;
  double   TargetM;

  double   CompoundU;
  TLorentzVector CompoundP;
  double   CompoundM;
  
  int NFragments;
  TClonesArray * Fragments; //->

public:
  precoreaction();
  precoreaction(const int n);
  virtual ~precoreaction() { Fragments->Delete(); }


  void SetReactionNumber(const int n) { ReactionNumber = n; }

  void SetCrossSection(const double xs) { CrossSection = xs; }

  // Projectile 
  void SetProjectileA(const int a) { ProjectileA = a; }
  void SetProjectileZ(const int z) { ProjectileZ = z; }
  void SetProjectileMomentum(const TLorentzVector & p) { ProjectileP = p; }
  void SetProjectileMass(const double m) { ProjectileM = m; }
  // Target
  void SetTargetA(const int a) { TargetA = a; }
  void SetTargetZ(const int z) { TargetZ = z; }
  void SetTargetMass(const double m) { TargetM = m; }
  // Compopund
  void SetCompoundU(const double u) { CompoundU = u; }
  void SetCompoundMomentum(const TLorentzVector& p) { CompoundP = p; }
  void SetCompoundMass(const double m) { CompoundM = m; }

  void Clear()
  {
    if ( Fragments ) Fragments->Delete();
    NFragments = 0;
    return;
  }

  void AddFragment(const std::string& frag_name,
		   const int a,
		   const int z,
		   const double m,
		   const TLorentzVector& p,
		   const std::string& process);

  void AllocateFor(const int n)
  {
    Fragments->Expand(n);
    return;
  }


  int GetNumOfFragments() const
  { return Fragments->GetLast()+1; }
  
  int GetReactionNumber() const { return ReactionNumber; }
  double GetCrossSection() const { return CrossSection; }

  int GetProjectileA() const { return ProjectileA; }
  int GetProjectileZ() const { return ProjectileZ; }
  double GetProjectileKineticE() const { return ProjectileP.E()-ProjectileM; }
  TLorentzVector GetProjectileMomentum() const { return ProjectileP; }
  double GetProjectileMass() const { return ProjectileM; }

  int GetTargetA() const { return TargetA; }
  int GetTargetZ() const { return TargetZ; }
  double GetTargetMass() const { return TargetM; }
  int GetCompoundA() const { return this->GetProjectileA() + this->GetTargetA(); }
  int GetCompoundZ() const { return this->GetProjectileZ() + this->GetTargetZ(); }
  double GetCompoundU() const { return CompoundU; }
  TLorentzVector GetCompoundMomentum() const { return CompoundP; }
  double GetCompoundE() const { return CompoundP.E(); }
  double GetCompoundMass() const { return CompoundM; }

  int CheckConservationA() const;
  int CheckConservationZ() const;
  double CheckConservationE() const;
  double CheckConservationPMag() const;
  TVector3 CheckConservationP() const;

  bool IsNatural() const { return (TargetM == 0); }

  TVector3 GetBoostToCM() const { return -CompoundP.BoostVector();}
  
  TClonesArray * GetFragmentsList() const { return Fragments; } 

  precofragment * GetFragment(const int n) const
  {
    if ( n > -1 && n < NFragments )
      return (precofragment*)Fragments->At(n);
    else 
      return (precofragment*)0;
  }

  ClassDef(precoreaction,1)

};
#endif
