#ifndef expdata_h
#define expdata_h  

#include <iostream>
#include <sstream>
#include "TObject.h"
#include "TString.h"

class expdata : public TObject
{
  // Base class 
 public:
  // Default constructor
  inline expdata();

  // Copy constructor
  expdata(const expdata &);
  // Destructor
  inline virtual ~expdata();
  // operator =
  const expdata & operator=(const expdata & rhs);
    
  inline void SetTargetA(const Int_t anA);
  inline Int_t GetTargetA() const;
    
  inline void SetTargetZ(const Int_t aZ);
  inline Int_t GetTargetZ() const;
    
  inline void SetTargetSymbol(const TString& s);
  inline void SetTargetSymbol(const Char_t * s);
  inline TString GetTargetSymbol() const;
  
  inline void SetProjectileA(const Int_t anA);
  inline Int_t GetProjectileA() const;
    
  inline void SetProjectileZ(const Int_t aZ);
  inline Int_t GetProjectileZ() const;
    
  inline void SetProjectileSymbol(const TString& s);
  inline void SetProjectileSymbol(const Char_t * s);
  inline TString GetProjectileSymbol() const;
    
  inline void SetProjectileEnergy(const Double_t e);
  inline Double_t GetProjectileEnergy() const;
    
  inline void SetExforEntryCode(const TString& code);
  inline void SetExforEntryCode(const Char_t * code);
  inline TString GetExforEntryCode() const;
    
  inline void SetExforReaction(const TString& r );
  inline void SetExforReaction(const Char_t * r );
  inline TString GetExforReaction() const;
    
  inline void SetParticleA(const Int_t anA);
  inline Int_t GetParticleA() const;
    
  inline void SetParticleZ(const Int_t aZ);
  inline Int_t GetParticleZ() const;
    
  inline void SetParticleSymbol(const TString& s);
  inline void SetParticleSymbol(const Char_t * s);
  inline TString GetParticleSymbol() const;
    
  inline void SetCM();
  inline void UnsetCM();
  inline bool IsCM() const;

  inline void SetCutoff(const double cut);
  inline double GetCutoff() const;
  inline bool ThereIsCutoff() const;
    
  inline void ShowYourSelf() const;
  void ShowYourSelf(std::ostringstream & os) const;
  
private:
  
  Int_t targetA;
  Int_t targetZ;
  TString targetSymbol;
  Int_t projectileA;
  Int_t projectileZ;
  Double_t projectileEnergy;
  TString projectileSymbol;
  TString entry;
  TString reaction;
  Int_t particleA;
  Int_t particleZ;
  TString particleSymbol;
  bool CM;
  double cutoff;
  
  ClassDef(expdata,1)
};

#include "expdata.icc"

#endif
