#ifndef PCTCompositeNucleus_hh 
#define PCTCompositeNucleus_hh 1

#include "PCTProjectile.hh"
#include "PCTTarget.hh"

class G4Fragment;

class PCTCompositeNucleus
{
public:
  inline PCTCompositeNucleus(PCTProjectile* projectile, PCTTarget* target);

  PCTCompositeNucleus(PCTProjectile* projectile, PCTTarget* target, 
		      const G4int h, const G4int p, const G4int c);

  inline void SetExcitons(const G4int h, const G4int p, const G4int c);
  inline void RandomizeExcitons(const G4bool op = true);
  inline G4bool AreExcitonsRandomized() const;

  const G4Fragment * GetNewCNucleus();
  inline const G4Fragment * GetLastCNucleus() const;

  inline const PCTProjectile * GetProjectile() const;
  inline const PCTTarget * GetTarget() const;

private:
  inline PCTCompositeNucleus();
  inline PCTCompositeNucleus(const PCTCompositeNucleus& right);
  inline PCTCompositeNucleus& operator=(const PCTCompositeNucleus& right);
  
private:

  PCTProjectile * theProjectile;
  PCTTarget *     theTarget;
  G4bool          randExcitons;
  G4int           holes;
  G4int           particles;
  G4int           charged;
  G4Fragment *    lastCNucleus;
};

#include "PCTCompositeNucleus.icc"

#endif // PCTCompositeNucleus


