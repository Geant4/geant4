#ifndef PCTProjectile_hh
#define PCTProjectile_hh 1

#include "globals.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleDefinition.hh"

#include "PCTProjectileDirection.hh"


class PCTProjectile
{
public:
    PCTProjectile(const G4String partName);
    PCTProjectile(const G4int a, const G4int z);
    inline ~PCTProjectile();

    G4bool operator==(const PCTProjectile& right) const;
    G4bool operator!=(const PCTProjectile& right) const;


    void SetDirection(const PCTProjectileDirection aDirection);

    G4bool SetKineticEnergy(const G4double KineticEnergy);

    inline G4bool SetDirAndKineticE(const G4double KineticEnergy, 
				    const PCTProjectileDirection& aDirection);

    inline G4bool SetDirAndKineticE(const PCTProjectileDirection& aDirection, 
				    const G4double KineticEnergy);

    const G4LorentzVector& GetMomentum() const;

    inline PCTProjectileDirection GetDirection() const;

    inline G4ParticleDefinition * GetParticle() const;

    inline G4double GetMass() const;

    inline G4int GetZ() const;

    inline G4int GetA() const;
    
    inline G4ParticleDefinition * GetDefinition() const;

private:
    PCTProjectile() {};
    inline PCTProjectile(const PCTProjectile & p);
    inline const PCTProjectile& operator=(const PCTProjectile &right);

    inline void ResetMomentum();

private:
    G4ParticleDefinition * theParticle;
    PCTProjectileDirection dir;
    G4LorentzVector momentum;
};

#include "PCTProjectile.icc"

#endif // PCTProjectile_hh


