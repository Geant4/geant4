#ifndef PCTCrossSection_hh 
#define PCTCrossSection_hh 1

#include "globals.hh"
#include "G4VCrossSectionDataSet.hh"

class PCTTarget;
class PCTProjectile;


class PCTCrossSection
{
public:
    inline PCTCrossSection(PCTProjectile * projectile, PCTTarget * target);
    inline ~PCTCrossSection();

    G4double CalculateCrossSection(PCTProjectile * projectile, PCTTarget * target);
    inline G4double GetCrossSection() const;

private:
    inline PCTCrossSection();
    inline PCTCrossSection(const PCTCrossSection& right);
    inline PCTCrossSection& operator=(const PCTCrossSection& right);

    void CalculateNucleusXS(PCTProjectile * projectile, PCTTarget * target);
    void CalculateElementXS(PCTProjectile * projectile, PCTTarget * target);

    G4VCrossSectionDataSet * DataSetAdapter(const G4ParticleDefinition * particle);

private:
    G4double CrossSection;
};

#include "PCTCrossSection.icc"

#endif   // PCTCrossSection_hh 
