#ifndef G4IonsKoxCrossSection_h
#define G4IonsKoxCrossSection_h
//
// Class Description
// Implementation of Kox formulas 
// Kox et al. Phys. Rev. C 35 1678 (1987); 
// Total Reaction Cross Section for Nucleus-nucles reactions.
//
// Class Description - End
// 18-Sep-2003 First version is written by T. Koi
// 12-Nov-2003 Set upper limit at 10 GeV/n
// 12-Nov-2003 Insted of the lower limit, 
//             0 is returned to a partilce with energy lowae than 10 MeV/n 

#include "globals.hh"
#include "G4Proton.hh"

#include "G4VCrossSectionDataSet.hh"

class G4IonsKoxCrossSection : public G4VCrossSectionDataSet
{
   public:
      G4IonsKoxCrossSection():
         upperLimit ( 10 * GeV ), 
         lowerLimit ( 10 * MeV ), 
         r0 ( 1.1 * fermi ),
         rc ( 1.3 * fermi )
      {
      }

   
   virtual
   G4bool IsApplicable(const G4DynamicParticle* aDP, const G4Element*)
   {
      G4int baryonNumber = aDP->GetDefinition()->GetBaryonNumber();
      G4double kineticEnergy = aDP->GetKineticEnergy(); 
      if ( kineticEnergy / baryonNumber <= upperLimit ) 
         return true;
      return false;
   }

   virtual
   G4double GetCrossSection(const G4DynamicParticle*, 
                            const G4Element*, G4double aTemperature);

   virtual
   void BuildPhysicsTable(const G4ParticleDefinition&)
   {}

   virtual
   void DumpPhysicsTable(const G4ParticleDefinition&) 
   {G4cout << "G4IonsKoxCrossSection: uses Kox formula"<<G4endl;}

   private:
      const G4double upperLimit; 
      const G4double lowerLimit; 
      const G4double r0;
      const G4double rc;

      G4double calEcm ( G4double , G4double , G4double ); 
      G4double calCeValue ( G4double ); 
};

#endif
