//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

#ifndef G4IonsShenCrossSection_h
#define G4IonsShenCrossSection_h
//
// Class Description
// Implementation of formulas 
// Shen et al. Nuc. Phys. A 491 130 (1989); 
// Total Reaction Cross Section for Heavy-Ion Collisions 
//
// Class Description - End
// 18-Sep-2003 First version is written by T. Koi
// 12-Nov-2003 Set upper limit at 10 GeV/n
// 12-Nov-2003 Insted of the lower limit,
//             0 is returned to a partilce with energy lowae than 10 MeV/n

#include "globals.hh"
#include "G4Proton.hh"

#include "G4VCrossSectionDataSet.hh"

class G4IonsShenCrossSection : public G4VCrossSectionDataSet
{
   public:
      G4IonsShenCrossSection ():
         upperLimit ( 10 * GeV ),
         lowerLimit ( 10 * MeV ),
         r0 ( 1.1 )
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
   {G4cout << "G4IonsShenCrossSection: uses Shen formula"<<G4endl;}

   private:
      const G4double upperLimit;
      const G4double lowerLimit; 
      const G4double r0;

      G4double calEcmValue ( const G4double , const G4double , const G4double ); 
      G4double calCeValue ( const G4double ); 
};

#endif
