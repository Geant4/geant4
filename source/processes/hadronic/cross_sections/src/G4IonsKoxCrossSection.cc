// 18-Sep-2003 First version is written by T. Koi

#include "G4IonsKoxCrossSection.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"

G4double G4IonsKoxCrossSection::
GetCrossSection(const G4DynamicParticle* aParticle, 
                const G4Element* anElement, G4double )
{
   G4double xsection = 0.0;

   G4int At = int ( anElement->GetN() + 0.5 );
   G4int Zt = int ( anElement->GetZ() + 0.5 );  

   G4int Ap = aParticle->GetDefinition()->GetBaryonNumber();
   G4int Zp = int ( aParticle->GetDefinition()->GetPDGCharge() / eplus + 0.5); 
 
   G4double one_third = 1.0 / 3.0;

   G4double cubicrAt = pow ( At , one_third );  
   G4double cubicrAp = pow ( Ap , one_third );  


   G4double Bc = Zt * Zp / ( ( rc / fermi ) * (  cubicrAp + cubicrAt ) );   // rc divide fermi

   G4double Rvol = r0 * (  cubicrAp + cubicrAt );

   G4double ke_per_N = aParticle->GetKineticEnergy() / At;
   G4double c = calCeValue ( ke_per_N / MeV  );  

   G4double a = 1.85;
   G4double Rsurf = r0 * ( a * cubicrAp * cubicrAt / ( cubicrAp + cubicrAt ) - c); 
   G4double D = 5 * ( At - 2 * Zt ) * Zp / ( Ap * At );
   Rsurf = Rsurf + D * fermi;  // multiply D by fermi 

   G4double Rint = Rvol + Rsurf;

   G4double targ_mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass( Zt , At ); 
   G4double proj_mass = aParticle->GetMass(); 
   G4double proj_momentum = aParticle->GetMomentum().mag(); 

   G4double Ecm = calEcm ( proj_mass , targ_mass , proj_momentum ); 

   xsection = pi * Rint * Rint * ( 1 - Bc / ( Ecm / MeV ) );
  
   return xsection; 
}

G4double G4IonsKoxCrossSection::calEcm ( G4double mp , G4double mt , G4double Plab )
{
   G4double Elab = sqrt ( mp * mp + Plab * Plab );
   G4double Ecm = sqrt ( mp * mp + mt * mt + 2 * Elab * mt );
   G4double Pcm = Plab * mt / Ecm;
   G4double KEcm = sqrt ( Pcm * Pcm + mp * mp ) - mp;
   return KEcm;
}

G4double G4IonsKoxCrossSection::calCeValue( const G4double ke )
{
   // Calculate c value 
   // This value is indepenent from projectile and target particle 
   // ke is projectile kinetic energy per nucleon in the Lab system with MeV unit 
   // fitting function is made by T. Koi 
   // There are no data below 30 MeV/n in Kox et al., 

   G4double Ce; 
   G4double log10_ke = log10 ( ke );   
   if ( log10_ke > 1.5 ) 
   {
      Ce = - 10.0 / pow ( log10_ke , 5 ) + 2.0;
   }
   else
   {
      Ce = ( - 10.0 / pow ( 1.5 , 5 ) + 2.0 ) / pow ( 1.5 , 3 ) * pow ( log10_ke , 3 );

   }
   return Ce;
}
