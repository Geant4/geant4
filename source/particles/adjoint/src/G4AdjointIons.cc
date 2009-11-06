#include <fstream>
#include <iomanip>

#include "G4AdjointIons.hh"

 // ######################################################################
// ###                           ADJOINT Ions                          ###
// #######################################################################

G4AdjointIons::G4AdjointIons(
       const G4String&     aName,        G4double            mass,
       G4double            width,        G4double            charge,   
       G4int               iSpin,        G4int               iParity,    
       G4int               iConjugation, G4int               iIsospin,   
       G4int               iIsospin3,    G4int               gParity,
       const G4String&     pType,        G4int               lepton,      
       G4int               baryon,       G4int               encoding,
       G4bool              stable,       G4double            lifetime,
       G4DecayTable        *decaytable , G4bool              shortlived,
       const G4String&     subType,
       G4int               anti_encoding,
       G4double            excitation      )
  : G4ParticleDefinition( aName,mass,width,charge,iSpin,iParity,
           iConjugation,iIsospin,iIsospin3,gParity,pType,
           lepton,baryon,encoding,stable,lifetime,decaytable,
           shortlived, subType, anti_encoding)
{
  // initialize excitation energy/level
   theExcitationEnergy = excitation;

   SetAtomicNumber( G4int(-GetPDGCharge()/eplus) );
   SetAtomicMass( GetBaryonNumber() );

   //G4cout << "G4AdjointIons::" << GetParticleName() << G4endl; 
}


G4AdjointIons::~G4AdjointIons()
{

  //G4cout << "G4AdjointIons::" << GetParticleName() << G4endl; 
}


G4AdjointIons* G4AdjointIons::IonsDefinition()
{
  return this;
}



