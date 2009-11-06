#include "G4AdjointDeuteron.hh"
#include "G4ParticleTable.hh"

// ######################################################################
// ###                           ADJOINT DEUTERON                     ###
// ######################################################################

G4AdjointDeuteron* G4AdjointDeuteron::theInstance = 0;

G4AdjointDeuteron* G4AdjointDeuteron::Definition()
{
  if (theInstance !=0) return theInstance;
 
  const G4String name = "adj_deuteron";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4AdjointIons* anInstance =  reinterpret_cast<G4AdjointIons*>(pTable->FindParticle(name));
  if (anInstance ==0)
  {
  // create particle
  //
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived      subType    anti_encoding
  //             excitation 
    anInstance = new G4AdjointIons(
                 name,   1.875613*GeV,       0.0*MeV,  -1.0*eplus,
                    2,              +1,             0,
                    0,               0,             0,
            "adjoint_nucleus",               0,            +2, 1000010020,
                 true,            -1.0,          NULL,
		false,           "static",          0,
                  0.0
              );

    // Magnetic Moment
    G4double mN = eplus*hbar_Planck/2./(proton_mass_c2 /c_squared);
    anInstance->SetPDGMagneticMoment( 0.857438230 * mN);

  }
  //No Anti particle registered
  anInstance->SetAntiPDGEncoding(0);

  theInstance = reinterpret_cast<G4AdjointDeuteron*>(anInstance);
  return theInstance;
}

G4AdjointDeuteron*  G4AdjointDeuteron::DeuteronDefinition()
{
  return Definition();
}

G4AdjointDeuteron*  G4AdjointDeuteron::Deuteron()
{
  return Definition();
}


