#include "G4AdjointHe3.hh"
#include "G4ParticleTable.hh"

// ######################################################################
// ###                           ADJOINT He3                          ###
// ######################################################################

G4AdjointHe3* G4AdjointHe3::theInstance = 0;

G4AdjointHe3* G4AdjointHe3::Definition()
{
 
  if (theInstance !=0) return theInstance;
  const G4String name = "adj_He3";
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
                 name,   2.80923*GeV,       0.0*MeV,  -2.0*eplus,
                    1,              +1,             0,
                    0,               0,             0,
            "adjoint_nucleus",               0,            +3, 1000020030,
                 true,            -1.0,          NULL,
		false,        "static",             0,
                  0.0
              );
 
    // Magnetic Moment
    G4double mN = eplus*hbar_Planck/2./(proton_mass_c2 /c_squared);
    anInstance->SetPDGMagneticMoment( -2.12762485 * mN);

  }
  //No Anti particle registered
  anInstance->SetAntiPDGEncoding(0);

  theInstance = reinterpret_cast<G4AdjointHe3*>(anInstance);
  return theInstance;
}

G4AdjointHe3*  G4AdjointHe3::He3Definition()
{
  return Definition();
}

G4AdjointHe3*  G4AdjointHe3::He3()
{
  return Definition();
}


