
#include "G4AdjointAlpha.hh"
#include "G4ParticleTable.hh"

// ######################################################################
// ###                           ADJOINT ALPHA                                ###
// ######################################################################

G4AdjointAlpha* G4AdjointAlpha::theInstance = 0;

G4AdjointAlpha* G4AdjointAlpha::Definition()
{
  if (theInstance !=0) return theInstance;
  const G4String name = "adj_alpha";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4AdjointIons* anInstance = reinterpret_cast<G4AdjointIons*>(pTable->FindParticle(name));
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
                 name,   3.727417*GeV,       0.0*MeV,  -2.0*eplus,
                    0,              +1,             0,
                    0,               0,             0,
            "adjoint_nucleus",               0,            +4,  1000020040,
                 true,            -1.0,          NULL,
		 false,       "static",          0,
                 0.0
               );

  }

  //No Anti particle registered
  anInstance->SetAntiPDGEncoding(0);

  theInstance = reinterpret_cast<G4AdjointAlpha*>(anInstance);
  return theInstance;
}

G4AdjointAlpha*  G4AdjointAlpha::AlphaDefinition()
{
  return Definition();
}

G4AdjointAlpha*  G4AdjointAlpha::Alpha()
{
  return Definition();
}


