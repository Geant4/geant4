#include "G4VMolecularDecayDisplacer.hh"
#include "G4Molecule.hh"

DisplacementType G4VMolecularDecayDisplacer::Last = 0;
const DisplacementType G4VMolecularDecayDisplacer::NoDisplacement = G4VMolecularDecayDisplacer::AddDisplacement();

G4VMolecularDecayDisplacer::G4VMolecularDecayDisplacer()
{
    fVerbose = 0;
}

G4VMolecularDecayDisplacer::~G4VMolecularDecayDisplacer()
{
    ;
}

DisplacementType G4VMolecularDecayDisplacer::AddDisplacement()
{
    DisplacementType output = Last;
    Last++;
    return output;
}
