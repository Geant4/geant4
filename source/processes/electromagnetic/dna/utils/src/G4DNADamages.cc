#include "G4DNADamages.hh"
#include "G4UnitsTable.hh"

G4DNADamages* G4DNADamages::fInstance(0);

G4DNAIndirectHit::G4DNAIndirectHit(const G4String& baseName,
                                   const G4Molecule* molecule,
                                   const G4ThreeVector& position,
                                   G4double time) : G4VDNAHit(),
    fpMolecule(molecule)
{
    fBaseName = baseName;
    fPosition = position;
    fTime = time;
}


G4DNAIndirectHit::~G4DNAIndirectHit()
{
    ;
}

void G4DNAIndirectHit::Print()
{
    G4cout << "Reaction : " << fpMolecule->GetName() << " + " << fBaseName
           << " at position : " << G4BestUnit(fPosition,"Length")
           << " and time : " << G4BestUnit(fTime,"Time") << G4endl;
}


G4DNADamages* G4DNADamages::Instance()
{
    if(!fInstance) new G4DNADamages();

    return fInstance;
}

G4DNADamages::G4DNADamages() : G4VStateDependent()
{
    fInstance = this;
}


G4DNADamages::~G4DNADamages()
{
    for(int i = 0 ; i <(int) fIndirectHits.size() ; i++)
    {
        if(fIndirectHits[i])
            delete fIndirectHits[i];
    }
    fIndirectHits.clear();
}

G4bool G4DNADamages::Notify(G4ApplicationState /*requestedState*/)
{
    return true;
}

void G4DNADamages::Reset()
{
    for(int i = 0 ; i <(int) fIndirectHits.size() ; i++)
    {
        if(fIndirectHits[i])
            delete fIndirectHits[i];
    }
    fIndirectHits.clear();
}

void G4DNADamages::AddIndirectDamage(const G4String& baseName,
                                     const G4Molecule* molecule,
                                     const G4ThreeVector& position,
                                     G4double time)
{
    G4DNAIndirectHit* indirectHit  = 0;
    std::map<const G4Molecule, const G4Molecule*>::iterator it = fMolMap.find(*molecule);

    if(it == fMolMap.end())
    {
        G4Molecule* mol(0);
        fMolMap[*molecule] = (mol = new G4Molecule(*molecule));
        indirectHit = new G4DNAIndirectHit(baseName, mol, position, time);
    }
    else
    {
        indirectHit = new G4DNAIndirectHit(baseName, it->second, position, time);
    }
    fIndirectHits.push_back(indirectHit);
}
