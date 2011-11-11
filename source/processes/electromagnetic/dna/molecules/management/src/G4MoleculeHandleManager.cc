#include "G4MoleculeHandleManager.hh"
#include "G4Molecule.hh"

using namespace std;

G4MoleculeHandleManager* G4MoleculeHandleManager::fInstance (0);

G4MoleculeHandleManager::G4MoleculeHandleManager()
{
//    G4cout << "G4MoleculeHandleManager::G4MoleculeHandleManager()" << G4endl;
}

G4bool G4MoleculeHandleManager::CompMoleculePointer::operator()(const G4Molecule* mol1, const G4Molecule* mol2) const
{
    return (*mol1) < (*mol2);
}

G4MoleculeHandleManager::~G4MoleculeHandleManager()
{
    if(!fMoleculeHandle.empty())
    {
        MoleculeHandleMap::iterator it = fMoleculeHandle.begin();
        for( ; it != fMoleculeHandle.end() ; it++)
        {
            it->second.reset();
        }
    }
}

void G4MoleculeHandleManager::DeleteInstance()
{
    if(fInstance)
    {
        delete fInstance;
        fInstance = 0;
    }
}

G4MoleculeHandleManager* G4MoleculeHandleManager::Instance()
{
    if(!fInstance)
    {
        fInstance = new G4MoleculeHandleManager;
    }
    return fInstance;
}

G4MoleculeHandle G4MoleculeHandleManager::GetMoleculeHandle(const G4Molecule* molecule)
{
    MoleculeHandleMap::iterator it = fMoleculeHandle.find(molecule);
    G4MoleculeHandle molHandle;

    if(it != fMoleculeHandle.end())
    {
        molHandle = G4MoleculeHandle(it->second);
    }
    else
    {
        molHandle = G4MoleculeHandle(molecule);
        fMoleculeHandle.insert(make_pair(molecule, G4MoleculeHandle(molHandle)));
    }

    return molHandle;
}
