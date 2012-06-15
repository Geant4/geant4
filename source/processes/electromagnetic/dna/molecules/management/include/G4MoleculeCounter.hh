#ifndef G4MoleculeCounter_h
#define G4MoleculeCounter_h

#include "G4Molecule.hh"
#include <map>
#include "G4ParticleDefinition.hh"

using namespace std;

struct compDoubleWithPrecision
{
    bool operator() (const double& a, const double& b) const
    {
        if(fabs(a - b) < fPrecision)
        {
            return false;
        }
        else
        {
            return a < b;
        }
    }

    static double fPrecision ;
};

typedef map<G4double, G4int, compDoubleWithPrecision> NbMoleculeAgainstTime;

class G4MoleculeCounter
{
private:
    G4MoleculeCounter();
    ~G4MoleculeCounter(){;}
    static G4MoleculeCounter* fpInstance;
    typedef std::map<const G4Molecule, NbMoleculeAgainstTime> CounterMapType;

    CounterMapType fCounterMap;

    std::map<const G4MoleculeDefinition*, G4bool> fDontRegister ;
    G4bool fUse;

    G4int fVerbose ;

public:
    static void DeleteInstance();

    static  G4MoleculeCounter* GetMoleculeCounter();
    inline const NbMoleculeAgainstTime& GetNbMoleculeAgainstTime(const G4Molecule &molecule);
    vector<G4Molecule> RecordMolecules();
    void AddAMoleculeAtTime(const G4Molecule&, G4double);
    void RemoveAMoleculeAtTime(const G4Molecule&, G4double);

    // inline void DontRegister(G4MoleculeID);
    inline void DontRegister(const G4MoleculeDefinition*);
    inline void ResetDontRegister();

    void Use(G4bool flag)
    {
        fUse=flag;
    }
    G4bool InUse()
    {
        return fUse;
    }

    inline void SetVerbose(G4int);
    inline G4int GetVerbose();

    void ResetCounter();
};

inline void G4MoleculeCounter::ResetCounter()
{
    fCounterMap.clear();
}

inline const NbMoleculeAgainstTime& G4MoleculeCounter::GetNbMoleculeAgainstTime(const G4Molecule& molecule)
{
    return fCounterMap[molecule];
}

inline void G4MoleculeCounter::SetVerbose(G4int level)
{
    fVerbose = level;
}

inline G4int G4MoleculeCounter::GetVerbose()
{
    return fVerbose ;
}

inline void G4MoleculeCounter::DontRegister(const G4MoleculeDefinition* molDef)
{
    fDontRegister[molDef] = true ;
}

inline void G4MoleculeCounter::ResetDontRegister()
{
    fDontRegister.clear();
}

#endif
