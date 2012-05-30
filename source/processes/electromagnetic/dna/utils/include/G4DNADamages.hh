#ifndef G4DNADAMAGES_HH
#define G4DNADAMAGES_HH 1

#include "G4VStateDependent.hh"
#include "G4Molecule.hh"

class G4VDNAHit
{
public :
    G4VDNAHit(){;}
    virtual ~G4VDNAHit(){;}
};

class G4DNAIndirectHit : public G4VDNAHit
{
public :
    G4DNAIndirectHit(const G4String& baseName, const G4Molecule* molecule,
                     const G4ThreeVector& position, G4double time);
    virtual ~G4DNAIndirectHit();

    inline const G4Molecule* GetMolecule() {return fpMolecule;}
    inline const G4ThreeVector& GetPosition() {return fPosition;}
    inline const G4String& GetBaseName() {return fBaseName;}
    inline double GetTime() {return fTime;}

    void Print();

protected :
    const G4Molecule* fpMolecule;
    G4ThreeVector fPosition;
    G4double fTime;
    G4String fBaseName;
};

class G4DNADamages : public G4VStateDependent
{
    static G4DNADamages* fInstance;
public:
    G4DNADamages();
    static G4DNADamages* Instance();
    virtual ~G4DNADamages();

    virtual G4bool Notify(G4ApplicationState requestedState) ;

    void Reset();

    //void AddDirectDamage();
    void AddIndirectDamage(const G4String& baseName,const G4Molecule* molecule,
                           const G4ThreeVector& position, double time);

    inline const std::vector<G4DNAIndirectHit*>* GetIndirectHits();

protected :
    std::vector<G4DNAIndirectHit*> fIndirectHits;
    std::map<const G4Molecule, const G4Molecule*> fMolMap;
};

inline const std::vector<G4DNAIndirectHit*>* G4DNADamages::GetIndirectHits()
{
    return &fIndirectHits;
}

#endif // G4DNADAMAGES_HH
