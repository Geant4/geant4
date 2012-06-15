#ifndef G4DNAMolecularMaterial_HH
#define G4DNAMolecularMaterial_HH

#include "globals.hh"
#include "G4ios.hh"
#include <map>
#include <vector>
#include "G4VStateDependent.hh"

class G4Material;

struct CompareMaterial
{
    // If the materials derives from a base material,
    // it should be able to find the derived material using the base material.
    bool operator() (const G4Material* mat1, const G4Material* mat2) const;
};

typedef std::map<const G4Material*, double,CompareMaterial> ComponentMap;

// G4DNAMolecularMaterial is initialized when G4ApplicationState == G4State_Idle

class G4DNAMolecularMaterial : public G4VStateDependent
{
public:
    static G4DNAMolecularMaterial* Instance();
    void DeleteInstance();
    void Initialize();

    virtual G4bool Notify(G4ApplicationState requestedState) ;

    inline const std::vector<ComponentMap>* GetMassFractionTable() const;
    inline const std::vector<ComponentMap>* GetDensityTable() const;
//    const std::vector<double>* GetMassFractionTableFor(const G4Material*) const;
    const std::vector<double>* GetDensityTableFor(const G4Material*) const;
    const std::vector<double>* GetNumMolPerVolTableFor(const G4Material*) const;

protected :
    static G4DNAMolecularMaterial* fInstance;
    G4DNAMolecularMaterial();
    G4DNAMolecularMaterial(const G4DNAMolecularMaterial& right);
    G4DNAMolecularMaterial& operator=(const G4DNAMolecularMaterial&);
    virtual ~G4DNAMolecularMaterial();
    void Create();
    void InitializeNumMolPerVol();
    void InitializeDensity();
    void RecordMolecularMaterial(G4Material* parentMaterial, G4Material* molecularMaterial, G4double fraction);
    void SearchMolecularMaterial(G4Material* parentMaterial, G4Material* material, double currentFraction);

    void AddMaterial(const G4Material*, double fraction);

    void PrintNotAMolecularMaterial(const char* methodName, const G4Material* lookForMaterial) const;

    std::vector<ComponentMap>* fpCompFractionTable;
    std::vector<ComponentMap>* fpCompDensityTable;
    std::vector<ComponentMap>* fpCompNumMolPerVolTable;

    mutable std::map<const G4Material*,std::vector<double>*,CompareMaterial> fAskedDensityTable;
    mutable std::map<const G4Material*,std::vector<double>*,CompareMaterial> fAskedNumPerVolTable;
    mutable std::map<const G4Material*,bool,CompareMaterial> fWarningPrinted;

    G4bool fIsInitialized;
};

inline const std::vector<ComponentMap> *G4DNAMolecularMaterial::GetMassFractionTable() const
{
    return fpCompFractionTable;
}

inline const std::vector<ComponentMap>* G4DNAMolecularMaterial::GetDensityTable() const
{
    return fpCompDensityTable;
}

#endif // G4DNAMolecularMaterial_HH
