

#ifndef FluoDataSet_hh
#define FluoDataSet_hh 1

#include "globals.hh"
#include "G4DataVector.hh"
#include "G4VEMDataSet.hh"

class G4VDataSetAlgorithm;

class XrayFluoDataSet : public G4VEMDataSet {

public:

  XrayFluoDataSet(G4int Z,
		  G4DataVector* points, 
	      G4DataVector* values,
	      const G4VDataSetAlgorithm* interpolation,
	      G4double unitE = MeV, G4double unitData = barn);

  XrayFluoDataSet(G4int Z,
		  const G4String& dataFile,
	      const G4VDataSetAlgorithm* interpolation,
	      G4double unitE = MeV, G4double unitData = barn);

  ~XrayFluoDataSet();
 
  G4double FindValue(G4double e, G4int id = 0) const;
  
  void PrintData() const;

  const G4DataVector& GetEnergies(G4int i) const { return *energies; }
  const G4DataVector& GetData(G4int i) const { return *data; }

private:

  void LoadData(const G4String& dataFile);
  G4int z;
  G4int FindBinLocation(G4double energy) const;

  G4DataVector* energies; // Owned pointer
  G4DataVector* data;     // Owned pointer

  const G4VDataSetAlgorithm* algorithm; // Not owned pointer 
  
  G4double unit1;
  G4double unit2;

  size_t numberOfBins;

};
 
#endif
