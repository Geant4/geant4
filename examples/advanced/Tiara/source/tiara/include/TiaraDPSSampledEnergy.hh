#ifndef TiaraDPSSampledEnergy_hh
#define TiaraDPSSampledEnergy_hh TiaraDPSSampledEnergy_hh

#ifdef G4ANALYSIS_USE

#include "TiaraVSourceEnergyGenerator.hh"
#include <memory>
#include "g4std/map"

namespace AIDA {
  class IDataPointSet;
  class ITree;
}

class TiaraDPSSampledEnergy : public TiaraVSourceEnergyGenerator {
public:
  TiaraDPSSampledEnergy(const G4String &eng,
		     G4double minEnergyCut,
		     const G4String &sourceTree,
		     const G4String &nameExt = G4String(""));
  ~TiaraDPSSampledEnergy();
  TiaraDPSSampledEnergy(const TiaraDPSSampledEnergy& rhs);


  virtual G4double GetEnergy();
  virtual TiaraVSourceEnergyGenerator *Clone() const;

  TiaraDPSSampledEnergy& operator=(const TiaraDPSSampledEnergy& rhs);
  
private:

  void getBounds(G4int &cL, G4int &cH, G4double v);

  std::auto_ptr<AIDA::ITree> fTree;
  std::auto_ptr<AIDA::IDataPointSet> fSampleDPS;
  G4double fMinEnergyCut;
  G4std::map<int, double> fEnergy_Flux;
  G4double fMaxProb;
  G4double fMinE;
  G4double fMaxE;
};

#endif
#endif
