#ifndef TiaraSampledEnergy_hh
#define TiaraSampledEnergy_hh TiaraSampledEnergy_hh

#include "TiaraVSourceEnergyGenerator.hh"
#include <memory>

namespace AIDA {
  class IHistogram1D;
  class ITree;
}

class TiaraSampledEnergy : public TiaraVSourceEnergyGenerator {
public:
  TiaraSampledEnergy(const G4String &eng,
		     G4double minEnergyCut,
		     const G4String &sourceTree,
		     const G4String &nameExt);
  ~TiaraSampledEnergy();
  TiaraSampledEnergy(const TiaraSampledEnergy& rhs);


  virtual G4double GetEnergy();
  virtual TiaraVSourceEnergyGenerator *Clone() const;

  TiaraSampledEnergy& operator=(const TiaraSampledEnergy& rhs);
  
private:
  std::auto_ptr<AIDA::ITree> fTree;
  std::auto_ptr<AIDA::IHistogram1D> fSampleHisto;
  G4double fScale;
  G4double fMinEnergyCut;
};


#endif
