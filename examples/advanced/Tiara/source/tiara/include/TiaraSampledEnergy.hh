// $Id: TiaraSampledEnergy.hh,v 1.3 2003-06-16 17:06:46 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
//
// Class TiaraSampledEnergy
//

#ifndef TiaraSampledEnergy_hh
#define TiaraSampledEnergy_hh TiaraSampledEnergy_hh

#ifdef G4ANALYSIS_USE

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
#endif
