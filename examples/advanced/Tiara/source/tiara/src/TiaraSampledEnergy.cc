#include "TiaraSampledEnergy.hh"
#include "Randomize.hh"
#include "TiaraFileAcess.hh"


#include "AIDA/AIDA.h"



TiaraSampledEnergy::TiaraSampledEnergy(const G4String &eng,
				       G4double minEnergyCut,
				       const G4String &sourceTree,
				       const G4String &nameExt)
  :
  TiaraVSourceEnergyGenerator(),
  fTree((checkFileIsReadable(sourceTree, "TiaraSampledEnergy::TiaraSampledEnergy"),
	 AIDA_createAnalysisFactory()
	 ->createTreeFactory()
	 ->create(sourceTree, "xml",true,false))),
  fSampleHisto(dynamic_cast<AIDA::IHistogram1D *>(fTree->find(G4String("h" + eng + "int" + nameExt)))),
  fScale(fSampleHisto->axis().bins()),
  fMinEnergyCut(minEnergyCut)
{}

TiaraSampledEnergy::TiaraSampledEnergy(const TiaraSampledEnergy& rhs) :
  TiaraVSourceEnergyGenerator()
{
  *this = rhs;
}

TiaraSampledEnergy& TiaraSampledEnergy::
operator=(const TiaraSampledEnergy& rhs) {
  if (this!=&rhs) {
    TiaraSampledEnergy &nonConstRhs = static_cast<TiaraSampledEnergy &>(rhs);
    fTree = nonConstRhs.fTree;
    fSampleHisto = nonConstRhs.fSampleHisto;
    fScale = nonConstRhs.fScale;
  }
  return *this;
}

TiaraSampledEnergy::~TiaraSampledEnergy(){
}


TiaraVSourceEnergyGenerator *TiaraSampledEnergy::Clone() const {
  return new TiaraSampledEnergy(*this);
}

G4double TiaraSampledEnergy::GetEnergy() {

  G4double e(0.0);
  while (!(e>fMinEnergyCut)) { 
    G4double r(G4UniformRand() * fScale);
    G4int bin(static_cast<int>(r));
    e = fSampleHisto->binHeight(bin) * MeV;
  }
  return e;
}

