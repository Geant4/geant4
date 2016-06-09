//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: TiaraSampledEnergy.cc,v 1.4 2003/06/25 09:13:11 gunter Exp $
// GEANT4 tag $Name: geant4-06-00 $
//

#include "TiaraSampledEnergy.hh"

#ifdef G4ANALYSIS_USE

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

#endif
