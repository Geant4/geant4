//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: TiaraSampledEnergy.cc,v 1.5 2006/06/29 15:45:29 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
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
