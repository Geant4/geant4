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
// $Id: TiaraSampledEnergy.hh,v 1.5 2006/06/29 15:44:12 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
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
