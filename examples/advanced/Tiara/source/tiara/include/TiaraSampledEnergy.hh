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
// $Id: TiaraSampledEnergy.hh,v 1.4 2003/06/25 09:12:50 gunter Exp $
// GEANT4 tag $Name: geant4-05-02-patch-01 $
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
