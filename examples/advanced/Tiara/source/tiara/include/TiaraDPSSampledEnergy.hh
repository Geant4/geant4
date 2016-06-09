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
// $Id: TiaraDPSSampledEnergy.hh,v 1.5 2003/06/25 09:12:37 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// ----------------------------------------------------------------------
//
// Class TiaraDPSSampledEnergy
//

#ifndef TiaraDPSSampledEnergy_hh
#define TiaraDPSSampledEnergy_hh TiaraDPSSampledEnergy_hh

#ifdef G4ANALYSIS_USE

#include "TiaraVSourceEnergyGenerator.hh"
#include <memory>
#include <map>

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
  std::map<int, double> fEnergy_Flux;
  G4double fMaxProb;
  G4double fMinE;
  G4double fMaxE;
};

#endif
#endif
