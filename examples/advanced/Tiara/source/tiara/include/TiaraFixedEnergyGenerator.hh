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
// $Id: TiaraFixedEnergyGenerator.hh,v 1.3 2003/06/25 09:12:41 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// ----------------------------------------------------------------------
//
// Class TiaraFixedEnergyGenerator
//

#ifndef TiaraFixedEnergyGenerator_hh
#define TiaraFixedEnergyGenerator_hh TiaraFixedEnergyGenerator_hh

#include "TiaraVSourceEnergyGenerator.hh"

class TiaraFixedEnergyGenerator : public TiaraVSourceEnergyGenerator{
public:
  explicit TiaraFixedEnergyGenerator(G4double energy);
  ~TiaraFixedEnergyGenerator();
  virtual G4double GetEnergy();
  virtual TiaraVSourceEnergyGenerator *Clone() const;
private:
  G4double fEnergy;
};

#endif
