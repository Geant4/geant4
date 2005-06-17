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
#ifndef G4_CASCADE_CHANNEL_HH
#define G4_CASCADE_CHANNEL_HH

#include "globals.hh"
#include <vector>

class G4CascadeChannel {

public:

  G4CascadeChannel();
  virtual ~G4CascadeChannel();
 
protected: 

  std::pair<G4int, G4double> interpolateEnergy(G4double ke) const;
  G4int sampleFlat(std::vector<G4double> sigma) const;

  std::vector<G4int> getQnums(G4int type) const;
private:

  static const double energyScale[31];
};        

#endif
