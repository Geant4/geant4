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
// by J.P Wellisch, Sun Sep 15 2002.

#include "g4std/vector"
#include "g4std/algorithm"
#include "globals.hh"

class G4PiData : public G4std::vector<G4std::pair<G4double, G4std::pair<G4double, G4double > > > 
{
  public:
    G4PiData(const G4double * aTotal, const G4double * aInelastic, const G4double * anEnergy, G4int nPoints);
    struct Delete{void operator()(G4PiData * aP){delete aP;} };
    G4bool AppliesTo(G4double kineticEnergy);
    G4double ReactionXSection(G4double kineticEnergy);
    G4double ElasticXSection(G4double kineticEnergy);
    
  private:
    
};

#endif
