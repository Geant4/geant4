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
// by J.P Wellisch, Sun Sep 15 2002.

#ifndef G4PiData_h
#define G4PiData_h

#include <vector>
#include <algorithm>
#include "globals.hh"

class G4PiData : public std::vector<std::pair<G4double, std::pair<G4double, G4double > > > 
{
  public:
    G4PiData(const G4double * aTotal, const G4double * aInelastic, 
             const G4double * anEnergy, G4int nPoints);

    struct Delete{void operator()(G4PiData * aP){delete aP;} };

    G4bool AppliesTo(G4double kineticEnergy);

    G4double ReactionXSection(G4double kineticEnergy);
    G4double ElasticXSection(G4double kineticEnergy);
    G4double TotalXSection(G4double kineticEnergy);
    
  private:
    
};

#endif
