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

#pragma once
#include "G4VUserParallelWorld.hh"
#include "globals.hh"

<<<<<<< HEAD:source/processes/hadronic/models/parton_string/qgsm/include/G4QGSMParameters.hh
class G4QGSMParameters 
    {
public:
      G4QGSMParameters();
      ~G4QGSMParameters();
      
private:
      G4QGSMParameters(const G4QGSMParameters &right);
      int operator==(const G4QGSMParameters &right) const;
      int operator!=(const G4QGSMParameters &right) const;
     };

#endif     
=======
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ParallelWorld: public G4VUserParallelWorld
{
public:
    ParallelWorld(G4String);
    ~ParallelWorld() override;
    void Construct() override;

private:
    void SetupChemistryGeometry();
};
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c:examples/extended/medical/dna/dnadamage1/include/ParallelWorld.hh
