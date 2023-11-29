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
#ifndef G4VDNAHitModel_hh
#define G4VDNAHitModel_hh
#include "globals.hh"
#include <vector>
#include <variant>
class G4DNAComponentNode;
class G4VPhysicalVolume;
class G4Track;

class G4VDNAHitModel
{
  using DNANode =
    std::variant<const G4DNAComponentNode*, /*for dnadamage chain*/
                 const G4VPhysicalVolume* /*for molecularDNA chain*/>;

 public:
  explicit G4VDNAHitModel(const G4String& name);
  virtual ~G4VDNAHitModel() = default;
  //  delete assignment operator
  G4VDNAHitModel& operator=(const G4VDNAHitModel& right) = delete;
  G4VDNAHitModel(const G4VDNAHitModel&)                  = delete;
  virtual G4double CalculateReactionTime(const G4Track& trackA, DNANode&) = 0;
  virtual G4bool DoReaction(const G4Track& track, const G4double&,
                            const DNANode&)                               = 0;

 private:
  const G4String fName;
};
#endif
