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
//
//
#ifndef G4StringModel_h
#define G4StringModel_h 1

#include "G4VHighEnergyGenerator.hh"
#include "G4KineticTrackVector.hh"

class G4V3DNucleus;
class G4VStringFragmentation;


class G4StringModel : public G4VHighEnergyGenerator 
{
  public:
    G4StringModel();
    ~G4StringModel() override;

    G4StringModel(const G4StringModel &right) = delete;
    const G4StringModel & operator=(const G4StringModel &right) = delete;
    G4bool operator==(const G4StringModel &right) const = delete;
    G4bool operator!=(const G4StringModel &right) const = delete;

    void Set3DNucleus(G4V3DNucleus *const  value);
    void SetStringFragmentationModel(G4VStringFragmentation *const  value);

  private: 
    G4V3DNucleus *the3DNucleus;
    G4VStringFragmentation *theStringFragmentationModel;
};

inline void G4StringModel::Set3DNucleus(G4V3DNucleus *const  value)
{
  the3DNucleus = value;
}

inline void G4StringModel::SetStringFragmentationModel(G4VStringFragmentation *const  value)
{
  theStringFragmentationModel = value;
}

#endif

