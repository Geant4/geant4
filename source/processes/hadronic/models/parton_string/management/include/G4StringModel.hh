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
// $Id: G4StringModel.hh 100828 2016-11-02 15:25:59Z gcosmo $
//
#ifndef G4StringModel_h
#define G4StringModel_h 1

#include "G4VHighEnergyGenerator.hh"
#include "G4EventGenerator.hh"
#include "G4KineticTrackVector.hh"
class G4V3DNucleus;
class G4VStringFragmentation;


class G4StringModel : public G4VHighEnergyGenerator 
{
  public:
    G4StringModel();
    ~G4StringModel();

  private:
    G4StringModel(const G4StringModel &right);
    const G4StringModel & operator=(const G4StringModel &right);
    int operator==(const G4StringModel &right) const;
    int operator!=(const G4StringModel &right) const;

  public:
    void Set3DNucleus(G4V3DNucleus *const  value);
    void SetStringFragmentationModel(G4VStringFragmentation *const  value);
    void SetGenerator(G4EventGenerator *const  value);

  private:
    const G4V3DNucleus * Get3DNucleus() const;
    const G4VStringFragmentation * GetStringFragmentationModel() const;
    const G4EventGenerator * GetGenerator() const;

  private: 
    G4V3DNucleus *the3DNucleus;
    G4VStringFragmentation *theStringFragmentationModel;
    G4EventGenerator *theGenerator;
};

inline const G4V3DNucleus * G4StringModel::Get3DNucleus() const
{
  return the3DNucleus;
}

inline void G4StringModel::Set3DNucleus(G4V3DNucleus *const  value)
{
  the3DNucleus = value;
}

inline const G4VStringFragmentation * G4StringModel::GetStringFragmentationModel() const
{
  return theStringFragmentationModel;
}

inline void G4StringModel::SetStringFragmentationModel(G4VStringFragmentation *const  value)
{
  theStringFragmentationModel = value;
}

inline const G4EventGenerator * G4StringModel::GetGenerator() const
{
  return theGenerator;
}

inline void G4StringModel::SetGenerator(G4EventGenerator *const  value)
{
  theGenerator = value;
}

#endif

