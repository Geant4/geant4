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
//
// $Id: G4StringModel.hh,v 1.1 2003/10/07 11:26:00 hpw Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
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
