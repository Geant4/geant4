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
// John Allison, 30th August 2020

#ifndef G4QT3DQENTITY_HH
#define G4QT3DQENTITY_HH

#include <Qt3DCore>

#include "G4PhysicalVolumeModel.hh"

class G4Qt3DQEntity: public Qt3DCore::QEntity {
public:
  G4Qt3DQEntity(Qt3DCore::QNode *parent = nullptr)
  : Qt3DCore::QEntity(parent) {}
  const G4PhysicalVolumeModel::G4PhysicalVolumeNodeID& GetPVNodeID () const
  {return fPVNodeID;}
  void SetPVNodeID(const G4PhysicalVolumeModel::G4PhysicalVolumeNodeID& id)
  {fPVNodeID = id;}
private:
  G4PhysicalVolumeModel::G4PhysicalVolumeNodeID fPVNodeID;
};

#endif
