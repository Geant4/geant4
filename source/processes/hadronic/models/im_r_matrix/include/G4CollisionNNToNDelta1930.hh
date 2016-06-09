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
#ifndef G4CollisionNNToNDelta1930_h
#define G4CollisionNNToNDelta1930_h

#include "globals.hh"
#include "G4GeneralNNCollision.hh"
#include "G4VCrossSectionSource.hh"
#include "G4VAngularDistribution.hh"
#include "G4KineticTrackVector.hh"
#include <vector>

class G4CollisionNNToNDelta1930 : public G4GeneralNNCollision
{

public:

  G4CollisionNNToNDelta1930();

  virtual ~G4CollisionNNToNDelta1930();


  virtual G4String GetName() const { return "NN -> N Delta(1930) Collision"; }
  virtual const std::vector<G4String>& GetListOfColliders(G4int ) const
  {
    throw G4HadronicException(__FILE__, __LINE__, "Tried to call G4CollisionNNToNDelta1930::GetListOfColliders. Please find out why!");
    std::vector<G4String> * aList = new std::vector<G4String>;
    return *aList;
  } 
  
private:
  G4CollisionNNToNDelta1930(const G4CollisionNNToNDelta1930 &);
  G4CollisionNNToNDelta1930 & operator= (const G4CollisionNNToNDelta1930 &);

protected:
  
  virtual const G4CollisionVector* GetComponents() const { return components; } 

private:  

  G4CollisionVector* components;

};

#endif
