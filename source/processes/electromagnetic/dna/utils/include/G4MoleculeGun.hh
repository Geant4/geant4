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
// Author: Mathieu Karamitros

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef MOLECULEGUN_HH_
#define MOLECULEGUN_HH_

#include "G4ITGun.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4memory.hh"

#include <vector>
#include <map>

class G4Track;
class G4MoleculeGunMessenger;
class G4MoleculeShootMessenger;
class G4MoleculeGun;

using G4ContinuousMedium = G4int;

//------------------------------------------------------------------------------

/*
 * Define a specific species shoot
 * Multiple shoots maybe be defined
 * TODO: make the shoots fully time dependent with
 *       user access in the new framework
 */

class G4MoleculeShoot : public G4enable_shared_from_this<G4MoleculeShoot>
{
public:

  G4MoleculeShoot();
  virtual ~G4MoleculeShoot();
  virtual void Shoot(G4MoleculeGun*) = 0;

  template<typename TYPE> G4shared_ptr<G4MoleculeShoot> ChangeType();

  G4String fMoleculeName;
  G4ThreeVector fPosition;
  G4double fTime;
  G4int fNumber;
  G4ThreeVector* fBoxSize;

  static void RandomPosInBox(const G4ThreeVector& boxSize,
                             G4ThreeVector& output);
};

//------------------------------------------------------------------------------

/*
 * Define a shoot type =
 *   track (used by the Smoluchowski code)
 *   or continuous medium (used by the gillespie code)
 */

template<typename TYPE>
class TG4MoleculeShoot : public G4MoleculeShoot
{
public:
  TG4MoleculeShoot() : G4MoleculeShoot(){;}
  virtual ~TG4MoleculeShoot(){;}
  void Shoot(G4MoleculeGun*){}

protected:
  void ShootAtRandomPosition(G4MoleculeGun*){}
  void ShootAtFixedPosition(G4MoleculeGun*){}
};

template<typename TYPE>
G4shared_ptr<G4MoleculeShoot> G4MoleculeShoot::ChangeType()
{
  G4shared_ptr<G4MoleculeShoot> output(new TG4MoleculeShoot<TYPE>);
  output->fMoleculeName = fMoleculeName;
  output->fPosition = fPosition;
  output->fTime = fTime;
  output->fNumber = fNumber;
  output->fBoxSize = fBoxSize;
  return output;
}


//------------------------------------------------------------------------------

class G4MoleculeGun : public G4ITGun
{
public:
  G4MoleculeGun();
  virtual ~G4MoleculeGun();

  virtual void DefineTracks();

  /*
   * Create a single molecule
   * @param moleculeName name of the molecule such as recorded in molecule table
   * @param position position where the molecule should pop up
   * @param time time at which the molecule should pop up
   */
  void AddMolecule(const G4String& moleculeName,
                   const G4ThreeVector& position,
                   G4double time = 0);

  /*
   * Create N molecules at a single point
   * @param n number of molecules to create
   * @param moleculeName name of the molecules such as recorded in molecule table
   * @param position position where the molecules should pop up
   * @param time time at which the molecules should pop up
   */
  void AddNMolecules(std::size_t n,
                     const G4String& moleculeName,
                     const G4ThreeVector& position,
                     G4double time = 0);

  /*
   * Create N molecules in a box
   * @param n number of molecules to create
   * @param moleculeName name of the molecules such as recorded in molecule table
   * @param boxCenter center of the box
   * @param boxExtension size of the box
   * @param time time at which the molecules should pop up
   */
  void AddMoleculesRandomPositionInBox(std::size_t n,
                                       const G4String& moleculeName,
                                       const G4ThreeVector& boxCenter,
                                       const G4ThreeVector& boxExtension,
                                       G4double time = 0);

  /*
   * Create N molecules as component of the continuous medium in a box
   * @param n number of molecules to create
   * @param moleculeName name of the molecules such as recorded in molecule table
   * @param boxCenter center of the box
   * @param boxExtension size of the box
   * @param time time at which the molecules should pop up
   */
//  void AddMoleculeInCMRepresentation(std::size_t n,
//                                     const G4String& moleculeName,
//                                     const G4ThreeVector& boxCenter,
//                                     const G4ThreeVector& boxExtension,
//                                     G4double time = 0);

  void AddMoleculeInCMRepresentation(std::size_t n,
                                     const G4String& moleculeName,
                                     G4double time = 0);

  const std::vector<G4shared_ptr<G4MoleculeShoot> >&
      GetMoleculeShoot() {
    return fShoots;
  }

  typedef std::map<G4String, G4int> NameNumber;
  void GetNameAndNumber(NameNumber&);


  void AddMoleculeShoot(G4shared_ptr<G4MoleculeShoot>);

protected:
  void BuildAndPushTrack(const G4String& name,
                         const G4ThreeVector& position,
                         G4double time = 0);
  G4MoleculeGunMessenger* fpMessenger;

  std::vector<G4shared_ptr<G4MoleculeShoot> > fShoots;
  friend class G4MoleculeShoot;
  template<class T> friend class TG4MoleculeShoot;
};

#endif /* MOLECULEGUN_HH_ */
