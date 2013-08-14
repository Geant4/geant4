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
// $Id: G4EvaporationInuclCollider.hh 71942 2013-06-28 19:08:11Z mkelsey $
//
// 20100315  M. Kelsey -- Remove "using" directive and unnecessary #includes.
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class, make other colliders
//		simple data members
// 20100714  M. Kelsey -- Switch to new G4CascadeColliderBase class
// 20110728  M. Kelsey -- Fix Coverity #23843, add destructor.
// 20130620  Address Coverity complaint about missing copy actions
// 20130622  Inherit from G4CascadeDeexciteBase, move to deExcite() interface
//		with G4Fragment

#ifndef G4EVAPORATIONINUCL_COLLIDER_HH
#define G4EVAPORATIONINUCL_COLLIDER_HH
 
#include "G4CascadeDeexciteBase.hh"

class G4InuclParticle;
class G4CollisionOutput;
class G4EquilibriumEvaporator;
class G4BigBanger;

class G4EvaporationInuclCollider : public G4CascadeDeexciteBase {
public:
  G4EvaporationInuclCollider();
  ~G4EvaporationInuclCollider();

  virtual void deExcite(const G4Fragment& target, G4CollisionOutput& output);
  
private: 
  G4EquilibriumEvaporator* theEquilibriumEvaporator;

private:
  // Copying of modules is forbidden
  G4EvaporationInuclCollider(const G4EvaporationInuclCollider&);
  G4EvaporationInuclCollider& operator=(const G4EvaporationInuclCollider&);
};        

#endif /* G4EVAPORATIONINUCL_COLLIDER_HH */


