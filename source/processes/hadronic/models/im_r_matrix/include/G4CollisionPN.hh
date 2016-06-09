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
// $Id: G4CollisionPN.hh,v 1.3 2006-06-29 20:34:41 gunter Exp $ //
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//      For information related to this code contact:
//
//      File name:     G4CollisionNN
//
//      Author:        Maria Grazia Pia
// 
//      Creation date: 15 April 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4COLLISIONPN_HH
#define G4COLLISIONPN_HH

#include "globals.hh"
#include "G4CollisionComposite.hh"
#include "G4CollisionVector.hh"
#include "G4VCrossSectionSource.hh"
#include <vector>

class G4KineticTrack;

class G4CollisionPN : public G4CollisionComposite
{

public:

  G4CollisionPN();

  virtual ~G4CollisionPN();

  G4bool operator==(const G4CollisionPN &right) const;
  G4bool operator!=(const G4CollisionPN &right) const;

  virtual G4String GetName() const { return "PN CollisionComposite"; }

private:
  G4CollisionPN(const G4CollisionPN &);
  G4CollisionPN & operator= (const G4CollisionPN &);

protected:

  virtual const G4VCrossSectionSource* GetCrossSectionSource() const 
  { return crossSectionSource; }
  virtual const G4VAngularDistribution* GetAngularDistribution() const 
  { return 0; }

  virtual const std::vector<G4String>& GetListOfColliders(G4int whichOne) const;  

private:  

  G4VCrossSectionSource* crossSectionSource;

  std::vector<G4String> colliders1;
  std::vector<G4String> colliders2;
};

#endif
