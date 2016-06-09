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
// $Id: G4CollisionPN.hh,v 1.1 2003/10/07 12:37:27 hpw Exp $ //
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


protected:

  virtual const G4VCrossSectionSource* GetCrossSectionSource() const { return crossSectionSource; }
  virtual const G4VAngularDistribution* GetAngularDistribution() const { return 0; }

  virtual const G4CollisionVector* GetComponents() const { return components; } 

  virtual const std::vector<G4String>& GetListOfColliders(G4int whichOne) const;  

private:  

  G4CollisionVector* components;

  G4VCrossSectionSource* crossSectionSource;

  std::vector<G4String> colliders1;
  std::vector<G4String> colliders2;
};

#endif
