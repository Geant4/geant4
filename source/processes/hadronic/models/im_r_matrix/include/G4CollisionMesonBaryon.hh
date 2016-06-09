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
// $Id: G4CollisionMesonBaryon.hh,v 1.3 2006-06-29 20:32:29 gunter Exp $ //
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

#ifndef G4COLLISIONMESONBARYON_HH
#define G4COLLISIONMESONBARYON_HH

#include "globals.hh"
#include "G4CollisionComposite.hh"
#include "G4CollisionVector.hh"
#include "G4VCrossSectionSource.hh"
#include <vector>

class G4KineticTrack;

class G4CollisionMesonBaryon : public G4CollisionComposite
{

public:

  G4CollisionMesonBaryon();

  virtual ~G4CollisionMesonBaryon() {};

  virtual G4String GetName() const { return "Meson Baryon CollisionComposite"; }

private:
  G4CollisionMesonBaryon(const G4CollisionMesonBaryon &);
  G4CollisionMesonBaryon & operator= (const G4CollisionMesonBaryon &);

protected:

  std::vector<G4String> result;
  virtual const std::vector<G4String>& GetListOfColliders(G4int ) const
  {
    throw G4HadronicException(__FILE__, __LINE__, "Tried to call G4CollisionNNToDeltaDelta::GetListOfColliders. Please find out why!");
    return result;
  } 
};

#endif
