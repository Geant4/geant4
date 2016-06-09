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
// $Id: G4CollisionMesonBaryon.hh,v 1.7 2003/06/16 17:07:54 gunter Exp $ //
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


protected:

  std::vector<G4String> result;
  virtual const std::vector<G4String>& GetListOfColliders(G4int ) const
  {
    G4Exception("Tried to call G4CollisionNNToDeltaDelta::GetListOfColliders. Please find out why!");
    return result;
  } 
};

#endif
