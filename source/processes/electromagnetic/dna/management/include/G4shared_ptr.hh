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
// $Id: $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// -------------------------------------------------------------------

#ifndef G4SHARED_PTR_HH_
#define G4SHARED_PTR_HH_

#include "CLHEP/Utility/memory.h"

namespace G4
{
 using CLHEP::shared_ptr;
 using CLHEP::weak_ptr;
 using CLHEP::enable_shared_from_this;
 using CLHEP::static_pointer_cast;
 using CLHEP::const_pointer_cast;
 using CLHEP::dynamic_pointer_cast;
 //using CLHEP::polymorphic_cast_tag;
}

#define G4shared_ptr G4::shared_ptr

#endif /* G4SHARED_PTR_HH_ */
