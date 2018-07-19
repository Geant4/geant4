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
// $Id: G4LockcoutDestination.hh 103582 2017-04-18 17:24:45Z adotti $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file 
//
// Class Description:
//
// Implements a output destination to std::cout / std::cerr with a
// mutex lock access to the shared resource.

//      ---------------- G4LockcoutDestination ----------------
//
// Author: A.Dotti (SLAC), April 2017
// --------------------------------------------------------------------
#ifndef G4LOCKCOUTDESTINATION_HH_
#define G4LOCKCOUTDESTINATION_HH_

#include "G4coutDestination.hh"

class G4LockcoutDestination : public G4coutDestination
{
  public:

    G4LockcoutDestination() = default;
    virtual ~G4LockcoutDestination();
    virtual G4int ReceiveG4cout(const G4String& msg) override;
    virtual G4int ReceiveG4cerr(const G4String& msg) override;
};

#endif /* G4LOCKCOUTDESTINATION_HH_ */
