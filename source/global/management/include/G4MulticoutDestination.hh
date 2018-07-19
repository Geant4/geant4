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
// $Id: G4MulticoutDestination.hh 103582 2017-04-18 17:24:45Z adotti $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4MulticoutDestination.hh
//
// Class description:
//
// Extends G4coutDestination and allows multiple G4coutDestination objects
// to be chained in the same job.
// Note that formatters can be attached to this destination that will apply
// them before passing over the message to the child destination.
//
// Usage:
//   User should create and set-up G4coutDestination objects the usual way
//   these should then be added to an instance of this container class
//   this class should then be added to the streaming. E.g.
//           class MyCout1 : public G4coutDestination {};
//           class MyCout2 : public G4coutDestination {};
//           auto multi = new G4MulticoutDestination();
//           multi->push_back( G4coutDestinationUPtr( new MyCout1 ) );
//           multi->push_back( G4coutDestinationUPtr( new MyCout2 ) );
//           G4coutbuf.SetDestination( multi ); // or G4cerrbuf

//      ---------------- G4MulticoutDestination ----------------
//
// Author: A.Dotti (SLAC), April 2017
// --------------------------------------------------------------------
#ifndef G4MULTICOUTDESTINATION_HH_
#define G4MULTICOUTDESTINATION_HH_

#include <memory>
#include <vector>

#include "G4coutDestination.hh"

using G4coutDestinationUPtr = std::unique_ptr<G4coutDestination>;
using G4coutDestinationVector = std::vector<G4coutDestinationUPtr>;

class G4MulticoutDestination : public G4coutDestination,
                               public G4coutDestinationVector
{
  public:

    G4MulticoutDestination() = default;
    virtual ~G4MulticoutDestination() {}

    // Forward call to contained destination. Note that the message may have
    // been modified by formatters attached to this
    virtual G4int ReceiveG4cout(const G4String& msg) override
    {
      G4bool result = true;
      std::for_each( begin(), end(),
        [&](G4coutDestinationUPtr& e) { result &= (e->ReceiveG4cout_(msg)==0); }
      );
      return ( result ? 0 : -1);
    }

    virtual G4int ReceiveG4cerr(const G4String& msg) override
    {
      G4bool result = true;
      std::for_each( begin(), end(),
        [&](G4coutDestinationUPtr& e) { result &= (e->ReceiveG4cerr_(msg)==0); }
      );
      return ( result ? 0 : -1);
    }
};

#endif /* G4MULTICOUTDESTINATION_HH_ */
