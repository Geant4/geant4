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
// G4coutDestination class implementation
//
// Authors: H.Yoshida, M.Nagamatu, 1997
// --------------------------------------------------------------------

#include "G4coutDestination.hh"

#include <algorithm>

G4coutDestination* G4coutDestination::masterG4coutDestination = nullptr;

// --------------------------------------------------------------------
void G4coutDestination::ResetTransformers()
{
  transformersCout.clear();
  transformersCerr.clear();
}

// --------------------------------------------------------------------
G4int G4coutDestination::ReceiveG4cout(const G4String& msg)
{
  std::cout << msg << std::flush;
  return 0;
}

// --------------------------------------------------------------------
G4int G4coutDestination::ReceiveG4cerr(const G4String& msg)
{
  std::cerr << msg << std::flush;
  return 0;
}

// --------------------------------------------------------------------
G4int G4coutDestination::ReceiveG4cout_(const G4String& msg)
{
  // Avoid copy of string if not necessary
  if(!transformersCout.empty())
  {
    G4String m    = msg;
    G4bool result = true;
    for(const auto& el : transformersCout)
    {
      result &= el(m);
      if(!result)
      {
        break;
      }
    }
    return (result ? ReceiveG4cout(m) : 0);
  }
  
  return ReceiveG4cout(msg);
}

// --------------------------------------------------------------------
G4int G4coutDestination::ReceiveG4cerr_(const G4String& msg)
{
  if(!transformersCout.empty())
  {
    G4String m = msg;
    std::for_each(transformersCerr.begin(), transformersCerr.end(),
                  [&m](const Transformer& t) { t(m); }
                  // Call transforming function on message
    );
    return ReceiveG4cerr(m);
  }
  
  return ReceiveG4cerr(msg);
}
