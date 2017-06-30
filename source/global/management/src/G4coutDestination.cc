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
// $Id: G4coutDestination.cc 103661 2017-04-20 14:57:11Z gcosmo $
//
// ----------------------------------------------------------------------
// G4coutDestination
// ----------------------------------------------------------------------

#include "G4coutDestination.hh"

#include <algorithm>

G4coutDestination* G4coutDestination::masterG4coutDestination = 0;

G4coutDestination::~G4coutDestination()
{
}

void G4coutDestination::ResetTransformers()
{
  transformersCout.clear(); transformersCerr.clear();
}

G4int G4coutDestination::ReceiveG4cout(const G4String& msg)
{
  std::cout<<msg<<std::flush;
  return 0;
}

G4int G4coutDestination::ReceiveG4cerr(const G4String& msg)
{
  std::cerr<<msg<<std::flush;
  return 0;
}

G4int G4coutDestination::ReceiveG4cout_(const G4String& msg)
{
  // Avoid copy of string if not necessary
  if( transformersCout.size() > 0 )
  {
    G4String m = msg;
    G4bool result = true;
    for ( const auto& el : transformersCout )
    {
      result &= el(m);
      if ( ! result ) break;
    }
    return ( result ? ReceiveG4cout(m) : 0 );
  }
  else
  {
    return ReceiveG4cout(msg);
  }
}

G4int G4coutDestination::ReceiveG4cerr_(const G4String& msg)
{
  if( transformersCout.size() > 0 )
  {
    G4String m = msg;
    std::for_each( transformersCerr.begin() , transformersCerr.end() ,
                   [&m](const Transformer& t) { t(m); }
                   // Call transforming function on message
                 );
    return ReceiveG4cerr(m);
  }
  else
  {
    return ReceiveG4cerr(msg);
  }
}
