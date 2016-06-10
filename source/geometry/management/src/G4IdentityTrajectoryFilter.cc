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
// $Id: G4IdentityTrajectoryFilter.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
// ------------------------------------------------------------------------

#include "G4IdentityTrajectoryFilter.hh"

G4IdentityTrajectoryFilter::G4IdentityTrajectoryFilter()
{
}

G4IdentityTrajectoryFilter::~G4IdentityTrajectoryFilter()
{
}

void
G4IdentityTrajectoryFilter::TakeIntermediatePoint( G4ThreeVector newPoint )
{
  // Just store every single point, initially. (jacek 30/10/2002)
  // Implement more sophisticated filters later.  Copy by value into
  // the vector; the vector will never be copied by value itself. In
  // the final version, will probably want to create the intermediate
  // points at this stage.
  //
  fpFilteredPoints->push_back( newPoint );
}
