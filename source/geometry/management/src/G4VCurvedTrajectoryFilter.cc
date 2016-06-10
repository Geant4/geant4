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
// $Id: G4VCurvedTrajectoryFilter.cc 66356 2012-12-18 09:02:32Z gcosmo $
// --------------------------------------------------------------------

#include "G4VCurvedTrajectoryFilter.hh"

G4VCurvedTrajectoryFilter::G4VCurvedTrajectoryFilter()
  : fpFilteredPoints(0)
{
}

G4VCurvedTrajectoryFilter::~G4VCurvedTrajectoryFilter()
{
}

std::vector<G4ThreeVector>* 
G4VCurvedTrajectoryFilter::GimmeThePointsAndForgetThem()
{
  std::vector<G4ThreeVector>* tmp = fpFilteredPoints;
  // ParticleChangeForTransport invokes this method (via
  // PropagatorInField) at every Step, even if the step did not
  // involve PropagatorInField. Must, therefore, ensure that points
  // submitted by previous invocations of PIF are not
  // copied. Therefore the points must be cleared. (Note that the
  // responsibility for deleting the vector lies with the
  // SmoothTrajctoryPoint, which is the vector's final destination.)
  // (jacek 08/11/2002)
  fpFilteredPoints = 0;
  return tmp;
}

void
G4VCurvedTrajectoryFilter::CreateNewTrajectorySegment( )
{
  if (fpFilteredPoints)
  {
    // GimmePoints has not been called (it would have set the
    // pointer to NULL), therefore nobody has taken charge of the
    // points and they will never be deleted!
    G4cout << "!!!!!!!! Filter: auxiliary points are being memory leaked !!!!!"
           << G4endl;
  }
  fpFilteredPoints = new std::vector<G4ThreeVector>;
}    
