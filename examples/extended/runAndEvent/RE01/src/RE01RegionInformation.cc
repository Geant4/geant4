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
// $Id: RE01RegionInformation.cc,v 1.1 2004/11/26 07:37:42 asaim Exp $
// GEANT4 tag $Name: geant4-07-01 $
//

#include "RE01RegionInformation.hh"
#include "G4ios.hh"

RE01RegionInformation::RE01RegionInformation()
:isWorld(false),isTracker(false),isCalorimeter(false)
{;}

RE01RegionInformation::~RE01RegionInformation()
{;}

void RE01RegionInformation::Print() const
{
 G4cout << "I'm ";
 if(isWorld) { G4cout << "World."; }
 else if(isTracker) { G4cout << "Tracker."; }
 else if(isCalorimeter) { G4cout << "Calorimeter."; }
 else { G4cout << "unknown."; }
 G4cout << G4endl;
}

