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
// $Id: MyMultiNavigator.hh,v 1.1 2006-11-11 01:35:38 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Class to help test Multi-navigator geometry 
//
// History:
// 7.11.06 J. Apostolakis  Adapted from MyNavigator in ApproxEqual

#ifndef MYMULTINAV_HH
#define MYMULTINAV_HH

#include "globals.hh"
#include "geomdefs.hh"

#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"
#include "G4MultiNavigator.hh"

// Derived multi-navigator class

class MyMultiNavigator : public G4MultiNavigator
{
  public:
  
    MyMultiNavigator() : G4MultiNavigator() { 
      G4ThreeVector zero(0.0, 0.0, 0.0); 
      // PrepareNewTrack( zero, zero );
    }
    ~MyMultiNavigator(){}

    G4ThreeVector CurrentLocalCoordinate() const
    {
      G4Navigator* fpMassNavigator;    
      fpMassNavigator= G4MultiNavigator::GetNavigator(0); 
      return fpMassNavigator->GetCurrentLocalCoordinate();
    }
    G4ThreeVector GetNetTranslation() const
    {
      G4Navigator* fpMassNavigator;    
      // G4ThreeVector nouse= G4Navigator::NetTranslation();
      fpMassNavigator= GetNavigator(0); 
      return fpMassNavigator->NetTranslation(); 
    }
    G4RotationMatrix GetNetRotation() const
    {
      G4Navigator* fpMassNavigator;    
      fpMassNavigator= GetNavigator(0); 
      return fpMassNavigator->NetRotation();
    }

  private:
    // G4Navigator* fpMassNavigator;    
};
#endif
