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
// $Id: RE01RegionInformation.hh,v 1.2 2006-06-29 17:43:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef RE01RegionInformation_H
#define RE01RegionInformation_H 1

#include "globals.hh"
#include "G4VUserRegionInformation.hh"

class RE01RegionInformation : public G4VUserRegionInformation
{
  public:
    RE01RegionInformation(); 
    ~RE01RegionInformation();
    void Print() const;

  private:
    G4bool isWorld;
    G4bool isTracker;
    G4bool isCalorimeter;

  public:
    inline void SetWorld(G4bool v=true) {isWorld = v;}
    inline void SetTracker(G4bool v=true) {isTracker = v;}
    inline void SetCalorimeter(G4bool v=true) {isCalorimeter = v;}
    inline G4bool IsWorld() const {return isWorld;}
    inline G4bool IsTracker() const {return isTracker;}
    inline G4bool IsCalorimeter() const {return isCalorimeter;}
};

#endif

