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
// $Id: G4RTRun.hh 73913 2013-09-16 21:51:09Z asaim $
//
//

// class description:
//

//////////////////////
//G4RTRun
/////////////////////

#ifndef G4RTRun_h
#define G4RTRun_h 1

#include "globals.hh"
#include "G4Run.hh"
#include "G4THitsMap.hh"
#include "G4Colour.hh"

class G4Event;
class G4TheMTRayTracer;
class G4VisAttributes;
class G4RayTrajectoryPoint;

class G4RTRun : public G4Run
{
  public:
    G4RTRun();
    virtual ~G4RTRun();

    virtual void RecordEvent(const G4Event*);
    virtual void Merge(const G4Run*);

  private:
    G4THitsMap<G4Colour>* colorMap;

  public:
    G4THitsMap<G4Colour>* GetMap() const { return colorMap; }

  private:
    G4Colour backgroundColour;
    G4ThreeVector lightDirection;
    G4double attenuationLength;

  private:
    G4Colour GetSurfaceColour(G4RayTrajectoryPoint*);
    G4Colour GetMixedColour(G4Colour,G4Colour,G4double);
    G4Colour Attenuate(G4RayTrajectoryPoint*,G4Colour);
    G4bool ValidColour(const G4VisAttributes*);
};

#endif
