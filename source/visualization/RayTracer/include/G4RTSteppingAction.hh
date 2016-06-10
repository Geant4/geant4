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
// $Id: G4RTSteppingAction.hh 74050 2013-09-20 09:38:19Z gcosmo $
//
//

// class description:
//
//  This is a concrete class of G4UserSteppingAction. This class is used
// by G4RayTracer for managing a ray tracked through volumes. An object
// of this class is constructed by G4RayTracer and set to G4SteppingManager
// with replacement of user defined stepping action during the period of
// ray tracing.
//

//////////////////////
//G4RTSteppingAction
/////////////////////


#ifndef G4RTSteppingAction_h
#define G4RTSteppingAction_h 1


#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4RTSteppingAction : public G4UserSteppingAction
{
  public:
    G4RTSteppingAction();
    virtual ~G4RTSteppingAction(){;}

    virtual void UserSteppingAction(const G4Step*);

  private:
    static G4bool ignoreTransparency;

  public:
    static void SetIgnoreTransparency(G4bool val);
    static G4bool GetIgnoreTransparency();
};

#endif
