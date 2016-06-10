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
// $Id: RE05PrimaryGeneratorAction.hh 69764 2013-05-14 09:59:36Z gcosmo $
//
/// \file RE05/include/RE05PrimaryGeneratorAction.hh
/// \brief Definition of the RE05PrimaryGeneratorAction class
//

#ifndef RE05PrimaryGeneratorAction_h
#define RE05PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4VPrimaryGenerator;
class G4Event;
class RE05PrimaryGeneratorMessenger;

class RE05PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    RE05PrimaryGeneratorAction();
    virtual ~RE05PrimaryGeneratorAction();

  public:
    virtual void GeneratePrimaries(G4Event* anEvent);

  private:
    static G4VPrimaryGenerator* HEPEvt;
    G4VPrimaryGenerator* particleGun;
    RE05PrimaryGeneratorMessenger* messenger;
    G4bool useHEPEvt;

  public:
    inline void SetHEPEvtGenerator(G4bool f)
    { useHEPEvt = f; }
    inline G4bool GetHEPEvtGenerator()
    { return useHEPEvt; }
};

#endif


