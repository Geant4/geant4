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
// $Id: RemSimPhysicsList.hh,v 1.1 2004-01-30 12:18:25 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
// ExRemSimPhysicsList
//  Construct/define particles and physics processes
//
//  Particle defined in ExampleRemSim
//    geantino
//
//  Process defined in ExampleRemSim
//    transportation
//

#ifndef RemSimPhysicsList_h
#define RemSimPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class RemSimPhysicsList: public G4VUserPhysicsList
{
public:
  RemSimPhysicsList();
  ~RemSimPhysicsList();

protected:
  // Construct particle and physics process
  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();

private:
  void ConstructEM();
};
#endif







