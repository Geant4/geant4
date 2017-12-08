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
// $Id: G4PhysListFactoryMessenger.hh 66704 2013-01-10 18:20:17Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4PhysListFactoryMessenger
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 09.11.2005 V.Ivanchenko edit to provide a standard
// 19.06.2006 V.Ivanchenko add mu-nuclear process
//
//----------------------------------------------------------------------------
//

#ifndef G4PhysListFactoryMessenger_h
#define G4PhysListFactoryMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"

class G4VModularPhysicsList;

class G4PhysListFactoryMessenger: public G4UImessenger
{
public:
  G4PhysListFactoryMessenger(G4VModularPhysicsList* pl);
  virtual ~G4PhysListFactoryMessenger();

  void SetNewValue(G4UIcommand* aComm, G4String aS);

private:
  G4VModularPhysicsList*   thePhysList;
  G4UIcommand*             theRadDecay;
  G4UIcommand*             theOptical;
  G4UIdirectory*           theDir;
};

#endif
