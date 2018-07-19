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
// $Id$
//
// 
//---------------------------------------------------------------
//
//  G4FastSimulationHelper.hh
//
//  Description:
//      Helper class to insert the G4FastSimulationManagerProcess
//      in a G4ProcessManager object
//
//  History:
//      Nov 2016: M. Verderi - First Implementation.
//
//---------------------------------------------------------------


#ifndef G4FastSimulationHelper_h
#define G4FastSimulationHelper_h 1

#include "globals.hh"

class G4ProcessManager;
class G4FastSmulationManagerProcess;


class G4FastSimulationHelper
{
public:
  // -- Activate fast simulation for particle with pmanager in the mass geometry:
  static void ActivateFastSimulation(G4ProcessManager* pmanager);
  // -- Activate fast simulation for particle with pmanager in the parallel geometry:
  static void ActivateFastSimulation(G4ProcessManager* pmanager, G4String parallelGeometryName );
  
};

#endif
