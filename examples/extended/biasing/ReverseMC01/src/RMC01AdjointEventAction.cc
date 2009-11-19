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
// $Id: RMC01AdjointEventAction.cc,v 1.1 2009-11-19 22:41:18 ldesorgh Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//////////////////////////////////////////////////////////////
//      Class Name:	RMC01AdjointEventAction
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
//////////////////////////////////////////////////////////////
#include "RMC01AdjointEventAction.hh"
#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryPoint.hh"



#include "G4UnitsTable.hh"
#include "G4SDManager.hh"
#include "Randomize.hh"
#include <iomanip>
#include "RMC01SD.hh"
#include "RMC01AnalysisManager.hh"
RMC01AdjointEventAction::RMC01AdjointEventAction()
{;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
//
RMC01AdjointEventAction::~RMC01AdjointEventAction()
{;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
//
void RMC01AdjointEventAction::BeginOfEventAction(const G4Event* anEvent )
{   G4int i_event=anEvent->GetEventID();
    G4bool print_nb=false;
    if (i_event <100) print_nb =true;
    else if (i_event<500 && (i_event/100)*100 == i_event) print_nb = true;
    else if (i_event<5000 && (i_event/500)*500 == i_event) print_nb = true;
    else if ((i_event/5000)*5000 == i_event) print_nb = true;
    if (print_nb) G4cout<<"nb event "<<i_event<<std::endl;
   
    
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
//
void RMC01AdjointEventAction::EndOfEventAction(const G4Event* )
{ ;
}  
