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
// $Id: G4VITStepModel.cc 81769 2014-06-05 08:30:21Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4VITStepModel.hh"

G4VITStepModel::G4VITStepModel(const G4String& aName)
{
    //ctor
    fpTimeStepper        = 0;
    fpReactionProcess    = 0;
    fpReactionTable      = 0;

    fType1              = -1;
    fType2              = -1;
    fName               = aName;
}

G4VITStepModel::~G4VITStepModel()
{
    //dtor
    if(fpTimeStepper)        delete fpTimeStepper;
    if(fpReactionProcess)    delete fpReactionProcess;
    //if(fReactionTable)      delete fReactionTable;
    // Let the concrete class delete the reactionTable
}

G4VITStepModel::G4VITStepModel(const G4VITStepModel& right)
{
    //copy ctor
    fName               = right.fName;
    fType1              = right.fType1;
    fType2              = right.fType2;
    fpReactionTable     = 0;
    fpTimeStepper     = right.fpTimeStepper->Clone();
    fpReactionProcess = right.fpReactionProcess->Clone();
}

// should not be used
G4VITStepModel& G4VITStepModel::operator=(const G4VITStepModel& right)
{
    if (this == &right) return *this; // handle self assignment

    fName               = right.fName;
    fType1              = right.fType1;
    fType2              = right.fType2;
    fpReactionTable      = 0;
    if(fpTimeStepper) delete fpTimeStepper;
    fpTimeStepper        = right.fpTimeStepper->Clone();
    if(fpReactionProcess) delete fpReactionProcess;
    fpReactionProcess    = right.fpReactionProcess->Clone();

    //assignment operator
    return *this;
}

void G4VITStepModel::IsApplicable(G4ITType& type1, G4ITType& type2)
{
    type1 = fType1;
    type2 = fType2;
    PrintInfo();
}

void G4VITStepModel::Initialize()
{
    fpReactionProcess->SetReactionTable(fpReactionTable);
    fpTimeStepper->SetReactionTable(fpReactionTable);
    fpTimeStepper->Initialize();
    fpReactionProcess->Initialize();
}
