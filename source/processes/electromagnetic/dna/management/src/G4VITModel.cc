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
// $Id: G4VITModel.cc 60427 2012-07-11 16:34:35Z matkara $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4VITModel.hh"

G4VITModel::G4VITModel(const G4String& aName)
{
    //ctor
    fpTimeStepper        = 0;
    fpReactionProcess    = 0;
    fpReactionTable      = 0;

    fType1              = -1;
    fType2              = -1;
    fName               = aName;
}

G4VITModel::~G4VITModel()
{
    //dtor
    if(fpTimeStepper)        delete fpTimeStepper;
    if(fpReactionProcess)    delete fpReactionProcess;
    //if(fReactionTable)      delete fReactionTable;
    // Let the concrete class delete the reactionTable
}

G4VITModel::G4VITModel(const G4VITModel& other)
{
    //copy ctor
    fName               = other.fName;
    fType1              = other.fType1;
    fType2              = other.fType2;
    fpReactionTable      = 0;
    fpTimeStepper     = other.fpTimeStepper->Clone();
    fpReactionProcess = other.fpReactionProcess->Clone();
}

// should not be used
G4VITModel& G4VITModel::operator=(const G4VITModel& rhs)
{
    if (this == &rhs) return *this; // handle self assignment

    fName               = rhs.fName;
    fType1              = rhs.fType1;
    fType2              = rhs.fType2;
    fpReactionTable      = 0;
    if(fpTimeStepper) delete fpTimeStepper;
    fpTimeStepper        = rhs.fpTimeStepper->Clone();
    if(fpReactionProcess) delete fpReactionProcess;
    fpReactionProcess    = rhs.fpReactionProcess->Clone();

    //assignment operator
    return *this;
}

void G4VITModel::IsApplicable(G4ITType& type1, G4ITType& type2)
{
    type1 = fType1;
    type2 = fType2;
    PrintInfo();
}

void G4VITModel::Initialize()
{
    fpReactionProcess->SetReactionTable(fpReactionTable);
    fpTimeStepper->SetReactionTable(fpReactionTable);
    fpTimeStepper->Initialize();
    fpReactionProcess->Initialize();
}
