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
/// \file analysis/A01/include/A01EventAction.hh
/// \brief Definition of the A01EventAction class
//
// $Id$
// --------------------------------------------------------------
//
#ifndef A01EventAction_h
#define A01EventAction_h 1


#include "G4UserEventAction.hh"
#include "globals.hh"

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
using namespace AIDA;
#endif // G4ANALYSIS_USE

class A01EventActionMessenger;

class A01EventAction : public G4UserEventAction
{
  public:
    A01EventAction();
    virtual ~A01EventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  private:
    G4int fHHC1ID;
    G4int fHHC2ID;
    G4int fDHC1ID;
    G4int fDHC2ID;
    G4int fECHCID;
    G4int fHCHCID;

    A01EventActionMessenger* fMessenger;
    G4int fVerboseLevel;

#ifdef G4ANALYSIS_USE
    IHistogram1D* fDc1Hits;
    IHistogram1D* fDc2Hits;
    ICloud2D* fDc1XY;
    ICloud2D* fDc2XY;
    ICloud2D* fEvstof;
    ITuple* fTuple;
    IPlotter* fPlotter;
#endif // G4ANALYSIS_USE

  public:
    inline void SetVerbose(G4int val) { fVerboseLevel = val; }
    inline G4int GetVerbose() const { return fVerboseLevel; }
};

#endif
