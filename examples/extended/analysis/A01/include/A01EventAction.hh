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
// $Id: A01EventAction.hh,v 1.6 2006-06-29 16:31:04 gunter Exp $
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

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  private:
    G4int HHC1ID;
    G4int HHC2ID;
    G4int DHC1ID;
    G4int DHC2ID;
    G4int ECHCID;
    G4int HCHCID;

    A01EventActionMessenger* messenger;
    G4int verboseLevel;

#ifdef G4ANALYSIS_USE
    IHistogram1D* dc1Hits;
	IHistogram1D* dc2Hits;
	ICloud2D* dc1XY;
	ICloud2D* dc2XY;
	ICloud2D* evstof;
	ITuple* tuple;
	IPlotter* plotter;
#endif // G4ANALYSIS_USE

  public:
    inline void SetVerbose(G4int val) { verboseLevel = val; }
    inline G4int GetVerbose() const { return verboseLevel; }
};

#endif
