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
// $Id: A01EventAction.hh,v 1.3 2002-12-13 11:34:28 gunter Exp $
// --------------------------------------------------------------
//
#ifndef A01EventAction_h
#define A01EventAction_h 1

#include "AIDA/AIDA.h"
#include "G4UserEventAction.hh"
#include "globals.hh"

using namespace AIDA;

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
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    G4int HCHCID;
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    A01EventActionMessenger* messenger;
    G4int verboseLevel;

    IHistogram1D* dc1Hits;
	IHistogram1D* dc2Hits;
	ICloud2D* dc1XY;
	ICloud2D* dc2XY;
	ICloud2D* evstof;
	ITuple* tuple;
	IPlotter* plotter;

  public:
    inline void SetVerbose(G4int val) { verboseLevel = val; }
    inline G4int GetVerbose() const { return verboseLevel; }
};

#endif
