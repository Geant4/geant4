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
// $Id: G4ParallelTransportConfigurator.hh,v 1.1 2002-10-10 13:25:30 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelTransportConfigurator
//
// Class description:
// This class builds and places  the G4ParallelTransport
// If the object is deleted the process is removed from the 
// process list.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4ParallelTransportConfigurator_hh
#define G4ParallelTransportConfigurator_hh G4ParallelTransportConfigurator_hh

#include "G4VSamplerConfigurator.hh"

#include "globals.hh"
#include "G4ParallelTransport.hh"
#include "G4ProcessPlacer.hh"


class G4ParallelWorld;

class G4ParallelTransportConfigurator : public G4VSamplerConfigurator{
public:
  G4ParallelTransportConfigurator(const G4String &particlename,
			     G4ParallelWorld &pworld);
  ~G4ParallelTransportConfigurator();
  void Configure(G4VSamplerConfigurator *preConf);
  G4VTrackTerminator *GetTrackTerminator();
  
private:
  G4ProcessPlacer fPlacer;
  G4ParallelWorld &fPWorld;
  G4ParallelTransport *fParallelTransport;
};


#endif
