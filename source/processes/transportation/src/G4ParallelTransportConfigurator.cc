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
// $Id: G4ParallelTransportConfigurator.cc,v 1.3 2002-11-04 10:47:56 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelTransportConfigurator
//

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#include "G4ParallelTransportConfigurator.hh"
#include "G4ParallelWorld.hh"

G4ParallelTransportConfigurator::
G4ParallelTransportConfigurator(const G4String &particlename,
			   G4ParallelWorld &pworld)
  :
  fPlacer(particlename),
  fPWorld(pworld),
  fParallelTransport(0)
{}

G4ParallelTransportConfigurator::
~G4ParallelTransportConfigurator(){
  if (fParallelTransport) {
    fPlacer.RemoveProcess(fParallelTransport);
    delete fParallelTransport;
  }
}

void G4ParallelTransportConfigurator::
Configure(G4VSamplerConfigurator *preConf){
  fParallelTransport = new 
    G4ParallelTransport(fPWorld.GetGeoDriver(), 
			fPWorld.GetParallelStepper());
  if (!fParallelTransport) {
    G4Exception("ERROR:G4ParallelTransportConfigurator::Configure new failed to create G4ParallelTransport!");
  }
  fPlacer.AddProcessAsSecondDoIt(fParallelTransport);
}

const G4VTrackTerminator *G4ParallelTransportConfigurator::
GetTrackTerminator() const {
  return 0;
}
