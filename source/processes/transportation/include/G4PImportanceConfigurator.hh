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
// $Id: G4PImportanceConfigurator.hh,v 1.6 2006/06/29 21:10:00 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// ----------------------------------------------------------------------
// Class G4PImportanceConfigurator
//
// Class description:
// This class builds and places the G4PImportanceProcess.
// If the object is deleted the process is removed from the 
// process list.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4PImportanceConfigurator_hh
#define G4PImportanceConfigurator_hh G4PImportanceConfigurator_hh

#include "G4Types.hh"
#include "G4VSamplerConfigurator.hh"
#include "G4ImportanceSplitExaminer.hh"
#include "G4ProcessPlacer.hh"

class G4ParallelWorld;
class G4VImportanceAlgorithm;
class G4ParallelImportanceProcess;

class G4PImportanceConfigurator : public G4VSamplerConfigurator
{

public:  // with description

  G4PImportanceConfigurator(const G4String &particlename,
                              G4ParallelWorld &parallelWorld,
                              G4VIStore &istore,
                              const G4VImportanceAlgorithm *ialg);
  virtual ~G4PImportanceConfigurator();
  virtual void Configure(G4VSamplerConfigurator *preConf);
  virtual const G4VTrackTerminator *GetTrackTerminator() const;
  
private:

  G4PImportanceConfigurator(const G4PImportanceConfigurator &);
  G4PImportanceConfigurator &
  operator=(const G4PImportanceConfigurator &);
  G4ProcessPlacer fPlacer;
  G4ParallelWorld &fPWorld;
  G4bool fDeleteIalg;
  const G4VImportanceAlgorithm *fIalgorithm;
  G4ImportanceSplitExaminer fExaminer;
  G4ParallelImportanceProcess *fParallelImportanceProcess;

};

#endif
