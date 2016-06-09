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
// $Id: G4PImportanceConfigurator.hh,v 1.5 2003/11/26 14:51:48 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
