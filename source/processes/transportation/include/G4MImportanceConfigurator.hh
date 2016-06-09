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
// $Id: G4MImportanceConfigurator.hh,v 1.5 2003/11/26 14:51:47 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// ----------------------------------------------------------------------
// Class G4MImportanceConfigurator
//
// Class description:
// This class builds and places the G4ImportanceProcess.
// If the object is deleted the process is removed from the 
// process list.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef G4MImportanceConfigurator_hh
#define G4MImportanceConfigurator_hh G4MImportanceConfigurator_hh

#include "G4Types.hh"
#include "G4ProcessPlacer.hh"
#include "G4VSamplerConfigurator.hh"

class G4MassImportanceProcess;
class G4VImportanceAlgorithm;
class G4VIStore;

class G4MImportanceConfigurator : public G4VSamplerConfigurator
{

public:  // with description

  G4MImportanceConfigurator(const G4String &particlename,
                            G4VIStore &istore,
                            const G4VImportanceAlgorithm *ialg);

  virtual ~G4MImportanceConfigurator();
  virtual void Configure(G4VSamplerConfigurator *preConf);
  virtual const G4VTrackTerminator *GetTrackTerminator() const;

private:

  G4MImportanceConfigurator(const G4MImportanceConfigurator &);
  G4MImportanceConfigurator &
  operator=(const G4MImportanceConfigurator &);

  G4ProcessPlacer fPlacer;
  G4VIStore &fIStore;
  G4bool fDeleteIalg;
  const G4VImportanceAlgorithm *fIalgorithm;
  G4MassImportanceProcess *fMassImportanceProcess;
};


#endif
