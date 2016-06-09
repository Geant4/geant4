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
// $Id: G4MWeightWindowConfigurator.cc,v 1.4 2003/11/26 14:51:49 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// ----------------------------------------------------------------------
// Class G4MWeightWindowConfigurator
//
// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#include "G4MWeightWindowConfigurator.hh"
#include "G4WeightWindowAlgorithm.hh"
#include "G4MassWeightWindowProcess.hh"

G4MWeightWindowConfigurator::
G4MWeightWindowConfigurator(const G4String &particlename,
                            G4VWeightWindowStore &wwstore,
                            const G4VWeightWindowAlgorithm *wwAlg,
                            G4PlaceOfAction placeOfAction)
  : fPlacer(particlename),
    fWeightWindowStore(wwstore),
    fDeleteWWalg( ( ! wwAlg) ),
    fWWalgorithm(( (fDeleteWWalg) ? 
                   new G4WeightWindowAlgorithm(5,3,5) : wwAlg)),
    fMassWeightWindowProcess(0),
    fPlaceOfAction(placeOfAction)
{
}

G4MWeightWindowConfigurator::~G4MWeightWindowConfigurator()
{  
  if (fMassWeightWindowProcess)
  {
    fPlacer.RemoveProcess(fMassWeightWindowProcess);
    delete fMassWeightWindowProcess;
  }
  if (fDeleteWWalg)
  {
    delete fWWalgorithm;
  }
}

void
G4MWeightWindowConfigurator::Configure(G4VSamplerConfigurator *preConf)
{
  const G4VTrackTerminator *terminator = 0;
  if (preConf)
  {
    terminator = preConf->GetTrackTerminator();
  };

  fMassWeightWindowProcess = 
    new G4MassWeightWindowProcess(*fWWalgorithm, 
                                  fWeightWindowStore, 
                                  terminator,
                                  fPlaceOfAction);
  fPlacer.AddProcessAsSecondDoIt(fMassWeightWindowProcess);
}

const G4VTrackTerminator *
G4MWeightWindowConfigurator::GetTrackTerminator() const 
{
  return fMassWeightWindowProcess;
}

