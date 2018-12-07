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
/// \file hadronic/Hadr02/include/IonUrQMDPhysics.hh
/// \brief Definition of the IonUrQMDPhysics class
//
// GRAS tag Name:
//
//---------------------------------------------------------------------------
//
// Header:    G4IonUrQMDPhysics
//
// Author:    2012 Andrea Dotti
//
// 
//
// Modified:     
//
// ------------------------------------------------------------
//

#ifndef G4IonUrQMDPhysics_h
#define G4IonUrQMDPhysics_h 1

#include "G4VHadronPhysics.hh"
#include "globals.hh"

class G4HadronicInteraction;
class G4VCrossSectionDataSet;
class G4UrQMD1_3Model;

class IonUrQMDPhysics : public G4VHadronPhysics
{
public:

  IonUrQMDPhysics(G4int ver = 0);
  virtual ~IonUrQMDPhysics();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  virtual void ConstructProcess();

private:

  void AddProcess(const G4String&, G4ParticleDefinition*, G4bool isIon);

  G4VCrossSectionDataSet*   fTripathi;
  G4VCrossSectionDataSet*   fTripathiLight;
  G4VCrossSectionDataSet*   fShen;
  G4VCrossSectionDataSet*   fIonH;
  
  G4UrQMD1_3Model * fModel;

  G4int  fVerbose;
  G4bool fWasActivated;
};


#endif








