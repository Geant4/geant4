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
/// \file hadronic/Hadr02/include/IonDPMJETPhysics.hh
/// \brief Definition of the IonDPMJETPhysics class
//
// $Id$
// GRAS tag $Name: gras-02-05-02 $
//
//---------------------------------------------------------------------------
//
// Header:    IonDPMJETPhysics
//
// Author:    copy from P.Truscott manuel DPMJET2.5 
//
// 
// Customer:          
// Contract:          
//
// Modifications are provided according to
//
// Organisation:        
// Customer:            
// Contract:            
//
// Modified:     26.08.2010
//
// ------------------------------------------------------------
//

#ifndef IonDPMJETPhysics_h
#define IonDPMJETPhysics_h 1

#include "G4VHadronPhysics.hh"
#include "globals.hh"

class G4BinaryLightIonReaction;
class G4DPMJET2_5Model;
class G4DPMJET2_5CrossSection;
class G4VCrossSectionDataSet;

class IonDPMJETPhysics : public G4VHadronPhysics
{
public:

  IonDPMJETPhysics(G4bool val);
  virtual ~IonDPMJETPhysics();

  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type
  void ConstructProcess();

private:

  void AddProcess(const G4String& name, G4ParticleDefinition* part,
                  G4bool isIon);

  G4VCrossSectionDataSet* fTripathi;
  G4VCrossSectionDataSet* fTripathiLight;
  G4VCrossSectionDataSet* fShen;
  G4VCrossSectionDataSet* fIonH;
  G4BinaryLightIonReaction*  fIonBC;
  G4DPMJET2_5Model*          fDPM;
  G4DPMJET2_5CrossSection*   fDpmXS;
  G4bool                  fUseDPMJETXS;
};

#endif
