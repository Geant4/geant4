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
// 05-11-21 NeutronHP or Low Energy Parameterization Models 
//          Implemented by T. Koi (SLAC/SCCS)
//          If NeutronHP data do not available for an element, then Low Energy 
//          Parameterization models handle the interactions of the element.
//

#ifndef G4NeutronHPorLEInelasticModel_h
#define G4NeutronHPorLEInelasticModel_h 1

#include "G4HadronicInteraction.hh"
//#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPorLEInelastic.hh"
#include "G4NeutronHPNames.hh"
#include "G4LENeutronInelastic.hh"

class G4NeutronHPorLEInelasticModel : public G4HadronicInteraction
{
   public:
      G4NeutronHPorLEInelasticModel();
      ~G4NeutronHPorLEInelasticModel();

      G4HadFinalState * ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus);
      G4VCrossSectionDataSet* GiveHPXSectionDataSet() { return theHPInelastic->GiveXSectionDataSet(); } 

   private: 
      //G4NeutronHPElastic* theHPElastic;
      G4NeutronHPorLEInelastic* theHPInelastic;
      G4LENeutronInelastic* theLEInelastic;

      G4NeutronHPNames* theHPNames; 
};
#endif
