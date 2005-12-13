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
// 05-11-21 NeutronHP or Low Energy Parameterization Models 
//          Implemented by T. Koi (SLAC/SCCS)
//          If NeutronHP data do not available for an element, then Low Energy 
//          Parameterization models handle the interactions of the element.
//

#ifndef G4NeutronHPorLFissionModel_h
#define G4NeutronHPorLFissionModel_h 1

#include "G4HadronicInteraction.hh"
#include "G4NeutronHPorLFission.hh"
#include "G4NeutronHPNames.hh"
#include "G4LFission.hh"

class G4NeutronHPorLFissionModel : public G4HadronicInteraction
{
   public:
      G4NeutronHPorLFissionModel();
      ~G4NeutronHPorLFissionModel();

      G4HadFinalState * ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& aTargetNucleus);
      G4VCrossSectionDataSet* GiveHPXSectionDataSet() { return theHPFission->GiveXSectionDataSet(); } 

   private: 
      G4NeutronHPorLFission* theHPFission;
      G4LFission* theLFission;

      G4NeutronHPNames* theHPNames; 
};
#endif
