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
// $Id: G4HadronCrossSectionPlugin.hh,v 1.5 2006-12-13 15:44:54 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Plug-in for G4CrossSectionDataTest
// F.W. Jones, TRIUMF, 30-MAR-98
//  


#ifndef G4HadronCrossSectionPlugin_h
#define G4HadronCrossSectionPlugin_h 1
 
#include "G4VCrossSectionDataSet.hh"

#include "globals.hh"
#include "G4Element.hh"
#include "G4PionPlus.hh"
#include "G4PionZero.hh"
#include "G4PionMinus.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Neutron.hh"


class G4HadronCrossSectionPlugin : public G4VCrossSectionDataSet
{
public:

   G4HadronCrossSectionPlugin()
   {
   }

   ~G4HadronCrossSectionPlugin()
   {
   }

   G4double
   GetCrossSection(const G4DynamicParticle* aParticle,
                   const G4Element* anElement, G4double aTemp)
   {
      return 999.*millibarn;
   }

   G4bool
   IsApplicable(const G4DynamicParticle* aParticle,
                const G4Element* anElement)
   {
      if (verboseLevel > 1) {
         G4cout << "G4HadronCrossSectionPlugin::IsApplicable:" << G4endl;
         G4cout << "  Particle: " <<
                 aParticle->GetDefinition()->GetParticleName() << G4endl;
         G4cout << "  Energy:   " << aParticle->GetKineticEnergy()/GeV << G4endl;
         G4cout << "  Element:  " << anElement->GetName() << G4endl;
      }
      if (aParticle->GetDefinition() != G4Proton::Proton()) return 0;
      G4double ekin = aParticle->GetKineticEnergy()/GeV;
      return (ekin > 0.95 && ekin < 1.05);
   }

   void BuildPhysicsTable(const G4ParticleDefinition&)
   {
   }

   void DumpPhysicsTable(const G4ParticleDefinition&)
   {
   }

private:

};
#endif
