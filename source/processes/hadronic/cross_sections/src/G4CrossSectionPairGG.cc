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
// $Id: G4CrossSectionPairGG.cc 94961 2016-01-08 16:31:48Z gcosmo $
// $ GEANT4 tag $Name: not supported by cvs2svn $
//
//   Class G4CrossSectionPairGG
//
//     smoothly join two cross section sets by scaling the second at a given 
//       transition energy to match the first.
//
//  Author:  Gunter Folger
//           November 2009
//

#include "G4CrossSectionPairGG.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadTmpUtil.hh"
#include "G4NistManager.hh"
#include "G4ThreeVector.hh"
#include "G4NistManager.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"


G4CrossSectionPairGG::G4CrossSectionPairGG(G4VCrossSectionDataSet* low,
      G4double Etransit) :
      G4VCrossSectionDataSet("G4CrossSectionPairGG"), theLowX(low), ETransition(
            Etransit) {
   NistMan = G4NistManager::Instance();
   theHighX = new G4ComponentGGHadronNucleusXsc();
   verboseLevel = 0;
}

G4CrossSectionPairGG::~G4CrossSectionPairGG() {
   delete theHighX;
   // The cross section registry will delete theLowX
}

void G4CrossSectionPairGG::CrossSectionDescription(
      std::ostream& outFile) const {
   outFile << "G4CrossSectionPairGG is used to add the relativistic rise to\n"
         << "hadronic cross section data sets above a given energy.  In this\n"
         << "case, the Glauber-Gribov cross section is used above 91 GeV.\n"
         << "At this energy the low energy cross section is smoothly joined\n"
         << "to the high energy cross section.  Below 91 GeV the Barashenkov\n"
         << "cross section is used for pions (G4PiNuclearCrossSection), the\n"
         << "Axen-Wellisch cross section is used for protons\n"
         << "(G4ProtonInelasticCrossSection), and the Wellisch-Laidlaw cross\n"
         << "section is used for neutrons (G4NeutronInelasticCrossSection).\n";
}

G4bool G4CrossSectionPairGG::IsElementApplicable(
      const G4DynamicParticle* aParticle, G4int Z, const G4Material* mat) {
   G4bool isApplicable(false);
   G4double Ekin = aParticle->GetKineticEnergy();
   if (Ekin <= ETransition) {
      isApplicable = theLowX->IsElementApplicable(aParticle, Z, mat);
   } else if (Z > 1) {
      isApplicable = true;
   }
   return isApplicable;
}

G4double G4CrossSectionPairGG::GetElementCrossSection(
      const G4DynamicParticle* aParticle, G4int ZZ, const G4Material* mat)
{
   G4double Xsec(0.);

   if (aParticle->GetKineticEnergy() < ETransition)
   {
      Xsec = theLowX->GetElementCrossSection(aParticle, ZZ, mat);
   } else {

      std::vector<ParticleXScale>::iterator iter = scale_factors.begin();
      const G4ParticleDefinition * pDef = aParticle->GetDefinition();
      while (iter != scale_factors.end() && (*iter).first != pDef)  /* Loop checking, 08.01.2016, W. Pokorski */
      {
         ++iter;
      }
      if (iter != scale_factors.end() )
      {
         G4int AA = G4lrint(NistMan->GetAtomicMassAmu(ZZ));
         Xsec = theHighX->GetInelasticGlauberGribov(aParticle, ZZ, AA)
				            * (*iter).second[ZZ];
         if (verboseLevel > 2)
         {
            G4cout << " scaling .." << ZZ << " " << AA << " "
                  << (*iter).second[ZZ] << " "
                  << theHighX->GetInelasticGlauberGribov(aParticle, ZZ, AA)
                  << "  " << Xsec << G4endl;
         }
      } else {
         // BuildPhysicsTable not done for pDef=aParticle->GetDefinition
         //  build table, and recurse
         BuildPhysicsTable(*pDef);
         Xsec=GetElementCrossSection(aParticle, ZZ, mat);
      }
   }

   return Xsec;
}

void G4CrossSectionPairGG::BuildPhysicsTable(const G4ParticleDefinition& pDef) {
   theLowX->BuildPhysicsTable(pDef);
   theHighX->BuildPhysicsTable(pDef);

   if (verboseLevel > 0) {
      G4cout << "G4CrossSectionPairGG::BuildPhysicsTable "
            << theLowX->GetName() << "  " << theHighX->GetName() << G4endl;
   }

   const G4ParticleDefinition * myDef = &pDef;
   std::vector<ParticleXScale>::iterator iter;
   iter = scale_factors.begin();
   while (iter != scale_factors.end() && (*iter).first != myDef)  /* Loop checking, 08.01.2016, W. Pokorski */
     {
       ++iter;
     }

   //  new particle, initialise

   G4Material* mat = 0;

   if (iter == scale_factors.end()) {
      XS_factors factors(93);
      G4ThreeVector mom(0.0, 0.0, 1.0);
      G4DynamicParticle DynPart(myDef, mom, ETransition); // last is kinetic Energy

      if (verboseLevel > 0) {
         G4cout << "G4CrossSectionPairGG::BuildPhysicsTable for particle "
               << pDef.GetParticleName() << G4endl;
      }
      for (G4int aZ = 1; aZ < 93; ++aZ) {
         factors[aZ] = 1.; // default, to give reasonable value if only high is applicable
         G4int AA = G4lrint(NistMan->GetAtomicMassAmu(aZ));
         G4bool isApplicable = theLowX->IsElementApplicable(&DynPart, aZ,
               mat) && (aZ > 1);

         if (isApplicable) {
            factors[aZ] = theLowX->GetElementCrossSection(&DynPart, aZ, mat)
						      / theHighX->GetInelasticGlauberGribov(&DynPart, aZ, AA);

         }
         if (verboseLevel > 0) {
            G4cout << "Z=" << aZ << ",  A=" << AA << ", scale="
                  << factors[aZ];
            if (verboseLevel == 1) {
               G4cout << G4endl;
            } else {
               if (isApplicable) {
                  G4cout << ",  low / high "
                        << theLowX->GetElementCrossSection(&DynPart, aZ,
                              mat) << "  "
                              << theHighX->GetInelasticGlauberGribov(&DynPart,
                                    aZ, AA) << G4endl;
               } else {
                  G4cout << ",   N/A" << G4endl;
               }
            }
         }
      }
      ParticleXScale forPart(myDef, factors);
      scale_factors.push_back(forPart);
   }
}

/*
 void G4CrossSectionPairGG::DumpHtml(const G4ParticleDefinition&,
 std::ofstream outFile)
 {
 outFile << "         <li><b>"
 << " G4CrossSectionPairGG: " << theLowX->GetName() << " cross sections \n";
 outFile << "below " << ETransition/GeV << " GeV, Glauber-Gribov above \n";
 }
 */

void G4CrossSectionPairGG::DumpPhysicsTable(const G4ParticleDefinition&) {
   G4cout << std::setw(24) << " " << " G4CrossSectionPairGG: "
         << theLowX->GetName() << " cross sections " << G4endl;
   G4cout << std::setw(27) << " " << "below " << ETransition / GeV
         << " GeV, Glauber-Gribov above " << G4endl;
}
