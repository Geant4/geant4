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
#include "G4LENDCombinedCrossSection.hh"
#include "G4LENDElasticCrossSection.hh"
#include "G4LENDInelasticCrossSection.hh"
#include "G4LENDCaptureCrossSection.hh"
#include "G4LENDFissionCrossSection.hh"
#include "Randomize.hh"


G4LENDCombinedCrossSection::G4LENDCombinedCrossSection( G4ParticleDefinition* pd )
:G4LENDCrossSection("LENDCombinedCrossSection")
{
   proj = pd;
   elasticXS = new G4LENDElasticCrossSection( pd );   
   inelasticXS = new G4LENDInelasticCrossSection( pd );   
   captureXS = new G4LENDCaptureCrossSection( pd );   
   fissionXS = new G4LENDFissionCrossSection( pd );   
}

void G4LENDCombinedCrossSection::BuildPhysicsTable( const G4ParticleDefinition& pd )
{
   elasticXS->BuildPhysicsTable( pd );
   inelasticXS->BuildPhysicsTable( pd );
   captureXS->BuildPhysicsTable( pd );
   fissionXS->BuildPhysicsTable( pd );
   create_used_target_map();
}

G4double G4LENDCombinedCrossSection::GetIsoCrossSection( const G4DynamicParticle* dp , G4int iZ , G4int iA ,
                             const G4Isotope* isotope , const G4Element* /*elment*/ , const G4Material* material )
{
   G4double XS = 0.0;
   XS += elasticXS->GetIsoCrossSection( dp, iZ, iA, isotope, NULL , material );
   XS += inelasticXS->GetIsoCrossSection( dp, iZ, iA, isotope, NULL , material );
   XS += captureXS->GetIsoCrossSection( dp, iZ, iA, isotope, NULL , material );
   XS += fissionXS->GetIsoCrossSection( dp, iZ, iA, isotope, NULL , material );
   //G4cout << "G4LENDCombinedCrossSection::GetIsoCrossSection " << XS/CLHEP::barn << " [barn]" << G4endl;
   return XS;
}

G4int G4LENDCombinedCrossSection::SelectChannel( const G4DynamicParticle* dp , G4int iZ , G4int iA ,
                             const G4Isotope* isotope , const G4Element* /*elment*/ , const G4Material* material )
{
   G4int ichannel=-1;
   G4double XSs[4];
   XSs[0] = elasticXS->GetIsoCrossSection( dp, iZ, iA, isotope, NULL , material );
   XSs[1] = XSs[0] + inelasticXS->GetIsoCrossSection( dp, iZ, iA, isotope, NULL , material );
   XSs[2] = XSs[1] + captureXS->GetIsoCrossSection( dp, iZ, iA, isotope, NULL , material );
   XSs[3] = XSs[2] + fissionXS->GetIsoCrossSection( dp, iZ, iA, isotope, NULL , material );

   G4double total = XSs[3];

   G4double random = G4UniformRand();
   for ( G4int i = 0 ; i != 4 ; i++ ) {
      if ( random*total <= XSs[i] ) { 
         ichannel = i;
         break; 
      }
   }
   
   return ichannel;
}
