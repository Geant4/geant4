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
//  The main program for G4IntegralHadrNucleus class
#include "G4IntegralHadrNucleus.hh"
#include <iostream>
#include "G4Proton.hh"
#include "G4DynamicParticle.hh"

  main()

 {
   FILE* inFile;

   G4double CrossSection, ElasticCrossSec, InelasticCrossSec, Energy, Momentum;
   G4double px=0;
   G4double py=0;
   G4double pz=100*GeV;
   G4ThreeVector aVector(px, py, pz);
   G4Proton* aProton = G4Proton::Proton();
   int inew;

   G4DynamicParticle* aParticle = new G4DynamicParticle;
   
   aParticle->SetDefinition(aProton);

   G4IntegralHadrNucleus*  IntegrDebug = new G4IntegralHadrNucleus();

 //  inFile = fopen("Elastic.dat", "rw");

   for(G4int i1=10; i1<14;i1+=2)
       {
         inew = i1;
         
         for(G4int i2=10;i2<10000;i2+=10)
           {
              Momentum = i2*GeV;
              aVector.setZ(Momentum);
              aParticle->SetMomentum(aVector);
	      Energy = aParticle->GetTotalEnergy();
   CrossSection=IntegrDebug->GetTotalCrossSection(aParticle, i1);
   ElasticCrossSec = IntegrDebug->GetElasticCrossSection(aParticle, i1);
   InelasticCrossSec = IntegrDebug->GetInelasticCrossSection(aParticle, i1);
//       CrossSection = Energy;
       G4cout << i1  <<" "<<  Energy<<" "  <<  CrossSection << " "
       << ElasticCrossSec<<" "<<InelasticCrossSec<< " new "<< G4endl;

//        fprintf(inFile,"  %d  %f8.2, %f8.2  \n", inew, Energy,
//                          CrossSection);

            }
       }
     delete IntegrDebug;
     delete aParticle;
     delete aProton;
//   fclose(inFile);
  }           

