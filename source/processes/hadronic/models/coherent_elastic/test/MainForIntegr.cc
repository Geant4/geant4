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
 //  The main program for G4IntegrHadrNucleus class

#include "G4IntegrHadrNucleus.hh"
#include "G4Proton.hh"
#include "G4PionPlus.hh"
#include "G4DynamicParticle.hh"

  int  main()

 {

   G4double CrossSection, ElasticCrossSec, InelasticCrossSec, 
            Energy, Momentum;


   G4double       px         = 0;
   G4double       py         = 0;
   G4double       pz         = 1000;
   G4ThreeVector  aVector(px, py, pz);
   G4Proton *     aProton    = G4Proton::Proton();
   G4PionPlus *   aPionPlus  = G4PionPlus::PionPlus();
   G4Nucleus *    aNucleus   = new G4Nucleus;

   G4DynamicParticle * aParticle = new G4DynamicParticle;
   
   aParticle->SetDefinition(aProton);
   aParticle->SetDefinition(aPionPlus);

   G4IntegrHadrNucleus *  IntegrDebug = new G4IntegrHadrNucleus;

       G4double N=12, Z=6;   

   for(G4int i1=12; i1<14;i1+=2)
       {
      
         N = (double) i1;
         Z = N/2;

         aNucleus->SetParameters(N, Z);
   
         for(G4int i2=1000; i2<10000; i2+=1000)
           {
              Momentum  = i2;
              aVector.setZ(Momentum);
              aParticle->SetMomentum(aVector);
	      Energy    = aParticle->GetTotalEnergy();

   CrossSection         =   IntegrDebug->GetTotalCrossSection(
                                                 aParticle, aNucleus);

   ElasticCrossSec      =   IntegrDebug->GetElasticCrossSection(
                                                 aParticle, aNucleus);

   InelasticCrossSec    =   IntegrDebug->
                    GetInelasticCrossSection(aParticle, aNucleus);

       G4cout << i1  <<" "<<  Energy<<" "  <<  CrossSection << " "
       << ElasticCrossSec<<" "<<InelasticCrossSec<< " new "<< endl;

            }   //  i2
       }        //  i1

     delete aParticle;

     return 0;
  }           

