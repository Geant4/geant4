//  The main program for G4IntegralHadrNucleus class
#include "G4IntegralHadrNucleus.hh"
#include "g4std/iostream"
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

