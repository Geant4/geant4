//  The main program for "G4ElasticHadrNucleusHE" class

#include "G4Proton.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4Nucleus.hh"
#include "G4IonConstructor.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "g4std/fstream"

 int main()


 {
	G4IonConstructor Ion;
	Ion.ConstructParticle();

        G4double CrossSection, ElasticCrossSec, Q2, dE, Tkin,
                 InelasticCrossSec, Energy, Momentum;

        G4double        px = 0;
        G4double        py = 0;
        G4double        pz = 1000;
        G4ThreeVector   inVector(px, py, pz);
        G4ThreeVector   outVector, aPosition(0., 0., 0.);

        G4PionPlus        *  aPionP    = G4PionPlus::
                                            PionPlusDefinition();
        G4PionMinus       *  aPionM    = G4PionMinus::
                                            PionMinusDefinition();
        G4KaonPlus        *  aKaonP    = G4KaonPlus::
                                            KaonPlusDefinition();
        G4KaonMinus       *  aKaonM    = G4KaonMinus::
                                            KaonMinusDefinition();
        G4Proton          *  aProton   = G4Proton::Proton();
        G4AntiProton      *  aProtonA  = G4AntiProton::AntiProton();

        G4DynamicParticle *  aParticle = new G4DynamicParticle;
   
        G4VParticleChange *  aChange   = new G4VParticleChange;
        G4Nucleus            aNucleus, * pNucl;

  const G4double N = 90, Z = 45;

        aParticle->SetDefinition(aProton);

        aNucleus.SetParameters( N, Z);
        pNucl   = &aNucleus;

        G4double aTime = 0.1;

                 Momentum = 5000;
                 Tkin     = sqrt(Momentum*Momentum+938.*938.)-938.;

                 inVector.setZ(Momentum);
                 aParticle->SetMomentum(inVector);

  G4cout<<" Befor constructor N, Z : "<<N<<"  "<<Z<<G4endl;

       G4ElasticHadrNucleusHE    aElasticRandom(aParticle, pNucl,
                                  2000., 10000., 5);

  G4cout<<G4endl<< " The array is created !!! "<<G4endl;

       G4Track *    secTrack, aTrack;

           for(G4int i1=1; i1<=2; i1+=5)
             {
                 Momentum = i1*100000;
                 Tkin     = sqrt(Momentum*Momentum+938.*938.)-938.;

                 inVector.setZ(Momentum);
                 aParticle->SetMomentum(inVector);
        aTrack = G4Track( aParticle, aTime, aPosition);

                inVector =  aParticle->GetMomentum();

            G4std::ofstream TestFile("q2from4.dat", G4std::ios::out);
            TestFile.precision(9);
            TestFile.setf(G4std::ios::scientific);

               for(G4int i2=1; i2<500001; i2++)
                  {

          Q2   = aElasticRandom.RandomElastic1(
                                  aParticle,  pNucl);

  if(Q2<0 || Q2 > 250000) {
  G4cout << " i2, Q2   "<<i2<<"    "<<Q2<<endl;   // For Random
//        Q2 = 20000; 
                           }
          TestFile<<"  "<<Q2<<G4endl;
         
                   }  //  i2
      G4cout<<" end cicle of Q2 "<<endl;
               }      //  i1
 
//     delete aParticle;
//     delete aChange;
//     delete aElasticRandom; 

     return 0;
  }           

