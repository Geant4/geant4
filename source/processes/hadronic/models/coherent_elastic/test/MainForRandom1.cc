//  The main program for "G4ElasticHadrNucleusHE" class


#include "G4Proton.hh"
#include "G4DynamicParticle.hh"
//#include "G4ParticleChange.hh"
//#include "G4Track.hh"
//#include "G4ThreeVector.hh"
#include "G4Nucleus.hh"
//#include "G4IonConstructor.hh"
#include "G4ElasticHadrNucleusHE.hh"
//#include "iostream.h"

 int main()

 {
//	G4IonConstructor Ion;
//	Ion.ConstructParticle();

        G4double CrossSection, ElasticCrossSec, Q2, dE, Tkin,
                 InelasticCrossSec, Energy, Momentum;

        G4double        px = 0;
        G4double        py = 1500;
        G4double        pz = 1000;
        G4ThreeVector   inVector(px, py, pz);
        G4ThreeVector   outVector, aPosition(0., 0., 0.);

        G4Proton          *    aProton   = G4Proton::Proton();
        G4DynamicParticle *    aParticle = new G4DynamicParticle;
   
        G4VParticleChange *    aChange   = new G4VParticleChange;
        G4Nucleus              aNucleus, * pNucl;

        G4Track           *    secTrack;

  const G4double N = 12, Z = 6;

        aParticle->SetDefinition(aProton);

        aNucleus.SetParameters( N, Z);
        pNucl   = &aNucleus;

        G4double aTime = 0.1;

           for(G4int i1=1; i1<=1; i1+=1)
             {
                 Momentum = i1*5000;
                 Tkin     = sqrt(Momentum*Momentum+938.*938.)-938.;

                 inVector.setZ(Momentum);
                 aParticle->SetMomentum(inVector);

//    G4cout << " Tkin "<< Tkin <<endl;

        G4ElasticHadrNucleusHE    aElasticRandom(aParticle, pNucl);

        G4cout <<  "  inMom befor " <<inVector.x()<<" " << 
inVector.y()<<" "<<  inVector.z()<<" " << " Tkin " << Tkin << endl;

                inVector =  aParticle->GetMomentum();

//        G4cout <<  "  after inMom " <<inVector.x()<<" " << 
//inVector.y()<<" "<<  inVector.z()<<" " << " Tkin " << Tkin << endl;

               for(G4int i2=1; i2<10; i2++)
                  {
          Q2   = aElasticRandom.RandomElastic0(aParticle,  pNucl);

  G4cout << " i2, Q2   "<<i2<<"    "<<Q2<<endl;   // For Random
//     G4cout << Q2/1000/1000 <<endl;         

                   }  //  i2

//      G4cout<<" end cicle of Q2 "<< Q2<<endl;

               }      //  i1
 
     return 0;
  }           

