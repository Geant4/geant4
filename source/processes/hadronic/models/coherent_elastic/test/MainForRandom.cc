       //  The main program for "G4ElasticHadrNucleusHE" class

#include "G4Proton.hh"
#include "G4DynamicParticle.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4IonConstructor.hh"
 int main()

 {
	G4IonConstructor Ion;
	Ion.ConstructParticle();

        G4double CrossSection, ElasticCrossSec, Q2, dE, Tkin,
                 InelasticCrossSec, Energy, Momentum;

        G4double        px = 0;
        G4double        py = 500;
        G4double        pz = 100;
        G4ThreeVector   inVector(px, py, pz);
        G4ThreeVector   outVector, aPosition(0., 0., 0.);

        G4Proton *                 aProton   = G4Proton::Proton();
        G4DynamicParticle *        aParticle = new G4DynamicParticle;
                                               
        G4ParticleChange *         aChange   = new G4ParticleChange;
        G4Nucleus                  aNucleus;
        G4ElasticHadrNucleusHE *   aElasticRandom = 
                                       new G4ElasticHadrNucleusHE;
        G4Track *                  secTrack;
  
  const G4double N = 6, Z = 3;

        aParticle->SetDefinition(aProton);

        aNucleus.SetParameters( N, Z);

        G4double aTime = 0.1;

    G4Track   aTrack = G4Track( aParticle, aTime, aPosition);

           for(G4int i1=1; i1<=1; i1+=1)
             {
                 Momentum = i1*1000;
                 Tkin     = sqrt(Momentum*Momentum+938.*938.)-938.;

                 inVector.setZ(Momentum);
                 aParticle->SetMomentum(inVector);

//        G4cout <<  "  inMom befor " <<inVector.x()<<" " << 
//inVector.y()<<" "<<  inVector.z()<<" " << " Tkin " << Tkin << G4endl;

                inVector =  aParticle->GetMomentum();
//        G4cout <<  "  after inMom " <<inVector.x()<<" " << 
//inVector.y()<<" "<<  inVector.z()<<" " << " Tkin " << Tkin << G4endl;

         G4double      Tkin2 = aTrack.GetKineticEnergy();
         G4double      TotE  = aTrack.GetTotalEnergy();
         G4ThreeVector DirM  = aTrack.GetMomentumDirection();
         G4ThreeVector Mom   = aTrack.GetMomentum();

//         G4cout << " Tkin "<<Tkin2<<" TotE "<<TotE<< G4endl;

//         G4cout << " Direction " << DirM.x() << " " 
//                 <<  DirM.y() << " " << DirM.z() << G4endl;

//         G4cout  << " Direction Mom " << Mom.x() << " " 
//                 <<  Mom.y() << " " << Mom.z() << G4endl;

               for(G4int i2=1; i2<1000; i2++)
                  {

             aChange   = aElasticRandom->ApplyYourself(
                                    aTrack, aNucleus);

             outVector = *aChange->GetMomentumChange();
   G4double  outE      =  aChange->GetEnergyChange();
  
   G4double  outP      = sqrt(outE*outE-938.*938.);
 

             secTrack  =  aChange->GetSecondary(0);
   G4double  secEnergy =  secTrack->GetKineticEnergy();

       G4cout<<"         "<<secEnergy<<"         "<< outE<<G4endl;

//             G4cout<< " outVector^2 "<<pow(outVector.x(),2)+
//                      pow(outVector.y(),2)+
//                      pow(outVector.z(),2)<<G4endl;

             Q2        =  pow(outE-TotE,2)-
                        pow((outVector.x()*outP-inVector.x()),2)-
                        pow((outVector.y()*outP-inVector.y()),2)-
                        pow((outVector.z()*outP-inVector.z()),2);

//    G4cout <<" Common i1 i2"<< i1 <<i2<<" dE "<<(outE-TotE)/1000 <<
//    " outVector " <<sqrt(pow(outVector.x()*outP,2)+
//                         pow(outVector.y()*outP,2)+
//                         pow(outVector.z()*outP,2))/1000
//   <<" Q2 " <<Q2<<G4endl<<G4endl;

                   }  //  i2
               }      //  i1
 
//     delete aParticle;
//     delete aChange;
//     delete aElasticRandom; 

     return 0;
  }           

