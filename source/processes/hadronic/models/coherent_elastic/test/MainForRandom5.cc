	//  The main program for "G4ElasticHadrNucleusHE" class

#include "G4Proton.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4Nucleus.hh"
#include "G4IonConstructor.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "fstream.h"
#include "iostream.h"
#include "g4std/strstream"

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

        G4PionPlus   *       aPionP    = G4PionPlus::
                                            PionPlusDefinition();
        G4PionMinus  *       aPionM    = G4PionMinus::
                                            PionMinusDefinition();
        G4KaonPlus   *       aKaonP    = G4KaonPlus::
                                            KaonPlusDefinition();
        G4KaonMinus  *       aKaonM    = G4KaonMinus::
                                            KaonMinusDefinition();
        G4Proton     *       aProton   = G4Proton::Proton();
        G4AntiProton *       aProtonA  = G4AntiProton::AntiProton();

        G4DynamicParticle *  aParticle = new G4DynamicParticle;
   
        G4VParticleChange *  aChange   = new G4VParticleChange;
        G4Nucleus            aNucleus, * pNucl;

        G4Track *                  secTrack;

  const G4double N = 8, Z = 4;

        aParticle->SetDefinition(aProton);

        aNucleus.SetParameters( N, Z);
        pNucl   = &aNucleus;
  
       G4std::ostrstream osTest;
       osTest<<" Atomic Number "<<N<<ends;
       G4cout<<osTest.str()<<endl;

  G4DynamicParticle * dpNucleus =
                        aNucleus.ReturnTargetParticle();

      G4double aTime = 0.1;

//  -----  Declaration of the reading file class (constr. 5)  -------

       G4ElasticHadrNucleusHE    
                aElasticRandom(aParticle, pNucl, 137);

   cout<<"  The array is created "<<endl;
//    G4Track   aTrack = G4Track( aParticle, aTime, aPosition);

           for(G4int i1=1; i1<=2; i1+=5)
             {
                 Momentum = i1*100;
                 Tkin     = sqrt(Momentum*Momentum+938.*938.)-938.;

                 inVector.setZ(Momentum);
                 inVector.setY(100);
                 aParticle->SetMomentum(inVector);
    G4Track   aTrack = G4Track( aParticle, aTime, aPosition);

        G4cout <<  "  inMom befor " <<inVector.x()<<" " << 
inVector.y()<<" "<<  inVector.z()<<" " << " Tkin " << Tkin << endl;

                inVector =  aParticle->GetMomentum();
//        G4cout <<  "  after inMom " <<inVector.x()<<" " << 
//inVector.y()<<" "<<  inVector.z()<<" " << " Tkin " << Tkin << endl;


         G4double      Tkin2 = aTrack.GetKineticEnergy();
         G4double      TotE  = aTrack.GetTotalEnergy();
         G4ThreeVector DirM  = aTrack.GetMomentumDirection();
         G4ThreeVector Mom   = aTrack.GetMomentum();

         G4cout << " Tkin "<<Tkin2<<" TotE "<<TotE<< endl;

         G4cout << " Direction M" << DirM<< endl;

         G4cout  << " Direction Momentum " <<Mom<< endl;

//  ----------  This is for randomizi of q^2 -------------
             ofstream TestCount("q2from5.dat");
                      TestCount.setf(ios::scientific);

               for(G4int i2=1; i2<15; i2++)
                  {

////   G4cout<<" Before Applay "<<endl;

             aChange   = aElasticRandom.ApplyYourself(
                                    aTrack, aNucleus);

//++          Q2   = aElasticRandom.RandomElastic1(
//++                                  aParticle,  pNucl);

////   G4cout<<" After Apply "<<endl;

//++   G4cout << " i2, Q2 "<<i2<<"  "<<Q2<<endl;   // For Random
         
//++            TestCount<<Q2<<endl;

/*                   outVector  =  aTrack.GetMomentum();

         G4double  outE       =  aTrack.GetTotalEnergy();
         G4double  outP       =  sqrt(outE*outE-938.*938.);
 
                   secTrack   =  aChange->GetSecondary(0);
         G4double  secEnergy  =  secTrack->GetKineticEnergy();

   G4cout<<" SecEnergy "<<secEnergy<<" outE "<< outE<<endl;
   G4cout<<" outVector  "<< outVector <<endl;    
   G4cout<< " outVector^2 "<<pow(outVector.x(),2)+
                      pow(outVector.y(),2)+
                      pow(outVector.z(),2)<<endl;

             Q2      =  pow(outE-TotE,2)-
                        pow((outVector.x()*outP-inVector.x()),2)-
                        pow((outVector.y()*outP-inVector.y()),2)-
                        pow((outVector.z()*outP-inVector.z()),2);

  G4cout <<" Common i1 i2 "<< i1<<" " <<i2<<" dE "<<(outE-TotE)/1000 <<
  " outVector " <<  sqrt(pow(outVector.x()*outP,2)+
                         pow(outVector.y()*outP,2)+
                         pow(outVector.z()*outP,2))/1000
   <<" Q2 " <<Q2<<endl<<endl;    
*/

                   }  //  i2

      G4cout<<" end cicle of Q2 "<<endl;

               }      //  i1
//++   TestCount.close();

     return 0;
  }           

