    //  The main program for G4DiffElasticHadrNucleus class


#include "G4Proton.hh"
#include "G4DynamicParticle.hh"
#include "G4DiffElasticHadrNucleus.hh"

 int main()

 {

        G4double CrossSection, ElasticCrossSec, 
                 InelasticCrossSec, Energy, Momentum,
                 Tkin, aS, Ecm, Pcm, TetCM, Pi1;

                      Pi1 = 3.1416;

        G4double      px = 0;
        G4double      py = 0;
        G4double      pz = 1000;
        G4ThreeVector aVector(px, py, pz);

        G4Proton* aProton = G4Proton::Proton();
        G4DiffElasticHadrNucleus * aCrossSection = new G4DiffElasticHadrNucleus;
        G4DynamicParticle *  aParticle = new G4DynamicParticle;
   
        aParticle->SetDefinition(aProton);

        G4double N = 28, Z = 14;

        G4Nucleus *   aNucleus = new G4Nucleus;

        aNucleus->SetParameters( N, Z);
        G4double hMass = 938.;

           for(G4int i1=1; i1<10; i1+=10)
             {
                 Tkin = i1*1000;
                 Momentum = sqrt(pow(Tkin+hMass,2)-
                   pow(hMass,2));

                 aVector.setZ(Momentum);
                 aParticle->SetMomentum(aVector);

       G4double  MomGeV   = Momentum/1000;
       G4double  TkinGeV  = Tkin/1000;
       G4double  hMassGeV = hMass/1000;

                 aS = 2*hMassGeV*N*(TkinGeV+hMassGeV)+
                      hMassGeV*hMassGeV*(1+N*N);
                 Ecm = (aS-(N*N-1)*hMassGeV*hMassGeV)/2/sqrt(aS);  
                 Pcm = sqrt(Ecm*Ecm-hMassGeV*hMassGeV);

//  G4cout << " Tkin "<<TkinGeV << " Mom " << MomGeV <<
//            " Pcm " << PcmGeV << endl;

               for(G4int i2=1; i2<100; i2++)
                  {
                    G4double Q2 = 0.005*(i2-1)*1000*1000;

             CrossSection = aCrossSection->HadrNuclDifferCrSec(
                          aParticle, aNucleus, Q2);

              TetCM = acos(1-Q2/1000/1000/2/Pcm/Pcm)*180/Pi1;

//  G4cout << i1  <<" "<<  MomGeV <<" "  <<  Q2/1000 << " "
//       << CrossSection*Pcm*Pcm/Pi1 << endl;

    G4cout <<" Mom "<<MomGeV<<" Ecm "<<Ecm<<" Pcm "<<Pcm<<" Q2 "<<
              Q2/1000/1000 << " Tet " << 
               TetCM << " Cr.-Sec. "<< CrossSection   << endl;

//   G4cout<< CrossSection<<endl;

                   }  //  i2
               }      //  i1
 
     delete aParticle;
     return 0;
  }           

