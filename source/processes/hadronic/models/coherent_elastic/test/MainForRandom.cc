//  pfnhfns
//  The main program for "G4ElasticHadrNucleusHE" class

#include "iostream"
#include  "fstream"

#include "G4Proton.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4Nucleus.hh"
#include "G4IonConstructor.hh"
#include "G4HadFinalState.hh"

#include "G4ElasticHadrNucleusHE.hh"

 int main()

 {
	G4IonConstructor Ion;
	Ion.ConstructParticle();

        G4double CrossSection, ElasticCrossSec, Q2, dE, Tkin,
                 InelasticCrossSec, Energy, Momentum;

        G4double            px = 0;
        G4double            py = 0;
        G4double            pz = 2000;
        G4ThreeVector       inVector(px, py, pz);
        G4ThreeVector       outVector, aPosition(0., 0., 0.);

        G4Proton          * aProton   = G4Proton::Proton();
        G4PionPlus        * aPionP    = G4PionPlus::
                                           PionPlusDefinition();
        G4PionMinus       * aPionM    = G4PionMinus::
                                            PionMinusDefinition();
        G4KaonPlus        * aKaonP    = G4KaonPlus::
                                            KaonPlusDefinition();
        G4KaonMinus       * aKaonM    = G4KaonMinus::
                                            KaonMinusDefinition();
        G4AntiProton      * aProtonA  = G4AntiProton::AntiProton();
        G4Lambda          * aLambda   = G4Lambda::Lambda();

        G4DynamicParticle   aParticle;
        G4HadFinalState   * HadSta    = new G4HadFinalState();   

        G4Nucleus           aNucleus, * pNucl;

        G4double        setOfNuclei[8]=
               {12.0, 16.0, 28.0, 40.0, 58.0, 64.0, 90.0, 208.0};

  const G4double N1 = 208, Z1 = 82;

        aParticle.SetDefinition(aProton);

        aNucleus.SetParameters( N1, Z1);
        pNucl   = &aNucleus;

        G4double aTime = 0.1;
        G4int    NumbPart, iq;

    G4ElasticHadrNucleusHE    aElasticRandom;

  std::ofstream TestFile("q2_for208.dat", std::ios::out);
             TestFile.setf(std::ios::scientific);
  std::ofstream TestHist("q2hist_for208.dat", std::ios::out);
             TestFile.setf(std::ios::scientific);
   
  G4double qhist[50], Hist[50];

  for(iq=1; iq<50; iq++)
    {
      Hist[iq]=0;
      qhist[iq] = iq*0.07/50*1000*1000;
    }

//  ------------ The cicle on energy of incident  -------------
        for(G4int i1=1; i1<=1; i1+=1)
        {
  const G4double N = setOfNuclei[i1];
  const G4double Z = N/2.0;
          aNucleus.SetParameters( N, Z);

          Momentum = 2000;
          Tkin     = sqrt(Momentum*Momentum+938.27*938.27)-938.27;

           inVector.setZ(Momentum);
           aParticle.SetMomentum(inVector);
           aParticle.SetKineticEnergy(Tkin);
//  -----------------------------------------------------
 const  G4HadProjectile    HadProj   = G4HadProjectile(aParticle);

//    G4ElasticHadrNucleusHE    aElasticRandom;

       G4cout <<  "  inVector " <<inVector.x()<<" " << 
  inVector.y()<<" "<<  inVector.z()<<" " << " Tkin " << Tkin <<G4endl;

                inVector =  aParticle.GetMomentum();
		inVector =  HadProj.Get4Momentum().vect();
                Tkin     =  HadProj.GetKineticEnergy();

         G4double      Tkin2 = HadProj.GetKineticEnergy();
         G4double      TotE  = HadProj.GetTotalEnergy();
         G4ThreeVector Mom   = HadProj.Get4Momentum().vect();

  G4cout<<" Tkin "<<Tkin2<<" TotE "<<TotE;
  G4cout<<" Direction Momentum "<<Mom<<G4endl;



	 for(G4int i2=1; i2<500000; i2++)  
        {

//  G4cout<<" Before Applay: NumbPart "<<NumbPart<<G4endl;
// ================================================
      HadSta   = aElasticRandom.ApplyYourself(HadProj, aNucleus);
// ================================================
      NumbPart = HadSta->GetNumberOfSecondaries();

//  --------------- Hadron's kinematics ---------------
            outVector  =  HadSta->GetMomentumChange();
   G4double       outE =  HadSta->GetEnergyChange();
   G4double       outP =  sqrt(outE*outE-938.27*938.27);
   G4double  secEnergy; // =  Had->GetKineticEnergy();
/*
   G4cout<<"Hadron outE "<<outE<<"  outP  "<<outP<<G4endl;
   G4cout<<" outVector  "<< outVector;    
   G4cout<<" outVector^2 "<<pow(outVector.x(),2)+
                            pow(outVector.y(),2)+
                            pow(outVector.z(),2)<<G4endl;
*/
//  --------------- Nucleus kinematics ----------------

   G4HadSecondary    * NuclSecond = HadSta->GetSecondary(0);
   G4DynamicParticle * SecPartNuc = NuclSecond->GetParticle();
/*
   G4cout<<" NucleusName "
         <<SecPartNuc->GetDefinition()->GetParticleName()<<G4endl;
*/
   G4ThreeVector outVectorN =  SecPartNuc->GetMomentum();

   G4double        outENucl =  SecPartNuc->GetTotalEnergy();
   G4double        outPNucl =  sqrt(outENucl*outENucl-938.*938.*N*N);
   G4double   secEnergyNucl =  SecPartNuc->GetKineticEnergy();
/*
   G4cout<<"Nucleus KinEnergy "<<secEnergyNucl<<" outE "
         <<outENucl<<"  NuclMass "<<outENucl-secEnergyNucl<<G4endl;
   G4cout<<" outVector  "<< outVectorN<<G4endl;    
   G4cout<< " outVector^2 "<<pow(outVectorN.x(),2)+
                      pow(outVectorN.y(),2)+
                      pow(outVectorN.z(),2)<<G4endl;
*/
 //  ---------------------------------------------------
 //   G4cout<<"  outVector "<<outVector<<G4endl;
             Q2      =  pow(outE-TotE,2)-
                        pow((outVector.x()*outP-inVector.x()),2)-
                        pow((outVector.y()*outP-inVector.y()),2)-
                        pow((outVector.z()*outP-inVector.z()),2);
/* 
G4cout<<" dE "<<outE-TotE
       <<" dx "<<outVector.x()*outP-inVector.x()
       <<" dy "<<outVector.y()*outP-inVector.y()
       <<" dz "<<outVector.z()*outP-inVector.z()<<G4endl;
*/
 G4cout<<" outVector " <<  sqrt(pow(outVector.x()*outP,2)+
                         pow(outVector.y()*outP,2)+
                         pow(outVector.z()*outP,2))/1000
       <<"    i2  "<<i2<<"    Q2 " <<abs(Q2)/1000/1000<<G4endl;    

//      G4cout<<" end cicle of Q2 "<< Q2<<G4endl;

      TestFile<<abs(Q2)/1000/1000<<G4endl;
/*
      iq=0;
	for(G4int ip=0; ip<50; ip++)
	{
 G4cout<<"1 iq qhist[iq] "<<iq<<"  "<<qhist[iq]<<"  "<<abs(Q2)<<G4endl;
          if(abs(Q2)>qhist[iq]) { iq++; }
        }

  iq=0; 
 G4cout<<"1 iq qhist[iq] "<<iq<<"  "<<qhist[iq]<<"  "<<abs(Q2)<<G4endl;
 
  while(abs(Q2)<qhist[iq]) 
           {
 G4cout<<"2 iq qhist[iq] "<<iq<<"  "<<qhist[iq]<<"  "<<abs(Q2)<<G4endl;
        iq++;
           }
 */

   Hist[iq]++;

//   G4cout<<"3 iq qhist[iq] "<<iq<<"  "<<qhist[iq]<<G4endl;
                   }  //  i2
      G4cout<<" end cicle of Momentum "<< Q2<<G4endl;
               }      //  i1
 
	for(iq=0; iq<50; iq++)
	  TestHist<<qhist[iq]<<"  "<<Hist[iq]<<"   "<<G4endl;

//     delete aParticle;
//     delete aChange;
//     delete aElasticRandom; 
      G4cout<<" end before return "<<G4endl;
     return 0;
  }           

