  // G4HadronValues.cc

#include "G4HadronValues.hh"
//#include "G4
 void 
 G4HadronValues::GetHadronValues(const G4DynamicParticle* aHadron)
  {

       G4ParticleDefinition * dHadron = aHadron->GetDefinition();

       G4int iHadron;

         if(dHadron == G4Proton::Proton()||
            dHadron == G4Neutron::Neutron())         iHadron=0;

   else  if(dHadron == G4AntiProton::AntiProton()||
            dHadron == G4AntiNeutron::AntiNeutron()) iHadron=1;

   else  if(dHadron == G4PionPlus::PionPlus())       iHadron=2;

   else  if(dHadron == G4PionMinus::PionMinus())     iHadron=3;

   else  if(dHadron == G4KaonPlus::KaonPlus())       iHadron=4;

   else  if(dHadron == G4KaonMinus::KaonMinus())     iHadron=5;

   else { 
         G4cout << " There is not method for the hadron " <<
                  dHadron->GetParticleName() << endl;
        }

       G4double mHadr      = aHadron->GetMass()/1000.;         // In GeV
       G4double HadrEnergy = aHadron->GetTotalEnergy()/1000.;  // In GeV
       G4double sHadr = 2*HadrEnergy*0.938+0.938*0.938+mHadr*mHadr;
       G4double sqrS  =    sqrt(sHadr);

//  G4cout<<" From hadr. Cr.-Sec. E: "<<HadrEnergy<<" Mass "<<
//                   mHadr<<" Nh "<<iHadron<<endl;

        switch (iHadron)
        {

         case 0:                  //  proton

              HadrTot   =
              -4+4.85*log(HadrEnergy)+50*pow(HadrEnergy,-0.19); //  mb
              HadrSlope = 6.44+0.88*log(sHadr);               //  GeV-2 
              HadrReIm  = 0.112*log(sHadr/300)*exp(-0.0001*sHadr);
              DDSect2   = 11;                                    //mb*GeV-2
              DDSect3   = 3;                                     //mb*GeV-2
            break;

         case 1:              //   antiproton

              sqrS      = sqrt(sHadr);
              HadrTot   =
              -4+4.85*log(HadrEnergy)+50*pow(HadrEnergy,-0.19)
                 +57.32*pow(HadrEnergy,-0.664); //  mb
              HadrSlope = 6.44+0.88*log(sHadr)+10*pow(sHadr,-0.423); //GeV-2 
              HadrReIm  =0.06*(sqrS-2.236)*(sqrS-14.14)*pow(sHadr,-1.01);
              DDSect2   = 11;                                       //mb*GeV-2
              DDSect3   = 3;                                        //mb*GeV-2
            break;

         case 2:             //   pi plus

              HadrTot   =10.6+2.*log(HadrEnergy)+25*pow(HadrEnergy,-0.43); // mb 
              HadrSlope = 7.28+0.245*log(sHadr);  //This is GeV-2 
              HadrReIm  = 0.12*log(sHadr/100)*exp(-0.001*sHadr);
              DDSect2   = 4.6;                                      //mb*GeV-2
              DDSect3   = 1.33;                                     //mb*GeV-2
            break;

         case 3:             //   pi minus

              HadrTot   = 10.6+2*log(HadrEnergy)+30*pow(HadrEnergy,-0.43);
//This is mb 
              HadrSlope = 7.28+0.245*log(sHadr);
//This is GeV-2 
              HadrReIm  = 0.12*log(sHadr/100)*exp(-0.001*sHadr);
              DDSect2   = 4.6;                                 //mb*GeV-2
              DDSect3   = 1.33;                                //mb*GeV-2
            break;

         case 4:            //  K plus

              HadrTot   = 10+1.8*log(HadrEnergy)+8*pow(HadrEnergy,-0.5);
//This is mb 
              HadrSlope  = 5.28+1.76*log(sHadr)-2.84*pow(sHadr,-0.5);
//This is GeV-2
              HadrReIm  = 7.2*(sHadr-20)*(sHadr-150)*pow(sHadr+75,-2.6);
              DDSect2   = 3.5;                                    //mb*GeV-2
              DDSect3   = 1.03;                                   //mb*GeV-2
            break;

         case 5:              //   K minus

              HadrTot   = 10+1.8*log(HadrEnergy)+25*pow(HadrEnergy,-0.5);
//This is mb 

              HadrSlope = 6.98+0.127*log(sHadr);                 
//This is GeV-2 
              HadrReIm  = 7.2*(sHadr-20)*(sHadr-150)*pow(sHadr+7,-2.6);
              DDSect2   = 3.5;                                    //mb*GeV-2
              DDSect3   = 1.03;                                   //mb*GeV-2
            break;
      }
       
  }


/* end of file */
