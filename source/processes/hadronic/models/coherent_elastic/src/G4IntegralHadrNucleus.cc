// IntegralHadrNucleus.cc

#include "G4IntegralHadrNucleus.hh"

 void 
 G4IntegralHadrNucleus::GetHadronValues(G4int    iHadron,
                                        G4double HadrEnergy)
  {
       G4double  Seta, Seps;

       G4double  sHadr = 2*HadrEnergy*0.938+0.938*0.938+mHadr*mHadr;
       G4double  sqrS  = sqrt(sHadr);

        switch (iHadron)
        {

         case 0:                  //  proton

              HadrTot   =
              -4+4.85*log(HadrEnergy)+50*pow(HadrEnergy,-0.19); //  mb
              HadrSlope = 6.44+0.88*log(sHadr);               //  GeV-2 
              HadrReIm  = 0.112*log(sHadr/300)*exp(-0.0001*sHadr);
              DDSect2   = 11;                                    //mb*GeV-2
              DDSect3   = 3;                                           //mb*GeV-2
            break;

         case 1:              //   antiproton

              HadrTot   =
              -4+4.85*log(HadrEnergy)+50*pow(HadrEnergy,-0.19)
                 +57.32*pow(HadrEnergy,-0.664); //  mb
              HadrSlope = 6.44+0.88*log(sHadr)+10*pow(sHadr,-0.423); //GeV-2 
              HadrReIm  =0.06*(sqrS-2.236)*(sqrS-14.14)*pow(sHadr,-1.01);
              DDSect2   = 11;                                           //mb*GeV-2
              DDSect3   = 3;                                            //mb*GeV-2
            break;

         case 2:             //   pi plus

              HadrTot   =10.6+2.*log(HadrEnergy)+25*pow(HadrEnergy,-0.43); // mb 
              HadrSlope = 7.28+0.245*log(sHadr);  //This is GeV-2 
              HadrReIm  = 0.12*log(sHadr/100)*exp(-0.001*sHadr);
              DDSect2   = 4.6;                                         //mb*GeV-2
              DDSect3   = 1.33;                                        //mb*GeV-2
            break;

         case 3:             //   pi minus

              HadrTot   = 10.6+2*log(HadrEnergy)+30*pow(HadrEnergy,-0.43);
//This is mb 
              HadrSlope = 7.28+0.245*log(sHadr);
//This is GeV-2 
              HadrReIm  = 0.12*log(sHadr/100)*exp(-0.001*sHadr);
              DDSect2   = 4.6;                                         //mb*GeV-2
              DDSect3   = 1.33;                                //mb*GeV-2
            break;

         case 4:            //  K plus

              HadrTot   = 10+1.8*log(HadrEnergy)+8*pow(HadrEnergy,-0.5);
//This is mb 
              HadrSlope  = 5.28+1.76*log(sHadr)-2.84*pow(sHadr,-0.5);
//This is GeV-2
              HadrReIm  = 7.2*(sHadr-20)*(sHadr-150)*pow(s+75,-2.6);
              DDSect2   = 3.5;                                         //mb*GeV-2
              DDSect3   = 1.03;                                        //mb*GeV-2
            break;

         case 5:              //   K minus

              HadrTot   = 10+1.8*log(HadrEnergy)+25*pow(HadrEnergy,-0.5);
//This is mb 

              HadrSlope = 6.98+0.127*log(sHadr);                 
//This is GeV-2 
              HadrReIm  = 7.2*(sHadr-20)*(sHadr-150)*pow(sHadr+7,-2.6);
              DDSect2   = 3.5;                                         //mb*GeV-2
              DDSect3   = 1.03;                                        //mb*GeV-2
            break;
      }
       
  }

 void   G4IntegralHadrNucleus::GetIntegralCrSec(G4int   Anucleus)
 {

    G4int      i, k, l, m;
    G4double   N, N1, N2, N3, N4, Delta, Inel1, Inel3;
    G4double   Tot0, Inel0, Inel2, Prod0, Prod1, ak, Delt, Delt2, Delt3;
    G4double   Rnucl, R0, Stot, Bhad, Asq, MbToB, Pi1;
    G4double   Dtot, Dinel, Dprod, Rnuc2, RB, R2B, bk, bd;

          MbToB   = 2.568;
          Pi1     = 3.1416;

          Stot    = HadrTot*MbToB;                       //{In GeV-2}
          Bhad    = HadrSlope;                           //{In GeV-2}
          Asq     = 1+HadrReIm*HadrReIm;

          R0      = sqrt(0.99);                        //{ This is fermi}
          if (Anucleus >10)  R0 = sqrt(0.84);          
          if (Anucleus >20)  R0 = sqrt((35.34+0.5*Anucleus)/(40.97+Anucleus));

          Rnucl   = R0*pow(Anucleus,0.3333);            //{In Fermi }
          Rnuc2   = Rnucl*Rnucl*MbToB*10;                 //{ In GeV-4}
          RB      = Rnuc2+Bhad;
          R2B     = RB+Bhad;
          Delta   = Stot/R2B/2/Pi1;
          Delt    = Delta*Asq*0.5;
          Delt2   = Delta*2;
          Delt3   = Stot/RB/Bhad/16/Pi1*Asq*R2B;

          Tot0    = 0; Inel0 = 0; Inel1=0;
          N = N1  = -1/Delta; 
          N3      = -1/Delt2;

          for (i=1; i<= Anucleus; i++)

             {
               N     = -N*Delta*(Anucleus-i+1)/i;
               N1    = -N1*Delta*(2*Anucleus-i+1)/i;
               N3    = -N3*Delt2*(Anucleus-i+1)/i;
               Tot0  = Tot0+N/i;
               Inel0 = Inel0+N1/i;
               N2    = 1;
               N4    = 1;
               Inel1 = 0;
               Prod1 = 0;

                  for (l=0; l<= i; l++) 

                     {
//                       Inel1   = Inel1+N2/(i+l);
                         Prod1 = Prod1+N4*RB/(i*RB+l*Bhad);
//                       N2      = -N2*Delt*(i-l)/(l+1);
                         N4    = -N4*Delt3*(i-l)/(l+1);
                     } //  l

//                 Inel2      = Inel2+Inel1*N3;
                 Prod0    = Prod0+Prod1*N3;
          if(N1/i/Inel0<0.0001)  break;
              }  // i

             Tot0         = Tot0*HadrTot;
             Inel0        = Inel0*HadrTot*0.5;
//             Inel2      = Inel2*HadrTot;
             Prod0        = Prod0*HadrTot;

             ak           = (Rnuc2*2*Pi1/Stot);
             DDSect1      = (DDSect2+
                              DDSect3*log(1.06*2*HadrEnergy/Rnucl/sqrt(25.68)/4));
             Dtot         =  8*Pi1*ak/HadrTot*(1-(1+Anucleus/ak)*exp(-Anucleus/ak))*
                                    DDSect1/MbToB;

             bk           = (1-1/ak)/Stot/(1-1/ak/4);
             bd               = bk*bk*DDSect1*(1-(1+Anucleus/ak*(1-1/ak/4))*
                                    exp(-Anucleus/ak*(1-1/4/ak)))*Rnuc2;
             Dprod        = bd*4*Pi1*Pi1*MbToB;

             TotalCrSec   = Tot0-Dtot;
             InelCrSec    = Inel0-Dprod;
//             InelCrSec1 = Inel2;
             ProdCrSec    = Prod0-Dprod;
             ElasticCrSec = TotalCrSec-InelCrSec;
             QuasyElasticCrSec = InelCrSec-ProdCrSec;

   }

  G4double 
  G4IntegralHadrNucleus::GetElasticCrossSection(const G4DynamicParticle * aHadron,
                                                      G4int               iNucleus)
         {       
                TransformValues(aHadron, iNucleus);
                GetHadronValues(iHadron, HadrEnergy);
                GetIntegralCrSec(Anucleus);
                return(ElasticCrSec);
         }

  G4double 
  G4IntegralHadrNucleus::GetTotalCrossSection(const G4DynamicParticle * aHadron,
                                                    G4int               iNucleus)
        {
                TransformValues(aHadron, iNucleus);
                GetHadronValues(iHadron, HadrEnergy);
                GetIntegralCrSec(Anucleus);
                return(TotalCrSec);
        }

  G4double 
  G4IntegralHadrNucleus::GetInelasticCrossSection(const G4DynamicParticle * aHadron,
                                                        G4int               iNucleus)
       {
                TransformValues(aHadron, iNucleus);
                GetHadronValues(iHadron, HadrEnergy);
                GetIntegralCrSec(Anucleus);
                return(InelCrSec);
       }

  G4double 
  G4IntegralHadrNucleus::GetProductionCrossSection(const G4DynamicParticle * aHadron,
                                                         G4int               iNucleus)
       {
                TransformValues(aHadron, iNucleus);
                GetHadronValues(iHadron, HadrEnergy);
                GetIntegralCrSec(Anucleus);
                return(ProdCrSec);
       }

  G4double 
  G4IntegralHadrNucleus::GetQuasyElasticCrossSection(const G4DynamicParticle * aHadron,
                                                           G4int               iNucleus)
       {
                TransformValues(aHadron, iNucleus);
                GetHadronValues(iHadron, HadrEnergy);
                GetIntegralCrSec(Anucleus);
                return(QuasyElasticCrSec);
       }



  void 
  G4IntegralHadrNucleus::TransformValues(const G4DynamicParticle* aHadron,
                                               G4int              aNucleus)
         {
                Anucleus   = aNucleus;
                TypeHadron = "proton"; //aHadron->GetParticleName();
                mHadr      = aHadron->GetMass();
                HadrEnergy = aHadron->GetTotalEnergy();

            for(iHadron=0;iHadron<6;iHadron++)  {
                   if(HadronNames[iHadron]==TypeHadron) break;
                                                }

//                   G4cout << " Hadron :  " << TypeHadron, iHadron;
         }

   const G4String G4IntegralHadrNucleus::HadronNames[6] 
          = {"proton","anti_proton","pi+","pi-", "kaon+","kaon-"};

//

/* end of file */
