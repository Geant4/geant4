// G4DiffElasticHadrNucleus.cc

#include "G4DiffElasticHadrNucleus.hh"

   G4double   G4DiffElasticHadrNucleus::HadrNuclDifferCrSec(const
                                            G4DynamicParticle* aHadron, 
                                            G4Nucleus*         aNucleus,
                                            G4double           aQ2)
    {
      
      G4double HadrEnergy = aHadron->GetTotalEnergy()/1000;

      G4HadronValues::GetHadronValues(aHadron);

      G4double Q2 = aQ2/1000/1000;

      G4int    Nucleus      = (int)aNucleus->GetN();

        if(Nucleus>42) {
      G4cout << " This nucleus is rather heavy for this model" << G4endl;
                       }

      G4double  R0    = sqrt(0.99);                 // This is in fermi
          if( Nucleus > 10 ) R0 = sqrt (0.84);
          if( Nucleus > 20 ) R0 = sqrt((35.34+0.65*Nucleus)/(40.97+Nucleus));

      G4double    MbToB    = 2.568;      // Change 0.5 -> 0.7 for N =40 !!
      G4double    Pi1      = 3.1416;
      G4double    Stot     = HadrTot*MbToB;          //In GeV-2
      G4double    Bhad     = HadrSlope;              //In GeV-2
      G4double    Asq      = 1+HadrReIm*HadrReIm;
      G4double    Rho2     = 1+Asq;
      G4double    Rnucl    = R0*pow(Nucleus,0.3333); //  In Fermi 
      G4double    Rnuc2    = Rnucl*Rnucl*MbToB*10;      //  In GeV-2
      G4double    R1       = Rnucl*sqrt(MbToB*10);


// G4cout<<" Hadr mass (elastic) "<<aHadron->GetMass()<<G4endl;

      G4double    R2       = 7.5;
      G4double    Pnucl    = 0.8;
      G4double    R12      = R1*R1;
      G4double    R22      = R2*R2;
      G4double    R12B     = R12+2*Bhad;
      G4double    R22B     = R22+2*Bhad;
      G4double    R12Bp    = R12+20;
      G4double    R22Bp    = R22+20;
      G4double    R13Bp    = R12*R1/R12Bp;
      G4double    R23Bp    = R22*R2/R22Bp;
      G4double    R12Ap    = R12+20;
      G4double    R22Ap    = R22+20;
      G4double    R13Ap    = R12*R1/R12Ap;
      G4double    R23Ap    = R22*R2/R22Ap*Pnucl;
      G4double    R23dR13  = R23Ap/R13Ap;
      G4double    R12Apd   = 2/R12Ap;
      G4double    R22Apd   = 2/R22Ap;
      G4double R12ApdR22Ap = 0.5*(R12Apd+R22Apd);

      G4double DDSec1p  = (DDSect2+DDSect3*log(1.06*2*HadrEnergy/R1/
                             sqrt(25.68)/4));
      G4double DDSec2p  = (DDSect2+DDSect3*log(1.06*2*HadrEnergy/
                             sqrt((R12+R22)/2/25.68)/4));
      G4double DDSec3p  = (DDSect2+DDSect3*log(1.06*2*HadrEnergy/R2/
                             sqrt(25.68)/4));

      G4double    Norm     = R12*R1-Pnucl*R22*R2;
      G4double    R13      = R12*R1/R12B;
      G4double    R23      = Pnucl*R22*R2/R22B;
      G4double    Unucl    = Stot/2/Pi1/Norm*R13;
      G4double    UnuclScr= Stot/2/Pi1/Norm*R13Ap;
      G4double    UnucBUF  = Stot/2/Pi1/R12B;
      G4double    N        = -1;
      G4double    N2       = R23/R13;
      G4double    NBUF     = -1/UnucBUF;

      G4double    ImElasticAmpl0  = 0;
      G4double    ReElasticAmpl0  = 0;

          for(G4int i=1; i<=Nucleus; i++)
             {
                  N     = -N*Unucl*(Nucleus-i+1)/i;
      G4double    N4    = 1;
      G4double    Prod1 = exp(-Q2/i*R12/4)/i*R12;

                  for(G4int l=1; l<=i; l++)
                    {
                      G4double exp1  = l/R22B+(i-l)/R12B;
                               N4    = -N4*(i-l+1)/l*N2;
                               Prod1 = Prod1+N4/exp1*exp(-Q2/exp1/4);
                    }  // end l

               G4double          SinFi  = HadrReIm/sqrt(Rho2);
               G4double          FiH    = asin(SinFi);
                        ReElasticAmpl0  = ReElasticAmpl0+Prod1*N*sin(FiH*i);
                        ImElasticAmpl0  = ImElasticAmpl0+Prod1*N*cos(FiH*i);
             if(Prod1*N/ImElasticAmpl0 < 0.001) break;
               }  // i

      ImElasticAmpl0 = ImElasticAmpl0*Pi1/sqrt(Rho2*2.568); 
      ReElasticAmpl0 = ReElasticAmpl0*Pi1/sqrt(Rho2*2.568);  

      G4double N1p  = 1;
      G4double Din1 = 0.5*(R13Ap*R13Ap*exp(-Q2/8*R12Ap)/2*R12Ap/2*DDSec1p-
                      2*R23Ap*R13Ap/2/R12ApdR22Ap*exp(-Q2/4*R12ApdR22Ap)*DDSec2p+
                      R23Ap*R23Ap/2*R22Ap/2*exp(-Q2/8*R22Ap)*DDSec3p); //  at i=0

               for(i = 1; i<= Nucleus-2; i++)
                   {
                              N1p  = -N1p*UnuclScr*(Nucleus-i-1)/i;
                     G4double N2p  = 1;
                     G4double Din2 = 0;

                        for(G4int l = 0; l<=i; l++) 
                          {
                            G4double exp1  = l/R22B+(i-l)/R12B;
                            G4double exp1p = exp1+R12Apd;
                            G4double exp2p = exp1+R12ApdR22Ap;
                            G4double exp3p = exp1+R22Apd;
               Din2 = Din2 + N2p*binom(i,l)*
                     (R13Ap*R13Ap/2/exp1p*exp(-Q2/4/exp1p)*DDSec1p-
                    2*R13Ap*R23Ap/2/exp2p*exp(-Q2/4/exp2p)*DDSec2p+
                      R23Ap*R23Ap/2/exp3p*exp(-Q2/4/exp3p)*DDSec3p);
               N2p  = -N2p*R23dR13;
                           } // l

               Din1 = Din1+Din2*N1p*Mnoj[i]/(i+2)/(i+1);
             if(Din2*N1p/Din1 < 0.001) break;
                    }  //  i

             Din1 = 1*Din1*Nucleus*(Nucleus-1)
                     /2/Pi1/Norm/2/Pi1/Norm/sqrt(2.568);
             Din1 = 0;
       G4double DiffCrSec2 = (ReElasticAmpl0*ReElasticAmpl0+
                             (ImElasticAmpl0+Din1)*(ImElasticAmpl0+Din1))*
                                 2.568/4/Pi1;
//   G4cout<<" R1 "<<R1<<" Tot "<<HadrTot<<" B "<<Bhad<<
//            " Re/Im "<< HadrReIm<<" E "<<HadrEnergy/1000<<
//         " Q2 "<<Q2<<" diffSec "<<DiffCrSec2<< "  (elastic)"<<G4endl;
      
             return DiffCrSec2;
          }   // function

/*  End of file  */

