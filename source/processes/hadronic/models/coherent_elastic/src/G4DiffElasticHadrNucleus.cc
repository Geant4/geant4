//G4DiffElasticHadrNucleus.cc

#include "globals.hh"
#include "G4DiffElasticHadrNucleus.hh"

//  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   void  G4DiffElasticHadrNucleus::
         GetNucleusParameters(G4Nucleus   * aNucleus)
    {
      G4int       Nucleus  = (int)aNucleus->GetN();
                  R1       = 4.1408*pow(Nucleus,0.3018);
                  R2       = 3.806*pow(Nucleus-10.068,0.2685);
                  Pnucl    = 0.95;

     }

//  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   G4double   G4DiffElasticHadrNucleus::HadrNuclDifferCrSec(
                                  const   G4DynamicParticle * aHadron, 
                                          G4Nucleus         * aNucleus,
                                          G4double           aQ2)
    {
//   ------ All external kinematical variables are in MeV -------
//            ------ but internal in GeV !!!!  ------

      G4double HadrEnergy = aHadron->GetTotalEnergy()/1000;  //  GeV

      G4int    Nucleus      = (int)aNucleus->GetN();

 if(Nucleus<11)
         {
  G4Exception(" This nucleus is very light for this model !!!");
         }

  if(Nucleus>208)
         {
  G4Exception(" This nucleus is very heavy for this model !!!");
         }

  if(HadrEnergy < 2)
  G4Exception(" The hadron energy is very low for this model !!!");

      G4HadronValues::GetHadronValues(aHadron);

      GetNucleusParameters(aNucleus);

      G4double Q2 = aQ2/1000/1000;                           //  GeV

      G4double MomentumCMN, S, MassH, MassN, EcmH;

               MassH       = aHadron->GetMass()/1000;
               MassN       = Nucleus*0.938;
               S           = 2*MassN*HadrEnergy+MassN*MassN+MassH*MassH;
               EcmH        = (S-MassN*MassN+MassH*MassH)/2/sqrt(S);
               MomentumCMN = sqrt(EcmH*EcmH-MassH*MassH);

      G4double    MbToB    = 2.568;        //  from mb to GeV^-2
      G4double    Pi1      = 3.1416;
      G4double    Stot     = HadrTot*MbToB;          //In GeV-2
      G4double    Bhad     = HadrSlope;              //In GeV-2
      G4double    Asq      = 1+HadrReIm*HadrReIm;
      G4double    Rho2     = sqrt(Asq);
      G4double    Pnuclp   = 0.001;
                  Pnuclp   = Pnucl;
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
      G4double    R23Ap    = R22*R2/R22Ap*Pnuclp;
      G4double    R23App   = R22*R2/R22Ap*Pnuclp;
      G4double    R23dR13  = R23Ap/R13Ap;
      G4double    R12Apd   = 2/R12Ap;
      G4double    R22Apd   = 2/R22Ap;
      G4double R12ApdR22Ap = 0.5*(R12Apd+R22Apd);

      G4double DDSec1p  = (DDSect2+DDSect3*log(1.06*2*HadrEnergy/R1/4));
      G4double DDSec2p  = (DDSect2+DDSect3*log(1.06*2*HadrEnergy/
                             sqrt((R12+R22)/2)/4));
      G4double DDSec3p  = (DDSect2+DDSect3*log(1.06*2*HadrEnergy/R2/4));

      G4double    Norm     = R12*R1-Pnucl*R22*R2;
      G4double    Normp    = R12*R1-Pnuclp*R22*R2;
      G4double    R13      = R12*R1/R12B;
      G4double    R23      = Pnucl*R22*R2/R22B;
      G4double    Unucl    = Stot/2/Pi1/Norm*R13;
      G4double    UnuclScr = Stot/2/Pi1/Normp*R13Ap;
      G4double    SinFi    = HadrReIm/Rho2;
      G4double    FiH      = asin(SinFi);
      G4double    N        = -1;
      G4double    N2       = R23/R13;

      G4double    ImElasticAmpl0  = 0;
      G4double    ReElasticAmpl0  = 0;

      G4double    exp1;
      G4double    N4;
      G4double    Prod1;
      G4int       i;

        for( i=1; i<=Nucleus; i++)
             {
                  N     = -N*Unucl*(Nucleus-i+1)/i*Rho2;
                  N4    = 1;
                  Prod1 = exp(-Q2/i*R12B/4)/i*R12B;

                  for(G4int l=1; l<=i; l++)
                    {
                         exp1  = l/R22B+(i-l)/R12B;
                         N4    = -N4*(i-l+1)/l*N2;
                         Prod1 = Prod1+N4/exp1*exp(-Q2/exp1/4);
                    }  // end l

                  ReElasticAmpl0  = ReElasticAmpl0+Prod1*N*sin(FiH*i);
                  ImElasticAmpl0  = ImElasticAmpl0+Prod1*N*cos(FiH*i);

             if(abs(Prod1*N/ImElasticAmpl0) < 0.000001) break;
               }      // i

    ImElasticAmpl0 = ImElasticAmpl0*Pi1/2.568;   // The amplitude in mB
    ReElasticAmpl0 = ReElasticAmpl0*Pi1/2.568;   // The amplitude in mB

      G4double N1p  = 1;
      G4double Din1 = 0.5*(R13Ap*R13Ap*exp(-Q2/8*R12Ap)/2*R12Ap/2*DDSec1p-
              2*R23Ap*R13Ap/2/R12ApdR22Ap*exp(-Q2/4/R12ApdR22Ap)*DDSec2p+
              R23Ap*R23Ap/2*R22Ap/2*exp(-Q2/8*R22Ap)*DDSec3p);   // at i=0

      if(Nucleus>95 && HadrEnergy>1000) {Din1=0; goto tam;}

               G4double exp1p;
               G4double exp2p;
               G4double exp3p;
               G4double N2p;
               G4double Din2, BinCoeff = 1;

               for( i = 1; i<= Nucleus-2; i++)
                   {
                     N1p  = -N1p*UnuclScr*(Nucleus-i-1)/i*Rho2;
                     N2p  = 1;
                     Din2 = 0;

                        for(G4int l = 0; l<=i; l++) 
                          {

                          if(l == 0)      BinCoeff = 1;
                          else if(l !=0 ) BinCoeff = BinCoeff*(i-l+1)/l;

                             exp1  = l/R22B+(i-l)/R12B;
                             exp1p = exp1+R12Apd;
                             exp2p = exp1+R12ApdR22Ap;
                             exp3p = exp1+R22Apd;

               Din2 = Din2 + N2p*BinCoeff*
                     (R13Ap*R13Ap/2/exp1p*exp(-Q2/4/exp1p)*DDSec1p-
                    2*R13Ap*R23Ap/2/exp2p*exp(-Q2/4/exp2p)*DDSec2p+
                      R23Ap*R23Ap/2/exp3p*exp(-Q2/4/exp3p)*DDSec3p);

               N2p  = -N2p*R23dR13;
                           }     // l

               Din1 = Din1+Din2*N1p*Mnoj[i]/(i+2)/(i+1)*cos(FiH*i);
             if(abs(Din2*N1p/Din1) < 0.000001) break;
                    }           //  i

                   Din1 = -1*Din1*Nucleus*(Nucleus-1)
                         /2/Pi1/Normp/2/Pi1/Normp*16*Pi1*Pi1;
  
//  ----------------  dSigma/dOmegaCM,  mb/Ster  -----------------

tam:
    G4double DiffCrSec2 = (ReElasticAmpl0*ReElasticAmpl0+
                          (ImElasticAmpl0+Din1)*
                          (ImElasticAmpl0+Din1))*
                          2.568/4/Pi1/Pi1
                          *MomentumCMN*MomentumCMN;

       return DiffCrSec2;  //  dSig/dOmegaCM, mb/Ster
   }   // function

/*  End of file  */
