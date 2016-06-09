//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4DiffElasticHadrNucleus.cc,v 1.15 2005/06/10 13:23:42 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $

//G4DiffElasticHadrNucleus.cc

#include "globals.hh"
#include "G4DiffElasticHadrNucleus.hh"

//  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   void  G4DiffElasticHadrNucleus::
         GetNucleusParameters(G4Nucleus   * aNucleus)
    {
      G4int       Nucleus  = (int)aNucleus->GetN();
  
                  if(Nucleus == 208)
                    {  //  R1 = 20.73; R2 = 15.74.
                  R1       = 4.1408*std::pow(static_cast<double>(Nucleus),0.3018);
                  R2       = 3.806*std::pow(Nucleus-10.068,0.2685);
                  Pnucl    = 0.9;
                  Aeff     = 1.1;
                  R1       = 19.5;
                  R2       = 15.74;
                  Pnucl    = 0.4;
                  Aeff     = 0.7;
                    }

                  else if(Nucleus == 90) 
                    {
                      R1    = 16.50;
                      R2    = 11.62;
                      Pnucl = 0.4;
                      Aeff  = 0.9;
                      R1    = 16.5;
                      R2    = 11.62;
                      Pnucl = 0.4;
                      Aeff  = 0.7;
                    }


                  else if(Nucleus == 58) 
                    {
                      R1    = 15.0;
                      R2    = 9.9;
                      Pnucl = 0.45;
                      Aeff  = 0.85;
                    }
                  else if(Nucleus == 16)
                    {
                      R1    = 10.50;
                      R2    = 5.5;
                      Pnucl = 0.7;
                      Aeff  = 0.98;
//                      R1    = 11.3;
//                      R2    = 2.5;
//                      Pnucl = 0.75;
//                      Aeff  = 0.9;

                    }
                  else
                    {
                      R1    = 4.45*std::pow(static_cast<double>(Nucleus-1),0.309);
                 if(Nucleus == 28)
                      R1    = 4.25*std::pow(static_cast<double>(Nucleus-1),0.309);
                      R2    = 2.3*std::pow(static_cast<double>(Nucleus),0.36);
                      Pnucl = 0.176+0.00167*Nucleus+
                                 8.69E-6*Nucleus*Nucleus;
                      Aeff  = 0.9;
                    }

               if(Nucleus == 9) 
               {
                  R1    = 9.0;
                  R2    = 7.0;
                  Pnucl = 0.190;
                  Aeff  = 0.9;
               }
/*
               if(Nucleus == 12) 
               {
                  R1    = 9.336;
                  R2    = 5.63;
                  Pnucl = 0.197;
                  Aeff  = 01.0;
               }
*/
/*
               if(Nucleus == 11) 
               {
                  R1    = 10.8;
                  R2    = 7.5;
                  Pnucl = 0.85;
                  Aeff  = 1.2;
               }
*/
               if(Nucleus == 4)
               {
                  R1    = 5.5;
                  R2    = 3.7;
                  Pnucl = 0.4;   
                  Aeff  = 0.87;
               }
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

 if(Nucleus<4)
         {
  G4Exception(" This nucleus is very light for this model !!!");
         }

  if(Nucleus>208)
         {
  G4Exception(" This nucleus is very heavy for this model !!!");
         }

  if(HadrEnergy < 1.4999)
  {
  G4cout<<HadrEnergy<<G4endl;
  G4Exception(" The hadron energy is very low for this model: E= ");
  }
      G4HadronValues::GetHadronValues(aHadron);
      GetTotalCrossSection(aHadron, aNucleus);

//   G4cout<<"  Tot00 "<<Tot00<<" DTot00 "<<DTot00<<endl;

      GetNucleusParameters(aNucleus);

      G4double Q2 = aQ2/1000/1000;                           //  GeV

      G4double MomentumCMN, S, MassH, MassN, EcmH;

               MassH       = aHadron->GetMass()/1000;
               MassN       = Nucleus*0.938;
               S           = 2*MassN*HadrEnergy+MassN*MassN+MassH*MassH;
               EcmH        = (S-MassN*MassN+MassH*MassH)/2/std::sqrt(S);
               MomentumCMN = std::sqrt(EcmH*EcmH-MassH*MassH);

      G4double    MbToB    = 2.568;        //  from mb to GeV^-2
      G4double    Pi1      = 3.1416;
      G4double    Stot     = HadrTot*MbToB;          //In GeV-2
      G4double    Bhad     = HadrSlope;              //In GeV-2
      G4double    Asq      = 1+HadrReIm*HadrReIm;
      G4double    Rho2     = std::sqrt(Asq);
      G4double    Pnuclp   = 0.001;
                  Pnuclp   = Pnucl;
      G4double    R12      = R1*R1;
      G4double    R22      = R2*R2;
      G4double    R12B     = R12+2*Bhad;
      G4double    R22B     = R22+2*Bhad;
//      G4double    R12Bp    = R12+20;
//      G4double    R22Bp    = R22+20;
//      G4double    R13Bp    = R12*R1/R12Bp;
//      G4double    R23Bp    = R22*R2/R22Bp;
      G4double    R12Ap    = R12+20;
      G4double    R22Ap    = R22+20;
      G4double    R13Ap    = R12*R1/R12Ap;
      G4double    R23Ap    = R22*R2/R22Ap*Pnuclp;
//      G4double    R23App   = R22*R2/R22Ap*Pnuclp;
      G4double    R23dR13  = R23Ap/R13Ap;
      G4double    R12Apd   = 2/R12Ap;
      G4double    R22Apd   = 2/R22Ap;
      G4double R12ApdR22Ap = 0.5*(R12Apd+R22Apd);

      G4double DDSec1p  = (DDSect2+DDSect3*std::log(1.06*2*HadrEnergy/R1/4));
      G4double DDSec2p  = (DDSect2+DDSect3*std::log(1.06*2*HadrEnergy/
                             std::sqrt((R12+R22)/2)/4));
      G4double DDSec3p  = (DDSect2+DDSect3*std::log(1.06*2*HadrEnergy/R2/4));

      G4double    Norm     = (R12*R1-Pnucl*R22*R2)*Aeff;
      G4double    Normp    = (R12*R1-Pnuclp*R22*R2)*Aeff;
      G4double    R13      = R12*R1/R12B;
      G4double    R23      = Pnucl*R22*R2/R22B;
      G4double    Unucl    = Stot/2/Pi1/Norm*R13;
      G4double    UnuclScr = Stot/2/Pi1/Normp*R13Ap;
      G4double    SinFi    = HadrReIm/Rho2;
      G4double    FiH      = std::asin(SinFi);
      G4double    N        = -1;
      G4double    N2       = R23/R13;

      G4double    ImElasticAmpl0  = 0;
      G4double    ReElasticAmpl0  = 0;

      G4double    exp1;
      G4double    N4;
      G4double    Prod1, Tot1=0, medTot, DTot1, DmedTot;
      G4int       i;

        for( i=1; i<=Nucleus; i++)
             {
                  N       = -N*Unucl*(Nucleus-i+1)/i*Rho2;
                  N4      = 1;
                  Prod1   = std::exp(-Q2/i*R12B/4)/i*R12B;
                  medTot  = R12B/i;

                  for(G4int l=1; l<=i; l++)
                    {
                     exp1    = l/R22B+(i-l)/R12B;
                     N4      = -N4*(i-l+1)/l*N2;
                     Prod1   = Prod1+N4/exp1*std::exp(-Q2/exp1/4);
                     medTot  = medTot+N4/exp1;
                    }  // end l

                  ReElasticAmpl0  = ReElasticAmpl0+Prod1*N*std::sin(FiH*i);
                  ImElasticAmpl0  = ImElasticAmpl0+Prod1*N*std::cos(FiH*i);
                  Tot1            = Tot1+medTot*N*std::cos(FiH*i);
             if(std::fabs(Prod1*N/ImElasticAmpl0) < 0.000001) break;
               }      // i

    ImElasticAmpl0 = ImElasticAmpl0*Pi1/2.568;   // The amplitude in mB
    ReElasticAmpl0 = ReElasticAmpl0*Pi1/2.568;   // The amplitude in mB
    Tot1           = Tot1*Pi1*2.0/2.568;

      G4double N1p  = 1;
      G4double Din1 = 0.5*(R13Ap*R13Ap*std::exp(-Q2/8*R12Ap)/2*R12Ap/2*DDSec1p-
              2*R23Ap*R13Ap/2/R12ApdR22Ap*std::exp(-Q2/4/R12ApdR22Ap)*DDSec2p+
              R23Ap*R23Ap/2*R22Ap/2*std::exp(-Q2/8*R22Ap)*DDSec3p);   // at i=0

           DTot1 = 0.5*(R13Ap*R13Ap/2*R12Ap/2*DDSec1p-
                       2*R23Ap*R13Ap/2/R12ApdR22Ap*DDSec2p+
                       R23Ap*R23Ap/2*R22Ap/2*DDSec3p);   // at i=0

//    G4cout<<" DTot00 "<<DTot00<<" Din1 "<<Din1<<" DTot1 "<<DTot1<<endl;

/////      if(Nucleus>95 && HadrEnergy>1000) {Din1=0; goto tam;}

               G4double exp1p;
               G4double exp2p;
               G4double exp3p;
               G4double N2p;
               G4double Din2, BinCoeff;

               BinCoeff = 1;

               for( i = 1; i<= Nucleus-2; i++)
               {
                 N1p     = -N1p*UnuclScr*(Nucleus-i-1)/i*Rho2;
                 N2p     = 1;
                 Din2    = 0;
                 DmedTot = 0;
                   for(G4int l = 0; l<=i; l++) 
                   {
                    if(l == 0)      BinCoeff = 1;
                     else if(l !=0 ) BinCoeff = BinCoeff*(i-l+1)/l;

                     exp1  = l/R22B+(i-l)/R12B;
                     exp1p = exp1+R12Apd;
                     exp2p = exp1+R12ApdR22Ap;
                     exp3p = exp1+R22Apd;

                Din2     = Din2 + N2p*BinCoeff*
                     (R13Ap*R13Ap/2/exp1p*std::exp(-Q2/4/exp1p)*DDSec1p-
                    2*R13Ap*R23Ap/2/exp2p*std::exp(-Q2/4/exp2p)*DDSec2p+
                      R23Ap*R23Ap/2/exp3p*std::exp(-Q2/4/exp3p)*DDSec3p);

                DmedTot = DmedTot + N2p*BinCoeff*
                           (R13Ap*R13Ap/2/exp1p*DDSec1p-
                          2*R13Ap*R23Ap/2/exp2p*DDSec2p+
                            R23Ap*R23Ap/2/exp3p*DDSec3p);
            
                N2p   = -N2p*R23dR13;
                           }     // l

               Din1  = Din1+Din2*N1p*Mnoj[i]/(i+2)/(i+1)*std::cos(FiH*i);
               DTot1 = DTot1+DmedTot*N1p*Mnoj[i]/(i+2)/(i+1)*std::cos(FiH*i);
             if(std::fabs(Din2*N1p/Din1) < 0.000001) break;
                    }           //  i

                   Din1 = -1*Din1*Nucleus*(Nucleus-1)
                         /2/Pi1/Normp/2/Pi1/Normp*16*Pi1*Pi1;

                   DTot1 = 1*DTot1*Nucleus*(Nucleus-1)
                         /2/Pi1/Normp/2/Pi1/Normp*16*Pi1*Pi1;

///      Din1 *= 0;
//  ----------------  dSigma/dOmegaCM,  mb/Ster  -----------------

// tam:
   G4double Corr0=Tot00/Tot1/1.0; // Corr1=DTot00/DTot1/1.0;

             ImElasticAmpl0 *= Corr0;
//             Din1           *= Corr1;

    G4double DiffCrSec2 = (ReElasticAmpl0*ReElasticAmpl0+
                          (ImElasticAmpl0+Din1)*
                          (ImElasticAmpl0+Din1))*
                          2.568/4/Pi1/Pi1
                          *MomentumCMN*MomentumCMN;

    AIm = ImElasticAmpl0;
    ARe = ReElasticAmpl0;
    DIm = Din1;

// if(Q2<0.001)
 {
// G4cout<<" DTot00 "<<DTot00<<" Din1 "<<Din1<<" DTot1 "<<DTot1
//       <<" Tot00 "<<Tot00<<" Tot1 "<<Tot1
//       <<" Ampl "<<ImElasticAmpl0<<" dSdT "<<CrSecT<<G4endl;
 }

     return DiffCrSec2;  //  dSig/dOmegaCM, mb/Ster
   }   // function
//  ##############################################
   G4double G4DiffElasticHadrNucleus::
               Thickness(G4int  A,  G4double b)
   {
       G4double An=A, Dn, Ct, Tn, Bn, En;
       G4int    kk;
       G4double dr, r, WSo, WS, SumZ=0, SumN;

       dr = rAmax*rAfm/(NpointsB-1);

       G4double  Norm = 
       3.0/4.0/3.1416/rAfm/rAfm/rAfm/(1+0.53*0.53*3.14*3.14/rAfm/rAfm);

//       Norm = 1/(SqrtPi*Rlight*Rlight*Rlight*3.1416)/
//                (1+1.5*alpha*Rlight*Rlight);

       SumN = 0.0;
       for(kk=0; kk<NpointsB; kk++)
       {
         r     = kk*dr;
         WS    = 1/(1+std::exp((std::sqrt(b*b+r*r)-rAfm)/0.53));
         WSo   = 1/(1+std::exp((r-rAfm)/0.53));

//         WS     = (1.0+alpha*(b*b+r*r))*
//	           std::exp(-(b*b+r*r)/Rlight/Rlight);
//         WSo    = (1.0+alpha*r*r)*
//                  std::exp(-r*r/Rlight/Rlight);

	 SumZ += WS*dr*Norm*A*2.0;
	 SumN += 4*3.1416*r*r*WSo*dr*Norm;
// G4cout<<" ThickI  r  "<<r<<"  WS  "<<WS<<"  Wo  "<<WSo
//       <<"  SumZ  "<<SumZ<<"  SumN  "<<SumN<<G4endl;
       }

       Dn = 0.42*std::pow(An,-0.26);  //  fm^-2
       Bn = 8.0e-4*An*An;
       Ct = An*Dn/3.1416/std::log(1+Bn);
       En = std::exp(-Dn*b*b);
       Tn = Ct*Bn*En/(1+Bn*En);  //  fm^-2
/*
 G4cout<<"  Th: b "<<"  "<<b<<"  Sz  "
       <<SumZ
       <<"   SumN  "<<SumN<<"  Norm  "<<Norm
       <<"   Tn  "<<Tn<<G4endl;
*/
     return SumZ;
//     return Tn;
   }
//  ##############################################
    G4double G4DiffElasticHadrNucleus::
               GetIntegrandS(G4int    Anucleus,
                             G4double ImpactPar,
                             G4int    Kind)
   {
     G4int     iInt;
     G4double ValS, Integ=0;
     G4double FunS, ExpS, I0bs;

     for(iInt=1; iInt<NpointsB; iInt++)
     {
       ValS  = iInt*stepB;     //  GeV^(-1)
       FunS = ImpactPar*ValS/HadrSlope;  

       if(Kind == 0)  Thick[iInt] = 
              Thickness(Anucleus, ValS/std::sqrt(25.68))/25.68;

       I0bs = 0.0;
       if(FunS > 320) continue;

       ExpS  = std::exp(-(ValS*ValS+ImpactPar*ImpactPar)/2.0/HadrSlope);
       I0bs  = MyI0(FunS);

       FunS = ValS*ExpS*I0bs*Thick[iInt];
       Integ += FunS;
/*
  if(ImpactPar*ValS/Slope>200 && ImpactPar>59 && ImpactPar<60 )
  G4cout<<" Int. S: S  "<<ValS<<"  Exp "<<ExpS
        <<" I0  "<<I0bs<<"   Th(fm)  "<<Thick[iInt]*25.68
        <<"   F  "<<FunS<<"  Int  "<<Integ<<G4endl;
*/
     } 
     return Integ*stepB;
   }
//  #############################################
   void   G4DiffElasticHadrNucleus::
               GetIntegrandB(G4int  Anucleus)
   {
     G4double ValB;
     G4int    iInt;
     G4double expB, IntegS, InExp;
     G4double Sigm=HadrTot*2.568;

// G4cout<<G4endl<<" Nucleus  (r0,r01):"<<r0<<"  "<<r01
//       <<"  "<<Anucleus<<"  r  "<<rAfm
//       <<" fm ("<<rAGeV<<" GeV^1)"<<G4endl<<G4endl;

     for(iInt=0; iInt<NpointsB; iInt++)
     {
       ValB   = iInt*stepB;
       IntegS = GetIntegrandS(Anucleus, ValB, iInt);
       InExp  = -Sigm/2/HadrSlope*IntegS;
       expB   = std::exp(InExp);

       ReIntegrand[iInt]  = (1.0-expB*std::cos(HadrReIm*InExp));
       ImIntegrand[iInt]  =     expB*std::sin(HadrReIm*InExp);

// G4cout<<" b  Int(b) :  "<<ValB<<"  InExp "<<InExp<<"  "
//       <<ReIntegrand[iInt]<<"  1-Re(B)  "<<ReIntegrand[iInt]
//       <<"  Im(B)  "<<ImIntegrand[iInt]<<G4endl;
     } 
   }
//  #############################################
   G4double G4DiffElasticHadrNucleus::
  HadrNuclDifferCrSecT(const  G4DynamicParticle *  aHadron,
                              G4Nucleus         *  aNucleus, 
                              G4double             aQ2,
                              G4int                Kind)
   {
     G4double J0qb, Re;
     G4double HadrEnergy = aHadron->GetTotalEnergy()/1000;  //  GeV
     G4int    Nucleus    = (int)aNucleus->GetN();

  if(Nucleus<4)
    G4Exception(" This nucleus is very light for this model !!!");

  if(Nucleus>208)
    G4Exception(" This nucleus is very heavy for this model !!!");

//  G4cout<<" Energy (T) "<<HadrEnergy<<G4endl;
  if(HadrEnergy < 1.4)
//    G4cout<<" The hadron energy is very low for this model !!!"<<G4endl;
    G4Exception(" The hadron energy is very low for this model !!!");

      G4HadronValues::GetHadronValues(aHadron);

      G4double Q2 = aQ2/1000/1000;             //  GeV

      G4double S, MassH, MassN, EcmH;

      switch(Nucleus)
      {
      case(208) :   r0      = 1.125;
                    r01     = 1.16*0.0;
                    break;
      case(90) :    r0       = 1.12;
                    r01      = 1.16*0.0;
                    break;
      case(58) :    r0       = 1.09;
                    r01      = 1.16*0.0;
                    break;
      case(48) :    r0       = 1.07;
                    r01      = 1.16*0.0;
                    break;

      case(40) :    r0       = 1.15;
                    r01      = 1.16*0.0;
                    break;
      case(28) :    r0       = 0.93;
                    r01      = 1.16*0.0;
                    break;
      case(16) :    r0       = 0.92;
                    r01      = 1.16*0.0;
                    break;
      case(12) :    r0       = 0.8;
                    r01      = 1.16*0.0;
                    break;
      case(64) :    r0       = 1.1;
                    r01      = 1.16*0.0;
                    break;
      default  :    Re = std::pow(static_cast<double>(Nucleus), 0.3333);
                    r0 = 1.16*(1.0-1.16/Re/Re);
        }

               MassH       = aHadron->GetMass()/1000;
               MassN       = Nucleus*0.938;
               S           = 2*MassN*HadrEnergy+MassN*MassN+MassH*MassH;
               EcmH        = (S-MassN*MassN+MassH*MassH)/2/std::sqrt(S);
               MomentumCMN = std::sqrt(EcmH*EcmH-MassH*MassH);

      G4int kk;      
      G4double ValB;
      G4double ReIntSumm=0.0, ImIntSumm=0.0, Section;

      if(Kind == 0)
      {
        MassH       = aHadron->GetMass()/1000;
        MassN       = Nucleus*0.938;
        S           = 2*MassN*HadrEnergy+MassN*MassN+MassH*MassH;
        EcmH        = (S-MassN*MassN+MassH*MassH)/2/std::sqrt(S);
        MomentumCMN = std::sqrt(EcmH*EcmH-MassH*MassH);

        rAfm        = r0*std::pow(static_cast<double>(Nucleus), 0.3333)*
                      (1-r01/std::pow(static_cast<double>(Nucleus),0.666));
        rAGeV       = rAfm*std::sqrt(25.68);
        stepB       = rAmax*rAGeV/(NpointsB-1);

        GetIntegrandB(Nucleus);

// G4cout<<G4endl<<" Nucleus  "<<Nucleus<<"  r  "<<rAfm
//       <<" fm ("<<rAGeV<<" GeV^1)  "<<Kind<<G4endl<<G4endl;
      }

     if(Kind==0) InCohI = 0.0;
     for(kk=0; kk<NpointsB; kk++)
      {
        ValB = stepB*kk;

        if(Kind==0) InCohI += Thick[kk]*ValB*
		      std::exp(-HadrTot*2.568*Thick[kk])
                       *2.0*3.1416*stepB
                       /16/3.1416*std::pow(HadrTot*2.568,2)*
                       (1+std::pow(HadrReIm,2))/2.56;

        J0qb = MyJ0(std::sqrt(Q2)*ValB)*ValB;
        ReIntSumm += J0qb*ReIntegrand[kk];
        ImIntSumm += J0qb*ImIntegrand[kk];
/*
 G4cout<<" Q2  "<<Q2<<" B  "<<ValB<<" Re   "
       <<ReIntegrand[kk]<<" Im  "<<ImIntegrand[kk]
       <<"   J0qb  "<<J0qb<<"  ReF "<<ReIntSumm
       <<"  ImF  "<<ImIntSumm<<G4endl;
*/
      }

     InCoh  = InCohI*
//   2.0*3.1416*stepB
// /16/3.1416*std::pow(HadrTot*2.568,2)*
//      (1+std::pow(HadrReIm,2))/2.568*
     std::exp(-HadrSlope*Q2);

     Section = (ReIntSumm*ReIntSumm+ImIntSumm*ImIntSumm)
       /2.0/3.1416*2.568*MomentumCMN*MomentumCMN 
                *stepB*stepB;

//     Section = (ReIntSumm*ReIntSumm+ImIntSumm*ImIntSumm)
//                /2.0/Pi1*2.568*MomentumCMN*MomentumCMN;

   return Section;
   }  

/*  End of file  */
