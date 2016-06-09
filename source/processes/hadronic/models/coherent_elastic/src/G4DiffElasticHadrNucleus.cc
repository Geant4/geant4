//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4DiffElasticHadrNucleus.cc,v 1.23 2006/06/29 20:09:23 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
//  Geant4 class G4DiffElasticHadrNucleus  
//
//  High energy hadron-nucleus elastic scattering
//  Kinetic energy T > 1 GeV
//  N.  Starkov 2003.
//
//  14.11.05 The HE elastic scattering on proton is added (N.Starkov)
//  30.05.06 New version, which not uses elastic data (N.Starkov)
//  30.05.06 cleanup (int -> G4int) (V.Ivanchenko)
//

#include "globals.hh"
#include "G4DiffElasticHadrNucleus.hh"

//  ##########################################################
void NucleusParameters::GetNucleusParameters(G4int Nucleus)
  {
    
    if(Nucleus == 208)
    {  
//      R1 = 20.73; R2 = 15.74.
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
//      R1    = 11.3;
//      R2    = 2.5;
//      Pnucl = 0.75;
//      Aeff  = 0.9;
    }
    else
    {
      if(Nucleus == 28)
	R1    = 4.25*std::pow(static_cast<double>(Nucleus-1),0.309);
      else
	R1    = 4.45*std::pow(static_cast<double>(Nucleus-1),0.309);
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

//   G4cout<<" Nucl.Par. "<<Nucleus<<"  R1  "<<R1<<G4endl;
  }

//  ##########################################################
void  G4DiffElasticHadrNucleus::
         GetNucleusParameters(G4Nucleus   * aNucleus)
  {
    G4int Nucleus  = (G4int)aNucleus->GetN();
  
//    G4cout<<" Nucl.Par. "<<Nucleus<<G4endl;

    if(Nucleus == 208)
    {  
//      R1 = 20.73; R2 = 15.74.
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
//      R1    = 11.3;
//      R2    = 2.5;
//      Pnucl = 0.75;
//      Aeff  = 0.9;
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

//   G4cout<<" Nucl.Par. "<<Nucleus<<"  R1  "<<R1<<G4endl;
 }

//  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   G4double G4DiffElasticHadrNucleus::
                HadrNuclDifferCrSec(
                            const   G4DynamicParticle * aHadron, 
                                    G4Nucleus         * aNucleus,
                                    G4double           aQ2)
{
//   ------ All external kinematical variables are in MeV -------
//            ------ but internal in GeV !!!!  ------

  G4double energy = aHadron->GetTotalEnergy()/GeV;  //  GeV
  G4double MassH  = aHadron->GetMass()/GeV;

  G4int    Nucleus      = (G4int)aNucleus->GetN();

  if(Nucleus==2 || Nucleus == 3)
    G4Exception(" This model does not work for nuclei with A=3 0r A= 4");

  if(Nucleus>238)
    G4Exception(" This nucleus is very heavy for this model !!!");

  if(energy-MassH < 1.0)
    G4Exception(" The hadron kinetic energy is low for this model (< 1 GeV");

    G4HadronValues::GetHadronValues(aHadron);
    GetTotalCrossSection(aHadron, aNucleus);

    GetNucleusParameters(aNucleus);

    G4double theQ2 = aQ2/GeV/GeV;                           //  GeV

    G4double MassN = 
      aNucleus->AtomicMass(aNucleus->GetN(),aNucleus->GetZ())/GeV;
    G4double S     = 2*MassN*energy+MassN*MassN+MassH*MassH;
    G4double ecmH  = (S-MassN*MassN+MassH*MassH)/2/std::sqrt(S);
    G4double momentumCMN = std::sqrt(ecmH*ecmH-MassH*MassH);

    G4double    MbToB    = 2.568;        //  from mb to GeV^-2
    G4double    Pi1      = pi;
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

    G4double DDSec1p  = (DDSect2+DDSect3*std::log(1.06*2*energy/R1/4));
    G4double DDSec2p  = (DDSect2+DDSect3*std::log(1.06*2*energy/
                             std::sqrt((R12+R22)/2)/4));
    G4double DDSec3p  = (DDSect2+DDSect3*std::log(1.06*2*energy/R2/4));

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
      Prod1   = std::exp(-theQ2/i*R12B/4)/i*R12B;
      medTot  = R12B/i;

       for(G4int l=1; l<=i; l++)
       {
         exp1    = l/R22B+(i-l)/R12B;
         N4      = -N4*(i-l+1)/l*N2;
         Prod1   = Prod1+N4/exp1*std::exp(-theQ2/exp1/4);
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
    G4double Din1 = 0.5*(R13Ap*R13Ap*std::exp(-theQ2/8*R12Ap)/2*R12Ap/2*DDSec1p-
           2*R23Ap*R13Ap/2/R12ApdR22Ap*std::exp(-theQ2/4/R12ApdR22Ap)*DDSec2p+
           R23Ap*R23Ap/2*R22Ap/2*std::exp(-theQ2/8*R22Ap)*DDSec3p);   // at i=0

    DTot1 = 0.5*(R13Ap*R13Ap/2*R12Ap/2*DDSec1p-
                   2*R23Ap*R13Ap/2/R12ApdR22Ap*DDSec2p+
                     R23Ap*R23Ap/2*R22Ap/2*DDSec3p);   // at i=0

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

          Din2  = Din2 + N2p*BinCoeff*
                  (R13Ap*R13Ap/2/exp1p*std::exp(-theQ2/4/exp1p)*DDSec1p-
                   2*R13Ap*R23Ap/2/exp2p*std::exp(-theQ2/4/exp2p)*DDSec2p+
                   R23Ap*R23Ap/2/exp3p*std::exp(-theQ2/4/exp3p)*DDSec3p);

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
                          *momentumCMN*momentumCMN;

     AIm = ImElasticAmpl0;
     ARe = ReElasticAmpl0;
     DIm = Din1;

     return DiffCrSec2;  //  dSig/dOmegaCM, mb/Ster
   }   // function
//  ##############################################
   G4double G4DiffElasticHadrNucleus::
               Thickness(G4int  A,  G4double b)
  {
    G4double An=A, Dn, Ct, Tn, Bn, En, A1, A2;
    G4int    kk;
    G4double dr, r, WSo, WS, SumZ=0, SumN;

    dr = rAmax*rAfm/(NpointsB-1);

    G4double  Norm = 3.0/4.0/3.1416/rAfm/rAfm/rAfm/
                      (1+0.53*0.53*3.14*3.14/rAfm/rAfm);

//       Norm = 1/(SqrtPi*Rlight*Rlight*Rlight*3.1416)/
//                (1+1.5*alpha*Rlight*Rlight);

    A1   = 2.0*A*dr*Norm;
    A2   = 4.0*3.1416*dr*Norm;

    SumN = 0.0;
    for(kk=0; kk<NpointsB; kk++)
    {
      r     = kk*dr;
      WS    = 1/(1+std::exp((std::sqrt(b*b+r*r)-rAfm)/0.53));
      WSo   = 1/(1+std::exp((r-rAfm)/0.53));

//      WS     = (1.0+alpha*(b*b+r*r))*
//	         std::exp(-(b*b+r*r)/Rlight/Rlight);
//      WSo    = (1.0+alpha*r*r)*
//                  std::exp(-r*r/Rlight/Rlight);

//      SumZ += WS*dr*Norm*A*2.0;
//      SumN += 4*3.1416*r*r*WSo*dr*Norm;
      SumZ += WS;
      SumN += r*r*WSo;
     }

     SumZ *= A1;
     SumN *= A2;

     Dn = 0.42*std::pow(An,-0.26);  //  fm^-2
     Bn = 8.0e-4*An*An;
     Ct = An*Dn/pi/std::log(1+Bn);
     En = std::exp(-Dn*b*b);
     Tn = Ct*Bn*En/(1+Bn*En);  //  fm^-2

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
              Thickness(Anucleus, ValS/FmToGeV)/Fm2ToGeV2;

       I0bs = 0.0;
       if(FunS > 320) continue;

       ExpS  = std::exp(-(ValS*ValS+ImpactPar*ImpactPar)/2.0/HadrSlope);
       I0bs  = MyI0(FunS);

       FunS = ValS*ExpS*I0bs*Thick[iInt];
       Integ += FunS;
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
     G4double Sigm=HadrTot/2/HadrSlope*MbToGeV2;

     for(iInt=0; iInt<NpointsB; iInt++)
     {
       ValB   = iInt*stepB;
       IntegS = GetIntegrandS(Anucleus, ValB, iInt);
       InExp  = -Sigm*IntegS;
       expB   = std::exp(InExp);

       ReIntegrand[iInt]  = (1.0-expB*std::cos(HadrReIm*InExp));
       ImIntegrand[iInt]  =     expB*std::sin(HadrReIm*InExp);
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
  G4double energy     = aHadron->GetTotalEnergy()/GeV;  //  GeV
  G4int    Nucleus    = (G4int)aNucleus->GetN();
  G4double MassH      = aHadron->GetMass()/GeV;

  if(Nucleus==2 || Nucleus==3)
    G4Exception(" This nucleus is very light for this model !!!");

  if(Nucleus>238)
    G4Exception(" This nucleus is very heavy for this model !!!");

  if(energy-MassH < 1.0)
    G4Exception(" The hadron energy is very low for this model !!!");

  G4HadronValues::GetHadronValues(aHadron);

  G4double theQ2 = aQ2/GeV/GeV;             //  GeV

  G4int kk;      
  G4double ValB;
  G4double ReIntSumm=0.0, ImIntSumm=0.0;
  G4double Mnojit = 0.0;

  if(Kind == 0)
    {
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

        Mnojit      =    stepB/8.*
                         std::pow(HadrTot*2.568,2)*
                         (1+std::pow(HadrReIm,2))/2.56;


        rAfm        = r0*std::pow(static_cast<double>(Nucleus), 0.3333)*
                      (1-r01/std::pow(static_cast<double>(Nucleus),0.666));
        rAGeV       = rAfm*std::sqrt(25.68);
        stepB       = rAmax*rAGeV/(NpointsB-1);

        GetIntegrandB(Nucleus);
        InCohI = 0.0;
    }

  G4double Simps = 4;

  for(kk=0; kk<NpointsB; kk++)
    {
      ValB = stepB*kk;

      if(Kind==0) InCohI += Thick[kk]*ValB*Simps*Mnojit*
		    std::exp(-HadrTot*2.568*Thick[kk]);


      J0qb = MyJ0(std::sqrt(theQ2)*ValB)*ValB;
      ReIntSumm += J0qb*ReIntegrand[kk];
      ImIntSumm += J0qb*ImIntegrand[kk];
      if(Simps>3) Simps = 2;
    }

  InCoh  = InCohI*
//   2.0*3.1416*stepB
// /16/3.1416*std::pow(HadrTot*2.568,2)*
//      (1+std::pow(HadrReIm,2))/2.568*
    std::exp(-HadrSlope*theQ2)/3.;

  G4double MassN = 
    aNucleus->AtomicMass(aNucleus->GetN(),aNucleus->GetZ())/GeV;
  G4double S     = 2*MassN*energy+MassN*MassN+MassH*MassH;
  G4double ecmH  = (S-MassN*MassN+MassH*MassH)/2/std::sqrt(S);
  G4double momentumCMN = std::sqrt(ecmH*ecmH-MassH*MassH);

  G4double Section = (ReIntSumm*ReIntSumm+ImIntSumm*ImIntSumm)
    /twopi*2.568*momentumCMN*momentumCMN*stepB*stepB/3.;

//     Section = (ReIntSumm*ReIntSumm+ImIntSumm*ImIntSumm)
//                /2.0/Pi1*2.568*MomentumCMN*MomentumCMN;

  return Section;
}  
//  =====================================================
   void  G4DiffElasticHadrNucleus::
       GetKinematics(const G4DynamicParticle * aHadron)
   {
     HadronName = aHadron->GetDefinition()->GetParticleName();

     GetHadronValues(aHadron);

     SigTot = HadrTot;
     ReOnIm = HadrReIm;

     if(Slope==0) Slope  = HadrSlope;

     IntConst = (1-Coeff1)/Slope;

     ProtonM   = proton_mass_c2/GeV;
     HadronE   = aHadron->GetTotalEnergy()/GeV; // GeV
     HdrE      = HadronE;
     MomentumH = aHadron->GetTotalMomentum()/GeV; // GeV
     HadronM   = aHadron->GetMass()/GeV;      // GeV

     PM2    = ProtonM*ProtonM;
     HM2    = HadronM*HadronM;
     Sh     = 2.0*ProtonM*HadronE+PM2+HM2;    // GeV
     SqrtS  = std::sqrt(Sh);
     ConstU = 2*PM2+2*HM2-Sh;

     EcmH   = (Sh+HM2-PM2)/2/SqrtS;
     EcmP   = (Sh-HM2+PM2)/2/SqrtS;

     Kcm    = std::sqrt(EcmH*EcmH-HM2);
     MaxT   = 4*Kcm*Kcm;

     BoundaryP[0]=9.0; BoundaryTG[0]=5.0;BoundaryTL[0]=MaxT/2.0;
     BoundaryP[1]=20.0;BoundaryTG[1]=1.5;BoundaryTL[1]=MaxT;
     BoundaryP[2]=4.0; BoundaryTG[2]=1.0;BoundaryTL[2]=1.5;
     BoundaryP[3]=4.0; BoundaryTG[3]=3.0;BoundaryTL[3]=MaxT;
     BoundaryP[4]=4.0; BoundaryTG[4]=3.0;BoundaryTL[4]=MaxT;
     BoundaryP[5]=5.0; BoundaryTG[5]=2.0;BoundaryTL[5]=MaxT;
     BoundaryP[6]=4.0; BoundaryTG[6]=1.5;BoundaryTL[6]=3.0;
     
     HadrCodes[0] =  2212; HadrCodes[1] =  2112;
     HadrCodes[2] = -2212; HadrCodes[3] =  211;
     HadrCodes[4] = -211;  HadrCodes[5] =  321;
     HadrCodes[6] = -321;

     HadrCode = aHadron->GetDefinition()->GetPDGEncoding();

     G4int NumberH=0;     
   
     while(HadrCode!=HadrCodes[NumberH]) NumberH++;

     if(MomentumH<BoundaryP[NumberH]) MaxTR = BoundaryTL[NumberH];
     else MaxTR = BoundaryTG[NumberH];

     GetParametersHP(aHadron);   
   }
//  +++++++++++++++++++++++++++++++++++++++
  G4double  G4DiffElasticHadrNucleus::
         HadronProtonDiffCrSec(G4double aQ2)
  {
    G4double dSigPodT;

    dSigPodT = SigTot*SigTot*(1+ReOnIm*ReOnIm)*(
                  Coeff1*std::exp(-Slope1*std::sqrt(aQ2))+
                  Coeff2*std::exp( Slope2*(ConstU+aQ2))+
                  (1-Coeff1-Coeff0)*std::exp(-Slope*aQ2)+
                 +Coeff0*std::exp(-Slope0*aQ2)
//                +0.1*(1-std::fabs(CosTh))
                  )/16/pi*2.568;

    return dSigPodT;
  }
//  ++++++++++++++++++++++++++++++++++++++++
  void  G4DiffElasticHadrNucleus::   
       GetParametersHP(const G4DynamicParticle * aHadron)
  {
    G4double HadrP;
    G4int  iStep;

//    HadrCode = aHadron->GetDefinition()->GetPDGEncoding();
    HadrP    = aHadron->GetTotalMomentum()/GeV;

    switch(HadrCode)
    {
      case 2212 :
      case 3122 : 
      case 3222 :
      case 3112 : 
      case 3212 :
      case 3312 :
      case 3322 : 
      case 3334 :
        {
	  G4double EnP[6]={1.5,3.0,5.0,9.0,14.0,19.0};
          G4double C0P[6]={0.15,0.02,0.06,0.08,0.0003,0.0002};
          G4double C1P[6]={0.05,0.02,0.03,0.025,0.0,0.0};
          G4double B0P[6]={1.5,2.5,3.0,4.5,1.4,1.25};
          G4double B1P[6]={5.0,1.0,3.5,4.0,4.8,4.8};
          
          if(BoundaryP[0]<HadrP) 
          {
            Coeff0=Coeff1=Coeff2=0.0;
            Slope0=Slope1=Slope2=0.0;
          }
          else
          {
            iStep = 0;  while(HadrP>EnP[iStep])  iStep++;

            Coeff0 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  C0P[iStep], C0P[iStep-1],  
				  HadrP);            

            Coeff1 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  C1P[iStep], C1P[iStep-1],  
				  HadrP);

            Slope0 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  B0P[iStep], B0P[iStep-1],  
				  HadrP);

            Slope1 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  B1P[iStep], B1P[iStep-1],  
				  HadrP);
          }       //  else
        } break;  //  case 2212

      case  2112 : 
      case -2112 :
        {
	  G4double EnN[5]={1.5,5.0,10.0,14.0,20.0};
          G4double C0N[5]={0.0,0.0,0.02,0.02,0.01};
          G4double C1N[5]={0.06,0.008,0.0015,0.001,0.0003};
          G4double B0N[5]={1.5,2.5,3.8,3.8,3.5};
          G4double B1N[5]={1.5,2.2,3.6,4.5,4.8};
          
          if(BoundaryP[1]<HadrP) 
          {
            Coeff0=Coeff1=0.0;
            Slope0=Slope1=0.0;
          }
          else
          {
            iStep = 0;  while(HadrP>EnN[iStep])   iStep++;

            Coeff0 = LineInterpol(EnN[iStep], EnN[iStep-1],
                                  C0N[iStep], C0N[iStep-1],  
				  HadrP);            

            Coeff1 = LineInterpol(EnN[iStep], EnN[iStep-1],
                                  C1N[iStep], C1N[iStep-1],  
				  HadrP);

            Slope0 = LineInterpol(EnN[iStep], EnN[iStep-1],
                                  B0N[iStep], B0N[iStep-1],  
				  HadrP);

            Slope1 = LineInterpol(EnN[iStep], EnN[iStep-1],
                                  B1N[iStep], B1N[iStep-1],  
				  HadrP);
            Coeff2 = 0.8/HadrP/HadrP;
            Slope2 = 5;
          }       //  else
        } break;  //  case 2112

      case -2212 : 
      case -3122 : 
      case -3222 :
      case -3112 : 
      case -3212 :
      case -3312 :
      case -3322 : 
      case -3334 :
        {
	  G4double EnP[2]={1.5,4.0};
          G4double C0P[2]={0.001,0.0005};
          G4double C1P[2]={0.003,0.001};
          G4double B0P[2]={2.5,4.5};
          G4double B1P[2]={1.0,4.0};
          
          if(BoundaryP[2]<HadrP) 
          {
            Coeff0=Coeff1=Coeff2=0.0;
            Slope0=Slope1=Slope2=0.0;
          }

          else
          {
            iStep = 0;  while(HadrP>EnP[iStep])  iStep++;

            Coeff0 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  C0P[iStep], C0P[iStep-1],  
				  HadrP);            

            Coeff1 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  C1P[iStep], C1P[iStep-1],  
				  HadrP);

            Slope0 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  B0P[iStep], B0P[iStep-1],  
				  HadrP);

            Slope1 = LineInterpol(EnP[iStep], EnP[iStep-1],
                                  B1P[iStep], B1P[iStep-1],  
				  HadrP);
          }       //  else

            Coeff2=0.0;
            Slope2=5.0;
        } break;  //  case -2212

      case 211 : 
        {
          G4double EnPP[4]={1.0,2.0,3.0,4.0};
          G4double C0PP[4]={0.0,0.0,0.0,0.0};
          G4double C1PP[4]={0.15,0.08,0.02,0.01};
          G4double B0PP[4]={1.5,2.8,3.8,3.8};
          G4double B1PP[4]={0.8,1.6,3.6,4.6};

          if(BoundaryP[3]<HadrP) 
          {
            Coeff0=Coeff1=Coeff2=0.0;
            Slope0=Slope1=Slope2=0.0;
          }
          else
          {
            iStep = 0;  while(HadrP>EnPP[iStep])   iStep++;

            Coeff0 = LineInterpol(EnPP[iStep], EnPP[iStep-1],
                                  C0PP[iStep], C0PP[iStep-1],  
				  HadrP);            

            Coeff1 = LineInterpol(EnPP[iStep], EnPP[iStep-1],
                                  C1PP[iStep], C1PP[iStep-1],  
				  HadrP);

            Slope0 = LineInterpol(EnPP[iStep], EnPP[iStep-1],
                                  B0PP[iStep], B0PP[iStep-1],  
				  HadrP);

            Slope1 = LineInterpol(EnPP[iStep], EnPP[iStep-1],
                                  B1PP[iStep], B1PP[iStep-1],  
				  HadrP);
            Coeff2 = 0.02/HadrP;
            Slope2 = 5;
          }       //  else
        } break;  //  case 211

      case -211 : 
        {
          G4double EnPP[4]={1.0,2.0,3.0,4.0};
          G4double C0PP[4]={0.0,0.0,0.0,0.0};
          G4double C1PP[4]={0.15,0.08,0.02,0.01};
          G4double B0PP[4]={1.5,2.8,3.8,3.8};
          G4double B1PP[4]={0.8,1.6,3.6,4.6};

          if(BoundaryP[4]<HadrP) 
          {
            Coeff0=Coeff1=Coeff2=0.0;
            Slope0=Slope1=Slope2=0.0;
          }
          else
          {
            iStep = 0;  while(HadrP>EnPP[iStep])    iStep++;

            Coeff0 = LineInterpol(EnPP[iStep], EnPP[iStep-1],
                                  C0PP[iStep], C0PP[iStep-1],  
				  HadrP);            

            Coeff1 = LineInterpol(EnPP[iStep], EnPP[iStep-1],
                                  C1PP[iStep], C1PP[iStep-1],  
				  HadrP);

            Slope0 = LineInterpol(EnPP[iStep], EnPP[iStep-1],
                                  B0PP[iStep], B0PP[iStep-1],  
				  HadrP);

            Slope1 = LineInterpol(EnPP[iStep], EnPP[iStep-1],
                                  B1PP[iStep], B1PP[iStep-1],  
				  HadrP);
            Coeff2 = 0.02/HadrP;
            Slope2 = 5;
          }       //  else
        } break;  //  case -211

      case 321 : 
        {
	  G4double EnK[4]={1.4,2.33,3.0,5.0};
          G4double C0K[4]={0.0,0.0,0.0,0.0};
          G4double C1K[4]={0.01,0.007,0.005,0.003};
          G4double B0K[4]={1.5,2.0,3.8,3.8};
          G4double B1K[4]={1.6,1.6,1.6,1.6};
          
          if(BoundaryP[5]<HadrP) 
          {
            Coeff0=Coeff1=Coeff2=0.0;
            Slope0=Slope1=Slope2=0.0;
          }
          else
          {
            iStep = 0;  while(HadrP>EnK[iStep])   iStep++;

            Coeff0 = LineInterpol(EnK[iStep], EnK[iStep-1],
                                  C0K[iStep], C0K[iStep-1],  
                                  HadrP);            

            Coeff1 = LineInterpol(EnK[iStep], EnK[iStep-1],
                                  C1K[iStep], C1K[iStep-1],  
                                  HadrP);

            Slope0 = LineInterpol(EnK[iStep], EnK[iStep-1],
                                  B0K[iStep], B0K[iStep-1],  
                                  HadrP);

            Slope1 = LineInterpol(EnK[iStep], EnK[iStep-1],
                                  B1K[iStep], B1K[iStep-1],  
                                  HadrP);
            Coeff2 = 0.34/HadrP/HadrP/HadrP;
            Slope2 = 5;
          }       //  else
        } break;  //  case 321

      case -321 : 
        {
	  G4double EnKM[2]={1.4,4.0};
          G4double C0KM[2]={0.006,0.002};
          G4double C1KM[2]={0.00,0.00};
          G4double B0KM[2]={2.5,3.5};
          G4double B1KM[2]={1.6,1.6};
          
          if(BoundaryP[6]<HadrP) 
          {
            Coeff0=Coeff1=Coeff2=0.0;
            Slope0=Slope1=Slope2=0.0;
          }
          else
          {
            iStep = 0;  while(HadrP>EnKM[iStep])   iStep++;

            Coeff0 = LineInterpol(EnKM[iStep], EnKM[iStep-1],
                                  C0KM[iStep], C0KM[iStep-1],  
				  HadrP);            

            Coeff1 = LineInterpol(EnKM[iStep], EnKM[iStep-1],
                                  C1KM[iStep], C1KM[iStep-1],  
				  HadrP);

            Slope0 = LineInterpol(EnKM[iStep], EnKM[iStep-1],
                                  B0KM[iStep], B0KM[iStep-1],  
				  HadrP);

            Slope1 = LineInterpol(EnKM[iStep], EnKM[iStep-1],
                                  B1KM[iStep], B1KM[iStep-1],  
				  HadrP);
            Coeff2 = 0.01/HadrP/HadrP/HadrP;
            Slope2 = 5;
          }       //  else
        } break;  //  case -321
    }
  }
//  ++++++++++++++++++++++++++++++++++++++++
  G4double G4DiffElasticHadrNucleus::
               LineInterpol(G4double p1, G4double p2,
                            G4double c1, G4double c2, 
                            G4double p)
  {
    G4double c;

    c = c1+(p-p1)/(p2-p1)*(c2-c1);
    return c;
  }
//  ++++++++++++++++++++++++++++++++++++++++
/*  End of file  */
