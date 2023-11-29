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
// ABLAXX statistical de-excitation model
// Jose Luis Rodriguez, GSI (translation from ABLA07 and contact person)
// Pekka Kaitaniemi, HIP (initial translation of ablav3p)
// Aleksandra Kelic, GSI (ABLA07 code)
// Davide Mancusi, CEA (contact person INCL)
// Aatos Heikkinen, HIP (project coordination)
//

#define ABLAXX_IN_GEANT4_MODE 1

#include "globals.hh"
#include <time.h>
#include <cmath>

#include "G4Abla.hh"
#include "G4AblaDataFile.hh"
#include "G4AblaRandom.hh"
#ifdef ABLAXX_IN_GEANT4_MODE
G4Abla::G4Abla(G4Volant *aVolant, G4VarNtp *aVarntp)
#else
G4Abla::G4Abla(G4INCL::Config *config, G4Volant *aVolant, G4VarNtp *aVarntp)
#endif
{
#ifndef ABLAXX_IN_GEANT4_MODE
  theConfig = config;
#endif
  verboseLevel = 0;
  ilast = 0;
  volant = aVolant; // ABLA internal particle data
  volant->iv = 0;
  varntp = aVarntp; // Output data structure
  varntp->ntrack = 0;
 
  verboseLevel = 0;
  gammaemission= 0;// 0 presaddle, 1 postsaddle
  T_freeze_out = 0.;
  Ainit=0;
  Zinit=0;
  Sinit=0;

  pace = new G4Pace();
  ald = new G4Ald();
  eenuc = new G4Eenuc();
  ec2sub = new G4Ec2sub();
  ecld = new G4Ecld();
  masses = new G4Mexp();
  fb = new G4Fb();
  fiss = new G4Fiss();
  opt = new G4Opt();
}

void G4Abla::setVerboseLevel(G4int level)
{
  verboseLevel = level;
}

G4Abla::~G4Abla()
{
  delete pace;
  delete ald;
  delete eenuc;
  delete ec2sub;
  delete ecld;
  delete masses;
  delete fb;
  delete fiss;
  delete opt;
}

// Main interface to the evaporation without lambda evaporation
void G4Abla::DeexcitationAblaxx(G4int nucleusA, G4int nucleusZ, G4double excitationEnergy, G4double angularMomentum, G4double momX, G4double momY, G4double momZ, G4int eventnumber)
{
 DeexcitationAblaxx(nucleusA,nucleusZ,excitationEnergy,angularMomentum,momX,momY,momZ,eventnumber,0);
}

// Main interface to the evaporation with lambda emission
void G4Abla::DeexcitationAblaxx(G4int nucleusA, G4int nucleusZ, G4double excitationEnergy, G4double angularMomentum, G4double momX, G4double momY, G4double momZ, G4int eventnumber, G4int nucleusS)
{

  const G4double amu = 931.4940; //  MeV/C^2
  const G4double C = 29.9792458; // cm/ns

  SetParametersG4(nucleusZ, nucleusA);

  mult10:
  G4int IS = 0;

  if(nucleusS>0)nucleusS=0;// S=1 from INCL ????

  G4int NbLam0 = std::abs(nucleusS);

  Ainit=-1*nucleusA;
  Zinit=-1*nucleusZ;
  Sinit=-1*nucleusS;

  G4double aff = 0.0;
  G4double zff = 0.0;
  G4int ZFP1 = 0, AFP1 = 0, AFPIMF = 0, ZFPIMF = 0, ZFP2 = 0, AFP2 = 0, SFP1 = 0, SFP2 = 0, SFPIMF = 0;
  G4double vx_eva = 0.0, vy_eva = 0.0, vz_eva = 0.0;  
  G4double VX_PREF=0.,VY_PREF=0.,VZ_PREF=00,VP1X,VP1Y,VP1Z,VXOUT,VYOUT,VZOUT,V_CM[3],VFP1_CM[3],VFP2_CM[3],VIMF_CM[3],VX2OUT,VY2OUT,VZ2OUT;
  G4double zf = 0.0, af = 0.0, mtota = 0.0, tkeimf = 0.0, jprf0=0.;
  G4int ff = 0,afpnew=0,zfpnew=0,aprfp=0,zprfp=0,IOUNSTABLE=0,ILOOP=0,IEV_TAB=0,IEV_TAB_TEMP=0;
  G4int fimf = 0,INMIN=0,INMAX=0;
  G4int ftype=0;//,ftype1=0;
  G4int inum = eventnumber;
  G4int inttype = 0;
  opt->optimfallowed=1;

  if(fiss->zt>56){
  fiss->ifis = 1;
  }else {
  fiss->ifis = 0;
  }

  if(NbLam0>0){
   opt->nblan0 = NbLam0;
  }
  
  G4double aprf = (G4double) nucleusA;
  G4double zprf = (G4double) nucleusZ;
  G4double ee = excitationEnergy;
  G4double jprf = angularMomentum; // actually root-mean-squared

  G4double pxrem = momX;
  G4double pyrem = momY;
  G4double pzrem = momZ;
  G4double zimf,aimf;

  volant->clear(); // Clean up an initialize ABLA output.
  varntp->clear(); // Clean up an initialize ABLA output.
  varntp->ntrack = 0;
  varntp->kfis = 0;
  volant->iv = 0;
  gammaemission=0;
  G4double T_init=0.,T_diff=0.,a_tilda=0.,a_tilda_BU=0., EE_diff=0., EINCL=0., A_FINAL=0., Z_FINAL=0., E_FINAL=0.;

  G4double A_diff=0.,ASLOPE1,ASLOPE2,A_ACC,ABU_SLOPE, ABU_SUM=0., AMEM=0., ZMEM=0., EMEM=0., JMEM=0., PX_BU_SUM = 0.0, PY_BU_SUM = 0.0, PZ_BU_SUM = 0.0, ETOT_SUM=0., P_BU_SUM=0., ZBU_SUM=0.,Z_Breakup_sum=0.,A_Breakup,Z_Breakup,N_Breakup,G_SYMM,CZ,Sigma_Z,Z_Breakup_Mean,ZTEMP=0.,ATEMP=0.;

  G4double ETOT_PRF=0.0,PXPRFP=0.,PYPRFP=0.,PZPRFP=0.,PPRFP=0., VX1_BU=0., VY1_BU=0., VZ1_BU=0., VBU2=0., GAMMA_REL=1.0, Eexc_BU_SUM=0., VX_BU_SUM = 0., VY_BU_SUM =0.,VZ_BU_SUM =0., E_tot_BU=0.,EKIN_BU=0.,ZIMFBU=0., AIMFBU=0., ZFFBU=0., AFFBU=0., AFBU=0., ZFBU=0., EEBU=0.,TKEIMFBU=0.,vx_evabu=0.,vy_evabu=0.,vz_evabu=0., Bvalue_BU=0.,P_BU=0.,ETOT_BU=1.,PX_BU=0.,PY_BU=0.,PZ_BU=0.,VX2_BU=0.,VY2_BU=0.,VZ2_BU=0.;

  G4int ABU_DIFF,ZBU_DIFF,NBU_DIFF;
  G4int INEWLOOP = 0, ILOOPBU=0;

  G4double BU_TAB_TEMP[200][6], BU_TAB_TEMP1[200][6]; 
  G4double EV_TAB_TEMP[200][6],EV_TEMP[200][6];
  G4int IMEM_BU[200], IMEM=0;

 if(nucleusA<1){
  std::cout << "Error - Remnant with a mass number A below 1." << std::endl;
 //INCL_ERROR("Remnant with a mass number A below 1.");
 return;
 }

  for(G4int j=0;j<3;j++){
   V_CM[j]=0.;
   VFP1_CM[j]=0.;
   VFP2_CM[j]=0.;
   VIMF_CM[j]=0.;
  }

   for(G4int I1=0;I1<200;I1++){
       for(G4int I2 = 0;I2<12;I2++)
        BU_TAB[I1][I2] = 0.0;
       for(G4int I2 = 0;I2<6;I2++){
        BU_TAB_TEMP[I1][I2] = 0.0;
        BU_TAB_TEMP1[I1][I2] = 0.0;
        EV_TAB_TEMP[I1][I2] = 0.0;
        EV_TAB[I1][I2] = 0.0;
        EV_TAB_SSC[I1][I2] = 0.0;
        EV_TEMP[I1][I2] = 0.0;
       }
   }

  G4int idebug = 0;
  if(idebug == 1) {
    zprf =   81.;
    aprf =   201.;
//    ee =   86.5877686;
    ee = 100.0;
    jprf =   10.;
    zf =   0.;
    af =   0.;
    mtota =   0.;
    ff =  1;
    inttype =  0;
    //inum =  2;
  }
//
      G4double AAINCL = aprf;
      G4double ZAINCL = zprf;
      EINCL = ee;
//
// Velocity after the first stage of reaction (INCL)
// For coupling with INCL, comment the lines below, and use output
// of INCL as pxincl, pyincl,pzincl
//
      G4double pincl = std::sqrt(pxrem*pxrem + pyrem*pyrem + pzrem*pzrem);
// PPRFP is in MeV/c
      G4double ETOT_incl = std::sqrt(pincl*pincl + (AAINCL * amu)*(AAINCL * amu));
      G4double VX_incl = C * pxrem / ETOT_incl;
      G4double VY_incl = C * pyrem / ETOT_incl;
      G4double VZ_incl = C * pzrem / ETOT_incl;
//
// Multiplicity in the break-up event
      G4int  IMULTBU = 0;
      G4int  IMULTIFR = 0;
      G4int  I_Breakup=0;
      G4int  NbLamprf= 0;
           IEV_TAB = 0;
/*
C     Set maximum temperature for sequential decay (evaporation)
C     Remove additional energy by simultaneous break up
C                          (vaporisation or multi-fragmentation)

C     Idea: If the temperature of the projectile spectator exceeds
c           the limiting temperature T_freeze_out, the additional
C           energy which is present in the spectator is used for
C           a stage of simultaneous break up. It is either the
C           simultaneous emission of a gaseous phase or the simultaneous
C           emission of several intermediate-mass fragments. Only one
C           piece of the projectile spectator (assumed to be the largest
C           one) is kept track.

C        MVR, KHS, October 2001
C        KHS, AK 2007 - Masses from the power low; slope parameter dependent on
C                      energy  per nucleon; symmtery-energy coeff. dependent on
C                      energy per nucleon.

c       Clear BU_TAB (array of multifragmentation products)
*/
        if(T_freeze_out_in >= 0.0){
          T_freeze_out = T_freeze_out_in;
        }else{
         T_freeze_out = max(9.33*std::exp(-0.00282*AAINCL),5.5);
//         ! See: J. Natowitz et al, PRC65 (2002) 034618
//        T_freeze_out=DMAX1(9.0D0*DEXP(-0.001D0*AAABRA),
//     &                     5.5D0)
        }
//
        a_tilda = ald->av*aprf + ald->as*std::pow(aprf,2.0/3.0) + ald->ak*std::pow(aprf,1.0/3.0);

        T_init = std::sqrt(EINCL/a_tilda);

        T_diff = T_init - T_freeze_out;

        if(T_diff>0.1 && zprf>2. && (aprf-zprf)>0.){
        // T_Diff is set to be larger than 0.1 MeV in order to avoid strange cases for which
        // T_Diff is of the order of 1.e-3 and less.
        varntp->kfis = 10;

        for(G4int i=0;i<5;i++){
            EE_diff = EINCL - a_tilda * T_freeze_out*T_freeze_out;
//            Energy removed 10*5/T_init per nucleon removed in simultaneous breakup
//            adjusted to frag. xsections 238U (1AGeV) + Pb data, KHS Dec. 2005
// This should maybe be re-checked, in a meanwhile several things in break-up description
// have changed (AK).

            A_diff = dint(EE_diff / (8.0 * 5.0 / T_freeze_out));

            if(A_diff>AAINCL) A_diff = AAINCL;

            A_FINAL = AAINCL - A_diff;

            a_tilda = ald->av*A_FINAL + ald->as*std::pow(A_FINAL,2.0/3.0) + ald->ak*std::pow(A_FINAL,1.0/3.0);
            E_FINAL = a_tilda * T_freeze_out*T_freeze_out;

            if(A_FINAL<4.0){  // To avoid numerical problems
              EE_diff = EINCL - E_FINAL;
              A_FINAL = 1.0;
              Z_FINAL = 1.0;
              E_FINAL = 0.0;
              goto mul4325;
            }
        }
        mul4325:
// The idea is similar to Z determination of multifragment - Z of "heavy" partner is not
// fixed by the A/Z of the prefragment, but randomly picked from Gaussian
         // Z_FINAL_MEAN = dint(zprf * A_FINAL / (aprf));

          Z_FINAL = dint(zprf * A_FINAL / (aprf));

          if(E_FINAL<0.0) E_FINAL = 0.0;

          aprf = A_FINAL;
          zprf = Z_FINAL;
          ee = E_FINAL;

          A_diff = AAINCL - aprf;

// Creation of multifragmentation products by breakup
          if(A_diff<=1.0){
           aprf = AAINCL;
           zprf = ZAINCL;
           ee = EINCL;
           IMULTIFR = 0;
           goto mult7777;
          }else if(A_diff>1.0){

          A_ACC = 0.0;
// Energy-dependence of the slope parameter, acc. to A. Botvina, fits also to exp. data (see
// e.g. Sfienti et al, NPA 2007)
          ASLOPE1 = -2.400;  // e*/a=7   -2.4
          ASLOPE2 = -1.200;  // e*/a=3   -1.2

          a_tilda = ald->av*AAINCL + ald->as*std::pow(AAINCL,2.0/3.0) + ald->ak*std::pow(AAINCL,1.0/3.0);

          E_FINAL = a_tilda * T_freeze_out*T_freeze_out;

          ABU_SLOPE = (ASLOPE1-ASLOPE2)/4.0*(E_FINAL/AAINCL)+
                     ASLOPE1-(ASLOPE1-ASLOPE2)*7.0/4.0;

// Botvina et al, PRC 74 (2006) 044609, fig. 5 for B0=18 MeV
//          ABU_SLOPE = 5.57489D0-2.08149D0*(E_FINAL/AAABRA)+
//     &    0.3552D0*(E_FINAL/AAABRA)**2-0.024927D0*(E_FINAL/AAABRA)**3+
//     &    7.268D-4*(E_FINAL/AAABRA)**4
// They fit with A**(-tau) and here is done A**(tau)
//          ABU_SLOPE = ABU_SLOPE*(-1.D0)

//           ABU_SLOPE = -2.60D0
//          print*,ABU_SLOPE,(E_FINAL/AAABRA)

          if(ABU_SLOPE > -1.01) ABU_SLOPE = -1.01;

          I_Breakup = 0;
          Z_Breakup_sum = Z_FINAL;
          ABU_SUM = 0.0;
          ZBU_SUM = 0.0;

          for(G4int i=0;i<100;i++){
             IS = 0;
             mult4326:
             A_Breakup = dint(G4double(IPOWERLIMHAZ(ABU_SLOPE,1,idnint(A_diff))));
                  // Power law with exponent ABU_SLOPE
             IS = IS +1;
             if(IS>100){
             std::cout << "WARNING: IPOWERLIMHAZ CALLED MORE THAN 100 TIMES WHEN CALCULATING A_BREAKUP IN Rn07.FOR. NEW EVENT WILL BE DICED: " << A_Breakup << std::endl;
             goto mult10;
             }

             if(A_Breakup>AAINCL) goto mult4326;

             if(A_Breakup<=0.0){
              std::cout << "A_BREAKUP <= 0 " << std::endl;
              goto mult10;
             }

            A_ACC = A_ACC + A_Breakup;

            if(A_ACC<=A_diff){

              Z_Breakup_Mean = dint(A_Breakup * ZAINCL / AAINCL);

              Z_Breakup_sum = Z_Breakup_sum + Z_Breakup_Mean;
//
// See G.A. Souliotis et al, PRC 75 (2007) 011601R (Fig. 2)
              G_SYMM = 34.2281 - 5.14037 * E_FINAL/AAINCL;
              if(E_FINAL/AAINCL < 2.0) G_SYMM = 25.0;
              if(E_FINAL/AAINCL > 4.0) G_SYMM = 15.0;

//             G_SYMM = 23.6;

              G_SYMM = 25.0;      //25
              CZ = 2.0 * G_SYMM * 4.0 / A_Breakup;
              // 2*CZ=d^2(Esym)/dZ^2, Esym=Gamma*(A-2Z)**2/A
               // gamma = 23.6D0 is the symmetry-energy coefficient
              G4int IIS = 0;
              Sigma_Z = std::sqrt(T_freeze_out/CZ);

              IS = 0;
              mult4333:
              Z_Breakup =  dint( G4double(gausshaz(1,Z_Breakup_Mean,Sigma_Z)));
              IS = IS +1;
//
              if(IS>100){
               std::cout << "WARNING: GAUSSHAZ CALLED MORE THAN 100 TIMES WHEN CALCULATING Z_BREAKUP IN Rn07.FOR. NEW EVENT WILL BE DICED: " << A_Breakup << " " << Z_Breakup << std::endl;
               goto mult10;
              }

             if(Z_Breakup<0.0 ) goto mult4333;
             if((A_Breakup-Z_Breakup)<0.0) goto mult4333;
             if((A_Breakup-Z_Breakup)==0.0 && Z_Breakup!=1.0) goto mult4333;

             if(Z_Breakup>=ZAINCL){
               IIS = IIS + 1;
                 if(IIS > 10){
                   std::cout << "Z_BREAKUP RESAMPLED MORE THAN 10 TIMES; EVENT WILL BE RESAMPLED AGAIN " << std::endl;
                   goto mult10;
                 }
               goto mult4333;
             }

//     *** Find the limits that fragment is bound :
        isostab_lim(idnint(Z_Breakup),&INMIN,&INMAX);
//        INMIN = MAX(1,INMIN-2)
        if(Z_Breakup > 2.0){
          if(idnint(A_Breakup-Z_Breakup)<INMIN || idnint(A_Breakup-Z_Breakup)>(INMAX+5)){
//             PRINT*,'N_Breakup >< NMAX',
//     &      IDNINT(Z_Breakup),IDNINT(A_Breakup-Z_Breakup),INMIN,INMAX
            goto mult4343;
          }
        }

  mult4343:

// We consider all products, also nucleons created in the break-up
//               I_Breakup = I_Breakup + 1;// moved below

               N_Breakup = A_Breakup - Z_Breakup;
               BU_TAB[I_Breakup][0] = dint(Z_Breakup);    // Mass of break-up product
               BU_TAB[I_Breakup][1] = dint(A_Breakup);    // Z of break-up product
               ABU_SUM = ABU_SUM + BU_TAB[i][1];
               ZBU_SUM = ZBU_SUM + BU_TAB[i][0];
//
// Break-up products are given zero angular momentum (simplification)
               BU_TAB[I_Breakup][3] = 0.0;
               I_Breakup = I_Breakup + 1;
               IMULTBU = IMULTBU + 1;
           }else{
//     There are A_DIFF - A_ACC nucleons lost by breakup, but they do not end up in multifragmentation products.
//     This is a deficiency of the Monte-Carlo method applied above to determine the sizes of the fragments
//     according to the power law.
//            print*,'Deficiency',IDNINT(A_DIFF-A_ACC)

             goto mult4327;
           }// if(A_ACC<=A_diff)
          }//for
          //mult4327:
          //IMULTIFR = 1;
          } //  if(A_diff>1.0)
          mult4327:
          IMULTIFR = 1;

// "Missing" A and Z picked from the power law:
        ABU_DIFF = idnint(ABU_SUM+aprf-AAINCL);
        ZBU_DIFF = idnint(ZBU_SUM+zprf-ZAINCL);
        NBU_DIFF = idnint((ABU_SUM-ZBU_SUM)+(aprf-zprf)-(AAINCL-ZAINCL));
//
        if(IMULTBU > 200)
        std::cout << "WARNING - MORE THAN 200 BU " << IMULTBU  << std::endl;

        if(IMULTBU < 1)
        std::cout << "WARNING - LESS THAN 1 BU " << IMULTBU << std::endl; 
         //,AABRA,ZABRA,IDNINT(APRF),IDNINT(ZPRF),ABU_DIFF,ZBU_DIFF

        G4int IPROBA = 0;
        for(G4int i=0;i<IMULTBU;i++)
        IMEM_BU[i] = 0;

        while(NBU_DIFF!=0 && ZBU_DIFF!=0){
// (APRF,ZPRF) is also inlcuded in this game, as from time to time the program
// is entering into endless loop, as it can not find proper nucleus for adapting A and Z.
         IS = 0;
         mult5555:    
         G4double RHAZ = G4AblaRandom::flat()*G4double(IMULTBU);
         IPROBA = IPROBA + 1;
         IS = IS + 1;
         if(IS>100){
          std::cout << "WARNING: HAZ CALLED MORE THAN 100 TIMES WHEN CALCULATING N_BREAKUP IN Rn07.FOR. NEW EVENT WILL BE DICED." << std::endl; 
          goto mult10;
         }
          G4int IEL = G4int(RHAZ);
          if(IMEM_BU[IEL]==1) goto mult5555;
	  if(!(IEL<200))std::cout << "5555:" << IEL << RHAZ << IMULTBU << std::endl; 
           if(IEL<0)std::cout << "5555:"<< IEL << RHAZ << IMULTBU << std::endl; 
           if(IEL<=IMULTBU){
            N_Breakup = dint(BU_TAB[IEL][1]-BU_TAB[IEL][0] - DSIGN(1.0,G4double(NBU_DIFF)));
            }else if(IEL>IMULTBU){
            N_Breakup = dint(aprf - zprf - DSIGN(1.0,G4double(NBU_DIFF)));
            }
            if(N_Breakup<0.0){
             IMEM_BU[IEL] = 1;
             goto mult5555;
            }
             if(IEL<=IMULTBU){
             ZTEMP = dint(BU_TAB[IEL][0] - DSIGN(1.0,G4double(ZBU_DIFF)));
             }else if(IEL>IMULTBU){
             ZTEMP = dint(zprf - DSIGN(1.0,G4double(ZBU_DIFF)));
             }
              if(ZTEMP<0.0){
               IMEM_BU[IEL] = 1;
               goto mult5555;
              }
              if(ZTEMP<1.0 && N_Breakup<1.0){
               IMEM_BU[IEL] = 1;
               goto mult5555;
              }
// Nuclei with A=Z and Z>1 are allowed in this stage, as otherwise,
// for more central collisions there is not enough mass which can be
// shufeled in order to conserve A and Z. These are mostly nuclei with
// Z=2 and in less extent 3, 4 or 5.
//             IF(ZTEMP.GT.1.D0 .AND. N_Breakup.EQ.0.D0) THEN
//              GOTO 5555
//             ENDIF
            if(IEL<=IMULTBU){
            BU_TAB[IEL][0] = dint(ZTEMP);
            BU_TAB[IEL][1] = dint(ZTEMP + N_Breakup);
            }else if(IEL>IMULTBU){
            zprf = dint(ZTEMP);
            aprf = dint(ZTEMP + N_Breakup);
            }
          NBU_DIFF = NBU_DIFF - ISIGN(1,NBU_DIFF);
          ZBU_DIFF = ZBU_DIFF - ISIGN(1,ZBU_DIFF);
        }// while

        IPROBA = 0;
        for(G4int i=0;i<IMULTBU;i++)
        IMEM_BU[i] = 0;

        if(NBU_DIFF != 0 && ZBU_DIFF == 0){
         while(NBU_DIFF > 0 || NBU_DIFF < 0){
         IS = 0;
         mult5556:    
         G4double RHAZ = G4AblaRandom::flat()*G4double(IMULTBU);
         IS = IS + 1;
         if(IS>100){
          std::cout << "WARNING: HAZ CALLED MORE THAN 100 TIMES WHEN CALCULATING N_BREAKUP IN Rn07.FOR. NEW EVENT WILL BE DICED." << std::endl; 
          goto mult10;
         }
         G4int IEL = G4int(RHAZ);
         if(IMEM_BU[IEL]==1) goto mult5556;
//         IPROBA = IPROBA + 1;
         if(IPROBA>IMULTBU+1 && NBU_DIFF>0){
         std::cout << "###',IPROBA,IMULTBU,NBU_DIFF,ZBU_DIFF,T_freeze_out" << std::endl; 
         IPROBA = IPROBA + 1;
           if(IEL<=IMULTBU){
            BU_TAB[IEL][1] = dint(BU_TAB[IEL][1]-G4double(NBU_DIFF));
           }else{ if(IEL>IMULTBU)
            aprf = dint(aprf - G4double(NBU_DIFF));
           }
         goto mult5432;
         }
	 if(!(IEL<200))std::cout << "5556:" << IEL << RHAZ << IMULTBU << std::endl; 
           if(IEL<0)std::cout << "5556:"<< IEL << RHAZ << IMULTBU << std::endl; 
           if(IEL<=IMULTBU){
            N_Breakup = dint(BU_TAB[IEL][1]-BU_TAB[IEL][0] - DSIGN(1.0,G4double(NBU_DIFF)));
            }else if(IEL>IMULTBU){
            N_Breakup = dint(aprf - zprf - DSIGN(1.0,G4double(NBU_DIFF)));
            }
            if(N_Breakup<0.0){
             IMEM_BU[IEL] = 1;
             goto mult5556;
            }
            if(IEL<=IMULTBU){
             ATEMP = dint(BU_TAB[IEL][0] + N_Breakup);
             }else if(IEL>IMULTBU){
             ATEMP = dint(zprf + N_Breakup);
             }
             if((ATEMP - N_Breakup)<1.0 && N_Breakup<1.0){
              IMEM_BU[IEL] = 1;
              goto mult5556;
             }
//             IF((ATEMP - N_Breakup).GT.1.D0 .AND.
//     &        N_Breakup.EQ.0.D0) THEN
//              IMEM_BU(IEL) = 1
//              GOTO 5556
//             ENDIF
             if(IEL<=IMULTBU)
             BU_TAB[IEL][1] = dint(BU_TAB[IEL][0] + N_Breakup);
             else if(IEL>IMULTBU)
             aprf = dint(zprf + N_Breakup);
//
           NBU_DIFF = NBU_DIFF - ISIGN(1,NBU_DIFF);
         }//while(NBU_DIFF > 0 || NBU_DIFF < 0)

        IPROBA = 0;
        for(G4int i=0;i<IMULTBU;i++)
        IMEM_BU[i] = 0;

        }else{// if(NBU_DIFF != 0 && ZBU_DIFF == 0)
          if(ZBU_DIFF != 0 && NBU_DIFF == 0){
            while(ZBU_DIFF > 0 || ZBU_DIFF < 0){
             IS = 0;
             mult5557:    
             G4double RHAZ = G4AblaRandom::flat()*G4double(IMULTBU);
             IS = IS + 1;
             if(IS>100){
              std::cout << "WARNING: HAZ CALLED MORE THAN 100 TIMES WHEN CALCULATING N_BREAKUP IN Rn07.FOR. NEW EVENT WILL BE DICED." << std::endl; 
              goto mult10;
             }
             G4int IEL = G4int(RHAZ);
             if(IMEM_BU[IEL]==1) goto mult5557;
             //IPROBA = IPROBA + 1;
             if(IPROBA>IMULTBU+1 && ZBU_DIFF>0){
             std::cout << "###',IPROBA,IMULTBU,NBU_DIFF,ZBU_DIFF,T_freeze_out" << std::endl; 
             IPROBA = IPROBA + 1;
             if(IEL<=IMULTBU){
               N_Breakup = dint(BU_TAB[IEL][1]-BU_TAB[IEL][0]);
               BU_TAB[IEL][0] = dint(BU_TAB[IEL][0] - G4double(ZBU_DIFF));
               BU_TAB[IEL][1] = dint(BU_TAB[IEL][0] + N_Breakup);
             }else{ 
                 if(IEL>IMULTBU){
                  N_Breakup = aprf - zprf;
                  zprf = dint(zprf - G4double(ZBU_DIFF));
                  aprf = dint(zprf + N_Breakup);
                 }
             }
          goto mult5432;
          }
	   if(!(IEL<200))std::cout << "5557:" << IEL << RHAZ << IMULTBU << std::endl; 
           if(IEL<0)std::cout << "5557:"<< IEL << RHAZ << IMULTBU << std::endl; 
           if(IEL<=IMULTBU){
            N_Breakup = dint(BU_TAB[IEL][1]-BU_TAB[IEL][0]);
            ZTEMP = dint(BU_TAB[IEL][0] - DSIGN(1.0,G4double(ZBU_DIFF)));
           }else if(IEL>IMULTBU){
            N_Breakup = dint(aprf - zprf);
            ZTEMP = dint(zprf - DSIGN(1.0,G4double(ZBU_DIFF)));
           }
            ATEMP = dint(ZTEMP + N_Breakup);
            if(ZTEMP<0.0){
             IMEM_BU[IEL] = 1;
             goto mult5557;
            }
            if((ATEMP-ZTEMP)<0.0){
             IMEM_BU[IEL] = 1;
             goto mult5557;
            }
            if((ATEMP-ZTEMP)<1.0 && ZTEMP<1.0){
             IMEM_BU[IEL] = 1;
             goto mult5557;
            }
             if(IEL<=IMULTBU){
               BU_TAB[IEL][0] = dint(ZTEMP);
               BU_TAB[IEL][1] = dint(ZTEMP + N_Breakup);
             }else{ 
                 if(IEL>IMULTBU){
                  zprf = dint(ZTEMP);
                  aprf = dint(ZTEMP + N_Breakup);
                 }
             } 
            ZBU_DIFF = ZBU_DIFF - ISIGN(1,ZBU_DIFF);
            }//while
          }//if(ZBU_DIFF != 0 && NBU_DIFF == 0)
        }// if(NBU_DIFF != 0 && ZBU_DIFF == 0)

        mult5432:
// Looking for the heaviest fragment among all multifragmentation events, and
// "giving" excitation energy to fragments
         ZMEM = 0.0;

         for(G4int i =0;i<IMULTBU;i++){
//For particles with Z>2 we calculate excitation energy from freeze-out temperature.
// For particels with Z<3 we assume that they form a gas, and that temperature results
// in kinetic energy (which is sampled from Maxwell distribution with T=Tfreeze-out)
// and not excitation energy.
           if(BU_TAB[i][0]>2.0){
            a_tilda_BU = ald->av*BU_TAB[i][1] + ald->as*std::pow(BU_TAB[i][1],2.0/3.0) + ald->ak*std::pow(BU_TAB[i][1],1.0/3.0);
            BU_TAB[i][2] = a_tilda_BU * T_freeze_out*T_freeze_out; // E* of break-up product
           }else{
            BU_TAB[i][2] = 0.0;
           }
//
           if(BU_TAB[i][0] > ZMEM){
            IMEM = i;
            ZMEM = BU_TAB[i][0];
            AMEM = BU_TAB[i][1];
            EMEM = BU_TAB[i][2];
            JMEM = BU_TAB[i][3];
           }
         }//for IMULTBU

         if(zprf < ZMEM){
          BU_TAB[IMEM][0] = zprf;
          BU_TAB[IMEM][1] = aprf;
          BU_TAB[IMEM][2] = ee;
          BU_TAB[IMEM][3] = jprf;
          zprf =  ZMEM;
          aprf =  AMEM;
          aprfp = idnint(aprf);
          zprfp = idnint(zprf);
          ee   =  EMEM;
          jprf =  JMEM;
         }

//     Just for checking:
        ABU_SUM = aprf;
        ZBU_SUM = zprf;
        for(G4int i = 0;i<IMULTBU;i++){
         ABU_SUM = ABU_SUM + BU_TAB[i][1];
         ZBU_SUM = ZBU_SUM + BU_TAB[i][0];
        }
        ABU_DIFF = idnint(ABU_SUM-AAINCL);
        ZBU_DIFF = idnint(ZBU_SUM-ZAINCL);
//
        if(ABU_DIFF!=0 || ZBU_DIFF!=0)
         std::cout << "Problem of mass in BU " << ABU_DIFF << " " << ZBU_DIFF << std::endl;
        PX_BU_SUM = 0.0;
        PY_BU_SUM = 0.0;
        PZ_BU_SUM = 0.0;
// Momenta of break-up products are calculated. They are all given in the rest frame
// of the primary prefragment (i.e. after incl):
// Goldhaber model ****************************************
// "Heavy" residue
        AMOMENT(AAINCL,aprf,1,&PXPRFP,&PYPRFP,&PZPRFP);
        PPRFP = std::sqrt(PXPRFP*PXPRFP + PYPRFP*PYPRFP + PZPRFP*PZPRFP);
// ********************************************************
// PPRFP is in MeV/c
        ETOT_PRF = std::sqrt(PPRFP*PPRFP + (aprf * amu)*(aprf * amu));
        VX_PREF = C * PXPRFP / ETOT_PRF;
        VY_PREF = C * PYPRFP / ETOT_PRF;
        VZ_PREF = C * PZPRFP / ETOT_PRF;

// Contribution from Coulomb repulsion ********************
        tke_bu(zprf,aprf,ZAINCL,AAINCL,&VX1_BU,&VY1_BU,&VZ1_BU);

// Lorentz kinematics
//        VX_PREF = VX_PREF + VX1_BU
//        VY_PREF = VY_PREF + VY1_BU
//        VZ_PREF = VZ_PREF + VZ1_BU
// Lorentz transformation
        lorentz_boost(VX1_BU,VY1_BU,VZ1_BU,
                VX_PREF,VY_PREF,VZ_PREF,
                &VXOUT,&VYOUT,&VZOUT);

        VX_PREF = VXOUT;
        VY_PREF = VYOUT;
        VZ_PREF = VZOUT;

// Total momentum: Goldhaber + Coulomb
        VBU2 = VX_PREF*VX_PREF + VY_PREF*VY_PREF + VZ_PREF*VZ_PREF;
        GAMMA_REL = std::sqrt(1.0 - VBU2 / (C*C));
        ETOT_PRF = aprf * amu / GAMMA_REL;
        PXPRFP = ETOT_PRF * VX_PREF / C;
        PYPRFP = ETOT_PRF * VY_PREF / C;
        PZPRFP = ETOT_PRF * VZ_PREF / C;

// ********************************************************
//  Momentum: Total width of abrasion and breakup assumed to be given 
//  by Fermi momenta of nucleons
// *****************************************

        PX_BU_SUM = PXPRFP;
        PY_BU_SUM = PYPRFP;
        PZ_BU_SUM = PZPRFP;

        Eexc_BU_SUM = ee;
        Bvalue_BU = eflmac(idnint(aprf),idnint(zprf),1,0);

        for(I_Breakup=0;I_Breakup<IMULTBU;I_Breakup++){
//       For bu products:
          Bvalue_BU = Bvalue_BU + eflmac(idnint(BU_TAB[I_Breakup][1]), idnint(BU_TAB[I_Breakup][0]),1,0);
          Eexc_BU_SUM = Eexc_BU_SUM + BU_TAB[I_Breakup][2];

          AMOMENT(AAINCL,BU_TAB[I_Breakup][1],1,&PX_BU,&PY_BU,&PZ_BU);
          P_BU = std::sqrt(PX_BU*PX_BU + PY_BU*PY_BU + PZ_BU*PZ_BU);
// *******************************************************
//        PPRFP is in MeV/c
          ETOT_BU = std::sqrt(P_BU*P_BU + (BU_TAB[I_Breakup][1]*amu)*(BU_TAB[I_Breakup][1]*amu));
          BU_TAB[I_Breakup][4] = C * PX_BU / ETOT_BU;    // Velocity in x
          BU_TAB[I_Breakup][5] = C * PY_BU / ETOT_BU;    // Velocity in y
          BU_TAB[I_Breakup][6] = C * PZ_BU / ETOT_BU;    // Velocity in z
//        Contribution from Coulomb repulsion:
          tke_bu(BU_TAB[I_Breakup][0],BU_TAB[I_Breakup][1],ZAINCL,AAINCL,&VX2_BU,&VY2_BU,&VZ2_BU);
// Lorentz kinematics
//          BU_TAB(I_Breakup,5) = BU_TAB(I_Breakup,5) + VX2_BU ! velocity change by Coulomb repulsion
//          BU_TAB(I_Breakup,6) = BU_TAB(I_Breakup,6) + VY2_BU
//          BU_TAB(I_Breakup,7) = BU_TAB(I_Breakup,7) + VZ2_BU
// Lorentz transformation
          lorentz_boost(VX2_BU,VY2_BU,VZ2_BU,
                BU_TAB[I_Breakup][4],BU_TAB[I_Breakup][5],BU_TAB[I_Breakup][6],
                &VXOUT,&VYOUT,&VZOUT);

          BU_TAB[I_Breakup][4] = VXOUT;
          BU_TAB[I_Breakup][5] = VYOUT;
          BU_TAB[I_Breakup][6] = VZOUT;

// Total momentum: Goldhaber + Coulomb
          VBU2 = BU_TAB[I_Breakup][4]*BU_TAB[I_Breakup][4] +
                 BU_TAB[I_Breakup][5]*BU_TAB[I_Breakup][5] +
                 BU_TAB[I_Breakup][6]*BU_TAB[I_Breakup][6];
          GAMMA_REL = std::sqrt(1.0 - VBU2 / (C*C));
          ETOT_BU = BU_TAB[I_Breakup][1]*amu/GAMMA_REL;
          PX_BU = ETOT_BU * BU_TAB[I_Breakup][4] / C;
          PY_BU = ETOT_BU * BU_TAB[I_Breakup][5] / C;
          PZ_BU = ETOT_BU * BU_TAB[I_Breakup][6] / C;

          PX_BU_SUM = PX_BU_SUM + PX_BU;
          PY_BU_SUM = PY_BU_SUM + PY_BU;
          PZ_BU_SUM = PZ_BU_SUM + PZ_BU;

         }//for I_Breakup

//   In the frame of source (i.e. prefragment after abrasion or INCL)
        P_BU_SUM = std::sqrt(PX_BU_SUM*PX_BU_SUM + PY_BU_SUM*PY_BU_SUM +
                   PZ_BU_SUM*PZ_BU_SUM);
// ********************************************************
// PPRFP is in MeV/c
        ETOT_SUM = std::sqrt(P_BU_SUM*P_BU_SUM +
                   (AAINCL * amu)*(AAINCL * amu));

        VX_BU_SUM = C * PX_BU_SUM / ETOT_SUM;
        VY_BU_SUM = C * PY_BU_SUM / ETOT_SUM;
        VZ_BU_SUM = C * PZ_BU_SUM / ETOT_SUM;

// Lorentz kinematics - DM 17/5/2010
//        VX_PREF = VX_PREF - VX_BU_SUM
//        VY_PREF = VY_PREF - VY_BU_SUM
//        VZ_PREF = VZ_PREF - VZ_BU_SUM
// Lorentz transformation
        lorentz_boost(-VX_BU_SUM,-VY_BU_SUM,-VZ_BU_SUM,
                VX_PREF,VY_PREF,VZ_PREF,
                &VXOUT,&VYOUT,&VZOUT);

        VX_PREF = VXOUT;
        VY_PREF = VYOUT;
        VZ_PREF = VZOUT;

        VBU2 = VX_PREF*VX_PREF + VY_PREF*VY_PREF + VZ_PREF*VZ_PREF;
        GAMMA_REL = std::sqrt(1.0 - VBU2 / (C*C));
        ETOT_PRF = aprf * amu / GAMMA_REL;
        PXPRFP = ETOT_PRF * VX_PREF / C;
        PYPRFP = ETOT_PRF * VY_PREF / C;
        PZPRFP = ETOT_PRF * VZ_PREF / C;

        PX_BU_SUM = 0.0;
        PY_BU_SUM = 0.0;
        PZ_BU_SUM = 0.0;

        PX_BU_SUM = PXPRFP;
        PY_BU_SUM = PYPRFP;
        PZ_BU_SUM = PZPRFP;
        E_tot_BU = ETOT_PRF;

        EKIN_BU = aprf * amu / GAMMA_REL - aprf * amu;

        for(I_Breakup=0;I_Breakup<IMULTBU;I_Breakup++){
// Lorentz kinematics - DM 17/5/2010
//         BU_TAB(I_Breakup,5) = BU_TAB(I_Breakup,5) - VX_BU_SUM
//         BU_TAB(I_Breakup,6) = BU_TAB(I_Breakup,6) - VY_BU_SUM
//         BU_TAB(I_Breakup,7) = BU_TAB(I_Breakup,7) - VZ_BU_SUM
// Lorentz transformation
          lorentz_boost(-VX_BU_SUM,-VY_BU_SUM,-VZ_BU_SUM,
                BU_TAB[I_Breakup][4],BU_TAB[I_Breakup][5],BU_TAB[I_Breakup][6],
                &VXOUT,&VYOUT,&VZOUT);

          BU_TAB[I_Breakup][4] = VXOUT;
          BU_TAB[I_Breakup][5] = VYOUT;
          BU_TAB[I_Breakup][6] = VZOUT;

          VBU2 = BU_TAB[I_Breakup][4]*BU_TAB[I_Breakup][4] +
                 BU_TAB[I_Breakup][5]*BU_TAB[I_Breakup][5] +
                 BU_TAB[I_Breakup][6]*BU_TAB[I_Breakup][6];
          GAMMA_REL = std::sqrt(1.0 - VBU2 / (C*C));

          ETOT_BU = BU_TAB[I_Breakup][1]*amu/GAMMA_REL;

          EKIN_BU = EKIN_BU + BU_TAB[I_Breakup][1] * amu /
                    GAMMA_REL - BU_TAB[I_Breakup][1] * amu;

          PX_BU = ETOT_BU * BU_TAB[I_Breakup][4] / C;
          PY_BU = ETOT_BU * BU_TAB[I_Breakup][5] / C;
          PZ_BU = ETOT_BU * BU_TAB[I_Breakup][6] / C;
          E_tot_BU = E_tot_BU + ETOT_BU;

          PX_BU_SUM = PX_BU_SUM + PX_BU;
          PY_BU_SUM = PY_BU_SUM + PY_BU;
          PZ_BU_SUM = PZ_BU_SUM + PZ_BU;
        }// for I_Breakup

        if(std::abs(PX_BU_SUM)>10. || std::abs(PY_BU_SUM)>10. ||
         std::abs(PZ_BU_SUM)>10.){

//   In the frame of source (i.e. prefragment after INCL)
        P_BU_SUM = std::sqrt(PX_BU_SUM*PX_BU_SUM + PY_BU_SUM*PY_BU_SUM +
                  PZ_BU_SUM*PZ_BU_SUM);
// ********************************************************
// PPRFP is in MeV/c
        ETOT_SUM = std::sqrt(P_BU_SUM*P_BU_SUM +
                  (AAINCL * amu)*(AAINCL * amu));

        VX_BU_SUM = C * PX_BU_SUM / ETOT_SUM;
        VY_BU_SUM = C * PY_BU_SUM / ETOT_SUM;
        VZ_BU_SUM = C * PZ_BU_SUM / ETOT_SUM;

// Lorentz kinematics
//        VX_PREF = VX_PREF - VX_BU_SUM
//        VY_PREF = VY_PREF - VY_BU_SUM
//        VZ_PREF = VZ_PREF - VZ_BU_SUM
// Lorentz transformation
        lorentz_boost(-VX_BU_SUM,-VY_BU_SUM,-VZ_BU_SUM,
                VX_PREF,VY_PREF,VZ_PREF,
                &VXOUT,&VYOUT,&VZOUT);

        VX_PREF = VXOUT;
        VY_PREF = VYOUT;
        VZ_PREF = VZOUT;

        VBU2 = VX_PREF*VX_PREF + VY_PREF*VY_PREF + VZ_PREF*VZ_PREF;
        GAMMA_REL = std::sqrt(1.0 - VBU2 / (C*C));
        ETOT_PRF = aprf * amu / GAMMA_REL;
        PXPRFP = ETOT_PRF * VX_PREF / C;
        PYPRFP = ETOT_PRF * VY_PREF / C;
        PZPRFP = ETOT_PRF * VZ_PREF / C;

        PX_BU_SUM = 0.0;
        PY_BU_SUM = 0.0;
        PZ_BU_SUM = 0.0;

        PX_BU_SUM = PXPRFP;
        PY_BU_SUM = PYPRFP;
        PZ_BU_SUM = PZPRFP;
        E_tot_BU = ETOT_PRF;

        EKIN_BU = aprf * amu / GAMMA_REL - aprf * amu;

         for(I_Breakup=0;I_Breakup<IMULTBU;I_Breakup++){
// Lorentz kinematics - DM 17/5/2010
//         BU_TAB(I_Breakup,5) = BU_TAB(I_Breakup,5) - VX_BU_SUM
//         BU_TAB(I_Breakup,6) = BU_TAB(I_Breakup,6) - VY_BU_SUM
//         BU_TAB(I_Breakup,7) = BU_TAB(I_Breakup,7) - VZ_BU_SUM
// Lorentz transformation
          lorentz_boost(-VX_BU_SUM,-VY_BU_SUM,-VZ_BU_SUM,
                BU_TAB[I_Breakup][4],BU_TAB[I_Breakup][5],BU_TAB[I_Breakup][6],
                &VXOUT,&VYOUT,&VZOUT);

          BU_TAB[I_Breakup][4] = VXOUT;
          BU_TAB[I_Breakup][5] = VYOUT;
          BU_TAB[I_Breakup][6] = VZOUT;

          VBU2 = BU_TAB[I_Breakup][4]*BU_TAB[I_Breakup][4] +
                 BU_TAB[I_Breakup][5]*BU_TAB[I_Breakup][5] +
                 BU_TAB[I_Breakup][6]*BU_TAB[I_Breakup][6];
          GAMMA_REL = std::sqrt(1.0 - VBU2 / (C*C));

          ETOT_BU = BU_TAB[I_Breakup][1]*amu/GAMMA_REL;

          EKIN_BU = EKIN_BU + BU_TAB[I_Breakup][1] * amu /
                    GAMMA_REL - BU_TAB[I_Breakup][1] * amu;

          PX_BU = ETOT_BU * BU_TAB[I_Breakup][4] / C;
          PY_BU = ETOT_BU * BU_TAB[I_Breakup][5] / C;
          PZ_BU = ETOT_BU * BU_TAB[I_Breakup][6] / C;
          E_tot_BU = E_tot_BU + ETOT_BU;

          PX_BU_SUM = PX_BU_SUM + PX_BU;
          PY_BU_SUM = PY_BU_SUM + PY_BU;
          PZ_BU_SUM = PZ_BU_SUM + PZ_BU;
         }// for I_Breakup
        }// if DABS(PX_BU_SUM).GT.10.d0
//
//      Find the limits that fragment is bound - only done for neutrons and LCPs and for
//      nuclei with A=Z, for other nuclei it will be done after decay:

         INEWLOOP = 0;
         for(G4int i=0;i<IMULTBU;i++){
          if(BU_TAB[i][0]<3.0 || BU_TAB[i][0]==BU_TAB[i][1]){
            unstable_nuclei(idnint(BU_TAB[i][1]),idnint(BU_TAB[i][0]), &afpnew,&zfpnew,IOUNSTABLE,
            BU_TAB[i][4], BU_TAB[i][5], BU_TAB[i][6],
            &VP1X,&VP1Y,&VP1Z,BU_TAB_TEMP,&ILOOP);

            if(IOUNSTABLE>0){
// Properties of "heavy fragment":
             BU_TAB[i][1] = G4double(afpnew);
             BU_TAB[i][0] = G4double(zfpnew);
             BU_TAB[i][4] = VP1X;
             BU_TAB[i][5] = VP1Y;
             BU_TAB[i][6] = VP1Z;

//Properties of "light" fragments:
             for(int IJ=0;IJ<ILOOP;IJ++){
              BU_TAB[IMULTBU+INEWLOOP+IJ][0] = BU_TAB_TEMP[IJ][0];
              BU_TAB[IMULTBU+INEWLOOP+IJ][1] = BU_TAB_TEMP[IJ][1];
              BU_TAB[IMULTBU+INEWLOOP+IJ][4] = BU_TAB_TEMP[IJ][2];
              BU_TAB[IMULTBU+INEWLOOP+IJ][5] = BU_TAB_TEMP[IJ][3];
              BU_TAB[IMULTBU+INEWLOOP+IJ][6] = BU_TAB_TEMP[IJ][4];
              BU_TAB[IMULTBU+INEWLOOP+IJ][2] = 0.0;
              BU_TAB[IMULTBU+INEWLOOP+IJ][3] = 0.0;
             }// for ILOOP

             INEWLOOP = INEWLOOP + ILOOP;

            }// if IOUNSTABLE.GT.0
          }//if BU_TAB[I_Breakup][0]<3.0
         }// for IMULTBU

// Increased array of BU_TAB
        IMULTBU = IMULTBU + INEWLOOP;
// Evaporation from multifragmentation products
        opt->optimfallowed = 1;  //  IMF is allowed
        fiss->ifis = 0;          //  fission is not allowed
        gammaemission=0;
        ILOOPBU = 0;

//  Arrays for lambda emission from breakup fragments
         G4double * problamb;
         problamb = new G4double[IMULTBU];
         G4double sumN = aprf - zprf;
         for(G4int i=0;i<IMULTBU;i++)sumN=sumN+BU_TAB[i][1]-BU_TAB[i][0];

         for(G4int i=0;i<IMULTBU;i++){
         problamb[i] = (BU_TAB[i][1]-BU_TAB[i][0])/sumN;
         }
         G4int * Nblamb;
         Nblamb = new G4int[IMULTBU];
         for(G4int i=0;i<IMULTBU;i++)Nblamb[i] = 0;
         for(G4int j=0;j<NbLam0;){
          G4double probtotal = (aprf - zprf)/sumN;
          G4double ran =  G4AblaRandom::flat();
//   Lambdas in the heavy breakup fragment
          if(ran <= probtotal){
           NbLamprf++;
           goto directlamb0;
          }
          for(G4int i=0;i<IMULTBU;i++){
//   Lambdas in the light breakup residues
           if(probtotal < ran && ran <= probtotal+problamb[i]){
           Nblamb[i] = Nblamb[i] + 1;
           goto directlamb0;
           }
           probtotal = probtotal + problamb[i];
          }
          directlamb0:
          j++;
         }
//
         for(G4int i=0;i<IMULTBU;i++){
          EEBU = BU_TAB[i][2];
          BU_TAB[i][10] = BU_TAB[i][6];
          G4double jprfbu = BU_TAB[i][9];
          if(BU_TAB[i][0]>2.0){
           G4int nbl = Nblamb[i];
           evapora(BU_TAB[i][0],BU_TAB[i][1],&EEBU,0.0, &ZFBU, &AFBU, &mtota, &vz_evabu, &vx_evabu,&vy_evabu, &ff, &fimf, &ZIMFBU, &AIMFBU,&TKEIMFBU, &jprfbu, &inttype, &inum,EV_TEMP,&IEV_TAB_TEMP,&nbl);

           Nblamb[i] = nbl;
           BU_TAB[i][9] = jprfbu;

//Velocities of evaporated particles (in the frame of the primary prefragment)
               for(G4int IJ = 0; IJ< IEV_TAB_TEMP;IJ++){ 
               EV_TAB[IJ+IEV_TAB][0] = EV_TEMP[IJ][0];
               EV_TAB[IJ+IEV_TAB][1] = EV_TEMP[IJ][1];
               EV_TAB[IJ+IEV_TAB][5] = EV_TEMP[IJ][5];
//Lorentz kinematics
//                 DO IK = 3, 5, 1
//                 EV_TAB(IJ+IEV_TAB,IK) = EV_TEMP(IJ,IK) + BU_TAB(I,IK+2)
//                 ENDDO
// Lorentz transformation
               lorentz_boost(BU_TAB[i][4],BU_TAB[i][5],BU_TAB[i][6],
                EV_TEMP[IJ][2],EV_TEMP[IJ][3],EV_TEMP[IJ][4],
                &VXOUT,&VYOUT,&VZOUT);
               EV_TAB[IJ+IEV_TAB][2] = VXOUT;
               EV_TAB[IJ+IEV_TAB][3] = VYOUT;
               EV_TAB[IJ+IEV_TAB][4] = VZOUT;
               }
               IEV_TAB = IEV_TAB + IEV_TAB_TEMP;

//All velocities in the frame of the "primary" prefragment (after INC)
// Lorentz kinematics
//                BU_TAB(I,5) = BU_TAB(I,5) + VX_EVABU
//                BU_TAB(I,6) = BU_TAB(I,6) + VY_EVABU
//                BU_TAB(I,7) = BU_TAB(I,7) + VZ_EVABU
// Lorentz transformation
               lorentz_boost(vx_evabu,vy_evabu,vz_evabu,
                BU_TAB[i][4],BU_TAB[i][5],BU_TAB[i][6],
                &VXOUT,&VYOUT,&VZOUT);
               BU_TAB[i][4] = VXOUT;
               BU_TAB[i][5] = VYOUT;
               BU_TAB[i][6] = VZOUT;

               if(fimf==0){
                 BU_TAB[i][7] = dint(ZFBU);
                 BU_TAB[i][8] = dint(AFBU);
                 BU_TAB[i][11]= nbl;
               }// if fimf==0

               if(fimf==1){
//            PRINT*,'IMF EMISSION FROM BU PRODUCTS'
// IMF emission: Heavy partner is not allowed to fission or to emitt IMF.
               //double FEE = EEBU;
               G4int FFBU1 = 0;
               G4int FIMFBU1 = 0;
               opt->optimfallowed = 0;  //  IMF is not allowed
               fiss->ifis = 0;          //  fission is not allowed
// Velocities of IMF and partner: 1 denotes partner, 2 denotes IMF
               G4double EkinR1 = TKEIMFBU * AIMFBU / (AFBU+AIMFBU);
               G4double EkinR2 = TKEIMFBU * AFBU / (AFBU+AIMFBU);
               G4double V1 = std::sqrt(EkinR1/AFBU) * 1.3887;
               G4double V2 = std::sqrt(EkinR2/AIMFBU) * 1.3887;
               G4double VZ1_IMF = (2.0 * G4AblaRandom::flat() - 1.0) * V1;
               G4double VPERP1 = std::sqrt(V1*V1 - VZ1_IMF*VZ1_IMF);
               G4double ALPHA1 = G4AblaRandom::flat() * 2. * 3.142;
               G4double VX1_IMF = VPERP1 * std::sin(ALPHA1);
               G4double VY1_IMF = VPERP1 * std::cos(ALPHA1);
               G4double VX2_IMF = - VX1_IMF / V1 * V2;
               G4double VY2_IMF = - VY1_IMF / V1 * V2;
               G4double VZ2_IMF = - VZ1_IMF / V1 * V2;

               G4double EEIMFP = EEBU * AFBU /(AFBU + AIMFBU);
               G4double EEIMF = EEBU * AIMFBU /(AFBU + AIMFBU);

// Decay of heavy partner
     G4double IINERTTOT = 0.40 * 931.490 * 1.160*1.160 *( std::pow(AIMFBU,5.0/3.0) + std::pow(AFBU,5.0/3.0)) + 931.490 * 1.160*1.160*AIMFBU*AFBU/(AIMFBU+AFBU)*(std::pow(AIMFBU,1./3.) + std::pow(AFBU,1./3.))*(std::pow(AIMFBU,1./3.) + std::pow(AFBU,1./3.));

     G4double JPRFHEAVY = BU_TAB[i][9] * 0.4 * 931.49 * 1.16*1.16 * std::pow(AFBU,5.0/3.0) / IINERTTOT;
     G4double JPRFLIGHT = BU_TAB[i][9] * 0.4 * 931.49 * 1.16*1.16 * std::pow(AIMFBU,5.0/3.0) / IINERTTOT;

// Lorentz kinematics
//           BU_TAB(I,5) = BU_TAB(I,5) + VX1_IMF
//           BU_TAB(I,6) = BU_TAB(I,6) + VY1_IMF
//           BU_TAB(I,7) = BU_TAB(I,7) + VZ1_IMF
// Lorentz transformation
               lorentz_boost(VX1_IMF,VY1_IMF,VZ1_IMF,
                BU_TAB[i][4],BU_TAB[i][5],BU_TAB[i][6],
                &VXOUT,&VYOUT,&VZOUT);
               BU_TAB[i][4] = VXOUT;
               BU_TAB[i][5] = VYOUT;
               BU_TAB[i][6] = VZOUT;

     G4double vx1ev_imf=0., vy1ev_imf=0., vz1ev_imf=0., zdummy=0., adummy=0., tkedummy=0.,jprf1=0.;

 //  Lambda particles 
      G4int NbLamH=0;
      G4int NbLamimf=0;
      G4double pbH = (AFBU-ZFBU) / (AFBU-ZFBU+AIMFBU-ZIMFBU); 
      for(G4int j=0;j<nbl;j++){
       if(G4AblaRandom::flat()<pbH){
        NbLamH++;
       }else{
        NbLamimf++;
       }
      }
// Decay of IMF's partner:
               evapora(ZFBU,AFBU,&EEIMFP,JPRFHEAVY, &ZFFBU, &AFFBU, &mtota,  &vz1ev_imf, &vx1ev_imf,&vy1ev_imf, &FFBU1, &FIMFBU1, &zdummy, &adummy,&tkedummy, &jprf1, &inttype, &inum,EV_TEMP,&IEV_TAB_TEMP,&NbLamH);

               for(G4int IJ = 0; IJ< IEV_TAB_TEMP;IJ++){ 
               EV_TAB[IJ+IEV_TAB][0] = EV_TEMP[IJ][0];
               EV_TAB[IJ+IEV_TAB][1] = EV_TEMP[IJ][1];
               EV_TAB[IJ+IEV_TAB][5] = EV_TEMP[IJ][5];
//Lorentz kinematics
//                 DO IK = 3, 5, 1
//                 EV_TAB(IJ+IEV_TAB,IK) = EV_TEMP(IJ,IK) + BU_TAB(I,IK+2)
//                 ENDDO
// Lorentz transformation
               lorentz_boost(BU_TAB[i][4],BU_TAB[i][5],BU_TAB[i][6],
                EV_TEMP[IJ][2],EV_TEMP[IJ][3],EV_TEMP[IJ][4],
                &VXOUT,&VYOUT,&VZOUT);
               EV_TAB[IJ+IEV_TAB][2] = VXOUT;
               EV_TAB[IJ+IEV_TAB][3] = VYOUT;
               EV_TAB[IJ+IEV_TAB][4] = VZOUT;
               }
               IEV_TAB = IEV_TAB + IEV_TAB_TEMP;

               BU_TAB[i][7] = dint(ZFFBU);
               BU_TAB[i][8] = dint(AFFBU);
               BU_TAB[i][11]= NbLamH;
//Lorentz kinematics
//           BU_TAB(I,5) = BU_TAB(I,5) + vx1ev_imf
//           BU_TAB(I,6) = BU_TAB(I,6) + vy1ev_imf
//           BU_TAB(I,7) = BU_TAB(I,7) + vz1ev_imf
               lorentz_boost(vx1ev_imf,vy1ev_imf,vz1ev_imf,
                BU_TAB[i][4],BU_TAB[i][5],BU_TAB[i][6],
                &VXOUT,&VYOUT,&VZOUT);
               BU_TAB[i][4] = VXOUT;
               BU_TAB[i][5] = VYOUT;
               BU_TAB[i][6] = VZOUT;
// For IMF - fission and IMF emission are not allowed
              G4int FFBU2 = 0;
              G4int FIMFBU2 = 0;
              opt->optimfallowed = 0;  //  IMF is not allowed
              fiss->ifis = 0;          //  fission is not allowed
// Decay of IMF
              G4double zffimf, affimf,zdummy1, adummy1, tkedummy1, jprf2, vx2ev_imf, vy2ev_imf, vz2ev_imf;

              evapora(ZIMFBU,AIMFBU,&EEIMF,JPRFLIGHT, &zffimf, &affimf, &mtota, &vz2ev_imf, &vx2ev_imf,&vy2ev_imf, &FFBU2, &FIMFBU2, &zdummy1, &adummy1,&tkedummy1, &jprf2, &inttype, &inum,EV_TEMP,&IEV_TAB_TEMP,&NbLamimf);

               for(G4int IJ = 0; IJ< IEV_TAB_TEMP;IJ++){ 
               EV_TAB[IJ+IEV_TAB][0] = EV_TEMP[IJ][0];
               EV_TAB[IJ+IEV_TAB][1] = EV_TEMP[IJ][1];
               EV_TAB[IJ+IEV_TAB][5] = EV_TEMP[IJ][5];
//Lorentz kinematics
//            EV_TAB(IJ+IEV_TAB,3) = EV_TEMP(IJ,3) + BU_TAB(I,5) +VX2_IMF
//            EV_TAB(IJ+IEV_TAB,4) = EV_TEMP(IJ,4) + BU_TAB(I,6) +VY2_IMF
//            EV_TAB(IJ+IEV_TAB,5) = EV_TEMP(IJ,5) + BU_TAB(I,7) +VZ2_IMF
// Lorentz transformation
               lorentz_boost(BU_TAB[i][4],BU_TAB[i][5],BU_TAB[i][6],
                EV_TEMP[IJ][2],EV_TEMP[IJ][3],EV_TEMP[IJ][4],
                &VXOUT,&VYOUT,&VZOUT);
               lorentz_boost(VX2_IMF,VY2_IMF,VZ2_IMF,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
               EV_TAB[IJ+IEV_TAB][2] = VX2OUT;
               EV_TAB[IJ+IEV_TAB][3] = VY2OUT;
               EV_TAB[IJ+IEV_TAB][4] = VZ2OUT;
               }
               IEV_TAB = IEV_TAB + IEV_TAB_TEMP;

               BU_TAB[IMULTBU+ILOOPBU][0] = BU_TAB[i][0];
               BU_TAB[IMULTBU+ILOOPBU][1] = BU_TAB[i][1];
               BU_TAB[IMULTBU+ILOOPBU][2] = BU_TAB[i][2];
               BU_TAB[IMULTBU+ILOOPBU][3] = BU_TAB[i][3];
               BU_TAB[IMULTBU+ILOOPBU][7] = dint(zffimf);
               BU_TAB[IMULTBU+ILOOPBU][8] = dint(affimf);
               BU_TAB[IMULTBU+ILOOPBU][11]= NbLamimf;
// Lorentz transformation
               lorentz_boost(VX2_IMF,VY2_IMF,VZ2_IMF,
                BU_TAB[i][4],BU_TAB[i][5],BU_TAB[i][6],
                &VXOUT,&VYOUT,&VZOUT);
               lorentz_boost(vx2ev_imf,vy2ev_imf,vz2ev_imf,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
               BU_TAB[IMULTBU+ILOOPBU][4] = VX2OUT;
               BU_TAB[IMULTBU+ILOOPBU][5] = VY2OUT;
               BU_TAB[IMULTBU+ILOOPBU][6] = VZ2OUT;
               ILOOPBU = ILOOPBU + 1;
               }// if fimf==1

          } else {// if BU_TAB(I,1).GT.2.D0
	    //BU_TAB[i][0] = BU_TAB[i][0];
	    //BU_TAB[i][1] = BU_TAB[i][1];
	    //BU_TAB[i][2] = BU_TAB[i][2];
	    //BU_TAB[i][3] = BU_TAB[i][3];
           BU_TAB[i][7] = BU_TAB[i][0];
           BU_TAB[i][8] = BU_TAB[i][1];
           //BU_TAB[i][4] = BU_TAB[i][4];
           //BU_TAB[i][5] = BU_TAB[i][5];
           //BU_TAB[i][6] = BU_TAB[i][6];
           BU_TAB[i][11]= Nblamb[i];
          }// if BU_TAB(I,1).GT.2.D0
         }// for IMULTBU

         IMULTBU = IMULTBU + ILOOPBU;
//
// RESOLVE UNSTABLE NUCLEI
//
      INEWLOOP = 0;
      ABU_SUM = 0.0;
      ZBU_SUM = 0.0;
//
      for(G4int i=0;i<IMULTBU;i++){
       ABU_SUM = ABU_SUM + BU_TAB[i][8];
       ZBU_SUM = ZBU_SUM + BU_TAB[i][7];
       unstable_nuclei(idnint(BU_TAB[i][8]),idnint(BU_TAB[i][7]), &afpnew,&zfpnew,IOUNSTABLE,
            BU_TAB[i][4], BU_TAB[i][5], BU_TAB[i][6],
            &VP1X,&VP1Y,&VP1Z,BU_TAB_TEMP1,&ILOOP);

//From now on, all neutrons and LCP created in above subroutine are part of the
// BU_TAB array (see below - Properties of "light" fragments). Therefore,
// NEVA, PEVA ... are not needed any more in the break-up stage.

         if(IOUNSTABLE>0){
// Properties of "heavy fragment":
            ABU_SUM = ABU_SUM + G4double(afpnew) - BU_TAB[i][8];
            ZBU_SUM = ZBU_SUM + G4double(zfpnew) - BU_TAB[i][7];
             BU_TAB[i][8] = G4double(afpnew);
             BU_TAB[i][7] = G4double(zfpnew);
             BU_TAB[i][4] = VP1X;
             BU_TAB[i][5] = VP1Y;
             BU_TAB[i][6] = VP1Z;

//Properties of "light" fragments:
             for(G4int IJ=0;IJ<ILOOP;IJ++){
              BU_TAB[IMULTBU+INEWLOOP+IJ][7] = BU_TAB_TEMP1[IJ][0];
              BU_TAB[IMULTBU+INEWLOOP+IJ][8] = BU_TAB_TEMP1[IJ][1];
              BU_TAB[IMULTBU+INEWLOOP+IJ][4] = BU_TAB_TEMP1[IJ][2];
              BU_TAB[IMULTBU+INEWLOOP+IJ][5] = BU_TAB_TEMP1[IJ][3];
              BU_TAB[IMULTBU+INEWLOOP+IJ][6] = BU_TAB_TEMP1[IJ][4];
              BU_TAB[IMULTBU+INEWLOOP+IJ][2] = 0.0;
              BU_TAB[IMULTBU+INEWLOOP+IJ][3] = 0.0;
              BU_TAB[IMULTBU+INEWLOOP+IJ][0] = BU_TAB[i][0];
              BU_TAB[IMULTBU+INEWLOOP+IJ][1] = BU_TAB[i][1];
              BU_TAB[IMULTBU+INEWLOOP+IJ][11] = BU_TAB[i][11];
              ABU_SUM = ABU_SUM + BU_TAB[IMULTBU+INEWLOOP+IJ][8];
              ZBU_SUM = ZBU_SUM + BU_TAB[IMULTBU+INEWLOOP+IJ][7];
             }// for ILOOP

             INEWLOOP = INEWLOOP + ILOOP;
         }// if(IOUNSTABLE>0)
      }// for IMULTBU unstable

// Increased array of BU_TAB
        IMULTBU = IMULTBU + INEWLOOP;

// Transform all velocities into the rest frame of the projectile
        lorentz_boost(VX_incl,VY_incl,VZ_incl,
                VX_PREF,VY_PREF,VZ_PREF,
                &VXOUT,&VYOUT,&VZOUT);
        VX_PREF = VXOUT;
        VY_PREF = VYOUT;
        VZ_PREF = VZOUT;

        for(G4int i=0;i<IMULTBU;i++){
         lorentz_boost(VX_incl,VY_incl,VZ_incl,
                BU_TAB[i][4],BU_TAB[i][5],BU_TAB[i][6],
                &VXOUT,&VYOUT,&VZOUT);
         BU_TAB[i][4] = VXOUT;
         BU_TAB[i][5] = VYOUT;
         BU_TAB[i][6] = VZOUT;
        }
        for(G4int i=0;i<IEV_TAB;i++){
         lorentz_boost(VX_incl,VY_incl,VZ_incl,
                EV_TAB[i][2],EV_TAB[i][3],EV_TAB[i][4],
                &VXOUT,&VYOUT,&VZOUT);
         EV_TAB[i][2] = VXOUT;
         EV_TAB[i][3] = VYOUT;
         EV_TAB[i][4] = VZOUT;
        }
        if(IMULTBU>200)std::cout << "IMULTBU>200 " << IMULTBU << std::endl;
        delete[] problamb;
        delete[] Nblamb;
        }// if(T_diff>0.1)
        // End of multi-fragmentation
      mult7777:

// Start basic de-excitation of fragments
      aprfp = idnint(aprf);
      zprfp = idnint(zprf);
     
      if(IMULTIFR == 0){
// These momenta are in the frame of the projectile (or target in case of direct kinematics)
       VX_PREF = VX_incl;
       VY_PREF = VY_incl;
       VZ_PREF = VZ_incl;
      }
// Lambdas after multi-fragmentation
      if(IMULTIFR == 1){
       NbLam0 = NbLamprf;
      }
//
// CALL THE EVAPORATION SUBROUTINE
//
      opt->optimfallowed = 1; //  IMF is allowed
      fiss->ifis = 1;         //  fission is allowed
      fimf=0;
      ff=0;

// To spare computing time; these events in any case cannot decay
//      IF(ZPRFP.LE.2.AND.ZPRFP.LT.APRFP)THEN FIXME: <= or <
      if(zprfp<=2 && zprfp<aprfp){
       zf = zprf;
       af = aprf;
       ee = 0.0;
       ff = 0;
       fimf = 0;
       ftype = 0;
       aimf = 0.0;
       zimf = 0.0;
       tkeimf = 0.0;
       vx_eva = 0.0;
       vy_eva = 0.0;
       vz_eva = 0.0;
       jprf0 = jprf;
       goto a1972;
       }

//      if(ZPRFP.LE.2.AND.ZPRFP.EQ.APRFP)
      if(zprfp<=2 && zprfp==aprfp){
       unstable_nuclei(aprfp,zprfp,&afpnew,&zfpnew,IOUNSTABLE,
       VX_PREF, VY_PREF, VZ_PREF,
       &VP1X,&VP1Y,&VP1Z,EV_TAB_TEMP,&ILOOP);
         af = G4double(afpnew);
         zf = G4double(zfpnew);
         VX_PREF = VP1X;
         VY_PREF = VP1Y;
         VZ_PREF = VP1Z;
         for(G4int I = 0;I<ILOOP;I++){
          for(G4int IJ = 0; IJ<6; IJ++)
           EV_TAB[I+IEV_TAB][IJ] = EV_TAB_TEMP[I][IJ];
         }
        IEV_TAB = IEV_TAB + ILOOP;
        ee = 0.0;
        ff = 0;
        fimf = 0;
        ftype = 0;
        aimf = 0.0;
        zimf = 0.0;
        tkeimf = 0.0;
        vx_eva = 0.0;
        vy_eva = 0.0;
        vz_eva = 0.0;
        jprf0 = jprf;
       goto a1972;
       }

//      IF(ZPRFP.EQ.APRFP)THEN
      if(zprfp==aprfp){
       unstable_nuclei(aprfp,zprfp,&afpnew,&zfpnew,IOUNSTABLE,
       VX_PREF, VY_PREF, VZ_PREF,
       &VP1X,&VP1Y,&VP1Z,EV_TAB_TEMP,&ILOOP);
         af = G4double(afpnew);
         zf = G4double(zfpnew);
         VX_PREF = VP1X;
         VY_PREF = VP1Y;
         VZ_PREF = VP1Z;
         for(G4int I = 0;I<ILOOP;I++){
          for(G4int IJ = 0; IJ<6; IJ++)
           EV_TAB[I+IEV_TAB][IJ] = EV_TAB_TEMP[I][IJ];
         }
        IEV_TAB = IEV_TAB + ILOOP;
        ee = 0.0;
        ff = 0;
        fimf = 0;
        ftype = 0;
        aimf = 0.0;
        zimf = 0.0;
        tkeimf = 0.0;
        vx_eva = 0.0;
        vy_eva = 0.0;
        vz_eva = 0.0;
        jprf0 = jprf;
       goto a1972;
      }
//
      evapora(zprf,aprf,&ee,jprf, &zf, &af, &mtota, &vz_eva, &vx_eva, &vy_eva, &ff, &fimf, &zimf, &aimf,&tkeimf, &jprf0, &inttype, &inum,EV_TEMP,&IEV_TAB_TEMP,&NbLam0);
//
               for(G4int IJ = 0; IJ< IEV_TAB_TEMP;IJ++){ 
               EV_TAB[IJ+IEV_TAB][0] = EV_TEMP[IJ][0];
               EV_TAB[IJ+IEV_TAB][1] = EV_TEMP[IJ][1];
               EV_TAB[IJ+IEV_TAB][5] = EV_TEMP[IJ][5];
//
//               EV_TAB(IJ+IEV_TAB,3) = EV_TEMP(IJ,3) + VX_PREF
//               EV_TAB(IJ+IEV_TAB,4) = EV_TEMP(IJ,4) + VY_PREF
//               EV_TAB(IJ+IEV_TAB,5) = EV_TEMP(IJ,5) + VZ_PREF
// Lorentz transformation
               lorentz_boost(VX_PREF,VY_PREF,VZ_PREF,
                EV_TEMP[IJ][2],EV_TEMP[IJ][3],EV_TEMP[IJ][4],
                &VXOUT,&VYOUT,&VZOUT);
               EV_TAB[IJ+IEV_TAB][2] = VXOUT;
               EV_TAB[IJ+IEV_TAB][3] = VYOUT;
               EV_TAB[IJ+IEV_TAB][4] = VZOUT;
               }
               IEV_TAB = IEV_TAB + IEV_TAB_TEMP;

      a1972:

// vi_pref - velocity of the prefragment; vi_eva - recoil due to evaporation
      lorentz_boost(VX_PREF,VY_PREF,VZ_PREF,
      vx_eva,vy_eva,vz_eva,
      &VXOUT,&VYOUT,&VZOUT);
      V_CM[0] = VXOUT;
      V_CM[1] = VYOUT;
      V_CM[2] = VZOUT;
//
      if(ff == 0 && fimf == 0){
// Evaporation of neutrons and LCP; no IMF, no fission
      ftype = 0;
      ZFP1 = idnint(zf);
      AFP1 = idnint(af);
      SFP1 = NbLam0;
      AFPIMF = 0;
      ZFPIMF = 0;
      SFPIMF = 0;
      ZFP2 = 0;
      AFP2 = 0;
      SFP2 = 0;
      VFP1_CM[0] = V_CM[0];
      VFP1_CM[1] = V_CM[1];
      VFP1_CM[2] = V_CM[2];
       for(G4int j=0;j<3;j++){
        VIMF_CM[j] = 0.0;
        VFP2_CM[j] = 0.0;
       }
      }
//
      if(ff == 1 && fimf == 0) ftype = 1;     // fission
      if(ff == 0 && fimf == 1) ftype = 2;     // IMF emission
//
// AFP,ZFP IS THE FINAL FRAGMENT IF NO FISSION OR IMF EMISSION OCCURS
// IN CASE OF FISSION IT IS THE NUCLEUS THAT UNDERGOES FISSION OR IMF
//

//***************** FISSION ***************************************
//
    if(ftype == 1){
    varntp->kfis = 1;
    if(NbLam0>0)varntp->kfis = 20;
   //   ftype1=0;

      G4int IEV_TAB_FIS = 0,imode=0;

      G4double vx1_fission=0.,vy1_fission=0.,vz1_fission=0.;
      G4double vx2_fission=0.,vy2_fission=0.,vz2_fission=0.;
      G4double vx_eva_sc=0.,vy_eva_sc=0.,vz_eva_sc=0.;

      fission(af,zf,ee,jprf0,
          &vx1_fission,&vy1_fission,&vz1_fission,
          &vx2_fission,&vy2_fission,&vz2_fission,
          &ZFP1,&AFP1,&SFP1,&ZFP2,&AFP2,&SFP2,&imode,
          &vx_eva_sc,&vy_eva_sc,&vz_eva_sc,EV_TEMP,&IEV_TAB_FIS,&NbLam0);

               for(G4int IJ = 0; IJ< IEV_TAB_FIS;IJ++){ 
               EV_TAB[IJ+IEV_TAB][0] = EV_TEMP[IJ][0];
               EV_TAB[IJ+IEV_TAB][1] = EV_TEMP[IJ][1];
               EV_TAB[IJ+IEV_TAB][5] = EV_TEMP[IJ][5];
// Lorentz kinematics
//               EV_TAB(IJ+IEV_TAB,3) = EV_TEMP(IJ,3) + VX_PREF
//               EV_TAB(IJ+IEV_TAB,4) = EV_TEMP(IJ,4) + VY_PREF
//               EV_TAB(IJ+IEV_TAB,5) = EV_TEMP(IJ,5) + VZ_PREF
// Lorentz transformation
               lorentz_boost(V_CM[0],V_CM[1],V_CM[2],
                EV_TEMP[IJ][2],EV_TEMP[IJ][3],EV_TEMP[IJ][4],
                &VXOUT,&VYOUT,&VZOUT);
               EV_TAB[IJ+IEV_TAB][2] = VXOUT;
               EV_TAB[IJ+IEV_TAB][3] = VYOUT;
               EV_TAB[IJ+IEV_TAB][4] = VZOUT;
               }
               IEV_TAB = IEV_TAB + IEV_TAB_FIS;

    //  if(imode==1) ftype1 = 1;    // S1 mode
    //  if(imode==2) ftype1 = 2;    // S2 mode

      AFPIMF = 0;
      ZFPIMF = 0;
      SFPIMF = 0;

// VX_EVA_SC,VY_EVA_SC,VZ_EVA_SC - recoil due to particle emisison
// between saddle and scission
// Lorentz kinematics
//        VFP1_CM(1) = V_CM(1) + VX1_FISSION + VX_EVA_SC ! Velocity of FF1 in x
//        VFP1_CM(2) = V_CM(2) + VY1_FISSION + VY_EVA_SC ! Velocity of FF1 in y
//        VFP1_CM(3) = V_CM(3) + VZ1_FISSION + VZ_EVA_SC ! Velocity of FF1 in x
        lorentz_boost(vx1_fission,vy1_fission,vz1_fission,
                V_CM[0],V_CM[1],V_CM[2],
                &VXOUT,&VYOUT,&VZOUT);
        lorentz_boost(vx_eva_sc,vy_eva_sc,vz_eva_sc,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
        VFP1_CM[0] = VX2OUT;
        VFP1_CM[1] = VY2OUT;
        VFP1_CM[2] = VZ2OUT;

// Lorentz kinematics
//        VFP2_CM(1) = V_CM(1) + VX2_FISSION + VX_EVA_SC ! Velocity of FF2 in x
//        VFP2_CM(2) = V_CM(2) + VY2_FISSION + VY_EVA_SC ! Velocity of FF2 in y
//        VFP2_CM(3) = V_CM(3) + VZ2_FISSION + VZ_EVA_SC ! Velocity of FF2 in x
        lorentz_boost(vx2_fission,vy2_fission,vz2_fission,
                V_CM[0],V_CM[1],V_CM[2],
                &VXOUT,&VYOUT,&VZOUT);
        lorentz_boost(vx_eva_sc,vy_eva_sc,vz_eva_sc,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
        VFP2_CM[0] = VX2OUT;
        VFP2_CM[1] = VY2OUT;
        VFP2_CM[2] = VZ2OUT;

//************** IMF EMISSION ************************************************
//
    }else if(ftype == 2){
// IMF emission: Heavy partner is allowed to fission and to emitt IMF, but ONLY once.
      G4int FF11 = 0;
      G4int FIMF11 = 0;
      opt->optimfallowed = 1;  //  IMF is allowed
      fiss->ifis = 1;          //  fission is allowed
//  Lambda particles 
      G4int NbLamH=0;
      G4int NbLamimf=0;
      G4double pbH = (af-zf) / (af-zf+aimf-zimf);
      //double pbL = aimf / (af+aimf);  
      for(G4int i=0;i<NbLam0;i++){
       if(G4AblaRandom::flat()<pbH){
        NbLamH++;
       }else{
        NbLamimf++;
       }
      }
//
//  Velocities of IMF and partner: 1 denotes partner, 2 denotes IMF
      G4double EkinR1 = tkeimf * aimf / (af+aimf);
      G4double EkinR2 = tkeimf * af / (af+aimf);
      G4double V1 = std::sqrt(EkinR1/af) * 1.3887;
      G4double V2 = std::sqrt(EkinR2/aimf) * 1.3887;
      G4double VZ1_IMF = (2.0 * G4AblaRandom::flat() - 1.0) * V1;
      G4double VPERP1 = std::sqrt(V1*V1 - VZ1_IMF*VZ1_IMF);
      G4double ALPHA1 = G4AblaRandom::flat() * 2. * 3.142;
      G4double VX1_IMF = VPERP1 * std::sin(ALPHA1);
      G4double VY1_IMF = VPERP1 * std::cos(ALPHA1);
      G4double VX2_IMF = - VX1_IMF / V1 * V2;
      G4double VY2_IMF = - VY1_IMF / V1 * V2;
      G4double VZ2_IMF = - VZ1_IMF / V1 * V2;

      G4double EEIMFP = ee * af /(af + aimf);
      G4double EEIMF = ee * aimf /(af + aimf);

// Decay of heavy partner
     G4double IINERTTOT = 0.40 * 931.490 * 1.160*1.160 *( std::pow(aimf,5.0/3.0) + std::pow(af,5.0/3.0)) + 931.490 * 1.160*1.160*aimf*af/(aimf+af)*(std::pow(aimf,1./3.) + std::pow(af,1./3.))*(std::pow(aimf,1./3.) + std::pow(af,1./3.));

     G4double JPRFHEAVY = jprf0 * 0.4 * 931.49 * 1.16*1.16 * std::pow(af,5.0/3.0) / IINERTTOT;
     G4double JPRFLIGHT = jprf0 * 0.4 * 931.49 * 1.16*1.16 * std::pow(aimf,5.0/3.0) / IINERTTOT;
     if(af<2.0) std::cout << "RN117-4,AF,ZF,EE,JPRFheavy" << std::endl;

     G4double vx1ev_imf=0., vy1ev_imf=0., vz1ev_imf=0., zdummy=0., adummy=0., tkedummy=0.,jprf1=0.;

     evapora(zf,af,&EEIMFP,JPRFHEAVY, &zff, &aff, &mtota, &vz1ev_imf, &vx1ev_imf,&vy1ev_imf, &FF11, &FIMF11, &zdummy, &adummy,&tkedummy, &jprf1, &inttype, &inum,EV_TEMP,&IEV_TAB_TEMP,&NbLamH);

               for(G4int IJ = 0; IJ< IEV_TAB_TEMP;IJ++){ 
               EV_TAB[IJ+IEV_TAB][0] = EV_TEMP[IJ][0];
               EV_TAB[IJ+IEV_TAB][1] = EV_TEMP[IJ][1];
               EV_TAB[IJ+IEV_TAB][5] = EV_TEMP[IJ][5];
//
//               EV_TAB(IJ+IEV_TAB,3) = EV_TEMP(IJ,3) + VX_PREF
//               EV_TAB(IJ+IEV_TAB,4) = EV_TEMP(IJ,4) + VY_PREF
//               EV_TAB(IJ+IEV_TAB,5) = EV_TEMP(IJ,5) + VZ_PREF
// Lorentz transformation
               lorentz_boost(V_CM[0],V_CM[1],V_CM[2],
                EV_TEMP[IJ][2],EV_TEMP[IJ][3],EV_TEMP[IJ][4],
                &VXOUT,&VYOUT,&VZOUT);
               lorentz_boost(vx1ev_imf,vy1ev_imf,vz1ev_imf,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
               EV_TAB[IJ+IEV_TAB][2] = VX2OUT;
               EV_TAB[IJ+IEV_TAB][3] = VY2OUT;
               EV_TAB[IJ+IEV_TAB][4] = VZ2OUT;
               }
               IEV_TAB = IEV_TAB + IEV_TAB_TEMP;

// For IMF - fission and IMF emission are not allowed
     G4int FF22 = 0;
     G4int FIMF22 = 0;
      opt->optimfallowed = 0; //  IMF is not allowed
      fiss->ifis = 0;         //  fission is not allowed

// Decay of IMF
     G4double zffimf, affimf,zdummy1=0., adummy1=0., tkedummy1=0.,jprf2,vx2ev_imf,vy2ev_imf,
    vz2ev_imf;

     evapora(zimf,aimf,&EEIMF,JPRFLIGHT, &zffimf, &affimf, &mtota, &vz2ev_imf, &vx2ev_imf,&vy2ev_imf, &FF22, &FIMF22, &zdummy1, &adummy1,&tkedummy1, &jprf2, &inttype, &inum,EV_TEMP,&IEV_TAB_TEMP,&NbLamimf);

               for(G4int IJ = 0; IJ< IEV_TAB_TEMP;IJ++){ 
               EV_TAB[IJ+IEV_TAB][0] = EV_TEMP[IJ][0];
               EV_TAB[IJ+IEV_TAB][1] = EV_TEMP[IJ][1];
               EV_TAB[IJ+IEV_TAB][5] = EV_TEMP[IJ][5];
//
//               EV_TAB(IJ+IEV_TAB,3) = EV_TEMP(IJ,3) + VX_PREF
//               EV_TAB(IJ+IEV_TAB,4) = EV_TEMP(IJ,4) + VY_PREF
//               EV_TAB(IJ+IEV_TAB,5) = EV_TEMP(IJ,5) + VZ_PREF
// Lorentz transformation
               lorentz_boost(V_CM[0],V_CM[1],V_CM[2],
                EV_TEMP[IJ][2],EV_TEMP[IJ][3],EV_TEMP[IJ][4],
                &VXOUT,&VYOUT,&VZOUT);
               lorentz_boost(VX2_IMF,VY2_IMF,VZ2_IMF,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
               EV_TAB[IJ+IEV_TAB][2] = VX2OUT;
               EV_TAB[IJ+IEV_TAB][3] = VY2OUT;
               EV_TAB[IJ+IEV_TAB][4] = VZ2OUT;
               }
               IEV_TAB = IEV_TAB + IEV_TAB_TEMP;
// As IMF is not allowed to emit IMF, adummy1=zdummy1=0

      AFPIMF = idnint(affimf);
      ZFPIMF = idnint(zffimf);
      SFPIMF = NbLamimf;

// vi1_imf, vi2_imf - velocities of imf and partner from TKE;
// vi1ev_imf, vi2_imf - recoil of partner and imf due to evaporation
// Lorentz kinematics - DM 18/5/2010
//        VIMF_CM(1) = V_CM(1) + VX2_IMF + VX2EV_IMF
//        VIMF_CM(2) = V_CM(2) + VY2_IMF + VY2EV_IMF
//        VIMF_CM(3) = V_CM(3) + VZ2_IMF + VZ2EV_IMF
        lorentz_boost(VX2_IMF,VY2_IMF,VZ2_IMF,
                V_CM[0],V_CM[1],V_CM[2],
                &VXOUT,&VYOUT,&VZOUT);
        lorentz_boost(vx2ev_imf,vy2ev_imf,vz2ev_imf,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
        VIMF_CM[0] = VX2OUT;
        VIMF_CM[1] = VY2OUT;
        VIMF_CM[2] = VZ2OUT;
// Lorentz kinematics 
//       VFP1_CM(1) = V_CM(1) + VX1_IMF + VX1EV_IMF
//       VFP1_CM(2) = V_CM(2) + VY1_IMF + VY1EV_IMF
//       VFP1_CM(3) = V_CM(3) + VZ1_IMF + VZ1EV_IMF
        lorentz_boost(VX1_IMF,VY1_IMF,VZ1_IMF,
                V_CM[0],V_CM[1],V_CM[2],
                &VXOUT,&VYOUT,&VZOUT);
        lorentz_boost(vx1ev_imf,vy1ev_imf,vz1ev_imf,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
        VFP1_CM[0] = VX2OUT;
        VFP1_CM[1] = VY2OUT;
        VFP1_CM[2] = VZ2OUT;

      if(FF11==0 && FIMF11==0){
// heavy partner deexcites by emission of light particles
      AFP1 = idnint(aff);
      ZFP1 = idnint(zff);
      SFP1 = NbLamH;
      ZFP2 = 0;
      AFP2 = 0;
      SFP2 = 0;
      ftype = 2;
      AFPIMF = idnint(affimf);
      ZFPIMF = idnint(zffimf);
      SFPIMF = NbLamimf;
        for(G4int I=0;I<3;I++)
         VFP2_CM[I] = 0.0;


      } else if(FF11==1 && FIMF11==0){
// Heavy partner fissions
      varntp->kfis = 1;
      if(NbLam0>0)varntp->kfis = 20;
//
      opt->optimfallowed = 0; //  IMF is not allowed
      fiss->ifis = 0;         //  fission is not allowed
//
      zf = zff;
      af = aff;
      ee = EEIMFP;
    //  ftype1=0;
      ftype=21;

      G4int IEV_TAB_FIS = 0,imode=0;

      G4double vx1_fission=0.,vy1_fission=0.,vz1_fission=0.;
      G4double vx2_fission=0.,vy2_fission=0.,vz2_fission=0.;
      G4double vx_eva_sc=0.,vy_eva_sc=0.,vz_eva_sc=0.;

      fission(af,zf,ee,jprf1,
          &vx1_fission,&vy1_fission,&vz1_fission,
          &vx2_fission,&vy2_fission,&vz2_fission,
          &ZFP1,&AFP1,&SFP1,&ZFP2,&AFP2,&SFP2,&imode,
          &vx_eva_sc,&vy_eva_sc,&vz_eva_sc,EV_TEMP,&IEV_TAB_FIS,&NbLamH);

               for(int IJ = 0; IJ< IEV_TAB_FIS;IJ++){ 
               EV_TAB[IJ+IEV_TAB][0] = EV_TEMP[IJ][0];
               EV_TAB[IJ+IEV_TAB][1] = EV_TEMP[IJ][1];
               EV_TAB[IJ+IEV_TAB][5] = EV_TEMP[IJ][5];
// Lorentz kinematics
//               EV_TAB(IJ+IEV_TAB,3) = EV_TEMP(IJ,3) + VX_PREF
//               EV_TAB(IJ+IEV_TAB,4) = EV_TEMP(IJ,4) + VY_PREF
//               EV_TAB(IJ+IEV_TAB,5) = EV_TEMP(IJ,5) + VZ_PREF
// Lorentz transformation
               lorentz_boost(VFP1_CM[0],VFP1_CM[1],VFP1_CM[2],
                EV_TEMP[IJ][2],EV_TEMP[IJ][3],EV_TEMP[IJ][4],
                &VXOUT,&VYOUT,&VZOUT);
               EV_TAB[IJ+IEV_TAB][2] = VXOUT;
               EV_TAB[IJ+IEV_TAB][3] = VYOUT;
               EV_TAB[IJ+IEV_TAB][4] = VZOUT;
               }
               IEV_TAB = IEV_TAB + IEV_TAB_FIS;

    //  if(imode==1) ftype1 = 1;    // S1 mode
    //  if(imode==2) ftype1 = 2;    // S2 mode

// Lorentz kinematics
//        VFP1_CM(1) = V_CM(1) + VX1_IMF + VX1EV_IMF + VX1_FISSION +
//     &               VX_EVA_SC ! Velocity of FF1 in x
//        VFP1_CM(2) = V_CM(2) + VY1_IMF + VY1EV_IMF + VY1_FISSION +
//     &               VY_EVA_SC ! Velocity of FF1 in y
//        VFP1_CM(3) = V_CM(3) + VZ1_IMF + VZ1EV_IMF + VZ1_FISSION +
//     &               VZ_EVA_SC ! Velocity of FF1 in x
        lorentz_boost(VX1_IMF,VY1_IMF,VZ1_IMF,
                V_CM[0],V_CM[1],V_CM[2],
                &VXOUT,&VYOUT,&VZOUT);
        lorentz_boost(vx1ev_imf,vy1ev_imf,vz1ev_imf,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
        lorentz_boost(vx1_fission,vy1_fission,vz1_fission,
                VX2OUT,VY2OUT,VZ2OUT,
                &VXOUT,&VYOUT,&VZOUT);
        lorentz_boost(vx_eva_sc,vy_eva_sc,vz_eva_sc,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
        VFP1_CM[0] = VX2OUT;
        VFP1_CM[1] = VY2OUT;
        VFP1_CM[2] = VZ2OUT;

// Lorentz kinematics
//        VFP2_CM(1) = V_CM(1) + VX1_IMF + VX1EV_IMF + VX2_FISSION +
//     &               VX_EVA_SC ! Velocity of FF2 in x
//        VFP2_CM(2) = V_CM(2) + VY1_IMF + VY1EV_IMF + VY2_FISSION +
//     &               VY_EVA_SC ! Velocity of FF2 in y
//        VFP2_CM(3) = V_CM(3) + VZ1_IMF + VZ1EV_IMF + VZ2_FISSION +
//     &               VZ_EVA_SC ! Velocity of FF2 in x
        lorentz_boost(VX1_IMF,VY1_IMF,VZ1_IMF,
                V_CM[0],V_CM[1],V_CM[2],
                &VXOUT,&VYOUT,&VZOUT);
        lorentz_boost(vx1ev_imf,vy1ev_imf,vz1ev_imf,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
        lorentz_boost(vx2_fission,vy2_fission,vz2_fission,
                VX2OUT,VY2OUT,VZ2OUT,
                &VXOUT,&VYOUT,&VZOUT);
        lorentz_boost(vx_eva_sc,vy_eva_sc,vz_eva_sc,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
        VFP2_CM[0] = VX2OUT;
        VFP2_CM[1] = VY2OUT;
        VFP2_CM[2] = VZ2OUT;

      } else if(FF11==0 && FIMF11==1){
// Heavy partner emits imf, consequtive imf emission or fission is not allowed
      opt->optimfallowed = 0; //  IMF is not allowed
      fiss->ifis = 0;         //  fission is not allowed
//
      zf = zff;
      af = aff;
      ee = EEIMFP;
      aimf = adummy;
      zimf = zdummy;
      tkeimf = tkedummy;
      FF11 = 0;
      FIMF11 = 0;
      ftype = 22;
//  Lambda particles 
      G4int NbLamH1=0;
      G4int NbLamimf1=0;
      G4double pbH1 = (af-zf) / (af-zf+aimf-zimf); 
      for(G4int i=0;i<NbLamH;i++){
       if(G4AblaRandom::flat()<pbH1){
        NbLamH1++;
       }else{
        NbLamimf1++;
       }
      }
//
// Velocities of IMF and partner: 1 denotes partner, 2 denotes IMF
      EkinR1 = tkeimf * aimf / (af+aimf);
      EkinR2 = tkeimf * af / (af+aimf);
      V1 = std::sqrt(EkinR1/af) * 1.3887;
      V2 = std::sqrt(EkinR2/aimf) * 1.3887;
      G4double VZ1_IMFS = (2.0 * G4AblaRandom::flat() - 1.0) * V1;
             VPERP1 = std::sqrt(V1*V1 - VZ1_IMFS*VZ1_IMFS);
             ALPHA1 = G4AblaRandom::flat() * 2. * 3.142;
      G4double VX1_IMFS = VPERP1 * std::sin(ALPHA1);
      G4double VY1_IMFS = VPERP1 * std::cos(ALPHA1);
      G4double VX2_IMFS = - VX1_IMFS / V1 * V2;
      G4double VY2_IMFS = - VY1_IMFS / V1 * V2;
      G4double VZ2_IMFS = - VZ1_IMFS / V1 * V2;

             EEIMFP = ee * af /(af + aimf);
             EEIMF = ee * aimf /(af + aimf);

// Decay of heavy partner
      IINERTTOT = 0.40 * 931.490 * 1.160*1.160 *( std::pow(aimf,5.0/3.0) + std::pow(af,5.0/3.0)) + 931.490 * 1.160*1.160*aimf*af/(aimf+af)*(std::pow(aimf,1./3.) + std::pow(af,1./3.))*(std::pow(aimf,1./3.) + std::pow(af,1./3.));

      JPRFHEAVY = jprf1 * 0.4 * 931.49 * 1.16*1.16 * std::pow(af,5.0/3.0) / IINERTTOT;
      JPRFLIGHT = jprf1 * 0.4 * 931.49 * 1.16*1.16 * std::pow(aimf,5.0/3.0) / IINERTTOT;

     G4double zffs=0.,affs=0.,vx1ev_imfs=0.,vy1ev_imfs=0.,vz1ev_imfs=0.,jprf3=0.;

     evapora(zf,af,&EEIMFP,JPRFHEAVY, &zffs, &affs, &mtota, &vz1ev_imfs, &vx1ev_imfs,&vy1ev_imfs, &FF11, &FIMF11, &zdummy, &adummy,&tkedummy, &jprf3, &inttype, &inum,EV_TEMP,&IEV_TAB_TEMP,&NbLamH1);

               for(G4int IJ = 0; IJ< IEV_TAB_TEMP;IJ++){ 
               EV_TAB[IJ+IEV_TAB][0] = EV_TEMP[IJ][0];
               EV_TAB[IJ+IEV_TAB][1] = EV_TEMP[IJ][1];
               EV_TAB[IJ+IEV_TAB][5] = EV_TEMP[IJ][5];
//
//               EV_TAB(IJ+IEV_TAB,3) = EV_TEMP(IJ,3) + VX_PREF
//               EV_TAB(IJ+IEV_TAB,4) = EV_TEMP(IJ,4) + VY_PREF
//               EV_TAB(IJ+IEV_TAB,5) = EV_TEMP(IJ,5) + VZ_PREF
// Lorentz transformation
               lorentz_boost(VFP1_CM[0],VFP1_CM[1],VFP1_CM[2],
                EV_TEMP[IJ][2],EV_TEMP[IJ][3],EV_TEMP[IJ][4],
                &VXOUT,&VYOUT,&VZOUT);
               lorentz_boost(vx1ev_imfs,vy1ev_imfs,vz1ev_imfs,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
               EV_TAB[IJ+IEV_TAB][2] = VX2OUT;
               EV_TAB[IJ+IEV_TAB][3] = VY2OUT;
               EV_TAB[IJ+IEV_TAB][4] = VZ2OUT;
               }
               IEV_TAB = IEV_TAB + IEV_TAB_TEMP;

// For IMF - fission and IMF emission are not allowed
      opt->optimfallowed = 0; //  IMF is not allowed
      fiss->ifis = 0;         //  fission is not allowed
//
      FF22 = 0;
      FIMF22 = 0;
// Decay of "second" IMF
     G4double zffimfs=0.,affimfs=0.,vx2ev_imfs=0.,vy2ev_imfs=0.,vz2ev_imfs=0.,jprf4=0.;

     evapora(zimf,aimf,&EEIMF,JPRFLIGHT, &zffimfs, &affimfs, &mtota, &vz2ev_imfs, &vx2ev_imfs,&vy2ev_imfs, &FF22, &FIMF22, &zdummy1, &adummy1,&tkedummy1, &jprf4, &inttype, &inum,EV_TEMP,&IEV_TAB_TEMP,&NbLamimf1);

               for(G4int IJ = 0; IJ< IEV_TAB_TEMP;IJ++){ 
               EV_TAB[IJ+IEV_TAB][0] = EV_TEMP[IJ][0];
               EV_TAB[IJ+IEV_TAB][1] = EV_TEMP[IJ][1];
               EV_TAB[IJ+IEV_TAB][5] = EV_TEMP[IJ][5];
//
//               EV_TAB(IJ+IEV_TAB,3) = EV_TEMP(IJ,3) + VX_PREF
//               EV_TAB(IJ+IEV_TAB,4) = EV_TEMP(IJ,4) + VY_PREF
//               EV_TAB(IJ+IEV_TAB,5) = EV_TEMP(IJ,5) + VZ_PREF
// Lorentz transformation
               lorentz_boost(VFP1_CM[0],VFP1_CM[1],VFP1_CM[2],
                EV_TEMP[IJ][2],EV_TEMP[IJ][3],EV_TEMP[IJ][4],
                &VXOUT,&VYOUT,&VZOUT);
               lorentz_boost(vx2ev_imfs,vy2ev_imfs,vz2ev_imfs,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
               EV_TAB[IJ+IEV_TAB][2] = VX2OUT;
               EV_TAB[IJ+IEV_TAB][3] = VY2OUT;
               EV_TAB[IJ+IEV_TAB][4] = VZ2OUT;
               }
               IEV_TAB = IEV_TAB + IEV_TAB_TEMP;

      AFP1 = idnint(affs);
      ZFP1 = idnint(zffs);
      SFP1 = NbLamH1;
      ZFP2 = idnint(zffimfs);
      AFP2 = idnint(affimfs);
      SFP2 = NbLamimf1;

// Velocity of final heavy residue
// Lorentz kinematics 
//       VFP1_CM(1) = V_CM(1) + VX1_IMF + VX1EV_IMF
//       VFP1_CM(2) = V_CM(2) + VY1_IMF + VY1EV_IMF
//       VFP1_CM(3) = V_CM(3) + VZ1_IMF + VZ1EV_IMF
        lorentz_boost(VX1_IMF,VY1_IMF,VZ1_IMF,
                V_CM[0],V_CM[1],V_CM[2],
                &VXOUT,&VYOUT,&VZOUT);
               lorentz_boost(vx1ev_imf,vy1ev_imf,vz1ev_imf,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
               lorentz_boost(VX1_IMFS,VY1_IMFS,VZ1_IMFS,
                VX2OUT,VY2OUT,VZ2OUT,
                &VXOUT,&VYOUT,&VZOUT);
               lorentz_boost(vx1ev_imfs,vy1ev_imfs,vz1ev_imfs,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
        VFP1_CM[0] = VX2OUT;
        VFP1_CM[1] = VY2OUT;
        VFP1_CM[2] = VZ2OUT;

// Velocity of the second IMF
// Lorentz kinematics 
//       VFP1_CM(1) = V_CM(1) + VX1_IMF + VX1EV_IMF
//       VFP1_CM(2) = V_CM(2) + VY1_IMF + VY1EV_IMF
//       VFP1_CM(3) = V_CM(3) + VZ1_IMF + VZ1EV_IMF
        lorentz_boost(VX1_IMF,VY1_IMF,VZ1_IMF,
                V_CM[0],V_CM[1],V_CM[2],
                &VXOUT,&VYOUT,&VZOUT);
               lorentz_boost(vx1ev_imf,vy1ev_imf,vz1ev_imf,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
               lorentz_boost(VX2_IMFS,VY2_IMFS,VZ2_IMFS,
                VX2OUT,VY2OUT,VZ2OUT,
                &VXOUT,&VYOUT,&VZOUT);
               lorentz_boost(vx2ev_imfs,vy2ev_imfs,vz2ev_imfs,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
        VFP2_CM[0] = VX2OUT;
        VFP2_CM[1] = VY2OUT;
        VFP2_CM[2] = VZ2OUT;
      }//second decay
    }// if(ftype == 2)

// Only evaporation of light particles
    if(ftype!=1 && ftype!=21){

// ----------- RESOLVE UNSTABLE NUCLEI
      IOUNSTABLE=0;

      unstable_nuclei(AFP1,ZFP1,&afpnew,&zfpnew,IOUNSTABLE,
       VFP1_CM[0],VFP1_CM[1],VFP1_CM[2],
       &VP1X,&VP1Y,&VP1Z,EV_TAB_TEMP,&ILOOP);

      if(IOUNSTABLE==1){
      AFP1 = afpnew;
      ZFP1 = zfpnew;
      VFP1_CM[0] = VP1X;
      VFP1_CM[1] = VP1Y;
      VFP1_CM[2] = VP1Z;
         for(G4int I = 0;I<ILOOP;I++){
          for(G4int IJ = 0; IJ<5; IJ++)
           EV_TAB[I+IEV_TAB][IJ] = EV_TAB_TEMP[I][IJ];
         }
        IEV_TAB = IEV_TAB + ILOOP;
      }

      if(ftype>1){
       IOUNSTABLE=0;

       unstable_nuclei(AFPIMF,ZFPIMF,&afpnew,&zfpnew,IOUNSTABLE,
        VIMF_CM[0],VIMF_CM[1],VIMF_CM[2],
        &VP1X,&VP1Y,&VP1Z,EV_TAB_TEMP,&ILOOP);

        if(IOUNSTABLE==1){
        AFPIMF = afpnew;
        ZFPIMF = zfpnew;
        VIMF_CM[0] = VP1X;
        VIMF_CM[1] = VP1Y;
        VIMF_CM[2] = VP1Z;
         for(G4int I = 0;I<ILOOP;I++){
          for(G4int IJ = 0; IJ<5; IJ++)
           EV_TAB[I+IEV_TAB][IJ] = EV_TAB_TEMP[I][IJ];
         }
        IEV_TAB = IEV_TAB + ILOOP;
        }

       if(ftype>2){
        IOUNSTABLE=0;

        unstable_nuclei(AFP2,ZFP2,&afpnew,&zfpnew,IOUNSTABLE,
        VFP2_CM[0],VFP2_CM[1],VFP2_CM[2],
        &VP1X,&VP1Y,&VP1Z,EV_TAB_TEMP,&ILOOP);

        if(IOUNSTABLE==1){
        AFP2 = afpnew;
        ZFP2 = zfpnew;
        VFP2_CM[0] = VP1X;
        VFP2_CM[1] = VP1Y;
        VFP2_CM[2] = VP1Z;
         for(G4int I = 0;I<ILOOP;I++){
          for(G4int IJ = 0; IJ<5; IJ++)
           EV_TAB[I+IEV_TAB][IJ] = EV_TAB_TEMP[I][IJ];
         }
         IEV_TAB = IEV_TAB + ILOOP;
        }
       }// ftype>2
      }// ftype>1
    }


// For the case of fission:
    if(ftype==1 || ftype==21){
// ----------- RESOLVE UNSTABLE NUCLEI
      IOUNSTABLE=0;
// ----------- Fragment 1
      unstable_nuclei(AFP1,ZFP1,&afpnew,&zfpnew,IOUNSTABLE,
       VFP1_CM[0],VFP1_CM[1],VFP1_CM[2],
       &VP1X,&VP1Y,&VP1Z,EV_TAB_TEMP,&ILOOP);

      if(IOUNSTABLE==1){
      AFP1 = afpnew;
      ZFP1 = zfpnew;
      VFP1_CM[0] = VP1X;
      VFP1_CM[1] = VP1Y;
      VFP1_CM[2] = VP1Z;
         for(G4int I = 0;I<ILOOP;I++){
          for(G4int IJ = 0; IJ<5; IJ++)
           EV_TAB[I+IEV_TAB][IJ] = EV_TAB_TEMP[I][IJ];
         }
        IEV_TAB = IEV_TAB + ILOOP;
      }

      IOUNSTABLE=0;
// ----------- Fragment 2
      unstable_nuclei(AFP2,ZFP2,&afpnew,&zfpnew,IOUNSTABLE,
       VFP2_CM[0],VFP2_CM[1],VFP2_CM[2],
       &VP1X,&VP1Y,&VP1Z,EV_TAB_TEMP,&ILOOP);

      if(IOUNSTABLE==1){
      AFP2 = afpnew;
      ZFP2 = zfpnew;
      VFP2_CM[0] = VP1X;
      VFP2_CM[1] = VP1Y;
      VFP2_CM[2] = VP1Z;
         for(G4int I = 0;I<ILOOP;I++){
          for(G4int IJ = 0; IJ<5; IJ++)
           EV_TAB[I+IEV_TAB][IJ] = EV_TAB_TEMP[I][IJ];
         }
        IEV_TAB = IEV_TAB + ILOOP;
      }

       if(ftype==21){
       IOUNSTABLE=0;
// ----------- Fragment IMF
       unstable_nuclei(AFPIMF,ZFPIMF,&afpnew,&zfpnew,IOUNSTABLE,
        VIMF_CM[0],VIMF_CM[1],VIMF_CM[2],
        &VP1X,&VP1Y,&VP1Z,EV_TAB_TEMP,&ILOOP);

        if(IOUNSTABLE==1){
        AFPIMF = afpnew;
        ZFPIMF = zfpnew;
        VIMF_CM[0] = VP1X;
        VIMF_CM[1] = VP1Y;
        VIMF_CM[2] = VP1Z;
         for(G4int I = 0;I<ILOOP;I++){
          for(G4int IJ = 0; IJ<5; IJ++)
           EV_TAB[I+IEV_TAB][IJ] = EV_TAB_TEMP[I][IJ];
         }
        IEV_TAB = IEV_TAB + ILOOP;
        }
       }// ftype=21
    }

// Cross check
      if((ftype == 1 ||  ftype == 21) && (AFP2<=0 || AFP1<=0 || ZFP2<=0 || ZFP1<=0)){
       std::cout << "ZFP1:" << ZFP1 << std::endl;
       std::cout << "AFP1:" << AFP1 << std::endl;
       std::cout << "ZFP2:" << ZFP2 << std::endl;
       std::cout << "AFP2:" << AFP2 << std::endl;
      }

//     Put heavy residues in the EV_TAB array
       EV_TAB[IEV_TAB][0] = ZFP1;
       EV_TAB[IEV_TAB][1] = AFP1;
       EV_TAB[IEV_TAB][5] = SFP1;
       EV_TAB[IEV_TAB][2] = VFP1_CM[0];
       EV_TAB[IEV_TAB][3] = VFP1_CM[1];
       EV_TAB[IEV_TAB][4] = VFP1_CM[2];
       IEV_TAB = IEV_TAB + 1;

       if(AFP2>0){
       EV_TAB[IEV_TAB][0] = ZFP2;
       EV_TAB[IEV_TAB][1] = AFP2;
       EV_TAB[IEV_TAB][5] = SFP2;
       EV_TAB[IEV_TAB][2] = VFP2_CM[0];
       EV_TAB[IEV_TAB][3] = VFP2_CM[1];
       EV_TAB[IEV_TAB][4] = VFP2_CM[2];
       IEV_TAB = IEV_TAB + 1;
       }

       if(AFPIMF>0){
       EV_TAB[IEV_TAB][0] = ZFPIMF;
       EV_TAB[IEV_TAB][1] = AFPIMF;
       EV_TAB[IEV_TAB][5] = SFPIMF;
       EV_TAB[IEV_TAB][2] = VIMF_CM[0];
       EV_TAB[IEV_TAB][3] = VIMF_CM[1];
       EV_TAB[IEV_TAB][4] = VIMF_CM[2];
       IEV_TAB = IEV_TAB + 1;
       }
// Put the array of particles in the root file of INCL
   FillData(IMULTBU,IEV_TAB);
   return;
}

// Evaporation code
void G4Abla::initEvapora()
{

 //     40	C                       BFPRO,SNPRO,SPPRO,SHELL                         
  //     41	C                                                                       
  //     42	C     AP,ZP,AT,ZT   - PROJECTILE AND TARGET MASSES                      
  //     43	C     EAP,BETA      - BEAM ENERGY PER NUCLEON, V/C                      
  //     44	C     BMAXNUC       - MAX. IMPACT PARAMETER FOR NUCL. REAC.             
  //     45	C     CRTOT,CRNUC   - TOTAL AND NUCLEAR REACTION CROSS SECTION          
  //     46	C     R_0,R_P,R_T,  - RADIUS PARAMETER, PROJECTILE+ TARGET RADII        
  //     47	C     IMAX,IRNDM,PI - MAXIMUM NUMBER OF EVENTS, DUMMY, 3.141...         
  //     48	C     BFPRO         - FISSION BARRIER OF THE PROJECTILE                 
  //     49	C     SNPRO         - NEUTRON SEPARATION ENERGY OF THE PROJECTILE       
  //     50	C     SPPRO         - PROTON    "           "   "    "   "              
  //     51	C     SHELL         - GROUND STATE SHELL CORRECTION                     
  //     52	C---------------------------------------------------------------------  
  //     53	C                                                                       
  //     54	C     ENERGIES WIDTHS AND CROSS SECTIONS FOR EM EXCITATION              
  //     55	C     COMMON /EMDPAR/ EGDR,EGQR,FWHMGDR,FWHMGQR,CREMDE1,CREMDE2,        
  //     56	C                     AE1,BE1,CE1,AE2,BE2,CE2,SR1,SR2,XR                
  //     57	C                                                                       
  //     58	C     EGDR,EGQR       - MEAN ENERGY OF GDR AND GQR                      
  //     59	C     FWHMGDR,FWHMGQR - FWHM OF GDR, GQR                                
  //     60	C     CREMDE1,CREMDE2 - EM CROSS SECTION FOR E1 AND E2                  
  //     61	C     AE1,BE1,CE1     - ARRAYS TO CALCULATE                             
  //     62	C     AE2,BE2,CE2     - THE EXCITATION ENERGY AFTER E.M. EXC.           
  //     63	C     SR1,SR2,XR      - WITH MONTE CARLO                                
  //     64	C---------------------------------------------------------------------  
  //     65	C                                                                       
  //     66	C     DEFORMATIONS AND G.S. SHELL EFFECTS                               
  //     67	C     COMMON /ECLD/   ECGNZ,ECFNZ,VGSLD,ALPHA                           
  //     68	C                                                                       
  //     69	C     ECGNZ - GROUND STATE SHELL CORR. FRLDM FOR A SPHERICAL G.S.       
  //     70	C     ECFNZ - SHELL CORRECTION FOR THE SADDLE POINT (NOW: == 0)         
  //     71	C     VGSLD - DIFFERENCE BETWEEN DEFORMED G.S. AND LDM VALUE            
  //     72	C     ALPHA - ALPHA GROUND STATE DEFORMATION (THIS IS NOT BETA2!)       
  //     73	C             BETA2 = SQRT(5/(4PI)) * ALPHA                             
  //     74	C---------------------------------------------------------------------  
  //     75	C                                                                       
  //     76	C     ARRAYS FOR EXCITATION ENERGY BY STATISTICAL HOLE ENERY MODEL      
  //     77	C     COMMON /EENUC/  SHE, XHE                                          
  //     78	C                                                                       
  //     79	C     SHE, XHE - ARRAYS TO CALCULATE THE EXC. ENERGY AFTER              
  //     80	C                ABRASION BY THE STATISTICAL HOLE ENERGY MODEL          
  //     81	C---------------------------------------------------------------------  
  //     82	C                                                                       
  //     83	C     G.S. SHELL EFFECT                                                 
  //     84	C     COMMON /EC2SUB/ ECNZ                                              
  //     85	C                                                                       
  //     86	C     ECNZ G.S. SHELL EFFECT FOR THE MASSES (IDENTICAL TO ECGNZ)        
  //     87	C---------------------------------------------------------------------  
  //       

  G4double MN = 939.5653301;   
  G4double MP = 938.7829835;                             

#ifdef ABLAXX_IN_GEANT4_MODE
  G4AblaDataFile *dataInterface = new G4AblaDataFile();
#else
  G4AblaDataFile *dataInterface = new G4AblaDataFile(theConfig);
#endif
  if(dataInterface->readData() == true) {
    if(verboseLevel > 0) {
      // G4cout <<"G4Abla: Datafiles read successfully." << G4endl;
    }
  }
  else {
    //    G4Exception("ERROR: Failed to read datafiles.");
  }
  
  for(G4int z = 0; z < 99; z++) { //do 30  z = 0,98,1                                                 
    for(G4int n = 0; n < 154; n++) { //do 31  n = 0,153,1                                              
      ecld->ecfnz[n][z] = 0.e0;
      ec2sub->ecnz[n][z] = dataInterface->getEcnz(n,z);
      ecld->ecgnz[n][z] = dataInterface->getEcnz(n,z);
      ecld->alpha[n][z] = dataInterface->getAlpha(n,z);
      ecld->vgsld[n][z] = dataInterface->getVgsld(n,z);
      ecld->rms[n][z] = dataInterface->getRms(n,z);
    }
  }

  for(G4int z = 0; z < 137; z++){                                                  
    for(G4int n = 0; n < 251; n++){  
      ecld->beta2[n][z] = dataInterface->getBeta2(n,z);
      ecld->beta4[n][z] = dataInterface->getBeta4(n,z);
    }
  }

  for(G4int z = 0; z < 500; z++) {
    for(G4int a = 0; a < 500; a++) {
      pace->dm[z][a] = dataInterface->getPace2(z,a);
    }
  }



  G4double mfrldm[154][13];
// For 2 < Z < 12 we take "experimental" shell corrections instead of calculated
// Read FRLDM tables
  for(G4int i=1;i<13;i++){
   for(G4int j=1;j<154;j++){
      if(dataInterface->getMexpID(j,i)==1){
       masses->mexpiop[j][i]=1;  
      }else{
       masses->mexpiop[j][i]=0;
      } 
// LD masses (even-odd effect is later considered according to Ignatyuk)
      if(i==0 && j==0)
       mfrldm[j][i] = 0.;
      else
       mfrldm[j][i] = MP*i+MN*j+eflmac(i+j,i,1,0);
   }
  }

  G4double e0=0.;
  for(G4int i=1;i<13;i++){
   for(G4int j=1;j<154;j++){
      masses->bind[j][i]=0.;
      if(masses->mexpiop[j][i]==1){
        if(j<3){

          ec2sub->ecnz[j][i] = 0.0;
          ecld->ecgnz[j][i] = ec2sub->ecnz[j][i];
          masses->bind[j][i] = dataInterface->getMexp(j,i)-MP*i -MN*j;
          ecld->vgsld[j][i]=0.;

          e0=0.;
        }else{
// For these nuclei, we take "experimental" ground-state shell corrections
//
// Parametrization of CT model by Ignatyuk; note that E0 is shifted to correspond
// to pairing shift in Fermi-gas model (there, energy is shifted taking odd-odd nuclei as bassis)
              G4double para=0.;
              parite(j+i,&para);
                if(para<0.0){
// e-o, o-e
                 e0 =  0.285+11.17*std::pow(j+i,-0.464) -0.390-0.00058*(j+i);
                }else{
                  G4double parz=0.;
                  parite(i,&parz);
                  if (parz>0.0){
// e-e
                   e0 = 22.34*std::pow(j+i,-0.464)-0.235;
                  }else{
// o-o
                   e0 = 0.0;
                  }
                }
//
                if((j==i)&&mod(j,2)==1&&mod(i,2)==1){
                 e0 = e0 - 30.0*(1.0/G4double(j+i));
                }

            G4double delta_tot = ec2sub->ecnz[j][i] - ecld->vgsld[j][i];
            ec2sub->ecnz[j][i] = dataInterface->getMexp(j,i) - (mfrldm[j][i] - e0);

            ecld->vgsld[j][i] = max(0.0,ec2sub->ecnz[j][i] - delta_tot);
            ecld->ecgnz[j][i] = ec2sub->ecnz[j][i];

        }//if j
     }//if mexpiop
   }
  }
//
  delete dataInterface;
}

void G4Abla::SetParametersG4(G4int z, G4int a)
{
  //A and Z for the target
  fiss->at = a;
  fiss->zt = z;

  // shell+pairing.0-1-2-3 for IMFs
  opt->optshpimf = 0;

  //collective enhancement switched on 1 or off 0 in densn (qr=val or =1.)
  fiss->optcol = 1;
  if(fiss->zt<83 && fiss->zt>56){
  fiss->optshp = 1;
  }
  if(fiss->zt<=56){  
  fiss->optcol = 0;
  fiss->optshp = 3;
  }
}

void G4Abla::SetParameters()
{
/*
C     IFIS =   INTEGER SWITCH FOR FISSION
C     OPTSHP = INTEGER SWITCH FOR SHELL CORRECTION IN MASSES/ENERGY
C            =0 NO MICROSCOPIC CORRECTIONS IN MASSES AND ENERGY
C            =1 SHELL , NO PAIRING CORRECTION
C            =2 PAIRING, NO SHELL CORRECTION
C            =3 SHELL AND PAIRING CORRECTION IN MASSES AND ENERGY
C     OPTCOL =0,1 COLLECTIVE ENHANCEMENT SWITCHED ON 1 OR OFF 0 IN DENSN
C     OPTAFAN=0,1 SWITCH FOR AF/AN = 1 IN DENSNIV 0 AF/AN>1 1 AF/AN=1
C     BET  =  REAL    REDUCED FRICTION COEFFICIENT / 10**(+21) S**(-1)
C     OPTXFIS= INTEGER 0,1,2 FOR MYERS & SWIATECKI, DAHLINGER, ANDREYEV
C              FISSILITY PARAMETER.
C
C     NUCLEAR LEVEL DENSITIES:
C     AV     = REAL KOEFFICIENTS FOR CALCULATION OF A(TILDE)
C     AS     = REAL LEVEL DENSITY PARAMETER
C     AK     = REAL
*/

  // switch-fission.1=on.0=off
  fiss->ifis = 1;

  // shell+pairing.0-1-2-3
  fiss->optshp = 3; 
  if(fiss->zt<84 && fiss->zt>56)
  fiss->optshp = 1;

  // optemd =0,1  0 no emd, 1 incl. emd                                
  opt->optemd = 1;
  // read(10,*,iostat=io) dum(10),optcha                               
  opt->optcha = 1;

  // shell+pairing.0-1-2-3 for IMFs
  opt->optshpimf = 0;
  opt->optimfallowed = 1;

  // nuclear.viscosity.(beta)
  fiss->bet = 4.5;

  //collective enhancement switched on 1 or off 0 in densn (qr=val or =1.)
  fiss->optcol = 1;
  if(fiss->zt<=56){  
  fiss->optcol = 0;
  fiss->optshp = 3;
  }
  //collective enhancement parameters
  fiss->ucr = 40.;
  fiss->dcr = 10.;

  // switch for temperature constant model (CTM)
  fiss->optct = 1;

  ald->optafan = 0;

  ald->av = 0.0730;
  ald->as = 0.0950;
  ald->ak = 0.0000;

  fiss->optxfis = 3;

// Multi-fragmentation
  T_freeze_out_in = -6.5;

}

void G4Abla::mglw(G4double a, G4double z, G4double *el)
{
  // MODEL DE LA GOUTTE LIQUIDE DE C. F. WEIZSACKER.
  // USUALLY AN OBSOLETE OPTION

  G4double xv = 0.0, xs = 0.0, xc = 0.0, xa = 0.0;                                   

  if ((a <= 0.01) || (z < 0.01)) {
    (*el) = 1.0e38;
  }
  else {
    xv = -15.56*a;
    xs = 17.23*std::pow(a,(2.0/3.0));

    if (a > 1.0) {
      xc = 0.7*z*(z-1.0)*std::pow((a-1.0),(-1.e0/3.e0));
    }
    else {
      xc = 0.0;
    }
  }

  xa = 23.6*(std::pow((a-2.0*z),2)/a);
  (*el) = xv+xs+xc+xa;
  return;	
}

void G4Abla::mglms(G4double a, G4double z, G4int refopt4, G4double *el)
{
  // USING FUNCTION EFLMAC(IA,IZ,0)                                    
  // 
  // REFOPT4 = 0 : WITHOUT MICROSCOPIC CORRECTIONS                     
  // REFOPT4 = 1 : WITH SHELL CORRECTION                               
  // REFOPT4 = 2 : WITH PAIRING CORRECTION                             
  // REFOPT4 = 3 : WITH SHELL- AND PAIRING CORRECTION                  

  //   1839	C-----------------------------------------------------------------------
  //   1840	C     A1       LOCAL    MASS NUMBER (INTEGER VARIABLE OF A)             
  //   1841	C     Z1       LOCAL    NUCLEAR CHARGE (INTEGER VARIABLE OF Z)          
  //   1842	C     REFOPT4           OPTION, SPECIFYING THE MASS FORMULA (SEE ABOVE) 
  //   1843	C     A                 MASS NUMBER                                     
  //   1844	C     Z                 NUCLEAR CHARGE                                  
  //   1845	C     DEL               PAIRING CORRECTION                              
  //   1846	C     EL                BINDING ENERGY                                  
  //   1847	C     ECNZ( , )         TABLE OF SHELL CORRECTIONS                      
  //   1848	C-----------------------------------------------------------------------
  //   1849	C                                                                       
  G4int a1 = idnint(a);
  G4int z1 = idnint(z);
  G4int n1 = a1-z1;

  if ( (a1 <= 0) || (z1 <= 0) || ((a1-z1) <= 0) )  { //then 
    // modif pour recuperer une masse p et n correcte:
    (*el) = 1.e38;
    return;
    //    goto mglms50;
  }
  else {
    // binding energy incl. pairing contr. is calculated from                
    // function eflmac                                                       
    (*el) = eflmac(a1,z1,0,refopt4);

    if (refopt4 > 0) {
      if (refopt4 != 2) {
	(*el) = (*el) + ec2sub->ecnz[a1-z1][z1];
      }
    }
    
    if(z1>=90){
      if(n1<=145){
         (*el) = (*el) + (12.552-0.1436*z1);
      }else{
        if(n1>145&&n1<=152){
         (*el) = (*el) + ((152.4-1.77*z1)+(-0.972+0.0113*z1)*n1);
        }
      } 
    }

  }
  return;
}

G4double G4Abla::spdef(G4int a, G4int z, G4int optxfis)
{

  // INPUT:  A,Z,OPTXFIS MASS AND CHARGE OF A NUCLEUS,                     
  // OPTION FOR FISSILITY                                          
  // OUTPUT: SPDEF                                                         

  // ALPHA2 SADDLE POINT DEF. COHEN&SWIATECKI ANN.PHYS. 22 (1963) 406      
  // RANGING FROM FISSILITY X=0.30 TO X=1.00 IN STEPS OF 0.02              

  G4int index = 0;
  G4double x = 0.0, v = 0.0, dx = 0.0;

  const G4int alpha2Size = 37;
  // The value 0.0 at alpha2[0] added by PK.
  G4double alpha2[alpha2Size] = {0.0, 2.5464e0, 2.4944e0, 2.4410e0, 2.3915e0, 2.3482e0,
				 2.3014e0, 2.2479e0, 2.1982e0, 2.1432e0, 2.0807e0, 2.0142e0, 1.9419e0,
				 1.8714e0, 1.8010e0, 1.7272e0, 1.6473e0, 1.5601e0, 1.4526e0, 1.3164e0,
				 1.1391e0, 0.9662e0, 0.8295e0, 0.7231e0, 0.6360e0, 0.5615e0, 0.4953e0,
				 0.4354e0, 0.3799e0, 0.3274e0, 0.2779e0, 0.2298e0, 0.1827e0, 0.1373e0,
				 0.0901e0, 0.0430e0, 0.0000e0};

  dx = 0.02;
  x  = fissility(a,z,0,0.,0.,optxfis);

  v  = (x - 0.3)/dx + 1.0;
  index = idnint(v);

  if (index < 1) {
    return(alpha2[1]);
  }

  if (index == 36) {                                             
    return(alpha2[36]); 
  }
  else {
    return(alpha2[index] + (alpha2[index+1] - alpha2[index]) / dx * ( x - (0.3e0 + dx*(index-1))));
  }                                                       

  return alpha2[0]; // The algorithm is not supposed to reach this point.
}

G4double G4Abla::fissility(G4int a, G4int z, G4int ny, G4double sn, G4double slam, G4int optxfis)
{
  // CALCULATION OF FISSILITY PARAMETER                                 
  // 
  // INPUT: A,Z INTEGER MASS & CHARGE OF NUCLEUS                        
  // OPTXFIS = 0 : MYERS, SWIATECKI                              
  //           1 : DAHLINGER                                     
  //           2 : ANDREYEV                                      

  G4double aa = 0.0, zz = 0.0, i = 0.0,z2a,C_S,R,W,G,G1,G2,A_CC;
  G4double fissilityResult = 0.0;

  aa = G4double(a);
  zz = G4double(z);
  i  = G4double(a-2*z) / aa;
  z2a= zz*zz/aa-ny*(1115.-939.+sn-slam)/(0.7053*std::pow(a,2./3.));

  // myers & swiatecki droplet modell                        
  if (optxfis == 0) { //then                                            
    fissilityResult = std::pow(zz,2) / aa /50.8830e0 / (1.0e0 - 1.7826e0 * std::pow(i,2));
  }

  if (optxfis == 1) {
    // dahlinger fit:                                          
    fissilityResult = std::pow(zz,2) / aa * std::pow((49.22e0*(1.e0 - 0.3803e0*std::pow(i,2) - 20.489e0*std::pow(i,4))),(-1));
  }

  if (optxfis == 2) {
    // dubna fit:                                              
    fissilityResult = std::pow(zz,2) / aa  /(48.e0*(1.e0 - 17.22e0*std::pow(i,4)));
  }

  if (optxfis == 3) {
//  Fissiilty is calculated according to FRLDM, see Sierk, PRC 1984.
         C_S = 21.13 * (1.0 - 2.3*i*i);
         R = 1.16 * std::pow(aa,1.0/3.0);
         W = 0.704/R;
         G1 = 1.0 - 15.0/8.0*W+21.0/8.0*W*W*W;
         G2 = 1.0 + 9.0/2.0*W + 7.0*W*W + 7.0/2.0*W*W*W;
         G = 1.0 - 5.0*W*W*(G1 - 3.0/4.0*G2*std::exp(-2.0/W));
         A_CC = 3.0/5.0 * 1.44 * G / 1.16;
         fissilityResult = z2a * A_CC/(2.0*C_S);
  }

  if (fissilityResult > 1.0) {
    fissilityResult = 1.0;
  }

  if (fissilityResult < 0.0) {
    fissilityResult = 0.0;
  }

  return fissilityResult;
}

void G4Abla::evapora(G4double zprf, G4double aprf, G4double *ee_par, G4double jprf_par,G4double *zf_par, G4double *af_par, G4double *mtota_par,G4double *vleva_par, G4double *vxeva_par, G4double *vyeva_par,
G4int *ff_par,G4int *fimf_par, G4double *fzimf, G4double *faimf,G4double *tkeimf_par,G4double *jprfout, G4int *inttype_par, G4int *inum_par,G4double EV_TEMP[200][6],G4int *iev_tab_temp_par, G4int *NbLam0_par)
{
  G4double zf = zprf;
  G4double af = aprf;
  G4double ee = (*ee_par);
  G4double jprf = dint(jprf_par);
  G4double mtota = (*mtota_par);
  G4double vleva = 0.;
  G4double vxeva = 0.;
  G4double vyeva = 0.;
  G4int ff = (*ff_par);
  G4int fimf = (*fimf_par);
  G4double tkeimf = (*tkeimf_par);
  G4int inttype = (*inttype_par);
  G4int inum = (*inum_par);
  G4int NbLam0 = (*NbLam0_par);

  //    533	C                                                                       
  //    534	C     INPUT:                                                            
  //    535	C                                                                       
  //    536	C     ZPRF, APRF, EE(EE IS MODIFIED!), JPRF                             
  //    537	C                                                                       
  //    538	C     PROJECTILE AND TARGET PARAMETERS + CROSS SECTIONS                 
  //    539	C     COMMON /ABRAMAIN/ AP,ZP,AT,ZT,EAP,BETA,BMAXNUC,CRTOT,CRNUC,       
  //    540	C                       R_0,R_P,R_T, IMAX,IRNDM,PI,                     
  //    541	C                       BFPRO,SNPRO,SPPRO,SHELL                         
  //    542	C                                                                       
  //    543	C     AP,ZP,AT,ZT   - PROJECTILE AND TARGET MASSES                      
  //    544	C     EAP,BETA      - BEAM ENERGY PER NUCLEON, V/C                      
  //    545	C     BMAXNUC       - MAX. IMPACT PARAMETER FOR NUCL. REAC.             
  //    546	C     CRTOT,CRNUC   - TOTAL AND NUCLEAR REACTION CROSS SECTION          
  //    547	C     R_0,R_P,R_T,  - RADIUS PARAMETER, PROJECTILE+ TARGET RADII        
  //    548	C     IMAX,IRNDM,PI - MAXIMUM NUMBER OF EVENTS, DUMMY, 3.141...         
  //    549	C     BFPRO         - FISSION BARRIER OF THE PROJECTILE                 
  //    550	C     SNPRO         - NEUTRON SEPARATION ENERGY OF THE PROJECTILE       
  //    551	C     SPPRO         - PROTON    "           "   "    "   "              
  //    552	C     SHELL         - GROUND STATE SHELL CORRECTION                     
  //    553	C                                                                       
  //    554	C---------------------------------------------------------------------  
  //    555	C     FISSION BARRIERS                                                  
  //    556	C     COMMON /FB/     EFA                                               
  //    557	C     EFA    - ARRAY OF FISSION BARRIERS                                
  //    558	C---------------------------------------------------------------------  
  //    559	C     OUTPUT:                                                           
  //    560	C              ZF, AF, MTOTA, PLEVA, PTEVA, FF, INTTYPE, INUM           
  //    561	C                                                                       
  //    562	C     ZF,AF - CHARGE AND MASS OF FINAL FRAGMENT AFTER EVAPORATION       
  //    563	C     MTOTA _ NUMBER OF EVAPORATED ALPHAS                               
  //    564	C     PLEVA,PXEVA,PYEVA - MOMENTUM RECOIL BY EVAPORATION               
  //    565	C     INTTYPE - TYPE OF REACTION 0/1 NUCLEAR OR ELECTROMAGNETIC         
  //    566	C     FF      - 0/1 NO FISSION / FISSION EVENT                          
  //    567	C     INUM    - EVENTNUMBER                                             
  //    568	C   ____________________________________________________________________
  //    569	C  /                                                                    
  //    570	C  /  CALCUL DE LA MASSE ET CHARGE FINALES D'UNE CHAINE D'EVAPORATION   
  //    571	C  /                                                                    
  //    572	C  /  PROCEDURE FOR CALCULATING THE FINAL MASS AND CHARGE VALUES OF A   
  //    573	C  /  SPECIFIC EVAPORATION CHAIN, STARTING POINT DEFINED BY (APRF, ZPRF,
  //    574	C  /  EE)                                                               
  //    575	C  /  On ajoute les 3 composantes de l'impulsion (PXEVA,PYEVA,PLEVA)
  //    576	C  /    (actuellement PTEVA n'est pas correct; mauvaise norme...)                                               
  //    577	C  /____________________________________________________________________
  //    578	C                                                                       
  //    612	C                                                                       
  //    613	C-----------------------------------------------------------------------
  //    614	C     IRNDM             DUMMY ARGUMENT FOR RANDOM-NUMBER FUNCTION       
  //    615	C     SORTIE   LOCAL    HELP VARIABLE TO END THE EVAPORATION CHAIN      
  //    616	C     ZF                NUCLEAR CHARGE OF THE FRAGMENT                  
  //    617	C     ZPRF              NUCLEAR CHARGE OF THE PREFRAGMENT               
  //    618	C     AF                MASS NUMBER OF THE FRAGMENT                     
  //    619	C     APRF              MASS NUMBER OF THE PREFRAGMENT                  
  //    620	C     EPSILN            ENERGY BURNED IN EACH EVAPORATION STEP          
  //    621	C     MALPHA   LOCAL    MASS CONTRIBUTION TO MTOTA IN EACH EVAPORATION  
  //    622	C                        STEP                                           
  //    623	C     EE                EXCITATION ENERGY (VARIABLE)                    
  //    624	C     PROBP             PROTON EMISSION PROBABILITY                     
  //    625	C     PROBN             NEUTRON EMISSION PROBABILITY                    
  //    626	C     PROBA             ALPHA-PARTICLE EMISSION PROBABILITY             
  //    627	C     PTOTL             TOTAL EMISSION PROBABILITY                      
  //    628	C     E                 LOWEST PARTICLE-THRESHOLD ENERGY                
  //    629	C     SN                NEUTRON SEPARATION ENERGY                       
  //    630	C     SBP               PROTON SEPARATION ENERGY PLUS EFFECTIVE COULOMB 
  //    631	C                        BARRIER                                        
  //    632	C     SBA               ALPHA-PARTICLE SEPARATION ENERGY PLUS EFFECTIVE 
  //    633	C                        COULOMB BARRIER                                
  //    634	C     BP                EFFECTIVE PROTON COULOMB BARRIER                
  //    635	C     BA                EFFECTIVE ALPHA COULOMB BARRIER                 
  //    636	C     MTOTA             TOTAL MASS OF THE EVAPORATED ALPHA PARTICLES    
  //    637	C     X                 UNIFORM RANDOM NUMBER FOR NUCLEAR CHARGE        
  //    638	C     AMOINS   LOCAL    MASS NUMBER OF EVAPORATED PARTICLE              
  //    639	C     ZMOINS   LOCAL    NUCLEAR CHARGE OF EVAPORATED PARTICLE           
  //    640	C     ECP               KINETIC ENERGY OF PROTON WITHOUT COULOMB        
  //    641	C                        REPULSION                                      
  //    642	C     ECN               KINETIC ENERGY OF NEUTRON                       
  //    643	C     ECA               KINETIC ENERGY OF ALPHA PARTICLE WITHOUT COULOMB
  //    644	C                        REPULSION                                      
  //    645	C     PLEVA             TRANSVERSAL RECOIL MOMENTUM OF EVAPORATION      
  //    646	C     PTEVA             LONGITUDINAL RECOIL MOMENTUM OF EVAPORATION     
  //    647	C     FF                FISSION FLAG                                    
  //    648	C     INTTYPE           INTERACTION TYPE FLAG                           
  //    649	C     RNDX              RECOIL MOMENTUM IN X-DIRECTION IN A SINGLE STEP 
  //    650	C     RNDY              RECOIL MOMENTUM IN Y-DIRECTION IN A SINGLE STEP 
  //    651	C     RNDZ              RECOIL MOMENTUM IN Z-DIRECTION IN A SINGLE STEP 
  //    652	C     RNDN              NORMALIZATION OF RECOIL MOMENTUM FOR EACH STEP  
  //    653	C-----------------------------------------------------------------------
  //    654	C                                                                       
  //                        
  G4double epsiln = 0.0, probp = 0.0, probd = 0.0, probt = 0.0, probn = 0.0, probhe = 0.0, proba = 0.0, probg = 0.0, probimf=0.0, problamb0 = 0.0, ptotl = 0.0, e = 0.0, tcn = 0.0;  
  G4double sn = 0.0, sbp = 0.0, sbd = 0.0, sbt = 0.0, sbhe = 0.0, sba = 0.0, x = 0.0, amoins = 0.0, zmoins = 0.0,sp = 0.0, sd = 0.0, st = 0.0, she = 0.0, sa = 0.0, slamb0 = 0.0;
  G4double ecn = 0.0, ecp = 0.0, ecd = 0.0, ect = 0.0,eche = 0.0,eca = 0.0, ecg = 0.0, eclamb0 = 0.0, bp = 0.0, bd = 0.0, bt = 0.0, bhe = 0.0, ba = 0.0;
  G4double zimf= 0.0,aimf= 0.0,bimf= 0.0,sbimf= 0.0,timf= 0.0;
  G4int itest = 0, sortie=0;
  G4double probf = 0.0;
  G4double ctet1 = 0.0, stet1 = 0.0, phi1 = 0.0;
  G4double rnd = 0.0;
  G4double ef = 0.0;
  G4double ts1 = 0.0;
  G4int fgamma = 0, gammadecay = 0, flamb0decay=0;
  G4double pc = 0.0, malpha = 0.0;
  G4double jprfn=0.0, jprfp=0.0, jprfd=0.0, jprft=0.0, jprfhe=0.0, jprfa=0.0, jprflamb0 = 0.0;
  G4double tsum = 0.0;
  G4int twon;

  const G4double c = 29.9792458;
  const G4double mu = 931.494;
  const G4double mu2 = 931.494*931.494;

  G4double pleva = 0.0;
  G4double pxeva = 0.0;
  G4double pyeva = 0.0;
  G4int IEV_TAB_TEMP=0;

  for(G4int I1=0;I1<200;I1++)
  for(G4int I2=0;I2<6;I2++)
  EV_TEMP[I1][I2] = 0.0;
//
  ff = 0;
  itest = 0;
//
  evapora10:
  //
  // calculation of the probabilities for the different decay channels     
  // plus separation energies and kinetic energies of the particles  
  //
  if(ee<0.|| zf<3.)goto evapora100;
  direct(zf,af,ee,jprf,&probp,&probd,&probt,&probn,&probhe,&proba,&probg,&probimf,&probf,&problamb0,&ptotl,
	 &sn,&sbp,&sbd,&sbt,&sbhe,&sba,&slamb0,
         &ecn,&ecp,&ecd,&ect,&eche,&eca,&ecg,&eclamb0,
         &bp,&bd,&bt,&bhe,&ba,&sp,&sd,&st,&she,&sa,&ef,&ts1,inttype,inum,itest,&sortie,&tcn,
         &jprfn, &jprfp, &jprfd, &jprft, &jprfhe, &jprfa, &jprflamb0, &tsum, NbLam0);
//
// HERE THE FINAL STEPS OF THE EVAPORATION ARE CALCULATED
//
  if(ptotl==0.0) goto evapora100;

   e = dmin1(sba,sbhe,dmin1(sbt,sbhe,dmin1(sn,sbp,sbd)));

  if(e>1e30)std::cout << "ERROR AT THE EXIT OF EVAPORA,E>1.D30,AF="<< af << " ZF=" << zf << std::endl;

  if(sortie==1){
   if (probn!=0.0) {
    amoins = 1.0;
    zmoins = 0.0;
    epsiln = sn + ecn;
    pc = std::sqrt(std::pow((1.0 + (ecn)/9.3956e2),2.) - 1.0) * 9.3956e2;
    malpha = 0.0;
    fgamma = 0;
    fimf = 0;
    flamb0decay=0;
    gammadecay = 0;
   }
   else if(probp!=0.0){
    amoins = 1.0;
    zmoins = 1.0;
    epsiln = sp + ecp;
    pc = std::sqrt(std::pow((1.0 + ecp/9.3827e2),2.) - 1.0) * 9.3827e2;
    malpha = 0.0;
    fgamma = 0;
    fimf = 0;
    flamb0decay=0;
    gammadecay = 0;
   }
   else if(probd!=0.0){
    amoins = 2.0;
    zmoins = 1.0;
    epsiln = sd + ecd;
    pc = std::sqrt(std::pow((1.0 + ecd/1.875358e3),2) - 1.0) * 1.875358e3;
    malpha = 0.0;
    fgamma = 0;
    fimf = 0;
    flamb0decay=0;
    gammadecay = 0;
   }
   else if(probt!=0.0){
    amoins = 3.0;
    zmoins = 1.0;
    epsiln = st + ect;
    pc = std::sqrt(std::pow((1.0 + ect/2.80828e3),2) - 1.0) * 2.80828e3;
    malpha = 0.0;
    fgamma = 0;
    fimf = 0;
    flamb0decay=0;
    gammadecay = 0;
   }
   else if(probhe!=0.0){
    amoins = 3.0;
    zmoins = 2.0;
    epsiln = she + eche;
    pc = std::sqrt(std::pow((1.0 + eche/2.80826e3),2) - 1.0) * 2.80826e3;
    malpha = 0.0;
    fgamma = 0;
    fimf = 0;
    flamb0decay=0;
    gammadecay = 0;
   }
   else{ if(proba!=0.0){
    amoins = 4.0;
    zmoins = 2.0;
    epsiln = sa + eca;
    pc = std::sqrt(std::pow((1.0 + eca/3.72834e3),2) - 1.0) * 3.72834e3;
    malpha = 4.0;
    fgamma = 0;
    fimf = 0;
    flamb0decay=0;
    gammadecay = 0;
    }
   }
  goto direct99;
  }

  // here the normal evaporation cascade starts                            

  // random number for the evaporation
  x = G4AblaRandom::flat() * ptotl;

  itest = 0;
  if (x < proba) {
    // alpha evaporation                                                     
    amoins = 4.0;
    zmoins = 2.0;
    epsiln = sa + eca;
    pc = std::sqrt(std::pow((1.0 + eca/3.72834e3),2) - 1.0) * 3.72834e3;
    malpha = 4.0;
    fgamma = 0;
    fimf = 0;
    ff = 0;
    flamb0decay=0;
    gammadecay = 0;
    jprf=jprfa;
  }
  else if (x < proba+probhe) {
    // He3 evaporation                                                    
    amoins = 3.0;
    zmoins = 2.0;
    epsiln = she + eche;
    pc = std::sqrt(std::pow((1.0 + eche/2.80826e3),2) - 1.0) * 2.80826e3;
    malpha = 0.0;
    fgamma = 0;
    fimf = 0;
    ff = 0;
    flamb0decay=0;
    gammadecay = 0;
    jprf=jprfhe;
  }
  else if (x < proba+probhe+probt) {
    // triton evaporation                                                    
    amoins = 3.0;
    zmoins = 1.0;
    epsiln = st + ect;
    pc = std::sqrt(std::pow((1.0 + ect/2.80828e3),2) - 1.0) * 2.80828e3;
    malpha = 0.0;
    fgamma = 0;
    fimf = 0;
    ff = 0;
    flamb0decay=0;
    gammadecay = 0;
    jprf=jprft;
  }
  else if (x < proba+probhe+probt+probd) {
    // deuteron evaporation                                                    
    amoins = 2.0;
    zmoins = 1.0;
    epsiln = sd + ecd;
    pc = std::sqrt(std::pow((1.0 + ecd/1.875358e3),2) - 1.0) * 1.875358e3;
    malpha = 0.0;
    fgamma = 0;
    fimf = 0;
    ff = 0;
    flamb0decay=0;
    gammadecay = 0;
    jprf=jprfd;
  }
  else if (x < proba+probhe+probt+probd+probp) {
    // proton evaporation                                                    
    amoins = 1.0;
    zmoins = 1.0;
    epsiln = sp + ecp;
    pc = std::sqrt(std::pow((1.0 + ecp/9.3827e2),2) - 1.0) * 9.3827e2;
    malpha = 0.0;
    fgamma = 0;
    fimf = 0;
    ff = 0;
    flamb0decay=0;
    gammadecay = 0;
    jprf=jprfp;
  }
  else if (x < proba+probhe+probt+probd+probp+probn) {
    // neutron evaporation                                                   
    amoins = 1.0;
    zmoins = 0.0;
    epsiln = sn + ecn;
    pc = std::sqrt(std::pow((1.0 + (ecn)/9.3956e2),2.) - 1.0) * 9.3956e2;
    malpha = 0.0;
    fgamma = 0;
    fimf = 0;
    ff = 0;
    flamb0decay=0;
    gammadecay = 0;
    jprf=jprfn;
  }
  else if (x < proba+probhe+probt+probd+probp+probn+problamb0) {
    // lambda0 evaporation  
    amoins = 1.0;
    zmoins = 0.0;
    epsiln = slamb0 + eclamb0;
    pc = std::sqrt(std::pow((1.0 + (eclamb0)/11.1568e2),2.) - 1.0) * 11.1568e2;
    malpha = 0.0;
    fgamma = 0;
    fimf = 0;
    ff = 0;
    flamb0decay = 1;
    opt->nblan0 = opt->nblan0 -1;
    NbLam0 = NbLam0 -1;
    gammadecay = 0;
    jprf=jprflamb0;
  }
  else if (x < proba+probhe+probt+probd+probp+probn+problamb0+probg) {
    // gamma evaporation                                                    
    amoins = 0.0;
    zmoins = 0.0;
    epsiln = ecg;
    pc = ecg;
    malpha = 0.0;
    flamb0decay = 0;
    gammadecay = 1;
    //Next IF command is to shorten the calculations when gamma-emission is the only
    //possible channel
    if(probp==0.0 && probn==0.0 && probd==0.0 && probt==0.0 && proba==0.0 && probhe==0.0 && problamb0==0.0 && probimf==0.0 && probf==0.0)fgamma = 1;
    fimf = 0;
    ff = 0;
  }
  else if (x < proba+probhe+probt+probd+probp+probn+problamb0+probg+probimf) {
    // imf evaporation                                                   
// AIMF and ZIMF obtained from complete procedure (integration over all
// possible Gamma(IMF) and then randomly picked

   G4int iloop=0;
   dir1973:
   imf(af,zf,tcn,ee,&zimf,&aimf,&bimf,&sbimf,&timf,jprf);
   iloop++;
   if(iloop>100)std::cout << "Problem in EVAPORA: IMF called > 100 times" << std::endl;
   if(zimf>=(zf-2.0)) goto dir1973;
   if(zimf>zf/2.0){
          zimf = zf - zimf;
          aimf = af - aimf;
   }
   // These cases should in principle never happen
   if(zimf==0.0 || aimf==0.0 || sbimf>ee)std::cout << "warning: Look in EVAPORA CALL IMF" << std::endl;

// I sample the total kinetic energy consumed by the system of two nuclei
// from the distribution determined with the temperature at saddle point
// TKEIMF is the kinetic energy in the centre of mass of IMF and its partner

   G4int ii=0;
   dir1235:
   tkeimf= fmaxhaz(timf);
   ii++;
   if(ii>100){
   tkeimf=min(2.0*timf,ee-sbimf);
   goto dir1000;
   }
   if(tkeimf<=0.0)goto dir1235;
   if(tkeimf>(ee-sbimf) && timf>0.5)goto dir1235;
   dir1000:
   tkeimf = tkeimf + bimf;

    amoins = aimf;
    zmoins = zimf;
    epsiln = (sbimf-bimf) + tkeimf;
    pc = 0.0;
    malpha = 0.0;
    fgamma = 0; 
    fimf = 1; 
    ff = 0;
    flamb0decay = 0;
    gammadecay = 0;
  }
  else {
    // fission                                                               
    // in case of fission-events the fragment nucleus is the mother nucleus  
    // before fission occurs with excitation energy above the fis.- barrier. 
    // fission fragment mass distribution is calulated in subroutine fisdis  

    amoins = 0.0;
    zmoins = 0.0;
    epsiln = ef;
//
    malpha = 0.0;
    pc = 0.0;
    ff = 1;
    fimf = 0;
    fgamma = 0;
    flamb0decay = 0;
    gammadecay = 0;
  }
//
  direct99:
  if (ee <= 0.01)ee = 0.01;
// Davide Mancusi (DM) - 2010
      if(gammadecay==1 && ee<(epsiln+0.010)){
        epsiln = ee - 0.010;
       // fgamma = 1;
      }

      if(epsiln<0.0){
       std::cout << "***WARNING epsilon<0***" << std::endl;
       //epsiln=0.;
       //PRINT*,IDECAYMODE,IDNINT(AF),IDNINT(ZF),EE,EPSILN
      }
  // calculation of the daughter nucleus                                   
  af = af - amoins;
  zf = zf - zmoins;
  ee = ee - epsiln;
  if (ee <= 0.01)ee = 0.01;
  mtota = mtota + malpha;


  //if(amoins==2 && zmoins==0)std::cout << ee << std::endl;
 

      secondneutron:
      if(amoins==2 && zmoins==0){twon=1;amoins=1;}else{ twon=0;}


// Determination of x,y,z components of momentum from known emission momentum PC
      if(ff==0 && fimf==0){
        //
        if(flamb0decay==1){
        EV_TEMP[IEV_TAB_TEMP][0] = 0.;
        EV_TEMP[IEV_TAB_TEMP][1] = -2;
        EV_TEMP[IEV_TAB_TEMP][5] = 1.;
        }else{
        EV_TEMP[IEV_TAB_TEMP][0] = zmoins;
        EV_TEMP[IEV_TAB_TEMP][1] = amoins;
        EV_TEMP[IEV_TAB_TEMP][5] = 0.;
        }
        rnd = G4AblaRandom::flat();
        ctet1 = 2.0*rnd - 1.0;            // z component: uniform probability between -1 and 1
        stet1 = std::sqrt(1.0 - std::pow(ctet1,2)); // component perpendicular to z
        rnd = G4AblaRandom::flat();
        phi1 = rnd*2.0*3.141592654;       // angle in x-y plane: uniform probability between 0 and 2*pi
        G4double xcv = stet1*std::cos(phi1);// x component
        G4double ycv = stet1*std::sin(phi1);// y component
        G4double zcv = ctet1;               // z component
// In the CM system
        if(gammadecay==0){
// Light particle
           G4double ETOT_LP = std::sqrt(pc*pc + amoins*amoins * mu2);
           if(flamb0decay==1)ETOT_LP = std::sqrt(pc*pc + 1115.683*1115.683);
           EV_TEMP[IEV_TAB_TEMP][2] = c * pc * xcv / ETOT_LP;
           EV_TEMP[IEV_TAB_TEMP][3] = c * pc * ycv / ETOT_LP;
           EV_TEMP[IEV_TAB_TEMP][4] = c * pc * zcv / ETOT_LP;
        }else{
// gamma ray
           EV_TEMP[IEV_TAB_TEMP][2] = pc * xcv;
           EV_TEMP[IEV_TAB_TEMP][3] = pc * ycv;
           EV_TEMP[IEV_TAB_TEMP][4] = pc * zcv;
        }
        G4double VXOUT=0.,VYOUT=0.,VZOUT=0.;
        lorentz_boost(vxeva,vyeva,vleva,
            EV_TEMP[IEV_TAB_TEMP][2],EV_TEMP[IEV_TAB_TEMP][3],
            EV_TEMP[IEV_TAB_TEMP][4],
            &VXOUT,&VYOUT,&VZOUT);
        EV_TEMP[IEV_TAB_TEMP][2] = VXOUT;
        EV_TEMP[IEV_TAB_TEMP][3] = VYOUT;
        EV_TEMP[IEV_TAB_TEMP][4] = VZOUT;
// Heavy residue
        if(gammadecay==0){
        G4double v2 = std::pow(EV_TEMP[IEV_TAB_TEMP][2],2.) +
             std::pow(EV_TEMP[IEV_TAB_TEMP][3],2.) +
             std::pow(EV_TEMP[IEV_TAB_TEMP][4],2.);
        G4double gamma = 1.0/std::sqrt(1.0 - v2 / (c*c));
        G4double etot_lp = amoins*mu * gamma;
        pxeva = pxeva - EV_TEMP[IEV_TAB_TEMP][2] * etot_lp / c;
        pyeva = pyeva - EV_TEMP[IEV_TAB_TEMP][3] * etot_lp / c;
        pleva = pleva - EV_TEMP[IEV_TAB_TEMP][4] * etot_lp / c;
        }else{
// in case of gammas, EV_TEMP contains momentum components and not velocity
        pxeva = pxeva - EV_TEMP[IEV_TAB_TEMP][2];
        pyeva = pyeva - EV_TEMP[IEV_TAB_TEMP][3];
        pleva = pleva - EV_TEMP[IEV_TAB_TEMP][4];
        }
        G4double pteva = std::sqrt(pxeva*pxeva + pyeva*pyeva);
// To be checked:
        G4double etot = std::sqrt ( pleva*pleva + pteva*pteva + af*af * mu2 );
        vxeva = c * pxeva / etot;  // recoil velocity components of residue due to evaporation
        vyeva = c * pyeva / etot;
        vleva = c * pleva / etot;
        IEV_TAB_TEMP = IEV_TAB_TEMP + 1;
      }

  if(twon==1){goto secondneutron;}

  // condition for end of evaporation                                   
  if (zf < 3. || (ff == 1) || (fgamma == 1) || (fimf==1)) {
    goto evapora100;
  }
  goto evapora10;

  evapora100:
  (*zf_par) = zf;
  (*af_par) = af;
  (*ee_par) = ee;
  (*faimf) = aimf;
  (*fzimf) = zimf;
  (*jprfout) = jprf;
  (*tkeimf_par) = tkeimf;
  (*mtota_par) = mtota;
  (*vleva_par) = vleva;
  (*vxeva_par) = vxeva;
  (*vyeva_par) = vyeva;
  (*ff_par) = ff;
  (*fimf_par) = fimf;
  (*inttype_par) = inttype;  
  (*iev_tab_temp_par)= IEV_TAB_TEMP;                                       
  (*inum_par) = inum;
  (*NbLam0_par) = NbLam0;
  return;
}

void G4Abla::direct(G4double zprf, G4double a, G4double ee, G4double jprf, G4double *probp_par, G4double *probd_par, G4double *probt_par, G4double *probn_par, G4double *probhe_par, G4double *proba_par, G4double *probg_par,G4double *probimf_par,G4double *probf_par,G4double *problamb0_par, G4double *ptotl_par, G4double *sn_par, G4double *sbp_par, G4double *sbd_par, G4double *sbt_par, G4double *sbhe_par, G4double *sba_par,G4double *slamb0_par, G4double *ecn_par, G4double *ecp_par, G4double *ecd_par, G4double *ect_par,G4double *eche_par,G4double *eca_par, G4double *ecg_par,  G4double *eclamb0_par, G4double *bp_par, G4double *bd_par, G4double *bt_par, G4double *bhe_par, G4double *ba_par,G4double *sp_par,G4double *sd_par,G4double *st_par,G4double *she_par,G4double *sa_par, G4double *ef_par,G4double *ts1_par, G4int, G4int inum, G4int itest, G4int *sortie, G4double *tcn,G4double *jprfn_par, G4double *jprfp_par, G4double *jprfd_par, G4double *jprft_par, G4double *jprfhe_par, G4double *jprfa_par, G4double *jprflamb0_par, G4double *tsum_par, G4int NbLam0)
{
  G4double probp = (*probp_par);
  G4double probd = (*probd_par);
  G4double probt = (*probt_par);
  G4double probn = (*probn_par);
  G4double probhe = (*probhe_par);
  G4double proba = (*proba_par);
  G4double probg = (*probg_par);
  G4double probimf = (*probimf_par);
  G4double probf = (*probf_par);
  G4double problamb0 = (*problamb0_par);
  G4double ptotl = (*ptotl_par);
  G4double sn = (*sn_par);
  G4double sp = (*sp_par);
  G4double sd = (*sd_par);
  G4double st = (*st_par);
  G4double she = (*she_par);
  G4double sa = (*sa_par);
  G4double slamb0 = 0.0;
  G4double sbp = (*sbp_par);
  G4double sbd = (*sbd_par);
  G4double sbt = (*sbt_par);
  G4double sbhe = (*sbhe_par);
  G4double sba = (*sba_par);
  G4double ecn = (*ecn_par);
  G4double ecp = (*ecp_par);
  G4double ecd = (*ecd_par);
  G4double ect = (*ect_par);
  G4double eche = (*eche_par);
  G4double eca = (*eca_par);
  G4double ecg = (*ecg_par);
  G4double eclamb0 = (*eclamb0_par);
  G4double bp = (*bp_par);
  G4double bd = (*bd_par);
  G4double bt = (*bt_par);
  G4double bhe = (*bhe_par);
  G4double ba = (*ba_par);
  G4double tsum = (*tsum_par);

  // CALCULATION OF PARTICLE-EMISSION PROBABILITIES & FISSION     / 
  // BASED ON THE SIMPLIFIED FORMULAS FOR THE DECAY WIDTH BY      / 
  // MORETTO, ROCHESTER MEETING TO AVOID COMPUTING TIME           / 
  // INTENSIVE INTEGRATION OF THE LEVEL DENSITIES                 / 
  // USES EFFECTIVE COULOMB BARRIERS AND AN AVERAGE KINETIC ENERGY/ 
  // OF THE EVAPORATED PARTICLES                                  / 
  // COLLECTIVE ENHANCMENT OF THE LEVEL DENSITY IS INCLUDED       / 
  // DYNAMICAL HINDRANCE OF FISSION IS INCLUDED BY A STEP FUNCTION/ 
  // APPROXIMATION. SEE A.R. JUNGHANS DIPLOMA THESIS              / 
  // SHELL AND PAIRING STRUCTURES IN THE LEVEL DENSITY IS INCLUDED/ 

  // INPUT:                                                            
  // ZPRF,A,EE  CHARGE, MASS, EXCITATION ENERGY OF COMPOUND     
  // NUCLEUS                                         
  // JPRF       ROOT-MEAN-SQUARED ANGULAR MOMENTUM                           

  // DEFORMATIONS AND G.S. SHELL EFFECTS                               
  // COMMON /ECLD/   ECGNZ,ECFNZ,VGSLD,ALPHA                           

  // ECGNZ - GROUND STATE SHELL CORR. FRLDM FOR A SPHERICAL G.S.       
  // ECFNZ - SHELL CORRECTION FOR THE SADDLE POINT (NOW: == 0)         
  // VGSLD - DIFFERENCE BETWEEN DEFORMED G.S. AND LDM VALUE            
  // ALPHA - ALPHA GROUND STATE DEFORMATION (THIS IS NOT BETA2!)       
  // BETA2 = SQRT((4PI)/5) * ALPHA                             

  //OPTIONS AND PARAMETERS FOR FISSION CHANNEL                        
  //COMMON /FISS/    AKAP,BET,HOMEGA,KOEFF,IFIS,                       
  //                 OPTSHP,OPTXFIS,OPTLES,OPTCOL               
  //
  // AKAP   - HBAR**2/(2* MN * R_0**2) = 10 MEV, R_0 = 1.4 FM          
  // BET    - REDUCED NUCLEAR FRICTION COEFFICIENT IN (10**21 S**-1)   
  // HOMEGA - CURVATURE OF THE FISSION BARRIER = 1 MEV                 
  // KOEFF  - COEFFICIENT FOR THE LD FISSION BARRIER == 1.0            
  // IFIS   - 0/1 FISSION CHANNEL OFF/ON                               
  // OPTSHP - INTEGER SWITCH FOR SHELL CORRECTION IN MASSES/ENERGY     
  //          = 0 NO MICROSCOPIC CORRECTIONS IN MASSES AND ENERGY        
  //          = 1 SHELL ,  NO PAIRING                                    
  //          = 2 PAIRING, NO SHELL                                      
  //          = 3 SHELL AND PAIRING                                      
  // OPTCOL - 0/1 COLLECTIVE ENHANCEMENT SWITCHED ON/OFF               
  // OPTXFIS- 0,1,2 FOR MYERS & SWIATECKI, DAHLINGER, ANDREYEV         
  //                FISSILITY PARAMETER.                                     
  // OPTLES - CONSTANT TEMPERATURE LEVEL DENSITY FOR A,Z > TH-224      
  // OPTCOL - 0/1 COLLECTIVE ENHANCEMENT OFF/ON                        

  // LEVEL DENSITY PARAMETERS                                          
  // COMMON /ALD/    AV,AS,AK,OPTAFAN                                  
  //                 AV,AS,AK - VOLUME,SURFACE,CURVATURE DEPENDENCE OF THE             
  //                            LEVEL DENSITY PARAMETER                                
  // OPTAFAN - 0/1  AF/AN >=1 OR AF/AN ==1                             
  //           RECOMMENDED IS OPTAFAN = 0                              

  // FISSION BARRIERS                                                  
  // COMMON /FB/     EFA                                               
  // EFA    - ARRAY OF FISSION BARRIERS                                


  // OUTPUT: PROBN,PROBP,PROBA,PROBF,PTOTL:                            
  // - EMISSION PROBABILITIES FOR N EUTRON, P ROTON,  A LPHA     
  // PARTICLES, F ISSION AND NORMALISATION                     
  // SN,SBP,SBA: SEPARATION ENERGIES N P A                     
  // INCLUDING EFFECTIVE BARRIERS                              
  // ECN,ECP,ECA,BP,BA                                         
  // - AVERAGE KINETIC ENERGIES (2*T) AND EFFECTIVE BARRIERS     

  G4double bk = 0.0;
  G4double bksp = 0.0;
  G4double bc = 0.0;
  G4int afp = 0;
  G4double het = 0.0;
  G4double at = 0.0;
  G4double bs = 0.0;
  G4double bssp = 0.0;
  G4double bshell = 0.0;
  G4double cf = 0.0;
  G4double defbet = 0.0;
  G4double densa = 0.0;
  G4double denshe = 0.0;
  G4double densg = 0.0;
  G4double densn = 0.0;
  G4double densp = 0.0;
  G4double densd = 0.0;
  G4double denst = 0.0;
  G4double denslamb0 = 0.0;
  G4double eer = 0.0;
  G4double ecor = 0.0;
  G4double ef = 0.0;
  G4double ft = 0.0;
  G4double timf = 0.0;
  G4double qr = 0.0;
  G4double qrcn = 0.0;
  G4double omegap=0.0;
  G4double omegad=0.0;
  G4double omegat=0.0;
  G4double omegahe=0.0;
  G4double omegaa=0.0;
  G4double ga = 0.0;
  G4double ghe = 0.0;
  G4double gf = 0.0;
  G4double gff = 0.0;
  G4double gn = 0.0;
  G4double gp = 0.0;
  G4double gd = 0.0;
  G4double gt = 0.0;
  G4double gg = 0.0;
  G4double glamb0 = 0.0;
  G4double gimf = 0.0;
  G4double gimf3 = 0.0;
  G4double gimf5 = 0.0;
  G4double bimf = 0.0;
  G4double bsimf = 0.0;
  G4double sbimf = 0.0;
  G4double densimf = 0.0;
  G4double defbetimf = 0.0;
  G4double b_imf = 0.0;
  G4double a_imf = 0.0;
  G4double omegaimf = 0.0;
  G4int izimf = 0;
  G4double zimf = 0.0;
  G4double gsum = 0.0;
  G4double gtotal=0.0;
  G4double hbar = 6.582122e-22;
  G4double emin = 0.0;
  G4int il = 0;
  G4int choice_fisspart = 0;
  G4double t_lapse=0.0;
  G4int imaxwell = 0;
  G4int in = 0;
  G4int iz = 0;
  G4int ind = 0;
  G4int izd = 0;
  G4int j = 0;
  G4int k = 0;
  G4double ma1z = 0.0;
  G4double mazz = 0.0;
  G4double ma2z = 0.0;
  G4double ma1z1 = 0.0;
  G4double ma2z1 = 0.0;
  G4double ma3z1 = 0.0;
  G4double ma3z2 = 0.0;
  G4double ma4z2 = 0.0;
  G4double maz = 0.0;
  G4double nt = 0.0;
  G4double pi = 3.1415926535;
  G4double pt = 0.0;
  G4double dt = 0.0;
  G4double tt = 0.0;
  G4double lamb0t = 0.0;
  G4double gtemp = 0.0;
  G4double rdt = 0.0;
  G4double rtt = 0.0;
  G4double rat = 0.0;
  G4double rhet = 0.0;
  G4double refmod = 0.0;
  G4double rnt = 0.0;
  G4double rpt = 0.0;
  G4double rlamb0t = 0.0;
  G4double sbfis = 1.e40;
  G4double segs = 0.0;
  G4double selmax = 0.0;
  G4double tauc = 0.0;
  G4double temp = 0.0;
  G4double ts1 = 0.0;
  G4double xx = 0.0;
  G4double y = 0.0;
  G4double k1 = 0.0;
  G4double omegasp=0.0;
  G4double homegasp=0.0;
  G4double omegags=0.0;
  G4double homegags=0.0;
  G4double pa = 0.0;
  G4double gamma = 0.0;
  G4double gfactor = 0.0;
  G4double bscn;
  G4double bkcn;
  G4double bccn;
  G4double ftcn=0.0;
  G4double mfcd;
  G4double jprfn=jprf;
  G4double jprfp=jprf;
  G4double jprfd=jprf;
  G4double jprft=jprf;
  G4double jprfhe=jprf;
  G4double jprfa=jprf;
  G4double jprflamb0=jprf;
  G4double djprf=0.0;
  G4double dlout=0.0;
  G4double sdlout=0.0;
  G4double iinert=0.0;
  G4double erot=0.0;
  G4double erotn=0.0;
  G4double erotp=0.0;
  G4double erotd=0.0;
  G4double erott=0.0;
  G4double erothe=0.0;
  G4double erota=0.0;
  G4double erotlamb0=0.0;
  G4double erotcn=0.0;
 // G4double ecorcn=0.0;
  G4double imfarg=0.0;
  G4double width_imf=0.0;
  G4int IDjprf=0;
  G4int fimf_allowed=opt->optimfallowed;

  if(itest==1){

  }
  // Switch to calculate Maxwellian distribution of kinetic energies
  imaxwell = 1;
  *sortie = 0;

  // just a change of name until the end of this subroutine                
  eer = ee;
  if (inum == 1) {
    ilast = 1;
  }
  // calculation of masses                                           
  // refmod = 1 ==> myers,swiatecki model                              
  // refmod = 0 ==> weizsaecker model                                  
  refmod = 1;  // Default = 1
//
  if (refmod == 1) {
    mglms(a,zprf,fiss->optshp,&maz);
    mglms(a-1.0,zprf,fiss->optshp,&ma1z);
    mglms(a-2.0,zprf,fiss->optshp,&ma2z);
    mglms(a-1.0,zprf-1.0,fiss->optshp,&ma1z1);
    mglms(a-2.0,zprf-1.0,fiss->optshp,&ma2z1);
    mglms(a-3.0,zprf-1.0,fiss->optshp,&ma3z1);
    mglms(a-3.0,zprf-2.0,fiss->optshp,&ma3z2);
    mglms(a-4.0,zprf-2.0,fiss->optshp,&ma4z2);
  }
  else {
    mglw(a,zprf,&maz);
    mglw(a-1.0,zprf,&ma1z);
    mglw(a-1.0,zprf-1.0,&ma1z1);
    mglw(a-2.0,zprf-1.0,&ma2z1);
    mglw(a-3.0,zprf-1.0,&ma3z1);
    mglw(a-3.0,zprf-2.0,&ma3z2);
    mglw(a-4.0,zprf-2.0,&ma4z2);
  }

  if((a-1.)==3.0 && (zprf-1.0)==2.0) ma1z1=-7.7181660;
  if((a-1.)==4.0 && (zprf-1.0)==2.0) ma1z1=-28.295992;
  
  // separation energies                     
  sn = ma1z - maz;
  sp = ma1z1 - maz;
  sd = ma2z1 - maz - 2.2246;
  st = ma3z1 - maz - 8.481977;
  she = ma3z2 - maz - 7.7181660;
  sa = ma4z2 - maz - 28.295992;
  //
  if(NbLam0>1){
   sn = gethyperbinding(a,zprf,NbLam0)-gethyperbinding(a-1.,zprf,NbLam0);
   sp = gethyperbinding(a,zprf,NbLam0)-gethyperbinding(a-1.,zprf-1.,NbLam0);
   sd = gethyperbinding(a,zprf,NbLam0)-gethyperbinding(a-2.,zprf-1.,NbLam0);
   st = gethyperbinding(a,zprf,NbLam0)-gethyperbinding(a-3.,zprf-1.,NbLam0);
   she = gethyperbinding(a,zprf,NbLam0)-gethyperbinding(a-3.,zprf-2.,NbLam0);
   sa =  gethyperbinding(a,zprf,NbLam0)-gethyperbinding(a-4.,zprf-2.,NbLam0);
   slamb0 = gethyperbinding(a,zprf,NbLam0)-gethyperbinding(a-1.,zprf,NbLam0-1);
  }
   if(NbLam0==1){
   G4double deltasn = sn - (gethyperbinding(a,zprf,0)-gethyperbinding(a-1.,zprf,0));
   G4double deltasp = sp - (gethyperbinding(a,zprf,0)-gethyperbinding(a-1.,zprf-1,0));
   G4double deltasd = sd - (gethyperbinding(a,zprf,0)-gethyperbinding(a-2.,zprf-1,0));
   G4double deltast = st - (gethyperbinding(a,zprf,0)-gethyperbinding(a-3.,zprf-1,0));
   G4double deltashe = she - (gethyperbinding(a,zprf,0)-gethyperbinding(a-3.,zprf-2,0));
   G4double deltasa = sa - (gethyperbinding(a,zprf,0)-gethyperbinding(a-4.,zprf-2,0));

   sn = deltasn + gethyperbinding(a,zprf,NbLam0)-gethyperbinding(a-1.,zprf,NbLam0);
   sp = deltasp + gethyperbinding(a,zprf,NbLam0)-gethyperbinding(a-1.,zprf-1.,NbLam0);
   sd = deltasd + gethyperbinding(a,zprf,NbLam0)-gethyperbinding(a-2.,zprf-1.,NbLam0);
   st = deltast + gethyperbinding(a,zprf,NbLam0)-gethyperbinding(a-3.,zprf-1.,NbLam0);
   she = deltashe + gethyperbinding(a,zprf,NbLam0)-gethyperbinding(a-3.,zprf-2.,NbLam0);
   sa = deltasa + gethyperbinding(a,zprf,NbLam0)-gethyperbinding(a-4.,zprf-2.,NbLam0);
   slamb0 = gethyperseparation(a,zprf,NbLam0);
  }

// coulomb barriers
//Proton
  if (zprf <= 1.0e0 || a <= 1.0e0 || (a-zprf) < 0.0) {
    sbp = 1.0e75;
    bp = 1.0e75;
  }else{
  barrs(idnint(zprf-1.),idnint(a-1.),1,1,&bp,&omegap);
  bp = max(bp,0.1);
  sbp = sp + bp;
  }

//Deuteron
  if (zprf <= 1.0e0 || a <= 2.0e0 || (a-zprf) < 1.0) {
    sbd = 1.0e75;
    bd = 1.0e75;
  }else{
  barrs(idnint(zprf-1.),idnint(a-2.),1,2,&bd,&omegad);
  bd = max(bd,0.1);
  sbd = sd + bd;
  }

//Triton
  if (zprf <= 1.0e0 || a <= 3.0e0 || (a-zprf) < 2.0) {
    sbt = 1.0e75;
    bt = 1.0e75;
  }else{
  barrs(idnint(zprf-1.),idnint(a-3.),1,3,&bt,&omegat);
  bt = max(bt,0.1);
  sbt = st + bt;
  }

//Alpha
  if (a-4.0<=0.0 || zprf<=2.0 || (a-zprf)<2.0) {
    sba = 1.0e+75;
    ba = 1.0e+75;
  }else{
  barrs(idnint(zprf-2.),idnint(a-4.),2,4,&ba,&omegaa);
  ba = max(ba,0.1);
  sba = sa + ba;
  }

//He3
  if (a-3.0 <= 0.0 || zprf<=2.0 || (a-zprf)<1.0) {
    sbhe = 1.0e+75;
    bhe = 1.0e+75;
  }else{
  barrs(idnint(zprf-2.),idnint(a-3.),2,3,&bhe,&omegahe);
  bhe = max(bhe,0.1);
  sbhe = she + bhe;
  }

// Dealing with particle-unbound systems
   emin = dmin1(sba,sbhe,dmin1(sbt,sbhe,dmin1(sn,sbp,sbd)));

   if(emin<=0.0){
   *sortie = 1;
   unbound(sn,sp,sd,st,she,sa,bp,bd,bt,bhe,ba,&probf,&probn,&probp,&probd,&probt,&probhe,&proba,&probimf,&probg,&ecn,&ecp,&ecd,&ect,&eche,&eca);
   goto direct70;
   }
//
  k = idnint(zprf);
  j = idnint(a - zprf);
  if (fiss->ifis > 0) {
    // now ef is calculated from efa that depends on the subroutine
    // barfit which takes into account the modification on the ang. mom.
    // note *** shell correction (ecgnz)
    il = idnint(jprf);
    barfit(k,k+j,il,&sbfis,&segs,&selmax);
    if ((fiss->optshp == 1) || (fiss->optshp == 3)) {
       ef = G4double(sbfis) -  ecld->ecgnz[j][k];
// JLRS - Nov 2016 - Corrected values of fission barriers for actinides
     if(k==90){
      if(mod(j,2)==1){
         ef = ef*(4.5114-2.2687*(a-zprf)/zprf);
      }else{
         ef = ef*(3.3931-1.5338*(a-zprf)/zprf);
      }
     }
     if(k==92){
     if((a-zprf)/zprf>1.52)ef=ef*(1.1222-0.10886*(a-zprf)/zprf)-0.1;
     }
     if(k>=94&&k<=98&&j<158){// Data in this range have been tested
// e-e
        if(mod(j,2)==0&&mod(k,2)==0){
        if(k>=94){ef = ef-(11.54108*(a-zprf)/zprf-18.074);}
        }
// O-O
        if(mod(j,2)==1&&mod(k,2)==1){
        if(k>=95){ef = ef-(14.567*(a-zprf)/zprf-23.266);}
        }
// Odd A
        if(mod(j,2)==0&&mod(k,2)==1){
        if(j>=144){ef = ef-(13.662*(a-zprf)/zprf-21.656);}
        }

        if(mod(j,2)==1&&mod(k,2)==0){
        if(j>=144){ef = ef-(13.662*(a-zprf)/zprf-21.656);}
        }
     }
    } 
    else {
      ef = G4double(sbfis);
    }
//
// TO AVOID NEGATIVE VALUES FOR IMPOSSIBLE NUCLEI
// THE FISSION BARRIER IS SET TO ZERO IF SMALLER THAN ZERO.
//                                    
    if (ef < 0.0)ef = 0.0;
    fb->efa[j][k]=ef;
//
// Hyper-fission barrier
//
    if(NbLam0>0){
     ef = ef + 0.51*(1115.-938.+sn-slamb0)/std::pow(a,2./3.);
    }
//
// Set fission barrier
//
    (*ef_par) = ef;
//
  // calculation of surface and curvature integrals needed to      
  // to calculate the level density parameter at the saddle point
    xx = fissility((k+j),k,NbLam0,sn,slamb0,fiss->optxfis);
    y = 1.00 - xx;
    if(y<0.0) y = 0.0;
    if(y>1.0) y = 1.0;
    bssp = bipol(1,y);
    bksp = bipol(2,y);
  }
  else {
    ef = 1.0e40;
    sbfis = 1.0e40;
    bssp = 1.0;
    bksp = 1.0;
  }
//
// COMPOUND NUCLEUS LEVEL DENSITY
//
//  AK 2007 - Now DENSNIV called with correct BS, BK

  afp = idnint(a);
  iz = idnint(zprf);
  in = afp - iz;
  bshell = ecld->ecgnz[in][iz]- ecld->vgsld[in][iz];
  defbet = ecld->beta2[in][iz];

  iinert = 0.4 * 931.49 * 1.16*1.16 * std::pow(a,5.0/3.0)*(1.0 + 0.5*std::sqrt(5./(4.*pi))*defbet);
  erot = jprf * jprf * 197.328 * 197.328 /(2. * iinert);
  erotcn = erot;

  bsbkbc(a,zprf,&bscn,&bkcn,&bccn);

 // if(ee > erot+emin){
  densniv(a,zprf,ee,0.0,&densg,bshell,bscn,bkcn,&temp,fiss->optshp,fiss->optcol,defbet,&ecor,jprf,0,&qrcn);
  ftcn = temp;
/*
  //ecorcn = ecor;
  }else{
// If EE < EROT, only gamma emission can take place
         probf = 0.0;
   	 probp = 0.0;
   	 probd = 0.0;
   	 probt = 0.0;
  	 probn = 0.0;
   	 probhe = 0.0;
   	 proba = 0.0;
   	 probg = 1.0;
         probimf = 0.0;
//c JLRS 03/2017 - Added this calculation
//C According to A. Ignatyuk, GG :
//C Here BS=BK=1, as this was assumed in the parameterization
         pa = (ald->av)*a + (ald->as)*std::pow(a,2./3.) + (ald->ak)*std::pow(a,1./3.);
         gamma = 2.5 * pa * std::pow(a,-4./3.);
         gfactor = 1.+gamma*ecld->ecgnz[in][iz];
         if(gfactor<=0.){
          gfactor = 0.0;
         }
//
         gtemp = 17.60/(std::pow(a,0.699) * std::sqrt(gfactor));
         ecg = 4.0 * gtemp;
//
         goto direct70;
  }
*/

//  ---------------------------------------------------------------
//        LEVEL DENSITIES AND TEMPERATURES OF THE FINAL STATES
//  ---------------------------------------------------------------
//
//  MVR - in case of charged particle emission temperature
//  comes from random kinetic energy from a Maxwelliam distribution
//  if option imaxwell = 1 (otherwise E=2T)
//
//  AK - LEVEL DENSITY AND TEMPERATURE AT THE SADDLE POINT -> now calculated in the subroutine FISSION_WIDTH
//
//
// LEVEL DENSITY AND TEMPERATURE IN THE NEUTRON DAUGHTER
//
// KHS, AK 2007 - Reduction of angular momentum due to orbital angular momentum of emitted fragment
// JLRS Nov-2016 - Added these caculations in abla++

  if (in >= 2) {
    ind=idnint(a)-idnint(zprf)-1;
    izd=idnint(zprf);
    if(jprf>0.10){
         lorb(a,a-1.,jprf,ee-sn,&dlout,&sdlout);
         djprf = gausshaz(1,dlout,sdlout);
         if(IDjprf==1) djprf = 0.0;
         jprfn = jprf + djprf;
         jprfn = dint(std::abs(jprfn));   // The nucleus just turns the other way around
    }
    bshell = ecld->ecgnz[ind][izd] - ecld->vgsld[ind][izd];
    defbet = ecld->beta2[ind][izd];

    iinert = 0.4 * 931.49 * 1.16*1.16 * std::pow(a-1.,5.0/3.0)*(1.0 + 0.5*std::sqrt(5./(4.*pi))*defbet);
    erotn = jprfn * jprfn * 197.328 * 197.328 /(2. * iinert);
    bsbkbc(a-1.,zprf,&bs,&bk,&bc);   

    // level density and temperature in the neutron daughter                 
    densniv(a-1.0,zprf,ee,sn,&densn,bshell, bs,bk,&temp,fiss->optshp,fiss->optcol,defbet,&ecor,jprfn,0,&qr);
    nt = temp;
    ecn=0.0;
    if(densn>0.){
     G4int IS=0;
     if(imaxwell == 1){
      rnt = nt;
      dir1234:
      ecn=fvmaxhaz_neut(rnt);
      IS++;
      if(IS>100){std::cout << "WARNING: FVMAXHAZ_NEUT CALLED MORE THAN 100 TIMES" << std::endl;
      goto exi1000;
      }
       if(ecn>(ee-sn)){
        if((ee-sn)<rnt)
           ecn = ee-sn;
        else
           goto dir1234;
           }
        if(ecn<=0.0) goto dir1234;
     }else{
      ecn = 2.0 * nt;
     }
    }
  } 
  else {
    densn = 0.0;
    ecn = 0.0;
    nt = 0.0;
  }
  exi1000:

// LEVEL DENSITY AND TEMPERATURE IN THE PROTON DAUGHTER
//
// Reduction of angular momentum due to orbital angular momentum of emitted fragment
  if (iz >= 2) {
    ind=idnint(a)-idnint(zprf);
    izd=idnint(zprf)-1;
    if(jprf>0.10){
         lorb(a,a-1.,jprf,ee-sbp,&dlout,&sdlout);
         djprf = gausshaz(1,dlout,sdlout);
         if(IDjprf==1) djprf = 0.0;
         jprfp = jprf + djprf;
         jprfp = dint(std::abs(jprfp));   // The nucleus just turns the other way around
    }
    bshell = ecld->ecgnz[ind][izd] - ecld->vgsld[ind][izd];
    defbet =ecld->beta2[ind][izd];

    iinert = 0.4 * 931.49 * 1.16*1.16 * std::pow(a-1.,5.0/3.0)*(1.0 + 0.5*std::sqrt(5./(4.*pi))*defbet);
    erotp = jprfp * jprfp * 197.328 * 197.328 /(2. * iinert);

    bsbkbc(a-1.,zprf-1.,&bs,&bk,&bc);
   
    // level density and temperature in the proton daughter                  
    densniv(a-1.0,zprf-1.0,ee,sbp,&densp,bshell,bs,bk,&temp,fiss->optshp,fiss->optcol,defbet,&ecor,jprfp,0,&qr);
    pt = temp;
    ecp = 0.;
    if(densp>0.){
     G4int IS=0;
     if(imaxwell == 1){
      rpt = pt;
      dir1235:
      ecp=fvmaxhaz(rpt);
      IS++;
      if(IS>100){std::cout << "WARNING: FVMAXHAZ CALLED MORE THAN 100 TIMES" << std::endl;
      goto exi1001;
      }
       if(ecp>(ee-sbp)){
        if((ee-sbp)<rpt)
           ecp = ee-sbp;
        else
           goto dir1235;
           }
        if(ecp<=0.0) goto dir1235;
      ecp = ecp + bp;
     }else{
      ecp = 2.0 * pt + bp;
     }
    }
  }
  else {
    densp = 0.0;
    ecp = 0.0;
    pt = 0.0;
  }
  exi1001:

//  FINAL LEVEL DENSITY AND TEMPERATURE AFTER DEUTERON EMISSION
//
// Reduction of angular momentum due to orbital angular momentum of emitted fragment
  if ((in >= 2) && (iz >= 2)) {
    ind=idnint(a)-idnint(zprf)-1;
    izd=idnint(zprf)-1;
    if(jprf>0.10){
         lorb(a,a-2.,jprf,ee-sbd,&dlout,&sdlout);
         djprf = gausshaz(1,dlout,sdlout);
         if(IDjprf==1) djprf = 0.0;
         jprfd = jprf + djprf;
         jprfd = dint(std::abs(jprfd));   // The nucleus just turns the other way around
    }
    bshell = ecld->ecgnz[ind][izd] - ecld->vgsld[ind][izd];
    defbet = ecld->beta2[ind][izd];

    iinert = 0.4 * 931.49 * 1.16*1.16 * std::pow(a-2.,5.0/3.0)*(1.0 + 0.5*std::sqrt(5./(4.*pi))*defbet);
    erotd = jprfd * jprfd * 197.328 * 197.328 /(2. * iinert);

    bsbkbc(a-2.,zprf-1.,&bs,&bk,&bc);

    // level density and temperature in the deuteron daughter                   
    densniv(a-2.0,zprf-1.0e0,ee,sbd,&densd,bshell,bs,bk,&temp,fiss->optshp,fiss->optcol,defbet,&ecor,jprfd,0,&qr);

    dt = temp;
    ecd = 0.0;
    if(densd>0.){
     G4int IS=0;
     if(imaxwell == 1){
      rdt = dt;
      dir1236:
      ecd=fvmaxhaz(rdt);
      IS++;
      if(IS>100){std::cout << "WARNING: FVMAXHAZ CALLED MORE THAN 100 TIMES" << std::endl;
      goto exi1002;
      }
       if(ecd>(ee-sbd)){
        if((ee-sbd)<rdt)
           ecd = ee-sbd;
        else
           goto dir1236;
           }
        if(ecd<=0.0) goto dir1236;
      ecd = ecd + bd;
     }else{
      ecd = 2.0 * dt + bd;
     }
    }
  }
  else {
    densd = 0.0;
    ecd = 0.0;
    dt = 0.0;
  }
  exi1002:

//  FINAL LEVEL DENSITY AND TEMPERATURE AFTER TRITON EMISSION
//
// Reduction of angular momentum due to orbital angular momentum of emitted fragment
  if ((in >= 3) && (iz >= 2)) {
    ind=idnint(a)-idnint(zprf)-2;
    izd=idnint(zprf)-1;
    if(jprf>0.10){
         lorb(a,a-3.,jprf,ee-sbt,&dlout,&sdlout);
         djprf = gausshaz(1,dlout,sdlout);
         if(IDjprf==1) djprf = 0.0;
         jprft = jprf + djprf;
         jprft = dint(std::abs(jprft));   // The nucleus just turns the other way around
    }
    bshell = ecld->ecgnz[ind][izd] - ecld->vgsld[ind][izd];
    defbet = ecld->beta2[ind][izd];

    iinert = 0.4 * 931.49 * 1.16*1.16 * std::pow(a-3.,5.0/3.0)*(1.0 + 0.5*std::sqrt(5./(4.*pi))*defbet);
    erott = jprft * jprft * 197.328 * 197.328 /(2. * iinert);

    bsbkbc(a-3.,zprf-1.,&bs,&bk,&bc);

    // level density and temperature in the triton daughter                   
    densniv(a-3.0,zprf-1.0,ee,sbt,&denst,bshell,bs,bk,&temp,fiss->optshp,fiss->optcol,defbet,&ecor,jprft,0,&qr);

    tt = temp;
    ect=0.;
    if(denst>0.){
     G4int IS=0;
     if(imaxwell == 1){
      rtt = tt;
      dir1237:
      ect=fvmaxhaz(rtt);
      IS++;
      if(IS>100){std::cout << "WARNING: FVMAXHAZ CALLED MORE THAN 100 TIMES" << std::endl;
      goto exi1003;
      }
       if(ect>(ee-sbt)){
        if((ee-sbt)<rtt)
           ect = ee-sbt;
        else
           goto dir1237;
           }
        if(ect<=0.0) goto dir1237;
      ect = ect + bt;
     }else{
      ect = 2.0 * tt + bt;
     }
    }
  }
  else {
    denst = 0.0;
    ect = 0.0;
    tt = 0.0;
  }
  exi1003: 

// LEVEL DENSITY AND TEMPERATURE IN THE ALPHA DAUGHTER
//
// Reduction of angular momentum due to orbital angular momentum of emitted fragment
  if ((in >= 3) && (iz >= 3)) {
    ind=idnint(a)-idnint(zprf)-2;
    izd=idnint(zprf)-2;
    if(jprf>0.10){
         lorb(a,a-4.,jprf,ee-sba,&dlout,&sdlout);
         djprf = gausshaz(1,dlout,sdlout);
         if(IDjprf==1) djprf = 0.0;
         jprfa = jprf + djprf;
         jprfa = dint(std::abs(jprfa));   // The nucleus just turns the other way around
    }
    bshell = ecld->ecgnz[ind][izd] - ecld->vgsld[ind][izd];
    defbet = ecld->beta2[ind][izd];

    iinert = 0.4 * 931.49 * 1.16*1.16 * std::pow(a-4.,5.0/3.0)*(1.0 + 0.5*std::sqrt(5./(4.*pi))*defbet);
    erota = jprfa * jprfa * 197.328 * 197.328 /(2. * iinert);

    bsbkbc(a-4.,zprf-2.,&bs,&bk,&bc);

    // level density and temperature in the alpha daughter                   
    densniv(a-4.0,zprf-2.0,ee,sba,&densa,bshell,bs,bk,&temp,fiss->optshp,fiss->optcol,defbet,&ecor,jprfa,0,&qr);

    at = temp;
    eca = 0.0;
    if(densa>0.){
     G4int IS=0;
     if(imaxwell == 1){
      rat = at;
      dir1238:
      eca=fvmaxhaz(rat);
      IS++;
      if(IS>100){std::cout << "WARNING: FVMAXHAZ CALLED MORE THAN 100 TIMES" << std::endl;
      goto exi1004;
      }
       if(eca>(ee-sba)){
        if((ee-sba)<rat)
           eca = ee-sba;
        else
           goto dir1238;
           }
        if(eca<=0.0) goto dir1238;
      eca = eca + ba;
     }else{
      eca = 2.0 * at + ba;
     }
    }
  }
  else {
    densa = 0.0;
    eca = 0.0;
    at = 0.0;
  }
  exi1004: 

//  FINAL LEVEL DENSITY AND TEMPERATURE AFTER 3HE EMISSION
//
// Reduction of angular momentum due to orbital angular momentum of emitted fragment
  if ((in >= 2) && (iz >= 3)) {
    ind=idnint(a)-idnint(zprf)-1;
    izd=idnint(zprf)-2;
    if(jprf>0.10){
         lorb(a,a-3.,jprf,ee-sbhe,&dlout,&sdlout);
         djprf = gausshaz(1,dlout,sdlout);
         if(IDjprf==1) djprf = 0.0;
         jprfhe = jprf + djprf;
         jprfhe = dint(std::abs(jprfhe));   // The nucleus just turns the other way around
    }
    bshell = ecld->ecgnz[ind][izd] - ecld->vgsld[ind][izd];
    defbet = ecld->beta2[ind][izd];

    iinert = 0.4 * 931.49 * 1.16*1.16 * std::pow(a-3.,5.0/3.0)*(1.0 + 0.5*std::sqrt(5./(4.*pi))*defbet);
    erothe = jprfhe * jprfhe * 197.328 * 197.328 /(2. * iinert);

    bsbkbc(a-3.,zprf-2.,&bs,&bk,&bc);

    // level density and temperature in the he3 daughter                   
    densniv(a-3.0,zprf-2.0,ee,sbhe,&denshe,bshell,bs,bk,&temp,fiss->optshp,fiss->optcol,defbet,&ecor,jprfhe,0,&qr);

    het = temp;
    eche = 0.0;
    if(denshe>0.){
     G4int IS=0;
     if(imaxwell == 1){
      rhet = het;
      dir1239:
      eche=fvmaxhaz(rhet);
      IS++;
      if(IS>100){std::cout << "WARNING: FVMAXHAZ CALLED MORE THAN 100 TIMES" << std::endl;
      goto exi1005;
      }
       if(eche>(ee-sbhe)){
        if((ee-sbhe)<rhet)
           eche = ee-sbhe;
        else
           goto dir1239;
           }
        if(eche<=0.0) goto dir1239;
      eche = eche + bhe;
     }else{
      eche = 2.0 * het + bhe;
     }
    }
  }
  else {
    denshe = 0.0;
    eche = 0.0;
    het = 0.0;
  }
  exi1005:

// LEVEL DENSITY AND TEMPERATURE IN THE LAMBDA0 DAUGHTER
//
// - Reduction of angular momentum due to orbital angular momentum of emitted fragment
// JLRS Jun-2017 - Added these caculations in abla++

  if (in >= 2 && NbLam0>0) {
    ind=idnint(a)-idnint(zprf)-1;
    izd=idnint(zprf);
    if(jprf>0.10){
         lorb(a,a-1.,jprf,ee-slamb0,&dlout,&sdlout);
         djprf = gausshaz(1,dlout,sdlout);
         if(IDjprf==1) djprf = 0.0;
         jprflamb0 = jprf + djprf;
         jprflamb0 = dint(std::abs(jprflamb0));   // The nucleus just turns the other way around
    }
    bshell = ecld->ecgnz[ind][izd] - ecld->vgsld[ind][izd];
    defbet = ecld->beta2[ind][izd];

    iinert = 0.4 * 931.49 * 1.16*1.16 * std::pow(a-1.,5.0/3.0)*(1.0 + 0.5*std::sqrt(5./(4.*pi))*defbet);
    erotlamb0 = jprflamb0 * jprflamb0 * 197.328 * 197.328 /(2. * iinert);
    bsbkbc(a-1.,zprf,&bs,&bk,&bc);   

    // level density and temperature in the neutron daughter                 
    densniv(a-1.0,zprf,ee,slamb0,&denslamb0,bshell, bs,bk,&temp,fiss->optshp,fiss->optcol,defbet,&ecor,jprflamb0,0,&qr);
    lamb0t = temp;
    eclamb0=0.0;
    if(denslamb0>0.){
     G4int IS=0;
     if(imaxwell == 1){
      rlamb0t = lamb0t;
      dir1240:
      eclamb0=fvmaxhaz_neut(rlamb0t);
      IS++;
      if(IS>100){std::cout << "WARNING: FVMAXHAZ_NEUT CALLED MORE THAN 100 TIMES" << std::endl;
      goto exi1006;
      }
       if(eclamb0>(ee-slamb0)){
        if((ee-slamb0)<rlamb0t)
           eclamb0 = ee-slamb0;
        else
           goto dir1240;
           }
        if(eclamb0<=0.0) goto dir1240;
     }else{
      eclamb0 = 2.0 * lamb0t;
     }
    }
  } 
  else {
    denslamb0 = 0.0;
    eclamb0 = 0.0;
    lamb0t = 0.0;
  }
  exi1006:



// Decay widths for particles
  if ( densg > 0.) {
//
// CALCULATION OF THE PARTIAL DECAY WIDTH
// USED FOR BOTH THE TIME SCALE AND THE EVAPORATION DECAY WIDTH
//
//      AKAP = HBAR**2/(2* MN * R_0**2) = 10 MEV    *** input param ***
//
// AK, KHS 2005 - Energy-dependen inverse cross sections included, influence of
//                Coulomb barrier for LCP, tunnelling for LCP
// JLRS 2017 - Implementation in abla++ 

       if(densn<=0.0){
        gn = 0.0;
       }else{
        gn  = width(a,zprf,1.0,0.0,nt,0.0,sn,ee-erotn)* densn/densg;
       }
       if(densp<=0.0){
        gp = 0.0;
       }else{
        gp  = width(a,zprf,1.0,1.0,pt,bp,sbp,ee-erotp)*densp/densg* pen(a, 1.0, omegap, pt);
       }
       if(densd<=0.0){
        gd = 0.0;
       }else{
        gd  = width(a,zprf,2.0,1.0,dt,bd,sbd,ee-erotd)*densd/densg* pen(a, 2.0, omegad, dt);
       }
       if(denst<=0.0){
        gt = 0.0;
       }else{
        gt  = width(a,zprf,3.0,1.0,tt,bt,sbt,ee-erott)*denst/densg* pen(a, 3.0, omegat, tt);
       }
       if(denshe<=0.0){
        ghe = 0.0;
       }else{
        ghe =width(a,zprf,3.0,2.0,het,bhe,sbhe,ee-erothe)  * denshe/densg* pen(a, 3.0, omegahe, het);
       }
       if(densa<=0.0){
        ga = 0.0;
       }else{
        ga  = width(a,zprf,4.0,2.0,at,ba,sba,ee-erota)  * densa/densg* pen(a, 4.0, omegaa, at);
       }
       if(denslamb0<=0.0){
        glamb0 = 0.0;
       }else{
        glamb0  = width(a,zprf,1.0,-2.0,lamb0t,0.0,slamb0,ee-erotlamb0)* denslamb0/densg;
       }

//     **************************
//     *  Treatment of IMFs     *
//     * KHS, AK, MVR 2005-2006 *
//     **************************

       G4int izcn=0,incn=0,inmin=0,inmax=0,inmi=0,inma=0;
       G4double aimf,mares,maimf;

       if(fimf_allowed==0 || zprf<=5.0 || a<=7.0){
        gimf = 0.0;
       }else{
//      Estimate the total decay width for IMFs (Z >= 3)
//      By using the logarithmic slope between GIMF3 and GIMF5

        mglms(a,zprf,opt->optshpimf,&mazz);

        gimf3 = 0.0;
        zimf = 3.0;
        izimf = 3;
//      *** Find the limits that both IMF and partner are bound :
        izcn = idnint(zprf);                  // Z of CN
        incn = idnint(a) - izcn;              // N of CN

        isostab_lim(izimf,&inmin,&inmax);     // Bound isotopes for IZIMF from INMIN to INIMFMA
        isostab_lim(izcn-izimf,&inmi,&inma);  // Daughter nucleus after IMF emission,
                                           //     limits of bound isotopes
        inmin = max(inmin,incn-inma);      //     Both IMF and daughter must be bound
        inmax = min(inmax,incn-inmi);      //        "

        inmax = max(inmax,inmin);          // In order to keep the variables below

        for(G4int iaimf=izimf+inmin;iaimf<=izimf+inmax;iaimf++){
         aimf=G4double(iaimf);
         if(aimf>=a || zimf>=zprf){
          width_imf = 0.0;
         }else{
          // Q-values
          mglms(a-aimf,zprf-zimf,opt->optshpimf,&mares);
          mglms(aimf,zimf,opt->optshpimf,&maimf);
          // Bass barrier
          barrs(idnint(zprf-zimf),idnint(a-aimf),izimf,idnint(aimf),&bimf,&omegaimf);
          sbimf = maimf+mares-mazz+bimf+getdeltabinding(a,NbLam0);
          // Rotation energy
          defbetimf = ecld->beta2[idnint(aimf-zimf)][idnint(zimf)]+ecld->beta2[idnint(a-aimf-zprf+zimf)][idnint(zprf-zimf)];

          iinert= 0.40 * 931.490 * 1.160*1.160 * std::pow(a,5.0/3.0)*(std::pow(aimf,5.0/3.0) + std::pow(a - aimf,5.0/3.0)) + 931.490 * 1.160*1.160 * aimf * (a-aimf) / a *(std::pow(aimf,1.0/3.0) + std::pow(a - aimf,1.0/3.0))*(std::pow(aimf,1.0/3.0) + std::pow(a - aimf,1.0/3.0));

          erot = jprf * jprf * 197.328 * 197.328 /(2.0 * iinert);

          // Width
          if(densg==0.0 || ee < (sbimf + erot)){
           width_imf = 0.0;
          }else{
          // To take into account that at the barrier the system is deformed:
          //      BSIMF = ((A-AIMF)**(2.D0/3.D0) + AIMF**(2.D0/3.D0))/A**(2.D0/3.D0)
           bsimf = bscn;
           densniv(a,zprf,ee,sbimf,&densimf,0.0,bsimf,1.0,&timf,0,0,defbetimf,&ecor,jprf,2,&qr);

           imfarg = (sbimf+erotcn-erot)/timf;
           if(imfarg > 200.0) imfarg = 200.0;

// For IMF - The available phase space is given by the level densities in CN at the
// barrier; applaying MOrretto -> G=WIDTH*ro_CN(E-SBIMF)/ro_CN(E).
// Constant temperature approximation: ro(E+dE)/ro(E)=exp(dE/T)
// Ratio  DENSIMF/DENSCN is included to take into account that at the barrier system
// is deformed. If (above) BSIMF = 1 no deformation is considered and this ratio
// is equal to 1.
           width_imf = 0.0;
         //
           width_imf = width(a,zprf,aimf,zimf,timf,bimf,sbimf,ee-erot)*std::exp(-imfarg)*qr/qrcn;
          }// if densg
         }// if aimf
         gimf3 = gimf3 + width_imf;
        }// for IAIMF

//   zimf = 5       
        gimf5 = 0.0;
        zimf = 5.0;
        izimf = 5;
//      *** Find the limits that both IMF and partner are bound :
        izcn = idnint(zprf);                  // Z of CN
        incn = idnint(a) - izcn;              // N of CN

        isostab_lim(izimf,&inmin,&inmax);     // Bound isotopes for IZIMF from INMIN to INIMFMA
        isostab_lim(izcn-izimf,&inmi,&inma);  // Daughter nucleus after IMF emission,
                                           //     limits of bound isotopes
        inmin = max(inmin,incn-inma);      //     Both IMF and daughter must be bound
        inmax = min(inmax,incn-inmi);      //        "

        inmax = max(inmax,inmin);          // In order to keep the variables below

        for(G4int iaimf=izimf+inmin;iaimf<=izimf+inmax;iaimf++){
         aimf=G4double(iaimf);
         if(aimf>=a || zimf>=zprf){
          width_imf = 0.0;
         }else{
          // Q-values
          mglms(a-aimf,zprf-zimf,opt->optshpimf,&mares);
          mglms(aimf,zimf,opt->optshpimf,&maimf);
          // Bass barrier
          barrs(idnint(zprf-zimf),idnint(a-aimf),izimf,idnint(aimf),&bimf,&omegaimf);
          sbimf = maimf+mares-mazz+bimf+getdeltabinding(a,NbLam0);
          // Rotation energy
          defbetimf = ecld->beta2[idnint(aimf-zimf)][idnint(zimf)]+ecld->beta2[idnint(a-aimf-zprf+zimf)][idnint(zprf-zimf)];

          iinert= 0.40 * 931.490 * 1.160*1.160 * std::pow(a,5.0/3.0)*(std::pow(aimf,5.0/3.0) + std::pow(a - aimf,5.0/3.0)) + 931.490 * 1.160*1.160 * aimf * (a-aimf) / a *(std::pow(aimf,1.0/3.0) + std::pow(a - aimf,1.0/3.0))*(std::pow(aimf,1.0/3.0) + std::pow(a - aimf,1.0/3.0));

          erot = jprf * jprf * 197.328 * 197.328 /(2.0 * iinert);
//
          // Width
          if(densg==0.0 || ee < (sbimf + erot)){
           width_imf = 0.0;
          }else{
          // To take into account that at the barrier the system is deformed:
          //      BSIMF = ((A-AIMF)**(2.D0/3.D0) + AIMF**(2.D0/3.D0))/A**(2.D0/3.D0)
           bsimf = bscn;
           densniv(a,zprf,ee,sbimf,&densimf,0.0,bsimf,1.0,&timf,0,0,defbetimf,&ecor,jprf,2,&qr);
//
           imfarg = (sbimf+erotcn-erot)/timf;
           if(imfarg > 200.0) imfarg = 200.0;
//
// For IMF - The available phase space is given by the level densities in CN at the
// barrier; applaying MOrretto -> G=WIDTH*ro_CN(E-SBIMF)/ro_CN(E).
// Constant temperature approximation: ro(E+dE)/ro(E)=exp(dE/T)
// Ratio  DENSIMF/DENSCN is included to take into account that at the barrier system
// is deformed. If (above) BSIMF = 1 no deformation is considered and this ratio
// is equal to 1.
           width_imf = 0.0;
           width_imf = width(a,zprf,aimf,zimf,timf,bimf,sbimf,ee-erot)*std::exp(-imfarg)*qr/qrcn;//*densimf/densg;
          }// if densg
         }// if aimf
         gimf5 = gimf5 + width_imf;
        }// for IAIMF
// It is assumed that GIMFi = A_IMF*ZIMF**B_IMF; to get the total GIMF one integrates
// Int(A_IMF*ZIMF**B_IMF)(3->ZPRF)

        if(gimf3<=0.0 || gimf5<=0.0){
         gimf = 0.0;
         b_imf = -100.0;
         a_imf = 0.0;
        }else{
//
        b_imf = (std::log10(gimf3) - std::log10(gimf5))/(std::log10(3.0)-std::log10(5.0));
//
         if(b_imf >= -1.01) b_imf = -1.01;
         if(b_imf <= -100.0) {
          b_imf = -100.0;
          a_imf = 0.0;
          gimf = 0.0;
          goto direct2007;
         }
//
         a_imf = gimf3 / std::pow(3.0,b_imf);
         gimf = a_imf * ( std::pow(zprf,b_imf+1.0) - std::pow(3.0,b_imf+1.0)) /(b_imf + 1.0);
        }

       direct2007:
       if(gimf < 1.e-10) gimf = 0.0;
       }// if fimf_allowed
//
//c JLRS 2016 - Added this calculation
//C AK 2004 - Gamma width
//C According to A. Ignatyuk, GG :
//C Here BS=BK=1, as this was assumed in the parameterization
      pa = (ald->av)*a + (ald->as)*std::pow(a,2./3.) + (ald->ak)*std::pow(a,1./3.);
      gamma = 2.5 * pa * std::pow(a,-4./3.);
      gfactor = 1.+gamma*ecld->ecgnz[in][iz];
      if(gfactor<=0.){
       gfactor = 0.0;
      }
//
      gtemp = 17.60/(std::pow(a,0.699) * std::sqrt(gfactor));
//
//C If one switches gammas off, one should also switch off tunneling through the fission barrier.
      gg = 0.624e-9*std::pow(a,1.6)*std::pow(gtemp,5.);
//gammaemission==1
//C For fission fragments, GG is ~ 2 times larger than for
//c "oridnary" nuclei (A. Ignatyuk, private communication).
      if(gammaemission==1){
      gg = 2.0 * gg;
      }
      ecg = 4.0 * gtemp;
//
//
  gsum = ga + ghe + gd + gt + gp + gn + gimf + gg + glamb0;

  //std::cout << gn << " " << gd << " " << gp << std::endl;

  if (gsum > 0.0) {
    ts1  = hbar / gsum;
  }
  else {
    ts1  = 1.0e99;
    goto direct69;
  }
//
//Case of nuclei below Businaro-Gallone mass asymmetry point
    if(fiss->ifis==0 || (zprf*zprf/a<=22.74 && zprf<60.)){
     goto direct69;
    }
//
// Calculation of the fission decay width
// Deformation is calculated using the fissility  
//
    defbet = y;
    fission_width(zprf,a,ee,bssp,bksp,ef,y,&gf,&temp,jprf,0,1,fiss->optcol,fiss->optshp,densg);
    ft=temp;
//
// Case of very heavy nuclei that have no fission barrier
// For them fission is the only decay channel available
       if(ef<=0.0){
         probf = 1.0;
   	 probp = 0.0;
   	 probd = 0.0;
   	 probt = 0.0;
  	 probn = 0.0;
   	 probhe = 0.0;
   	 proba = 0.0;
   	 probg = 0.0;
         probimf = 0.0;
         problamb0 = 0.0;
         goto direct70;
       }

       if(fiss->bet<=0.){
        gtotal = ga + ghe + gp + gd + gt + gn + gg +gimf + gf + glamb0;
        if(gtotal<=0.0){
         probf = 0.0;
   	 probp = 0.0;
   	 probd = 0.0;
   	 probt = 0.0;
  	 probn = 0.0;
   	 probhe = 0.0;
   	 proba = 0.0;
   	 probg = 0.0;
         probimf = 0.0;
         problamb0 = 0.0;
         goto direct70;
        }else{
         probf = gf/gtotal;
         probn = gn/gtotal;
         probp = gp/gtotal;
         probd = gd/gtotal;
         probt = gt/gtotal;
         probhe = ghe/gtotal;
         proba = ga/gtotal;
         probg = gg/gtotal;
         probimf = gimf/gtotal;
         problamb0 = glamb0/gtotal;
         goto direct70;
        }
       }
  }else{
   goto direct69;
  }
//
  if (inum > ilast) {  // new event means reset the time scale
    tsum = 0.;
  }
//
// kramers factor for the dynamical hindrances of fission
  fomega_sp(a,y,&mfcd,&omegasp,&homegasp);
  cf = cram(fiss->bet,homegasp);
//
// We calculate the transient time
  fomega_gs(a,zprf,&k1,&omegags,&homegags);
  tauc=tau(fiss->bet,homegags,ef,ft);
  gf=gf*cf;
//
/*
c The subroutine part_fiss calculates the fission width GFF that corresponds to the time
c dependence of the probability distribution obtained by solving the FOKKER-PLANCK eq
c using a nucleus potential that is approximated by a parabola. It also gives the
c decay time for this step T_LAPSE that includes all particle decay channels and the
c fission channel. And it decides whether the nucleus decays by particle evaporation
c CHOICE_FISSPART = 1 or fission CHOICE_FISSPART = 2
*/
//
 part_fiss(fiss->bet,gsum,gf,y,tauc,ts1,tsum, &choice_fisspart,zprf,a,ft,&t_lapse,&gff);
 gf = gff;
//
// We accumulate in TSUM the mean decay for this step including all particle decay channels and fission
 tsum = tsum + t_lapse;

//   If fission occurs
    if(choice_fisspart==2){
	  probf = 1.0;
	  probp = 0.0;
	  probd = 0.0;
	  probt = 0.0;
	  probn = 0.0;
	  probhe = 0.0;
	  proba = 0.0;
	  probg = 0.0;
          probimf = 0.0;
          problamb0 = 0.0;
          goto direct70;
    }else{
// If particle evaporation occurs
// The probabilities for the different decays are calculated taking into account the fission width GFF that corresponds to this step

       gtotal=ga + ghe + gp + gd + gt + gn + gimf + gg + glamb0;
        if(gtotal<=0.0){
         probf = 0.0;
   	 probp = 0.0;
   	 probd = 0.0;
   	 probt = 0.0;
  	 probn = 0.0;
   	 probhe = 0.0;
   	 proba = 0.0;
   	 probg = 0.0;
         probimf = 0.0;
         problamb0 = 0.0;
         goto direct70;
       }else{
         probf = 0.0;
         probn = gn/gtotal;
         probp = gp/gtotal;
         probd = gd/gtotal;
         probt = gt/gtotal;
         probhe = ghe/gtotal;
         proba = ga/gtotal;
         probg = gg/gtotal;
         probimf = gimf/gtotal;
         problamb0 = glamb0/gtotal;
         goto direct70;
        }
    }
//
  direct69:
        gtotal = ga + ghe + gp + gd + gt + gn + gg + gimf + glamb0;
        if(gtotal<=0.0){
         probf = 0.0;
   	 probp = 0.0;
   	 probd = 0.0;
   	 probt = 0.0;
  	 probn = 0.0;
   	 probhe = 0.0;
   	 proba = 0.0;
   	 probg = 0.0;
         probimf = 0.0;
         problamb0 = 0.0;
       }else{
         probf = 0.0;
         probn = gn/gtotal;
         probp = gp/gtotal;
         probd = gd/gtotal;
         probt = gt/gtotal;
         probhe = ghe/gtotal;
         proba = ga/gtotal;
         probg = gg/gtotal;
         probimf = gimf/gtotal;
         problamb0 = glamb0/gtotal;
        }

  direct70:
  ptotl = probp+probd+probt+probn+probhe+proba+probg+probimf+probf+problamb0;  
  //
  ee = eer;
  ilast = inum;

  // Return values:
  (*probp_par) = probp;
  (*probd_par) = probd;
  (*probt_par) = probt;
  (*probn_par) = probn;
  (*probhe_par) = probhe;
  (*proba_par) = proba;
  (*probg_par) = probg;
  (*probimf_par) = probimf;
  (*problamb0_par) = problamb0;
  (*probf_par) = probf;
  (*ptotl_par) = ptotl;
  (*sn_par) = sn;
  (*sp_par) = sp;
  (*sd_par) = sd;
  (*st_par) = st;
  (*she_par) = she;
  (*sa_par) = sa;
  (*slamb0_par) = slamb0;
  (*sbp_par) = sbp;
  (*sbd_par) = sbd;
  (*sbt_par) = sbt;
  (*sbhe_par) = sbhe;
  (*sba_par) = sba;
  (*ecn_par) = ecn;
  (*ecp_par) = ecp;
  (*ecd_par) = ecd;
  (*ect_par) = ect;
  (*eche_par) = eche;
  (*eca_par) = eca;
  (*ecg_par) = ecg;
  (*eclamb0_par) = eclamb0;
  (*bp_par) = bp;
  (*bd_par) = bd;
  (*bt_par) = bt;
  (*bhe_par) = bhe;
  (*ba_par) = ba;
  (*tcn) = ftcn;
  (*ts1_par) = ts1;
  (*jprfn_par) = jprfn;
  (*jprfp_par) = jprfp;
  (*jprfd_par) = jprfd;
  (*jprft_par) = jprft;
  (*jprfhe_par) = jprfhe;
  (*jprfa_par) = jprfa;
  (*jprflamb0_par) = jprflamb0;
  (*tsum_par) = tsum;
  return;
}

void G4Abla::densniv(G4double a, G4double z, G4double ee, G4double esous, G4double *dens, G4double bshell, G4double bsin, G4double bkin, G4double *temp, G4int optshp, G4int optcol, G4double defbet, G4double *ecor, G4double jprf, G4int ifis,G4double *qr)
{
  //   1498	C                                                                       
  //   1499	C     INPUT:                                                            
  //   1500	C             A,EE,ESOUS,OPTSHP,BS,BK,BSHELL,DEFBET                     
  //   1501	C                                                                       
  //   1502	C     LEVEL DENSITY PARAMETERS                                          
  //   1503	C     COMMON /ALD/    AV,AS,AK,OPTAFAN                                  
  //   1504	C     AV,AS,AK - VOLUME,SURFACE,CURVATURE DEPENDENCE OF THE             
  //   1505	C                LEVEL DENSITY PARAMETER                                
  //   1506	C     OPTAFAN - 0/1  AF/AN >=1 OR AF/AN ==1                             
  //   1507	C               RECOMMENDED IS OPTAFAN = 0                              
  //   1508	C---------------------------------------------------------------------  
  //   1509	C     OUTPUT: DENS,TEMP                                                 
  //   1510	C                                                                       
  //   1511	C   ____________________________________________________________________
  //   1512	C  /                                                                    
  //   1513	C  /  PROCEDURE FOR CALCULATING THE STATE DENSITY OF A COMPOUND NUCLEUS 
  //   1514	C  /____________________________________________________________________
  //   1515	C                                                                       
  //   1516	      INTEGER AFP,IZ,OPTSHP,OPTCOL,J,OPTAFAN                            
  //   1517	      REAL*8 A,EE,ESOUS,DENS,E,Y0,Y1,Y2,Y01,Y11,Y21,PA,BS,BK,TEMP       
  //   1518	C=====INSERTED BY KUDYAEV===============================================
  //   1519	      COMMON /ALD/ AV,AS,AK,OPTAFAN                                     
  //   1520	      REAL*8 ECR,ER,DELTAU,Z,DELTPP,PARA,PARZ,FE,HE,ECOR,ECOR1,Pi6      
  //   1521	      REAL*8 BSHELL,DELTA0,AV,AK,AS,PONNIV,PONFE,DEFBET,QR,SIG,FP       
  //   1522	C=======================================================================
  //   1523	C                                                                       
  //   1524	C                                                                       
  //   1525	C-----------------------------------------------------------------------
  //   1526	C     A                 MASS NUMBER OF THE DAUGHTER NUCLEUS             
  //   1527	C     EE                EXCITATION ENERGY OF THE MOTHER NUCLEUS         
  //   1528	C     ESOUS             SEPARATION ENERGY PLUS EFFECTIVE COULOMB BARRIER
  //   1529	C     DENS              STATE DENSITY OF DAUGHTER NUCLEUS AT EE-ESOUS-EC
  //   1530	C     BSHELL            SHELL CORRECTION                                
  //   1531	C     TEMP              NUCLEAR TEMPERATURE                             
  //   1532	C     E        LOCAL    EXCITATION ENERGY OF THE DAUGHTER NUCLEUS       
  //   1533	C     E1       LOCAL    HELP VARIABLE                                   
  //   1534	C     Y0,Y1,Y2,Y01,Y11,Y21                                              
  //   1535	C              LOCAL    HELP VARIABLES                                  
  //   1536	C     PA       LOCAL    STATE-DENSITY PARAMETER                         
  //   1537	C     EC                KINETIC ENERGY OF EMITTED PARTICLE WITHOUT      
  //   1538	C                        COULOMB REPULSION                              
  //   1539	C     IDEN              FAKTOR FOR SUBSTRACTING KINETIC ENERGY IDEN*TEMP
  //   1540	C     DELTA0            PAIRING GAP 12 FOR GROUND STATE                 
  //   1541	C                       14 FOR SADDLE POINT                             
  //   1542	C     EITERA            HELP VARIABLE FOR TEMPERATURE ITERATION         
  //   1543	C-----------------------------------------------------------------------
  //   1544	C                                                                       
  //   1545	C                                                                       
  G4double delta0 = 0.0;
  G4double deltau = 0.0;
  G4double deltpp = 0.0;
  G4double e = 0.0;
  G4double e0 = 0.0;
  G4double ecor1 = 0.0;
  G4double ecr = 10.0;
  G4double fe = 0.0;
  G4double he = 0.0;
  G4double pa = 0.0;
  G4double para = 0.0;
  G4double parz = 0.0;
  G4double ponfe = 0.0;
  G4double ponniv = 0.0;
  G4double fqr = 1.0;
  G4double y01 = 0.0;
  G4double y11 = 0.0;
  G4double y2 = 0.0;
  G4double y21 = 0.0;
  G4double y1 = 0.0;
  G4double y0 = 0.0;
  G4double fnorm=0.0;
  G4double fp_per=0.;
  G4double fp_par=0.;
  G4double sig_per=0.;
  G4double sig_par=0.;
  G4double sigma2;
  G4double jfact=1.;
  G4double erot=0.;
  G4double fdens=0.;
  G4double fecor=0.;
  G4double BSHELLCT=0.;
  G4double gamma=0.;
  G4double ftemp=0.0;
  G4double tempct=0.0;
  G4double densfm = 0.0;
  G4double densct = 0.0;
  G4double ein=0.;
  G4double elim;
  G4double tfm;
  G4double bs=bsin;
  G4double bk=bkin;
  G4int IPARITE;
  G4int IOPTCT=fiss->optct;
//
  G4double pi6 = std::pow(3.1415926535,2) / 6.0;
  G4double pi = 3.1415926535;
//
  G4int afp=idnint(a);
  G4int iz=idnint(z);
  G4int in=afp-iz;
//
  if(ifis!=1){
     BSHELLCT = ecld->ecgnz[in][iz];
  }else{
     BSHELLCT = 0.0;
  }
  if(afp<=20) BSHELLCT = 0.0;
  //
   parite(a,&para);
   if (para < 0.0){
// Odd A
	IPARITE=1;
   }else{   
// Even A                                                      
      parite(z,&parz);
      if(parz > 0.0){
// Even Z, even N
        IPARITE=2;
      }else{
// Odd Z, odd N 
        IPARITE=0;
      }
   }
//
   ein = ee - esous;
//
   if(ein>1.e30){
      fdens = 0.0;
      ftemp = 0.5;
      goto densniv100;
   }
//
   e = ee - esous;
//
   if(e<0.0&&ifis!=1){  // TUNNELING
        fdens = 0.0;
        densfm = 0.0;
        densct = 0.0;
        if(ald->optafan == 1) {
         pa = (ald->av)*a + (ald->as)*std::pow(a,(2.e0/3.e0)) + (ald->ak)*std::pow(a,(1.e0/3.e0));
        }else {
         pa = (ald->av)*a + (ald->as)*bsin*std::pow(a,(2.e0/3.e0)) + (ald->ak)*bkin*std::pow(a,(1.e0/3.e0));
        }
        gamma = 2.5 * pa * std::pow(a,-4.0/3.0);
        fecor=0.0;
        goto densniv100;
   }
//
   if(ifis==0&&bs!=1.0){
// - With increasing excitation energy system in getting less and less deformed:
     G4double ponq = (e-100.0)/5.0;
     if(ponq>700.0) ponq = 700.0;
       bs = 1.0/(1.0+std::exp(-ponq)) + 1.0/(1.0+std::exp(ponq)) * bsin;
       bk = 1.0/(1.0+std::exp(-ponq)) + 1.0/(1.0+std::exp(ponq)) * bkin;
   }
//
  // level density parameter                                               
  if(ald->optafan == 1) {
    pa = (ald->av)*a + (ald->as)*std::pow(a,(2.e0/3.e0)) + (ald->ak)*std::pow(a,(1.e0/3.e0));
  }
  else {
    pa = (ald->av)*a + (ald->as)*bs*std::pow(a,(2.e0/3.e0)) + (ald->ak)*bk*std::pow(a,(1.e0/3.e0));
  }
//
  gamma = 2.5 * pa * std::pow(a,-4.0/3.0);
//
// AK - 2009 - trial, in order to have transition to constant-temperature approach
// Idea - at the phase transition superfluid-normal fluid, TCT = TEMP, and this
// determines critical energy for pairing.
  if(a>0.0){
   ecr = pa*17.60/(std::pow(a,0.699) * std::sqrt(1.0+gamma*BSHELLCT))*17.60/(std::pow(a,0.699) * std::sqrt(1.0+gamma*BSHELLCT));
  }

  // pairing corrections                                                   
  if (ifis == 1) {
    delta0 = 14;
  }
  else {
    delta0 = 12;
  }

  // shell corrections                                                     
  if (optshp > 0) {
    deltau = bshell;
    if (optshp == 2) {
      deltau = 0.0;
    }
    if (optshp >= 2) {
      // pairing energy shift with condensation energy a.r.j. 10.03.97        
    //deltpp = -0.25e0* (delta0/pow(sqrt(a),2)) * pa /pi6 + 2.e0*delta0/sqrt(a);
      deltpp = -0.25e0* std::pow((delta0/std::sqrt(a)),2) * pa /pi6 + 22.34e0*std::pow(a,-0.464)-0.235;
      // Odd A
      if (IPARITE == 1) {
	//e = e - delta0/sqrt(a);
        e=e-(0.285+11.17*std::pow(a,-0.464)-0.390-0.00058*a);//-30./a;//FIXME
      }
      // Even Z, even N
      if(IPARITE==2){ 
           e=e-(22.34*std::pow(a,-0.464)-0.235);//-30./a;//FIXME
      }
      // Odd Z, odd N
      if(IPARITE==0){
        if(in==iz){
         //  e = e;
        }else{
         //  e = e-30./a;
        }
      }
    } else {                                                          
      deltpp = 0.0;
    }
  }else {
    deltau = 0.0;
    deltpp = 0.0;
  }

  if(e < 0.0){
    e = 0.0;
    ftemp = 0.5;
  }

  // washing out is made stronger                          
  ponfe = -2.5*pa*e*std::pow(a,(-4.0/3.0));

  if (ponfe < -700.0)  {
    ponfe = -700.0;
  }
  fe = 1.0 - std::exp(ponfe);
  if (e < ecr) {
    // priv. comm. k.-h. schmidt                                         
    he = 1.0 - std::pow((1.0 - e/ecr),2);
  }
  else {
    he = 1.0;
  }
  // Excitation energy corrected for pairing and shell effects             
  // washing out with excitation energy is included.                        
  fecor = e + deltau*fe + deltpp*he;
  if (fecor <= 0.1) {
    fecor = 0.1;
  }                                          
  // iterative procedure according to grossjean and feldmeier              
  // to avoid the singularity e = 0                                        
  if (ee < 5.0) {
    y1 = std::sqrt(pa*fecor);
    for(G4int j = 0; j < 5; j++) {
      y2 = pa*fecor*(1.e0-std::exp(-y1));
      y1 = std::sqrt(y2);
    }
    y0 = pa/y1;
    ftemp=1.0/y0;
    fdens = std::exp(y0*fecor)/ (std::pow((std::pow(fecor,3)*y0),0.5)*std::pow((1.0-0.5*y0*fecor*std::exp(-y1)),0.5))* std::exp(y1)*(1.0-std::exp(-y1))*0.1477045;
    if (fecor < 1.0) {
      ecor1=1.0;
      y11 = std::sqrt(pa*ecor1);
      for(G4int j = 0; j < 7; j++) {
	y21 = pa*ecor1*(1.0-std::exp(-y11));
	y11 = std::sqrt(y21);
      }

      y01 = pa/y11;
      fdens = fdens*std::pow((y01/y0),1.5);
      ftemp = ftemp*std::pow((y01/y0),1.5);
    }
  }
  else {
    ponniv = 2.0*std::sqrt(pa*fecor);
    if (ponniv > 700.0) {
      ponniv = 700.0;
    }
    // fermi gas state density                                               
    fdens = 0.1477045 * std::exp(ponniv)/(std::pow(pa,0.25)*std::pow(fecor,1.25));
    ftemp = std::sqrt(fecor/pa);
  }
//
  densfm = fdens;
  tfm = ftemp;
//
  if(IOPTCT==0) goto densniv100;
  tempct = 17.60/( std::pow(a,0.699) * std::sqrt(1.+gamma*BSHELLCT));
  //tempct = 1.0 / ( (0.0570 + 0.00193*BSHELLCT) * pow(a,0.6666667));  // from  PRC 80 (2009) 054310

// - CONSTANT-TEMPERATURE LEVEL DENSITY PARAMETER (ONLY AT LOW ENERGIES)
  if(e<30.){
   if(a>0.0){
     if(optshp>=2){
// Parametrization of CT model by Ignatyuk; note that E0 is shifted to correspond
// to pairing shift in Fermi-gas model (there, energy is shifted taking odd-odd nuclei
//  as bassis)
// e-o, o-e
            if (IPARITE == 1) { e0 = 0.285+11.17*std::pow(a,-0.464) - 0.390-0.00058*a;}
// e-e
            if (IPARITE == 2) { e0 = 22.34*std::pow(a,-0.464)-0.235;}
// o-o
            if (IPARITE == 0){ e0 = 0.0;}

        ponniv = (ein-e0)/tempct;
        if(ifis!=1) ponniv = max(0.0,(ein-e0)/tempct);
        if(ponniv>700.0){ ponniv = 700.0;}
        densct = std::exp(ponniv)/tempct*std::exp(0.079*BSHELLCT/tempct);

        elim = ein;

        if(elim>=ecr&&densfm<=densct){
          fdens = densfm;
        //  IREGCT = 0;
        }else{
          fdens = densct;
         // IREGCT = 1;
//         ecor = min(ein-e0,0.10);
        }
        if(elim>=ecr&&tfm>=tempct){
          ftemp = tfm;
        }else{
          ftemp = tempct;
        }
     }else{
// Case of no pairing considered
//        ETEST = PA * TEMPCT**2
        ponniv = (ein)/tempct;
        if(ponniv>700.0){ ponniv = 700.0;}
        densct = std::exp(ponniv)/tempct;

        if(ein>=ecr && densfm<=densct){
          fdens = densfm;
          ftemp = tfm;
        //  IREGCT = 0;
        }else{
          fdens = densct;
          ftemp = tempct;
//          ECOR = DMIN1(EIN,0.1D0)
        }

        if(ein>=ecr && tfm>=tempct){
          ftemp = tfm;
        }else{
          ftemp = tempct;
        }
     }
   }
  }


 densniv100:

 if(fdens==0.0){
    if(a>0.0){
// Parametrization of CT model by Ignatyuk done for masses > 20
         ftemp = 17.60/( std::pow(a,0.699) * std::sqrt(1.0+gamma*BSHELLCT));
       //  ftemp = 1.0 / ( (0.0570 + 0.00193*BSHELLCT) * pow(a,0.6666667));  // from  PRC 80 (2009) 054310
    }else{
         ftemp = 0.5;
    }
 }
//
// spin cutoff parameter
/*
C PERPENDICULAR AND PARALLEL MOMENT OF INERTIA
c fnorm = R0*M0/hbar**2 = 1.16fm*931.49MeV/c**2 /(6.582122e-22 MeVs)**2 and is
c in units 1/MeV
*/
 fnorm = std::pow(1.16,2)*931.49*1.e-2/(9.0* std::pow(6.582122,2));

 if(ifis==0 || ifis==2){
/*
C GROUND STATE:
C FP_PER ~ 1+0.5*alpha2, FP_PAR ~ 1-alpha2 (Hasse & Myers, Geom. relat. macr. nucl. phys.)
C alpha2 = sqrt(5/(4*pi))*beta2
*/
    fp_per = 0.4*std::pow(a,5.0/3.0)*fnorm*(1.0+0.50*defbet*std::sqrt(5.0/(4.0*pi)));
    fp_par = 0.40*std::pow(a,5.0/3.0)*fnorm*(1.0-defbet*std::sqrt(5.0/(4.0*pi)));

 }else{
  if(ifis==1){
/*
C SADDLE POINT
C See Hasse&Myer, p. 100
C Perpendicular moment of inertia
*/
    fp_per = 2.0/5.0*std::pow(a,5.0/3.0)*fnorm*(1.0+7.0/6.0*defbet*(1.0+1396.0/255.0*defbet));
// Parallel moment of inertia
    fp_par = 2.0/5.0*std::pow(a,5.0/3.0)*fnorm*(1.0-7.0/3.0*defbet*(1.0-389.0/255.0*defbet));
  }else{
   if(ifis==20){
// IMF - two fragments in contact; it is asumed that both are spherical.
// See Hasse&Myers, p.106
// Here, DEFBET = R1/R2, where R1 and R2 are radii of IMF and its partner
// Perpendicular moment of inertia
   fp_per = 0.4*std::pow(a,5.0/3.0)*fnorm*3.50*(1.0 + std::pow(defbet,5.))/std::pow(1.0 + defbet*defbet*defbet,5.0/3.0);
   fp_par = 0.4*std::pow(a,5.0/3.0)*fnorm*(1.0 + std::pow(defbet,5.0))/std::pow(1.0 + defbet*defbet*defbet,5.0/3.0);
   }
  }
 }
  if(fp_par<0.0)fp_par=0.0;
  if(fp_per<0.0)fp_per=0.0;
//
  sig_per = std::sqrt(fp_per * ftemp);
  sig_par = std::sqrt(fp_par * ftemp);
//
  sigma2 = sig_per*sig_per + sig_par*sig_par;
  jfact = (2.*jprf+1.)*std::exp(-1.*jprf*(jprf+1.0)/(2.0*sigma2))/(std::sqrt(8.0*3.1415)*std::pow(sigma2,1.5));
  erot = jprf*jprf/(2.0*std::sqrt(fp_par*fp_par+fp_per*fp_per));
//
  // collective enhancement                                            
  if (optcol == 1) {
    qrot(z,a,defbet,sig_per,fecor-erot,&fqr);
  }
  else {
    fqr   = 1.0;
  }
//
  fdens = fdens * fqr *jfact;
//
  if(fdens<1e-300)fdens=0.0;
//
  *dens =fdens;
  *ecor=fecor;
  *temp=ftemp;
  *qr=fqr;
}

void G4Abla::qrot(G4double z, G4double a, G4double bet, G4double sig, G4double u, G4double *qr)
{
/*
C QROT INCLUDING DAMPING
C
C INPUT: Z,A,DEFBET,SIG,U
C
C OUTPUT: QR - COLLECTIVE ENHANCEMENT FACTOR
C
C SEE  JUNGHANS ET AL., NUCL. PHYS. A 629 (1998) 635
C
C
C   FR(U)    EXPONENTIAL FUNCTION TO DEFINE DAMPING
C   UCR      CRITICAL ENERGY FOR DAMPING
C   DCR      WIDTH OF DAMPING
C   DEFBET   BETA-DEFORMATION !
C   SIG      PERPENDICULAR SPIN CUTOFF FACTOR
C     U      ENERGY
C    QR      COEFFICIENT OF COLLECTIVE ENHANCEMENT
C     A      MASS NUMBER
C     Z      CHARGE NUMBER
C
*/
// JLRS: July 2016: new values for the collective parameters
//

  G4double ucr = fiss->ucr; // Critical energy for damping.
  G4double dcr = fiss->dcr; // Width of damping.
  G4double ponq = 0.0, dn = 0.0, n = 0.0, dz = 0.0;
  G4int distn,distz,ndist, zdist;
  G4int nmn[8]= {2, 8, 14, 20, 28, 50, 82, 126};
  G4int nmz[8]= {2, 8, 14, 20, 28, 50, 82, 126};
//
  sig = sig*sig;
//
 if(std::abs(bet)<=0.15){
  goto qrot10;
 }else{
  goto qrot11;
 }
//
 qrot10:
  n = a - z;
  distn = 10000000;
  distz = 10000000;

  for(G4int i =0;i<8;i++){
   ndist = std::fabs(idnint(n) - nmn[i]);
   if(ndist < distn) distn = ndist;
   zdist = std::fabs(idnint(z) - nmz[i]);
   if(zdist < distz) distz = zdist;
  }

  dz = G4float(distz);
  dn = G4float(distn);

  bet = 0.022 + 0.003*dn + 0.002*dz;

  sig = 75.0*std::pow(bet,2.) * sig;

// NO VIBRATIONAL ENHANCEMENT
 qrot11:   
  ponq = (u - ucr)/dcr;

  if (ponq > 700.0) {
    ponq = 700.0;
  }
  if (sig < 1.0) {
    sig = 1.0;
  }
  (*qr) = 1.0/(1.0 + std::exp(ponq)) * (sig - 1.0) + 1.0;

  if ((*qr) < 1.0) {
    (*qr) = 1.0;
  }

  return;
}

void G4Abla::lpoly(G4double x, G4int n, G4double pl[])
{
  // THIS SUBROUTINE CALCULATES THE ORDINARY LEGENDRE POLYNOMIALS OF   
  // ORDER 0 TO N-1 OF ARGUMENT X AND STORES THEM IN THE VECTOR PL.    
  // THEY ARE CALCULATED BY RECURSION RELATION FROM THE FIRST TWO      
  // POLYNOMIALS.                                                      
  // WRITTEN BY A.J.SIERK  LANL  T-9  FEBRUARY, 1984            
  // NOTE: PL AND X MUST BE DOUBLE PRECISION ON 32-BIT COMPUTERS!      

  pl[0] = 1.0;
  pl[1] = x;

  for(G4int i = 2; i < n; i++) {
    pl[i] = ((2*G4double(i+1) - 3.0)*x*pl[i-1] - (G4double(i+1) - 2.0)*pl[i-2])/(G4double(i+1)-1.0);
  }
}

G4double G4Abla::eflmac(G4int ia, G4int iz, G4int flag, G4int optshp)
{
  // CHANGED TO CALCULATE TOTAL BINDING ENERGY INSTEAD OF MASS EXCESS.     
  // SWITCH FOR PAIRING INCLUDED AS WELL.                                  
  // BINDING = EFLMAC(IA,IZ,0,OPTSHP)                                      
  // FORTRAN TRANSCRIPT OF /U/GREWE/LANG/EEX/FRLDM.C                       
  // A.J. 15.07.96                                                         

  // this function will calculate the liquid-drop nuclear mass for spheri
  // configuration according to the preprint NUCLEAR GROUND-STATE        
  // MASSES and DEFORMATIONS by P. M"oller et al. from August 16, 1993 p.
  // All constants are taken from this publication for consistency.      

  // Parameters:                                                         
  // a:    nuclear mass number                                         
  // z:    nuclear charge                                              
  // flag:     0       - return mass excess                            
  //       otherwise   - return pairing (= -1/2 dpn + 1/2 (Dp + Dn))   

  G4double eflmacResult = 0.0;

  if(ia==0)return eflmacResult;

  G4int in = 0;
  G4double z = 0.0, n = 0.0, a = 0.0, av = 0.0, as = 0.0;
  G4double a0 = 0.0, c1 = 0.0, c4 = 0.0, b1 = 0.0, b3 = 0.0;
  G4double ff = 0.0, ca = 0.0, w = 0.0, efl = 0.0; 
  G4double r0 = 0.0, kf = 0.0, ks = 0.0;
  G4double kv = 0.0, rp = 0.0, ay = 0.0, aden = 0.0, x0 = 0.0, y0 = 0.0;
  G4double esq = 0.0, ael = 0.0, i = 0.0, e0 = 0.0;
  G4double pi = 3.141592653589793238e0;

  // fundamental constants
  // electronic charge squared
  esq = 1.4399764;

  // constants from considerations other than nucl. masses
  // electronic binding
  ael = 1.433e-5;

  // proton rms radius
  rp  = 0.8;

  // nuclear radius constant
  r0  = 1.16;

  // range of yukawa-plus-expon. potential
  ay  = 0.68;

  // range of yukawa function used to generate                          
  // nuclear charge distribution
  aden= 0.70;

  // wigner constant
  w   = 30.0;

  // adjusted parameters
  // volume energy
  av  = 16.00126;

  // volume asymmetry
  kv  =  1.92240;

  // surface energy
  as  = 21.18466;

  // surface asymmetry
  ks  =  2.345;
  // a^0 constant
  a0  =  2.615;

  // charge asymmetry
  ca  =  0.10289;

  z   = G4double(iz);
  a   = G4double(ia);
  in  = ia - iz;                                                       
  n   = G4double(in);

  if(flag==1){goto eflmac311;}

  if(iz<13&&in<3){
     if(masses->mexpiop[in][iz]==1){
         return masses->bind[in][iz];
     }
  }

  eflmac311:
  
  c1  = 3.0/5.0*esq/r0;
  c4  = 5.0/4.0*std::pow((3.0/(2.0*pi)),(2.0/3.0)) * c1;
  kf  = std::pow((9.0*pi*z/(4.0*a)),(1.0/3.0))/r0;
  
  ff = -1.0/8.0*rp*rp*esq/std::pow(r0,3) * (145.0/48.0 - 327.0/2880.0*std::pow(kf,2) * std::pow(rp,2) + 1527.0/1209600.0*std::pow(kf,4) * std::pow(rp,4));
  i   = (n-z)/a;

  x0  = r0 * std::pow(a,(1.0/3.0)) / ay;
  y0  = r0 * std::pow(a,(1.0/3.0)) / aden;

  b1  = 1.0 - 3.0/(std::pow(x0,2)) + (1.0 + x0) * (2.0 + 3.0/x0 + 3.0/std::pow(x0,2)) * std::exp(-2.0*x0);

  b3  = 1.0 - 5.0/std::pow(y0,2) * (1.0 - 15.0/(8.0*y0) + 21.0/(8.0 * std::pow(y0,3))
			       - 3.0/4.0 * (1.0 + 9.0/(2.0*y0) + 7.0/std::pow(y0,2)
					    + 7.0/(2.0 * std::pow(y0,3))) * std::exp(-2.0*y0));

  // now calculation of total binding energy a.j. 16.7.96                   

  efl = -1.0 * av*(1.0 - kv*i*i)*a + as*(1.0 - ks*i*i)*b1 * std::pow(a,(2.0/3.0)) + a0
    + c1*z*z*b3/std::pow(a,(1.0/3.0)) - c4*std::pow(z,(4.0/3.0))/std::pow(a,(1.e0/3.e0))
    + ff*std::pow(z,2)/a -ca*(n-z) - ael * std::pow(z,(2.39e0));

  efl = efl + w*std::abs(i);

  // pairing is made optional                                              
  if (optshp >= 2) {
    // average pairing
    if (in==iz && (mod(in,2) == 1) && (mod(iz,2) == 1) && in>0.) {
      efl = efl + w/a;
    }

// AK 2008 - Parametrization of CT model by Ignatyuk;
// The following part has been introduced  in order to have correspondance
// between pairing in masses and level densities;
// AK 2010  note that E0 is shifted to correspond to pairing shift in
// Fermi-gas model (there, energy is shifted taking odd-odd nuclei
// as bassis)

    G4double para=0.;
    parite(a,&para);

    if(para<0.0){
// e-o, o-e
      e0 =  0.285+11.17*std::pow(a,-0.464) -0.390-0.00058*(a);
    }else{
        G4double parz=0.;
        parite(z,&parz);
        if (parz>0.0){
// e-e
         e0 = 22.34*std::pow(a,-0.464)-0.235;
        }else{
// o-o
         e0 = 0.0;
        }
    }
    efl = efl - e0;
  // end if for pairing term                                  
  }

  eflmacResult = efl;

  return eflmacResult;
}

void G4Abla::appariem(G4double a, G4double z, G4double *del)
{
  // CALCUL DE LA CORRECTION, DUE A L'APPARIEMENT, DE L'ENERGIE DE     
  // LIAISON D'UN NOYAU                                                
  // PROCEDURE FOR CALCULATING THE PAIRING CORRECTION TO THE BINDING   
  // ENERGY OF A SPECIFIC NUCLEUS                                      

  G4double para = 0.0, parz = 0.0;
  // A                 MASS NUMBER                                     
  // Z                 NUCLEAR CHARGE                                  
  // PARA              HELP VARIABLE FOR PARITY OF A                   
  // PARZ              HELP VARIABLE FOR PARITY OF Z                   
  // DEL               PAIRING CORRECTION                              

  parite(a, &para);

  if (para < 0.0) {
    (*del) = 0.0;
  }
  else {
    parite(z, &parz);
    if (parz > 0.0) {
      (*del) = -12.0/std::sqrt(a);
    }
    else {
      (*del) = 12.0/std::sqrt(a);
    }
  }
}

void G4Abla::parite(G4double n, G4double *par)
{
  // CALCUL DE LA PARITE DU NOMBRE N                                   
  //
  // PROCEDURE FOR CALCULATING THE PARITY OF THE NUMBER N.             
  // RETURNS -1 IF N IS ODD AND +1 IF N IS EVEN                        

  G4double n1 = 0.0, n2 = 0.0, n3 = 0.0;

  // N                 NUMBER TO BE TESTED                             
  // N1,N2             HELP VARIABLES                                  
  // PAR               HELP VARIABLE FOR PARITY OF N                   

  n3 = G4double(idnint(n));
  n1 = n3/2.0;
  n2 = n1 - dint(n1);

  if (n2 > 0.0) {
    (*par) = -1.0;
  }
  else {
    (*par) = 1.0;
  }
}

G4double G4Abla::tau(G4double bet, G4double homega, G4double ef, G4double t)
{
  // INPUT : BET, HOMEGA, EF, T                                          
  // OUTPUT: TAU - RISE TIME IN WHICH THE FISSION WIDTH HAS REACHED      
  //               90 PERCENT OF ITS FINAL VALUE                               
  // 
  // BETA   - NUCLEAR VISCOSITY                                          
  // HOMEGA - CURVATURE OF POTENTIAL                                     
  // EF     - FISSION BARRIER                                            
  // T      - NUCLEAR TEMPERATURE                                        

  G4double tauResult = 0.0;

  G4double tlim = 8.e0 * ef;
  if (t > tlim) {
    t = tlim;
  }
  //
  if (bet/(std::sqrt(2.0)*10.0*(homega/6.582122)) <= 1.0) {
    tauResult = std::log(10.0*ef/t)/(bet*1.0e21);
  }
  else {
    tauResult = std::log(10.0*ef/t)/ (2.0*std::pow((10.0*homega/6.582122),2))*(bet*1.0e-21);
  } //end if                                                            

  return tauResult;
}

G4double G4Abla::cram(G4double bet, G4double homega)
{
  // INPUT : BET, HOMEGA  NUCLEAR VISCOSITY + CURVATURE OF POTENTIAL      
  // OUTPUT: KRAMERS FAKTOR  - REDUCTION OF THE FISSION PROBABILITY       
  //                           INDEPENDENT OF EXCITATION ENERGY                             

  G4double rel = bet/(20.0*homega/6.582122);
  G4double cramResult = std::sqrt(1.0 + std::pow(rel,2)) - rel;
  // limitation introduced   6.1.2000  by  khs

  if (cramResult > 1.0) {
    cramResult = 1.0;
  }

  return cramResult;
}

G4double G4Abla::bipol(G4int iflag, G4double y)
{
  // CALCULATION OF THE SURFACE BS OR CURVATURE BK OF A NUCLEUS        
  // RELATIVE TO THE SPHERICAL CONFIGURATION                           
  // BASED ON  MYERS, DROPLET MODEL FOR ARBITRARY SHAPES               

  // INPUT: IFLAG - 0/1 BK/BS CALCULATION                              
  //         Y    - (1 - X) COMPLEMENT OF THE FISSILITY                

  // LINEAR INTERPOLATION OF BS BK TABLE                               

  G4int i = 0;

  G4double bipolResult = 0.0;

  const G4int bsbkSize = 54;

  G4double bk[bsbkSize] = {0.0, 1.00000,1.00087,1.00352,1.00799,1.01433,1.02265,1.03306,
			   1.04576,1.06099,1.07910,1.10056,1.12603,1.15651,1.19348,
			   1.23915,1.29590,1.35951,1.41013,1.44103,1.46026,1.47339,
			   1.48308,1.49068,1.49692,1.50226,1.50694,1.51114,1.51502,
			   1.51864,1.52208,1.52539,1.52861,1.53177,1.53490,1.53803,
			   1.54117,1.54473,1.54762,1.55096,1.55440,1.55798,1.56173,
			   1.56567,1.56980,1.57413,1.57860,1.58301,1.58688,1.58688,
			   1.58688,1.58740,1.58740, 0.0}; //Zeroes at bk[0], and at the end added by PK

  G4double bs[bsbkSize] = {0.0, 1.00000,1.00086,1.00338,1.00750,1.01319,
			   1.02044,1.02927,1.03974,
			   1.05195,1.06604,1.08224,1.10085,1.12229,1.14717,1.17623,1.20963,
			   1.24296,1.26532,1.27619,1.28126,1.28362,1.28458,1.28477,1.28450,
			   1.28394,1.28320,1.28235,1.28141,1.28042,1.27941,1.27837,1.27732,
			   1.27627,1.27522,1.27418,1.27314,1.27210,1.27108,1.27006,1.26906,
			   1.26806,1.26707,1.26610,1.26514,1.26418,1.26325,1.26233,1.26147,
			   1.26147,1.26147,1.25992,1.25992, 0.0};

  i = idint(y/(2.0e-02)) + 1;
    
  if((i + 1) >= bsbkSize) {
    if(verboseLevel > 2) {
      // G4cout <<"G4Abla error: index " << i + 1 << " is greater than array size permits." << G4endl;
    }
    bipolResult = 0.0;
  }
  else {
    if (iflag == 1) {
      bipolResult = bs[i] + (bs[i+1] - bs[i])/2.0e-02 * (y - 2.0e-02*(i - 1));
    }
    else {
      bipolResult = bk[i] + (bk[i+1] - bk[i])/2.0e-02 * (y - 2.0e-02*(i - 1));
    }
  }
  
  return bipolResult;
}

void G4Abla::fomega_sp(G4double AF,G4double Y,G4double *MFCD,G4double *sOMEGA,G4double *sHOMEGA)
{
/*
c  Y                 1 - Fissility
c  OMEGA             Frequency at the ground state, in units 1.e-21 s
*/
  G4double OMEGA,HOMEGA,ES0,MR02;

      ES0 = 20.760*std::pow(AF,2.0/3.0);
// In units 1.e-42 MeVs**2; r0 = 1.175e-15 m, u=931.49MeV/c**2=103.4MeV*s**2/m**2
// divided by 1.e-4 to go from 1.e-46 to 1.e-42
      MR02 = std::pow(AF,5.0/3.0)*1.0340*0.010*1.175*1.175;
// Determination of the inertia of the fission collective degree of freedom
      (*MFCD) = MR02 * 3.0/10.0*(1.0+3.0*Y);
// Omega at saddle
      OMEGA = std::sqrt(ES0/MR02)*std::sqrt(8.0/3.0*Y*(1.0+304.0*Y/255.0));
//
      HOMEGA = 6.58122*OMEGA/10.0;
//
   (*sOMEGA)=OMEGA;
   (*sHOMEGA)=HOMEGA;
//
  return;
}


void G4Abla::fomega_gs(G4double AF,G4double ZF,G4double *K1,G4double *sOMEGA,G4double *sHOMEGA)
{
/*
c  Y                 1 - Fissility
c  OMEGA             Frequency at the ground state, in units 1.e-21 s
*/
  G4double OMEGA,HOMEGA,MR02,MINERT,C,fk1;
//
      MR02 = std::pow(AF,5.0/3.0)*1.0340*0.01*1.175*1.175;
      MINERT = 3.*MR02/10.0;
      C = 17.9439*(1.-1.7826*std::pow((AF-2.0*ZF)/AF,2));
      fk1 = 0.4*C*std::pow(AF,2.0/3.0)-0.1464*std::pow(ZF,2)/std::pow(AF,1./3.);
      OMEGA = std::sqrt(fk1/MINERT);
      HOMEGA = 6.58122*OMEGA/10.0;
//
   (*K1)=fk1;
   (*sOMEGA)=OMEGA;
   (*sHOMEGA)=HOMEGA;
//
  return;
}

void G4Abla::barrs(G4int Z1,G4int A1,G4int Z2,G4int A2,G4double *sBARR,G4double *sOMEGA)
{/*
C AK 2004 - Barriers for LCP and IMF are calculated now according to the
C           Bass model (Nucl. Phys. A (1974))
C KHS 2007 - To speed up, barriers are read from tabels; in case thermal
C            expansion is considered, barriers are calculated.
C INPUT:
C EA    - Excitation energy per nucleon
C Z11, A11 - Charge and mass of daughter nucleus
C Z22, A22 - Charge and mass of LCP or IMF
C
C OUTPUT:
C BARR - Barrier
C OMEGA - Curvature of the potential
C
C BASS MODEL NPA 1974 - used only if expansion is considered (OPTEXP=1)
C                        or one wants this model explicitly (OPTBAR=1)
C October 2011 - AK - new parametrization of the barrier and its position,
C                    see W.W. Qu et al., NPA 868 (2011) 1; this is now
C                    default option (OPTBAR=0)
c
c November 2016 - JLRS - Added this function from abla07v4
c
*/
     G4double BARR, OMEGA, RMAX;
     RMAX = 1.1 * (ecld->rms[A1-Z1][Z1]+ecld->rms[A2-Z2][Z2]) + 2.8;
     BARR = 1.345 * Z1 * Z2 / RMAX;
//C Omega according to Avishai:
     OMEGA = 4.5 / 197.3287;

    // if(Z1<60){
    //  if(Z2==1 && A2==2) BARR = BARR * 1.1;
    //  if(Z2==1 && A2==3) BARR = BARR * 1.1;
   //   if(Z2==2 && A2==3) BARR = BARR * 1.3;
    //  if(Z2==2 && A2==4) BARR = BARR * 1.1;
    // }

     (*sOMEGA)=OMEGA;
     (*sBARR)=BARR;
//
  return;
}

void G4Abla::barfit(G4int iz, G4int ia, G4int il, G4double *sbfis, G4double *segs, G4double *selmax)
{
  //   2223	C     VERSION FOR 32BIT COMPUTER                                        
  //   2224	C     THIS SUBROUTINE RETURNS THE BARRIER HEIGHT BFIS, THE              
  //   2225	C     GROUND-STATE ENERGY SEGS, IN MEV, AND THE ANGULAR MOMENTUM        
  //   2226	C     AT WHICH THE FISSION BARRIER DISAPPEARS, LMAX, IN UNITS OF        
  //   2227	C     H-BAR, WHEN CALLED WITH INTEGER AGUMENTS IZ, THE ATOMIC           
  //   2228	C     NUMBER, IA, THE ATOMIC MASS NUMBER, AND IL, THE ANGULAR           
  //   2229	C     MOMENTUM IN UNITS OF H-BAR. (PLANCK'S CONSTANT DIVIDED BY         
  //   2230	C     2*PI).                                                            
  //   2231	C                                                                       
  //   2232	C        THE FISSION BARRIER FO IL = 0 IS CALCULATED FROM A 7TH         
  //   2233	C     ORDER FIT IN TWO VARIABLES TO 638 CALCULATED FISSION              
  //   2234	C     BARRIERS FOR Z VALUES FROM 20 TO 110. THESE 638 BARRIERS ARE      
  //   2235	C     FIT WITH AN RMS DEVIATION OF 0.10 MEV BY THIS 49-PARAMETER        
  //   2236	C     FUNCTION.                                                         
  //   2237	C     IF BARFIT IS CALLED WITH (IZ,IA) VALUES OUTSIDE THE RANGE OF      
  //   2238	C     THE BARRIER HEIGHT IS SET TO 0.0, AND A MESSAGE IS PRINTED        
  //   2239	C     ON THE DEFAULT OUTPUT FILE.                                       
  //   2240	C                                                                       
  //   2241	C        FOR IL VALUES NOT EQUAL TO ZERO, THE VALUES OF L AT WHICH      
  //   2242	C     THE BARRIER IS 80% AND 20% OF THE L=0 VALUE ARE RESPECTIVELY      
  //   2243	C     FIT TO 20-PARAMETER FUNCTIONS OF Z AND A, OVER A MORE             
  //   2244	C     RESTRICTED RANGE OF A VALUES, THAN IS THE CASE FOR L = 0.         
  //   2245	C     THE VALUE OF L WHERE THE BARRIER DISAPPEARS, LMAX IS FIT TO       
  //   2246	C     A 24-PARAMETER FUNCTION OF Z AND A, WITH THE SAME RANGE OF        
  //   2247	C     Z AND A VALUES AS L-80 AND L-20.                                  
  //   2248	C        ONCE AGAIN, IF AN (IZ,IA) PAIR IS OUTSIDE OF THE RANGE OF      
  //   2249	C     VALIDITY OF THE FIT, THE BARRIER VALUE IS SET TO 0.0 AND A        
  //   2250	C     MESSAGE IS PRINTED. THESE THREE VALUES (BFIS(L=0),L-80, AND       
  //   2251	C     L-20) AND THE CONSTRINTS OF BFIS = 0 AND D(BFIS)/DL = 0 AT        
  //   2252	C     L = LMAX AND L=0 LEAD TO A FIFTH-ORDER FIT TO BFIS(L) FOR         
  //   2253	C     L>L-20. THE FIRST THREE CONSTRAINTS LEAD TO A THIRD-ORDER FIT     
  //   2254	C     FOR THE REGION L < L-20.                                          
  //   2255	C                                                                       
  //   2256	C        THE GROUND STATE ENERGIES ARE CALCULATED FROM A                
  //   2257	C     120-PARAMETER FIT IN Z, A, AND L TO 214 GROUND-STATE ENERGIES     
  //   2258	C     FOR 36 DIFFERENT Z AND A VALUES.                                  
  //   2259	C     (THE RANGE OF Z AND A IS THE SAME AS FOR L-80, L-20, AND          
  //   2260	C     L-MAX)                                                            
  //   2261	C                                                                       
  //   2262	C        THE CALCULATED BARRIERS FROM WHICH THE FITS WERE MADE WERE     
  //   2263	C     CALCULATED IN 1983-1984 BY A. J. SIERK OF LOS ALAMOS              
  //   2264	C     NATIONAL LABORATORY GROUP T-9, USING YUKAWA-PLUS-EXPONENTIAL      
  //   2265	C     G4DOUBLE FOLDED NUCLEAR ENERGY, EXACT COULOMB DIFFUSENESS           
  //   2266	C     CORRECTIONS, AND DIFFUSE-MATTER MOMENTS OF INERTIA.               
  //   2267	C     THE PARAMETERS OF THE MODEL R-0 = 1.16 FM, AS 21.13 MEV,          
  //   2268	C     KAPPA-S = 2.3, A = 0.68 FM.                                       
  //   2269	C     THE DIFFUSENESS OF THE MATTER AND CHARGE DISTRIBUTIONS USED       
  //   2270	C     CORRESPONDS TO A SURFACE DIFFUSENESS PARAMETER (DEFINED BY        
  //   2271	C     MYERS) OF 0.99 FM. THE CALCULATED BARRIERS FOR L = 0 ARE          
  //   2272	C     ACCURATE TO A LITTLE LESS THAN 0.1 MEV; THE OUTPUT FROM           
  //   2273	C     THIS SUBROUTINE IS A LITTLE LESS ACCURATE. WORST ERRORS MAY BE    
  //   2274	C     AS LARGE AS 0.5 MEV; CHARACTERISTIC UNCERTAINY IS IN THE RANGE    
  //   2275	C     OF 0.1-0.2 MEV. THE RMS DEVIATION OF THE GROUND-STATE FIT         
  //   2276	C     FROM THE 214 INPUT VALUES IS 0.20 MEV. THE MAXIMUM ERROR          
  //   2277	C     OCCURS FOR LIGHT NUCLEI IN THE REGION WHERE THE GROUND STATE      
  //   2278	C     IS PROLATE, AND MAY BE GREATER THAN 1.0 MEV FOR VERY NEUTRON      
  //   2279	C     DEFICIENT NUCLEI, WITH L NEAR LMAX. FOR MOST NUCLEI LIKELY TO     
  //   2280	C     BE ENCOUNTERED IN REAL EXPERIMENTS, THE MAXIMUM ERROR IS          
  //   2281	C     CLOSER TO 0.5 MEV, AGAIN FOR LIGHT NUCLEI AND L NEAR LMAX.        
  //   2282	C                                                                       
  //   2283	C     WRITTEN BY A. J. SIERK, LANL T-9                                  
  //   2284	C     VERSION 1.0 FEBRUARY, 1984                                        
  //   2285	C                                                                       
  //   2286	C     THE FOLLOWING IS NECESSARY FOR 32-BIT MACHINES LIKE DEC VAX,      
  //   2287	C     IBM, ETC                                                          

  G4double pa[7],pz[7],pl[10];
  for(G4int init_i = 0; init_i < 7; init_i++) {
    pa[init_i] = 0.0; 
    pz[init_i] = 0.0; 
  }
  for(G4int init_i = 0; init_i < 10; init_i++) {
    pl[init_i] = 0.0;
  }

  G4double a = 0.0, z = 0.0, amin = 0.0, amax = 0.0, amin2 = 0.0;
  G4double amax2 = 0.0, aa = 0.0, zz = 0.0, bfis = 0.0;
  G4double bfis0 = 0.0, ell = 0.0, el = 0.0, egs = 0.0, el80 = 0.0, el20 = 0.0;
  G4double elmax = 0.0, sel80 = 0.0, sel20 = 0.0, x = 0.0, y = 0.0, q = 0.0, qa = 0.0, qb = 0.0;
  G4double aj = 0.0, ak = 0.0, a1 = 0.0, a2 = 0.0;

  G4int i = 0, j = 0, k = 0, m = 0;
  G4int l = 0;

  G4double emncof[4][5] = {{-9.01100e+2,-1.40818e+3, 2.77000e+3,-7.06695e+2, 8.89867e+2}, 
			   {1.35355e+4,-2.03847e+4, 1.09384e+4,-4.86297e+3,-6.18603e+2},
			   {-3.26367e+3, 1.62447e+3, 1.36856e+3, 1.31731e+3, 1.53372e+2},
			   {7.48863e+3,-1.21581e+4, 5.50281e+3,-1.33630e+3, 5.05367e-2}};

  G4double elmcof[4][5] = {{1.84542e+3,-5.64002e+3, 5.66730e+3,-3.15150e+3, 9.54160e+2},
			   {-2.24577e+3, 8.56133e+3,-9.67348e+3, 5.81744e+3,-1.86997e+3},
			   {2.79772e+3,-8.73073e+3, 9.19706e+3,-4.91900e+3, 1.37283e+3},
			   {-3.01866e+1, 1.41161e+3,-2.85919e+3, 2.13016e+3,-6.49072e+2}};

  G4double emxcof[4][6] = {{9.43596e4,-2.241997e5,2.223237e5,-1.324408e5,4.68922e4,-8.83568e3},
			   {-1.655827e5,4.062365e5,-4.236128e5,2.66837e5,-9.93242e4,1.90644e4},
			   {1.705447e5,-4.032e5,3.970312e5,-2.313704e5,7.81147e4,-1.322775e4},
			   {-9.274555e4,2.278093e5,-2.422225e5,1.55431e5,-5.78742e4,9.97505e3}};

  G4double elzcof[7][7] = {{5.11819909e+5,-1.30303186e+6, 1.90119870e+6,-1.20628242e+6, 5.68208488e+5, 5.48346483e+4,-2.45883052e+4},
			   {-1.13269453e+6, 2.97764590e+6,-4.54326326e+6, 3.00464870e+6, -1.44989274e+6,-1.02026610e+5, 6.27959815e+4},
			   {1.37543304e+6,-3.65808988e+6, 5.47798999e+6,-3.78109283e+6, 1.84131765e+6, 1.53669695e+4,-6.96817834e+4},
			   {-8.56559835e+5, 2.48872266e+6,-4.07349128e+6, 3.12835899e+6, -1.62394090e+6, 1.19797378e+5, 4.25737058e+4},
			   {3.28723311e+5,-1.09892175e+6, 2.03997269e+6,-1.77185718e+6, 9.96051545e+5,-1.53305699e+5,-1.12982954e+4},
			   {4.15850238e+4, 7.29653408e+4,-4.93776346e+5, 6.01254680e+5, -4.01308292e+5, 9.65968391e+4,-3.49596027e+3},
			   {-1.82751044e+5, 3.91386300e+5,-3.03639248e+5, 1.15782417e+5, -4.24399280e+3,-6.11477247e+3, 3.66982647e+2}};

  const G4int sizex = 5;
  const G4int sizey = 6;
  const G4int sizez = 4;

  G4double egscof[sizey][sizey][sizez];

  G4double egs1[sizey][sizex] = {{1.927813e5, 7.666859e5, 6.628436e5, 1.586504e5,-7.786476e3},
				 {-4.499687e5,-1.784644e6,-1.546968e6,-4.020658e5,-3.929522e3},
				 {4.667741e5, 1.849838e6, 1.641313e6, 5.229787e5, 5.928137e4},
				 {-3.017927e5,-1.206483e6,-1.124685e6,-4.478641e5,-8.682323e4},
				 {1.226517e5, 5.015667e5, 5.032605e5, 2.404477e5, 5.603301e4},
				 {-1.752824e4,-7.411621e4,-7.989019e4,-4.175486e4,-1.024194e4}};

  G4double egs2[sizey][sizex] = {{-6.459162e5,-2.903581e6,-3.048551e6,-1.004411e6,-6.558220e4},
				 {1.469853e6, 6.564615e6, 6.843078e6, 2.280839e6, 1.802023e5},
				 {-1.435116e6,-6.322470e6,-6.531834e6,-2.298744e6,-2.639612e5},
				 {8.665296e5, 3.769159e6, 3.899685e6, 1.520520e6, 2.498728e5},      
				 {-3.302885e5,-1.429313e6,-1.512075e6,-6.744828e5,-1.398771e5},
				 {4.958167e4, 2.178202e5, 2.400617e5, 1.167815e5, 2.663901e4}};

  G4double egs3[sizey][sizex] = {{3.117030e5, 1.195474e6, 9.036289e5, 6.876190e4,-6.814556e4},
				 {-7.394913e5,-2.826468e6,-2.152757e6,-2.459553e5, 1.101414e5},
				 {7.918994e5, 3.030439e6, 2.412611e6, 5.228065e5, 8.542465e3},
				 {-5.421004e5,-2.102672e6,-1.813959e6,-6.251700e5,-1.184348e5},
				 {2.370771e5, 9.459043e5, 9.026235e5, 4.116799e5, 1.001348e5},
				 {-4.227664e4,-1.738756e5,-1.795906e5,-9.292141e4,-2.397528e4}};

  G4double egs4[sizey][sizex] = {{-1.072763e5,-5.973532e5,-6.151814e5, 7.371898e4, 1.255490e5},
				 {2.298769e5, 1.265001e6, 1.252798e6,-2.306276e5,-2.845824e5},
				 {-2.093664e5,-1.100874e6,-1.009313e6, 2.705945e5, 2.506562e5},
				 {1.274613e5, 6.190307e5, 5.262822e5,-1.336039e5,-1.115865e5},
				 {-5.715764e4,-2.560989e5,-2.228781e5,-3.222789e3, 1.575670e4},
				 {1.189447e4, 5.161815e4, 4.870290e4, 1.266808e4, 2.069603e3}};

  for(i = 0; i < sizey; i++) {
    for(j = 0; j < sizex; j++) {
      egscof[i][j][0] = egs1[i][j];
      egscof[i][j][1] = egs2[i][j];
      egscof[i][j][2] = egs3[i][j];
      egscof[i][j][3] = egs4[i][j];
    }
  }

  // the program starts here                                           
  if (iz < 19  ||  iz > 111) {
    goto barfit900;
  }

  if(iz > 102   &&  il > 0) {
    goto barfit902;
  }

  z=G4double(iz);
  a=G4double(ia);
  el=G4double(il);
  amin= 1.2e0*z + 0.01e0*z*z;
  amax= 5.8e0*z - 0.024e0*z*z;

  if(a  <  amin  ||  a  >  amax) {
    goto barfit910;
  }

  // angul.mom.zero barrier                 
  aa=2.5e-3*a;
  zz=1.0e-2*z;
  ell=1.0e-2*el;
  bfis0 = 0.0;
  lpoly(zz,7,pz);
  lpoly(aa,7,pa);

  for(i = 0; i < 7; i++) { //do 10 i=1,7                                                       
    for(j = 0; j < 7; j++) { //do 10 j=1,7                                                       
      bfis0=bfis0+elzcof[j][i]*pz[i]*pa[j];
    }
  }

  bfis=bfis0;
  
  (*sbfis)=bfis;
  egs=0.0;
  (*segs)=egs;

  // values of l at which the barrier        
  // is 20%(el20) and 80%(el80) of l=0 value    
  amin2 = 1.4e0*z + 0.009e0*z*z;
  amax2 = 20.e0 + 3.0e0*z;

  if((a < amin2-5.e0  ||  a > amax2+10.e0) &&  il > 0) {
    goto barfit920;
  }

  lpoly(zz,5,pz);
  lpoly(aa,4,pa);
  el80=0.0;
  el20=0.0;
  elmax=0.0;

  for(i = 0; i < 4; i++) {
    for(j = 0; j < 5; j++) {
            el80 = el80 + elmcof[i][j]*pz[j]*pa[i];
            el20 = el20 + emncof[i][j]*pz[j]*pa[i];
    }
  }

  sel80 = el80;
  sel20 = el20;

  // value of l (elmax) where barrier disapp.
  lpoly(zz,6,pz);
  lpoly(ell,9,pl);

  for(i = 0; i < 4; i++) { //do 30 i= 1,4                                                      
    for(j = 0; j < 6; j++) { //do 30 j=1,6
      elmax = elmax + emxcof[i][j]*pz[j]*pa[i];
    }
  }

  (*selmax)=elmax;

  // value of barrier at ang.mom.  l          
  if(il < 1){
    return;                                                
  }

  x = sel20/(*selmax);
  y = sel80/(*selmax);
  
  if(el <= sel20) {
    // low l              
    q = 0.2/(std::pow(sel20,2)*std::pow(sel80,2)*(sel20-sel80));
    qa = q*(4.0*std::pow(sel80,3) - std::pow(sel20,3));
    qb = -q*(4.0*std::pow(sel80,2) - std::pow(sel20,2));
    bfis = bfis*(1.0 + qa*std::pow(el,2) + qb*std::pow(el,3));
  }
  else {
    // high l             
    aj = (-20.0*std::pow(x,5) + 25.e0*std::pow(x,4) - 4.0)*std::pow((y-1.0),2)*y*y;
    ak = (-20.0*std::pow(y,5) + 25.0*std::pow(y,4) - 1.0) * std::pow((x-1.0),2)*x*x;
    q = 0.2/(std::pow((y-x)*((1.0-x)*(1.0-y)*x*y),2));
    qa = q*(aj*y - ak*x);
    qb = -q*(aj*(2.0*y + 1.0) - ak*(2.0*x + 1.0));
    z = el/(*selmax);
    a1 = 4.0*std::pow(z,5) - 5.0*std::pow(z,4) + 1.0;
    a2 = qa*(2.e0*z + 1.e0);
    bfis=bfis*(a1 + (z - 1.e0)*(a2 + qb*z)*z*z*(z - 1.e0));
  }
  
  if(bfis <= 0.0) {
    bfis=0.0;
  }

  if(el > (*selmax)) {
    bfis=0.0;
  }
  (*sbfis)=bfis;

  // now calculate rotating ground state energy                        
  if(el > (*selmax)) {
    return;                                           
  }

  for(k = 0; k < 4; k++) {
    for(l = 0; l < 6; l++) {
      for(m = 0; m < 5; m++) {
	egs = egs + egscof[l][m][k]*pz[l]*pa[k]*pl[2*m];
      }
    }
  }

  (*segs)=egs;
  if((*segs) < 0.0) {
    (*segs)=0.0;
  }

  return;                                                            

 barfit900:  //continue                                                          
  (*sbfis)=0.0;
  // for z<19 sbfis set to 1.0e3                                            
  if (iz < 19)  {
    (*sbfis) = 1.0e3;
  }
  (*segs)=0.0;
  (*selmax)=0.0;
  return;                                                            

 barfit902:
  (*sbfis)=0.0;
  (*segs)=0.0;
  (*selmax)=0.0;
  return;                                                            

 barfit910:
  (*sbfis)=0.0;
  (*segs)=0.0;
  (*selmax)=0.0;
  return;                                                           

 barfit920:
  (*sbfis)=0.0;
  (*segs)=0.0;
  (*selmax)=0.0;
  return;                                                            
}

G4double G4Abla::erf(G4double x)
{
 G4double ferf;

 if(x<0.){
   ferf=-gammp(0.5,x*x);
 }else{
   ferf=gammp(0.5,x*x);;
 }
 return ferf;
}

G4double G4Abla::gammp(G4double a, G4double x)
{
 G4double fgammp;
 G4double gammcf,gamser,gln=0.;

 if(x<0.0 || a<=0.0)std::cout << "G4Abla::gammp = bad arguments in gammp" << std::endl;
 if(x<a+1.){
  gser(&gamser,a,x,gln);
  fgammp=gamser;
 }else{
  gcf(&gammcf,a,x,gln);
  fgammp=1.-gammcf;
 }
 return fgammp;
}

void G4Abla::gcf(G4double *gammcf,G4double a,G4double x,G4double gln)
{
 G4double fgammcf,del;
 G4double eps=3e-7;
 G4double fpmin=1e-30;
 G4int itmax=100;
 G4double an,b,c,d,h;

 gln=gammln(a);
 b=x+1.-a;
 c=1./fpmin;
 d=1./b;
 h=d;
 for(G4int i=1;i<=itmax;i++){
  an=-i*(i-a);
  b=b+2.;
  d=an*d+b;
  if(std::fabs(d)<fpmin)d=fpmin;
  c=b+an/c;
  if(std::fabs(c)<fpmin)c=fpmin;
  d=1.0/d;
  del=d*c;
  h=h*del;
  if(std::fabs(del-1.)<eps)goto dir1;
 }
 std::cout << "a too large, ITMAX too small in gcf" << std::endl;
 dir1:
 fgammcf=std::exp(-x+a*std::log(x)-gln)*h;
 (*gammcf)=fgammcf;
 return;
}

void G4Abla::gser(G4double *gamser,G4double a,G4double x,G4double gln)
{
 G4double fgamser,ap,sum,del;
 G4double eps=3e-7;
 G4int itmax=100;

 gln=gammln(a);
 if(x<=0.){
   if(x<0.)std::cout << "G4Abla::gser = x < 0 in gser" << std::endl;
   (*gamser)=0.0;
   return;
 }
 ap=a;
 sum=1./a;
 del=sum;
 for(G4int n=0;n<itmax;n++){
 ap=ap+1.;
 del=del*x/ap;
 sum=sum+del;
 if(std::fabs(del)<std::fabs(sum)*eps)goto dir1;
 }
 std::cout << "a too large, ITMAX too small in gser" << std::endl;
 dir1:
 fgamser=sum*std::exp(-x+a*std::log(x)-gln);
 (*gamser)=fgamser;
 return;
}

G4double G4Abla::gammln(G4double xx)
{
 G4double fgammln,x,ser,tmp,y;
 G4double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,
-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
 G4double stp=2.5066282746310005;

 x=xx;
 y=x;
 tmp=x+5.5;
 tmp=(x+0.5)*std::log(tmp)-tmp;
 ser=1.000000000190015;
 for(G4int j=0;j<6;j++){
 y=y+1.;
 ser=ser+cof[j]/y;
 }

 return fgammln=tmp+std::log(stp*ser/x);
}


G4double G4Abla::fd(G4double E)
{
  // DISTRIBUTION DE MAXWELL

  return (E*std::exp(-E));
}

G4double G4Abla::f(G4double E)
{
  // FONCTION INTEGRALE DE FD(E)
  return (1.0 - (E + 1.0) * std::exp(-E));
}

G4double G4Abla::fmaxhaz(G4double x)
{
 return ( -x*std::log(G4AblaRandom::flat()) -x*std::log(G4AblaRandom::flat()) -x*std::log(G4AblaRandom::flat()) ) ;
}

G4double G4Abla::fmaxhaz_old(G4double T)
{
  // tirage aleatoire dans une maxwellienne
  // t : temperature
  //
  // declaration des variables
  //

  const G4int pSize = 101;
  G4double p[pSize];

  // ial generateur pour le cascade (et les iy pour eviter les correlations)
  G4int i = 0;
  G4int itest = 0;
  // programme principal

  // calcul des p(i) par approximation de newton
  p[pSize-1] = 8.0;
  G4double x = 0.1;
  G4double x1 = 0.0;
  G4double y = 0.0;

  if (itest == 1) {
    goto fmaxhaz120;
  }

  for(i = 1; i <= 99; i++) {
  fmaxhaz20:
    x1 = x - (f(x) - G4double(i)/100.0)/fd(x);
    x = x1;
    if (std::fabs(f(x) - G4double(i)/100.0) < 1e-5) {
      goto fmaxhaz100;
    }
    goto fmaxhaz20;
  fmaxhaz100:
    p[i] = x;
  } //end do

  //  itest = 1;
  itest = 0;
  // tirage aleatoire et calcul du x correspondant 
  // par regression lineaire
 fmaxhaz120:
    y = G4AblaRandom::flat();
  i = nint(y*100);

  //   2590	c ici on evite froidement les depassements de tableaux....(a.b. 3/9/99)        
  if(i == 0) {
    goto fmaxhaz120;
  }

  if (i == 1) {
    x = p[i]*y*100;
  }
  else {
    x = (p[i] - p[i-1])*(y*100 - i) + p[i];
  }

  return(x*T);
}

G4double G4Abla::pace2(G4double a, G4double z)
{
  // PACE2
  // Cette fonction retourne le defaut de masse du noyau A,Z en MeV
  // Revisee pour a, z flottants 25/4/2002	                       =

  G4double fpace2 = 0.0;

  G4int ii = idint(a+0.5);
  G4int jj = idint(z+0.5);

  if(ii <= 0 || jj < 0) {
    fpace2=0.;
    return fpace2;
  }

  if(jj > 300) {
    fpace2=0.0;
  }
  else {
    fpace2=pace->dm[ii][jj];
  }
  fpace2=fpace2/1000.;

  if(pace->dm[ii][jj] == 0.) {
    if(ii < 12) {
      fpace2=-500.;
    }
    else {
      guet(&a, &z, &fpace2);
      fpace2=fpace2-ii*931.5;
      fpace2=fpace2/1000.;
    }
  }

  return fpace2;
}

void G4Abla::guet(G4double *x_par, G4double *z_par, G4double *find_par)
{
  // TABLE DE MASSES ET FORMULE DE MASSE TIRE DU PAPIER DE BRACK-GUET
  // Gives the theoritical value for mass excess...
  // Revisee pour x, z flottants 25/4/2002

  //real*8 x,z
  //	dimension q(0:50,0:70)
  G4double x = (*x_par);
  G4double z = (*z_par);
  G4double find = (*find_par);

  const G4int qrows = 50;
  const G4int qcols = 70;
  G4double q[qrows][qcols];
  for(G4int init_i = 0; init_i < qrows; init_i++) {
    for(G4int init_j = 0; init_j < qcols; init_j++) {
      q[init_i][init_j] = 0.0;
    }
  }

  G4int ix=G4int(std::floor(x+0.5));
  G4int iz=G4int(std::floor(z+0.5));
  G4double zz = iz;
  G4double xx = ix;
  find = 0.0;
  G4double avol = 15.776;
  G4double asur = -17.22;
  G4double ac = -10.24;
  G4double azer = 8.0;
  G4double xjj = -30.03;
  G4double qq = -35.4;
  G4double c1 = -0.737;
  G4double c2 = 1.28;

  if(ix <= 7) {
    q[0][1]=939.50;
    q[1][1]=938.21;
    q[1][2]=1876.1;
    q[1][3]=2809.39;
    q[2][4]=3728.34;
    q[2][3]=2809.4;
    q[2][5]=4668.8;
    q[2][6]=5606.5;
    q[3][5]=4669.1;
    q[3][6]=5602.9;
    q[3][7]=6535.27;
    q[4][6]=5607.3;
    q[4][7]=6536.1;
    q[5][7]=6548.3;
    find=q[iz][ix];
  }
  else {
    G4double xneu=xx-zz;
    G4double si=(xneu-zz)/xx;
    G4double x13=std::pow(xx,.333);
    G4double ee1=c1*zz*zz/x13;
    G4double ee2=c2*zz*zz/xx;
    G4double aux=1.+(9.*xjj/4./qq/x13);
    G4double ee3=xjj*xx*si*si/aux;
    G4double ee4=avol*xx+asur*(std::pow(xx,.666))+ac*x13+azer;
    G4double tota = ee1 + ee2 + ee3 + ee4;
    find = 939.55*xneu+938.77*zz - tota;
  }

  (*x_par) = x;
  (*z_par) = z;
  (*find_par) = find;
}
//

void G4Abla::FillData(G4int IMULTBU,G4int IEV_TAB){

    const G4double c = 29.9792458;
    const G4double fmp = 938.27231,fmn=939.56563,fml=1115.683; 

    varntp->ntrack = IMULTBU + IEV_TAB;

    G4int intp=0;

    for(G4int i=0;i<IMULTBU;i++){
 
    G4int iz = nint(BU_TAB[i][7]);
    G4int ia = nint(BU_TAB[i][8]);
    G4int is = nint(BU_TAB[i][11]);

    Ainit = Ainit + ia;
    Zinit = Zinit + iz;
    Sinit = Sinit - is;

    varntp->zvv[intp] = iz;
    varntp->avv[intp] = ia;
    varntp->svv[intp] = -1*is;
    varntp->itypcasc[intp] = 0;

    G4double v2 = BU_TAB[i][4]*BU_TAB[i][4]+BU_TAB[i][5]*BU_TAB[i][5]+BU_TAB[i][6]*BU_TAB[i][6];
    G4double gamma = std::sqrt(1.0 - v2 / (c*c));
    G4double avvmass = iz*fmp + (ia-iz-is)*fmn + is*fml + eflmac(ia,iz,0,3);
    G4double etot = avvmass / gamma;
    varntp->pxlab[intp] = etot * BU_TAB[i][4] / c;
    varntp->pylab[intp] = etot * BU_TAB[i][5] / c;
    varntp->pzlab[intp] = etot * BU_TAB[i][6] / c;
    varntp->enerj[intp] = etot - avvmass;
    intp++;
    }


    for(G4int i=0;i<IEV_TAB;i++){
 
    G4int iz = nint(EV_TAB[i][0]);
    G4int ia = nint(EV_TAB[i][1]);
    G4int is = EV_TAB[i][5];

    varntp->itypcasc[intp] = 0;

     if(ia>0){// normal particles
     varntp->zvv[intp] = iz;
     varntp->avv[intp] = ia;
     varntp->svv[intp] = -1*is;
     Ainit = Ainit + ia;
     Zinit = Zinit + iz;
     Sinit = Sinit - is;
     G4double v2 = EV_TAB[i][2]*EV_TAB[i][2]+EV_TAB[i][3]*EV_TAB[i][3]+EV_TAB[i][4]*EV_TAB[i][4];
     G4double gamma = std::sqrt(1.0 - v2 / (c*c));
     G4double avvmass = iz*fmp + (ia-iz-is)*fmn + is*fml + eflmac(ia,iz,0,3);
     G4double etot = avvmass / gamma;
     varntp->pxlab[intp] = etot * EV_TAB[i][2] / c;
     varntp->pylab[intp] = etot * EV_TAB[i][3] / c;
     varntp->pzlab[intp] = etot * EV_TAB[i][4] / c;
     varntp->enerj[intp] = etot - avvmass;
     }else if(ia==-2){// lambda0
     varntp->zvv[intp] = 0;
     varntp->avv[intp] = 1;
     varntp->svv[intp] = -1;
     Ainit = Ainit + 1;
     Sinit = Sinit - 1;
     G4double v2 = EV_TAB[i][2]*EV_TAB[i][2]+EV_TAB[i][3]*EV_TAB[i][3]+EV_TAB[i][4]*EV_TAB[i][4];
     G4double gamma = std::sqrt(1.0 - v2 / (c*c));
     G4double avvmass = fml;
     G4double etot = avvmass / gamma;
     varntp->pxlab[intp] = etot * EV_TAB[i][2] / c;
     varntp->pylab[intp] = etot * EV_TAB[i][3] / c;
     varntp->pzlab[intp] = etot * EV_TAB[i][4] / c;
     varntp->enerj[intp] = etot - avvmass;
     }else{// photons
     varntp->zvv[intp] = iz;
     varntp->avv[intp] = ia;
     varntp->svv[intp] = 0;
     Ainit = Ainit + ia;
     Zinit = Zinit + iz;
     Sinit = Sinit - is;
     varntp->pxlab[intp] = EV_TAB[i][2];
     varntp->pylab[intp] = EV_TAB[i][3];
     varntp->pzlab[intp] = EV_TAB[i][4];
     varntp->enerj[intp] = std::sqrt(EV_TAB[i][2]*EV_TAB[i][2]+EV_TAB[i][3]*EV_TAB[i][3]+EV_TAB[i][4]*EV_TAB[i][4]);
     }
    intp++;
    }
// 
 return;
}

// Utilities

G4double G4Abla::min(G4double a, G4double b)
{
  if(a < b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4Abla::min(G4int a, G4int b)
{
  if(a < b) {
    return a;
  }
  else {
    return b;
  }
}

G4double G4Abla::max(G4double a, G4double b)
{
  if(a > b) {
    return a;
  }
  else {
    return b;
  }
}

G4int G4Abla::max(G4int a, G4int b)
{
  if(a > b) {
    return a;
  }
  else {
    return b;
  }
}

G4double G4Abla::DSIGN(G4double a, G4double b){
// A function that assigns the sign of the second argument to the
// absolute value of the first

 if(b>=0){ 
  return std::abs(a);
 }else{
  return -1.0*std::abs(a);
 }
 return 0;
}

G4int G4Abla::ISIGN(G4int a, G4int b){
// A function that assigns the sign of the second argument to the
// absolute value of the first

 if(b>=0){ 
  return std::abs(a);
 }else{
  return -1*std::abs(a);
 }
 return 0;
}

G4int G4Abla::nint(G4double number)
{
  G4double intpart = 0.0;
  G4double fractpart = 0.0;
  fractpart = std::modf(number, &intpart);
  if(number == 0) {
    return 0;
  }
  if(number > 0) {
    if(fractpart < 0.5) {
      return G4int(std::floor(number));
    }
    else {
      return G4int(std::ceil(number));
    }
  }
  if(number < 0) {
    if(fractpart < -0.5) {
      return G4int(std::floor(number));
    }
    else {
      return G4int(std::ceil(number));
    }
  }

  return G4int(std::floor(number));
}

G4int G4Abla::secnds(G4int x)
{
  time_t mytime;
  tm *mylocaltime;

  time(&mytime);
  mylocaltime = localtime(&mytime);

  if(x == 0) {
    return(mylocaltime->tm_hour*60*60 + mylocaltime->tm_min*60 + mylocaltime->tm_sec);
  }
  else {
    return G4int(mytime - x);
  }
}

G4int G4Abla::mod(G4int a, G4int b)
{
  if(b != 0) {
    return a%b;
  }
  else {
    return 0;
  } 
}

G4double G4Abla::dint(G4double x)
{
  G4double value = 0.0;
/*
  if(a < 0.0) {
    value = double(std::ceil(a));
  }
  else {
    value = double(std::floor(a));
  }
*/
	if(x-std::floor(x) <= std::ceil(x)-x) 
         value = G4double(std::floor(x));
	else 
         value = G4double(std::ceil(x));

  return value;
}

G4int G4Abla::idint(G4double x)
{
  G4int value = 0;
	if(x-std::floor(x) <= std::ceil(x)-x) 
         value = G4int(std::floor(x));
	else 
         value = G4int(std::ceil(x));

  return value;
}

G4int G4Abla::idnint(G4double x)
{
	if(x-std::floor(x) <= std::ceil(x)-x) 
         return G4int(std::floor(x));
	else 
         return G4int(std::ceil(x));
}

G4double G4Abla::dmin1(G4double a, G4double b, G4double c)
{
  if(a < b && a < c) {
    return a;
  }
  if(b < a && b < c) {
    return b;
  }
  if(c < a && c < b) {
    return c;
  }
  return a;
}

G4double G4Abla::utilabs(G4double a)
{
  return std::abs(a);
}


G4double G4Abla::width(G4double AMOTHER,G4double ZMOTHER,G4double APART,G4double ZPART,G4double TEMP,G4double B1,G4double SB1,G4double EXC)
{
/*
* Implemented by JLRS for Abla c++: 06/11/2016
*
C  Last update:
C       28/10/13 - JLRS - from abrablav4 (AK)
*/
      G4int IZPART,IAPART,NMOTHER;
      G4double B,HBAR,PI,RGEOM,MPART,SB;
      G4double BKONST,C,C2,G,APARTNER,MU;
      G4double INT1,INT2,INT3,AKONST,EARG,R0,MPARTNER;
      G4double AEXP;
      G4double ARG;
      G4double PAR_A1=0.,PAR_B1=0.,FACT=1.;
      G4double fwidth=0.;
      G4int idlamb0=0;
      PI=3.141592654;

      if(ZPART==-2.){
       ZPART=0.;
       idlamb0=1;
      }

      IZPART = idnint(ZPART);
      IAPART = idnint(APART);

      B = B1;
      SB = SB1;
      NMOTHER = idnint(AMOTHER-ZMOTHER);

      PAR_A1 = 0.0;
      PAR_B1 = 0.0;

      if(SB>EXC){
       return fwidth=0.0;
      }else{
// in MeV*s
      HBAR = 6.582122e-22;
//      HBAR2 = HBAR * HBAR
// in m/s
      C = 2.99792458e8;
      C2 = C * C;
      APARTNER = AMOTHER - APART;
      MPARTNER = APARTNER * 931.49 / C2;

//           g=(2s+1)
      if(IAPART==1&&IZPART==0){
        G = 2.0;
        MPART =  939.56 / C2;
        if(idlamb0==1)MPART =  1115.683 / C2;
      }else{
       if(IAPART==1&&IZPART==1){
        G = 2.0;
        MPART = 938.27 / C2;
       }
       else{
        if(IAPART==2&&IZPART==0){
        G = 1.0;
        MPART = 2.*939.56 / C2;
        }else{
         if(IAPART==2&&IZPART==1){
         G = 3.0;
         MPART = 1876.10 / C2;
         }else{
          if(IAPART==3&&IZPART==1){
           G = 2.0;
           MPART = 2809.39 / C2;
          }else{
           if(IAPART==3&&IZPART==2){
            G = 2.0;
            MPART = 2809.37 / C2;
           }else{
            if(IAPART==4&&IZPART==2){
             G = 1.0;
             MPART = 3728.35 / C2;
            }else{
             // IMF
             G = 1.0;
             MPART = APART * 931.49 / C2;   
            }
           }
          }
         }
        }
       }
      }//end g

// Relative mass in MeV*s^2/m^2
      MU = MPARTNER * MPART / (MPARTNER + MPART);
// in m
      R0 = 1.16e-15;

      RGEOM = R0 * (std::pow(APART,1.0/3.0)+std::pow(AMOTHER-APART,1.0/3.0));

// in m*sqrt(MeV)
      AKONST = HBAR*std::sqrt(1.0 / MU);

// in  1/(MeV*m^2)
      BKONST =  MPART / ( PI * PI * HBAR * HBAR);
//
// USING ANALYTICAL APPROXIMATION

      INT1 = 2.0 * std::pow(TEMP,3.) / (2.0 * TEMP + B);

      ARG = std::sqrt(B/TEMP);
      EARG = (erf(ARG) - 1.0);
      if(std::abs(EARG)<1.e-9) EARG = 0.0;
      if(B==0.0){
        INT2 = 0.5 * std::sqrt(PI) * std::pow(TEMP,3.0/2.0);
      }else{
         AEXP = B/TEMP;
          if(AEXP>700.0) AEXP = 700.0;
         INT2 = (2.0*B*B +TEMP*B)/std::sqrt(B) + std::exp(AEXP) * std::sqrt(PI/(4.0*TEMP))*(4.0*B*B+4.0*B*TEMP - TEMP*TEMP) *EARG;
       if(INT2<0.0) INT2 = 0.0;
// For very low temperatures when EARG=0, INT2 get unreasonably high values
// comming from the first term. Therefore, for these cases INT2 is set to 0.
       if(EARG==0.0) INT2 = 0.0;
      }//if B

      INT3 = 2.0*TEMP*TEMP*TEMP / (2.0*TEMP*TEMP + 4.0*B*TEMP + B*B);

      if(IZPART<-1.0&&ZMOTHER<151.0){
//      IF(IZPART.LT.1)THEN
// For neutrons, the width is given by a mean value between geometrical and QM values;
// Only QM contribution (Rgeom -> Rgeom + Rlamda) seems to be too strong for neutrons
       fwidth = PI * BKONST *  G * std::sqrt((RGEOM * RGEOM * INT1 + 2.0 * AKONST * RGEOM * INT2 + AKONST * AKONST * INT3) * RGEOM * RGEOM * INT1);

      }else{
       fwidth = PI * BKONST *  G *(RGEOM * RGEOM * INT1 + 2.0 * AKONST * RGEOM * INT2 + AKONST * AKONST * INT3);
      }


// To correct for too high values of analytical width compared to
// numerical solution for energies close to the particle threshold:
       if(IZPART<3.0){
        if(AMOTHER<155.0){
         PAR_A1=std::exp(2.302585*0.2083*std::exp(-0.01548472*AMOTHER))-0.05;
         PAR_B1 = 0.59939389 + 0.00915657 * AMOTHER;
        }else{
         if(AMOTHER>154.0&&AMOTHER<195.0){
           PAR_A1=1.0086961-8.629e-5*AMOTHER;
           PAR_B1 = 1.5329331 + 0.00302074 * AMOTHER;
         }else{
          if(AMOTHER>194.0&&AMOTHER<208.0){
           PAR_A1=9.8356347-0.09294663*AMOTHER+2.441e-4*AMOTHER*AMOTHER;
           PAR_B1 = 7.7701987 - 0.02897401 * AMOTHER;
          }else{
           if(AMOTHER>207.0&&AMOTHER<228.0){
            PAR_A1=15.107385-0.12414415*AMOTHER+2.7222e-4*AMOTHER*AMOTHER;
            PAR_B1=-64.078009+0.56813179*AMOTHER-0.00121078*AMOTHER*AMOTHER; 
           }else{
             if(AMOTHER>227.0){
              if(mod(NMOTHER,2)==0&&NMOTHER>147.){
               PAR_A1 = 2.0*(0.9389118 + 6.4559e-5 * AMOTHER);
              }else{
               if(mod(NMOTHER,2)==1)PAR_A1 = 3.0*(0.9389118 + 6.4559e-5 * AMOTHER);
              }
              PAR_B1 = 2.1507177 + 0.00146119 * AMOTHER;
             }
           }
          }
         }
        }
       FACT = std::exp((2.302585*PAR_A1*std::exp(-PAR_B1*(EXC-SB))));
       if(FACT<1.0) FACT = 1.0;
       if(IZPART<-1.&&ZMOTHER<151.0){
//       IF(IZPART.LT.1)THEN
        fwidth = fwidth / std::sqrt(FACT);
       }else{
        fwidth = fwidth / FACT;
       }
       }//if IZPART<3.0

       if(fwidth<=0.0){
       std::cout <<"LOOK IN PARTICLE_WIDTH!" << std::endl;
       std::cout <<"ACN,APART :"<< AMOTHER << APART << std::endl;
       std::cout <<"EXC,TEMP,B,SB :" << EXC << " " << TEMP << " "  << B << " "  << SB << std::endl;
       std::cout <<"INTi, i=1-3 :" << INT1 << " "  << INT2 << " "  << INT3 << std::endl;
       std::cout <<" " << std::endl;
       }

      }//if SB>EXC
  return fwidth;
}

G4double G4Abla::pen(G4double A, G4double ap, G4double omega, G4double T)
{
// JLRS: 06/11/2016
// CORRECTIONS FOR BARRIER PENETRATION
// AK, KHS 2005 - Energy-dependen inverse cross sections included, influence of
//                Coulomb barrier for LCP, tunnelling for LCP

 G4double fpen=0., MU, HO;

// REDUCED MASSES (IN MeV/C**2)
       MU = (A - ap) * ap / A;

// ENERGY OF THE INVERSE PARABOLA AT THE POTENTIAL BARRIER (hbar*omega);
// HERE hbar = 197.3287 fm*MeV/c, omega is in c/fm
       HO = 197.3287 * omega;

     if(T<=0.0){
       fpen = 0.0;
     }else{
       fpen=std::pow(10.0,4.e-4*std::pow(T/(HO*HO*std::pow(MU,0.25)),-4.3/2.3026));
     }

 return fpen;
}

void G4Abla::bsbkbc(G4double A,G4double Z,G4double *BS,G4double *BK,G4double *BC)
{
// Calculate BS and BK needed for a level-density parameter:
// BETA2 and BETA4 = quadrupole and hexadecapole deformation

      G4double PI = 3.14159265;
      G4int IZ = idnint(Z);
      G4int IN = idnint(A - Z);
// alphaN = sqrt(2*N/(4*pi))*BetaN
      G4double ALPHA2 = std::sqrt(5.0/(4.0*PI))*ecld->beta2[IN][IZ];
      G4double ALPHA4 = std::sqrt(9.0/(4.0*PI))*ecld->beta4[IN][IZ];

      (*BS) = 1.0 + 0.4*ALPHA2*ALPHA2 - 4.0/105.0*ALPHA2*ALPHA2*ALPHA2 - 66.0/175.0*ALPHA2*ALPHA2*ALPHA2*ALPHA2 - 4.0/35.0*ALPHA2*ALPHA2*ALPHA4 + ALPHA4*ALPHA4;

      (*BK) = 1.0 + 0.4*ALPHA2*ALPHA2 + 16.0/105.0*ALPHA2*ALPHA2*ALPHA2 - 82.0/175.0*ALPHA2*ALPHA2*ALPHA2*ALPHA2 + 2.0/35.0*ALPHA2*ALPHA2*ALPHA4 + ALPHA4*ALPHA4;

      (*BC)=0.0;      

 return;
}

G4double G4Abla::fvmaxhaz( G4double T)
{
// Random generator according to a distribution similar to a
// Maxwell distribution with quantum-mech. x-section for charged particles according to KHS
//      Y = X**(1.5E0) / (B+X) * EXP(-X/T) (approximation:)

return (3.0 * T * std::pow(-1.*std::log(G4AblaRandom::flat()) * std::log(G4AblaRandom::flat())*std::log(G4AblaRandom::flat()),0.333333));
}

G4double G4Abla::func_trans(G4double TIME,G4double ZF,G4double AF,G4double bet,G4double Y,G4double FT,G4double T_0)
{
/*
c   This function determines the fission width as a function o time
c   according to the analytical solution of the FPE for the probability distribution
c   at the barrier when the nucleus potential is aproximated by a parabolic
c   potential. It is taken from S. Chandrasekhar, Rev. Mod. Phys. 15 (1943) 1
c                   
c***********************INPUT PARAMETERS*********************************
c  Time               Time at which we evaluate the fission width
c  ZF                 Z of nucleus
C  AF                 A of nucleus		
c  BET                Reduced dissipation coefficient
c  FT                 Nuclear temperature
C**************************************************************************
C********************************OUTPUT***********************************
C   Fission decay width at the corresponding time of the decay cascade
C*************************************************************************
c****************************OTHER VARIABLES******************************
C  SIGMA_SQR         Square of the width of the prob. distribution
C  XB                Deformation of the nucleus at the saddle point
c  NORM              Normalization factor of the probability distribution
c  W                 Probability distribution at the saddle deformation XB
c  W_INFIN           Probability distr. at XB at infinite time
c  MFCD              Mass of the fission collective degree of freedom
C*************************************************************************
*/
      G4double PI = 3.14159;
      G4double DEFO_INIT,OMEGA,HOMEGA,OMEGA_GS,HOMEGA_GS,K1,MFCD;
      G4double BET1,XACT,SIGMA_SQR,W_EXP,XB,NORM,SIGMA_SQR_INF,W_INFIN,W;
      G4double FUNC_TRANS,LOG_SLOPE_INF,LOG_SLOPE_ABS;
//
// Influence of initial deformation
// Initial alpha2 deformation (GS)
      DEFO_INIT = std::sqrt(5.0/(4.0*PI))*ecld->beta2[fiss->at-fiss->zt][fiss->zt];
//
      fomega_sp(AF,Y,&MFCD,&OMEGA,&HOMEGA);
      fomega_gs(AF,ZF,&K1,&OMEGA_GS,&HOMEGA_GS);
//
// Determination of the square of the width of the probability distribution
// For the overdamped regime BET**2 > 4*OMEGA**2
         if((bet*bet)>4.0*OMEGA_GS*OMEGA_GS){
          BET1=std::sqrt(bet*bet-4.0*OMEGA_GS*OMEGA_GS);
//
// REMEMBER THAT HOMEGA IS ACTUALLY HBAR*HOMEGA1=1MeV
// SO THAT HOMEGA1 = HOMEGA/HBAR
//
          SIGMA_SQR = (FT/K1)*(1.0 -((2.0*bet*bet/(BET1*BET1)* (0.5 * (std::exp(0.50*(BET1-bet)*1.e21*TIME) - std::exp(0.5*(-BET1-bet)*1.e21*TIME)))*(0.5 * (std::exp(0.50*(BET1-bet)*1.e21*TIME) - std::exp(0.5*(-BET1-bet)*1.e21*TIME)))) + (bet/BET1*0.50 * (std::exp((BET1-bet)*1.e21*TIME)-std::exp((-BET1-bet)*1.e21*TIME))) + 1. * std::exp(-bet*1.e21*TIME)));
//
// Evolution of the mean x-value (KHS March 2006)
          XACT = DEFO_INIT *std::exp(-0.5*(bet-BET1)*1.e21*(TIME-T_0));
//
         }else{
// For the underdamped regime BET**2 < 4*HOMEGA**2 BET1 becomes a complex number
// and the expression with sinh and cosh can be transformed in one with sin and cos
          BET1=std::sqrt(4.0*OMEGA_GS*OMEGA_GS-bet*bet);
          SIGMA_SQR = FT/K1*(1.-std::exp(-1.0*bet*1.e21*TIME)*(bet*bet/(BET1*BET1)*(1.-std::cos(BET1*1.e21*TIME)) + bet/BET1*std::sin(BET1*1.e21*TIME) + 1.0));
          XACT = DEFO_INIT*std::cos(0.5*BET1*1.e21*(TIME-T_0))*std::exp(-bet*1.e21*(TIME-T_0));
         }

// Determination of the deformation at the saddle point according to
// "Geometrical relationships of Macroscopic Nucl. Phys." from Hass and Myers page 100
// This corresponds to alpha2 deformation.
          XB = 7./3.*Y-938./765.*Y*Y+9.499768*Y*Y*Y-8.050944*Y*Y*Y*Y;
//
// Determination of the probability distribution at the saddle deformation
//
          if(SIGMA_SQR>0.0){
           NORM = 1./std::sqrt(2.*PI*SIGMA_SQR);
//
           W_EXP = -1.*(XB - XACT)*(XB - XACT)/(2.0 * SIGMA_SQR);
            if(W_EXP<(-708.0) ) W_EXP = -708.0;
           W = NORM * std::exp( W_EXP ) * FT / (K1 * SIGMA_SQR);
           }else{
           W = 0.0;
           }
//
// Determination of the fission decay width, we assume we are in the overdamped regime
//
              SIGMA_SQR_INF = FT/K1;
              W_EXP = -XB*XB/(2.0 * SIGMA_SQR_INF);
              if(W_EXP<(-708.0))W_EXP = -708.0;
              W_INFIN = std::exp(W_EXP)/std::sqrt(2.0*PI*SIGMA_SQR_INF);
              FUNC_TRANS = W / W_INFIN;
//
// Correction for the variation of the mean velocity at the fission barrier
//  (see B. Jurado et al, Nucl. Phys. A747, p. 14)
//
              LOG_SLOPE_INF = cram(bet,HOMEGA)*bet*MFCD*OMEGA/FT;
              LOG_SLOPE_ABS = (XB-XACT)/SIGMA_SQR-XB/SIGMA_SQR_INF+cram(bet,HOMEGA)*bet*MFCD*OMEGA/FT;
//
              FUNC_TRANS = FUNC_TRANS * LOG_SLOPE_ABS/LOG_SLOPE_INF;
//
 return FUNC_TRANS;
}


void G4Abla::part_fiss(G4double BET,G4double GP,G4double GF,G4double Y,G4double TAUF,G4double TS1,G4double TSUM,G4int *CHOICE,G4double ZF,G4double AF,G4double FT,G4double *T_LAPSE,G4double *GF_LOC)
{
/*
C     THIS SUBROUTINE IS AIMED TO CHOOSE BETWEEN PARTICLE EMISSION
C     AND FISSION
C     WE USE MONTE-CARLO METHODS AND SAMPLE TIME BETWEEN T=0 AND T=1.5*TAUF
c	TO SIMULATE THE TRANSIENT TIME WITH 30 STEPS (0.05*TAUF EACH)
C     FOR t>1.5*TAUF , GF=CONSTANT=ASYMPTOTICAL VALUE (INCLUDING KRAMERS FACTOR)
c------------------------------------------------------------------------
c    Modifications introduced by BEATRIZ JURADO 18/10/01:
c    1. Now this subrutine is included in the rutine direct
c    2. TSUM does not include the current particle decay time
C    3. T_LAPSE is the time until decay, taken as an output variable
C    4. GF_LOC is also taken as an output variable
C    5. BET (Diss. Coeff.) and HOMEGA (Frequency at the ground state
c       are included as input variables because they are needed for FUNC_TRANS
C-----------------------------------------------------------------------
C     ON INPUT:
C       GP                 Partial particle decay width
C       GF                 Asymptotic value of Gamma-f, including Kramers factor
C       AF                 Mass number of nucleus
C       TAUF               Transient time
C       TS1                Partial particle decay time for the next step
C       TSUM               Total sum of partial particle decay times, including
C                               the next expected one, which is in competition
C                               with fission now
C       ZF                 Z of nucleus
C       AF                 A of nucleus
C-----------------------------------------------------------------------
C     ON OUTPUT:
C       CHOICE             Key for decay mode: 0 = no decay (only internal)
C                                              1 = evaporation
C                                              2 = fission
C-----------------------------------------------------------------------
C     VARIABLES:
C       GP                 Partial particle decay width
C       GF                 Asymptotic value of Gamma-f, including Kramers factor
C       TAUF               Transient time
C       TS1                Partial particle decay time
C       TSUM               Total sum of partial particle decay times
C       CHOICE              Key for decay mode
C       ZF                 Z of nucleus
C       AF                 A of nucleus
C       FT                 Used for Fermi function in FUNC_TRANS
C       STEP_LENGTH        Step in time to sample different decays
C       BEGIN_TIME         Total sum of partial particle decay times, excluding
C                               the next expected one, which is in competition
C                               with fission now
C       LOC_TIME_BEGIN     Begin of time interval considered in one step
C       LOC_TIME_END       End of time interval considered in one step
C       GF_LOC             In-grow function for fission width,
c                                 normalized to asymptotic value
C       TS2                Effective partial fission decay time in one time step
C       HBAR               hbar
C       T_LAPSE            Effective decay time in one time step
C       REAC_PROB          Reaction probability in one time step
C       X                  Help variable for random generator
C------------------------------------------------------------------------
*/
  G4double K1,OMEGA,HOMEGA,t_0,STEP_LENGTH,LOC_TIME_BEGIN,LOC_TIME_END=0.,BEGIN_TIME=0.,FISS_PROB,X,TS2,LAMBDA,REAC_PROB;
  G4double HBAR=6.582122e-22;
  G4int fchoice=0;
  G4double fGF_LOC=0.,fT_LAPSE=0.;
//
	if(GF<=0.0){
          *CHOICE = 1;
          *T_LAPSE=TS1;
          *GF_LOC = 0.0;
  	  goto direct107;
	}
//
      fomega_gs(AF,ZF,&K1,&OMEGA,&HOMEGA);
//
// ****************************************************************
//    Calculation of the shift in time due to the initial conditions
//
//    Overdamped regime
      if(BET*BET>4.0*OMEGA*OMEGA){
//         REMEMBER THAT HOMEGA IS ACTUALLY HBAR*HOMEGA1=1MeV
//         SO THAT HOMEGA1 = HOMEGA/HBAR
//     Additional factor 1/16 proposed by KHS on 14/7/2010. Takes into
//     account the fact that the curvature of the potential is ~16 times
//     larger than what predicted by the liquid drop model, because of
//     shell effects.
          t_0 = BET*1.e21*HBAR*HBAR/(4.*HOMEGA*FT)/16.;
       }else{
//     Underdamped regime
         if(((2.*FT-HOMEGA/16.)>0.000001) && BET>0.0){
//     Additional factor 1/16 proposed by KHS on 14/7/2010. Takes into
//     account the fact that the curvature of the potential is ~16 times
//     larger than what predicted by the liquid drop model, because of
//     shell effects.
           t_0 = (std::log(2.*FT/(2.*FT-HOMEGA/16.)))/(BET*1.e21);
         }else{
//     Neglect fission transients if the time shift t_0 is too
//     large. Suppresses large, spurious fission cross section at very
//     low excitation energy in p+Ta.
//
              fchoice = 0;
              goto direct106;
         }
       }
// ********************************************************************+
      fchoice = 0;
      STEP_LENGTH = 1.5*TAUF/50.;
//
//  AT FIRST WE CACULATE THE REAL CURRENT TIME
//  TSUM includes only the time elapsed in the previous steps
//
      BEGIN_TIME = TSUM + t_0;
//
      if(BEGIN_TIME<0.0) std::cout << "CURRENT TIME < 0" << BEGIN_TIME << std::endl;
//
      if(BEGIN_TIME<1.50*TAUF){
        LOC_TIME_BEGIN = BEGIN_TIME;
//
        while((LOC_TIME_BEGIN<1.5*TAUF)&&fchoice==0){

         LOC_TIME_END = LOC_TIME_BEGIN + STEP_LENGTH;
//
// NOW WE ESTIMATE THE MEAN VALUE OF THE FISSION WIDTH WITHIN THE SMALL INTERVAL
         fGF_LOC=(func_trans(LOC_TIME_BEGIN,ZF,AF,BET,Y,FT,t_0)+func_trans(LOC_TIME_END,ZF,AF,BET,Y,FT,t_0))/2.0;
//
         fGF_LOC = fGF_LOC * GF;

// TS2 IS THE MEAN DECAY TIME OF THE FISSION CHANNEL
                 if(fGF_LOC>0.0){
                     TS2 = HBAR/fGF_LOC;
                 }else{
                     TS2 = 0.0;
                 }
//
                 if(TS2>0.0){
                     LAMBDA  = 1.0/TS1 + 1.0/TS2;
                 }else{
                     LAMBDA = 1.0/TS1;
                 }
//
// This is the probability to survive the decay at this step
                 REAC_PROB = std::exp(-1.0*STEP_LENGTH*LAMBDA);
// I GENERATE A RANDOM NUMBER
                 X = G4AblaRandom::flat();
                 if(X>REAC_PROB){
// THEN THE EVAPORATION OR FISSION HAS OCCURED
                        FISS_PROB = fGF_LOC / (fGF_LOC+GP);
                        X = G4AblaRandom::flat();
//                       WRITE(6,*)'X=',X
                        if(X<FISS_PROB){
// FISSION OCCURED
                           fchoice = 2;
                        }else{
// EVAPORATION OCCURED
                           fchoice = 1;
                        }
                  }// if x
                  LOC_TIME_BEGIN = LOC_TIME_END;
       }// while
// Take the real decay time of this decay step
      fT_LAPSE = LOC_TIME_END - BEGIN_TIME;
      }// if BEGIN_TIME
//
// NOW, IF NOTHING HAPPENED DURING TRANSIENT TIME
  direct106:
    if(fchoice==0){
            fGF_LOC=GF;
            FISS_PROB = GF / (GF+GP);
 
// Added for cases where already at the beginning BEGIN_TIME > 1.5d0*TAUF
            if(GF>0.0){
                 TS2 = HBAR/GF;
            }else{
                 TS2 = 0.0;
            }

            if(TS2>0.0){
                 LAMBDA  = 1./TS1 + 1./TS2;
            }else{
                 LAMBDA = 1./TS1;
            }
//
            X = G4AblaRandom::flat();

            if(X<FISS_PROB){
// FISSION OCCURED
                  fchoice = 2;
            }else{
// EVAPORATION OCCURED
                  fchoice = 1;
            }
//
//TIRAGE ALEATOIRE DANS UNE EXPONENTIELLLE : Y=EXP(-X/T)
//       EXPOHAZ=-T*LOG(HAZ(K))
       fT_LAPSE = fT_LAPSE -1.0/LAMBDA*std::log(G4AblaRandom::flat());
      }
//
  direct107:

  (*T_LAPSE)=fT_LAPSE;
  (*GF_LOC)=fGF_LOC;
  (*CHOICE)=fchoice;
  return;
}

G4double G4Abla::tunnelling(G4double A,G4double ZPRF,G4double Y,G4double EE,G4double EF,G4double TEMP,G4double DENSG,G4double DENSF,G4double ENH_FACT)
{
// Subroutine to caluclate fission width with included effects
// of tunnelling through the fission barrier

      G4double PI = 3.14159;
      G4int IZ, IN;
      G4double MFCD,OMEGA,HOMEGA1,HOMEGA2=0.,GFTUN;
      G4double E1,E2,EXP_FACT,CORR_FUNCT,FACT1,FACT2,FACT3;

      IZ = idnint(ZPRF);
      IN = idnint(A-ZPRF);

// For low energies system "sees" LD barrier
      fomega_sp(A,Y,&MFCD,&OMEGA,&HOMEGA1);

      if(mod(IN,2)==0&&mod(IZ,2)==0){    // e-e
// Due to pairing gap, even-even nuclei cannot tunnel for excitation energy lower
// than pairing gap (no levels at which system can be)
      EE = EE - 12.0/std::sqrt(A);
      HOMEGA2 = 1.04;
      }

      if(mod(IN,2)==1&&mod(IZ,2)==1){   // o-o
      HOMEGA2 = 0.65;
      }

      if(mod(IN,2)==1&&mod(IZ,2)==0){   // o-e
      HOMEGA2 = 0.8;
      }

      if(mod(IN,2)==0&&mod(IZ,2)==1){   // e-0
      HOMEGA2 = 0.8;
      }

      E1 = EF + HOMEGA1/2.0/PI*std::log(HOMEGA1*(2.0*PI+HOMEGA2)/4.0/PI/PI);

      E2 = EF + HOMEGA2/(2.0*PI)*std::log(1.0+2.0*PI/HOMEGA2);

// AKH May 2013 - Due to approximations in the analytical integration, at energies
// just above barrier Pf was to low, at energies below
// barrier it was somewhat higher. LInes below are supposed to correct for this.
// Factor 0.20 in EXP_FACT comes from the slope of the Pf(Eexc) (Gavron's data)
// around fission barrier.
      EXP_FACT = (EE-EF)/(HOMEGA2/(2.0*PI));
      if(EXP_FACT>700.0) EXP_FACT = 700.0;
      CORR_FUNCT = HOMEGA1 * (1.0-1.0/(1.0+std::exp(EXP_FACT)));
      if(mod(IN,2)==0&&mod(IZ,2)==0){
      CORR_FUNCT = HOMEGA1 * (1.0-1.0/(1.0+std::exp(EXP_FACT)));
      }

      FACT1 = HOMEGA1/(2.0*PI*TEMP+HOMEGA1);
      FACT2 = (2.0*PI/(2.0*PI+HOMEGA2)-HOMEGA1*(2.0*PI+HOMEGA2)/4.0/PI/PI)/(E2-E1);
      FACT3 = HOMEGA2/(2.0*PI*TEMP-HOMEGA2);

      if(EE<E1){
      GFTUN = FACT1*(std::exp(EE/TEMP)*std::exp(2.0*PI*(EE-EF)/HOMEGA1)-std::exp(-2.0*PI*EF/HOMEGA1));
      }else{
        if(EE>=E1&&EE<E2){
         GFTUN = std::exp(EE/TEMP)*(0.50+FACT2*(EE-EF-TEMP))-std::exp(E1/TEMP)*(0.5+FACT2*(E1-EF-TEMP))+FACT1*(std::exp(E1/TEMP)*std::exp(2.0*PI*(E1-EF)/HOMEGA1)-std::exp(-2.0*PI*EF/HOMEGA1));
        }else{
         GFTUN = std::exp(EE/TEMP)*(1.0+FACT3*std::exp(-2.0*PI*(EE-EF)/HOMEGA2))-std::exp(E2/TEMP)*(1.0+FACT3*std::exp(-2.0*PI*(E2-EF)/HOMEGA2))+std::exp(E2/TEMP)*(0.5+FACT2*(E2-EF-TEMP))-std::exp(E1/TEMP)*(0.5+FACT2*(E1-EF-TEMP))+FACT1*(std::exp(E1/TEMP)*std::exp(2.0*PI*(E1-EF)/HOMEGA1)-std::exp(-2.0*PI*EF/HOMEGA1));
        }
      }
      GFTUN = GFTUN/std::exp(EE/TEMP)*DENSF*ENH_FACT/DENSG/2.0/PI;
      GFTUN = GFTUN * CORR_FUNCT;
 return GFTUN;
}


void G4Abla::fission_width(G4double ZPRF,G4double A,G4double EE,G4double BS,G4double BK,G4double EF,G4double Y,G4double *GF,G4double *TEMP,G4double JPR,G4int IEROT,G4int FF_ALLOWED,G4int OPTCOL,G4int OPTSHP,G4double DENSG)
{
//
 G4double FNORM,MASS_ASYM_SADD_B,FP_PER,FP_PAR,SIG_PER_SP,SIG_PAR_SP;
 G4double Z2OVERA,ftemp,fgf,DENSF,ECOR,EROT,qr;
 G4double DCR,UCR,ENH_FACTA,ENH_FACTB,ENH_FACT,PONFE;
 G4double PI = 3.14159;

      DCR = fiss->dcr;
      UCR = fiss->ucr;
      Z2OVERA = ZPRF * ZPRF / A;

// Nuclei below Businaro-Gallone point do not go through fission
      if((ZPRF<=55.0) || (FF_ALLOWED==0)){
        (*GF) = 0.0;
        (*TEMP) = 0.5;
        return;
      }

// Level density above SP
// Saddle-point deformation is defbet as above. But, FP_PER and FP_PAR
// are calculated for fission in DENSNIV acc to Myers and Hasse, and their
// parametrization is done as function of y
      densniv(A,ZPRF,EE,EF,&DENSF,0.0,BS,BK,&ftemp,OPTSHP,0,Y,&ECOR,JPR,1,&qr);

      if(OPTCOL==0){
         fgf= DENSF/DENSG/PI/2.0*ftemp;
         (*TEMP)=ftemp;
         (*GF)= fgf;
         return;
      }

// FP = 2/5*M0*R0**2/HBAR**2 * A**(5/3) * (1 + DEFBET/3)
// FP is used to calculate the spin-cutoff parameter SIG=FP*TEMP/hbar**2; hbar**2
// is, therefore, included in FP in order to avoid problems with large exponents
// The factor fnorm inlcudes then R0, M0 and hbar**2 -
// fnorm = R0*M0/hbar**2 = 1.2fm*931.49MeV/c**2 /(6.582122e-22 MeVs)**2 and is
// in units 1/MeV
      FNORM = 1.2*1.2 * 931.49 * 1.e-2 / (9.0 * 6.582122*6.582122);
// FP_PER ~ 1+7*y/6, FP_PAR ~ 1-7*y/3 (Hasse & Myers, Geom. relat. macr. nucl. phys.)
// Perpendicular moment of inertia
      FP_PER = 2.0/5.0*std::pow(A,5.0/3.0)*FNORM*(1. + 7.0/6.0*Y*(1.0+1396.0/255.*Y));

// AK - Jan 2011 - following line is needed, as for these nuclei it seems that
// FP_PER calculated according to above formula has too large values, leading to too
// large ENH_FACT
      if(Z2OVERA<=30.0) FP_PER = 6.50;

// Parallel moment of inertia
      FP_PAR = 2.0/5.0*std::pow(A,5.0/3.0)*FNORM*(1.0 - 7.0/3.0*Y*(1.0-389.0/255.0*Y));
      if(FP_PAR<0.0) FP_PAR = 0.0;

      EROT = JPR * JPR / (2.0 * std::sqrt(FP_PAR*FP_PAR + FP_PER*FP_PER));
      if(IEROT==1) EROT = 0.0;

// Perpendicular spin cut-off parameter
      SIG_PER_SP = std::sqrt(FP_PER * ftemp);

      if(SIG_PER_SP<1.0) SIG_PER_SP = 1.0;

// Parallel spin cut-off parameter
      SIG_PAR_SP = std::sqrt(FP_PAR * ftemp);
      ENH_FACT = 1.0;
//
        if(A>223.0){
         MASS_ASYM_SADD_B = 2.0;
        }else{
         MASS_ASYM_SADD_B = 1.0;
        }

// actinides with low barriers
      if(Z2OVERA>35.&&Z2OVERA<=(110.*110./298.0)){
// Barrier A is axial asymmetric
       ENH_FACTA = std::sqrt(8.0*PI) * SIG_PER_SP*SIG_PER_SP * SIG_PAR_SP;
// Barrier B is axial symmetric
      ENH_FACTB = MASS_ASYM_SADD_B * SIG_PER_SP*SIG_PER_SP;
// Total enhancement
      ENH_FACT = ENH_FACTA * ENH_FACTB / (ENH_FACTA + ENH_FACTB);
      }else{
// nuclei with high fission barriers (only barrier B plays a role, axial symmetric)
        if(Z2OVERA<=35.){
          ENH_FACT = MASS_ASYM_SADD_B*SIG_PER_SP*SIG_PER_SP;
         }else{
// super-heavy nuclei  (only barrier A plays a role, axial asymmetric)
          ENH_FACT =  std::sqrt(8.0*PI) * SIG_PER_SP*SIG_PER_SP* SIG_PAR_SP;
        }
      }

// Fading-out with excitation energy above the saddle point:
      PONFE = (ECOR-UCR-EROT)/DCR;
      if(PONFE>700.) PONFE = 700.0;
// Fading-out according to Junghans:
      ENH_FACT = 1.0/(1.0+std::exp(PONFE))*ENH_FACT+1.0;

      if(ENH_FACT<1.0)ENH_FACT = 1.0;
      fgf= DENSF/DENSG/PI/2.0*ftemp*ENH_FACT;

// Tunneling
      if(EE<EF+1.){
      fgf=tunnelling(A,ZPRF,Y,EE,EF,ftemp,DENSG,DENSF,ENH_FACT);
      }
//
      (*GF)= fgf;
      (*TEMP)=ftemp;
 return;
}


void G4Abla::lorb(G4double AMOTHER,G4double ADAUGHTER,G4double LMOTHER,G4double EEFINAL,G4double *LORBITAL,G4double *SIGMA_LORBITAL)
{

    G4double AFRAGMENT,S4FINAL,ALEVDENS;
    G4double THETA_MOTHER,THETA_ORBITAL;

/*
C     Values on input:
C       AMOTHER          mass of mother nucleus
C       ADAUGHTER        mass of daughter fragment
C       LMOTHER          angular momentum of mother (may be real)
C       EEFINAL          excitation energy after emission
C                          (sum of daughter and fragment)
C
C     Values on output:
C       LORBITAL         mean value of orbital angular momentum
C                           (assumed to be fully aligned with LMOTHER)
C       SIGMA_LORBITAL   standard deviation of the orbital angular momentum
*/
    if (EEFINAL<=0.01) EEFINAL = 0.01;
        AFRAGMENT = AMOTHER - ADAUGHTER;
        ALEVDENS = 0.073*AMOTHER + 0.095*std::pow(AMOTHER,2.0/3.0);
        S4FINAL = ALEVDENS * EEFINAL;
        if(S4FINAL <= 0.0 || S4FINAL > 100000.){
            std::cout<< "S4FINAL:" << S4FINAL << ALEVDENS << EEFINAL  << idnint(AMOTHER) << idnint(AFRAGMENT) << std::endl;
        }
        THETA_MOTHER = 0.0111 * std::pow(AMOTHER,1.66667);
        THETA_ORBITAL = 0.0323 / std::pow(AMOTHER,2.) *std::pow(std::pow(AFRAGMENT,0.33333) + std::pow(ADAUGHTER,0.33333),2.) * AFRAGMENT*ADAUGHTER*(AFRAGMENT+ADAUGHTER);

        *LORBITAL = -1.* THETA_ORBITAL * (LMOTHER / THETA_MOTHER + std::sqrt(S4FINAL) /(ALEVDENS*LMOTHER));

        *SIGMA_LORBITAL = std::sqrt(std::sqrt(S4FINAL) * THETA_ORBITAL / ALEVDENS);

 return;
}

// Random generator according to a distribution similar to a
// Maxwell distribution with quantum-mech. x-section for neutrons according to KHS
//      Y = SQRT(X) * EXP(-X/T) (approximation:)
G4double G4Abla::fvmaxhaz_neut(G4double x){

 return (2.0 * x * std::sqrt(std::log(G4AblaRandom::flat()) * std::log(G4AblaRandom::flat())));
}

void G4Abla::imf(G4double ACN,G4double ZCN,G4double TEMP,G4double EE,G4double *ZIMF,G4double *AIMF,G4double *BIMF,G4double *SBIMF,G4double *TIMF,G4double JPRF)
{
//     input variables (compound nucleus) Acn, Zcn, Temp, EE
//     output variable (IMF) Zimf,Aimf,Bimf,Sbimf,IRNDM
//
//     SBIMF = separation energy + coulomb barrier
//
//     SDW(Z) is the sum over all isotopes for a given Z of the decay widths
//     DW(Z,A) is the decay width of a certain nuclide
//
//  Last update:
//             28/10/13 - JLRS - from abrablav4 (AK)
//             13/11/16 - JLRS - Included this function in Abla++

   G4int IZIMFMAX=0;
   G4int iz=0,in=0,IZIMF=0,INMI=0,INMA=0,IZCN=0,INCN=0,INIMFMI=0,INIMFMA=0,ILIMMAX=0,INNMAX=0,INMIN=0,IAIMF=0,IZSTOP=3,IZMEM=0,IA=0,INMINMEM=0,INMAXMEM=0,IIA=0;
   G4double BS=0,BK=0,BC=0,BSHELL=0,DEFBET=0,DEFBETIMF=0,EROT=0,MAIMF=0,MAZ=0,MARES=0,AIMF_1,OMEGAP=0,fBIMF=0.0,BSIMF=0,A1PAR=0,A2PAR=0,SUM_A,EEDAUG;
   G4double DENSCN=0,TEMPCN=0,ECOR=0,IINERT=0,EROTCN=0,WIDTH_IMF=0.0,WIDTH1=0,IMFARG=0,QR=0,QRCN=0,DENSIMF=0,fTIMF=0,fZIMF=0,fAIMF=0.0,NIMF=0,fSBIMF=0;
   G4double PI = 3.141592653589793238;
   G4double ZIMF_1=0.0;
   G4double SDWprevious=0,SUMDW_TOT=0,SUM_Z=0,X=0,SUMDW_N_TOT=0,XX=0;
   G4double SDW[98];
   G4double DW[98][251];
   G4double BBIMF[98][251];
   G4double SSBIMF[98][251];
   G4int OPTSHPIMF=opt->optshpimf;
     
   // Initialization
   for (G4int ia = 0; ia < 98; ia++)
    for (G4int ib = 0; ib < 251; ib++) {
      BBIMF[ia][ib] = 0.0;
      SSBIMF[ia][ib] = 0.0;
    }

   // take the half of the CN and transform it in integer (floor it)
   IZIMFMAX = idnint(ZCN / 2.0);

   if(IZIMFMAX<3){
         std::cout << "CHARGE_IMF line 46" << std::endl;
         std::cout << "Problem: IZIMFMAX < 3 " << std::endl;
         std::cout << "ZCN,IZIMFMAX," << ZCN << "," << IZIMFMAX << std::endl;
   }

  iz = idnint(ZCN);
  in = idnint(ACN) - iz;
  BSHELL = ecld->ecgnz[in][iz]- ecld->vgsld[in][iz];
  DEFBET = ecld->beta2[in][iz];

  bsbkbc(ACN,ZCN,&BS,&BK,&BC);

  densniv(ACN,ZCN,EE,0.0,&DENSCN,BSHELL,BS,BK,&TEMPCN,0,0,DEFBET,&ECOR,JPRF,0,&QRCN);

  IINERT = 0.4 * 931.49 * 1.16*1.16 * std::pow(ACN,5.0/3.0)*(1.0 + 0.5*std::sqrt(5./(4.*PI))*DEFBET);
  EROTCN = JPRF * JPRF * 197.328 * 197.328 /(2. * IINERT);
//
  for(IZIMF=3;IZIMF<=IZIMFMAX;IZIMF++){
  
     SDW[IZIMF] = 0.0;
     ZIMF_1 = 1.0*IZIMF;

//     *** Find the limits that both IMF and partner are bound :

     isostab_lim(IZIMF,&INIMFMI,&INIMFMA);// Bound isotopes for IZIMF from INMIN to INIMFMA
// Idea - very proton-rich nuclei can live long enough to evaporate IMF before decaying:
     INIMFMI = max(1,INIMFMI-2);

     IZCN = idnint(ZCN);        //  Z of CN
     INCN = idnint(ACN) - IZCN; //  N of CN

     isostab_lim(IZCN-IZIMF,&INMI,&INMA); // Daughter nucleus after IMF emission,
                                         // limits of bound isotopes
     INMI = max(1,INMI-2);
     INMIN = max(INIMFMI,INCN-INMA);  //  Both IMF and daughter must be bound
     INNMAX = min(INIMFMA,INCN-INMI); //   "

     ILIMMAX = max(INNMAX,INMIN);     // In order to keep the variables below
//     ***

     for(G4int INIMF=INMIN;INIMF<=ILIMMAX;INIMF++){ // Range of possible IMF isotopes
          IAIMF = IZIMF + INIMF;
          DW[IZIMF][IAIMF] = 0.0;
          AIMF_1 = 1.0*(IAIMF);

//         Q-values
          mglms(ACN-AIMF_1,ZCN-ZIMF_1,OPTSHPIMF,&MARES);
          mglms(AIMF_1,ZIMF_1,OPTSHPIMF,&MAIMF);
          mglms(ACN,ZCN,OPTSHPIMF,&MAZ);

//         Barrier
          if(ACN<=AIMF_1){
            SSBIMF[IZIMF][IAIMF] = 1.e37;
          }else{
            barrs(idnint(ZCN-ZIMF_1),idnint(ACN-AIMF_1),idnint(ZIMF_1),idnint(AIMF_1),&fBIMF,&OMEGAP);
            SSBIMF[IZIMF][IAIMF] = MAIMF + MARES - MAZ + fBIMF;
            BBIMF[IZIMF][IAIMF] = fBIMF;
          }

// *****  Width *********************
          DEFBETIMF = ecld->beta2[idnint(AIMF_1-ZIMF_1)][idnint(ZIMF_1)]+ecld->beta2[idnint(ACN-AIMF_1-ZCN+ZIMF_1)][idnint(ZCN-ZIMF_1)];

          IINERT = 0.40 * 931.490 * 1.160*1.160 * std::pow(ACN,5.0/3.0)*(std::pow(AIMF_1,5.0/3.0) + std::pow(ACN - AIMF_1,5.0/3.0)) + 931.490 * 1.160*1.160 * AIMF_1 * (ACN-AIMF_1) / ACN *(std::pow(AIMF_1,1.0/3.0) + std::pow(ACN - AIMF_1,1.0/3.0))*(std::pow(AIMF_1,1.0/3.0) + std::pow(ACN - AIMF_1,1.0/3.0));

          EROT = JPRF * JPRF * 197.328 * 197.328 /(2.0 * IINERT);

 //      IF(IEROT.EQ.1) EROT = 0.D0
          if (EE<(SSBIMF[IZIMF][IAIMF]+EROT) || DENSCN<=0.0){
           WIDTH_IMF = 0.0;
//          PRINT*,IDNINT(ACN),IDNINT(ZCN),IZIMF,IAIMF
          }else{
//          here the temperature at "saddle point" is used
// Increase of the level densitiy at the barrier due to deformation; see comment in ABLA
//          BSIMF = ((ACN-AIMF_1)**(2.D0/3.D0) + AIMF_1**(2.D0/3.D0))/
//     &                ACN**(2.D0/3.D0)
           BSIMF = BS;
           densniv(ACN,ZCN,EE,SSBIMF[IZIMF][IAIMF],&DENSIMF,0.0,BSIMF,1.0,&fTIMF,0,0,DEFBETIMF,&ECOR,JPRF,2,&QR);
           IMFARG = (SSBIMF[IZIMF][IAIMF]+EROTCN-EROT)/fTIMF;
           if(IMFARG>200.0) IMFARG = 200.0;

           WIDTH1 = width(ACN,ZCN,AIMF_1,ZIMF_1,fTIMF,fBIMF,SSBIMF[IZIMF][IAIMF],EE-EROT);

           WIDTH_IMF = WIDTH1 * std::exp(-IMFARG) * QR / QRCN;

           if(WIDTH_IMF<=0.0){
            std::cout << "GAMMA_IMF=0 -> LOOK IN GAMMA_IMF CALCULATIONS!" << std::endl;
            std::cout << "ACN,ZCN,AIMF,ZIMF:" << idnint(ACN) << "," << idnint(ZCN) << "," << idnint(AIMF_1) << "," << idnint(ZIMF_1) << std::endl;
            std::cout << "SSBIMF,TIMF :" << SSBIMF[IZIMF][IAIMF] << "," << fTIMF << std::endl;
            std::cout << "DEXP(-IMFARG) = " << std::exp(-IMFARG) << std::endl;
            std::cout << "WIDTH1 =" << WIDTH1 << std::endl;
           }
          }// if ee

          SDW[IZIMF] = SDW[IZIMF] + WIDTH_IMF;

          DW[IZIMF][IAIMF] = WIDTH_IMF;

     }// for INIMF
  }// for IZIMF
//     End loop to calculate the decay widths ************************
//     ***************************************************************

//     Loop to calculate where the gamma of IMF has the minimum ******
      SDWprevious = 1.e20;
      IZSTOP = 0;

      for(G4int III_ZIMF=3;III_ZIMF<=IZIMFMAX;III_ZIMF++){

        if(SDW[III_ZIMF]==0.0){
          IZSTOP = III_ZIMF - 1;
            goto imfs30;
        }

        if(SDW[III_ZIMF]>SDWprevious){
            IZSTOP = III_ZIMF - 1;
            goto imfs30;
        }else{
            SDWprevious = SDW[III_ZIMF];
        }

      }// for III_ZIMF

      imfs30:

      if(IZSTOP<=6){
       IZSTOP = IZIMFMAX;
       goto imfs15;
      }

      A1PAR = std::log10(SDW[IZSTOP]/SDW[IZSTOP-2])/std::log10((1.0*IZSTOP)/(1.0*IZSTOP-2.0));
      A2PAR = std::log10(SDW[IZSTOP]) - A1PAR * std::log10(1.0*(IZSTOP));
      if(A2PAR>0.)A2PAR=-1.*A2PAR;
      if(A1PAR>0.)A1PAR=-1.*A1PAR;

//     End loop to calculate where gamma of IMF has the minimum

      for(G4int II_ZIMF = IZSTOP;II_ZIMF<=IZIMFMAX;II_ZIMF++){
       SDW[II_ZIMF] =  std::pow(10.0,A2PAR) * std::pow(1.0*II_ZIMF,A1PAR);  // Power-low
       if(SDW[II_ZIMF]<0.0) SDW[II_ZIMF] = 0.0;
      }

      imfs15:

//    Sum of all decay widths (for normalisation)
      SUMDW_TOT = 0.0;
      for(G4int I_ZIMF = 3;I_ZIMF<=IZIMFMAX;I_ZIMF++){
        SUMDW_TOT = SUMDW_TOT + SDW[I_ZIMF];
      }
      if(SUMDW_TOT<=0.0){
        std::cout << "*********************" << std::endl;
        std::cout <<  "IMF function" << std::endl;
        std::cout <<  "SUM of decay widths = " << SUMDW_TOT << " IZIMFMAX = " << IZIMFMAX << std::endl;
        std::cout <<  "IZSTOP = " << IZSTOP << std::endl;
      }

//    End of Sum of all decay widths (for normalisation)
  
//    Loop to sample the nuclide that is emitted ********************
//    ------- sample Z -----------
      imfs10:
      X = haz(1)*SUMDW_TOT;

//      IF(X.EQ.0.D0) PRINT*,'WARNING: X=0',XRNDM,SUMDW_TOT
      SUM_Z = 0.0;
      fZIMF = 0.0;
      IZMEM = 0;

      for(G4int IZ = 3;IZ<=IZIMFMAX;IZ++){
         SUM_Z = SUM_Z + SDW[IZ];
         if(X<SUM_Z){
            fZIMF = 1.0*IZ;
            IZMEM = IZ;
            goto imfs20;
         }
      }//for IZ

      imfs20:

//     ------- sample N -----------

      isostab_lim(IZMEM,&INMINMEM,&INMAXMEM);
      INMINMEM = max(1,INMINMEM-2);

      isostab_lim(IZCN-IZMEM,&INMI,&INMA);  // Daughter nucleus after IMF emission,
      INMI = max(1,INMI-2);
                                            // limits of bound isotopes

      INMINMEM = max(INMINMEM,INCN-INMA); // Both IMF and daughter must be bound
      INMAXMEM = min(INMAXMEM,INCN-INMI); //   "

      INMAXMEM = max(INMINMEM,INMAXMEM);

      IA = 0;
      SUMDW_N_TOT = 0.0;
      for(G4int IIINIMF = INMINMEM;IIINIMF<=INMAXMEM;IIINIMF++){
       IA = IZMEM + IIINIMF;
       if(IZMEM>=3&&IZMEM<=95&&IA>=4&&IA<=250){
        SUMDW_N_TOT = SUMDW_N_TOT + DW[IZMEM][IA];
       }else{
         std::cout << "CHARGE IMF OUT OF RANGE" << IZMEM << ", " << IA << ", " << idnint(ACN) << ", " << idnint(ZCN) << ", " << TEMP << std::endl;
       }
      }

      XX = haz(1)*SUMDW_N_TOT;
      IIA = 0;
      SUM_A = 0.0;
      for(G4int IINIMF = INMINMEM;IINIMF<=INMAXMEM; IINIMF++){
        IIA = IZMEM + IINIMF;
  //      SUM_A = SUM_A + DW[IZ][IIA]; //FIXME
        SUM_A = SUM_A + DW[IZMEM][IIA];
        if(XX<SUM_A){
          fAIMF = G4double(IIA);
          goto imfs25;
        }
      }

      imfs25:
//     CHECK POINT 1
      NIMF = fAIMF - fZIMF;

      if((ACN-ZCN-NIMF)<=0.0 || (ZCN-fZIMF) <= 0.0){
       std::cout << "IMF Partner unstable:" << std::endl;
       std::cout << "System: Acn,Zcn,NCN:" << std::endl;
       std::cout << idnint(ACN) << ", " << idnint(ZCN) << ", " << idnint(ACN-ZCN) << std::endl;
       std::cout << "IMF: A,Z,N:" << std::endl;
       std::cout << idnint(fAIMF) << ", " << idnint(fZIMF) << ", " << idnint(fAIMF-fZIMF) << std::endl;
       std::cout << "Partner: A,Z,N:" << std::endl;
       std::cout << idnint(ACN-fAIMF) << ", " << idnint(ZCN-fZIMF) << ", " << idnint(ACN-ZCN-NIMF) << std::endl;
        std::cout << "----nmin,nmax" << INMINMEM << ", " << INMAXMEM << std::endl;
        std::cout << "----- warning: Zimf=" << fZIMF << " Aimf=" << fAIMF << std::endl;
        std::cout << "----- look in subroutine IMF" << std::endl;
        std::cout << "ACN,ZCN,ZIMF,AIMF,temp,EE,JPRF::" << ACN << ", " << ZCN << ", " << fZIMF << ", " << fAIMF << ", " << TEMP << ", " << EE << ", " << JPRF << std::endl;
std::cout << "-IZSTOP,IZIMFMAX:" << IZSTOP << ", " << IZIMFMAX << std::endl;
std::cout << "----X,SUM_Z,SUMDW_TOT:" << X << ", " << SUM_Z << ", " << SUMDW_TOT << std::endl;
//for(int III_ZIMF=3;III_ZIMF<=IZIMFMAX;III_ZIMF++)
  //    std::cout << "-**Z,SDW:" << III_ZIMF << ", " << SDW[III_ZIMF] << std::endl;

      goto imfs10;
      }
      if(fZIMF>=ZCN || fAIMF>=ACN || fZIMF<=2 || fAIMF<=3){
        std::cout << "----nmin,nmax" << INMINMEM << ", " << INMAXMEM << std::endl;
        std::cout << "----- warning: Zimf=" << fZIMF << " Aimf=" << fAIMF << std::endl;
        std::cout << "----- look in subroutine IMF" << std::endl;
        std::cout << "ACN,ZCN,ZIMF,AIMF,temp,EE,JPRF:" << ACN << ", " << ZCN << ", " << fZIMF << ", " << fAIMF << ", " << TEMP << ", " << EE << ", " << JPRF << std::endl;
std::cout << "-IZSTOP,IZIMFMAX:" << IZSTOP << ", " << IZIMFMAX << std::endl;
std::cout << "----X,SUM_Z,SUMDW_TOT:" << X << ", " << SUM_Z << ", " << SUMDW_TOT << std::endl;
for(int III_ZIMF=3;III_ZIMF<=IZIMFMAX;III_ZIMF++)
      std::cout << "-**Z,SDW:" << III_ZIMF << ", " << SDW[III_ZIMF] << std::endl;

        fZIMF = 3.0;  // provisorisch AK
        fAIMF = 4.0;
      }

// Characteristics of selected IMF (AIMF, ZIMF, BIMF, SBIMF, TIMF)
      fSBIMF = SSBIMF[idnint(fZIMF)][idnint(fAIMF)];
      fBIMF = BBIMF[idnint(fZIMF)][idnint(fAIMF)];

      if((ZCN-fZIMF)<=0.0)std::cout << "CHARGE_IMF ZIMF > ZCN" << std::endl;
      if((ACN-fAIMF)<=0.0)std::cout << "CHARGE_IMF AIMF > ACN" << std::endl;

      BSHELL = ecld->ecgnz[idnint(ACN-ZCN-NIMF)][idnint(ZCN-fZIMF)] -ecld->vgsld[idnint(ACN-ZCN-NIMF)][idnint(ZCN-fZIMF)];

      DEFBET = ecld->beta2[idnint(ACN-ZCN-NIMF)][idnint(ZCN-fZIMF)];
      EEDAUG = (EE - fSBIMF) * (ACN - fAIMF) / ACN;
      bsbkbc(ACN - fAIMF,ZCN-fZIMF,&BS,&BK,&BC);
      densniv(ACN-fAIMF,ZCN-fZIMF,EEDAUG,0.0,&DENSIMF,BSHELL,BS,BK,&fTIMF,0,0,DEFBET,&ECOR,0.0,0,&QR);

      if(fSBIMF>EE){
        std::cout << "----- warning: EE=" << EE << "," << " S+Bimf=" << fSBIMF  << std::endl;
        std::cout << "----- look in subroutine IMF" << std::endl;
        std::cout << "IMF will be resampled" << std::endl;
       goto imfs10;
      }
   (*ZIMF) = fZIMF;
   (*AIMF) = fAIMF;
   (*SBIMF) = fSBIMF;
   (*BIMF) = fBIMF;
   (*TIMF) = fTIMF;
  return;
}

void G4Abla::isostab_lim(G4int z, G4int *nmin, G4int *nmax)
{

G4int VISOSTAB[191][2]={
           {0 ,           7 },
           {1 ,           8 },
           {1 ,           9 },
           {2 ,          12 },
           {2 ,          14 },
           {2 ,          16 },
           {3 ,          18 },
           {4 ,          22 },
           {6 ,          22 },
           {6 ,          28 },
           {7 ,          28 },
           {7 ,          30 },
           {8 ,          28 },
           {8 ,          36 },
          {10 ,          38 },
          {10 ,          40 },
          {11 ,          38 },
          {10 ,          42 },
          {13 ,          50 },
          {14 ,          50 },
          {15 ,          52 },
          {16 ,          52 },
          {17 ,          54 },
          {18 ,          54 },
          {19 ,          60 },
          {19 ,          62 },
          {21 ,          64 },
          {20 ,          66 },
          {23 ,          66 },
          {24 ,          70 },
          {25 ,          70 },
          {26 ,          74 },
          {27 ,          78 },
          {29 ,          82 },
          {33 ,          82 },
          {31 ,          82 },
          {35 ,          82 },
          {34 ,          84 },
          {40 ,          84 },
          {36 ,          86 },
          {40 ,          92 },
          {38 ,          96 },
          {42 ,         102 },
          {42 ,         102 },
          {44 ,         102 },
          {42 ,         106 },
          {47 ,         112 },
          {44 ,         114 },
          {49 ,         116 },
          {46 ,         118 },
          {52 ,         120 },
          {52 ,         124 },
          {55 ,         126 },
          {54 ,         126 },
          {57 ,         126 },
          {57 ,         126 },
          {60 ,         126 },
          {58 ,         130 },
        {  62 ,         132 },
        {  60 ,         140 },
        {  67 ,         138 },
        {  64 ,         142 },
        {  67 ,         144 },
        {  68 ,         146 },
        {  70 ,         148 },
        {  70 ,         152 },
        {  73 ,         152 },
        {  72 ,         154 },
        {  75 ,         156 },
        {  77 ,         162 },
        {  79 ,         164 },
        {  78 ,         164 },
        {  82 ,         166 },
        {  80 ,         166 },
        {  85 ,         168 },
        {  83 ,         176 },
        {  87 ,         178 },
        {  88 ,         178 },
        {  91 ,         182 },
        {  90 ,         184 },
        {  96 ,         184 },
        {  95 ,         184 },
        {  99 ,         184 },
        {  98 ,         184 },
        { 105 ,         194 },
        { 102 ,         194 },
        { 108 ,         196 },
        { 106 ,         198 },
        { 115 ,         204 },
        { 110 ,         206 },
        { 119 ,         210 },
        { 114 ,         210 },
        { 124 ,         210 },
        { 117 ,         212 },
        { 130 ,         212 }
        };

      if (z<0){
        *nmin = 0;
        *nmax = 0;
      }else{ 
        if(z==0){
         *nmin = 1;
         *nmax = 1;
// AK (Dez2010) - Just to avoid numerical problems
        }else{ 
         if(z>95){
          *nmin = 130;
          *nmax = 200;
         }else{
          *nmin = VISOSTAB[z-1][0];
          *nmax = VISOSTAB[z-1][1];
         }
        }
      }

  return;
}


void G4Abla::evap_postsaddle(G4double A, G4double Z, G4double EXC, G4double *E_scission_post, G4double *A_scission, G4double *Z_scission,G4double &vx_eva,G4double &vy_eva,G4double &vz_eva,G4int *NbLam0_par){

//  AK 2006 - Now in case of fission deexcitation between saddle and scission
//            is explicitly calculated. Langevin calculations made by P. Nadtochy
//            used to parametrise saddle-to-scission time

  G4double af,zf,ee;
  G4double epsiln = 0.0, probp = 0.0, probd = 0.0, probt = 0.0, probn = 0.0, probhe = 0.0, proba = 0.0, probg = 0.0, probimf=0.0, problamb0=0.0, ptotl = 0.0, tcn = 0.0;  
  G4double sn = 0.0, sbp = 0.0, sbd = 0.0, sbt = 0.0, sbhe = 0.0, sba = 0.0, x = 0.0, amoins = 0.0, zmoins = 0.0,sp= 0.0,sd= 0.0,st= 0.0,she= 0.0,sa= 0.0, slamb0 = 0.0;
  G4double ecn = 0.0, ecp = 0.0, ecd = 0.0, ect = 0.0,eche = 0.0,eca = 0.0, ecg = 0.0, eclamb0 = 0.0, bp = 0.0, bd = 0.0, bt = 0.0, bhe = 0.0, ba = 0.0;

  G4double xcv=0.,ycv=0.,zcv=0.,VXOUT=0.,VYOUT=0.,VZOUT=0.;

  G4double jprfn=0.0, jprfp=0.0, jprfd=0.0, jprft=0.0, jprfhe=0.0, jprfa=0.0, jprflamb0=0.0;
  G4double ctet1 = 0.0, stet1 = 0.0, phi1 = 0.0;
  G4double rnd = 0.0;

  G4int itest = 0, sortie=0;
  G4double probf = 0.0;

  G4double ef = 0.0;
  G4double pc = 0.0;

  G4double time,tauf,tau0,a0,a1,emin,ts1,tsum=0.;
  G4int inttype=0,inum=0,gammadecay = 0, flamb0decay = 0;
  G4double pleva = 0.0;
  G4double pxeva = 0.0;
  G4double pyeva = 0.0;
  G4double pteva = 0.0;
  G4double etot = 0.0;
  G4int NbLam0= (*NbLam0_par);

  const G4double c = 29.9792458;
  const G4double mu = 931.494;
  const G4double mu2 = 931.494*931.494;

      vx_eva=0.;
      vy_eva=0.;
      vz_eva=0.;
      IEV_TAB_SSC = 0;


      af = dint(A);
      zf = dint(Z);
      ee = EXC;

      fiss->ifis = 0;
      opt->optimfallowed = 0;
      gammaemission=0;
// Initialsation
      time = 0.0;
 
// in sec
      tau0 = 1.0e-21;
      a0 = 0.66482503 - 3.4678935 * std::exp(-0.0104002*ee);
      a1  = 5.6846e-04 + 0.00574515 * std::exp(-0.01114307*ee);
      tauf = (a0 + a1 * zf*zf/std::pow(af,0.3333333)) * tau0;
//
      post10:
      direct(zf,af,ee,0.,&probp,&probd,&probt,&probn,&probhe,&proba,&probg,&probimf,&probf,&problamb0,&ptotl,
	 &sn,&sbp,&sbd,&sbt,&sbhe,&sba,&slamb0,
         &ecn,&ecp,&ecd,&ect,&eche,&eca,&ecg,&eclamb0,
         &bp,&bd,&bt,&bhe,&ba,&sp,&sd,&st,&she,&sa,&ef,&ts1,inttype,inum,itest,&sortie,&tcn,
         &jprfn, &jprfp, &jprfd, &jprft, &jprfhe, &jprfa, &jprflamb0, &tsum, NbLam0); //:::FIXME::: Call
//
// HERE THE FINAL STEPS OF THE EVAPORATION ARE CALCULATED
//
  if(ptotl<=0.)goto post100;

  emin = dmin1(sba,sbhe,dmin1(sbt,sbhe,dmin1(sn,sbp,sbd)));

  if(emin>1e30)std::cout << "ERROR AT THE EXIT OF EVAPORA,E>1.D30,AF" << std::endl;

  if(sortie==1){
   if (probn!=0.0) {
    amoins = 1.0;
    zmoins = 0.0;
    epsiln = sn + ecn;
    pc = std::sqrt(std::pow((1.0 + ecn/9.3956e2),2.) - 1.0) * 9.3956e2;
    gammadecay = 0;
    flamb0decay = 0;
   }
   else if(probp!=0.0){
    amoins = 1.0;
    zmoins = 1.0;
    epsiln = sp + ecp;
    pc = std::sqrt(std::pow((1.0 + ecp/9.3827e2),2.) - 1.0) * 9.3827e2;
    gammadecay = 0;
    flamb0decay = 0;
   }
   else if(probd!=0.0){
    amoins = 2.0;
    zmoins = 1.0;
    epsiln = sd + ecd;
    pc = std::sqrt(std::pow((1.0 + ecd/1.875358e3),2) - 1.0) * 1.875358e3;
    gammadecay = 0;
    flamb0decay = 0;
   }
   else if(probt!=0.0){
    amoins = 3.0;
    zmoins = 1.0;
    epsiln = st + ect;
    pc = std::sqrt(std::pow((1.0 + ect/2.80828e3),2) - 1.0) * 2.80828e3;
    gammadecay = 0;
    flamb0decay = 0;
   }
   else if(probhe!=0.0){
    amoins = 3.0;
    zmoins = 2.0;
    epsiln = she + eche;
    pc = std::sqrt(std::pow((1.0 + eche/2.80826e3),2) - 1.0) * 2.80826e3;
    gammadecay = 0;
    flamb0decay = 0;
   }
   else{ if(proba!=0.0){
    amoins = 4.0;
    zmoins = 2.0;
    epsiln = sa + eca;
    pc = std::sqrt(std::pow((1.0 + eca/3.72834e3),2) - 1.0) * 3.72834e3;
    gammadecay = 0;
    flamb0decay = 0;
    }
   }
  goto post99;
  }

  //    IRNDM = IRNDM+1;
//
// HERE THE NORMAL EVAPORATION CASCADE STARTS
// RANDOM NUMBER FOR THE EVAPORATION


  // random number for the evaporation
  x = G4AblaRandom::flat() * ptotl;

  itest = 0;
  if (x < proba) {
    // alpha evaporation                                                     
    amoins = 4.0;
    zmoins = 2.0;
    epsiln = sa + eca;
    pc = std::sqrt(std::pow((1.0 + eca/3.72834e3),2) - 1.0) * 3.72834e3;
    gammadecay = 0;
    flamb0decay = 0;
  }  
  else if (x < proba+probhe) {
    // He3 evaporation                                                    
    amoins = 3.0;
    zmoins = 2.0;
    epsiln = she + eche;
    pc = std::sqrt(std::pow((1.0 + eche/2.80826e3),2) - 1.0) * 2.80826e3;
    gammadecay = 0;
    flamb0decay = 0;
  }
  else if (x < proba+probhe+probt) {
    // triton evaporation                                                    
    amoins = 3.0;
    zmoins = 1.0;
    epsiln = st + ect;
    pc = std::sqrt(std::pow((1.0 + ect/2.80828e3),2) - 1.0) * 2.80828e3;
    gammadecay = 0;
    flamb0decay = 0;
  }
  else if (x < proba+probhe+probt+probd) {
    // deuteron evaporation                                                    
    amoins = 2.0;
    zmoins = 1.0;
    epsiln = sd + ecd;
    pc = std::sqrt(std::pow((1.0 + ecd/1.875358e3),2) - 1.0) * 1.875358e3;
    gammadecay = 0;
    flamb0decay = 0;
  }
  else if (x < proba+probhe+probt+probd+probp) {
    // proton evaporation                                                    
    amoins = 1.0;
    zmoins = 1.0;
    epsiln = sp + ecp;
    pc = std::sqrt(std::pow((1.0 + ecp/9.3827e2),2) - 1.0) * 9.3827e2;
    gammadecay = 0;
    flamb0decay = 0;
  }
  else if (x < proba+probhe+probt+probd+probp+probn) {
    // neutron evaporation                                                   
    amoins = 1.0;
    zmoins = 0.0;
    epsiln = sn + ecn;
    pc = std::sqrt(std::pow((1.0 + ecn/9.3956e2),2.) - 1.0) * 9.3956e2;
    gammadecay = 0;
    flamb0decay = 0;
  }
  else if (x < proba+probhe+probt+probd+probp+probn+problamb0) {
    // lambda0 evaporation  
    amoins = 1.0;
    zmoins = 0.0;
    epsiln = slamb0 + eclamb0;
    pc = std::sqrt(std::pow((1.0 + (eclamb0)/11.1568e2),2.) - 1.0) * 11.1568e2;
    opt->nblan0 = opt->nblan0 -1;
    NbLam0 = NbLam0 -1;
    gammadecay = 0;
    flamb0decay = 1;
  }
  else if (x < proba+probhe+probt+probd+probp+probn+problamb0+probg) {
    // gamma evaporation                                                    
    amoins = 0.0;
    zmoins = 0.0;
    epsiln = ecg;
    pc = ecg;
    gammadecay = 1;
    flamb0decay = 0;
    if(probp==0.0 && probn==0.0 && probd==0.0 && probt==0.0 && proba==0.0 && probhe==0.0 && problamb0==0.0 && probimf==0.0 && probf==0.0){
    //ee = ee-epsiln;
    //if(ee<=0.01) ee = 0.010;
    goto post100;
    }
  }

// CALCULATION OF THE DAUGHTER NUCLEUS
//
      post99:

      if(gammadecay==1 && ee<=0.01+epsiln){
       epsiln = ee-0.01;
       time = tauf + 1.;
      }

      af = af-amoins;
      zf = zf-zmoins;
      ee = ee-epsiln;

      if(ee<=0.01) ee = 0.010;

      if(af<2.5) goto post100;

      time = time + ts1;

// Determination of x,y,z components of momentum from known emission momentum
        if(flamb0decay==1){
        EV_TAB_SSC[IEV_TAB_SSC][0] = 0.;
        EV_TAB_SSC[IEV_TAB_SSC][1] = -2.;
        EV_TAB_SSC[IEV_TAB_SSC][5] = 1.;
        }else{
        EV_TAB_SSC[IEV_TAB_SSC][0] = zmoins;
        EV_TAB_SSC[IEV_TAB_SSC][1] = amoins;
        EV_TAB_SSC[IEV_TAB_SSC][5] = 0.;
        }

        rnd = G4AblaRandom::flat();
        ctet1 = 2.0*rnd - 1.0;           // z component: uniform probability between -1 and 1
        stet1 = std::sqrt(1.0 - std::pow(ctet1,2));// component perpendicular to z
        rnd = G4AblaRandom::flat();
        phi1 = rnd*2.0*3.141592654;   // angle in x-y plane: uniform probability between 0 and 2*pi
        xcv = stet1*std::cos(phi1);   // x component
        ycv = stet1*std::sin(phi1);   // y component
        zcv = ctet1;                  // z component
// In the CM system
        if(gammadecay==0){
// Light particle
           G4double ETOT_LP = std::sqrt(pc*pc + amoins*amoins * mu2);
           if(flamb0decay==1)ETOT_LP = std::sqrt(pc*pc + 1115.683*1115.683);
           EV_TAB_SSC[IEV_TAB_SSC][2] = c * pc * xcv / ETOT_LP;
           EV_TAB_SSC[IEV_TAB_SSC][3] = c * pc * ycv / ETOT_LP;
           EV_TAB_SSC[IEV_TAB_SSC][4] = c * pc * zcv / ETOT_LP;
        }else{
// gamma ray
           EV_TAB_SSC[IEV_TAB_SSC][2] = pc * xcv;
           EV_TAB_SSC[IEV_TAB_SSC][3] = pc * ycv;
           EV_TAB_SSC[IEV_TAB_SSC][4] = pc * zcv;
        }
        lorentz_boost(vx_eva,vy_eva,vz_eva,
            EV_TAB_SSC[IEV_TAB_SSC][2],EV_TAB_SSC[IEV_TAB_SSC][3],
            EV_TAB_SSC[IEV_TAB_SSC][4],
            &VXOUT,&VYOUT,&VZOUT);
        EV_TAB_SSC[IEV_TAB_SSC][2] = VXOUT;
        EV_TAB_SSC[IEV_TAB_SSC][3] = VYOUT;
        EV_TAB_SSC[IEV_TAB_SSC][4] = VZOUT;

// Heavy residue
        if(gammadecay==0){
        G4double v2 = std::pow(EV_TAB_SSC[IEV_TAB_SSC][2],2.) +
             std::pow(EV_TAB_SSC[IEV_TAB_SSC][3],2.) +
             std::pow(EV_TAB_SSC[IEV_TAB_SSC][4],2.);
        G4double gamma = 1.0/std::sqrt(1.0 - v2 / (c*c));
        G4double etot_lp = amoins*mu * gamma;
        pxeva = pxeva - EV_TAB_SSC[IEV_TAB_SSC][2] * etot_lp / c;
        pyeva = pyeva - EV_TAB_SSC[IEV_TAB_SSC][3] * etot_lp / c;
        pleva = pleva - EV_TAB_SSC[IEV_TAB_SSC][4] * etot_lp / c;
        }else{
// in case of gammas, EV_TEMP contains momentum components and not velocity
        pxeva = pxeva - EV_TAB_SSC[IEV_TAB_SSC][2];
        pyeva = pyeva - EV_TAB_SSC[IEV_TAB_SSC][3];
        pleva = pleva - EV_TAB_SSC[IEV_TAB_SSC][4];
        }
        pteva = std::sqrt(pxeva*pxeva + pyeva*pyeva);
// To be checked:
        etot = std::sqrt ( pleva*pleva + pteva*pteva + af*af * mu2 );
        vx_eva = c * pxeva / etot;  // recoil velocity components of residue due to evaporation
        vy_eva = c * pyeva / etot;
        vz_eva = c * pleva / etot;

        IEV_TAB_SSC = IEV_TAB_SSC +1;

      if(time<tauf)goto post10;
//
     post100:
//
      *A_scission= af;
      *Z_scission= zf;
      *E_scission_post = ee;
      *NbLam0_par = NbLam0;
 return;
}

G4double G4Abla::getdeltabinding(G4double A, G4int H){
 if(A<1.)return (1.*H)/A*(10.68*A-21.27*std::pow(A,2./3.))*10.;
 return (1.*H)/A*(10.68*A-21.27*std::pow(A,2./3.));
}


G4double G4Abla::gethyperseparation(G4double A, G4double Z, G4int ny){
 if(A<1.)return 1.e38;
// For light nuclei we take experimental values
// Journal of Physics G, Nucl Part Phys 32,363 (2006)
  if (ny == 1) {
    if (Z == 1 && A == 4)
      return 2.04;
    else if (Z == 2 && A == 4)
      return 2.39;
    else if (Z == 2 && A == 5)
      return 3.12;
    else if (Z == 2 && A == 6)
      return 4.18;
    else if (Z == 2 && A == 7)
      return 5.23;
    else if (Z == 2 && A == 8)
      return 7.16;
    else if (Z == 3 && A == 6)
      return 4.50;
    else if (Z == 3 && A == 7)
      return 5.58;
    else if (Z == 3 && A == 8)
      return 6.80;
    else if (Z == 3 && A == 9)
      return 8.50;
    else if (Z == 4 && A == 7)
      return 5.16;
    else if (Z == 4 && A == 8)
      return 6.84;
    else if (Z == 4 && A == 9)
      return 6.71;
    else if (Z == 4 && A == 10)
      return 9.11;
    else if (Z == 5 && A == 9)
      return 8.29;
    else if (Z == 5 && A == 10)
      return 9.01;
    else if (Z == 5 && A == 11)
      return 10.29;
    else if (Z == 5 && A == 12)
      return 11.43;
    else if (Z == 6 && A == 12)
      return 10.95;
    else if (Z == 6 && A == 13)
      return 11.81;
    else if (Z == 6 && A == 14)
      return 12.50;
    else if (Z == 7 && A == 14)
      return 12.17;
    else if (Z == 7 && A == 15)
      return 13.59;
    else if (Z == 8 && A == 16)
      return 12.50;
    else if (Z == 8 && A == 17)
      return 13.59;
    else if (Z == 14 && A == 28)
      return 16.0;
    else if (Z == 39 && A == 89)
      return 22.1;
    else if (Z == 57 && A == 139)
      return 23.8;
    else if (Z == 82 && A == 208)
      return 26.5;
 }//ny==1
// For other nuclei we take Bethe-Weizsacker mass formula
 return gethyperbinding(A, Z, ny)-gethyperbinding(A-1., Z, ny-1);
}

G4double G4Abla::gethyperbinding(G4double A, G4double Z, G4int ny){
//
// Bethe-Weizsacker mass formula
// Journal of Physics G, Nucl Part Phys 32,363 (2006)
//
 if(A<2 || Z<2)return 0.;
 G4double N = A-Z -1.*ny;
 G4double be=0., my = 1115.683,
 av = 15.77,
 as = 18.34,
 ac = 0.71,
 asym = 23.21,
 k = 17.,
 c = 30.,
 D = 0.;
 if(mod(N,2) == 1 && mod(Z,2) == 1)D = -12./std::sqrt(A);
 if(mod(N,2) == 0 && mod(Z,2) == 0)D = 12./std::sqrt(A);
//
 G4double deltanew = (1.-std::exp(-1.*A/c))*D;
//
 be= av*A-as*std::pow(A,2./3.)-ac*Z*(Z-1.)/std::pow(A,1./3.)-asym*(N-Z)*(N-Z)/((1.+std::exp(-1.*A/k))*A)+deltanew + ny*(0.0335*my-26.7-48.7/std::pow(A,2.0/3.0));
 return be;
}

void G4Abla::unbound(G4double SN,G4double SP,G4double  SD,G4double ST,G4double SHE,G4double SA,G4double BP,G4double BD,G4double BT,G4double BHE,G4double BA,G4double *PROBF,G4double *PROBN,G4double *PROBP,G4double *PROBD,G4double *PROBT,G4double *PROBHE,G4double *PROBA,G4double *PROBIMF,G4double *PROBG,G4double *ECN,G4double *ECP,G4double *ECD,G4double *ECT,G4double *ECHE,G4double *ECA)
{
 G4double SBP = SP + BP;
 G4double SBD = SD + BD;
 G4double SBT = ST + BT;
 G4double SBHE = SHE + BHE;
 G4double SBA = SA + BA;

 G4double e = dmin1(SBP,SBD,SBT);
        e = dmin1(SBHE,SN,e);
        e = dmin1(SBHE,SBA,e);
//
 if(SN==e){
     *ECN = (-1.0)*SN;
     *ECP = 0.0;
     *ECD = 0.0;
     *ECT = 0.0;
     *ECHE = 0.0;
     *ECA = 0.0;
     *PROBN = 1.0;
     *PROBP = 0.0;
     *PROBD = 0.0;
     *PROBT = 0.0;
     *PROBHE = 0.0;
     *PROBA = 0.0;
     *PROBIMF = 0.0;
     *PROBF = 0.0;
     *PROBG = 0.0;
 }
 else if(SBP==e){
     *ECN = 0.0;
     *ECP = (-1.0)*SP + BP;
     *ECD = 0.0;
     *ECT = 0.0;
     *ECHE = 0.0;
     *ECA = 0.0;
     *PROBN = 0.0;
     *PROBP = 1.0;
     *PROBD = 0.0;
     *PROBT = 0.0;
     *PROBHE = 0.0;
     *PROBA = 0.0;
     *PROBIMF = 0.0;
     *PROBF = 0.0;
     *PROBG = 0.0;
 }
 else if(SBD==e){
     *ECN = 0.0;
     *ECD = (-1.0)*SD + BD;
     *ECP = 0.0;
     *ECT = 0.0;
     *ECHE = 0.0;
     *ECA = 0.0;
     *PROBN = 0.0;
     *PROBP = 0.0;
     *PROBD = 1.0;
     *PROBT = 0.0;
     *PROBHE = 0.0;
     *PROBA = 0.0;
     *PROBIMF = 0.0;
     *PROBF = 0.0;
     *PROBG = 0.0;
 }
 else if(SBT==e){
     *ECN = 0.0;
     *ECT = (-1.0)*ST + BT;
     *ECD = 0.0;
     *ECP = 0.0;
     *ECHE = 0.0;
     *ECA = 0.0;
     *PROBN = 0.0;
     *PROBP = 0.0;
     *PROBD = 0.0;
     *PROBT = 1.0;
     *PROBHE = 0.0;
     *PROBA = 0.0;
     *PROBIMF = 0.0;
     *PROBF = 0.0;
     *PROBG = 0.0;
 }
 else if(SBHE==e){
     *ECN = 0.0;
     *ECHE= (-1.0)*SHE + BHE;
     *ECD = 0.0;
     *ECT = 0.0;
     *ECP = 0.0;
     *ECA = 0.0;
     *PROBN = 0.0;
     *PROBP = 0.0;
     *PROBD = 0.0;
     *PROBT = 0.0;
     *PROBHE = 1.0;
     *PROBA = 0.0;
     *PROBIMF = 0.0;
     *PROBF = 0.0;
     *PROBG = 0.0;
 }
 else{ 
    if(SBA==e){
     *ECN = 0.0;
     *ECA = (-1.0)*SA + BA;
     *ECD = 0.0;
     *ECT = 0.0;
     *ECHE = 0.0;
     *ECP = 0.0;
     *PROBN = 0.0;
     *PROBP = 0.0;
     *PROBD = 0.0;
     *PROBT = 0.0;
     *PROBHE = 0.0;
     *PROBA = 1.0;
     *PROBIMF = 0.0;
     *PROBF = 0.0;
     *PROBG = 0.0;
    }
 }

 return;
}

void G4Abla::fissionDistri(G4double &A,G4double &Z,G4double &E,
		           G4double &a1,G4double &z1,G4double &e1,G4double &v1,
		           G4double &a2,G4double &z2,G4double &e2,G4double &v2,
                           G4double &vx_eva_sc,G4double &vy_eva_sc,G4double &vz_eva_sc, 
                           G4int *NbLam0_par)
{

/*
  Last update:

  21/01/17 - J.L.R.S. - Implementation of this fission model in C++


  Authors: K.-H. Schmidt, A. Kelic, M. V. Ricciardi,J. Benlliure, and 
           J.L.Rodriguez-Sanchez(1995 - 2017)

  On input: A, Z, E (mass, atomic number and exc. energy of compound nucleus
                     before fission)
  On output: Ai, Zi, Ei (mass, atomic number and (absolute) exc. energy of
                         fragment 1 and 2 after fission)

*/
  /* This program calculates isotopic distributions of fission fragments    */
  /* with a semiempirical model                                             */
  /* The width and eventually a shift in N/Z (polarization) follows the     */
  /* following rules:                                                       */
  /*                                                                        */
  /* The line N/Z following UCD has an angle of atan(Zcn/Ncn)               */
  /* to the horizontal axis on a chart of nuclides.                         */
/*   (For 238U the angle is 32.2 deg.)                                      */
/*                                                                        */
/*   The following relations hold: (from Armbruster)
c
c    sigma(N) (A=const) = sigma(Z) (A=const)
c    sigma(A) (N=const) = sigma(Z) (N=const)
c    sigma(A) (Z=const) = sigma(N) (Z=const)
c
c   From this we get:
c    sigma(Z) (N=const) * N = sigma(N) (Z=const) * Z
c    sigma(A) (Z=const) = sigma(Z) (A=const) * A/Z
c    sigma(N) (Z=const) = sigma(Z) (A=const) * A/Z
c    Z*sigma(N) (Z=const) = N*sigma(Z) (N=const) = A*sigma(Z) (A=const)     */
//

/*   Model parameters: 
C     These parameters have been adjusted to the compound nucleus 238U.
c     For the fission of another compound nucleus, it might be
c     necessary to slightly adjust some parameter values.
c     The most important ones are
C      Delta_U1_shell_max and
c      Delta_u2_shell.
*/
      G4double Nheavy1_in;   //  'position of shell for Standard 1'
      Nheavy1_in = 83.0;

      G4double Zheavy1_in;   //  'position of shell for Standard 1'
      Zheavy1_in = 50.0;

      G4double Nheavy2;   //  'position of heavy peak valley 2'
      Nheavy2 = 89.0;

      G4double Delta_U1_shell_max;  //  'Shell effect for valley 1'
      Delta_U1_shell_max = -2.45;

      G4double U1NZ_SLOPE;  // Reduction of shell effect with distance to 132Sn
      U1NZ_SLOPE = 0.2;

      G4double Delta_U2_shell;  //  'Shell effect for valley 2'
      Delta_U2_shell = -2.45;

      G4double X_s2s;   //  'Ratio (C_sad/C_scis) of curvature of potential'
      X_s2s = 0.8;

      G4double hbom1,hbom2,hbom3;   //  'Curvature of potential at saddle'
      hbom1 = 0.2;  // hbom1 is hbar * omega1 / (2 pi) !!!
      hbom2 = 0.2;  // hbom2 is hbar * omega2 / (2 pi) !!!
      hbom3 = 0.2;  // hbom3 is hbar * omega3 / (2 pi) !!!

      G4double Fwidth_asymm1,Fwidth_asymm2,Fwidth_symm;
//         'Factors for widths of distr. valley 1 and 2'
      Fwidth_asymm1 = 0.65;
      Fwidth_asymm2 = 0.65;
      Fwidth_symm   = 1.16;

      G4double xLevdens;   // 'Parameter x: a = A/x'
      xLevdens = 10.75;
//     The value of 1/0.093 = 10.75 is consistent with the
//     systematics of the mass widths of Ref. (RuI97).

      G4double FGAMMA; // 'Factor to gamma'
      FGAMMA = 1.;  // Theoretical expectation, not adjusted to data.
//     Additional factor to attenuation coefficient of shell effects
//     with increasing excitation energy

      G4double FGAMMA1; // 'Factor to gamma_heavy1'
      FGAMMA1 = 2.;
//     Adjusted to reduce the weight of Standard 1 with increasing
//     excitation energies, as required by experimental data.

      G4double FREDSHELL;
      FREDSHELL = 0.;
//     Adjusted to the reduced attenuation of shells in the superfluid region.
//     If FGAMMA is modified,
//     FGAMMA * FREADSHELL should remain constant (0.65) to keep
//     the attenuation of the shell effects below the critical
//     pairing energy ECRIT unchanged, which has been carefully
//     adjusted to the mass yields of Vives and Zoeller in this
//     energy range. A high value of FGAMMA leads ot a stronger
//     attenuation of shell effects above the superfluid region.

      G4double Ecrit;
      Ecrit = 5.;
//     The value of ECRIT determines the transition from a weak
//     decrease of the shell effect below ECRIT to a stronger
//     decrease above the superfluid range.
      const G4double d = 2.0;   // 'Surface distance of scission configuration'
     // d = 2.0;
//    Charge polarisation from Wagemanns p. 397: 
      G4double cpol1; // Charge polarisation standard I
      cpol1 = 0.35;  // calculated internally with shells
      G4double cpol2;  // Charge polarisation standard II
      cpol2 = 0.;  // calculated internally from LDM
      G4double Friction_factor;
      Friction_factor = 1.0;
      G4double Nheavy1;    // position of valley St 1 in Z and N
      G4double Delta_U1,Delta_U2; // used shell effects
      G4double cN_asymm1_shell, cN_asymm2_shell;
      G4double gamma,gamma_heavy1,gamma_heavy2; // fading of shells
      G4double E_saddle_scission;  // friction from saddle to scission
      G4double Ysymm=0.;    // Yield of symmetric mode
      G4double Yasymm1=0.;  // Yield of asymmetric mode 1
      G4double Yasymm2=0.;  // Yield of asymmetric mode 2
      G4double Nheavy1_eff; // Effective position of valley 1
      G4double Nheavy2_eff; // Effective position of valley 2
      G4double eexc1_saddle; // Excitation energy above saddle 1
      G4double eexc2_saddle; // Excitation energy above saddle 2
      G4double EEXC_MAX;  // Excitation energy above lowest saddle
      G4double r_e_o;  // Even-odd effect in Z
      G4double cN_symm;  // Curvature of symmetric valley
      G4double CZ;  // Curvature of Z distribution for fixed A
      G4double Nheavy2_NZ;  // Position of Shell 2, combined N and Z
      G4double N;
      G4double Aheavy1,Aheavy2;
      G4double Sasymm1=0.,Sasymm2=0.,Ssymm=0.,Ysum=0.,Yasymm=0.;
      G4double Ssymm_mode1,Ssymm_mode2;
      G4double wNasymm1_saddle, wNasymm2_saddle, wNsymm_saddle;
      G4double wNasymm2_scission, wNsymm_scission;
      G4double wNasymm1, wNasymm2, wNsymm;
      G4int imode;
      G4double  rmode;
      G4double ZA1width;
      G4double N1r,N2r,A1r,N1,N2;
      G4double Zsymm,Nsymm;
      G4double N1mean, N1width;
      G4double dUeff;
      /* effective shell effect at lowest barrier */
      G4double Eld;
      /* Excitation energy with respect to ld barrier */
      G4double re1,re2,re3;
      G4double eps1,eps2;
      G4double Z1UCD,Z2UCD;
      G4double beta=0.,beta1=0.,beta2=0.;
     // double betacomplement;
      G4double DN1_POL;
      /* shift of most probable neutron number for given Z,
            according to polarization */
      G4int i_help;
      G4double A_levdens;
            /* level-density parameter */
     // double A_levdens_light1,A_levdens_light2;
      G4double A_levdens_heavy1,A_levdens_heavy2;

      G4double R0=1.16;

      G4double epsilon_1_saddle,epsilon0_1_saddle;
      G4double epsilon_2_saddle,epsilon0_2_saddle,epsilon_symm_saddle;
      G4double epsilon_1_scission;//,epsilon0_1_scission;
      G4double epsilon_2_scission;//,epsilon0_2_scission;
      G4double epsilon_symm_scission;
                                    /* modified energy */
      G4double E_eff1_saddle,E_eff2_saddle;
      G4double Epot0_mode1_saddle,Epot0_mode2_saddle,Epot0_symm_saddle;
      G4double Epot_mode1_saddle,Epot_mode2_saddle,Epot_symm_saddle;
      G4double E_defo,E_defo1,E_defo2,E_scission_pre=0.,E_scission_post;
      G4double E_asym;
      G4double E1exc=0.,E2exc=0.;
      G4double E1exc_sigma,E2exc_sigma;
      G4double TKER;
      G4double EkinR1,EkinR2;
      G4double MassCurv_scis, MassCurv_sadd;
      G4double cN_symm_sadd;
      G4double Nheavy1_shell,Nheavy2_shell;
      G4double wNasymm1_scission;
      G4double Aheavy1_eff,Aheavy2_eff;
      G4double Z1rr,Z1r;
      G4double E_HELP;
      G4double Z_scission,N_scission,A_scission;
      G4double Z2_over_A_eff;
      G4double beta1gs=0.,beta2gs=0.,betags=0.;
      G4double sigZmin;   // 'Minimum neutron width for constant Z'
      G4double DSN132,Delta_U1_shell,E_eff0_saddle;//,e_scission;
      G4int NbLam0= (*NbLam0_par);
      //
      sigZmin = 0.5;
      N = A - Z;  /*  neutron number of the fissioning nucleus  */
//
      cN_asymm1_shell = 0.700 * N/Z;
      cN_asymm2_shell = 0.040 * N/Z;

//*********************************************************************

      DSN132 = Nheavy1_in - N/Z * Zheavy1_in;
      Aheavy1 = Nheavy1_in + Zheavy1_in + 0.340 * DSN132;
     /* Neutron number of valley Standard 1 */
     /* It is assumed that the 82-neutron shell effect is stronger than
c         the 50-proton shell effect. Therefore, the deviation in N/Z of
c         the fissioning nucleus from the N/Z of 132Sn will
c         change the position of the combined shell in mass. For neutron-
c         deficient fissioning nuclei, the mass will increase and vice
c         versa.  */

      Delta_U1_shell = Delta_U1_shell_max + U1NZ_SLOPE * std::abs(DSN132);
      Delta_U1_shell = min(0.,Delta_U1_shell);
      /* Empirical reduction of shell effect with distance in N/Z of CN to 132Sn */
      /* Fits (239U,n)f and 226Th e.-m.-induced fission */

      Nheavy1 = N/A * Aheavy1;   /* UCD */
      Aheavy2 = Nheavy2 * A/N;

      Zsymm  = Z / 2.0; /* proton number in symmetric fission (centre) */
      Nsymm  = N / 2.0;
      A_levdens = A / xLevdens;
      gamma = A_levdens / (0.40 * std::pow(A,1.3333)) * FGAMMA;
      A_levdens_heavy1 = Aheavy1 / xLevdens;
      gamma_heavy1 = A_levdens_heavy1 / (0.40 * std::pow(Aheavy1,1.3333)) * FGAMMA * FGAMMA1;
      A_levdens_heavy2 = Aheavy2 / xLevdens;
      gamma_heavy2 = A_levdens_heavy2 / (0.40 * std::pow(Aheavy2,1.3333)) * FGAMMA;

//     Energy dissipated from saddle to scission
//     F. Rejmund et al., Nucl. Phys. A 678 (2000) 215, fig. 4 b    */
      E_saddle_scission = (-24. + 0.02227 * Z*Z/std::pow(A,0.33333))*Friction_factor;
      E_saddle_scission = max( 0.0, E_saddle_scission );

//     Fit to experimental result on curvature of potential at saddle
//     Parametrization of T. Enqvist according to Mulgin et al. 1998
//     MassCurv taken at scission.    */

       Z2_over_A_eff = Z*Z/A;

      if( Z2_over_A_eff< 34.0 )
        MassCurv_scis = std::pow(10., -1.093364 + 0.082933 * Z2_over_A_eff - 0.0002602 * Z2_over_A_eff*Z2_over_A_eff);
      else
        MassCurv_scis = std::pow(10., 3.053536 - 0.056477 * Z2_over_A_eff+ 0.0002454 * Z2_over_A_eff*Z2_over_A_eff );

//     to do:
//     fix the X with the channel intensities of 226Th (KHS at SEYSSINS,1998)
//     replace then (all) cN_symm by cN_symm_saddle (at least for Yields)
      MassCurv_sadd = X_s2s * MassCurv_scis;

      cN_symm      = 8.0 / std::pow(N,2.) * MassCurv_scis;
      cN_symm_sadd = 8.0 / std::pow(N,2.) * MassCurv_sadd;

      Nheavy1_shell = Nheavy1;

      if(E < 100.0)
      Nheavy1_eff = (cN_symm_sadd*Nsymm + cN_asymm1_shell * 
                     Uwash(E/A*Aheavy1,Ecrit,FREDSHELL,gamma_heavy1) *
                     Nheavy1_shell)
                     / (cN_symm_sadd +
                       cN_asymm1_shell *
                     Uwash(E/A*Aheavy1,Ecrit,FREDSHELL,gamma_heavy1));
      else
      Nheavy1_eff = (cN_symm_sadd*Nsymm +
                    cN_asymm1_shell*Nheavy1_shell)
                  / (cN_symm_sadd +
                    cN_asymm1_shell);
      
      /* Position of Standard II defined by neutron shell */
      Nheavy2_NZ = Nheavy2;
      Nheavy2_shell = Nheavy2_NZ;
      if (E < 100.)
      Nheavy2_eff = (cN_symm_sadd*Nsymm +
                    cN_asymm2_shell*
                    Uwash(E/A*Aheavy2,Ecrit,FREDSHELL,gamma_heavy2) *
                    Nheavy2_shell)
                  / (cN_symm_sadd +
                    cN_asymm2_shell*
                    Uwash(E/A*Aheavy2,Ecrit,FREDSHELL,gamma_heavy2));
      else
      Nheavy2_eff = (cN_symm_sadd*Nsymm +
                    cN_asymm2_shell*Nheavy2_shell)
                  / (cN_symm_sadd +
                    cN_asymm2_shell);
     
      Delta_U1 = Delta_U1_shell + (Nheavy1_shell - Nheavy1_eff)*(Nheavy1_shell - Nheavy1_eff) * cN_asymm1_shell; /* shell effect in valley of mode 1 */
      Delta_U1 = min(Delta_U1,0.0);
      Delta_U2 = Delta_U2_shell + (Nheavy2_shell - Nheavy2_eff)*(Nheavy2_shell - Nheavy2_eff) * cN_asymm2_shell; /* shell effect in valley of mode 2 */
      Delta_U2 = min(Delta_U2,0.0);

//    liquid drop energies at the centres of the different shell effects
//    with respect to liquid drop at symmetry
      Epot0_mode1_saddle = (Nheavy1_eff-Nsymm)*(Nheavy1_eff-Nsymm) * cN_symm_sadd;
      Epot0_mode2_saddle = (Nheavy2_eff-Nsymm)*(Nheavy2_eff-Nsymm) * cN_symm_sadd;
      Epot0_symm_saddle = 0.0;

//    energies including shell effects at the centres of the different
//    shell effects with respect to liquid drop at symmetry  */
      Epot_mode1_saddle = Epot0_mode1_saddle + Delta_U1;
      Epot_mode2_saddle = Epot0_mode2_saddle + Delta_U2;
      Epot_symm_saddle = Epot0_symm_saddle;

//    minimum of potential with respect to ld potential at symmetry
      dUeff = min( Epot_mode1_saddle, Epot_mode2_saddle);
      dUeff = min( dUeff, Epot_symm_saddle);
      dUeff = dUeff - Epot_symm_saddle;

      Eld = E + dUeff;
//     E   = energy above lowest effective barrier
//     Eld = energy above liquid-drop barrier
//     Due to this treatment the energy E on input means the excitation
//     energy above the lowest saddle.                                  */

//    excitation energies at saddle modes 1 and 2 without shell effect  */
      epsilon0_1_saddle = Eld - Epot0_mode1_saddle;
      epsilon0_2_saddle = Eld - Epot0_mode2_saddle;

//    excitation energies at saddle modes 1 and 2 with shell effect */
      epsilon_1_saddle = Eld - Epot_mode1_saddle;
      epsilon_2_saddle = Eld - Epot_mode2_saddle;

      epsilon_symm_saddle = Eld - Epot_symm_saddle;
//    epsilon_symm_saddle = Eld - dUeff;

      eexc1_saddle = epsilon_1_saddle;
      eexc2_saddle = epsilon_2_saddle;

//    EEXC_MAX is energy above the lowest saddle */
      EEXC_MAX = max( eexc1_saddle, eexc2_saddle);
      EEXC_MAX = max( EEXC_MAX, Eld);

//    excitation energy at scission */
      epsilon_1_scission = Eld + E_saddle_scission - Epot_mode1_saddle;
      epsilon_2_scission = Eld + E_saddle_scission - Epot_mode2_saddle;

//    excitation energy of symmetric fragment at scission  */
      epsilon_symm_scission = Eld + E_saddle_scission - Epot_symm_saddle;

//    calculate widhts at the saddle
      E_eff1_saddle = epsilon0_1_saddle - Delta_U1 *
      Uwash(epsilon_1_saddle/A*Aheavy1,Ecrit,FREDSHELL,gamma_heavy1);

      if( E_eff1_saddle < A_levdens * hbom1*hbom1)
          E_eff1_saddle = A_levdens * hbom1*hbom1;
      
      wNasymm1_saddle =
        std::sqrt(0.50 * std::sqrt(1.0/A_levdens*E_eff1_saddle) /
        (cN_asymm1_shell *
      Uwash(epsilon_1_saddle/A*Aheavy1,Ecrit,FREDSHELL,gamma_heavy1)+
        cN_symm_sadd));

      E_eff2_saddle = epsilon0_2_saddle -
      Delta_U2 *
      Uwash(epsilon_2_saddle/A*Aheavy2,Ecrit,FREDSHELL,gamma_heavy2);

      if(E_eff2_saddle < A_levdens * hbom2*hbom2)
          E_eff2_saddle = A_levdens * hbom2*hbom2;
      
      wNasymm2_saddle =
        std::sqrt(0.50 * std::sqrt(1.0/A_levdens*E_eff2_saddle) /
        (cN_asymm2_shell *
       Uwash(epsilon_2_saddle/A*Aheavy2,Ecrit,FREDSHELL,gamma_heavy2)+
        cN_symm_sadd));

      E_eff0_saddle = epsilon_symm_saddle;
      if(E_eff0_saddle < A_levdens * hbom3*hbom3)
          E_eff0_saddle = A_levdens * hbom3*hbom3;
      
      wNsymm_saddle =
          std::sqrt(0.50 * std::sqrt(1.0/A_levdens*E_eff0_saddle) /
                    cN_symm_sadd);

      if(epsilon_symm_scission > 0.0 ){
        E_HELP = max(E_saddle_scission,epsilon_symm_scission);
        wNsymm_scission =
        std::sqrt(0.50 * std::sqrt(1.0/A_levdens*(E_HELP)) /
                    cN_symm);
      }else{
        wNsymm_scission =
        std::sqrt(0.50 * std::sqrt(1.0/A_levdens*E_saddle_scission) /
                    cN_symm);
      }

//    Calculate widhts at the scission point: 
//    fits of ref. Beizin 1991 (Plots by Sergei Zhdanov) 

      if( E_saddle_scission == 0.0 ){
        wNasymm1_scission = wNasymm1_saddle;
        wNasymm2_scission = wNasymm2_saddle;
      }else{
        if( Nheavy1_eff > 75.0 ){
          wNasymm1_scission = std::sqrt(21.0)*N/A;
          wNasymm2_scission = max( 12.8 - 1.0 *(92.0 - Nheavy2_eff),1.0)*N/A;

        }else{
           wNasymm1_scission = wNasymm1_saddle;
           wNasymm2_scission = wNasymm2_saddle;
        }
      }

      wNasymm1_scission = max( wNasymm1_scission, wNasymm1_saddle );
      wNasymm2_scission = max( wNasymm2_scission, wNasymm2_saddle );

      wNasymm1 = wNasymm1_scission * Fwidth_asymm1;
      wNasymm2 = wNasymm2_scission * Fwidth_asymm2;
      wNsymm   = wNsymm_scission * Fwidth_symm;

//     mass and charge of fragments using UCD, needed for level densities
      Aheavy1_eff = Nheavy1_eff * A/N;
      Aheavy2_eff = Nheavy2_eff * A/N;

      A_levdens_heavy1 = Aheavy1_eff / xLevdens;
      A_levdens_heavy2 = Aheavy2_eff / xLevdens;
      gamma_heavy1 = A_levdens_heavy1 / (0.40 * std::pow(Aheavy1_eff,1.3333)) * FGAMMA * FGAMMA1;
      gamma_heavy2 = A_levdens_heavy2 / (0.40 * std::pow(Aheavy2_eff,1.3333)) * FGAMMA;

      if( epsilon_symm_saddle < A_levdens * hbom3*hbom3)
        Ssymm = 2.0 * std::sqrt(A_levdens*A_levdens * hbom3*hbom3) +
        (epsilon_symm_saddle - A_levdens * hbom3*hbom3)/hbom3;
      else
        Ssymm = 2.0 * std::sqrt(A_levdens*epsilon_symm_saddle);
      
      Ysymm = 1.0;

      if( epsilon0_1_saddle < A_levdens * hbom1*hbom1 )
        Ssymm_mode1 = 2.0 * std::sqrt(A_levdens*A_levdens * hbom1*hbom1) +
             (epsilon0_1_saddle - A_levdens * hbom1*hbom1)/hbom1;
      else
        Ssymm_mode1 = 2.0 * std::sqrt( A_levdens*epsilon0_1_saddle );

      if( epsilon0_2_saddle < A_levdens * hbom2*hbom2 )
        Ssymm_mode2 = 2.0 * std::sqrt(A_levdens*A_levdens * hbom2*hbom2) +
             (epsilon0_2_saddle - A_levdens * hbom2*hbom2)/hbom2;
      else
        Ssymm_mode2 = 2.0 * std::sqrt(A_levdens*epsilon0_2_saddle);


      if( epsilon0_1_saddle -
         Delta_U1*
       Uwash(epsilon_1_saddle/A*Aheavy1,Ecrit,FREDSHELL,gamma_heavy1)
           < A_levdens * hbom1*hbom1 )
        Sasymm1 = 2.0 * std::sqrt( A_levdens*A_levdens * hbom1*hbom1 ) +
         (epsilon0_1_saddle - Delta_U1 *
       Uwash(epsilon_1_saddle/A*Aheavy1,Ecrit,FREDSHELL,gamma_heavy1)
         - A_levdens * hbom1*hbom1)/hbom1;
      else
        Sasymm1 = 2.0 *std::sqrt( A_levdens*(epsilon0_1_saddle - Delta_U1 *
        Uwash(epsilon_1_saddle/A*Aheavy1,Ecrit,FREDSHELL,gamma_heavy1)));
      
      if( epsilon0_2_saddle -
        Delta_U2*
        Uwash(epsilon_2_saddle/A*Aheavy2,Ecrit,FREDSHELL,gamma_heavy2)
        < A_levdens * hbom2*hbom2 )
        Sasymm2 = 2.0 * std::sqrt( A_levdens*A_levdens * hbom2*hbom2 ) +
         (epsilon0_1_saddle-Delta_U1 *
         Uwash(epsilon_2_saddle/A*Aheavy2,Ecrit,FREDSHELL,gamma_heavy2)
         - A_levdens * hbom2*hbom2)/hbom2;
      else
        Sasymm2 = 2.0 *
             std::sqrt( A_levdens*(epsilon0_2_saddle - Delta_U2 *
       Uwash(epsilon_2_saddle/A*Aheavy2,Ecrit,FREDSHELL,gamma_heavy2)));

      Yasymm1 = ( std::exp(Sasymm1 - Ssymm) - std::exp(Ssymm_mode1 - Ssymm) ) *
                 wNasymm1_saddle / wNsymm_saddle * 2.0;

      Yasymm2 = ( std::exp(Sasymm2 - Ssymm) - std::exp(Ssymm_mode2 - Ssymm) ) *
                 wNasymm2_saddle / wNsymm_saddle * 2.0;

      Ysum = Ysymm + Yasymm1 + Yasymm2;  /* normalize */

      if( Ysum > 0.00 ){
        Ysymm = Ysymm / Ysum;
        Yasymm1 = Yasymm1 / Ysum;
        Yasymm2 = Yasymm2 / Ysum;
        Yasymm = Yasymm1 + Yasymm2;
      }else{
        Ysymm = 0.0;
        Yasymm1 = 0.0;
        Yasymm2 = 0.0;
//       search minimum threshold and attribute all events to this mode */
        if( (epsilon_symm_saddle < epsilon_1_saddle) &&
            (epsilon_symm_saddle < epsilon_2_saddle) )
          Ysymm = 1.0;
        else
          if( epsilon_1_saddle < epsilon_2_saddle )
            Yasymm1 = 1.0;
          else
            Yasymm2 = 1.0;
      }
  // even-odd effect
  // Parametrization from Rejmund et al. 
     if (mod(Z,2.0)== 0)
      r_e_o = std::pow(10.0,-0.0170 * (E_saddle_scission + Eld)*(E_saddle_scission + Eld));
     else
      r_e_o = 0.0;
     
/*     -------------------------------------------------------
c     selecting the fission mode using the yields at scission
c     -------------------------------------------------------
c     random decision: symmetric or asymmetric
c     IMODE = 1 means asymmetric fission, mode 1
c     IMODE = 2 means asymmetric fission, mode 2
c     IMODE = 3 means symmetric fission
c     testcase: 238U, E*= 6 MeV :    6467   8781   4752   (20000)
c                                  127798 176480  95722  (400000)
c                                  319919 440322 239759 (1000000)
c                     E*=12 MeV :  153407 293063 553530 (1000000) */

 fiss321:  // rmode = DBLE(HAZ(k))
        rmode = G4AblaRandom::flat();
      if( rmode < Yasymm1 )
         imode = 1;
      else
         if( (rmode > Yasymm1) && (rmode < Yasymm) )
            imode = 2;
         else
            imode = 3;
         
//    determine parameters of the neutron distribution of each mode
//    at scission

      if( imode == 1){
         N1mean = Nheavy1_eff;
         N1width = wNasymm1;
      }else{
       if( imode == 2 ){
         N1mean = Nheavy2_eff;
         N1width = wNasymm2;
       }else{
        //if( imode == 3 ) then
         N1mean = Nsymm;
         N1width = wNsymm;
       }
      }

//     N2mean needed by CZ below
    //  N2mean = N - N1mean;
      
//     fission mode found, then the determination of the
//     neutron numbers N1 and N2 at scission by randon decision
      N1r = 1.0;
      N2r = 1.0;
      while( N1r < 5.0 || N2r < 5.0 ){
       //  N1r = DBLE(GaussHaz(k,sngl(N1mean), sngl(N1width) ))
        // N1r = N1mean+G4AblaRandom::gaus(N1width);//
         N1r = gausshaz(0,N1mean,N1width);
         N2r = N - N1r;
      }

//     --------------------------------------------------
//     first approximation of fission fragments using UCD at saddle
//     --------------------------------------------------
      Z1UCD = Z/N * N1r;
      Z2UCD = Z/N * N2r;
      A1r = A/N * N1r;
//
//     --------------------------
//     deformations: starting ...
//     --------------------------  */
      if( imode == 1 ){
// ---   N = 82  */
        E_scission_pre = max( epsilon_1_scission, 1.0 );
//   ! Eexc at scission, neutron evaporation from saddle to scission not considered */
        if( N1mean > N*0.50 ){
          beta1  = 0.0;                   /*   1. fragment is spherical */
          beta2  = 0.55;                 /*   2. fragment is deformed  0.5*/
        }else{
          beta1  = 0.55;                 /*  1. fragment is deformed 0.5*/
          beta2  = 0.00;                  /*  2. fragment is spherical */
        }
       }
      if( imode == 2 ){
// ---   N appr. 86  */
        E_scission_pre = max( epsilon_2_scission, 1.0 );
        if( N1mean > N*0.50 ){
          beta1  = (N1r - 92.0) * 0.030 + 0.60;

          beta1gs = ecld->beta2[idint(N1r)][idint(Z1UCD)];
          beta2gs = ecld->beta2[idint(N2r)][idint(Z2UCD)];

          beta1 = max(beta1,beta1gs);
          beta2  = 1.0 - beta1;
          beta2 = max(beta2,beta2gs);
        }else{

          beta1gs = ecld->beta2[idint(N1r)][idint(Z1UCD)];
          beta2gs = ecld->beta2[idint(N2r)][idint(Z2UCD)];

          beta2  = (N2r -92.0) * 0.030 + 0.60;
          beta2 = max(beta2,beta2gs);
          beta1 = 1.0 - beta2;
          beta1 = max(beta1,beta1gs);
        }
      }
      beta = 0.0;
      if( imode == 3 ){
//      if( imode >0 ){
// ---   Symmetric fission channel
//       the fit function for beta is the deformation for optimum energy
//       at the scission point, d = 2
//       beta  : deformation of symmetric fragments
//       beta1 : deformation of first fragment
//       beta2 : deformation of second fragment
        betags = ecld->beta2[idint(Nsymm)][idint(Zsymm)];
        beta1gs = ecld->beta2[idint(N1r)][idint(Z1UCD)];
        beta2gs = ecld->beta2[idint(N2r)][idint(Z2UCD)];
        beta  = max(0.177963+0.0153241*Zsymm-1.62037e-4*Zsymm*Zsymm,betags);
        beta1 = max(0.177963+0.0153241*Z1UCD-1.62037e-4*Z1UCD*Z1UCD,beta1gs);
        beta2 = max(0.177963+0.0153241*Z2UCD-1.62037e-4*Z2UCD*Z2UCD,beta2gs);

        E_asym = frldm( Z1UCD, N1r, beta1 ) +
              frldm( Z2UCD, N2r, beta2 ) +
              ecoul( Z1UCD, N1r, beta1, Z2UCD, N2r, beta2, 2.0 ) -
              2.0 * frldm( Zsymm, Nsymm, beta ) -
              ecoul( Zsymm, Nsymm, beta, Zsymm, Nsymm, beta, 2.0 );
        E_scission_pre = max( epsilon_symm_scission - E_asym, 1. );
      }
//     -----------------------
//     ... end of deformations
//     -----------------------

//     ------------------------------------------
//     evaporation from saddle to scission ...
//     ------------------------------------------     
      if(E_scission_pre>5. && NbLam0<1){
       evap_postsaddle(A,Z,E_scission_pre,&E_scission_post,
        &A_scission,&Z_scission,vx_eva_sc,vy_eva_sc,vz_eva_sc,&NbLam0);
       N_scission = A_scission - Z_scission;
      }else{
       A_scission = A;
       Z_scission = Z;
       E_scission_post = E_scission_pre;
       N_scission = A_scission - Z_scission;
      }
//     ---------------------------------------------------
//     second approximation of fission fragments using UCD
//     --------------------------------------------------- */
//
      N1r = N1r * N_scission / N;
      N2r = N2r * N_scission / N;
      Z1UCD = Z1UCD * Z_scission / Z;
      Z2UCD = Z2UCD * Z_scission / Z;
      A1r = Z1UCD + N1r;

//     ---------------------------------------------------------
//     determination of the charge and mass of the fragments ...
//     ---------------------------------------------------------

//     - CZ is the curvature of charge distribution for fixed mass,
//       common to all modes, gives the width of the charge distribution.
//       The physics picture behind is that the division of the
//       fissioning nucleus in N and Z is slow when mass transport from
//       one nascent fragment to the other is concerned but fast when the
//       N/Z degree of freedom is concernded. In addition, the potential
//       minima in direction of mass transport are broad compared to the
//       potential minimum in N/Z direction.
//          The minima in direction of mass transport are calculated
//          by the liquid-drop (LD) potential (for superlong mode),
//          by LD + N=82 shell (for standard 1 mode) and
//          by LD + N=86 shell (for standard 2 mode).
//          Since the variation of N/Z is fast, it can quickly adjust to
//          the potential and is thus determined close to scission.
//          Thus, we calculate the mean N/Z and its width for fixed mass
//          at scission.
//          For the SL mode, the mean N/Z is calculated by the
//          minimum of the potential at scission as a function of N/Z for
//          fixed mass.
//          For the S1 and S2 modes, this correlation is imposed by the
//          empirical charge polarisation.
//          For the SL mode, the fluctuation in this width is calculated
//          from the curvature of the potential at scission as a function
//          of N/Z. This value is also used for the widths of S1 and S2.


//     Polarisation assumed for standard I and standard II:
//      Z - Zucd = cpol (for A = const);
//      from this we get (see remarks above)
//      Z - Zucd =  Acn/Ncn * cpol (for N = const)   */
//
      CZ = ( frldm( Z1UCD-1.0, N1r+1.0, beta1 ) +
             frldm( Z2UCD+1.0, N2r-1.0, beta2 ) +
             frldm( Z1UCD+1.0, N1r-1.0, beta1 ) +
             frldm( Z2UCD-1.0, N2r+1.0, beta2 ) +
             ecoul( Z1UCD-1.0, N1r+1.0, beta1,
                    Z2UCD+1.0, N2r-1.0, beta2, 2.0) +
             ecoul( Z1UCD+1.0, N1r-1.0, beta1,
                    Z2UCD-1.0, N2r+1.0, beta2, 2.0) -
         2.0*ecoul( Z1UCD, N1r, beta1, Z2UCD, N2r, beta2, 2.0) -
         2.0*frldm( Z1UCD, N1r, beta1 ) -
         2.0*frldm( Z2UCD, N2r, beta2) ) * 0.50;
//
      if(1.0/A_levdens*E_scission_post < 0.0)
        std::cout << "DSQRT 1 < 0" << A_levdens << " " << E_scission_post << std::endl;
      
      if(0.50 * std::sqrt(1.0/A_levdens*E_scission_post) / CZ < 0.0){
        std::cout << "DSQRT 2 < 0 " << CZ << std::endl;
        std::cout << "This event was not considered" << std::endl;
        goto fiss321;
      }

      ZA1width = std::sqrt(0.5*std::sqrt(1.0/A_levdens*E_scission_post)/CZ);

//     Minimum width in N/Z imposed.
//     Value of minimum width taken from 235U(nth,f) data
//     sigma_Z(A=const) = 0.4 to 0.5  (from Lang paper Nucl Phys. A345 (1980) 34)
//     sigma_N(Z=const) = 0.45 * A/Z  (= 1.16 for 238U)
//      therefore: SIGZMIN = 1.16                                              
//     Physics; variation in N/Z for fixed A assumed.
//      Thermal energy at scission is reduced by
//      pre-scission neutron evaporation"

       ZA1width = max(ZA1width,sigZmin);

      if(imode == 1 && cpol1 != 0.0){
//       --- asymmetric fission, mode 1 */
       G4int IS = 0;
       fiss2801:
       Z1rr = Z1UCD - cpol1 * A_scission/N_scission;
     // Z1r = DBLE(GaussHaz(k,sngl(Z1rr), sngl(ZA1width) ));
      // Z1r = Z1rr+G4AblaRandom::gaus(ZA1width);//
       Z1r =gausshaz(0,Z1rr,ZA1width);
       IS = IS +1;
       if(IS>100){
       std::cout << "WARNING: GAUSSHAZ CALLED MORE THAN 100 TIMES WHEN CALCULATING Z1R IN PROFI.FOR. A VALUE WILL BE FORCED" << std::endl;
         Z1r = Z1rr;
       }
       if ((utilabs(Z1rr - Z1r) > 3.0*ZA1width) || Z1r<1.0)goto fiss2801;
       N1r = A1r - Z1r;
      }else{
        if( imode == 2 && cpol2 != 0.0 ){
//       --- asymmetric fission, mode 2 */
        G4int IS = 0;
        fiss2802:
        Z1rr = Z1UCD - cpol2 * A_scission/N_scission;
        //Z1r = Z1rr+G4AblaRandom::gaus(ZA1width);//
        Z1r = gausshaz(0,Z1rr,ZA1width);
        IS = IS +1;
        if(IS>100){
        std::cout << "WARNING: GAUSSHAZ CALLED MORE THAN 100 TIMES WHEN CALCULATING Z1R IN PROFI.FOR. A VALUE WILL BE FORCED" << std::endl;
         Z1r = Z1rr;
        }
        if( (utilabs(Z1rr - Z1r) > 3.0*ZA1width) || Z1r < 1.0 ) goto fiss2802;
        N1r = A1r - Z1r;
        }else{
//      Otherwise do; /* Imode = 3 in any case; imode = 1 and 2 for CPOL = 0 */
//       and symmetric case     */
//         We treat a simultaneous split in Z and N to determine
//         polarisation  */

          re1 = frldm( Z1UCD-1.0, N1r+1.0, beta1 ) +
                frldm( Z2UCD+1.0, N2r-1.0, beta2 ) +
                ecoul( Z1UCD-1.0, N1r+1.0, beta1,
                      Z2UCD+1.0, N2r-1.0, beta2, d ); /* d = 2 fm */
          re2 = frldm( Z1UCD, N1r, beta1) +
                frldm( Z2UCD, N2r, beta2 ) +
                ecoul( Z1UCD, N1r, beta1,
                      Z2UCD, N2r, beta2, d );  /*  d = 2 fm */
          re3 = frldm( Z1UCD+1.0, N1r-1.0, beta1 ) +
                frldm( Z2UCD-1.0, N2r+1.0, beta2 ) +
                ecoul( Z1UCD+1.0, N1r-1.0, beta1,
                      Z2UCD-1.0, N2r+1.0, beta2, d ); /* d = 2 fm */
          eps2 = ( re1 - 2.0*re2 + re3 ) / 2.0;
          eps1 = ( re3 - re1 ) / 2.0;
          DN1_POL = -eps1 / ( 2.0 * eps2 );
//
          Z1rr = Z1UCD + DN1_POL;

//       Polarization of Standard 1 from shell effects around 132Sn
          if ( imode == 1 ){
            if ( Z1rr > 50.0 ){
              DN1_POL = DN1_POL - 0.6 * Uwash(E_scission_post,Ecrit,FREDSHELL,gamma);
              Z1rr = Z1UCD + DN1_POL;
              if ( Z1rr < 50. ) Z1rr = 50.0;
            }else{
              DN1_POL = DN1_POL + 0.60 * Uwash(E_scission_post,Ecrit,FREDSHELL,gamma);
              Z1rr = Z1UCD + DN1_POL;
              if ( Z1rr > 50.0 ) Z1rr = 50.0;
            }
          }

        G4int IS = 0;
        fiss2803:      
        //Z1r = Z1rr+G4AblaRandom::gaus(ZA1width);
        Z1r = gausshaz(0,Z1rr,ZA1width);
        IS = IS +1;
        if(IS>100){
        std::cout << "WARNING: GAUSSHAZ CALLED MORE THAN 100 TIMES WHEN CALCULATING Z1R IN PROFI.FOR. A VALUE WILL BE FORCED" << std::endl;
         Z1r = Z1rr;
        }

        if( (utilabs(Z1rr - Z1r) > 3.0*ZA1width) || (Z1r < 1.0) )goto fiss2803;
        N1r = A1r - Z1r;

        }
      }

//     ------------------------------------------
//     Integer proton number with even-odd effect
//     ------------------------------------------ 
      even_odd(Z1r, r_e_o, i_help);

      z1 = G4double(i_help);
      z2 = dint( Z_scission ) - z1;
      N1 = dint( N1r );
      N2 = dint( N_scission ) - N1;
      a1 = z1 + N1;
      a2 = z2 + N2;

      if( (z1 < 0) || (z2 < 0) || (a1 < 0) || (a2 < 0) ){
         std::cout << " -------------------------------" << std::endl;
         std::cout << " Z, A, N : " << Z  << " " << A << " " << N << std::endl;
         std::cout << z1 << " " << z2 << " " << a1 << " " << a2 << std::endl;
         std::cout << E_scission_post << " " << A_levdens << " " << CZ <<  std::endl;

         std::cout << " -------------------------------" << std::endl;
      }

//     -----------------------
//     excitation energies ...
//     -----------------------
//
      if( imode == 1 ){
// ----  N = 82
        if( N1mean > N*0.50 ){
//         (a) 1. fragment is spherical and  2. fragment is deformed */
          E_defo = 0.0;
          beta2gs = ecld->beta2[idint(N2)][idint(z2)];
          if(beta2< beta2gs) beta2 = beta2gs;
          E1exc = E_scission_pre * a1 / A + E_defo;
          E_defo = frldm( z2, N2, beta2 ) - frldm( z2, N2, beta2gs );
          E2exc = E_scission_pre * a2 / A + E_defo;
        }else{
//         (b) 1. fragment is deformed and  2. fragment is spherical */
          beta1gs = ecld->beta2[idint(N1)][idint(z1)];
          if(beta1< beta1gs) beta1 = beta1gs;
          E_defo = frldm( z1, N1, beta1 ) - frldm( z1, N1, beta1gs );
          E1exc = E_scission_pre * a1 / A + E_defo;
          E_defo = 0.0;
          E2exc = E_scission_pre * a2 / A + E_defo;
        }
      }


      if( imode == 2 ){
// ---   N appr. 86 */
        if( N1mean > N*0.5 ){       
          /*  2. fragment is spherical */
          beta1gs = ecld->beta2[idint(N1)][idint(z1)];
          if(beta1< beta1gs) beta1 = beta1gs;
          E_defo = frldm( z1, N1, beta1 ) - frldm( z1, N1, beta1gs );
          E1exc = E_scission_pre * a1 / A + E_defo;
          beta2gs = ecld->beta2[idint(N2)][idint(z2)];
          if(beta2< beta2gs) beta2 = beta2gs;
          E_defo = frldm( z2, N2, beta2 ) - frldm( z2, N2, beta2gs );
          E2exc = E_scission_pre * a2 / A + E_defo;
        }else{                           
          /*  1. fragment is spherical */
          beta2gs = ecld->beta2[idint(N2)][idint(z2)];
          if(beta2< beta2gs) beta2 = beta2gs;
          E_defo = frldm( z2, N2, beta2 ) - frldm( z2, N2, beta2gs );
          E2exc = E_scission_pre * a2 / A + E_defo;
          beta1gs = ecld->beta2[idint(N1)][idint(z1)];
          if(beta1< beta1gs) beta1 = beta1gs;
          E_defo = frldm( z1, N1, beta1 ) - frldm( z1, N1, beta1gs );
          E1exc = E_scission_pre * a1 / A + E_defo;
        }
      }

      if( imode == 3 ){
// ---   Symmetric fission channel
          beta1gs = ecld->beta2[idint(N1)][idint(z1)];
          if(beta1< beta1gs) beta1 = beta1gs;
          beta2gs = ecld->beta2[idint(N2)][idint(z2)];
          if(beta2< beta2gs) beta2 = beta2gs;
        E_defo1 = frldm( z1, N1, beta1 ) - frldm( z1, N1, beta1gs );
        E_defo2 = frldm( z2, N2, beta2 ) - frldm( z2, N2, beta2gs );
        E1exc = E_scission_pre * a1 / A + E_defo1;
        E2exc = E_scission_pre * a2 / A + E_defo2;
      }


//  pre-neutron-emission total kinetic energy */
    TKER = ( z1 * z2 * 1.440 ) /
           ( R0 * std::pow(a1,0.333330) * (1.0 + 2.0/3.0 * beta1 ) +
             R0 * std::pow(a2,0.333330) * (1.0 + 2.0/3.0 * beta2 ) + 2.0 );
//  Pre-neutron-emission kinetic energies of the fragments */
    EkinR1 = TKER * a2 / A;
    EkinR2 = TKER * a1 / A;
    v1 = std::sqrt(EkinR1/a1) * 1.3887;
    v2 = std::sqrt(EkinR2/a2) * 1.3887;

//  Extracted from Lang et al. Nucl. Phys. A 345 (1980) 34 */
    E1exc_sigma = 5.50;
    E2exc_sigma = 5.50;

    fis987:
    //e1 = E1exc+G4AblaRandom::gaus(E1exc_sigma);//
    e1 = gausshaz(0,E1exc,E1exc_sigma);
    if(e1<0.)goto fis987;
    fis988:
    //e2 = E2exc+G4AblaRandom::gaus(E2exc_sigma);//
    e2 = gausshaz(0,E2exc,E2exc_sigma);
    if(e2<0.)goto fis988;

    (*NbLam0_par) = NbLam0;
    return;
}


void G4Abla::even_odd(G4double r_origin,G4double r_even_odd,G4int &i_out)     
{
  // Procedure to calculate I_OUT from R_IN in a way that
  // on the average a flat distribution in R_IN results in a
  // fluctuating distribution in I_OUT with an even-odd effect as
  // given by R_EVEN_ODD

  //     /* ------------------------------------------------------------ */
  //     /* EXAMPLES :                                                   */
  //     /* ------------------------------------------------------------ */
  //     /*    If R_EVEN_ODD = 0 :                                       */
  //     /*           CEIL(R_IN)  ----                                   */
  //     /*                                                              */
  //     /*              R_IN ->                                         */
  //     /*            (somewhere in between CEIL(R_IN) and FLOOR(R_IN)) */                                            */
  //     /*                                                              */
  //     /*           FLOOR(R_IN) ----       --> I_OUT                   */
  //     /* ------------------------------------------------------------ */
  //     /*    If R_EVEN_ODD > 0 :                                       */
  //     /*      The interval for the above treatment is                 */
  //     /*         larger for FLOOR(R_IN) = even and                    */
  //     /*         smaller for FLOOR(R_IN) = odd                        */
  //     /*    For R_EVEN_ODD < 0 : just opposite treatment              */
  //     /* ------------------------------------------------------------ */

  //     /* ------------------------------------------------------------ */
  //     /* On input:   R_ORIGIN    nuclear charge (real number)         */
  //     /*             R_EVEN_ODD  requested even-odd effect            */
  //     /* Intermediate quantity: R_IN = R_ORIGIN + 0.5                 */
  //     /* On output:  I_OUT       nuclear charge (integer)             */
  //     /* ------------------------------------------------------------ */

  //      G4double R_ORIGIN,R_IN,R_EVEN_ODD,R_REST,R_HELP;
  G4double r_in = 0.0, r_rest = 0.0, r_help = 0.0;
  G4double r_floor = 0.0;
  G4double r_middle = 0.0;
  //      G4int I_OUT,N_FLOOR;
  G4int n_floor = 0;

  r_in = r_origin + 0.5;
  r_floor = (G4double)((G4int)(r_in));
  if (r_even_odd < 0.001) {
    i_out = (G4int)(r_floor);
  } 
  else {
    r_rest = r_in - r_floor;
    r_middle = r_floor + 0.5;
    n_floor = (G4int)(r_floor);
    if (n_floor%2 == 0) {
      // even before modif.
      r_help = r_middle + (r_rest - 0.5) * (1.0 - r_even_odd);
    } 
    else {
      // odd before modification
      r_help = r_middle + (r_rest - 0.5) * (1.0 + r_even_odd);
    }
    i_out = (G4int)(r_help);
  }
}

double G4Abla::umass(G4double z,G4double n,G4double beta)
{
  // liquid-drop mass, Myers & Swiatecki, Lysekil, 1967
  // pure liquid drop, without pairing and shell effects

  // On input:    Z     nuclear charge of nucleus
  //              N     number of neutrons in nucleus
  //              beta  deformation of nucleus
  // On output:   binding energy of nucleus

  G4double a = 0.0, fumass = 0.0;
  G4double alpha = 0.0;
  G4double xcom = 0.0, xvs = 0.0, xe = 0.0;
  const G4double pi = 3.1416;

  a = n + z;
  alpha = ( std::sqrt(5.0/(4.0*pi)) ) * beta;
  
  xcom = 1.0 - 1.7826 * ((a - 2.0*z)/a)*((a - 2.0*z)/a);
  // factor for asymmetry dependence of surface and volume term
  xvs = - xcom * ( 15.4941 * a - 
		   17.9439 * std::pow(a,2.0/3.0) * (1.0+0.4*alpha*alpha) );
  // sum of volume and surface energy
  xe = z*z * (0.7053/(std::pow(a,1.0/3.0)) * (1.0-0.2*alpha*alpha) - 1.1529/a);
  fumass = xvs + xe;
  
  return fumass;
}


double G4Abla::ecoul(G4double z1,G4double n1,G4double beta1,G4double z2,G4double n2,G4double beta2,G4double d)
{
  // Coulomb potential between two nuclei
  // surfaces are in a distance of d
  // in a tip to tip configuration

  // approximate formulation
  // On input: Z1      nuclear charge of first nucleus
  //           N1      number of neutrons in first nucleus
  //           beta1   deformation of first nucleus
  //           Z2      nuclear charge of second nucleus
  //           N2      number of neutrons in second nucleus
  //           beta2   deformation of second nucleus
  //           d       distance of surfaces of the nuclei

  //      G4double Z1,N1,beta1,Z2,N2,beta2,d,ecoul;
  G4double fecoul = 0;
  G4double dtot = 0;
  const G4double r0 = 1.16;

  dtot = r0 * ( std::pow((z1+n1),1.0/3.0) * (1.0+0.6666667*beta1)
		+ std::pow((z2+n2),1.0/3.0) * (1.0+0.6666667*beta2) ) + d;
  fecoul = z1 * z2 * 1.44 / dtot;

  return fecoul;
}


 G4double G4Abla::Uwash(G4double E, G4double Ecrit,G4double Freduction,G4double gamma){
        // E       excitation energy 
        // Ecrit   critical pairing energy 
        // Freduction  reduction factor for shell washing in superfluid region
        G4double R_wash,uwash;
        if(E < Ecrit)
          R_wash = std::exp(-E * Freduction * gamma);
        else
          R_wash = std::exp(- Ecrit * Freduction * gamma -(E-Ecrit) * gamma);
        
        uwash = R_wash;
 return uwash;
}


G4double G4Abla::frldm(G4double z,G4double n,G4double beta){

//     Liquid-drop mass, Myers & Swiatecki, Lysekil, 1967
//     pure liquid drop, without pairing and shell effects
//
//     On input:    Z     nuclear charge of nucleus
//                  N     number of neutrons in nucleus
//                  beta  deformation of nucleus
//     On output:   binding energy of nucleus
// The idea is to use FRLDM model for beta=0 and using Lysekil
// model to get the deformation energy

      G4double a;
      a = n + z;
      return eflmac_profi(a,z) + umass(z,n,beta) - umass(z,n,0.0);
}


//**********************************************************************
// *
// * this function will calculate the liquid-drop nuclear mass for spheri
// * configuration according to the preprint NUCLEAR GROUND-STATE
// * MASSES and DEFORMATIONS by P. M"oller et al. from August 16, 1993 p.
// * All constants are taken from this publication for consistency.
// *
// * Parameters:
// *   a:    nuclear mass number
// *   z:    nuclear charge
// **********************************************************************


G4double G4Abla::eflmac_profi(G4double ia, G4double iz)
{
  // CHANGED TO CALCULATE TOTAL BINDING ENERGY INSTEAD OF MASS EXCESS.     
  // SWITCH FOR PAIRING INCLUDED AS WELL.                                  
  // BINDING = EFLMAC(IA,IZ,0,OPTSHP)                                      
  // FORTRAN TRANSCRIPT OF /U/GREWE/LANG/EEX/FRLDM.C                       
  // A.J. 15.07.96                                                         

  // this function will calculate the liquid-drop nuclear mass for spheri
  // configuration according to the preprint NUCLEAR GROUND-STATE        
  // MASSES and DEFORMATIONS by P. M"oller et al. from August 16, 1993 p.
  // All constants are taken from this publication for consistency.      

  // Parameters:                                                         
  // a:    nuclear mass number                                         
  // z:    nuclear charge                                     

  G4double eflmacResult = 0.0;

  G4int in = 0;
  G4double z = 0.0, n = 0.0, a = 0.0, av = 0.0, as = 0.0;
  G4double a0 = 0.0, c1 = 0.0, c4 = 0.0, b1 = 0.0, b3 = 0.0;
  G4double ff = 0.0, ca = 0.0, w = 0.0, efl = 0.0; 
  G4double r0 = 0.0, kf = 0.0, ks = 0.0;
  G4double kv = 0.0, rp = 0.0, ay = 0.0, aden = 0.0, x0 = 0.0, y0 = 0.0;
  G4double esq = 0.0, ael = 0.0, i = 0.0;
  G4double pi = 3.141592653589793238e0;

  // fundamental constants
  // electronic charge squared
  esq = 1.4399764;

  // constants from considerations other than nucl. masses
  // electronic binding
  ael = 1.433e-5;

  // proton rms radius
  rp  = 0.8;

  // nuclear radius constant
  r0  = 1.16;

  // range of yukawa-plus-expon. potential
  ay  = 0.68;

  // range of yukawa function used to generate                          
  // nuclear charge distribution
  aden= 0.70;

  // wigner constant
  w   = 30.0;

  // adjusted parameters
  // volume energy
  av  = 16.00126;

  // volume asymmetry
  kv  =  1.92240;

  // surface energy
  as  = 21.18466;

  // surface asymmetry
  ks  =  2.345;
  // a^0 constant
  a0  =  2.615;

  // charge asymmetry
  ca  =  0.10289;

  z   = G4double(iz);
  a   = G4double(ia);
  in  = ia - iz;                                                       
  n   = G4double(in);

  
  c1  = 3.0/5.0*esq/r0;
  c4  = 5.0/4.0*std::pow((3.0/(2.0*pi)),(2.0/3.0)) * c1;
  kf  = std::pow((9.0*pi*z/(4.0*a)),(1.0/3.0))/r0;
  
  ff = -1.0/8.0*rp*rp*esq/std::pow(r0,3) * (145.0/48.0 - 327.0/2880.0*std::pow(kf,2) * std::pow(rp,2) + 1527.0/1209600.0*std::pow(kf,4) * std::pow(rp,4));

  i   = (n-z)/a;

  x0  = r0 * std::pow(a,(1.0/3.0)) / ay;
  y0  = r0 * std::pow(a,(1.0/3.0)) / aden;

  b1  = 1.0 - 3.0/(std::pow(x0,2)) + (1.0 + x0) * (2.0 + 3.0/x0 + 3.0/std::pow(x0,2)) * std::exp(-2.0*x0);

  b3  = 1.0 - 5.0/std::pow(y0,2) * (1.0 - 15.0/(8.0*y0) + 21.0/(8.0 * std::pow(y0,3))
			       - 3.0/4.0 * (1.0 + 9.0/(2.0*y0) + 7.0/std::pow(y0,2)
					    + 7.0/(2.0 * std::pow(y0,3))) * std::exp(-2.0*y0));

  // now calculation of total binding energy                  

  efl = -1.0 * av*(1.0 - kv*i*i)*a + as*(1.0 - ks*i*i)*b1 * std::pow(a,(2.0/3.0)) + a0
    + c1*z*z*b3/std::pow(a,(1.0/3.0)) - c4*std::pow(z,(4.0/3.0))/std::pow(a,(1.e0/3.e0))
    + ff*std::pow(z,2)/a -ca*(n-z) - ael * std::pow(z,(2.39e0));

  efl = efl + w*utilabs(i);

  eflmacResult = efl;

  return eflmacResult;
}
//
//
//
void G4Abla::unstable_nuclei(G4int AFP,G4int ZFP,G4int *AFPNEW,G4int *ZFPNEW,G4int &IOUNSTABLE,G4double VX,G4double VY,G4double VZ,G4double *VP1X,G4double *VP1Y,G4double *VP1Z,G4double BU_TAB_TEMP[200][6],G4int *ILOOP){
//
      G4int INMIN,INMAX,NDIF=0,IMEM;
      G4int NEVA=0,PEVA=0;
      G4double   VP2X,VP2Y,VP2Z;

      *AFPNEW = AFP;
      *ZFPNEW = ZFP;
      IOUNSTABLE = 0;
      *ILOOP = 0;
      IMEM = 0;
      for(G4int i=0;i<200;i++){
      BU_TAB_TEMP[i][0] = 0.0;
      BU_TAB_TEMP[i][1] = 0.0;
      BU_TAB_TEMP[i][2] = 0.0;
      BU_TAB_TEMP[i][3] = 0.0;
      BU_TAB_TEMP[i][4] = 0.0;
      //BU_TAB_TEMP[i][5] = 0.0;
      }
      *VP1X = 0.0;
      *VP1Y = 0.0;
      *VP1Z = 0.0;

      if(AFP==0 && ZFP==0){
//       PRINT*,'UNSTABLE NUCLEI, AFP=0, ZFP=0'
       return;
      }
      if((AFP==1 && ZFP==0) || 
         (AFP==1 && ZFP==1) || 
         (AFP==2 && ZFP==1) || 
         (AFP==3 && ZFP==1) || 
         (AFP==3 && ZFP==2) || 
         (AFP==4 && ZFP==2) || 
         (AFP==6 && ZFP==2) || 
         (AFP==8 && ZFP==2)
       ){
       *VP1X = VX;
       *VP1Y = VY;
       *VP1Z = VZ;
       return;
      }

      if ((AFP-ZFP)==0 && ZFP>1){
        for(G4int I = 0;I<=AFP-2;I++){
         unstable_tke(G4double(AFP-I),G4double(AFP-I),G4double(AFP-I-1),G4double(AFP-I-1),VX,VY,VZ,
           &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
        }
        // PEVA = PEVA + ZFP - 1;
         AFP = 1;
         ZFP = 1;
         IOUNSTABLE = 1;
      }
//
//*** Find the limits nucleus is bound :
      isostab_lim(ZFP,&INMIN,&INMAX);
      NDIF = AFP - ZFP;
      if(NDIF<INMIN){
// Proton unbound
        IOUNSTABLE = 1;
       for(G4int I = 1;I<=10; I++){
        isostab_lim(ZFP-I,&INMIN,&INMAX);
         if(INMIN<=NDIF){
         IMEM = I;
         ZFP = ZFP - I;
         AFP = ZFP + NDIF;
         PEVA = I;
         goto u10;
         }
       }
//
      u10:
       for(G4int I = 0;I< IMEM;I++){
         unstable_tke(G4double(NDIF+ZFP+IMEM-I),
         G4double(ZFP+IMEM-I),
         G4double(NDIF+ZFP+IMEM-I-1),
         G4double(ZFP+IMEM-I-1),
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[I+1+*ILOOP][0] = 1.0;
         BU_TAB_TEMP[I+1+*ILOOP][1] = 1.0;
         BU_TAB_TEMP[I+1+*ILOOP][2] = VP2X;
         BU_TAB_TEMP[I+1+*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[I+1+*ILOOP][4] = VP2Z;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
       }
         *ILOOP = *ILOOP + IMEM;

      }
      if(NDIF>INMAX){
// Neutron unbound
        NEVA = NDIF - INMAX;
        AFP = ZFP + INMAX;
        IOUNSTABLE = 1;
         for(G4int I = 0;I<NEVA;I++){
         unstable_tke(G4double(ZFP+NDIF-I),
         G4double(ZFP),
         G4double(ZFP+NDIF-I-1),
         G4double(ZFP),
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);

         BU_TAB_TEMP[*ILOOP][0] = 0.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }
      }

         if ((AFP>=2) && (ZFP==0)){
         for(G4int I = 0;I<= AFP-2;I++){
         unstable_tke(G4double(AFP-I),G4double(ZFP),
         G4double(AFP-I-1),G4double(ZFP),
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);

         BU_TAB_TEMP[*ILOOP][0] = 0.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }

	 //NEVA = NEVA + (AFP - 1);
          AFP = 1;
          ZFP = 0;
          IOUNSTABLE = 1;
         }
         if (AFP<ZFP){
          std::cout << "WARNING - BU UNSTABLE: AF < ZF" << std::endl;
          AFP = 0;
          ZFP = 0;
          IOUNSTABLE = 1;
         }
         if ((AFP>=4) && (ZFP==1)){
// Heavy residue is treated as 3H and the rest of mass is emitted as neutrons:
         for(G4int I = 0; I<AFP-3;I++){
         unstable_tke(G4double(AFP-I),G4double(ZFP),
         G4double(AFP-I-1),G4double(ZFP),
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);

         BU_TAB_TEMP[*ILOOP][0] = 0.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }

	 // NEVA = NEVA + (AFP - 3);
         AFP = 3;
         ZFP = 1;
         IOUNSTABLE = 1;
         }

         if ((AFP==4) && (ZFP==3)){
// 4Li -> 3He + p  ->
         AFP = 3;
         ZFP = 2;
         //PEVA = PEVA + 1;
         IOUNSTABLE = 1;
         unstable_tke(4.0,3.0,3.0,2.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);

         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         }
         if ((AFP==5) && (ZFP==2)){
// 5He -> 4He + n  ->
         AFP = 4;
         ZFP = 2;
         //NEVA = NEVA + 1;
         IOUNSTABLE = 1;
         unstable_tke(5.0,2.0,4.0,2.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 0.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         }

         if ((AFP==5) && (ZFP==3)){
// 5Li -> 4He + p
         AFP = 4;
         ZFP = 2;
         //PEVA = PEVA + 1;
         IOUNSTABLE = 1;
         unstable_tke(5.0,3.0,4.0,2.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         }

         if ((AFP==6) && (ZFP==4)){
// 6Be -> 4He + 2p (velocity in two steps: 6Be->5Li->4He)
         AFP = 4;
         ZFP = 2;
         //PEVA = PEVA + 2;
         IOUNSTABLE = 1;
// 6Be -> 5Li + p
         unstable_tke(6.0,4.0,5.0,3.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;

// 5Li -> 4He + p
         unstable_tke(5.0,3.0,4.0,2.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         }
         if ((AFP==7)&&(ZFP==2)){
// 7He -> 6He + n
         AFP = 6;
         ZFP = 2;
         //NEVA = NEVA + 1;
         IOUNSTABLE = 1;
         unstable_tke(7.0,2.0,6.0,2.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 0.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         }

         if ((AFP==7) && (ZFP==5)){
// 7B -> 6Be + p -> 4He + 3p
         for(int I = 0; I<= AFP-5;I++){
         unstable_tke(double(AFP-I),double(ZFP-I),
         double(AFP-I-1),double(ZFP-I-1),
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }

         AFP = 4;
         ZFP = 2;
         //PEVA = PEVA + 3;
         IOUNSTABLE = 1;
         }
         if ((AFP==8) && (ZFP==4)){
// 8Be  -> 4He + 4He
          AFP = 4;
          ZFP = 2;
         IOUNSTABLE = 1;
         unstable_tke(8.0,4.0,4.0,2.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 2.0;
         BU_TAB_TEMP[*ILOOP][1] = 4.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         }
         if ((AFP==8) && (ZFP==6)){
// 8C  -> 2p + 6Be
          AFP = 6;
          ZFP = 4;
          //PEVA = PEVA + 2;
         IOUNSTABLE = 1;

         unstable_tke(8.0,6.0,7.0,5.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;

         unstable_tke(7.0,5.0,6.0,4.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }

         if((AFP==9) && (ZFP==2)){
// 9He -> 8He + n
           AFP = 8;
           ZFP = 2;
           //NEVA = NEVA + 1;
         IOUNSTABLE = 1;

         unstable_tke(9.0,2.0,8.0,2.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 0.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }

         if((AFP==9) && (ZFP==5)){
// 9B -> 4He + 4He + p  ->
          AFP = 4;
          ZFP = 2;
          //PEVA = PEVA + 1;
         IOUNSTABLE = 1;
         unstable_tke(9.0,5.0,8.0,4.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;

         unstable_tke(8.0,4.0,4.0,2.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 2.0;
         BU_TAB_TEMP[*ILOOP][1] = 4.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }

         if((AFP==10) && (ZFP==2)){
// 10He -> 8He + 2n
           AFP = 8;
           ZFP = 2;
           //NEVA = NEVA + 2;
         IOUNSTABLE = 1;
// 10He -> 9He + n
         unstable_tke(10.0,2.0,9.0,2.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 0.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;

// 9He -> 8He + n
         unstable_tke(9.0,2.0,8.0,2.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 0.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }
         if ((AFP==10) && (ZFP==3)){
// 10Li -> 9Li + n  ->
          AFP = 9;
          ZFP = 3;
          //NEVA = NEVA + 1;
         IOUNSTABLE = 1;
         unstable_tke(10.0,3.0,9.0,3.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 0.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }
         if ((AFP==10) && (ZFP==7)){
// 10N -> 9C + p  ->
          AFP = 9;
          ZFP = 6;
          //PEVA = PEVA + 1;
         IOUNSTABLE = 1;
         unstable_tke(10.0,7.0,9.0,6.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }

         if((AFP==11) && (ZFP==7)){
// 11N -> 10C + p  ->
          AFP = 10;
          ZFP = 6;
          //PEVA = PEVA + 1;
         IOUNSTABLE = 1;
         unstable_tke(11.0,7.0,10.0,6.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }
         if ((AFP==12) && (ZFP==8)){
// 12O -> 10C + 2p  ->
          AFP = 10;
          ZFP = 6;
          //PEVA = PEVA + 2;
         IOUNSTABLE = 1;

         unstable_tke(12.0,8.0,11.0,7.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;

         unstable_tke(11.0,7.0,10.0,6.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }
         if ((AFP==15) && (ZFP==9)){
// 15F -> 14O + p  ->
          AFP = 14;
          ZFP = 8;
          //PEVA = PEVA + 1;
         IOUNSTABLE = 1;
         unstable_tke(15.0,9.0,14.0,8.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }

         if ((AFP==16) && (ZFP==9)){
// 16F -> 15O + p  ->
          AFP = 15;
          ZFP = 8;
          //PEVA = PEVA + 1;
         IOUNSTABLE = 1;
         unstable_tke(16.0,9.0,15.0,8.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }

         if ((AFP==16) && (ZFP==10)){
// 16Ne -> 14O + 2p  ->
          AFP = 14;
          ZFP = 8;
          //PEVA = PEVA + 2;
         IOUNSTABLE = 1;
         unstable_tke(16.0,10.0,15.0,9.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;

         unstable_tke(15.0,9.0,14.0,8.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }
         if((AFP==18) && (ZFP==11)){
// 18Na -> 17Ne + p  ->
          AFP = 17;
          ZFP = 10;
          //PEVA = PEVA + 1;
         IOUNSTABLE = 1;
         unstable_tke(18.0,11.0,17.0,10.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }
         if((AFP==19) && (ZFP==11)){
// 19Na -> 18Ne + p  ->
         AFP = 18;
         ZFP = 10;
         //PEVA = PEVA + 1;
         IOUNSTABLE = 1;
         unstable_tke(19.0,11.0,18.0,10.0,
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }
         if (ZFP>=4 && (AFP-ZFP)==1){
// Heavy residue is treated as 3He
           NEVA = AFP - 3;
           PEVA = ZFP - 2;

         for(G4int I = 0;I< NEVA;I++){
        unstable_tke(G4double(AFP-I),G4double(ZFP),
         G4double(AFP-I-1),G4double(ZFP),
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 0.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }
        for( G4int I = 0;I<PEVA;I++){
        unstable_tke(G4double(AFP-NEVA-I),G4double(ZFP-I),
         G4double(AFP-NEVA-I-1),G4double(ZFP-I-1),
         VX,VY,VZ,
         &(*VP1X),&(*VP1Y),&(*VP1Z),&VP2X,&VP2Y,&VP2Z);
         BU_TAB_TEMP[*ILOOP][0] = 1.0;
         BU_TAB_TEMP[*ILOOP][1] = 1.0;
         BU_TAB_TEMP[*ILOOP][2] = VP2X;
         BU_TAB_TEMP[*ILOOP][3] = VP2Y;
         BU_TAB_TEMP[*ILOOP][4] = VP2Z;
         *ILOOP = *ILOOP + 1;
         VX = *VP1X;
         VY = *VP1Y;
         VZ = *VP1Z;
         }

         AFP = 3;
         ZFP = 2;
         IOUNSTABLE = 1;
         }
//
      *AFPNEW = AFP;
      *ZFPNEW = ZFP;
      return;
}

//
//
void G4Abla::unstable_tke(G4double ain,G4double zin,G4double anew,G4double znew,G4double vxin,G4double vyin,G4double vzin,G4double *v1x,G4double *v1y,G4double *v1z,G4double *v2x,G4double *v2y,G4double *v2z){
//
      G4double EKIN_P1=0.,ekin_tot=0.;
      G4double PX1,PX2,PY1,PY2,PZ1,PZ2,PTOT;
      G4double RNDT,CTET1,STET1,RNDP,PHI1,ETOT_P1,ETOT_P2;
      G4double MASS,MASS1,MASS2;
      G4double vxout=0.,vyout=0.,vzout=0.;
      G4int iain,izin,ianew,iznew,inin,innew;
//
      G4double C = 29.97924580;//         cm/ns
      G4double AMU = 931.4940; //         MeV/C^2
//
      iain = idnint(ain);
      izin = idnint(zin);
      inin = iain - izin;
      ianew = idnint(anew);
      iznew = idnint(znew);
      innew = ianew - iznew;
      //
      if(ain==0)return;
      //
      if(izin>12){
      mglms(ain,zin,3,&MASS);
      mglms(anew,znew,3,&MASS1);
      mglms(ain-anew,zin-znew,3,&MASS2);
      ekin_tot = MASS-MASS1-MASS2;
      }else{
    //  ekin_tot = MEXP(ININ,IZIN)-(MEXP(INNEW,IZNEW)+MEXP(ININ-INNEW,IZIN-IZNEW));
      ekin_tot = masses->massexp[inin][izin]-(masses->massexp[innew][iznew]+masses->massexp[inin-innew][izin-iznew]);
      if(izin>12)std::cout << "*** ZIN > 12 ***" << izin << std::endl;
      }

      if( ekin_tot<0.00 ){
//         if( iain.ne.izin .and. izin.ne.0 ){
//            print *,"Negative Q-value in UNSTABLE_TKE"
//            print *,"ekin_tot=",ekin_tot
//            print *,"ain,zin=",ain,zin,MEXP(ININ,IZIN)
//            print *,"anew,znew=",anew,znew,MEXP(INNEW,IZNEW)
//            print *
//          }
        ekin_tot=0.0;
      }
//
      EKIN_P1 = ekin_tot*(ain-anew)/ ain;
      ETOT_P1 = EKIN_P1 + anew * AMU;
      PTOT = anew*AMU*std::sqrt((EKIN_P1/(anew*AMU)+1.0)*(EKIN_P1/(anew*AMU)+1.0)-1.0);  // MeV/C
//
      RNDT = G4AblaRandom::flat();
      CTET1 = 2.0*RNDT-1.0;
      STET1 = std::sqrt(1.0-CTET1*CTET1);
      RNDP = G4AblaRandom::flat();
      PHI1 = RNDP*2.0*3.141592654;
      PX1 = PTOT * STET1*std::cos(PHI1);
      PY1 = PTOT * STET1*std::sin(PHI1);
      PZ1 = PTOT * CTET1;
      *v1x = C * PX1 / ETOT_P1;
      *v1y = C * PY1 / ETOT_P1;
      *v1z = C * PZ1 / ETOT_P1;
      lorentz_boost(vxin,vyin,vzin,*v1x,*v1y,*v1z,&vxout,&vyout,&vzout);
      *v1x = vxout;
      *v1y = vyout;
      *v1z = vzout;
//
      PX2 = - PX1;
      PY2 = - PY1;
      PZ2 = - PZ1;
      ETOT_P2 = (ekin_tot - EKIN_P1) + (ain-anew) * AMU;
      *v2x = C * PX2 / ETOT_P2;
      *v2y = C * PY2 / ETOT_P2;
      *v2z = C * PZ2 / ETOT_P2;
      lorentz_boost(vxin,vyin,vzin,*v2x,*v2y,*v2z,&vxout,&vyout,&vzout);
      *v2x = vxout;
      *v2y = vyout;
      *v2z = vzout;
//
   return;
}
//
//**************************************************************************
//
void G4Abla::lorentz_boost(G4double VXRIN,G4double VYRIN,G4double VZRIN,G4double VXIN,G4double VYIN,G4double VZIN,G4double *VXOUT,G4double *VYOUT,G4double *VZOUT){
//
// Calculate velocities of a given fragment from frame 1 into frame 2.
// Frame 1 is moving with velocity v=(vxr,vyr,vzr) relative to frame 2.
// Velocity of the fragment in frame 1 -> vxin,vyin,vzin
// Velocity of the fragment in frame 2 -> vxout,vyout,vzout
//
      G4double  VXR,VYR,VZR;
      G4double  GAMMA,VR,C,CC,DENO,VXNOM,VYNOM,VZNOM;
//
      C = 29.9792458;        // cm/ns
      CC = C*C;
//
// VXR,VYR,VZR are velocities of frame 1 relative to frame 2; to go from 1 to 2
// we need to multiply them by -1
      VXR = -1.0 * VXRIN;
      VYR = -1.0 * VYRIN;
      VZR = -1.0 * VZRIN;
//
      VR = std::sqrt(VXR*VXR + VYR*VYR + VZR*VZR);
      if(VR<1e-9){
         *VXOUT = VXIN;
         *VYOUT = VYIN;
         *VZOUT = VZIN;
         return;
      }
      GAMMA = 1.0/std::sqrt(1.0 - VR*VR/CC);
      DENO = 1.0 - VXR*VXIN/CC - VYR*VYIN/CC - VZR*VZIN/CC;

// X component
      VXNOM = -GAMMA*VXR + (1.0+(GAMMA-1.0)*VXR*VXR/(VR*VR))*VXIN + (GAMMA-1.0)*VXR*VYR/(VR*VR)*VYIN + (GAMMA-1.0)*VXR*VZR/(VR*VR)*VZIN;

      *VXOUT = VXNOM / (GAMMA * DENO);

// Y component
      VYNOM = -GAMMA*VYR + (1.0+(GAMMA-1.0)*VYR*VYR/(VR*VR))*VYIN + (GAMMA-1.0)*VXR*VYR/(VR*VR)*VXIN + (GAMMA-1.0)*VYR*VZR/(VR*VR)*VZIN;

      *VYOUT = VYNOM / (GAMMA * DENO);

// Z component
      VZNOM = -GAMMA*VZR + (1.0+(GAMMA-1.0)*VZR*VZR/(VR*VR))*VZIN + (GAMMA-1.0)*VXR*VZR/(VR*VR)*VXIN + (GAMMA-1.0)*VYR*VZR/(VR*VR)*VYIN;

      *VZOUT = VZNOM / (GAMMA * DENO);

     return;
}

void G4Abla::fission(G4double AF,G4double ZF,G4double EE,G4double JPRF,
        G4double *VX1_FISSION_par,G4double *VY1_FISSION_par,G4double *VZ1_FISSION_par,
        G4double *VX2_FISSION_par,G4double *VY2_FISSION_par,G4double *VZ2_FISSION_par,
        G4int *ZFP1,G4int *AFP1,G4int *SFP1, G4int *ZFP2,G4int *AFP2,G4int *SFP2,G4int *imode_par, 
        G4double *VX_EVA_SC_par, G4double *VY_EVA_SC_par, G4double *VZ_EVA_SC_par,
        G4double EV_TEMP[200][6],G4int *IEV_TAB_FIS_par, G4int *NbLam0_par){
///
       G4double EFF1=0.,EFF2=0.,VFF1=0.,VFF2=0.,
                 AF1=0.,ZF1=0.,AF2=0.,ZF2=0.,
                 AFF1=0.,ZFF1=0.,AFF2=0.,ZFF2=0.,
                 vz1_eva=0., vx1_eva=0.,vy1_eva=0.,
                 vz2_eva=0., vx2_eva=0.,vy2_eva=0.,
                 vx_eva_sc=0.,vy_eva_sc=0.,vz_eva_sc=0.,
                 VXOUT=0.,VYOUT=0.,VZOUT=0.,
                 VX2OUT=0.,VY2OUT=0.,VZ2OUT=0.;
       G4int IEV_TAB_FIS=0,IEV_TAB_TEMP=0;
       G4double EV_TEMP1[200][6], EV_TEMP2[200][6],mtota=0.;
       G4int inttype = 0,inum=0;
       IEV_TAB_SSC=0;
       (*imode_par)=0;
       G4int NbLam0= (*NbLam0_par);

       for(G4int I1=0;I1<200;I1++)
       for(G4int I2=0;I2<6;I2++){
       EV_TEMP[I1][I2] = 0.0;
       EV_TEMP1[I1][I2] = 0.0;
       EV_TEMP2[I1][I2] = 0.0;
       }

       G4double et = EE - JPRF * JPRF * 197. * 197./(2.*0.4*931.*std::pow(AF,5.0/3.0)*1.16*1.16);

       fissionDistri(AF,ZF,et,AF1,ZF1,EFF1,VFF1,AF2,ZF2,EFF2,VFF2,
                     vx_eva_sc,vy_eva_sc,vz_eva_sc,&NbLam0);

//  Lambda particles 
      G4int NbLam1=0;
      G4int NbLam2=0;
      G4double pbH = (AF1 - ZF1) / (AF1 - ZF1 + AF2 - ZF2);  
      for(G4int i=0;i<NbLam0;i++){
       if(G4AblaRandom::flat()<pbH){
        NbLam1++;
       }else{
        NbLam2++;
       }
      }
//     Copy of the evaporated particles from saddle to scission
       for(G4int IJ = 0; IJ< IEV_TAB_SSC;IJ++){
             EV_TEMP[IJ][0] = EV_TAB_SSC[IJ][0];
             EV_TEMP[IJ][1] = EV_TAB_SSC[IJ][1];
             EV_TEMP[IJ][2] = EV_TAB_SSC[IJ][2];
             EV_TEMP[IJ][3] = EV_TAB_SSC[IJ][3];
             EV_TEMP[IJ][4] = EV_TAB_SSC[IJ][4];
             EV_TEMP[IJ][5] = EV_TAB_SSC[IJ][5];
       }
       IEV_TAB_FIS = IEV_TAB_FIS + IEV_TAB_SSC;

//    Velocities
      G4double VZ1_FISSION = (2.0 * G4AblaRandom::flat() - 1.0) * VFF1;
      G4double VPERP1 = std::sqrt(VFF1*VFF1 - VZ1_FISSION*VZ1_FISSION);
      G4double ALPHA1 = G4AblaRandom::flat() * 2. * 3.142;
      G4double VX1_FISSION = VPERP1 * std::sin(ALPHA1);
      G4double VY1_FISSION = VPERP1 * std::cos(ALPHA1);
      G4double VX2_FISSION = - VX1_FISSION / VFF1 * VFF2;
      G4double VY2_FISSION = - VY1_FISSION / VFF1 * VFF2;
      G4double VZ2_FISSION = - VZ1_FISSION / VFF1 * VFF2;
//
// Fission fragment 1
      if( (ZF1<=0.0) || (AF1<=0.0) || (AF1<ZF1) ){
       std::cout << "F1 unphysical: "<<ZF<< " "<<AF<< " "<<EE<< " "<<ZF1<< " "<<AF1 << std::endl;
      }else{
// fission and IMF emission are not allowed
     opt->optimfallowed = 0; //  IMF is not allowed
     fiss->ifis = 0;         //  fission is not allowed
     gammaemission=1;
     G4int FF11=0, FIMF11=0;
     G4double ZIMFF1=0., AIMFF1=0.,TKEIMF1=0.,JPRFOUT=0.;
//
     evapora(ZF1,AF1,&EFF1,0., &ZFF1, &AFF1, &mtota, &vz1_eva, &vx1_eva,&vy1_eva, &FF11, &FIMF11, &ZIMFF1, &AIMFF1,&TKEIMF1, &JPRFOUT, &inttype, &inum,EV_TEMP1,&IEV_TAB_TEMP,&NbLam1);

               for(G4int IJ = 0; IJ< IEV_TAB_TEMP;IJ++){ 
               EV_TEMP[IJ+IEV_TAB_FIS][0] = EV_TEMP1[IJ][0];
               EV_TEMP[IJ+IEV_TAB_FIS][1] = EV_TEMP1[IJ][1];
// Lorentz kinematics
//               EV_TEMP(IJ+IEV_TAB,3) = EV_TEMP(IJ,3) + VX_PREF
//               EV_TEMP(IJ+IEV_TAB,4) = EV_TEMP(IJ,4) + VY_PREF
//               EV_TEMP(IJ+IEV_TAB,5) = EV_TEMP(IJ,5) + VZ_PREF
// Lorentz transformation
               lorentz_boost(VX1_FISSION,VY1_FISSION,VZ1_FISSION,
                EV_TEMP1[IJ][2],EV_TEMP1[IJ][3],EV_TEMP1[IJ][4],
                &VXOUT,&VYOUT,&VZOUT);
               lorentz_boost(vx_eva_sc,vy_eva_sc,vz_eva_sc,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
               EV_TEMP[IJ+IEV_TAB_FIS][2] = VX2OUT;
               EV_TEMP[IJ+IEV_TAB_FIS][3] = VY2OUT;
               EV_TEMP[IJ+IEV_TAB_FIS][4] = VZ2OUT;
               //
               }
               IEV_TAB_FIS = IEV_TAB_FIS + IEV_TAB_TEMP;

      }
//
// Fission fragment 2
      if( (ZF2<=0.0) || (AF2<=0.0) || (AF2<ZF2) ){
       std::cout << "F2 unphysical: "<<ZF<< " "<<AF<< " "<<EE<< " "<<ZF2<< " "<<AF2 << std::endl;
      }else{
// fission and IMF emission are not allowed
     opt->optimfallowed = 0; //  IMF is not allowed
     fiss->ifis = 0;         //  fission is not allowed
     gammaemission=1;
     G4int FF22=0, FIMF22=0;
     G4double ZIMFF2=0., AIMFF2=0.,TKEIMF2=0.,JPRFOUT=0.;
//
     evapora(ZF2,AF2,&EFF2,0., &ZFF2, &AFF2, &mtota, &vz2_eva, &vx2_eva,&vy2_eva, &FF22, &FIMF22, &ZIMFF2, &AIMFF2,&TKEIMF2, &JPRFOUT, &inttype, &inum,EV_TEMP2,&IEV_TAB_TEMP,&NbLam2);

               for(G4int IJ = 0; IJ< IEV_TAB_TEMP;IJ++){
               EV_TEMP[IJ+IEV_TAB_FIS][0] = EV_TEMP2[IJ][0];
               EV_TEMP[IJ+IEV_TAB_FIS][1] = EV_TEMP2[IJ][1];
// Lorentz kinematics
//               EV_TEMP(IJ+IEV_TAB,3) = EV_TEMP(IJ,3) + VX_PREF
//               EV_TEMP(IJ+IEV_TAB,4) = EV_TEMP(IJ,4) + VY_PREF
//               EV_TEMP(IJ+IEV_TAB,5) = EV_TEMP(IJ,5) + VZ_PREF
// Lorentz transformation
               lorentz_boost(VX2_FISSION,VY2_FISSION,VZ2_FISSION,
                EV_TEMP2[IJ][2],EV_TEMP2[IJ][3],EV_TEMP2[IJ][4],
                &VXOUT,&VYOUT,&VZOUT);
               lorentz_boost(vx_eva_sc,vy_eva_sc,vz_eva_sc,
                VXOUT,VYOUT,VZOUT,
                &VX2OUT,&VY2OUT,&VZ2OUT);
               EV_TEMP[IJ+IEV_TAB_FIS][2] = VX2OUT;
               EV_TEMP[IJ+IEV_TAB_FIS][3] = VY2OUT;
               EV_TEMP[IJ+IEV_TAB_FIS][4] = VZ2OUT;
               //
               }
               IEV_TAB_FIS = IEV_TAB_FIS + IEV_TAB_TEMP;
      }
//
// Lorentz kinematics 
//      vx1_fission = vx1_fission + vx1_eva
//      vy1_fission = vy1_fission + vy1_eva
//      vz1_fission = vz1_fission + vz1_eva
//      vx2_fission = vx2_fission + vx2_eva
//      vy2_fission = vy2_fission + vy2_eva
//      vz2_fission = vz2_fission + vz2_eva
// The v_eva_sc contribution is considered in the calling subroutine
// Lorentz transformations
               lorentz_boost(vx1_eva,vy1_eva,vz1_eva,
                VX1_FISSION,VY1_FISSION,VZ1_FISSION,
                &VXOUT,&VYOUT,&VZOUT);
               VX1_FISSION = VXOUT;
               VY1_FISSION = VYOUT;
               VZ1_FISSION = VZOUT;
               lorentz_boost(vx2_eva,vy2_eva,vz2_eva,
                VX2_FISSION,VY2_FISSION,VZ2_FISSION,
                &VXOUT,&VYOUT,&VZOUT);
               VX2_FISSION = VXOUT;
               VY2_FISSION = VYOUT;
               VZ2_FISSION = VZOUT;
//
 (*ZFP1) = idnint(ZFF1);
 (*AFP1) = idnint(AFF1);
 (*SFP1) = NbLam1;
 (*VX1_FISSION_par) = VX1_FISSION;
 (*VY1_FISSION_par) = VY1_FISSION;
 (*VZ1_FISSION_par) = VZ1_FISSION;
 (*VX_EVA_SC_par)=vx_eva_sc;
 (*VY_EVA_SC_par)=vy_eva_sc;
 (*VZ_EVA_SC_par)=vz_eva_sc;
 (*ZFP2) = idnint(ZFF2);
 (*AFP2) = idnint(AFF2);
 (*SFP2) = NbLam2;
 (*VX2_FISSION_par) = VX2_FISSION;
 (*VY2_FISSION_par) = VY2_FISSION;
 (*VZ2_FISSION_par) = VZ2_FISSION;
 (*IEV_TAB_FIS_par) = IEV_TAB_FIS;
 (*NbLam0_par) = NbLam1 + NbLam2;
 if(NbLam0>(NbLam1 + NbLam2))varntp->kfis = 25;
 return;
}
//*************************************************************************
//
void G4Abla::tke_bu(G4double Z,G4double A,G4double ZALL,G4double AAL,G4double *VX,G4double *VY,G4double *VZ){

       G4double V_over_V0,R0,RALL,RHAZ,R,TKE,Ekin,V,VPERP,ALPHA1;

       V_over_V0 = 6.0;
       R0 = 1.16;

       if(Z < 1.0){
        *VX = 0.0;
        *VY = 0.0;
        *VZ = 0.0;
        return;
       }

       RALL = R0 * std::pow(V_over_V0,1.0/3.0) * std::pow(AAL,1.0/3.0);
       RHAZ = G4double(haz(1));
       R = std::pow(RHAZ,1.0/3.0) * RALL;
       TKE = 1.44 * Z * ZALL * R*R * (1.0 - A/AAL)*(1.0 - A/AAL) / std::pow(RALL,3.0);

       Ekin = TKE * (AAL - A) / AAL;
//       print*,'!!!',IDNINT(AAl),IDNINT(A),IDNINT(ZALL),IDNINT(Z)
       V = std::sqrt(Ekin/A) * 1.3887;
       *VZ = (2.0 * G4double(haz(1)) - 1.0) * V;
       VPERP = std::sqrt(V*V - (*VZ)*(*VZ));
       ALPHA1 = G4double(haz(1)) * 2.0 * 3.142;
       *VX = VPERP * std::sin(ALPHA1);
       *VY = VPERP * std::cos(ALPHA1);
 return;
}

G4double G4Abla::haz(G4int k)
{
 // const G4int pSize = 110;
 // static G4ThreadLocal G4double p[pSize];
  static G4ThreadLocal G4int ix = 0;
  static G4ThreadLocal G4double x = 0.0, y = 0.0;
  //  k =< -1 on initialise                                        
  //  k = -1 c'est reproductible                                   
  //  k < -1 || k > -1 ce n'est pas reproductible
/*
  // Zero is invalid random seed. Set proper value from our random seed collection:
  if(ix == 0) {
    //    ix = hazard->ial;
  }
*/
  if (k <= -1) { //then                                             
    if(k == -1) { //then                                            
      ix = 0;
    }
    else {
      x = 0.0;
      y = secnds(G4int(x));
      ix = G4int(y * 100 + 43543000);
      if(mod(ix,2) == 0) {
        ix = ix + 1;
      }
    }}

  return G4AblaRandom::flat();
}

//  Random generator according to the
//  powerfunction y = x**(lambda) in the range from xmin to xmax
//  xmin, xmax and y are integers.
//  lambda must be different from -1 !
G4int G4Abla::IPOWERLIMHAZ(G4double lambda,G4int xmin,G4int xmax){
       G4double y,l_plus,rxmin,rxmax;
         l_plus = lambda + 1.;
         rxmin = G4double(xmin) - 0.5;
         rxmax = G4double(xmax) + 0.5;
//       y=(HAZ(k)*(rxmax**l_plus-rxmin**l_plus)+ rxmin**l_plus)**(1.E0/l_plus)
         y=std::pow(G4AblaRandom::flat()*(std::pow(rxmax,l_plus)-std::pow(rxmin,l_plus))+ std::pow(rxmin,l_plus),1.0/l_plus);
         return nint(y);
}

void G4Abla::AMOMENT(G4double AABRA,G4double APRF, G4int IMULTIFR,G4double *PX,G4double *PY,G4double *PZ){

      G4int ISIGOPT = 0;
      G4double GOLDHA_BU=0.,GOLDHA=0.;
      G4double PI = 3.141592653589793;
      //nu = 1.d0

    //  G4double BETAP = sqrt(1.0 - 1.0/sqrt(1.0+EAP/931.494));
    //  G4double GAMMAP = 1.0 / sqrt(1. - BETAP*BETAP);
    //  G4double FACT_PROJ = (GAMMAP + 1.) / (BETAP * GAMMAP);

     // G4double R = 1.160 * pow(APRF,1.0/3.0);

    //  G4double RNDT = double(haz(1));
    //  G4double CTET = 2.0*RNDT-1.0;
    //  G4double TETA = acos(CTET);
    //  G4double RNDP = double(haz(1));
    //  G4double PHI = RNDP*2.0*PI;
    //  G4double STET = sqrt(1.0-CTET*CTET);
//      RX = R * STET * DCOS(PHI)
//      RY = R * STET * DSIN(PHI)
//      RZ = R * CTET

    //  G4double RZ = 0.0;
    //  G4double RY = R * sin(PHI);
    //  G4double RX = R * cos(PHI);

// In MeV/C
      G4double V0_over_VBU = 1.0 / 6.0;
      G4double SIGMA_0 = 118.50;
      G4double Efermi = 5.0 * SIGMA_0 * SIGMA_0 / (2.0 * 931.4940);

      if(IMULTIFR==1){
       if(ISIGOPT == 0){
// "Fermi model" picture:
// Influence of expansion:
        SIGMA_0 = SIGMA_0 * std::pow(V0_over_VBU,1.0/3.0);
// To take into account the influence of thermal motion of nucleons (see W. Bauer,
// PRC 51 (1995) 803)
//        Efermi = 5.D0 * SIGMA_0 * SIGMA_0 / (2.D0 * 931.49D0)

        GOLDHA_BU = SIGMA_0 * std::sqrt((APRF*(AABRA-APRF))/(AABRA-1.0));
        GOLDHA    = GOLDHA_BU*std::sqrt(1.0 +
                    5.0 * PI*PI / 12.0 * (T_freeze_out / Efermi)*(T_freeze_out / Efermi));
//       PRINT*,'AFTER BU fermi:',IDNINT(AABRA),IDNINT(APRF),GOLDHA,
//     &                          GOLDHA_BU
       }else{
// Thermal equilibrium picture (<=> to Boltzmann distribution in momentum with sigma2=M*T)
// The factor (AABRA-APRF)/AP comes from momentum conservation:
        GOLDHA_BU = std::sqrt(APRF *  T_freeze_out * 931.494 *
                   (AABRA - APRF) / AABRA);
        GOLDHA   = GOLDHA_BU;
//       PRINT*,'AFTER BU therm:',IDNINT(AABRA),IDNINT(APRF),GOLDHA,
//     &                          GOLDHA_BU
       }
      }else{
      GOLDHA = SIGMA_0 * std::sqrt((APRF*(AABRA-APRF))/(AABRA-1.0));
      }

      G4int IS = 0;
      mom123:  
      *PX = G4double(gausshaz(1,0.0,GOLDHA));
      IS = IS +1;
      if(IS>100){
      std::cout << "WARNING: GAUSSHAZ CALLED MORE THAN 100 TIMES WHEN CALCULATING PX IN Rn07.FOR. A VALUE WILL BE FORCED." << std::endl;
      *PX = (AABRA-1.0)*931.4940;
      }
      if(std::abs(*PX)>= AABRA*931.494){
//       PRINT*,'VX > C',PX,IDNINT(APRF)
       goto mom123;
      }
      IS = 0;
      mom456:  
      *PY = G4double(gausshaz(1,0.0,GOLDHA));
      IS = IS +1;
      if(IS>100){
      std::cout << "WARNING: GAUSSHAZ CALLED MORE THAN 100 TIMES WHEN CALCULATING PY IN Rn07.FOR. A VALUE WILL BE FORCED." << std::endl;
      *PY = (AABRA-1.0)*931.4940;
      }
      if(std::abs(*PY)>= AABRA*931.494){
//       PRINT*,'VX > C',PX,IDNINT(APRF)
       goto mom456;
      }
      IS = 0;
      mom789:  
      *PZ = G4double(gausshaz(1,0.0,GOLDHA));
      IS = IS +1;
      if(IS>100){
      std::cout << "WARNING: GAUSSHAZ CALLED MORE THAN 100 TIMES WHEN CALCULATING PZ IN Rn07.FOR. A VALUE WILL BE FORCED." << std::endl;
      *PZ = (AABRA-1.0)*931.4940;
      }
      if(std::abs(*PZ)>= AABRA*931.494){
//       PRINT*,'VX > C',PX,IDNINT(APRF)
       goto mom789;
      }
 return;
}

G4double G4Abla::gausshaz(G4int k, G4double xmoy, G4double sig)
{
  // Gaussian random numbers:

  //   1005       C*** TIRAGE ALEATOIRE DANS UNE GAUSSIENNE DE LARGEUR SIG ET MOYENNE XMOY
  static G4ThreadLocal G4int  iset = 0;
  static G4ThreadLocal G4double v1,v2,r,fac,gset,fgausshaz;

  if(iset == 0) { //then                                              
    do {
      v1 = 2.0*haz(k) - 1.0;
      v2 = 2.0*haz(k) - 1.0;
      r = std::pow(v1,2) + std::pow(v2,2);
    } while(r >= 1);

    fac = std::sqrt(-2.*std::log(r)/r);
    gset = v1*fac;
    fgausshaz = v2*fac*sig+xmoy;
    iset = 1;
  }
  else {
    fgausshaz=gset*sig+xmoy;
    iset=0;
  }
  return fgausshaz;                                                         
}
