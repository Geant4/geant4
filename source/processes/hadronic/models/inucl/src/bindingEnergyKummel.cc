//#define DEBUG

#include "InuclSpecialFunctions.h"

double InuclSpecialFunctions::bindingEnergyKummel(double A, double Z) {

// calculates the nuclei binding energy using Kummel mass formula

//     data for Kummel mass formula 

const double ASH[6] = { 20., 28., 50., 82., 126., 184. };

const double OMN[3] = { 2.511464E0, 7.640294E-2, 2.618317E-2 };

const double OMP[3] = { 2.511464E0, 7.640294E-2, 6.243892E-3 };

const double TMET[5] = { 9.477193E0, 7.166809E-1, 1.2169E-1, 4.07179E-2,
                              2.927847E-2 }; 
			      
const double TMET1[4] = { -4.775384E-2, 5.880883E-5, -1.98964E-3, 3.427961E-6 }; 

const double BET[5] = { 1.535628E4, 1.502538E4, 1.481831E4, 1.459232E4, 1.40664E4 };

const double SHN[5] = { 8., 22., 32., 44., 58. };

const double TNUN[5] = { 8.728778E2, 8.45789E1, 9.210738E1, 5.98E1, 3.724E1 }; 

const double TNUP[5] = { 8.728778E2, 8.45789E1, 3.74E1, 5.42E1, 3.724E1 };

const double TKSN[5] = { 3.34582E1, 3.858216E0, -1.489218E-1, -9.2E-1, -6.36E-1 };

const double TKSP[5] = { 3.34582E1, 3.858216E0, 7.2E-1, -8.2E-1, -6.36E-1 };

const double RON[5] = { -6.441735E3, -5.422058E3, -3.722645E3, -2.437172E4, -2.645867E4 };

const double ROP[5] = { -6.75601E3, -5.877588E3, -3.216382E4, -3.010688E4, -2.645867E4 }; 

const double SINGM[4] = { -3.578531E0, -2.6911E0, -7.487347E-1, 0. };

const double C[2] = { 6.1E3, 8.31E3 };

const double AKU = 6.04122E5;
const double US = 1.661835E4;
const double UC = 6.218614E2;
const double UT = 1.983151E4;
const double SIPG = -2.067547E0;
const double TAU = 2.231241E0;
const double AL0 = 0.151;
const double USB = 1.95114E4;
const double UCB = 753.3;
const double ALD = 15.4941;
const double ALD1 = 17.9439;
const double C3 = 0.7059;
const double C4 = 1.21129;
const double PKLD = 1.7826;
const double DMU = 3.132902E4;
        
double DM;
double AN = A - Z;

int INS;
double ANE;
double HNE;

for(int i = 1; i < 6; i++) {
  if(AN <= ASH[i]) {
    ANE = AN - ASH[i-1];
    HNE = ASH[i] - AN;
    INS = i - 1;
    break;
  };
};
 
int IPS;
double APR;
double HPR;

for(int i = 1; i < 6; i++) {
  if(Z <= ASH[i]) {
    APR = Z - ASH[i-1];
    HPR = ASH[i] - Z;
    IPS = i -1;
    break;
  };
};
 
double PPHN = ANE*HNE;
double PPHZ = APR*HPR;

//           omega terms
double OMT = 0.;
if(INS <= 2) OMT += OMN[INS]*PPHN;
if(IPS <= 2) OMT += OMP[IPS]*PPHZ;
#ifdef DEBUG
cout << " OMT " << OMT << endl;
#endif
//           theta term 1
double TET = 0.;
#ifdef DEBUG
cout << " PPHN " << PPHN << " PPHZ " << PPHZ << endl;
cout << " INS " << INS << " IPS " << IPS << endl;
#endif
if(PPHN > 0.5 && PPHZ > 0.5) {
  int IT = 0;
  if((INS+1)*(IPS+1) > 0) {
    if(INS - IPS == 1) {
      if(INS < 5) {
        IT = INS;
      }
       else {
        IT = -1; 
      };		 
    }
     else {
      IT = -1;
    }; 
  };
#ifdef DEBUG
  cout << " IT " << IT << endl; 
#endif
  if(IT >= 0) TET = TMET[IT] * PPHN * PPHZ;
};
#ifdef DEBUG
cout << " TET " << TET << endl;
#endif
//           theta term 2
double TET1 = 0.;

if(IPS == 2) {
  TET1 += TMET1[0] * PPHZ * PPHZ;
  if(INS == 3) TET1 += TMET1[1] * PPHN * PPHZ * HNE * HPR;
};
if(INS == 3) {
  double TVSP = PPHN * PPHN * HNE;
  TET1 += TMET1[2] * TVSP;
  TET1 += TMET1[3] * TVSP * PPHN;
};
#ifdef DEBUG
cout << " TET1 " << TET1 << endl;
#endif
//        betta, nu, ksy terms
double TBET = 0.;

if(INS != 0) {
  for(int i = 0; i <= INS - 1; i++) TBET += BET[i]*SHN[i];
};
TBET += BET[INS]*ANE;
	
if(IPS != 0) {
  for(int i = 0; i <= IPS - 1; i++) TBET += BET[i]*SHN[i];
};
TBET += BET[IPS]*APR;
#ifdef DEBUG
cout << " TBET " << TBET << endl;
#endif

double TBET1 = 0.;

if(PPHN > 0. || PPHZ > 0.) 
      TBET1 = 0.5*((TNUN[INS] + TKSN[INS]*ANE) * PPHN + 
	    (TNUP[IPS] + TKSP[IPS]*APR) * PPHZ);
TBET -= TBET1;
#ifdef DEBUG
cout << " TBET1 " << TBET1 << endl;
#endif
//           deformation
double TDEF = 0.;
double X = Z*Z/A;
double X1 = pow(A,0.3333333);
double X2 = X1 * X1;

if(IPS != INS && INS >= 3 && IPS >= 2) {
  double X3 = 2.*USB - UCB*X;
  double DNZ = 0.;
  if(!(INS < 3 || (INS - IPS) != 1)) DNZ = C[INS - 3];
  DNZ = TBET1 - DNZ;
  double X4 = 0.2*X2*AL0*AL0*X3;
  
  if(DNZ > X4) {
    double X5 = USB + X*UCB;
    double X6 = log(DNZ/X4);
    double X7 = sqrt(X6);
    double ALM = AL0*(X7 + 0.143*AL0*X5/X3);
    TDEF = -X4*(X6 + 1.) + DNZ + 0.038*X2*pow(AL0*X7,3.)*X5;
  };
};
#ifdef DEBUG
cout << " TDEF " << TDEF << endl;
#endif
//            pairng		
double TPE = 0.; 
double TPEN = 0.;
double TPEP = 0.;
double DN0 = 0.;
double DZ0 = 0.;
double AV = 2. * int(0.5*AN + 0.1);	
if(AN > AV) DN0 = 1.;
AV = 2. * int(0.5*Z + 0.1);
if(Z > AV) DZ0 = 1.;
if(DN0*DZ0 > 0.) TPE = DMU / A;
	
if(DN0 + DZ0 > 0.) {
  if(DN0 > 0.) {
    TPEN = RON[INS] / X1;
    if(INS >= 1) {
      if(AN >= 82.) {
	TPEN /= X1;
	if(INS == 3) {
          double PN = 0.;
	  double HN = 0.;
	  if(AN > 90.) PN = AN - 90.;
	  if(AN < 116.) HN = 116. - AN;
	  TPEN += TAU*PN*HN;
	};
       TPEN += SINGM[INS-1]*PPHN;
      }; 
    };
  };
  if(DZ0 > 0.) {
    TPEP = ROP[IPS] / X1;
    if(IPS == 1) TPEP += SIPG*PPHZ;
    if(IPS > 1) TPEP = TPEP / X1;
  };  
};
TPE += TPEN + TPEP;
#ifdef DEBUG
cout << " TPE " << TPE << " TPEN " << TPEN << " TPEP " << TPEP << endl;
#endif
//            collect everything
DM = (AKU - US*X2 + Z*(Z - 1.)/X1 * (OMT - UC) - UT / A *
	 (AN - Z)*(AN - Z) + TET + TET1 + TBET + TDEF + TPE) * 0.001;	  

#ifdef DEBUG
cout << " kummel " << endl; 
#endif
return DM;

}
