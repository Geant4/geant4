//#define DEBUG

#include "Fissioner.h"
#include "InuclNuclei.h"
#include "FissionStore.h"
#include "FissionConfiguration.h"

CollisionOutput Fissioner::collide(InuclParticle* bullet,
                     InuclParticle* target) {

const double one_third = 1./3.;
const double two_thirds = 2./3.;
const int itry_max = 1000;

CollisionOutput output;

if(InuclNuclei* nuclei_target = dynamic_cast<InuclNuclei*>(target)) {

#ifdef DEBUG
  cout << " Fissioner input " << endl;
  nuclei_target->printParticle();
#endif
 
  double A = nuclei_target->getA();
  double Z = nuclei_target->getZ();
  double EEXS = nuclei_target->getExitationEnergy();
  double mass_in = nuclei_target->getMass();
  double e_in = mass_in + 0.001*EEXS;
  double PARA = 0.055*pow(A,two_thirds)*(pow(A-Z,one_third) + pow(Z,one_third));
  double TEM = sqrt(EEXS/PARA);
  double TETA = 0.494*pow(A,one_third)*TEM;
  TETA = TETA/sinh(TETA);
  if(A < 246.) PARA += (nucleiLevelDensity(A) - PARA)*TETA;

  double A1 = int(A/2. + 1.1);
  double Z1;
  double A2 = A - A1;
  double ALMA = -1000.;
  double DM1 = bindingEnergy(A,Z);
  double EVV = EEXS - DM1;
  double DM2 = bindingEnergyAsymptotic(A,Z);
  double DTEM = (A < 220. ? 0.5 : 1.15);
  TEM += DTEM;
  
  vector<double> AL1(2,-0.15);
  vector<double> BET1(2,0.05);
  FissionStore fissionStore;
  double R12 = pow(A1,one_third) + pow(A2,one_third); 
  
  for(int i = 0; i < 50 && A1 > 30.; i++) {
    A1 -= 1.;
    A2 = A - A1;
    double X3 = 1./pow(A1,one_third);
    double X4 = 1./pow(A2,one_third);
    Z1 = int(getZopt(A1,A2,Z,X3,X4,R12)) - 1.;
    vector<double> EDEF1(2);
    double Z2 = Z - Z1;
    double VPOT, VCOUL;
    potentialMinimization(VPOT,EDEF1,VCOUL,A1,A2,Z1,Z2,AL1,BET1,R12);

    double DM3 = bindingEnergy(A1,Z1);
    double DM4 = bindingEnergyAsymptotic(A1,Z1);
    double DM5 = bindingEnergy(A2,Z2);
    double DM6 = bindingEnergyAsymptotic(A2,Z2);
    double DMT1 = DM4 + DM6 - DM2;
    double DMT = DM3 + DM5 - DM1;
    double EZL = EEXS + DMT - VPOT;
    
    if(EZL > 0.) { // generate fluctuations
//  faster, using randomGauss

      double C1 = sqrt(getC2(A1,A2,X3,X4,R12)/TEM);
      double DZ = randomGauss(C1);
      DZ = DZ > 0. ? int(DZ + 0.5) : -int(fabs(DZ - 0.5));
      Z1 += DZ;
      Z2 -= DZ;
     
      double DEfin = randomGauss(TEM);	
  
      double EZ = (DMT1+(DMT-DMT1)*TETA-VPOT+DEfin)/TEM;
      if(EZ >= ALMA) ALMA = EZ;
      double EK = VCOUL + DEfin + 0.5*TEM;
      double EV = EVV + bindingEnergy(A1,Z1) + bindingEnergy(A2,Z2) - EK;
       
      if(EV > 0.) fissionStore.addConfig(A1,Z1,EZ,EK,EV);
    };
  };

  int store_size = fissionStore.size();
  if(store_size > 0) {
    FissionConfiguration config = 
               fissionStore.generateConfiguration(ALMA,inuclRndm());
    A1 = config.afirst;
    A2 = A - A1;
    Z1 = config.zfirst;
    double Z2 = Z - Z1;
    InuclNuclei nuclei1(A1,Z1);
    InuclNuclei nuclei2(A2,Z2);        
    double mass1 = nuclei1.getMass();
    double mass2 = nuclei2.getMass();
    double EK = config.ekin;
    double pmod = sqrt(0.001*EK*mass1*mass2/mass_in);
    pair<double,double> COS_SIN = randomCOS_SIN();
    double Fi = randomPHI();
    double P1 = pmod*COS_SIN.second;
    vector<double> mom1(4);
    vector<double> mom2(4);
    mom1[1] = P1*cos(Fi);
    mom1[2] = P1*sin(Fi);
    mom1[3] = pmod*COS_SIN.first;
    for(int i = 1; i < 4; i++) mom2[i] = -mom1[i];
    double e_out = sqrt(pmod*pmod + mass1*mass1) + 
              sqrt(pmod*pmod + mass2*mass2);
    double EV = 1000.*(e_in - e_out)/A;
    if(EV > 0.) {
      double EEXS1 = EV*A1;
      double EEXS2 = EV*A2;
      InuclNuclei nuclei1(mom1,A1,Z1);        
      nuclei1.setExitationEnergy(EEXS1);
      nuclei1.setEnergy();
      output.addTargetFragment(nuclei1);
      InuclNuclei nuclei2(mom2,A2,Z2);        
      nuclei2.setExitationEnergy(EEXS2);
      nuclei2.setEnergy();
      output.addTargetFragment(nuclei2);
#ifdef DEBUG
      nuclei1.printParticle();
      nuclei2.printParticle();
#endif
    };
  };
}
 else {
  cout << " Fissioner -> target is not nuclei " << endl;    
}; 

  return output;
}

double Fissioner::getC2(double A1, double A2, double X3, double X4, double R12) const {
  double C2 = 124.57*(1./A1+1./A2) + 0.78*(X3+X4) - 176.9 *
     (pow(X3,4) + pow(X4,4)) + 219.36 * (1./(A1*A1) + 1./(A2*A2)) - 1.108/R12;
  return C2;   
}

double Fissioner::getZopt(double A1, double A2, double ZT, 
             double X3, double X4, double R12) const {
  double Zopt = (87.7*(X4-X3) * (1.-1.25*(X4+X3)) +
    ZT * ((124.57/A2+0.78*X4-176.9*pow(X4,4) + 219.36/(A2*A2))-0.554/R12))/
      getC2(A1,A2,X3,X4,R12);
  return Zopt;
}	     

void Fissioner::potentialMinimization(double& VP,vector<double> & ED,
   double& VC, double AF, double AS, double ZF, double ZS,
   vector<double>& AL1, vector<double>& BET1, double& R12) const {

const double huge = 2.e35;
const double one_third = 1./3.;
const double two_thirds = 2./3.;
const int itry_max = 2000;
const double DSOL1 = 1.e-6;
const double DS1 = 0.3;
const double DS2 = 1./DS1/DS1;

double A1[2];
A1[0] = AF;
A1[1] = AS;
double Z1[2];
Z1[0] = ZF;
Z1[1] = ZS;
double D = 1.01844*ZF*ZS;
double D0 = 1.E-3*D;
double R[2];
R12 = 0.;
double C[2];
double F[2];
double Y1;
double Y2;
for(int i = 0; i < 2; i++) {
  R[i] = pow(A1[i],one_third);
  Y1 = R[i]*R[i];
  Y2 = Z1[i]*Z1[i]/R[i];
  C[i] = 6.8*Y1 - 0.142*Y2;
  F[i] = 12.138*Y1 - 0.145*Y2; 
};

double SAL[2];
double SBE[2];
double X[2];
double X1[2];
double X2[2];
double RAL[2];
double RBE[2];
double A[4][4];
double B[4];

int itry = 0;
while(itry < itry_max) {
  itry++;
  double S = 0.;
  for(int i = 0; i < 2; i++) 
    S += R[i]*(1. + AL1[i] + BET1[i] - 0.257*AL1[i]*BET1[i]);

  R12 = 0.;
  Y1 = 0.;
  Y2 = 0.;
  for(int i = 0; i < 2; i++) {
    SAL[i] = R[i]*(1.-0.257*BET1[i]);
    SBE[i] = R[i]*(1.-0.257*AL1[i]);
    X[i] = R[i]/S;
    X1[i] = X[i]*X[i];
    X2[i] = X[i]*X1[i];
    Y1 += AL1[i]*X1[i];
    Y2 += BET1[i]*X2[i];
    R12 += R[i]*(1.- AL1[i]*(1. - 0.6*X[i])+BET1[i]*(1.- 0.429*X1[i]));
  }; 

  double Y3 = -0.6*Y1 + 0.857*Y2;
  double Y4 = (1.2*Y1 - 2.571*Y2)/S;
  double R2 = D0/(R12*R12);
  double R3 = 2.*R2/R12;
 
  for(int i = 0; i < 2; i++) {
    RAL[i] =-R[i]*(1.-0.6*X[i])+SAL[i]*Y3;
    RBE[i] = R[i]*(1.-0.429*X1[i])+SBE[i]*Y3;
  };

  double DX1;
  double DX2;
  
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      double DEL1 = i == j ? 1. : 0.;
      DX1 = 0.;
      DX2 = 0.;
      if(fabs(AL1[i]) >= DS1) {
        double XXX = AL1[i]*AL1[i]*DS2;
	double DEX = XXX > 100. ? huge : exp(XXX);
	DX1 = 2. * (1.+2.*AL1[i]*AL1[i]*DS2) * DEX * DS2;
      };
      if(fabs(BET1[i]) >= DS1) {
        double XXX = BET1[i]*BET1[i]*DS2;
	double DEX = XXX > 100. ? huge : exp(XXX);
	DX2 = 2. * (1.+2.*BET1[i]*BET1[i]*DS2) * DEX * DS2;
      };
      double DEL = 2.e-3*DEL1;
      A[i][j] = R3 * RBE[i] * RBE[j] - R2 * (-0.6 * 
       (X1[i]*SAL[j] + X1[j]*SAL[i]) + SAL[i]*SAL[j] * Y4) + 
        DEL * C[i] + DEL1 * DX1;
      int i1 = i + 2;
      int j1 = j + 2;
      A[i1][j1] = R3 * RBE[i] * RBE[j] - R2 * (0.857 * 
        (X2[i]*SBE[j] + X2[j]*SBE[i]) + SBE[i]*SBE[j] * Y4) +
	DEL * F[i] + DEL1 * DX2;
      A[i][j1] = R3 * RAL[i] * RBE[j] - R2 * (0.857 * 
        (X2[j]*SAL[i] - 0.6 * X1[i]*SBE[j]) + SBE[j]*SAL[i] * Y4 - 
	 0.257 * R[i] * Y3 * DEL1);
      A[j1][i] = A[i][j1];   	 
    };
  };
  for(int i = 0; i < 2; i++) {
    DX1 = 0.;
    DX2 = 0.;
    if(fabs(AL1[i]) >= DS1) DX1 = 2.*AL1[i] * DS2 * exp(AL1[i]*AL1[i]*DS2);
    if(fabs(BET1[i]) >= DS1) DX2 = 2.*BET1[i] * DS2 * exp(BET1[i]*BET1[i]*DS2);
    B[i] = R2*RAL[i] - 2.E-3*C[i]*AL1[i] + DX1;
    B[i+2] = R2*RBE[i] - 2.E-3*F[i]*BET1[i] + DX2;
  };

  double ST = 0.;
  double ST1 = 0.;
  for(int i = 0; i < 4; i++) {
    ST += B[i]*B[i];
    for(int j = 0; j < 4; j++) ST1 += A[i][j] * B[i] * B[j];
  };
  double STEP = ST/ST1;
  double DSOL = 0.;
  for(int i = 0; i < 2; i++) {
    AL1[i] += B[i]*STEP;
    BET1[i] += B[i+2]*STEP;
    DSOL += B[i]*B[i] + B[i+2]*B[i+2]; 
  };
  DSOL = sqrt(DSOL);
  if(DSOL < DSOL1) break;
};
if(itry == itry_max) 
   cout << " maximal number of iterations in potentialMinimization " << endl
        << " A1 " << AF << " Z1 " << ZF << endl; 

for(int i = 0; i < 2; i++) 
  ED[i] = F[i] * BET1[i] * BET1[i] + C[i] * AL1[i] * AL1[i]; 

VC = D/R12;
VP = VC + ED[0] + ED[1];

}
  
