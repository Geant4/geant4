#include "InuclSpecialFunctions.h"

pair<vector<double>,vector<double> > InuclSpecialFunctions::paraMaker(double Z) {
// calculates the coefficients for the phenomenological formulas for
// coulumb barier, c.s. etc needed for evaporators

const double Z1[5] = { 10., 20., 30., 50., 70. };
const double AP[5] = { 0.42, 0.58, 0.68, 0.77, 0.8 };
const double CP[5] = { 0.5, 0.28, 0.2, 0.15, 0.1 };
const double AA[5] = { 0.68, 0.82, 0.91, 0.97, 0.98 };
const double CA[5] = { 0.1, 0.1, 0.1, 0.08, 0.06 };

vector<double> AK(6);
vector<double> CPA(6);

AK[0] = CPA[0] = 0.;

double AK2, CP2, AK6, CP6;

if(Z < 10.) {
        AK2=0.42;
        CP2=0.5;
        AK6=0.68;
        CP6=0.1;
}
 else if(Z > 70.) {
          AK6=0.98;  // normal
          CP6=0.06;
          AK2=0.8;
          CP2=0.1;
//          AK6=1.1; // modified
//          CP6=0.;
}
 else {
  for(int i = 1; i < 5; i++) {
    if(Z <= Z1[i]) {
      double Z2 = 1./(Z1[i] - Z1[i-1]);
      AK2 = ((AP[i] - AP[i-1])*Z + AP[i-1]*Z1[i] - AP[i]*Z1[i-1])*Z2;
      CP2 = ((CP[i] - CP[i-1])*Z + CP[i-1]*Z1[i] - CP[i]*Z1[i-1])*Z2;
      AK6 = ((AA[i] - AA[i-1])*Z + AA[i-1]*Z1[i] - AA[i]*Z1[i-1])*Z2;
      CP6 = ((CA[i] - CA[i-1])*Z + CA[i-1]*Z1[i] - CA[i]*Z1[i-1])*Z2;
      break;
    };
  };
};

AK[1] = AK2;
AK[5] = AK6;
CPA[1] = CP2;
CPA[5] = CP6;
AK[2] = AK2 + 0.06;
CPA[2] = CP2 * 0.5;
AK[3] = AK2 + 0.12;
CPA[3] = CP2/3.;  
AK[4] = AK6 - 0.06;
CPA[4] = 4.*CP6/3.;

return pair<vector<double>,vector<double> >(AK,CPA);

}

pair<double,double> InuclSpecialFunctions::paraMakerTruncated(double Z) {
// truncated version of the previous one

const double Z1[5] = { 10., 20., 30., 50., 70. };
const double AP[5] = { 0.42, 0.58, 0.68, 0.77, 0.8 };
const double CP[5] = { 0.5, 0.28, 0.2, 0.15, 0.1 };

double AK2, CP2;

if(Z < 10.) {
        AK2=0.42;
        CP2=0.5;
}
 else if(Z > 70.) {
          AK2=0.8;
          CP2=0.1;
}
 else {
  for(int i = 1; i < 5; i++) {
    if(Z < Z1[i]) {
      double Z2 = 1./(Z1[i] - Z1[i-1]);
      AK2 = ((AP[i] - AP[i-1])*Z + AP[i-1]*Z1[i] - AP[i]*Z1[i-1])*Z2;
      CP2 = ((CP[i] - CP[i-1])*Z + CP[i-1]*Z1[i] - CP[i]*Z1[i-1])*Z2;
      break;
    };
  };
};


return pair<double,double>(AK2,CP2);

}
