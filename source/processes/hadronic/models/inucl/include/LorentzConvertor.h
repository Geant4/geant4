#ifndef LORENTZ_CONVERTOR_H
#define LORENTZ_CONVERTOR_H

#include "vector"

class LorentzConvertor {

public:

LorentzConvertor() : degenerated(false) {};

LorentzConvertor(const vector<double>& bmom, double bmass, 
                                 const vector<double>& tmom, double tmass) {

  setBullet(bmom,bmass);
  setTarget(tmom,tmass);
  degenerated = false;
  
}; 

void setBullet(const vector<double>& bmom, double bmass) {
  bullet_mom = bmom;
  bullet_mass = bmass;
//  cout << " bullet: e " << bmom[0] << " mass " << bmass << endl;
};

void setTarget(const vector<double>& tmom, double tmass) {
  target_mom = tmom;
  target_mass = tmass;
//  cout << " target: e " << tmom[0] << " mass " << tmass << endl;
};

void toTheCenterOfMass();
 
void toTheTargetRestFrame(); 

vector<double> backToTheLab(const vector<double>& mom) const; 

double getKinEnergyInTheTRS() const {
 double pv = bullet_mom[1]*target_mom[1] + bullet_mom[2]*target_mom[2] + 
        bullet_mom[3]*target_mom[3];  
 double ekin_trf = (target_mom[0]*bullet_mom[0] - pv)/target_mass - bullet_mass;
 return ekin_trf; 
};

double getTotalSCMEnergy() const { return ecm_tot; };

double getSCMMomentum() const { return pscm; };

double getTRSMomentum() const { return plab; };
 
vector<double> rotate(const vector<double> mom) const; 

vector<double> rotate(const vector<double> mom1,
                                   const vector<double> mom) const; 

bool reflectionNeeded() const; 

bool trivial() const { return degenerated; }; 

private: 

vector<double> bullet_mom;
double bullet_mass;

vector<double> target_mom;
double target_mass;

vector<double> velocity;
vector<double> scm_momentum;
double ecm_tot;
double pscm;
double plab;

double gamma;
double v2;
double ga;
double gb;
double gbpp;
double gapp;

bool degenerated;

};        

#endif // LORENTZ_CONVERTOR_H 
