#ifndef __NUCLMODEL_DEFINED__
#define __NUCLMODEL_DEFINED__

struct VECTOR;

double NormalizeDensity(double A,unsigned nSlices=50);
double GetDensity(double Norm,double Rad,double A);
double GetPotential(double Rad,double A,double Z);
enum ProbabilitiesTypes {EQUAL_TYPE,DIFF_TYPE};
enum ParticlesEnum {PROTON,NEUTRON};
double GetElasticData(double Ecm,unsigned nType,double& Total);
double GetAngle(double E1,double M1,double E2,double M2);
void ModifyMomentum(double costheta,double phi,VECTOR& Mom);
double GetXSection(double Elab,unsigned type);
#endif
