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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4Generator2BN
//
// Author:        Andreia Trindade (andreia@lip.pt)
//                Pedro Rodrigues  (psilva@lip.pt)
//                Luis Peralta     (luis@lip.pt)
//
// Creation date: 21 June 2003
//
// Modifications: 
// 21 Jun 2003                                 First implementation acording with new design
// 05 Nov 2003  MGP                            Fixed compilation warning
//
// Class Description: 
//
// Concrete base class for Bremsstrahlung Angular Distribution Generation - 2BN Distribution
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//    

#include "G4Generator2BN.hh"
#include "Randomize.hh"
//    

G4Generator2BN::G4Generator2BN(const G4String& name):G4VBremAngularDistribution(name)
{
  b = 1.2;
  index_min = -300;
  index_max = 20;

  // Set parameters minimum limits Ekelectron = 250 eV and kphoton = 100 eV
  kmin = 100*eV;
  Ekmin = 250*eV;
  kcut = 100*eV;

  // Increment Theta value for surface interpolation
  dtheta = 0.1*rad;

  // Majorant Surface Values
  Atab = new G4double[320];
  ctab = new G4double[320];

  // Construct Majorant Surface to 2BN cross-section
  ConstructMajorantSurface();
}

//    

G4Generator2BN::~G4Generator2BN() 
{

  delete[] Atab;
  delete[] ctab;

}

//

G4double G4Generator2BN::PolarAngle(const G4double initial_energy,
				    const G4double final_energy,
				    const G4int) // Z
{

  G4double theta = 0;

  G4double k = initial_energy - final_energy;

  theta = Generate2BN(initial_energy, k);

  return theta;
}
//

G4double G4Generator2BN::CalculateFkt(G4double k, G4double theta, G4double A, G4double c) const
{
  G4double Fkt_value = 0;
  Fkt_value = A*pow(k,-b)*theta/(1+c*theta*theta);
  return Fkt_value;
}

G4double G4Generator2BN::Calculatedsdkdt(G4double kout, G4double theta, G4double Eel) const
{

  G4double dsdkdt_value = 0.;
  G4double Z = 1;
  // classic radius (in cm)
  G4double r0 = 2.82E-13;
  //squared classic radius (in barn)
  G4double r02 = pow(r0,2)*1.0E+24;

  // Photon energy cannot be greater than electron kinetic energy
  if(kout > (Eel-electron_mass_c2)){
    dsdkdt_value = 0;
    return dsdkdt_value;
  }

     G4double E0 = Eel/electron_mass_c2;
     G4double k =  kout/electron_mass_c2;
     G4double E =  E0 - k;

     if(E <= 1*MeV ){                                 // Kinematic limit at 1 MeV !!
       dsdkdt_value = 0;
       return dsdkdt_value;
     }

     G4double p0 = sqrt(E0*E0-1);
     G4double p = sqrt(E*E-1);
     G4double L = log((E*E0-1+p*p0)/(E*E0-1-p*p0));
     G4double delta0 = E0 - p0*cos(theta);
     G4double epsilon = log((E+p)/(E-p));
     G4double Q = sqrt(pow(p0,2)+pow(k,2)-2*k*p0*cos(theta));
     G4double epsilonQ = log((Q+p)/(Q-p));

     dsdkdt_value = pow(Z,2) * (r02/(8*M_PI*137)) * (1/k) * (p/p0) *
       ( (8 * (pow(sin (theta),2)*(2*pow(E0,2)+1))/(pow(p0,2)*pow(delta0,4))) -
	 ((2*(5*pow(E0,2)+2*E*E0+3))/(pow(p0,2) * pow(delta0,2))) -
	 ((2*(pow(p0,2)-pow(k,2)))/((pow(Q,2)*pow(delta0,2)))) +
	 ((4*E)/(pow(p0,2)*delta0)) +
	 (L/(p*p0))*(
		 ((4*E0*pow(sin(theta),2)*(3*k-pow(p0,2)*E))/(pow(p0,2)*pow(delta0,4))) + 
		 ((4*pow(E0,2)*(pow(E0,2)+pow(E,2)))/(pow(p0,2)*pow(delta0,2))) +
		 ((2-2*(7*pow(E0,2)-3*E*E0+pow(E,2)))/(pow(p0,2)*pow(delta0,2))) +
		 (2*k*(pow(E0,2)+E*E0-1))/((pow(p0,2)*delta0)) 
		 ) - 
         ((4*epsilon)/(p*delta0)) +
         ((epsilonQ)/(p*Q))*
         (4/pow(delta0,2)-(6*k/delta0)-(2*k*(pow(p0,2)-pow(k,2)))/(pow(Q,2)*delta0))
	 );

     dsdkdt_value = dsdkdt_value*sin(theta);
     return dsdkdt_value;
}

void G4Generator2BN::ConstructMajorantSurface()
{

  G4double Eel;
  G4int vmax;
  G4double Ek;
  G4double k, theta, k0, theta0, thetamax;
  G4double ds, df, dsmax, dk, dt;
  G4double ratmin;
  G4double ratio = 0;
  G4double Vds, Vdf;
  G4double A, c;

  G4cout << "**** Constructing Majorant Surface for 2BN Distribution ****" << G4endl;

  if(kcut > kmin) kmin = kcut;

  G4int i = 0;
  // Loop on energy index
  for(G4int index = index_min; index < index_max; index++){

  G4double fraction = index/100.;
  Ek = pow(10.,fraction);
  Eel = Ek + electron_mass_c2;

  // find x-section maximum at k=kmin
  dsmax = 0.;
  thetamax = 0.;

  for(theta = 0.; theta < M_PI; theta = theta + dtheta){

    ds = Calculatedsdkdt(kmin, theta, Eel);

    if(ds > dsmax){
      dsmax = ds;
      thetamax = theta;
    }
  }

  // Compute surface paremeters at kmin
  if(Ek < kmin || thetamax == 0){
    c = 0;
    A = 0;
  }else{
    c = 1/pow(thetamax,2);
    A = 2*sqrt(c)*dsmax/(pow(kmin,-b));
  }

  // look for correction factor to normalization at kmin 
  ratmin = 1.;

  // Volume under surfaces
  Vds = 0.;
  Vdf = 0.;
  k0 = 0.;
  theta0 = 0.;

  vmax = G4int(100.*log10(Ek/kmin));

  for(G4int v = 0; v < vmax; v++){
    G4double fraction = (v/100.);
    k = pow(10.,fraction)*kmin;

    for(theta = 0.; theta < M_PI; theta = theta + dtheta){
      dk = k - k0;
      dt = theta - theta0;
      ds = Calculatedsdkdt(k,theta, Eel);
      Vds = Vds + ds*dk*dt;
      df = CalculateFkt(k,theta, A, c);
      Vdf = Vdf + df*dk*dt;

      if(ds != 0.){
	if(df != 0.) ratio = df/ds;
      }

      if(ratio < ratmin && ratio != 0.){
	ratmin = ratio;
      }
    }
  }


  // sampling efficiency
  Atab[i] = A/ratmin * 1.04;
  ctab[i] = c;

//  G4cout << Ek << " " << i << " " << index << " " << Atab[i] << " " << ctab[i] << " " << G4endl;
  i++;
  }

}

G4double G4Generator2BN::Generate2BN(G4double Ek, G4double k) const 
{

  G4double Eel;
  G4double kmin2; 
  G4double kmax, t;
  G4double cte2;
  G4double y, u;
  G4double fk, ft;
  G4double ds;
  G4double A2;
  G4double A, c;

  G4int trials, index;

  // find table index
  index = G4int(log10(Ek)*100) - index_min;
  Eel = Ek + electron_mass_c2;

  kmax = Ek;
  kmin2 = kcut;

  c = ctab[index];
  A = Atab[index];
  if(index < index_max){
    A2 = Atab[index+1];
    if(A2 > A) A = A2;
  }

  do{
  // generate k accordimg to pow(k,-b)
  trials++;

  // normalization constant 
//  cte1 = (1-b)/(pow(kmax,(1-b))-pow(kmin2,(1-b)));
//  y = G4UniformRand();
//  k = pow(((1-b)/cte1*y+pow(kmin2,(1-b))),(1/(1-b)));

  // generate theta accordimg to theta/(1+c*pow(theta,2)
  // Normalization constant
  cte2 = 2*c/log(1+c*pow(M_PI,2));

  y = G4UniformRand();
  t = sqrt((exp(2*c*y/cte2)-1)/c);
  u = G4UniformRand();

  // point acceptance algorithm
  fk = pow(k,-b);
  ft = t/(1+c*t*t);
  ds = Calculatedsdkdt(k,t,Eel);

  // violation point
  if(ds > (A*fk*ft)) G4cout << "WARNING IN G4Generator2BN !!!" << Ek << " " <<  (ds-A*fk*ft)/ds << G4endl;

  }while(u*A*fk*ft > ds);

  return t;

}

void G4Generator2BN::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Bremsstrahlung Angular Generator is 2BN Generator from 2BN Koch & Motz distribution (Rev Mod Phys 31(4), 920 (1959))" << G4endl;
  G4cout << "\n" << G4endl;
} 

