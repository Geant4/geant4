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
// File name:     G4PhotoElectricAngularGeneratorPolarized
//
// Author:       
// 
// Creation date:
//
// Modifications: 
// 20 January 2006       
//
// Class Description: 
//
// Concrete class for PhotoElectric Electron Angular Polarized Distribution Generation 
//
// Class Description: End 
//
// -------------------------------------------------------------------
//
//    

#include "G4PhotoElectricAngularGeneratorPolarized.hh"
#include "G4RotationMatrix.hh"
#include "Randomize.hh"

//    

G4PhotoElectricAngularGeneratorPolarized::G4PhotoElectricAngularGeneratorPolarized(const G4String& name):G4VPhotoElectricAngularDistribution(name)
{
    const G4int arraydim=980;

    //beta minimum
    betarray[0]=0.02;
    //beta step
    betarray[1]=0.001;
    //maximum index array for a and c tables
    betarray[2]=arraydim-1;

    for(G4int level = 0; level < 2; level++){

    char nameChar0[100] = "ftab0.dat";
    char nameChar1[100] = "ftab1.dat";

    G4String filename;
    if(level == 0) filename = nameChar0;
    if(level == 1) filename = nameChar1;

    char* path = getenv("G4LEDATA");
    if (!path)
      {
        G4String excep = "G4EMDataSet - G4LEDATA environment variable not set";
        G4Exception(excep);
      }

    G4String pathString(path);
    G4String dirFile = pathString + "/photoelectric_angular/" + filename;
    FILE *infile;
    infile = fopen(dirFile,"r"); // K-shell
    if (infile == NULL)
      {
      G4String excep = "G4PhotoElectricAngularGeneratorPolarized - data file: " + dirFile + " not found";
      G4Exception(excep);
      }

    G4float a_read,c_read, beta;
    for(G4int i=0 ; i<arraydim ;i++){
	fscanf(infile,"%f\t %e\t %e",&beta,&a_read,&c_read);
	a[i][level]=a_read;    
	c[i][level]=c_read;
    }
    fclose(infile);

    }
}

//    

G4PhotoElectricAngularGeneratorPolarized::~G4PhotoElectricAngularGeneratorPolarized() 
{;}

//

G4ThreeVector G4PhotoElectricAngularGeneratorPolarized::GetPhotoElectronDirection(G4ThreeVector direction, G4double eKineticEnergy,
										  G4ThreeVector polarization, G4int shellId)
{
  G4double gamma   = 1. + eKineticEnergy/electron_mass_c2;
  G4double beta  = std::sqrt(gamma*gamma-1.)/gamma;
  beta = 0.5;

  G4double theta = 0;
  G4double phi = 0;
  G4double a_beta = 0;
  G4double c_beta = 0;

  G4int level = 0;
  if(shellId <  3) level = 0; // K-shell
  if(shellId >= 3) level = 1; // L-shell

  PhotoElectronGetac(level,beta,&a_beta,&c_beta);
  PhotoElectronGenPhiTheta(level,beta,a_beta,c_beta,&phi,&theta);
  G4RotationMatrix rotation = PhotoElectronRotationMatrix(direction,polarization);
  G4ThreeVector final_direction = PhotoElectronGetPlab(rotation,theta,phi);
  return final_direction;
}

//

void G4PhotoElectricAngularGeneratorPolarized::PhotoElectronGenPhiTheta(G4int level,G4double beta, 
									G4double a_beta, G4double c_beta, 
									G4double *pphi, G4double *ptheta){
    G4double xi1,xi2;
    G4double phi = 0;
    G4double theta = 0;
    G4double u = 0;
    G4double xs = 0;
    G4double g = 0;
    G4double max_beta;

    do {

      xi1 = G4UniformRand();
      xi2 = G4UniformRand();
      u = G4UniformRand();
	
	phi=2*M_PI*xi1;

	if(level == 0){

	    theta=std::sqrt(((exp(xi2*std::log(1+c_beta*M_PI*M_PI)))-1)/c_beta);
	    g = G2Function(theta,c_beta);
	    xs = dsigma_k_shellGavrila1959(beta,theta,phi);

	} else {

	    theta = std::sqrt(((exp(xi2*std::log(1+c_beta*M_PI*M_PI)))-1)/c_beta);
	    g = G2Function(theta,c_beta);
	    xs = dsigma_L1_shellGavrila(beta,theta,phi);

	}

	max_beta=u*a_beta*g;

    }while(max_beta > xs);

    *pphi=phi;
    *ptheta=theta;
}

//

G4double G4PhotoElectricAngularGeneratorPolarized::G2Function(G4double theta,G4double c_beta){

  G4double g = 0; 
  g = theta/(1+c_beta*theta*theta);
  return g;

}

//

G4double G4PhotoElectricAngularGeneratorPolarized::dsigma_k_shellGavrila1959(G4double beta, G4double theta, G4double phi){

    //Double differential cross-section
    G4double beta2=beta*beta;
    G4double one_beta2 = 1-beta*beta;
    G4double sqrt_one_beta2 = std::sqrt(one_beta2);
    G4double one_beta2_to_3_2=std::pow(one_beta2,1.5);
    G4double costhe=std::cos(theta);
    G4double sinthe2=std::sin(theta)*std::sin(theta);
    G4double cosphi2=std::cos(phi)*std::cos(phi);
    G4double one_beta_costhe=1-beta*costhe;
    G4double dsigma = 0;
    G4double F = 0;
    G4double G = 0;
    G4double Z = 1;

    F = sinthe2*cosphi2/std::pow(one_beta_costhe,4)-(1 - sqrt_one_beta2)/(2*one_beta2) * 
	(sinthe2 * cosphi2)/std::pow(one_beta_costhe,3) + (1-sqrt_one_beta2)*
	(1-sqrt_one_beta2)/(4*one_beta2_to_3_2) * sinthe2/std::pow(one_beta_costhe,3);

    G = std::sqrt(1 - sqrt_one_beta2)/(std::pow(2,3.5)*beta2*std::pow(one_beta_costhe,2.5)) *
        (4*beta2/sqrt_one_beta2 * sinthe2*cosphi2/one_beta_costhe + 
        4*beta/one_beta2 * costhe * cosphi2
	- 4*(1-sqrt_one_beta2)/one_beta2 *(1+cosphi2)
        - beta2 * (1-sqrt_one_beta2)/one_beta2 * sinthe2/one_beta_costhe
        + 4*beta2*(1-sqrt_one_beta2)/one_beta2_to_3_2
        - 4*beta*(1-sqrt_one_beta2)*(1-sqrt_one_beta2)/one_beta2_to_3_2 * costhe)
        + (1-sqrt_one_beta2)/(4*beta2*one_beta_costhe*one_beta_costhe) *
        (beta/one_beta2 - 2/one_beta2 * costhe * cosphi2 + 
        (1-sqrt_one_beta2)/one_beta2_to_3_2 * costhe
	 - beta * (1-sqrt_one_beta2)/one_beta2_to_3_2);

    dsigma = ( F*(1-M_PI*fine_structure_const*Z/beta) + M_PI*fine_structure_const*Z*G);

    return dsigma;
}

//

G4double G4PhotoElectricAngularGeneratorPolarized::dsigma_L1_shellGavrila(G4double beta, G4double theta, G4double phi){

  //Double differential cross-section
  G4double beta2=beta*beta;
  G4double one_beta2 = 1-beta*beta;
  G4double sqrt_one_beta2 = std::sqrt(one_beta2);
  G4double one_beta2_to_3_2=std::pow(one_beta2,1.5);
  G4double costhe=std::cos(theta);
  G4double sinthe2=std::sin(theta)*std::sin(theta);
  G4double cosphi2=std::cos(phi)*std::cos(phi);
  G4double one_beta_costhe=1-beta*costhe;
	
  G4double dsigma = 0;
  G4double F = 0;
  G4double G = 0;
  G4double Z = 1;

  F = sinthe2*cosphi2/std::pow(one_beta_costhe,4)-(1 - sqrt_one_beta2)/(2*one_beta2)
     *  (sinthe2 * cosphi2)/std::pow(one_beta_costhe,3) + (1-sqrt_one_beta2)*
     (1-sqrt_one_beta2)/(4*one_beta2_to_3_2) * sinthe2/std::pow(one_beta_costhe,3);

  G = std::sqrt(1 - sqrt_one_beta2)/(std::pow(2,3.5)*beta2*std::pow(one_beta_costhe,2.5)) *
     (4*beta2/sqrt_one_beta2 * sinthe2*cosphi2/one_beta_costhe +       
      4*beta/one_beta2 * costhe * cosphi2
     - 4*(1-sqrt_one_beta2)/one_beta2 *(1+cosphi2)
     - beta2 * (1-sqrt_one_beta2)/one_beta2 * sinthe2/one_beta_costhe
     + 4*beta2*(1-sqrt_one_beta2)/one_beta2_to_3_2
     - 4*beta*(1-sqrt_one_beta2)*(1-sqrt_one_beta2)/one_beta2_to_3_2 * costhe)
     + (1-sqrt_one_beta2)/(4*beta2*one_beta_costhe*one_beta_costhe) *
     (beta/one_beta2 - 2/one_beta2 * costhe * cosphi2 + 
     (1-sqrt_one_beta2)/one_beta2_to_3_2 * costhe
      - beta * (1-sqrt_one_beta2)/one_beta2_to_3_2);

  dsigma = ( F*(1-M_PI*fine_structure_const*Z/beta) + M_PI*fine_structure_const*Z*G);

 return dsigma;

}

G4double G4PhotoElectricAngularGeneratorPolarized::getMax(G4double arg1, G4double arg2)
{
    if (arg1 > arg2)
	return arg1;
    else
	return arg2;
}

//

G4RotationMatrix G4PhotoElectricAngularGeneratorPolarized::PhotoElectronRotationMatrix(G4ThreeVector direction, 
										       G4ThreeVector polarization)
{
    G4double mK = direction.mag();
    G4double mS = polarization.mag();

    if(!(polarization.isOrthogonal(direction,1e-6)) || mS == 0){
      G4ThreeVector d0 = direction.unit();
      G4ThreeVector a1 = SetPerpendicularVector(d0); 
      G4ThreeVector a0 = a1.unit(); 
      G4double rand1 = G4UniformRand();
      G4double angle = twopi*rand1; 
      G4ThreeVector b0 = d0.cross(a0); 
      G4ThreeVector c;
      c.setX(std::cos(angle)*(a0.x())+std::sin(angle)*b0.x());
      c.setY(std::cos(angle)*(a0.y())+std::sin(angle)*b0.y());
      c.setZ(std::cos(angle)*(a0.z())+std::sin(angle)*b0.z());
      polarization = c.unit();
      mS = polarization.mag();
    }else
      {
	if ( polarization.howOrthogonal(direction) != 0)
	  {
	    polarization = polarization - polarization.dot(direction)/direction.dot(direction) * direction;
	  }
      }

    direction = direction/mK;
    polarization = polarization/mS;

    G4ThreeVector y = direction.cross(polarization);
    
    G4RotationMatrix R(polarization,y,direction);
    return R;
}

void G4PhotoElectricAngularGeneratorPolarized::PhotoElectronGetac(G4int level, G4double beta,G4double *pa, G4double *pc)
{
    G4int k;
    G4double a_beta,c_beta;
    G4double bmin,bstep;
    G4int indexmax;    
    if(level > 1) level = 1; // protection since we have only K and L1 polarized double differential cross-sections
    
    bmin=betarray[0];
    bstep=betarray[1];
    indexmax=(G4int)betarray[2];

    k=(G4int)((beta-bmin+1e-9)/bstep);    
    
    if(k<0)
	k=0;
    if(k>indexmax)
	k=indexmax; 
    
    
    if(k==0)
	a_beta=getMax(a[k][level],a[k+1][level]);
    else if(k==indexmax)
	a_beta=getMax(a[k-1][level],a[k][level]);
    else{
	a_beta=getMax(a[k-1][level],a[k][level]);
	a_beta=getMax(a_beta,a[k+1][level]);
	}   
    
    if(k==0)
	c_beta=getMax(c[k][level],c[k+1][level]);
    else if(k==indexmax)
	c_beta=getMax(c[k-1][level],c[k][level]);
    else{
	c_beta=getMax(c[k-1][level],c[k][level]);
	c_beta=getMax(c_beta,c[k+1][level]);
    }

    *pa=a_beta;
    *pc=c_beta;

}


//
G4ThreeVector G4PhotoElectricAngularGeneratorPolarized::PhotoElectronGetPlab(G4RotationMatrix rotation,G4double theta, G4double phi){

//computes the photoelectron momentum unitary vector 
  G4double px = std::cos(phi)*std::sin(theta);
  G4double py = std::sin(phi)*std::sin(theta);
  G4double pz = std::cos(theta);

  G4ThreeVector sampling_direction(px,py,pz);

  G4ThreeVector outgoing_direction = rotation*sampling_direction;
  return outgoing_direction;
}

//

void G4PhotoElectricAngularGeneratorPolarized::PrintGeneratorInformation() const
{
  G4cout << "\n" << G4endl;
  G4cout << "Polarized Photoelectric Angular Generator" << G4endl;
  G4cout << "(see Physics Reference Manual) \n" << G4endl;
} 

G4ThreeVector G4PhotoElectricAngularGeneratorPolarized::SetPerpendicularVector(G4ThreeVector& a)
{
  G4double dx = a.x();
  G4double dy = a.y();
  G4double dz = a.z();
  G4double x = dx < 0.0 ? -dx : dx;
  G4double y = dy < 0.0 ? -dy : dy;
  G4double z = dz < 0.0 ? -dz : dz;
  if (x < y) {
    return x < z ? G4ThreeVector(-dy,dx,0) : G4ThreeVector(0,-dz,dy);
  }else{
    return y < z ? G4ThreeVector(dz,0,-dx) : G4ThreeVector(-dy,dx,0);
  }
}
