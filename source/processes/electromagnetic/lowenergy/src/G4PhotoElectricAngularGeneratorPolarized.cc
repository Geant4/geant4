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
// Author: A. C. Farinha, L. Peralta, P. Rodrigues and A. Trindade
// 
// Creation date:
//
// Modifications: 
// 10 January 2006       
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
  const G4int arrayDim=980;

  //beta minimum
  betarray[0]=0.02;
  //beta step
  betarray[1]=0.001;
  //maximum index array for a and c tables
  betarray[2]=arrayDim-1;

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
    if (infile == 0)
      {
	G4String excep = "G4PhotoElectricAngularGeneratorPolarized - data file: " + dirFile + " not found";
	G4Exception(excep);
      }

    G4float aRead,cRead, beta;
    for(G4int i=0 ; i<arrayDim ;i++){
      fscanf(infile,"%f\t %e\t %e",&beta,&aRead,&cRead);
      a[i][level]=aRead;    
      c[i][level]=cRead;
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
  G4double aBeta = 0;
  G4double cBeta = 0;

  G4int level = 0;
  if(shellId <  3) level = 0; // K-shell
  if(shellId >= 3) level = 1; // L-shell

  PhotoElectronGetac(level,beta,&aBeta,&cBeta);
  PhotoElectronGenPhiTheta(level,beta,aBeta,cBeta,&phi,&theta);
  G4RotationMatrix rotation = PhotoElectronRotationMatrix(direction,polarization);
  G4ThreeVector final_direction = PhotoElectronGetPlab(rotation,theta,phi);
  return final_direction;
}

//

void G4PhotoElectricAngularGeneratorPolarized::PhotoElectronGenPhiTheta(G4int level,G4double beta, 
									G4double aBeta, G4double cBeta, 
									G4double *pphi, G4double *ptheta){
  G4double xi1,xi2;
  G4double phi = 0;
  G4double theta = 0;
  G4double u = 0;
  G4double xs = 0;
  G4double g = 0;
  G4double maxBeta;

  do {

    xi1 = G4UniformRand();
    xi2 = G4UniformRand();
    u = G4UniformRand();
	
    phi=2*pi*xi1;

    if(level == 0){

      theta=std::sqrt(((exp(xi2*std::log(1+cBeta*pi*pi)))-1)/cBeta);
      g = G2Function(theta,cBeta);
      xs = DSigmaKshellGavrila1959(beta,theta,phi);

    } else {

      theta = std::sqrt(((exp(xi2*std::log(1+cBeta*pi*pi)))-1)/cBeta);
      g = G2Function(theta,cBeta);
      xs = DSigmaL1shellGavrila(beta,theta,phi);

    }

    maxBeta=u*aBeta*g;

  }while(maxBeta > xs);

  *pphi=phi;
  *ptheta=theta;
}

//

G4double G4PhotoElectricAngularGeneratorPolarized::G2Function(G4double theta,G4double cBeta){

  G4double g = 0; 
  g = theta/(1+cBeta*theta*theta);
  return g;

}

//

G4double G4PhotoElectricAngularGeneratorPolarized::DSigmaKshellGavrila1959(G4double beta, G4double theta, G4double phi){

  //Double differential cross-section
  G4double beta2=beta*beta;
  G4double oneBeta2 = 1-beta*beta;
  G4double sqrt_oneBeta2 = std::sqrt(oneBeta2);
  G4double oneBeta2_to_3_2=std::pow(oneBeta2,1.5);
  G4double costhe=std::cos(theta);
  G4double sinthe2=std::sin(theta)*std::sin(theta);
  G4double cosphi2=std::cos(phi)*std::cos(phi);
  G4double oneBeta_costhe=1-beta*costhe;
  G4double dsigma = 0;
  G4double F = 0;
  G4double G = 0;
  G4double Z = 1;

  F = sinthe2*cosphi2/std::pow(oneBeta_costhe,4)-(1 - sqrt_oneBeta2)/(2*oneBeta2) * 
    (sinthe2 * cosphi2)/std::pow(oneBeta_costhe,3) + (1-sqrt_oneBeta2)*
    (1-sqrt_oneBeta2)/(4*oneBeta2_to_3_2) * sinthe2/std::pow(oneBeta_costhe,3);

  G = std::sqrt(1 - sqrt_oneBeta2)/(std::pow(2,3.5)*beta2*std::pow(oneBeta_costhe,2.5)) *
    (4*beta2/sqrt_oneBeta2 * sinthe2*cosphi2/oneBeta_costhe + 
     4*beta/oneBeta2 * costhe * cosphi2
     - 4*(1-sqrt_oneBeta2)/oneBeta2 *(1+cosphi2)
     - beta2 * (1-sqrt_oneBeta2)/oneBeta2 * sinthe2/oneBeta_costhe
     + 4*beta2*(1-sqrt_oneBeta2)/oneBeta2_to_3_2
     - 4*beta*(1-sqrt_oneBeta2)*(1-sqrt_oneBeta2)/oneBeta2_to_3_2 * costhe)
    + (1-sqrt_oneBeta2)/(4*beta2*oneBeta_costhe*oneBeta_costhe) *
    (beta/oneBeta2 - 2/oneBeta2 * costhe * cosphi2 + 
     (1-sqrt_oneBeta2)/oneBeta2_to_3_2 * costhe
     - beta * (1-sqrt_oneBeta2)/oneBeta2_to_3_2);

  dsigma = ( F*(1-pi*fine_structure_const*Z/beta) + pi*fine_structure_const*Z*G);

  return dsigma;
}

//

G4double G4PhotoElectricAngularGeneratorPolarized::DSigmaL1shellGavrila(G4double beta, G4double theta, G4double phi){

  //Double differential cross-section
  G4double beta2=beta*beta;
  G4double oneBeta2 = 1-beta*beta;
  G4double sqrt_oneBeta2 = std::sqrt(oneBeta2);
  G4double oneBeta2_to_3_2=std::pow(oneBeta2,1.5);
  G4double costhe=std::cos(theta);
  G4double sinthe2=std::sin(theta)*std::sin(theta);
  G4double cosphi2=std::cos(phi)*std::cos(phi);
  G4double oneBeta_costhe=1-beta*costhe;
	
  G4double dsigma = 0;
  G4double F = 0;
  G4double G = 0;
  G4double Z = 1;

  F = sinthe2*cosphi2/std::pow(oneBeta_costhe,4)-(1 - sqrt_oneBeta2)/(2*oneBeta2)
    *  (sinthe2 * cosphi2)/std::pow(oneBeta_costhe,3) + (1-sqrt_oneBeta2)*
    (1-sqrt_oneBeta2)/(4*oneBeta2_to_3_2) * sinthe2/std::pow(oneBeta_costhe,3);

  G = std::sqrt(1 - sqrt_oneBeta2)/(std::pow(2,3.5)*beta2*std::pow(oneBeta_costhe,2.5)) *
    (4*beta2/sqrt_oneBeta2 * sinthe2*cosphi2/oneBeta_costhe +       
     4*beta/oneBeta2 * costhe * cosphi2
     - 4*(1-sqrt_oneBeta2)/oneBeta2 *(1+cosphi2)
     - beta2 * (1-sqrt_oneBeta2)/oneBeta2 * sinthe2/oneBeta_costhe
     + 4*beta2*(1-sqrt_oneBeta2)/oneBeta2_to_3_2
     - 4*beta*(1-sqrt_oneBeta2)*(1-sqrt_oneBeta2)/oneBeta2_to_3_2 * costhe)
    + (1-sqrt_oneBeta2)/(4*beta2*oneBeta_costhe*oneBeta_costhe) *
    (beta/oneBeta2 - 2/oneBeta2 * costhe * cosphi2 + 
     (1-sqrt_oneBeta2)/oneBeta2_to_3_2 * costhe
     - beta * (1-sqrt_oneBeta2)/oneBeta2_to_3_2);

  dsigma = ( F*(1-pi*fine_structure_const*Z/beta) + pi*fine_structure_const*Z*G);

  return dsigma;

}

G4double G4PhotoElectricAngularGeneratorPolarized::GetMax(G4double arg1, G4double arg2)
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
  G4double aBeta,cBeta;
  G4double bMin,bstep;
  G4int indexMax;    
  if(level > 1) level = 1; // protection since we have only K and L1 polarized double differential cross-sections
    
  bMin=betarray[0];
  bstep=betarray[1];
  indexMax=(G4int)betarray[2];

  k=(G4int)((beta-bMin+1e-9)/bstep);    
    
  if(k<0)
    k=0;
  if(k>indexMax)
    k=indexMax; 
    
    
  if(k==0)
    aBeta=GetMax(a[k][level],a[k+1][level]);
  else if(k==indexMax)
    aBeta=GetMax(a[k-1][level],a[k][level]);
  else{
    aBeta=GetMax(a[k-1][level],a[k][level]);
    aBeta=GetMax(aBeta,a[k+1][level]);
  }   
    
  if(k==0)
    cBeta=GetMax(c[k][level],c[k+1][level]);
  else if(k==indexMax)
    cBeta=GetMax(c[k-1][level],c[k][level]);
  else{
    cBeta=GetMax(c[k-1][level],c[k][level]);
    cBeta=GetMax(cBeta,c[k+1][level]);
  }

  *pa=aBeta;
  *pc=cBeta;

}


//
G4ThreeVector G4PhotoElectricAngularGeneratorPolarized::PhotoElectronGetPlab(G4RotationMatrix rotation,G4double theta, G4double phi){

  //computes the photoelectron momentum unitary vector 
  G4double px = std::cos(phi)*std::sin(theta);
  G4double py = std::sin(phi)*std::sin(theta);
  G4double pz = std::cos(theta);

  G4ThreeVector samplingDirection(px,py,pz);

  G4ThreeVector outgoingDirection = rotation*samplingDirection;
  return outgoingDirection;
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
