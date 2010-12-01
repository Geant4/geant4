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
// $Id: G4PenelopeBremsstrahlungAngular.cc,v 1.10 2010-12-01 15:20:20 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// --------------------------------------------------------------
//
// File name:     G4PenelopeBremsstrahlungAngular
//
// Author:        Luciano Pandola
// 
// Creation date: February 2003
//
// History:
// -----------
// 04 Feb 2003  L. Pandola       1st implementation
// 19 Mar 2003  L. Pandola       Bugs fixed
// 07 Nov 2003  L. Pandola       Added GetAtomicNumber method for testing 
//                               purposes
//----------------------------------------------------------------

#include "G4PenelopeBremsstrahlungAngular.hh"
#include "G4PenelopeInterpolator.hh"
#include "Randomize.hh"
#include "globals.hh"

G4PenelopeBremsstrahlungAngular::G4PenelopeBremsstrahlungAngular (G4int Zed)
  : Zmat(Zed)
{
  InterpolationTableForZ();
  InterpolationForK();
}


G4PenelopeBremsstrahlungAngular::~G4PenelopeBremsstrahlungAngular()
{
}

G4int G4PenelopeBremsstrahlungAngular::GetAtomicNumber()
{
  return Zmat;
}

void G4PenelopeBremsstrahlungAngular::InterpolationTableForZ()
{
  G4double pZ[NumberofZPoints] = {2.0,8.0,13.0,47.0,79.0,92.0};
  G4double pX[NumberofZPoints],pY[NumberofZPoints];
  G4double QQ1[NumberofZPoints][NumberofEPoints][NumberofKPoints];
  G4double QQ2[NumberofZPoints][NumberofEPoints][NumberofKPoints];
 
  //Read information from DataBase file
  char* path = getenv("G4LEDATA");
  if (!path)
    {
      G4String excep = "G4PenelopeBremsstrahlungAngular - G4LEDATA environment variable not set!";
      G4Exception(excep);
      return;
    }
  G4String pathString(path);
  G4String pathFile = pathString + "/penelope/br-ang-pen.dat";
  std::ifstream file(pathFile);
  std::filebuf* lsdp = file.rdbuf();
  
  if (!(lsdp->is_open()))
    {
      G4String excep = "G4PenelopeBremsstrahlungAngular - data file " + pathFile + " not found!";
      G4Exception(excep);
    }
  G4int i=0,j=0,k=0; // i=index for Z, j=index for E, k=index for K 
  G4double a1,a2;
  while(i != -1) {
    file >> i >> j >> k >> a1 >> a2; 
    if (i > -1 && j > -1 && k >- 1)
      {
	QQ1[i][j][k]=a1;
	QQ2[i][j][k]=a2;
      }
  } 
  file.close();
  

  //Interpolation in Z
  for (i=0;i<NumberofEPoints;i++){
    for (j=0;j<NumberofKPoints;j++){
      for (k=0;k<NumberofZPoints;k++){
	pX[k]=std::log(QQ1[k][i][j]);
	pY[k]=QQ2[k][i][j];
      }
      G4PenelopeInterpolator* interpolator1 = new G4PenelopeInterpolator(pZ,pX,NumberofZPoints);
      Q1[i][j]=std::exp(interpolator1->CubicSplineInterpolation((G4double) Zmat));
      delete interpolator1;
      G4PenelopeInterpolator* interpolator2 = new G4PenelopeInterpolator(pZ,pY,NumberofZPoints);    
      Q2[i][j]=interpolator2->CubicSplineInterpolation((G4double) Zmat);
      delete interpolator2;
    }
  }
 
  
  //std::ofstream fil("matrice.dat",std::ios::app);
  //fil << "Numero atomico: " << Zmat << G4endl;
  //for (i=0;i<NumberofEPoints;i++)
  //{
  //  fil << Q1[i][0] << " " << Q1[i][1] << " " << Q1[i][2] << " " << Q1[i][3] << G4endl;
  //}
  //fil.close();
  
}

void G4PenelopeBremsstrahlungAngular::InterpolationForK()
{
  G4double pE[NumberofEPoints] = {1.0e-03,5.0e-03,1.0e-02,5.0e-02,1.0e-01,5.0e-01};
  G4double pK[NumberofKPoints] = {0.0,0.6,0.8,0.95};
  G4double ppK[reducedEnergyGrid];
  G4double pX[NumberofKPoints];
  G4int i,j;

  for(i=0;i<reducedEnergyGrid;i++){
    ppK[i]=((G4double) i) * 0.05;
  }

  for(i=0;i<NumberofEPoints;i++){
    betas[i]=std::sqrt(pE[i]*(pE[i]+2*electron_mass_c2))/(pE[i]+electron_mass_c2);
  }

  for (i=0;i<NumberofEPoints;i++){
    for (j=0;j<NumberofKPoints;j++){
      Q1[i][j]=Q1[i][j]/((G4double) Zmat);
    }
  }

  //Expanded table of distribution parameters
  for (i=0;i<NumberofEPoints;i++){
    for (j=0;j<NumberofKPoints;j++){
      pX[j]=std::log(Q1[i][j]); //logarithmic 
    }
    G4PenelopeInterpolator* interpolator = new G4PenelopeInterpolator(pK,pX,NumberofKPoints);
    for (j=0;j<reducedEnergyGrid;j++){
      Q1E[i][j]=interpolator->CubicSplineInterpolation(ppK[j]);
    }
    delete interpolator;
    for (j=0;j<NumberofKPoints;j++){
      pX[j]=Q2[i][j];
    }
    G4PenelopeInterpolator* interpolator2 = new G4PenelopeInterpolator(pK,pX,NumberofKPoints);
    for (j=0;j<reducedEnergyGrid;j++){
      Q2E[i][j]=interpolator2->CubicSplineInterpolation(ppK[j]);
    }
    delete interpolator2;
  }
}

G4double G4PenelopeBremsstrahlungAngular::ExtractCosTheta(G4double e1,G4double e2)
{
  //e1 = kinetic energy of the electron
  //e2 = energy of the bremsstrahlung photon

  G4double beta = std::sqrt(e1*(e1+2*electron_mass_c2))/(e1+electron_mass_c2);
  


  G4double RK=20.0*e2/e1;
  G4int ik=std::min((G4int) RK,19);
  
  G4double P10=0,P11=0,P1=0;
  G4double P20=0,P21=0,P2=0;
  G4double pX[NumberofEPoints];
  //First coefficient
  G4int i;
  G4int j = ik;
  for (i=0;i<NumberofEPoints;i++){
    pX[i]=Q1E[i][j];
  }
  G4PenelopeInterpolator* interpolator = new G4PenelopeInterpolator(betas,pX,NumberofEPoints);
  P10=interpolator->CubicSplineInterpolation(beta);
  delete interpolator;
  j++; //(j=ik+1)
  for (i=0;i<NumberofEPoints;i++){
    pX[i]=Q1E[i][j];
  }
  G4PenelopeInterpolator* interpolator2 = new G4PenelopeInterpolator(betas,pX,NumberofEPoints);
  P11=interpolator2->CubicSplineInterpolation(beta);
  delete interpolator2;
  P1=P10+(RK-(G4double) ik)*(P11-P10);
  
  //Second coefficient
  j = ik;
  for (i=0;i<NumberofEPoints;i++){
    pX[i]=Q2E[i][j];
  }
  G4PenelopeInterpolator* interpolator3 = new G4PenelopeInterpolator(betas,pX,NumberofEPoints);
  P20=interpolator3->CubicSplineInterpolation(beta);
  delete interpolator3;
  j++; //(j=ik+1)
  for (i=0;i<NumberofEPoints;i++){
    pX[i]=Q2E[i][j];
  }
  G4PenelopeInterpolator* interpolator4 = new G4PenelopeInterpolator(betas,pX,NumberofEPoints);
  P21=interpolator4->CubicSplineInterpolation(beta);
  delete interpolator4;
  P2=P20+(RK-(G4double) ik)*(P21-P20);
  
  //Sampling from the Lorenz-trasformed dipole distributions
  P1=std::min(std::exp(P1)/beta,1.0);
  G4double betap = std::min(std::max(beta*(1.0+P2/beta),0.0),0.9999);
  
  G4double cdt=0,testf=0;
  
  if (G4UniformRand() < P1){
    do{
      cdt = 2.0*G4UniformRand()-1.0;
      testf=2.0*G4UniformRand()-(1.0+cdt*cdt);
    }while(testf>0);
  }
  else{
    do{
      cdt = 2.0*G4UniformRand()-1.0;
      testf=G4UniformRand()-(1.0-cdt*cdt);
    }while(testf>0);
  }
  cdt = (cdt+betap)/(1.0+betap*cdt);
  return cdt;
}
