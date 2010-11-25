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
// $Id: G4PenelopeBremsstrahlungContinuous.cc,v 1.13 2010-11-25 09:43:47 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// --------------------------------------------------------------
//
// File name:     G4PenelopeBremsstrahlungContinuous
//
// Author:        Luciano Pandola
// 
// Creation date: February 2003
// History:
// -----------
// 20 Feb 2003  L. Pandola       1st implementation
// 17 Mar 2003  L. Pandola       Added the correction for positrons
// 19 Mar 2003  L. Pandola       Bugs fixed
// 17 Mar 2004  L. Pandola       Removed unnecessary calls to std::pow(a,b)
// 09 Dec 2008  L. Pandola       Update ReadFile() in a way to make explicit use of the units of 
//                               measurement. Also some cosmetics to improve readibility.
//----------------------------------------------------------------

#include "G4PenelopeBremsstrahlungContinuous.hh"
#include "G4PenelopeInterpolator.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include <sstream>

G4PenelopeBremsstrahlungContinuous::G4PenelopeBremsstrahlungContinuous (G4int Zed,G4double cut,
									G4double e1,
									G4double e2,
									const G4String name)  : 
  Zmat(Zed),tCut(cut),MinE(e1),MaxE(e2),partName(name)
{
  //Construct extended energy table      
  //200 bins between MinE and MaxE (logarithmic)
  G4double EL=0.99999*MinE;
  G4double EU=1.00001*MaxE;
  DLFC=std::log(EU/EL)/((G4double) (NumberofExtendedEGrid-1));
  ExtendedLogEnergy[0]=std::log(EL);
  for (size_t i=1;i<NumberofExtendedEGrid;i++){
    ExtendedLogEnergy[i]=ExtendedLogEnergy[i-1]+DLFC;
  }
  DLFC=1.0/DLFC;

  LoadFromFile();
  PrepareInterpolationTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopeBremsstrahlungContinuous::~G4PenelopeBremsstrahlungContinuous()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungContinuous::LoadFromFile()
{
  //Read information from DataBase File
 char* path = getenv("G4LEDATA");
 if (!path)
   {
     G4String excep = "G4PenelopeBremsstrahlungContinuous - G4LEDATA environment variable not set!";
     G4Exception(excep);
     return;
   }
 G4String pathString(path);
 G4String filename = "br-pen-cont-";
 std::ostringstream ost; 
 ost << filename << Zmat << ".dat";
 G4String name(ost.str());
 G4String dirFile = pathString + "/penelope/" + name;
 std::ifstream file(dirFile);
 if (!file.is_open())
     {
      G4String excep = "G4PenelopeBremsstrahlungContinuous - data file " + name + " not found!";
      G4Exception(excep);
     }
 G4double a1 = -1.;
 for (size_t i=0;i<NumberofEPoints;i++)
   {
     file >> a1;
     //1) reads energy in MeV
     Energies[i]=a1*MeV;
     //2) read 32 cross sections in cm2
     for (size_t j=0;j<NumberofKPoints;j++){
       file >> a1;
       ReducedCS[i][j]=a1*cm2; 
     }
     //3) read the total cross section, in cm2
     file >> a1;
     TotalCS[i]=a1*cm2; 
     // Check closing item
     file >> a1;
     if (a1 != ((G4double) -1))
       {
	 G4String excep = "G4PenelopeBremsstrahlungContinuous - Check the bremms data file " 
	   + name;
	 G4Exception(excep);
       }
   }
 file.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopeBremsstrahlungContinuous::PrepareInterpolationTable()
{
  //The energy-loss spectrum is re-normalized to reproduce the
  //total scaled cross-section of Berger and Seltzer

  G4double pK[NumberofKPoints] = {1.0e-12,0.05,0.075,0.1,0.125,0.15,0.2,0.25,
				  0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,
				  0.75,0.8,0.85,0.9,0.925,0.95,0.97,0.99,
				  0.995,0.999,0.9995,0.9999,0.99995,0.99999,1.0};
  G4double pY[NumberofKPoints];

  size_t i=0,j=0;
  for (i=0;i<NumberofEPoints;i++){
    for (j=0;j<NumberofKPoints;j++){
      pY[j]= ReducedCS[i][j];
    }
    G4PenelopeInterpolator* interpolator = new G4PenelopeInterpolator(pK,pY,NumberofKPoints);
    G4double Rsum = interpolator->CalculateMomentum(1.0,0);
   
    delete interpolator;
    G4double Fact = (Energies[i]+electron_mass_c2)*(1.0/fine_structure_const)*millibarn/
      (classic_electr_radius*classic_electr_radius*(Energies[i]+2.0*electron_mass_c2));
    G4double Normalization = TotalCS[i]/(Rsum*Fact);
    G4double TST = std::abs(Normalization-1.0);
    if (TST > 1.0) {
      G4String excep = "G4PenelopeBremsstrahlungContinuous - Check the bremms data file";
      G4Exception(excep);
    }
    for (j=0;j<NumberofKPoints;j++){
      ReducedCS[i][j] = ReducedCS[i][j]*Normalization;
    }
  }

  //Compute the scaled energy loss distribution and sampling parameters
  // for the energies in the simulation grid
  
  // Interpolation in E
  G4double pX[NumberofEPoints];
  G4double pYY[NumberofEPoints];
  for (i=0;i<NumberofEPoints;i++){
    pX[i] = std::log(Energies[i]);
  }
 
  for (j=0;j<NumberofKPoints;j++){
    for (i=0;i<NumberofEPoints;i++){
      pYY[i] = std::log(ReducedCS[i][j]);
    }
    G4PenelopeInterpolator* interpolator2 = new G4PenelopeInterpolator(pX,pYY,NumberofEPoints);
    for (i=0;i<NumberofExtendedEGrid;i++){
      G4double ELL = ExtendedLogEnergy[i];
      if (ELL >= pX[0]) {
	p0[i][j] = std::exp(interpolator2->CubicSplineInterpolation(ELL));
      }
      else
	{
	  G4double F1=interpolator2->CubicSplineInterpolation(pX[0]);
	  G4double FP1 = interpolator2->FirstDerivative(pX[0]);
	  p0[i][j] = std::exp(F1+FP1*(ELL-pX[0]));
	}
    }
    delete interpolator2;
  }
  
  //All this stuff might be useful later on 
  /*
  G4double PDF[NumberofKPoints];
  
  for (i=0;i<NumberofExtendedEGrid;i++){
    for (j=0;j<NumberofKPoints;j++){
      PDF[j]=p0[i][j];
    } 
    G4double Xc=0;
    if (i<(NumberofExtendedEGrid-1)){
      Xc=tCut/std::exp(ExtendedLogEnergy[i+1]);
    }
    else
      {
	Xc=tCut/std::exp(ExtendedLogEnergy[NumberofExtendedEGrid-1]);
      }
    
    G4PenelopeInterpolator* interpolator3 = new G4PenelopeInterpolator(pK,PDF,NumberofKPoints);
    Pbcut[i]=interpolator3->CalculateMomentum(Xc,-1);
    delete interpolator3;
  }
  */
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
		        
G4double G4PenelopeBremsstrahlungContinuous::CalculateStopping(G4double energy)
  //Stopping power expressed in MeV/mm*2
{
  G4double Xel=std::max(std::log(energy),ExtendedLogEnergy[0]);
  G4double Xe=1.0+(Xel-ExtendedLogEnergy[0])*DLFC;
  G4int Ke = (G4int) Xe; 
  G4double Xek = Xe-Ke; 
  
  //Global x-section factor
  G4double Fact=Zmat*Zmat*(energy+electron_mass_c2)*(energy+electron_mass_c2)/
    (energy*(energy+2.0*electron_mass_c2));
  Fact *= PositronCorrection(energy);

  //Moments of the scaled bremss x-section
  G4double wcre = tCut/energy;
  G4double pY[NumberofKPoints];
  G4double pK[NumberofKPoints] = {1.0e-12,0.05,0.075,0.1,0.125,0.15,0.2,0.25,
				  0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,
				  0.75,0.8,0.85,0.9,0.925,0.95,0.97,0.99,
				  0.995,0.999,0.9995,0.9999,0.99995,0.99999,1.0};

  for (size_t i=0;i<NumberofKPoints;i++){
    pY[i] = p0[Ke][i];
  }
  G4PenelopeInterpolator* interpolator1 = new G4PenelopeInterpolator(pK,pY,NumberofKPoints);
  G4double XS1A = interpolator1->CalculateMomentum(wcre,0);
  for (size_t k=0;k<NumberofKPoints;k++){
    pY[k] = p0[std::min(Ke+1,(G4int) NumberofExtendedEGrid-1)][k];
  }
  G4PenelopeInterpolator* interpolator2 = new G4PenelopeInterpolator (pK,pY,NumberofKPoints);
  G4double XS1B = interpolator2->CalculateMomentum(wcre,0);
 
  G4double XS1 = ((1.0-Xek)*XS1A+Xek*XS1B)*Fact*energy; //weighted mean between the energy bin of the grid

  //This is the 2nd momentum (straggling cross section). It might be useful later on
  /*
    G4double XS2A = interpolator1->CalculateMomentum(wcre,1);
    G4double XS2 = ((1.0-Xek)*XS2A+Xek*XS2B)*Fact*energy*energy; 
    G4double XS2B = interpolator2->CalculateMomentum(wcre,1);
  */

  delete interpolator1;
  delete interpolator2;

  return XS1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeBremsstrahlungContinuous::PositronCorrection(G4double energy)
{
  const G4double Coeff[7]={-1.2359e-01,6.1274e-2,-3.1516e-2,7.7446e-3,-1.0595e-3,
			   7.0568e-5,-1.8080e-6};
  G4double correct=0;
  if (partName == "e-") {
    return 1.0; //no correction for electrons
  }
  else if (partName == "e+"){
    G4double T=std::log(1+((1e6*energy)/(Zmat*Zmat*electron_mass_c2)));
    for (G4int i=0;i<7;i++){
      correct += Coeff[i]*std::pow(T,i+1);
    }
    correct = 1.0-std::exp(correct);
    return correct;
  }
  else //neither electrons nor positrons...exception
    {
      G4String excep = "G4PenelopeBremmstrahlungContinuous: the particle is not e- nor e+!";
      G4Exception(excep);
      return 0;
    }
}
