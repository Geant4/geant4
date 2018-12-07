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
//
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "XrayFluoNormalization.hh"
#include "XrayFluoDataSet.hh"
#include "G4SystemOfUnits.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"

XrayFluoNormalization::XrayFluoNormalization()

{ }

XrayFluoNormalization::~XrayFluoNormalization()

{ }

const XrayFluoDataSet* XrayFluoNormalization::Normalize(G4double minIntExtreme, G4double maxIntExtreme, G4int nBins, G4String fileName)
{
 
 G4VDataSetAlgorithm* interpolation = new G4LogLogInterpolation();

 XrayFluoDataSet* dataSet = 
   new XrayFluoDataSet(1,fileName,interpolation,keV,1);
  
 G4double integral = Integrate(minIntExtreme, maxIntExtreme, nBins, dataSet);
 
 G4VDataSetAlgorithm* interpolation2 = new G4LogLogInterpolation();
 
 XrayFluoDataSet* finalDataSet = new XrayFluoDataSet(1,fileName,interpolation2,keV,1/(integral/keV));
 return finalDataSet;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double XrayFluoNormalization::Integrate(G4double minIntExtreme, G4double maxIntExtreme, G4int nBins, XrayFluoDataSet* dataSet)
{
 G4double partialHeight = 0;
 G4int id = 0;
 G4double lenghtOfBins = (maxIntExtreme - minIntExtreme)/nBins;
 
 for (G4int i = 0; i<nBins; i++)
   {
     G4double midPoint = minIntExtreme + i*lenghtOfBins+lenghtOfBins/2;
    
     G4double heightOfRectangle = dataSet->FindValue(midPoint,id);
    
     partialHeight += heightOfRectangle;
   
   }

 G4double integral = lenghtOfBins * partialHeight;
 
 delete dataSet;
 dataSet = 0;
 return integral;
}




