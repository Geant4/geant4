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
// $Id: XrayFluoAnalysisMessenger.cc
// GEANT4 tag $Name: 
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "XrayFluoNormalization.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "XrayFluoDataSet.hh"

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

G4double XrayFluoNormalization::Integrate(G4double minIntExtreme, G4double maxIntExtreme, G4int nBins, XrayFluoDataSet* dataSet)
{
 G4double partialHeight = 0;
 G4double id = 0;
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




