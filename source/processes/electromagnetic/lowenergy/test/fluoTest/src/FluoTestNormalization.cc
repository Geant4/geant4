#include "FluoTestNormalization.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "FluoTestDataSet.hh"

FluoTestNormalization::FluoTestNormalization()

{ }

FluoTestNormalization::~FluoTestNormalization()

{ }

const FluoTestDataSet* FluoTestNormalization::Normalize(G4double minIntExtreme, G4double maxIntExtreme, G4int nBins)
{
 
 G4String fileName = "C_flare";
 G4VDataSetAlgorithm* interpolation = new G4LogLogInterpolation();
 FluoTestDataSet* dataSet = 
   new FluoTestDataSet(1,fileName,interpolation,keV,1);
 G4double integral = Integrate(minIntExtreme, maxIntExtreme, nBins, dataSet);

 G4VDataSetAlgorithm* interpolation2 = new G4LogLogInterpolation();
 FluoTestDataSet* finalDataSet = new FluoTestDataSet(1,fileName,interpolation2,keV,1/(integral/keV));

 return finalDataSet;
}

G4double FluoTestNormalization::Integrate(G4double minIntExtreme, G4double maxIntExtreme, G4int nBins, FluoTestDataSet* dataSet)
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




