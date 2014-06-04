#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <list>

#include <math.h>
#include <vector>

#include "Rtypes.h"
#include "TROOT.h"
#include "TRint.h"
#include "TObject.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"

const int NModels = 2;
// std::string ModelName[2]  = { "ftfp", "qgsp" };  
std::string ModelName[2]  = { "NuBeam", "ftfp_bert" };  
// std::string ModelName[2]  = { "qgsp_bert-with-decays", "ftfp_bert-with-decays" };  
int         ColorModel[5] = { kMagenta, 7, kRed, kBlack, 14 }; // 14 = grey, 7 = light "sky"-blue

static int isNA49UtilLoaded = 0;
static int isNA61UtilLoaded = 0;
static int isHARPUtilLoaded = 0;
static int isChi2UtilLoaded = 0;

void calcChi2DDiffXSecNA49( std::string beam, std::string target, std::string secondary )
{

   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadNA49Data.C");
      isNA49UtilLoaded = 1;
   }
   if ( isChi2UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/Chi2Calc.C");
      isChi2UtilLoaded = 1;
   }
   
   readDDiffSpectra( beam, target, secondary );
      
   for ( int m=0; m<NModels; ++m )
   {
      int NDF = 0;
      double chi2 = Chi2DDiffXSecNA49( beam, target, secondary, ModelName[m], NDF );
      std::cout << " chi2/NDF = " << chi2 << "/" << NDF << " for model " << ModelName[m] << std::endl;
   }   

   return;

}

void calcChi2MomSpectrumNA61( std::string beam, std::string target, std::string secondary )
{

   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadNA61Data.C");
      isNA61UtilLoaded = 1;
   }
   if ( isChi2UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/Chi2Calc.C");
      isChi2UtilLoaded = 1;
   }
         
   for ( int m=0; m<NModels; ++m )
   {
      double chi2 = Chi2MomSpectrumNA61( beam, target, secondary, ModelName[m] );
//      std::cout << " chi2 = " << chi2 << " for model " << ModelName[m] << std::endl;
   }   

   return;

}


void calcChi2MomSpectrumNA61( std::string beam, std::string target, std::string secondary )
{

   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadNA61Data.C");
      isNA61UtilLoaded = 1;
   }
   if ( isChi2UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/Chi2Calc.C");
      isChi2UtilLoaded = 1;
   }
   
   for ( int m=0; m<NModels; ++m )
   {
      double chi2 = 0.;
      int NDF = 0;
      readMomSpectrum( beam, target, secondary, "0", "20" );
      chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, secondary, 
                                           "0", "20",
				           ModelName[m],
					   NDF );
      readMomSpectrum( beam, target, secondary, "20", "40" );
      chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, secondary, 
                                           "20", "40",
				           ModelName[m],
				           NDF );
      readMomSpectrum( beam, target, secondary, "40", "60" );
      chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, secondary, 
                                           "40", "60",
				           ModelName[m],
				           NDF );
      readMomSpectrum( beam, target, secondary, "60", "100" );
      chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, secondary, 
                                           "60", "100",
				           ModelName[m],
					   NDF );
      readMomSpectrum( beam, target, secondary, "100", "140" );
      chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, secondary, 
                                           "100", "140",
				           ModelName[m],
					   NDF );
      readMomSpectrum( beam, target, secondary, "140", "180" );
      chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, secondary, 
                                           "140", "180",
				           ModelName[m],
					   NDF );
      readMomSpectrum( beam, target, secondary, "180", "240" );
      chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, secondary, 
                                           "180", "240",
				           ModelName[m],
					   NDF );
      readMomSpectrum( beam, target, secondary, "240", "300" );
      chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, secondary, 
                                           "240", "300",
				           ModelName[m],
					   NDF );
      readMomSpectrum( beam, target, secondary, "300", "360" );
      chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, secondary, 
                                           "300", "360",
				           ModelName[m],
					   NDF );
//
// Note: in the simulated results it's 360-420 instead of 360-400... don't know why !
//
//      chi2 += Chi2MomSpectrumNA61( beam, target, secondary, 
//                                   "360", "400",
//				   ModelName[m] );

      std::cout << " chi2/NDF = " << chi2 << "/" << NDF <<  " for model " << ModelName[m] << std::endl;
   }   

   return;

} 

void calcChi2MomSpectrumHARP( std::string beam, std::string target, std::string energy, std::string secondary )
{

   if ( isHARPUtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/ReadHARPData.C");
      isHARPUtilLoaded = 1;
   }
   if ( isChi2UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("./shared-root-macros/Chi2Calc.C");
      isChi2UtilLoaded = 1;
   }

   double chi2 = 0.; 
      
   for ( int m=0; m<NModels; ++m )
   {
      
      double chi2 = 0.;
      int NDF = 0;

      ReadHARPData( beam, target, energy, secondary, "FW" );   
      chi2 += Chi2MomSpectrumHARP( beam, target, energy, secondary, "FW", ModelName[m], NDF );
      ReadHARPData( beam, target, energy, secondary, "LA" );   
      chi2 += Chi2MomSpectrumHARP( beam, target, energy, secondary, "LA", ModelName[m], NDF );

      std::cout << " chi2/NDF = " << chi2 << "/" << NDF <<  " for model " << ModelName[m] << std::endl;

   }
   
   return;

}
