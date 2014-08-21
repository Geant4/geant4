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

const int NModels = 3;
std::string ModelName[3]  = { "ftfp_bert", "NuBeam", "NuBeam-7-10" };
// std::string ModelName[3]  = { "ftfp", "qgsp-g4lund-str-fragm", "qgsp" };  
// std::string ModelName[3]  = { "ftfp", "bertini", "qgsp_g4lund-str-fragm" };  
// std::string ModelName[2]  = { "NuBeam", "ftfp_bert" };  
// std::string ModelName[2]  = { "qgsp_bert-with-decays", "ftfp_bert-with-decays" };  
int         ColorModel[5] = { kMagenta, 7, kRed, kBlack, 14 }; // 14 = grey, 7 = light "sky"-blue

std::string energy[6] = { "3.0", "5.0", "8.0", "12.0", "31.0", "158.0" }
double table_pip[3][6] ;
double table_pim[3][6] ;
double table_p[3][6] ;


static int isNA49UtilLoaded = 0;
static int isNA61UtilLoaded = 0;
static int isHARPUtilLoaded = 0;
static int isChi2UtilLoaded = 0;


void calcChi2DDiffXSecNA49( std::string beam, std::string target, std::string secondary )
{

   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA49Data.C");
      isNA49UtilLoaded = 1;
   }
   if ( isChi2UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/Chi2Calc.C");
      isChi2UtilLoaded = 1;
   }
   
   readDDiffSpectra( beam, target, secondary );
      
   for ( int m=0; m<NModels; ++m )
   {
      int NDF = 0;
      double chi2 = Chi2DDiffXSecNA49( beam, target, secondary, ModelName[m], NDF );
      std::cout << " chi2/NDF = " << chi2 << "/" << NDF << " = " << chi2/NDF << " for model " << ModelName[m] << std::endl;
   }   

   return;

}


/*
void calcChi2MomSpectrumNA61( std::string beam, std::string target, std::string secondary )
{

   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA61Data.C");
      isNA61UtilLoaded = 1;
   }
   if ( isChi2UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/Chi2Calc.C");
      isChi2UtilLoaded = 1;
   }
         
   for ( int m=0; m<NModels; ++m )
   {
      double chi2 = Chi2MomSpectrumNA61( beam, target, secondary, ModelName[m] );
//      std::cout << " chi2 = " << chi2 << " for model " << ModelName[m] << std::endl;
   }   

   return;

}
*/

void calcChi2MomSpectrumNA61( std::string beam, std::string target, std::string secondary )
{

   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA61Data.C");
      isNA61UtilLoaded = 1;
   }
   if ( isChi2UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/Chi2Calc.C");
      isChi2UtilLoaded = 1;
   }
   
   for ( int m=0; m<NModels; ++m )
   {
      double chi2 = 0.;
      int NDF = 0;
//      readMomSpectrum( beam, target, secondary, "0", "20" );
//      chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, secondary, 
//                                           "0", "20",
//				           ModelName[m],
//					   NDF );
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
      gROOT->LoadMacro("../test23/shared-root-macros/ReadHARPData.C");
      isHARPUtilLoaded = 1;
   }
   if ( isChi2UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/Chi2Calc.C");
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

      std::cout << " chi2/NDF = " << chi2 << "/" << NDF << " = " << chi2/NDF << " for model " << ModelName[m] << std::endl;

   }
   
   return;

}

/* This does NOT work quite yet - naming conflicts (NPonts, etc...)

void Scan( std::string beam, std::string target, bool doHARP, bool doNA61, bool doNA49 )
{

   if ( isNA49UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA49Data.C");
      isNA49UtilLoaded = 1;
   }
   if ( isNA61UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadNA61Data.C");
      isNA61UtilLoaded = 1;
   }
   if ( isHARPUtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/ReadHARPData.C");
      isHARPUtilLoaded = 1;
   }

   if ( isChi2UtilLoaded <= 0 )
   {
      gROOT->LoadMacro("../test23/shared-root-macros/Chi2Calc.C");
      isChi2UtilLoaded = 1;
   }
  
   double chi2 = 0.;
   int NDF = 0;
   
   for ( int i=0; i<3; ++i )
   {
      for ( int j=0; j<6; ++j )
      {
         table_pip[i][j] = -1;
	 table_pim[i][j] = -1;
	 table_p[i][j]   = -1;
      }
   }

   for ( int im=0; im<NModels; ++im )
   {

      // first, do HARP data (3.0-12.0GeV)
      //  

      if ( doHARP )
      {

         if ( ModelName[im].find("qgs") == std::string::npos ) // check if NOT QGS (not valid in this E-range)
         {
            for ( int ie=0; ie<4; ++ie )
            {
         
	       NDF=0;
	       ReadHARPData( beam, target, energy[ie], "piplus", "FW" );   
               chi2 += Chi2MomSpectrumHARP( "proton", target, "3.0", "piplus", "FW", ModelName[im], NDF );
               ReadHARPData( beam, target, energy[ie], "piplus", "LA" );   
               chi2 += Chi2MomSpectrumHARP( beam, target, energy[ie], "piplus", "LA", ModelName[im], NDF );
	       table_pip[im][ie] = chi2/NDF;
	 
	       NDF=0;
	       ReadHARPData( beam, target, energy[ie], "piminus", "FW" );   
               chi2 += Chi2MomSpectrumHARP( beam, target, "3.0", "piminus", "FW", ModelName[im], NDF );
               ReadHARPData( beam, target, energy[ie], "piminus", "LA" );   
               chi2 += Chi2MomSpectrumHARP( beam, target, energy[ie], "piminus", "LA", ModelName[im], NDF );
	       table_pim[im][ie] = chi2/NDF;

	       NDF=0;
	       ReadHARPData( beam, target, energy[ie], "proton", "FW" );   
               chi2 += Chi2MomSpectrumHARP( beam, target, "3.0", "proton", "FW", ModelName[im], NDF );
               ReadHARPData( beam, target, energy[ie], "proton", "LA" );   
               chi2 += Chi2MomSpectrumHARP( beam, target, energy[ie], "proton", "LA", ModelName[im], NDF );
	       table_p[im][ie] = chi2/NDF;
	 
            }
         } // endif - check if NOT QGS

      } // endif doHARP

      if ( ModelName.find("bertini") != std::string::npos ) continue; // not valid above 10-12GeV anyway...

      // now do NA61 data (31.0GeV)
      //
      if ( doNA61 )
      {
         chi2 = 0.;
         NDF = 0;
         readMomSpectrum( beam, target, "piplus", "0", "20" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piplus", 
                                           "0", "20",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piplus", "20", "40" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piplus", 
                                           "20", "40",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piplus", "40", "60" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piplus", 
                                           "40", "60",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piplus", "60", "100" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piplus", 
                                           "60", "100",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piplus", "100", "140" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piplus", 
                                           "100", "140",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piplus", "140", "180" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piplus", 
                                           "140", "180",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piplus", "180", "240" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piplus", 
                                           "180", "240",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piplus", "240", "300" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piplus", 
                                           "240", "300",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piplus", "300", "360" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piplus", "300", "360",
				           ModelName[m], NDF );
         table_pip[im][4] = chi2/NDF;
      
         chi2 = 0.;
         NDF = 0;
         readMomSpectrum( beam, target, "piminus", "0", "20" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piminus", 
                                           "0", "20",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piminus", "20", "40" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piminus", 
                                           "20", "40",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piminus", "40", "60" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piminus", 
                                           "40", "60",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piminus", "60", "100" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piminus", 
                                           "60", "100",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piminus", "100", "140" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piminus", 
                                           "100", "140",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piminus", "140", "180" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piminus", 
                                           "140", "180",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piminus", "180", "240" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piminus", 
                                           "180", "240",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piminus", "240", "300" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piminus", 
                                           "240", "300",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "piminus", "300", "360" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "piminus", "300", "360",
				           ModelName[m], NDF );
         table_pim[im][4] = chi2/NDF;

         chi2 = 0.;
         NDF = 0;
         readMomSpectrum( beam, target, "proton", "0", "20" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "proton", 
                                           "0", "20",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "proton", "20", "40" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "proton", 
                                           "20", "40",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "proton", "40", "60" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "proton", 
                                           "40", "60",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "proton", "60", "100" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "proton", 
                                           "60", "100",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "proton", "100", "140" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "proton", 
                                           "100", "140",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "proton", "140", "180" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "proton", 
                                           "140", "180",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "proton", "180", "240" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "proton", 
                                           "180", "240",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "proton", "240", "300" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "proton", 
                                           "240", "300",
				           ModelName[m], NDF );
         readMomSpectrum( beam, target, "proton", "300", "360" );
         chi2 += Chi2MomSpectrumNA61ThetaBin( beam, target, "proton", "300", "360",
				           ModelName[m], NDF );
         table_p[im][4] = chi2/NDF;
            
      } // end-if doNA61

      if ( doNA49 )
      {

         NDF = 0;
         chi2 = 0.;
         chi2 = Chi2DDiffXSecNA49( beam, target, "piplus", ModelName[m], NDF );
         table_pip[im][5] = chi2/NDF;
      
         NDF = 0;
         chi2 = 0.;
         chi2 = Chi2DDiffXSecNA49( beam, target, "piminus", ModelName[m], NDF );
         table_pim[im][5] = chi2/NDF;

      }

   }
   
//   std::cout << "beam= " << beam << "+" << "target=" << target << " --> " << "piplus+X" << std::endl;    
   for ( int i=0; i<3; ++i )
   {
      for ( int j=0; j<6; ++j )
      {
         std::cout << " " << table_pip[i][j];
      }
      std::cout << std::endl;
   }
//   std::cout << "beam= " << beam << "+" << "target=" << target << " --> " << "piminus+X" << std::endl;    
   for ( int i=0; i<3; ++i )
   {
      for ( int j=0; j<6; ++j )
      {
         std::cout << " " << table_pim[i][j];
      }
      std::cout << std::endl;
   }
   
   return;

}
*/
