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

int    NPoints    = 0;

float* SecMom     = 0;

float* SecSigma   = 0; 
float* SecESigma  = 0;
float* SecESys    = 0;

float* K2PiRatio  = 0;
float* K2PiERatio = 0;

void readMomSpectrum( std::string beam, std::string target, 
                      std::string secondary, 
		      std::string sec_angle_min, std::string sec_angle_max )
{

   std::string dirname = "../test23/na61-exp-data/";
   std::string filename = beam + "_" + target + "_" + secondary;
   filename += "_" + sec_angle_min + "_" + sec_angle_max + ".dat";
   std::string file = dirname + filename;
   
   // std::cout << " file = " << file << std::endl;
   
   ifstream infile;
   infile.open( file.c_str() );
   
   infile >> NPoints;
   
   // std::cout << "NPoints = " << NPoints << std::endl;
   
   if ( SecMom ) delete [] SecMom;
   SecMom = new float[NPoints];
   
   if ( SecSigma ) delete [] SecSigma;
   SecSigma = new float[NPoints];
   
   if ( SecESigma ) delete [] SecESigma;
   SecESigma = new float[NPoints];
   
   if ( SecESys ) delete [] SecESys;
   SecESys = new float[NPoints];
   
   float pmin, pmax;
   for ( int i=0; i<NPoints; i++ )
   {  
      SecMom[i]    = 0.;
      SecSigma[i]  = 0.;
      SecESigma[i] = 0.;
      SecESys[i]   = 0.;
      
      if ( secondary != "proton" ) // FIXME !!! TMP STUFF untill I frigure out
                                   // about publication of the NA61's proton production data
      {
         infile >> pmin >> pmax >> SecSigma[i] >> SecESigma[i] >> SecESys[i];
      }
      else
      {
         infile >> pmin >> pmax >> SecSigma[i] >> SecESigma[i] ;
      }
      SecMom[i] = 0.5 * ( pmin + pmax );
      // std::cout << SecMom[i] << " " << SecSigma[i] << " " << SecESigma[i] << std::endl;
      if ( secondary != "proton" )
      {
         // scale to the total xsec (proton data come normalized)
	 SecSigma[i]  /= 229.3 ; // 229.3 is production sxec, while inelastic xsec = 257.2 
	                         // which one should be used for normalization ???        
         SecESigma[i] /= 229.3;
         SecESys[i]   /= 229.3 ;
      }
   }
      
   infile.close();
   
   return;

}

void readKPlus2PiPlusRatio( std::string beam, std::string target, 
		            std::string sec_angle_min, std::string sec_angle_max )
{

   std::string dirname = "../test23/na61-exp-data/";
   std::string filename = beam + "_" + target + "_kplus2piplus";
   filename += "_" + sec_angle_min + "_" + sec_angle_max + ".dat";
   std::string file = dirname + filename;
   
   // std::cout << " file = " << file << std::endl;
   
   ifstream infile;
   infile.open( file.c_str() );
   
   infile >> NPoints;
   
   if ( SecMom ) delete [] SecMom;
   SecMom = new float[NPoints];

   if ( K2PiRatio ) delete [] K2PiRatio;
   K2PiRatio = new float[NPoints];
   
   if ( K2PiERatio ) delete [] K2PiERatio;
   K2PiERatio = new float[NPoints];
   
   float pmin, pmax;
   for ( int i=0; i<NPoints; i++ )
   {  
      infile >> pmin >> pmax >> K2PiRatio[i] >> K2PiERatio[i] ;
      SecMom[i] = 0.5 * ( pmin + pmax );
      // std::cout << SecMom[i] << " " << K2PiRatio[i] << " " << K2PiERatio[i] << std::endl;
   }
   
   infile.close();

   return;

}
