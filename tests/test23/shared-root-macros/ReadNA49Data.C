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

// for all types of secondaries except pions
//

int    NPoints         = 0;

float* xF              = 0;

float* dXSecdxF     = 0;
float* err_dXSecdxF = 0;
float* dNdxF           = 0;
float* err_dNdxF       = 0;
float* pT              = 0;
float* err_pT          = 0;
float* pT2             = 0;
float* err_pT2         = 0;

float* XSec            = 0;
float* err_XSec        = 0;

// for pions (pi+ & pi-)
//
float* Y =0;
float* dNdY[2]     = { 0, 0 };
float* err_dNdY[2] = { 0, 0 };


void readIntegratedSpectra( std::string beam, std::string target, std::string secondary )
{

   std::string dirname = "./na49-exp-data/";
   
   std::string filename = beam + "_" + target + "_" + secondary;
   
   std::string filenam1;
   
   if ( secondary == "pion" )
   { 
      filenam1 = filename + "_dndy.dat";
   }
   else
   {
      filename1 = filename + "_integr.dat";
   }
   
   std::string file1 = dirname + filename1;
   
   ifstream infile1;
   infile1.open( file1.c_str() );

   infile1 >> NPoints;
   
   if ( secondary == "pion" )
   {
      if (Y) delete Y;
      Y = new float[NPoints];
      if ( dNdY[0] = 0 ) delete dNdY[0];
      if ( dNdY[1] = 0 ) delete dNdY[1];
      dNdY[0] = new float[NPoints]; // pi+
      dNdY[1] = new float[NPoints]; // pi-
      if ( err_dNdY[0] = 0 ) delete err_dNdY[0];
      if ( err_dNdY[1] = 0 ) delete err_dNdY[1];      
      err_dNdY[0] = new float[NPoints];
      err_dNdY[1] = new float[NPoints];
      for ( int i=0; i<NPoints; i++ )
      {
         infile1 >> Y[i] >> dNdY[0][i] >> dNdY[1][i];
	 err_dNdY[0][i] = dNdY[0][i]*0.02; // errors are said to be 2%
	 err_dNdY[1][i] = dNdY[1][i]*0.02;
      } 
      return;     
   }
   

   if (xF) delete [] xF;
   xF = new float[NPoints];
   
  
   if (dNdxF) delete dNdxF;
   dNdxF = new float[NPoints];
   
   if (err_dNdxF) delete err_dNdxF;
   err_dNdxF = new float[NPoints];
   
   if ( secondary != "neutron" )
   {
     if  (dXSecdxF) delete [] dXSecdxF;
     dXSecdxF = new float[NPoints];
     if  (err_dXSecdxF) delete [] err_dXSecdxF;
     err_dXSecdxF = new float[NPoints];
     if (pT) delete [] pT;
     pT = new float[NPoints];
     if ( err_pT ) delete [] err_pT;
     err_pT = new float[NPoints];
     if (pT2) delete [] pT2;
     pT2 = new float[NPoints];
     if ( err_pT2 ) delete [] err_pT2;
     err_pT2 = new float[NPoints];
   }

   for ( int i=0; i<NPoints; i++ )
   {
      if ( secondary == "neutron" )
      {
         infile1 >> xF[i] >> dNdxF[i] >> err_dNdxF[i];
	 err_dNdxF[i] *= 0.01;
      }
      else
      {
         infile1 >> xF[i] >> dXSecdxF[i] >> err_dXSecdxF[i]  >> dNdxF[i] >> err_dNdxF[i]  
	         >> pT[i] >> err_pT[i] >> pT2[i] >> err_pT2[i];
         err_dXSecdxF[i]    *= 0.01;
	 err_dNdxF[i]       *= 0.01;
	 err_pT[i]          *= 0.01;
	 err_pT2[i]         *= 0.01;
      }
   }
   
   return;

}

