#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
// #include <list>
#include <map>

#include <math.h>
#include <vector>

#include "Rtypes.h"
#include "TROOT.h"
#include "TRint.h"

// #include <boost/algorithm/string/split.hpp>
// #include <boost/algorithm/string/classification.hpp>

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

int NSubSets = 0;
double*     DDiff_xF = 0;
TObjArray*  DDiffDataHolder = 0;

void readDDiffSpectra( std::string beam, std::string target, std::string secondary )
{

// first of all, do the cleanups if necessary
   if ( DDiff_xF != 0 ) 
   {
      delete [] DDiff_xF;
      DDiff_xF = 0;
   }
   if ( DDiffDataHolder != 0 )
   {
      for ( int i=0; i<NSubSets; ++i )
      {
         DDiffDataHolder[i].Clear();
      }
      delete [] DDiffDataHolder;
      DDiffDataHolder = 0;
   }
   NSubSets = 0;
   
   // std::string dirname = "./na49-exp-data/";
   std::string dirname = "../test23/na49-exp-data/";
   
   std::string filename = beam + "_" + target + "_" + secondary;
   
   std::string filename1;
   
   filename1 = filename + "_ddiffxsec.dat";
   
   std::string file1 = dirname + filename1;

//   std::cout << " opening datafile " << file1 << std::endl;
         
   ifstream infile1;
   infile1.open( file1.c_str() );
   
   // std::string line = "";
   char line[256];
   for ( int i=0; i<256; ++i ) line[i] = ' ';
   std::vector<std::string> tokens;
   std::string del = " ";
   int counter = -1;
   
   while ( !infile1.eof()  )
   {

      for ( int i=0; i<256; ++i ) line[i] = ' '; // cleanup the read-in buffer before next iteration

      infile1.getline( line, 256 );

      std::string line1 = std::string( line ); // , 256 );
      
      if ( line1.find("#") == 0 ) // comment line
      {	 
	 continue;
      }
      if ( line1.find_first_not_of(del) == std::string::npos ) 
      {
	 continue;
      }
      if ( line1.find( "<NSubSets>" ) != std::string::npos )
      {
         infile1 >> NSubSets;
	 // Here I have to do some memory management !!!... or maybe not ???
	 DDiff_xF = new double[NSubSets];
	 DDiffDataHolder = new TObjArray[NSubSets];
	 continue; 
      }
      if ( line1.find("<nextXF>") != std::string::npos )
      {
	 for ( int i=0; i<line1.size(); ++i ) line[i] = ' ';
	 infile1.getline( line, 256  );
	 line1 = std::string( line );
	 // first line containing data in a subset
	 // now need to split into tokens
	 tokens.clear();
	 SplitString( line1, del, tokens );
	 if ( tokens.size() != 4 )
	 {
	    std::cout << " EMERGENCY !!!" << std::endl;
	 }
	 int sz = tokens.size();
	 counter++;
	 DDiff_xF[counter] = atof(tokens[0].c_str());
	 double x  = atof(tokens[1].c_str());
	 double y  = atof(tokens[2].c_str());
	 double ey =  atof(tokens[3].c_str()); // *y/100. ;
	 DDiffDataHolder[counter].Add( new TVector3( x, y, ey ) );
	 	 
	 continue;
      }
      // we get here for "normal" data line
      tokens.clear();
      SplitString( line1, del, tokens );
      DDiff_xF[counter] = atof(tokens[0].c_str());
      double x  = atof(tokens[1].c_str());
      double y  = atof(tokens[2].c_str());
      double ey =  atof(tokens[3].c_str()); // * y/100. ; // are errors given in percent ???
      DDiffDataHolder[counter].Add( new TVector3( x, y, ey ) );
   }

//   for ( int i=0; i<NSubSets; ++i )
//   {
//      std::cout << " size of subset " << i << " is " << DDiffDataHolder[i].GetEntries() << std::endl;
//   }
   
   return;

}

void SplitString( const std::string& str, const std::string del, std::vector<std::string>& tokens )
{
         
   string::size_type last = str.find_first_not_of(del, 0);
   string::size_type pos  = str.find_first_of(del, last);
   while (std::string::npos != pos || std::string::npos != last)
   {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(last, pos-last));
      // Skip delimiters.  Note the "not_of"
      last = str.find_first_not_of(del, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(del, last);
   }

   return;

}

void readIntegratedSpectra( std::string beam, std::string target, std::string secondary )
{

   // std::string dirname = "./na49-exp-data/";
   std::string dirname = "../test23/na49-exp-data/";
   
   std::string filename = beam + "_" + target + "_" + secondary;
      
   std::string filename1;
   
   if ( secondary == "pion" )
   { 
      filename1 = filename + "_dndy.dat";
   }
   else
   {
      filename1 = filename + "_integr.dat";
   }
   
   std::string file1 = dirname + filename1;
   
   std::cout << " input file = " << file1 << std::endl;
   
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
	 err_dNdY[0][i] = dNdY[0][i]*0.02; // common systematic errors are said to be 2%
	 err_dNdY[1][i] = dNdY[1][i]*0.02; // but that's for the dN/dy spectra only !!!
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

   float syserr = 0.;
   if ( beam == "H" )
   {
      syserr = 2; // in %
   }
   else
   {
      syserr = 2.5;
   }
   float err2 = 0.;
   for ( int i=0; i<NPoints; i++ )
   {
      if ( secondary == "neutron" )
      {
         infile1 >> xF[i] >> dNdxF[i] >> err_dNdxF[i];
	 err2 = err_dNdxF[i]*err_dNdxF[i] + syserr*syserr;
	 err_dNdxF[i]  = std::sqrt(err2)*0.01;
	 err_dNdxF[i] *= dNdxF[i];
      }
      else
      {
         infile1 >> xF[i] >> dXSecdxF[i] >> err_dXSecdxF[i]  >> dNdxF[i] >> err_dNdxF[i]  
	         >> pT[i] >> err_pT[i] >> pT2[i] >> err_pT2[i];
	 //
	 // here I also should account for the 2% (p+p) or 2.5% (p+C) systematics
	 //
	 err2 = err_dXSecdxF[i]*err_dXSecdxF[i] + syserr*syserr;
	 err_dXSecdxF[i]     = std::sqrt(err2)*0.01; 
	 err_dXSecdxF[i]    *= dXSecdxF[i];
         err2 = err_dNdxF[i]*err_dNdxF[i] + syserr*syserr;	 
	 err_dNdxF[i]        = std::sqrt(err2) * 0.01;
	 err_dNdxF[i]       *= dNdxF[i];
	 err2 = err_pT[i]*err_pT[i] + syserr*syserr;
	 err_pT[i]           = std::sqrt(err2)*0.01;
	 err_pT[i]          *= pT[i];
	 err2 = err_pT2[i]*err_pT2[i] + syserr*syserr;
	 err_pT2[i]          = std::sqrt(err2)*0.01;
	 err_pT2[i]         *= pT2[i];
      }
   }
   
   return;

}

