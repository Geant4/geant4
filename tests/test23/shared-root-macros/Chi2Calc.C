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


// FIXME: This assumes that either ReadHARP, or ReadNA61, or ReadNA49 utility is loaded,
//        along with the declaration of the variables and arrays it uses for data !!! 
//        Need to improve the infrastructure and maybe make it more modular.
//        

double Chi2DDiffXSecNA49( std::string beam, std::string target, std::string secondary, std::string model, int& NDF )
{
   
   double Chi2 = 0;
      
   std::string histofile = "./na49-histo/" + beam + target + "158.0GeV" + model + ".root"; 
   TFile* f = new TFile( histofile.c_str() );

   for ( int icount=0; icount<NSubSets; ++icount )
   {
      std::ostringstream cnt;
      cnt << (icount+7);
      std::string hname = "pT";
      if ( secondary == "piplus" )
      {
	    hname +="pip";
      }
      else if ( secondary == "piminus" )
      {
	    hname += "pim";
      }
      hname += cnt.str();
      TH1F* h = (TH1F*)f->Get( hname.c_str() );
      h->Scale(226.); // or should it be xsec=251. as in Geant4 ???   
      int NEntries = DDiffDataHolder[icount].GetEntries();
      int NX = h->GetNbinsX();
      for ( int k=0; k<=NX; ++k ) 
      {
         double xx1 = h->GetBinLowEdge(k);
	 double xx2 = h->GetBinWidth(k);
	 for (int kk=0; kk<NEntries; ++kk )
	 {
	    TVector3* tmpVec = DDiffDataHolder[icount].At(kk);
	    if ( xx1 < tmpVec->x() && xx1+xx2 > tmpVec->x() )
	    {
	       double yy1  = h->GetBinContent(k);
	       double eyy1 = h->GetBinError(k);
	       if ( ( eyy1*eyy1 + tmpVec->z()*tmpVec->z() ) > 0 )
	       {
	          Chi2 += ( yy1 - tmpVec->y() )*( yy1 - tmpVec->y() ) / ( eyy1*eyy1 + tmpVec->z()*tmpVec->z() );
	          // Chi2 += ( yy1 - tmpVec->y() )*( yy1 - tmpVec->y() ) / tmpVec->y();
	          ++NDF;
	       } 
	       break;
	    }
	 }
      }      
   } 
   
   // NDF -= 1;  
   // std::cout << " chi2/NDF = " << Chi2 << "/" << NDF << " for model " << model << std::endl;
   
   return Chi2;
   
}

double Chi2MomSpectrumNA61ThetaBin( std::string beam, std::string target, std::string secondary, 
                                    std::string sec_angle_min, std::string sec_angle_max,
			            std::string model,
				    int& NDF )
{

   double Chi2 = 0.;

   std::string histofile = "./na61-histo/" + beam + target + "31.0GeV" + model + ".root"; 
      
   TFile* f = new TFile( histofile.c_str() );
      
   std::string histoname = secondary + "Mult_" + sec_angle_min + "_" + sec_angle_max;
   
   TH1F* h = (TH1F*)f->Get( histoname.c_str() );
   int NX = h->GetNbinsX();

   for ( int k=0; k<=NX; k++ ) 
   {
      double xx1 = h->GetBinLowEdge(k);
      double xx2 = h->GetBinWidth(k);
      for (int kk=0; kk<NPoints; kk++ )
      {
	 if ( xx1 < SecMom[kk] && xx1+xx2 > SecMom[kk] )
	 {
	    double yy1  = h->GetBinContent(k);
	    double eyy1 = h->GetBinError(k);
	    // Chi2 += ( yy1 - SecSigma[kk] )*( yy1 - SecSigma[kk] ) / SecSigma[kk];
	    if ( ( eyy1*eyy1 + SecESigma[kk]*SecESigma[kk] ) > 0 )
	    {
	       Chi2 += ( yy1 - SecSigma[kk] )*( yy1 - SecSigma[kk] ) / ( eyy1*eyy1 + SecESigma[kk]*SecESigma[kk] );
	       ++NDF;
	    }
	    break; 
	 }
      }
   }   

   return Chi2; 

}

double Chi2MomSpectrumHARP( std::string beam, std::string target, std::string energy,
                            std::string secondary, 
                            std::string region, // must be "FW" or "LA"
			    std::string model,
			    int& NDF )
{

   double Chi2 = 0.;
      
   std::string histofile = "";
      
   histofile = "./harp-histo/";      
   histofile += ( beam + target + energy + "GeV" + model + ".root" ); 
            
   TFile* f = new TFile( histofile.c_str() );
   
   for ( int i=0; i<NSets; ++i )
   {
      
      char buf[5];
      sprintf( buf, "%i", i );      
      std::string histoname = secondary + "_" + region + "_";
      histoname.append( buf );
            
      TH1F* h = (TH1F*)f->Get( histoname.c_str() );
      
      int NX = h->GetNbinsX();
      
      for ( int k=0; k<=NX; k++ ) 
      {
         double xx1 = h->GetBinLowEdge(k);
         double xx2 = h->GetBinWidth(k);
         for (int kk=0; kk<NPoints[i]; kk++ )
         {
	    double xcenter = 0.5 * (XMin[i][kk]+XMax[i][kk]);
	    if ( xx1 < xcenter && xx1+xx2 > xcenter )
	    {
	       double yy1  = h->GetBinContent(k);
	       double eyy1 = h->GetBinError(k);
	       // Chi2 += ( yy1 - SecSigma[kk] )*( yy1 - SecSigma[kk] ) / SecSigma[kk];
	       if ( ( eyy1*eyy1 + EY[i][kk]*EY[i][kk] ) > 0 )
	       {
	          Chi2 += ( yy1 - Y[i][kk] )*( yy1 - Y[i][kk] ) / ( eyy1*eyy1 + EY[i][kk]*EY[i][kk] );
	          ++NDF;
	       }
	       break; 
	    }
         }
      }   
      
   }
   
   return Chi2;

}
