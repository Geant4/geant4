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
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TRefArray.h"
#include "TStyle.h"
#include "TGraph.h"

#include "TFileMerger.h"


// options for std::string sec_angle are:
// principle: "59.1", "119.0"  
// also available: ...

// beam energy should be given as "1.40", or "5.00", or "7.50"... 

// read-in exp.data business
//
int NPtKE = 0;
float BeamEn, SecAng;
// there's no particular reason for having size=30 
// usually, the number of exp.points is less than that, but... just to be on the safe side...
float KE[30], ExpValue[30], YY[30], Err1[30], Err2[30]; 

// regression test business
//
const int NVersions = 4;
std::string Versions[4] = { "geant4-09-06-ref07", "geant4-10-00-b01", "geant4-09-06-ref03", "geant4-09-06-p01" };
// std::string Versions[4] = { "geant4-09-04-ref10", "geant4-09-05", "geant4-09-05-p01", "geant4-09-05-ref02+MKtag"  };
// std::string Versions[4] = { "geant4-09-04-ref10", "geant4-09-05-ref01", "geant4-09-05-ref02+MKtag", "geant4-09-05-p01"  };
// int ColorVersion[4] = { kBlack, 7, kRed, kGreen }; // 7 = very light sky-blue
int ColorVersion[5] = { kRed, kGreen, 7, kBlack, 14 };


// model comparison business
//
 const int NModels = 3;
std::string Models[3] = { "bertini", "binary", "ftfp" };

//const int NModels = 2;
//std::string Models[2] = { "bertini", "ftfp" };

int ColorModel[6] = { 6, 3, 14 };

// --> General purspose exp.data read-in
//
void readKEData( std::string beam, std::string target, std::string energy, 
                 std::string secondary_type, std::string sec_angle )
{


   // read data
   std::string dirname = "./itep/" + beam + "/" + secondary_type + "/";
   std::string filename = target + energy + "GeV" + sec_angle + "deg.dat";
   std::string file = dirname + filename;
   std::cout << "reading from file: " << file << std::endl;
   
   ifstream infile;
   infile.open(file.c_str());
   infile >> BeamEn >> SecAng >> NPtKE;

   float staterr, syserr;
   for ( int ip=0; ip<NPtKE; ip++ )
   {
      infile >> KE[ip] >> ExpValue[ip] >> staterr >> syserr;
      syserr *= ExpValue[ip];
      YY[ip] = 1.;
      Err1[ip]  = sqrt(syserr*syserr+staterr*staterr);
      Err2[ip] = Err1[ip] / ExpValue[ip] ;
   }
   
   infile.close();

   return;

}

// --> Bertini regression business
//
void plotBertiniSummary( std::string beam, std::string target )
{
   

   TCanvas *myc1 = new TCanvas("myc1","",800,600);
   myc1->Divide(2,2);
   
   myc1->cd(1);
   plotMC2Data( beam, target, "1.40", "proton", "59.1", "bertini", 4 );
   
   myc1->cd(2);
   plotMC2Data( beam, target, "1.40", "neutron", "59.1", "bertini", 4 );
   
   myc1->cd(3);
   plotMC2Data( beam, target, "1.40", "proton", "119.0", "bertini",  4 );

   myc1->cd(4);
   plotMC2Data( beam, target, "1.40", "neutron", "119.0", "bertini", 4 );
   
   std::string en;
   if ( beam == "proton" )
   {
      en = "7.50";
   }
   else if ( beam == "piplus" || beam == "piminus" )
   {
      en = "5.00";
   }
   else
   {
      std::cout << " Nothing available for " << beam << " beam " << std::endl;
      return;
   }

   TCanvas *myc2 = new TCanvas("myc2","",800,600);
   myc2->Divide(2,2);
   
   
   myc2->cd(1);
   plotMC2Data( beam, target, en, "proton", "59.1", "bertini", 4 );
   
   myc2->cd(2);
   plotMC2Data( beam, target, en, "neutron", "59.1", "bertini", 4 );
   
   myc2->cd(3);
   plotMC2Data( beam, target, en, "proton", "119.0", "bertini", 4 );

   myc2->cd(4);
   plotMC2Data( beam, target, en, "neutron", "119.0", "bertini", 4 );

   myc2->cd();
   
   return;

}

// --> end of bertini regression business

// --> Binary regression business
//
void plotBinarySummary( std::string beam, std::string target )
{
   

   TCanvas *myc1 = new TCanvas("myc1","",800,600);
   myc1->Divide(2,2);
   
   myc1->cd(1);
   plotMC2Data( beam, target, "1.40", "proton", "59.1", "binary", 4 );
   
   myc1->cd(2);
   plotMC2Data( beam, target, "1.40", "neutron", "59.1", "binary", 4 );
   
   myc1->cd(3);
   plotMC2Data( beam, target, "1.40", "proton", "119.0", "binary",  4 );

   myc1->cd(4);
   plotMC2Data( beam, target, "1.40", "neutron", "119.0", "binary", 4 );
   
   myc1->cd();
   
   return;

}

// --> end of bertini regression business

// --> FTF (regression) business

void plotFTFSummary( std::string beam, std::string target )
{

   std::string en;
   if ( beam == "proton" )
   {
      en = "7.50";
   }
   else if ( beam == "piplus" || beam == "piminus" )
   {
      en = "5.00";
   }
   else
   {
      std::cout << " Nothing available for " << beam << " beam " << std::endl;
      return;
   }

   TCanvas *myc1 = new TCanvas("myc1","",800,600);
   myc1->Divide(2,2);

   myc1->cd(1);
   plotMC2Data( beam, target, en, "proton", "59.1", "ftfp" );
   myc1->cd(2);
   plotMC2Data( beam, target, en, "neutron", "59.1", "ftfp" );
   myc1->cd(3);
   plotMC2Data( beam, target, en, "proton", "119.0", "ftfp" );
   myc1->cd(4);
   plotMC2Data( beam, target, en, "neutron", "119.0", "ftfp" );

   myc1->cd();

   return;

}

// --> end of FTF business


// --> model comparison

void plotMC2Data( std::string beam, std::string target, std::string energy, 
                  std::string secondary_type, std::string sec_angle,
		  std::string model = "bertini",
	          int NVer=4 )
{
            
   readKEData( beam, target, energy, secondary_type, sec_angle );
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   for ( int ip=0; ip<NPtKE; ip++ )
   {
      if ( (YY[ip]+Err2[ip]) > ymax ) ymax = YY[ip]+Err2[ip];
      if ( (YY[ip]-Err2[ip]) < ymin ) ymin = YY[ip]-Err2[ip];
   }

   TGraph*  gr1 = new TGraphErrors(NPtKE,KE,YY,0,Err2);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   std::string xtitle = "Kinetic Energy of secondary " + secondary_type + " (GeV)";
   gr1->GetXaxis()->SetTitle( xtitle.c_str() );
   gr1->GetYaxis()->SetTitle("MC/Data");
   
   TGraph* gr[NVersions];
   std::string dir;
   if ( NVer > NVersions ) NVer = NVersions;
      
   for ( int iv=0; iv<NVer; iv++ )
   {

   // open G4 output root file (histo)
//   if ( iv == 0 ) 
//   { 
//      dir = "";
//   }
//   else
//   {
      dir = Versions[iv] + "/";
//   }
   std::string histofile = dir + beam + target + model + energy + "GeV.root";
   std::string histoname = "KE" + secondary_type + "0" + beam + target + model + energy + "GeV";
   int counter = 5 - sec_angle.length();
   for ( int il=0; il<counter; il++ )
   {
      histoname += " "; // pad up for the fact that initially sec.angle was supposed to be char[6] - no more, no less...
   }
   histoname += sec_angle;
   
   TFile* hfile = new TFile( histofile.c_str() );
   
   TH1F* hi = (TH1F*)hfile->Get( histoname.c_str() );
   
   //if ( iv == 1 || iv == 2 ) hi->Scale( (1./32.) );
   
   if ( hi == 0 || NPtKE <= 0 ) 
   {
      std::cout << " Invalid case for " << model << ": no exp.data or simulatiion, or both " << std::endl;
      return ;
   }
      
   
   float MC2DataX[30];
   float MC2DataY[30];
   float DX[30], DY[30];
   int np =0;
   
   int nx = hi->GetNbinsX();
   

   for ( int k=1; k <= nx; k++ )
   {        
        double xx1 = hi->GetBinLowEdge(k);
	double xx2 = hi->GetBinWidth(k);
	for (int kk=0; kk<NPtKE; kk++ )
	{
	   if ( xx1 < KE[kk] && xx1+xx2 > KE[kk] )
	   {
	      double yy = hi->GetBinContent(k);
	      MC2DataX[np] = KE[kk];
	      DX[np] = 0.;
	      MC2DataY[np] = yy / ExpValue[kk];
	      // also error calc here !...
	      DY[np]=Err1[kk]*MC2DataY[np]/ExpValue[kk];
	      if ( (MC2DataY[np]+DY[np]) > ymax ) ymax = MC2DataY[np]+DY[np];
	      if ( (MC2DataY[np]-DY[np]) < ymin ) ymin = MC2DataY[np]-DY[np];
	      np++;
	      break;
           }
	}
   }
   
   gr[iv] = new TGraphErrors( np, MC2DataX, MC2DataY, DX, DY );
   gr[iv]->SetTitle(hi->GetTitle());
   gr[iv]->GetXaxis()->SetTitle( xtitle.c_str() );
   gr[iv]->GetYaxis()->SetTitle("MC/Data");
   gr[iv]->SetMarkerColor(ColorVersion[iv]);  
   gr[iv]->SetMarkerStyle(21);
   gr[iv]->SetMarkerSize(1.6);
   if ( iv == 0 ) gr1->SetTitle(hi->GetTitle());
   
   } // end loop on versions
      
   
   if ( ymin > 0. ) ymin = 0;
   if ( ymax < 3. ) ymax = 3.;
   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->GetYaxis()->SetRangeUser(ymin, ymax);  
   
   gr1->Draw("apl");

   TLegend* leg = new TLegend(0.65, 0.75, 0.99, 0.99);   
   
   for ( int iv=0; iv<NVer; iv++ )
   {
      gr[iv]->Draw("lpsame");       
      leg->AddEntry( gr[iv], Versions[iv].c_str(), "p" );  
   }
    
   leg->AddEntry( gr1, "exp.data", "p");     

   leg->Draw();
   leg->SetFillColor(kWhite);
         
   return;

}


// --> models comparison


void plotForTalk()
{
   
   TCanvas *myc = new TCanvas("myc","",1200,800);
   myc->Divide(4,2);

   myc->cd(1);
   plotModelsMC2Data( "proton", "C", "7.50", "proton", "59.1" );
   myc->cd(2);
   plotModelsMC2Data( "proton", "U", "7.50", "proton", "59.1" );
   myc->cd(1);
   plotModelsMC2Data( "piminus", "C", "5.00", "proton", "59.1" );
   myc->cd(2);
   plotModelsMC2Data( "piminus", "U", "5.00", "proton", "59.1" );

   myc->cd(5);
   plotModelsMC2Data( "proton", "C", "7.50", "proton", "119.0" );
   myc->cd(6);
   plotModelsMC2Data( "proton", "U", "7.50", "proton", "119.0" );
   myc->cd(7);
   plotModelsMC2Data( "piminus", "C", "5.00", "proton", "119.0" );
   myc->cd(8);
   plotModelsMC2Data( "piminus", "U", "5.00", "proton", "119.0" );
  
   return;
}

// --> summary plots
//
void plotModelsMC2DataSummary( std::string beam, std::string target, std::string energy )
{

   TCanvas *myc1 = new TCanvas("myc1","",800,600);
   myc1->Divide(2,2);

   myc1->cd(1);
   plotModelsMC2Data( beam, target, energy, "proton", "59.1" );
   myc1->cd(2);
   plotModelsMC2Data( beam, target, energy, "neutron", "59.1" );
   myc1->cd(3);
   plotModelsMC2Data( beam, target, energy, "proton", "119.0" );
   myc1->cd(4);
   plotModelsMC2Data( beam, target, energy, "neutron", "119.0" );

   myc1->cd();

   return;

}

// -- > basic macro
//
void plotModelsMC2Data( std::string beam, std::string target, std::string energy, 
                        std::string secondary_type, std::string sec_angle )
{
            
   readKEData( beam, target, energy, secondary_type, sec_angle );
   
   double ymin = 10000.; // something big... don't know if I can use FLT_MAX
   double ymax = -1. ;
   for ( int ip=0; ip<NPtKE; ip++ )
   {
      if ( (YY[ip]+Err2[ip]) > ymax ) ymax = YY[ip]+Err2[ip];
      if ( (YY[ip]-Err2[ip]) < ymin ) ymin = YY[ip]-Err2[ip];
   }

   TGraph*  gr1 = new TGraphErrors(NPtKE,KE,YY,0,Err2);
   gr1->SetMarkerColor(4);  gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1.6);
   std::string xtitle = "Kinetic Energy of secondary " + secondary_type + " (GeV)";
   gr1->GetXaxis()->SetTitle( xtitle.c_str() );
   gr1->GetYaxis()->SetTitle("MC/Data");
   
   TGraph* gr[NModels];
   
//   int NMod = NModels;
//   if ( energy == "1.40" ) NMod -= 1;  
   
   std::string skip ="";
   if ( energy == "1.40" ) 
   {
      skip = "ftfp";
   }
   else 
   {
      skip = "binary";
   }  
   
   for ( int iv=0; iv<NModels; iv++ )
   {

      if ( Models[iv] == skip ) continue;
      
   // open G4 output root file (histo)
   std::string histofile = beam + target + Models[iv] + energy + "GeV.root";
   std::string histoname = "KE" + secondary_type + "0" + beam + target + Models[iv] + energy + "GeV";
   int counter = 5 - sec_angle.length();
   for ( int il=0; il<counter; il++ )
   {
      histoname += " "; // pad up for the fact that initially sec.angle was supposed to be char[6] - no more, no less...
   }
   histoname += sec_angle;
   
   TFile* hfile = new TFile( histofile.c_str() );
   
   std::cout << " Histo name : " << histoname << std::endl;
   TH1F* hi = (TH1F*)hfile->Get( histoname.c_str() );
   
   if ( hi == 0 || NPtKE <= 0 ) 
   {
      std::cout << " Invalid case: no exp.data or simulatiion, or both " << std::endl;
      return ;
   }
   
   //hi->Scale( (1./32.) );
      
   float MC2DataX[30];
   float MC2DataY[30];
   float DX[30], DY[30];
   int np =0;
   
   int nx = hi->GetNbinsX();
   

   for ( int k=1; k <= nx; k++ )
   {        
        double xx1 = hi->GetBinLowEdge(k);
	double xx2 = hi->GetBinWidth(k);
	for (int kk=0; kk<NPtKE; kk++ )
	{
	   if ( xx1 < KE[kk] && xx1+xx2 > KE[kk] )
	   {
	      double yy = hi->GetBinContent(k);
	      MC2DataX[np] = KE[kk];
	      DX[np] = 0.;
	      MC2DataY[np] = yy / ExpValue[kk];
	      // also error calc here !...
	      DY[np]=Err1[kk]*MC2DataY[np]/ExpValue[kk];
	      if ( (MC2DataY[np]+DY[np]) > ymax ) ymax = MC2DataY[np]+DY[np];
	      if ( (MC2DataY[np]-DY[np]) < ymin ) ymin = MC2DataY[np]-DY[np];
	      np++;
	      break;
           }
	}
   }
   
   // do NOT plot if less than 2 data points
   //
   if ( np < 2 ) return;
   
   gr[iv] = new TGraphErrors( np, MC2DataX, MC2DataY, DX, DY );
   gr[iv]->SetTitle(hi->GetTitle());
   gr[iv]->GetXaxis()->SetTitle( xtitle.c_str() );
   gr[iv]->GetYaxis()->SetTitle("MC/Data");
   gr[iv]->SetMarkerColor(ColorModel[iv]);  
   gr[iv]->SetMarkerStyle(21);
   gr[iv]->SetMarkerSize(1.0); // mrk size used to be 1.6
   if ( iv == 0 ) gr1->SetTitle(hi->GetTitle());
   
   } // end loop on models
      
   
   if ( ymin > 0. ) ymin = 0;
   if ( ymax < 3. ) ymax = 3.;
   // gr1->GetYaxis()->SetRangeUser(ymin-0.1, ymax+0.2);  
   gr1->GetYaxis()->SetRangeUser(ymin, ymax);  
   
   gr1->Draw("apl");
   
   TLegend* leg = new TLegend(0.75, 0.70, 0.99, 0.99);   
   
   for ( int iv=0; iv<NModels; iv++ )
   {
      if ( Models[iv] == skip ) continue;
      gr[iv]->Draw("lpsame");   
      leg->AddEntry( gr[iv], Models[iv].c_str(), "p" );
   } 
    
   leg->AddEntry( gr1, "exp.data", "p");     

   leg->Draw();
   leg->SetFillColor(kWhite);
   
   return;

}

void crudeMerge( std::string beam, std::string target, std::string energy, std::string model )
{

// Note: if one merges weighted/normilized histograms, in this case 
//       the resulting histogram(s) will be 32 times the statistics;
//       at the analysis stage they'd need to be scaled by 1/32 

   TFileMerger* fm = new TFileMerger();
   
   std::string output = beam + target + model + energy + "GeV.root";
   
   fm->OutputFile( output.c_str() );
   
   for ( int id=1; id<=32; id++ )
   {
      std::string filename = beam + target + model + energy + "GeV-";
      char buf[5];
      sprintf( buf, "%i", id );
      filename.append( buf ); 
      filename += ".root";    
      // std::cout << " file name = " << file << std::endl;           
      fm->AddFile( filename.c_str() );
   }
   
   fm->Merge();
   
   return;

}

void fancyMerge( std::string beam, std::string target, std::string energy, std::string model, bool doScale=false )
{
      
   std::string output = beam + target + model + energy + "GeV.root" ;
   
   targetFile = TFile::Open( output.c_str(), "RECREATE" );
   
   double scale = 1./32.;
   
   std::string input = beam + target + model + energy + "GeV-1.root";
   
   TFile* iFile1 = TFile::Open( input.c_str() );
   TIter  next( iFile1->GetListOfKeys() );
   TKey*  key = (TKey*)next();
   TH1* h = 0;
   while ( key )
   {   
         if ( !(TClass::GetClass(key->GetClassName())->InheritsFrom(TH1::Class())) ) continue;
         const char* kName = key->GetName();
         h = (TH1*)key->ReadObj();
         const char* hName = h->GetName();
         std::cout << " histoname = " << hName << std::endl;
	 TH1F* h1 = h->Clone();
	 for ( int id=2; id<=32; id++ )
	 {
	    std::string input_t = beam + target + model + energy + "GeV-" ;
            char buf[5];
            sprintf( buf, "%i", id );
            input_t.append( buf ); 
            input_t += ".root"; 
	    TFile* iFile_t = TFile::Open( input_t.c_str() );
	    TH1F* h_t = (TH1F*)iFile_t->Get( h->GetName() );
	    h1->Add( h_t );  
	    iFile_t->Close();
	 }
	 if ( doScale )
	 {
	    h1->Scale( scale );
	 }
	 targetFile->cd();
	 h1->Write();
         key = (TKey*)next();
   }
   
   targetFile->Close();
     
   return;

}




