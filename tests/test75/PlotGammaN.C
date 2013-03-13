#include "GammaData/gammaCu_300MeV.C"
#include "GammaData/gamma_400_750MeV.C"
#include "TLegend.h"
#include <string>


void PlotGammaN()
{

   PlotGamma300Cu();
      
   return;

}

void PlotGamma300Cu() 
{
   gStyle->SetOptStat(0);
   
   TCanvas* myc = new TCanvas( "myc", "", 1200, 400 );   
   myc->Divide(3,1);
   
   myc->cd(1);
   gPad->SetLogy();
   PlotGammaNuclear( "300", "Cu", "proton", "45" );

   myc->cd(2);
   gPad->SetLogy();
   PlotGammaNuclear( "300", "Cu", "proton", "90" );

   myc->cd(3);
   gPad->SetLogy();
   PlotGammaNuclear( "300", "Cu", "proton", "135" );

   return;
   
}

void PlotGamma668Cu()
{

   gStyle->SetOptStat(0);

   TCanvas* myc = new TCanvas( "myc", "", 800, 800 );   
   myc->Divide(2,2);
   
   myc->cd(1);
   gPad->SetLogy();
   PlotGammaNuclear( "668", "Cu", "pi-", "28" );

   myc->cd(2);
   gPad->SetLogy();
   PlotGammaNuclear( "668", "Cu", "pi+", "28" );   

   myc->cd(3);
   gPad->SetLogy();
   PlotGammaNuclear( "668", "Cu", "pi-", "44" );

   myc->cd(4);
   gPad->SetLogy();
   PlotGammaNuclear( "668", "Cu", "pi+", "44" );

   return;

}

void PlotGamma668Pb()
{

   gStyle->SetOptStat(0);

   TCanvas* myc = new TCanvas( "myc", "", 800, 400 );   
   myc->Divide(2,1);
   
   myc->cd(1);
   gPad->SetLogy();
   PlotGammaNuclear( "668", "Pb", "pi-", "44" );

   myc->cd(2);
   gPad->SetLogy();
   PlotGammaNuclear( "668", "Pb", "pi+", "44" );   

   return;

}

void PlotGammaNuclear( std::string energy, std::string target, std::string secondary, std::string angle )
{
   
// FIXME !!!
// Will need an update In case of several models
//

   std::string fName = "gamma"+ energy + "MeV" + target;
   
   std::string fName1 = fName + "Bertini.root";
   TFile* ifile1 = new TFile( fName1.c_str() );
   
// In case of CHIPS...
// -->    std::string fName2 = fName + "CHIPS.root";
// -->   TFile* ifile2 = new TFile( fName2.c_str() );
   
   TLegend* leg1 = new TLegend(0.6, 0.70, 0.9, 0.9);
   // TLegend* leg1 = new TLegend(0.2, 0.15, 0.5, 0.35);
   
   std::string hName = "";
   
   if ( secondary == "proton" )
   {
      hName = "p";
   } 
   else if ( secondary == "pi-" )
   {
      hName = "pim";
   }
   else if (secondary == "pi+" )
   {
      hName = "pip";
   }
   hName += angle ;
   hName +="deg";
      
   TH1F* histo1 = (TH1F*)ifile1->Get( hName.c_str() );   
   histo1->SetLineColor(kRed);
   histo1->SetLineWidth(3);
   TString hTit1 = "gamma + " + target;
   TString hTit2 = " X + " + secondary + " (" + angle + "deg)";
   histo1->SetTitle( hTit1 + " #rightarrow " + hTit2 );
   if ( secondary == "proton" )
   {
      histo1->GetXaxis()->SetTitle("kinetic energy of secondary proton (MeV)");
      histo1->GetYaxis()->SetTitle("d^{2}#sigma / dE d#Omega  [#mub MeV^{-1} sr^{-1}]");
   }
   else if ( secondary == "pi-" || secondary == "pi+" )
   {
      histo1->GetYaxis()->SetTitle("d^{2}#sigma / dp d#Omega  [#mub MeV^{-1} sr^{-1}]");   
   }
   histo1->GetYaxis()->SetTitleOffset(1.5);
   histo1->Draw();
   leg1->AddEntry( histo1, "Bertini", "L" );
   OverlayData( hName, target, leg1 );
   leg1->Draw();
   leg1->SetFillColor(kWhite);
   
   return;

}

TGraph* GetOverlayCu( const TString& hname ) 
{

    if (hname == "p45deg")  return gCup45();
    if (hname == "p90deg")  return gCup90();
    if (hname == "p135deg") return gCup135();
    
    if ( hname == "pim28deg" ) return gGam668Cu_PiMin_28();
    if ( hname == "pip28deg" ) return gGam668Cu_PiPls_28();
    if ( hname == "pim44deg" ) return gGam668Cu_PiMin_44();
    if ( hname == "pip44deg" ) return gGam668Cu_PiPls_44();

    return 0;
    
}

TGraph* GetOverlayPb( const TString& hname ) 
{

   if ( hname == "pim44deg" ) return gGam668Pb_PiMin_44();
   if ( hname == "pip44deg" ) return gGam668Pb_PiPls_44();

   return 0;

}

void OverlayData(const TString& hname, const TString& target, TLegend* leg=0 ) 
{

  TGraph* data = 0;
  if ( target == "Cu" )
  {
     data = GetOverlayCu(hname);
  }
  else if ( target == "Pb" )
  {
     data = GetOverlayPb(hname);
  }
  else if ( target == "C" ) // provision for future addition of gamma + C data
  {
     // data = GetOverlayC(hname);
  }
  
  if ( data )
  {
     data->SetMarkerStyle(22);
     data->SetMarkerColor(kBlue);
     data->SetMarkerSize(1.6);
     if ( leg ) leg->AddEntry( data, "exp.data", "p" );
     data->Draw("P same");
  }
  
  return;

}
