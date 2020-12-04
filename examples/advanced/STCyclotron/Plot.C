#include "TF1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH1.h"
#include "TMath.h"
#include <string.h>
#include "TGraph.h"
#include <map>

using namespace std;

void norm_th1_per_bin_width_per_primaries(TH1D* histo, int total_primaries)
{
  int xbins;
  double value,xbinwidth;
  
  //Normalize histogram per bin width and per incoming particle.
  xbins = histo->GetXaxis()->GetNbins();
  for(int i=1; i<xbins;i++)
    {
      xbinwidth = histo->GetBinWidth(i);
      value = histo->GetBinContent(i);
      value = value/xbinwidth/(total_primaries*1.);
      histo->SetBinContent(i,value);

      //Setting error bin content
      value = histo->GetBinError(i);
      value = value/xbinwidth/(total_primaries*1.);
      histo->SetBinError(i,value);
    }
}

void norm_th2_per_bin_width_per_primaries(TH2D* histo, int total_primaries)
{
  int xbins,ybins;
  double value,xbinwidth,ybinwidth;
  
  //Normalize histogram per bin width and per incoming particle.
  xbins = histo->GetXaxis()->GetNbins();
  ybins = histo->GetYaxis()->GetNbins();
  for(int i=1; i<xbins;i++)
    {
      xbinwidth = histo->GetXaxis()->GetBinWidth(i);

      for(int j=1; j<ybins;j++)
	{
	  ybinwidth = histo->GetYaxis()->GetBinWidth(j);
      
	  value = histo->GetBinContent(i,j);
	  value = value/xbinwidth/ybinwidth/(total_primaries*1.);
	  histo->SetBinContent(i,j,value);

	  //Setting error bin content
	  value = histo->GetBinError(i,j);
	  value = value/xbinwidth/ybinwidth/(total_primaries*1.);
	  histo->SetBinError(i,j,value);
	}
    }
}


void Plot(){

  
  //PARAMETERS
  double tMin = 0;
  double tMax = 30.;              //in hour(s)
  double halfLifeLimit = 1./60.;   //in hour(s)

  double tIrradiation;       //in hour(s)
  double beamCurrent;        //in µA
  int total_primaries;


  
  //VARIABLES
  string endLine;
  
  //Getting the parameters from the G4 output file.
  ifstream G4output;
  G4output.open("Output_General.txt");
  for(int i=0;i<6;i++)getline(G4output,endLine);
  G4output >> beamCurrent; getline(G4output,endLine);
  G4output >> tIrradiation; getline(G4output,endLine);
  for(int i=0;i<6;i++)getline(G4output,endLine);
  G4output >> total_primaries;
  G4output.close();
  
  beamCurrent*=1e6; //<--- convert from A to µA.

  /*
  cout << "Irradiation time = " << tIrradiation << " h." << endl;
  cout << "Beam current = " << beamCurrent << " µA." << endl;
  cout << "Total primaries = " << total_primaries << endl;
  getchar();*/

  
  system("rm -r Results");
  system("mkdir Results");
  system("mkdir Results/IsotopesProduction");
  system("mkdir Results/ParticlesEnergySpectra");
  system("mkdir Results/ParticlesEnergySpectra/Beam");
  system("mkdir Results/ParticlesEnergySpectra/Decay");
  system("mkdir Results/BeamData");

  
  ofstream results;
  results.open("Results.txt"); 
  results << "Parameters: " << endl;
  results << "Time of irradiation: " << tIrradiation << " hour(s)." << endl;
  results << "Beam current: " << beamCurrent << " µA." << endl;
  results << "Total number of simulated primaries: " << total_primaries << endl;
  results << "Please check they are the same as in the simulation. Otherwise change it by modifying the Plot.C file." << endl;

  //Opening root file.

  stringstream name_root_file;
  name_root_file << "./SolidTargetCyclotron.root";
  TFile *f = new TFile(name_root_file.str().c_str(),"open");

  //---------------------------------------------------------------//
  //       Energy Profile of the beam before/after the target      //
  //---------------------------------------------------------------//

  TCanvas *BeamEnergyTarget = new TCanvas("Beam energy profile before the target", "BeamTargetProfile");

  TH1D *energyProfileBeamTarget = (TH1D*)f->Get("H10;1");
  energyProfileBeamTarget->GetXaxis()->SetTitle("Energy (MeV)");
  energyProfileBeamTarget->GetYaxis()->SetTitle("P(E) (MeV^{-1}.particle^{-1})");
  energyProfileBeamTarget->SetTitle("Primary particles energy when reaching the target, per primary particle");
  
  energyProfileBeamTarget->GetXaxis()->SetMaxDigits(2);
  energyProfileBeamTarget->GetYaxis()->SetMaxDigits(3);

  //Normalize histogram per bin width and per incoming particle.
  norm_th1_per_bin_width_per_primaries(energyProfileBeamTarget,total_primaries);

  energyProfileBeamTarget->Draw("H");
  BeamEnergyTarget->Print("./Results/BeamData/BeamEnergyInTarget.pdf");

  

  TCanvas *BeamEnergyOutTarget = new TCanvas("Beam energy profile after the target", "BeamTargetOutProfile");

  TH1D *energyProfileBeamOutTarget = (TH1D*)f->Get("H12;1");
  energyProfileBeamOutTarget->GetXaxis()->SetTitle("Energy (MeV)");
  energyProfileBeamOutTarget->GetYaxis()->SetTitle("P(E) (MeV^{-1}.particle^{-1})");
  energyProfileBeamOutTarget->SetTitle("Primary particles energy when going out of the target, per primary particle");

  energyProfileBeamOutTarget->GetXaxis()->SetMaxDigits(2);
  energyProfileBeamOutTarget->GetYaxis()->SetMaxDigits(3);
  
  //Normalize histogram per bin width and per incoming particle.
  norm_th1_per_bin_width_per_primaries(energyProfileBeamOutTarget,total_primaries);
 
  energyProfileBeamOutTarget->Draw("H");
  BeamEnergyOutTarget->Print("./Results/BeamData/BeamEnergyOutTarget.pdf");



  //---------------------------------------------------------------//
  //        Energy Profile of the beam before/after the foil       //
  //---------------------------------------------------------------//

  TCanvas *BeamEnergyFoil = new TCanvas("Beam energy profile before the foil", "BeamFoilProfile");

  TH1D *energyProfileBeamFoil = (TH1D*)f->Get("H11;1");
  energyProfileBeamFoil->GetXaxis()->SetTitle("Energy (MeV)");
  energyProfileBeamFoil->GetYaxis()->SetTitle("P(E) (MeV^{-1}.particle^{-1})");
  energyProfileBeamFoil->SetTitle("Primary particles energy when reaching the foil, per primary particle");
  
  energyProfileBeamFoil->GetXaxis()->SetMaxDigits(2);
  energyProfileBeamFoil->GetYaxis()->SetMaxDigits(3);
  
  //Normalize histogram per bin width and per incoming particle.
  norm_th1_per_bin_width_per_primaries(energyProfileBeamFoil,total_primaries);
  
  energyProfileBeamFoil->Draw("H");
  BeamEnergyFoil->Print("./Results/BeamData/BeamEnergyInFoil.pdf");

  

  TCanvas *BeamEnergyOutFoil = new TCanvas("Beam energy profile after the foil", "BeamFoilOutProfile");
  TH1D *energyProfileBeamOutFoil = (TH1D*)f->Get("H13;1");
  energyProfileBeamOutFoil->GetXaxis()->SetTitle("Energy (MeV)");
  energyProfileBeamOutFoil->GetYaxis()->SetTitle("P(E) (MeV^{-1}.particle^{-1})");
  energyProfileBeamOutFoil->SetTitle("Primary particles energy when going out of the foil, per primary particle");
  
  energyProfileBeamOutFoil->GetXaxis()->SetMaxDigits(2);
  energyProfileBeamOutFoil->GetYaxis()->SetMaxDigits(3);
    
  //Normalize histogram per bin width and per incoming particle.
  norm_th1_per_bin_width_per_primaries(energyProfileBeamOutFoil,total_primaries);

  energyProfileBeamOutFoil->Draw("H");
  BeamEnergyOutFoil->Print("./Results/BeamData/BeamEnergyOutFoil.pdf");



  //---------------------------------------------------------------//
  //            Depth of isotope creation in the target            //
  //---------------------------------------------------------------//

  TCanvas *depthCreation = new TCanvas("DepthCreation", "Depth of isotope creation in the target per primary particle.");

  TH1D *hDepthCreation = (TH1D*) f->Get("H14;1");
  hDepthCreation->SetTitle("Depth of isotope creation in the target per primary particle.");
  hDepthCreation->GetXaxis()->SetTitle("Depth (mm)");
  hDepthCreation->GetYaxis()->SetTitle("N isotopes (mm^{-1}.particle^{-1})");

  hDepthCreation->GetYaxis()->SetMaxDigits(3);
   
  //Normalize histogram per bin width and per incoming particle.
  norm_th1_per_bin_width_per_primaries(hDepthCreation,total_primaries);

  hDepthCreation->SetMarkerStyle(4);
  hDepthCreation->SetMarkerSize(1);
  hDepthCreation->Draw("l");

  depthCreation->Print("./Results/IsotopesProduction/DepthCreation.pdf"); 


  //---------------------------------------------------------------//
  //                        Energy spectrum                        //
  //---------------------------------------------------------------//

  //----------------->> PARTICLES EMITTED DUE TO BEAM INTERACTIONS WITH THE TARGET

  //Positrons//
  TCanvas *PositronSpectrumBeam = new TCanvas("PositronSpectrumBeam", "Spectrum of the positrons created by the beam in the target");
  TH1D *hPositronSpectrumBeam = (TH1D*) f->Get("H15;1");
  if(hPositronSpectrumBeam->GetEntries()!=0)
    {
      hPositronSpectrumBeam->GetXaxis()->SetTitle("Energy (MeV)");
      hPositronSpectrumBeam->GetYaxis()->SetTitle("N positrons (MeV^{-1}.particle^{-1})");
      //Normalize histogram per bin width and per incoming particle.
      norm_th1_per_bin_width_per_primaries(hPositronSpectrumBeam,total_primaries);
      hPositronSpectrumBeam->GetYaxis()->SetMaxDigits(3);
      hPositronSpectrumBeam->GetYaxis()->SetTitleOffset(1.2);
      hPositronSpectrumBeam->Draw("H");
      PositronSpectrumBeam->SetLogy();
      PositronSpectrumBeam->Print("./Results/ParticlesEnergySpectra/Beam/PositronSpectrumBeam.pdf"); 
    }

 //Electrons//
  TCanvas *ElectronSpectrumBeam = new TCanvas("ElectronSpectrumBeam", "Spectrum of the electrons created by the beam in the target");
  TH1D *hElectronSpectrumBeam = (TH1D*) f->Get("H16;1");
  if(hElectronSpectrumBeam->GetEntries() !=0)
    {
      hElectronSpectrumBeam->GetXaxis()->SetTitle("Energy (MeV)");
      hElectronSpectrumBeam->GetYaxis()->SetTitle("N electrons (MeV^{-1}.particle^{-1})");
      hElectronSpectrumBeam->GetYaxis()->SetTitleOffset(1.2);
      //Normalize histogram per bin width and per incoming particle.
      norm_th1_per_bin_width_per_primaries(hElectronSpectrumBeam,total_primaries);
      hElectronSpectrumBeam->Draw("H");
      ElectronSpectrumBeam->SetLogy();
      ElectronSpectrumBeam->Print("./Results/ParticlesEnergySpectra/Beam/ElectronSpectrumBeam.pdf"); 
    }
  
  //Gammas//  
  TCanvas *GammaSpectrumBeam = new TCanvas("GammaSpectrumBeam", "Spectrum of the gammas created by the beam in the target");
  TH1D *hGammaSpectrumBeam = (TH1D*) f->Get("H17;1");
  if(hGammaSpectrumBeam->GetEntries() !=0)
    {
      hGammaSpectrumBeam->GetXaxis()->SetTitle("Energy (MeV)");
      hGammaSpectrumBeam->GetYaxis()->SetTitle("N Gammas (MeV^{-1}.particle^{-1})");
      hGammaSpectrumBeam->GetYaxis()->SetTitleOffset(1.2);
      hGammaSpectrumBeam->Draw("H");
      //Normalize histogram per bin width and per incoming particle.
      norm_th1_per_bin_width_per_primaries(hGammaSpectrumBeam,total_primaries);
      GammaSpectrumBeam->SetLogy();
      GammaSpectrumBeam->Print("./Results/ParticlesEnergySpectra/Beam/GammaSpectrumBeam.pdf"); 
    }

  //Neutrons//
  TCanvas *NeutronSpectrumBeam = new TCanvas("NeutronSpectrumBeam", "Spectrum of the neutrons created by the beam in the target");
  TH1D *hNeutronSpectrumBeam = (TH1D*) f->Get("H18;1");
  if(hNeutronSpectrumBeam->GetEntries() !=0)
    {
      hNeutronSpectrumBeam->GetXaxis()->SetTitle("Energy (MeV)");
      hNeutronSpectrumBeam->GetYaxis()->SetTitle("N neutrons (MeV^{-1}.particle^{-1})");
      hNeutronSpectrumBeam->GetYaxis()->SetTitleOffset(1.2);
      hNeutronSpectrumBeam->Draw("H");
      //Normalize histogram per bin width and per incoming particle.
      norm_th1_per_bin_width_per_primaries(hNeutronSpectrumBeam,total_primaries);
      NeutronSpectrumBeam->SetLogy();
      NeutronSpectrumBeam->Print("./Results/ParticlesEnergySpectra/Beam/NeutronSpectrumBeam.pdf"); 
    }

  
  //----------------->> PARTICLES EMITTED DUE TO ISOTOPE DECAY

  //Positrons//
  TCanvas *PositronSpectrumDecay = new TCanvas("PositronSpectrumDecay", "Spectrum of the positrons created by the decays in the target");
  TH1D *hPositronSpectrumDecay = (TH1D*) f->Get("H19;1");
  if(hPositronSpectrumDecay->GetEntries() !=0)
    {
      hPositronSpectrumDecay->GetXaxis()->SetTitle("Energy (MeV)");
      hPositronSpectrumDecay->GetYaxis()->SetTitle("N positrons (MeV^{-1}.particle^{-1})");
      hPositronSpectrumDecay->GetYaxis()->SetTitleOffset(1.2);
      //Normalize histogram per bin width and per incoming particle.
      norm_th1_per_bin_width_per_primaries(hPositronSpectrumDecay,total_primaries);
      hPositronSpectrumDecay->Draw("H");
      PositronSpectrumDecay->SetLogy();
      PositronSpectrumDecay->Print("./Results/ParticlesEnergySpectra/Decay/PositronSpectrumDecay.pdf"); 
    }
  
  //Electrons//
  TCanvas *ElectronSpectrumDecay = new TCanvas("ElectronSpectrumDecay", "Spectrum of the electrons created by the decays in the target");
  TH1D *hElectronSpectrumDecay = (TH1D*) f->Get("H110;1");
  if(hElectronSpectrumDecay->GetEntries() !=0)
    {
      hElectronSpectrumDecay->GetXaxis()->SetTitle("Energy (MeV)");
      hElectronSpectrumDecay->GetYaxis()->SetTitle("N electrons (MeV^{-1}.particle^{-1})");
      hElectronSpectrumDecay->GetYaxis()->SetTitleOffset(1.2);
      //Normalize histogram per bin width and per incoming particle.
      norm_th1_per_bin_width_per_primaries(hElectronSpectrumDecay,total_primaries);  
      hElectronSpectrumDecay->Draw("H");
      ElectronSpectrumDecay->SetLogy();
      ElectronSpectrumDecay->Print("./Results/ParticlesEnergySpectra/Decay/ElectronSpectrumDecay.pdf"); 
    }
  
  //Gammas//
  TCanvas *GammaSpectrumDecay = new TCanvas("GammaSpectrumDecay", "Spectrum of the gammas created by the decays in the target");
  TH1D *hGammaSpectrumDecay = (TH1D*) f->Get("H111;1");
  if(hGammaSpectrumDecay->GetEntries() !=0)
    {
      hGammaSpectrumDecay->GetXaxis()->SetTitle("Energy (MeV)");
      hGammaSpectrumDecay->GetYaxis()->SetTitle("N Gammas (MeV^{-1}.particle^{-1})");
      hGammaSpectrumDecay->GetYaxis()->SetTitleOffset(1.2);
      //Normalize histogram per bin width and per incoming particle.
      norm_th1_per_bin_width_per_primaries(hGammaSpectrumDecay,total_primaries);  
      hGammaSpectrumDecay->Draw("H");
      GammaSpectrumDecay->SetLogy();
      GammaSpectrumDecay->Print("./Results/ParticlesEnergySpectra/Decay/GammaSpectrumDecay.pdf"); 
    }

  //Neutrons//
   TCanvas *NeutronSpectrumDecay = new TCanvas("NeutronSpectrumDecay", "Spectrum of the neutrons created by the decays in the target");
   TH1D *hNeutronSpectrumDecay = (TH1D*) f->Get("H112;1");
   if(hNeutronSpectrumDecay->GetEntries() !=0)
     {
       hNeutronSpectrumDecay->GetXaxis()->SetTitle("Energy (MeV)");
       hNeutronSpectrumDecay->GetYaxis()->SetTitle("N Gammas (MeV^{-1}.particle^{-1})");
       hNeutronSpectrumDecay->GetYaxis()->SetTitleOffset(1.2);
       //Normalize histogram per bin width and per incoming particle.
       norm_th1_per_bin_width_per_primaries(hNeutronSpectrumDecay,total_primaries);  
       hNeutronSpectrumDecay->Draw("H");
       NeutronSpectrumDecay->SetLogy();
       NeutronSpectrumDecay->Print("./Results/ParticlesEnergySpectra/Decay/NeutronSpectrumBeam.pdf");
     }

  
  //Nu_e//
  TCanvas *NuESpectrumDecay = new TCanvas("NuESpectrumDecay", "Spectrum of the Nu_e created by the decays in the target");
  TH1D *hNuESpectrumDecay = (TH1D*) f->Get("H113;1");
  if(hNuESpectrumDecay->GetEntries() !=0)
    {
      hNuESpectrumDecay->GetXaxis()->SetTitle("Energy (MeV)");
      hNuESpectrumDecay->GetYaxis()->SetTitle("N Nu_{e} (MeV^{-1}.particle^{-1})");
      hNuESpectrumDecay->GetYaxis()->SetTitleOffset(1.2);
      //Normalize histogram per bin width and per incoming particle.
      norm_th1_per_bin_width_per_primaries(hNuESpectrumDecay,total_primaries);  
      hNuESpectrumDecay->Draw("H");
      NuESpectrumDecay->SetLogy();
      NuESpectrumDecay->Print("./Results/ParticlesEnergySpectra/Decay/NuESpectrumDecay.pdf");
    }
  
  //AntiNu_e//
  TCanvas *AntiNuESpectrumDecay = new TCanvas("AntiNuESpectrumDecay", "Spectrum of the Anti_Nu_e created by the decays in the target");  
  TH1D *hAntiNuESpectrumDecay = (TH1D*) f->Get("H114;1");
  if(hAntiNuESpectrumDecay->GetEntries() !=0)
    {
      hAntiNuESpectrumDecay->GetXaxis()->SetTitle("Energy (MeV)");
      hAntiNuESpectrumDecay->GetYaxis()->SetTitle("N AntiNu_{e} (MeV^{-1}.particle^{-1})");
      hAntiNuESpectrumDecay->GetYaxis()->SetTitleOffset(1.2);
      //Normalize histogram per bin width and per incoming particle.
      norm_th1_per_bin_width_per_primaries(hAntiNuESpectrumDecay,total_primaries);  
      hAntiNuESpectrumDecay->Draw("H");
      AntiNuESpectrumDecay->SetLogy();
      AntiNuESpectrumDecay->Print("./Results/ParticlesEnergySpectra/Decay/AntiNuESpectrumDecay.pdf");
    }

 

  
  /////////////////
  //2D histograms//
  /////////////////
  
  TH2D *hBeamIntensityTarget = (TH2D*) f->Get("H20;1");
  if(hBeamIntensityTarget->GetEntries()!=0)
    {
      TCanvas *BeamIntensityTarget = new TCanvas("BeamIntensityTarget", "Beam intensity (particle^{-1}.mm^{-2}) before hiting the target");
      hBeamIntensityTarget->GetXaxis()->SetTitle("X axis (mm)");
      hBeamIntensityTarget->GetYaxis()->SetTitle("Y axis (mm)");
      hBeamIntensityTarget->SetTitle("Beam intensity (particle^{-1}.mm^{-2}) before hiting the target");
      //Normalizing
      norm_th2_per_bin_width_per_primaries(hBeamIntensityTarget, total_primaries);
      hBeamIntensityTarget->GetXaxis()->SetMaxDigits(3);
      hBeamIntensityTarget->GetYaxis()->SetMaxDigits(3);
      hBeamIntensityTarget->GetZaxis()->SetMaxDigits(3);
      hBeamIntensityTarget->Draw("colz");
  
      gStyle->SetOptStat(0); 
      BeamIntensityTarget->Update();
      BeamIntensityTarget->Print("./Results/BeamData/BeamIntensityTarget.pdf");
      BeamIntensityTarget->Print("./Results/BeamData/BeamIntensityTarget.jpg");
    }
  
  TH2D *hBeamIntensityFoil = (TH2D*) f->Get("H21;1");
  if(hBeamIntensityFoil->GetEntries()!=0)
    {
      TCanvas *BeamIntensityFoil = new TCanvas("BeamIntensityFoil", "Beam intensity before hiting the foil");
      hBeamIntensityFoil->GetXaxis()->SetTitle("X axis (mm)");
      hBeamIntensityFoil->GetYaxis()->SetTitle("Y axis (mm)");
      hBeamIntensityFoil->SetTitle("Beam intensity (particle^{-1}.mm^{-2}) before hiting the foil");
      //Normalizing
      norm_th2_per_bin_width_per_primaries(hBeamIntensityFoil, total_primaries);
      hBeamIntensityFoil->GetXaxis()->SetMaxDigits(3);
      hBeamIntensityFoil->GetYaxis()->SetMaxDigits(3);
      hBeamIntensityFoil->GetZaxis()->SetMaxDigits(3);
      hBeamIntensityFoil->Draw("colz");
  
      BeamIntensityFoil->Print("./Results/BeamData/BeamIntensityFoil.pdf");
    }

  TH2D *hBeamIntensityOutTarget = (TH2D*) f->Get("H24;1");  
  if(hBeamIntensityOutTarget->GetEntries()!=0)
    {
      TCanvas *BeamIntensityOutTarget = new TCanvas("BeamIntensityOutTarget", "Beam intensity going out from the target");
      
      hBeamIntensityOutTarget->GetXaxis()->SetTitle("X axis (mm)");
      hBeamIntensityOutTarget->GetYaxis()->SetTitle("Y axis (mm)");
      hBeamIntensityOutTarget->SetTitle("Beam intensity (particle^{-1}.mm^{-2}) after hiting the target");
      hBeamIntensityOutTarget->Draw("colz");
     //Normalizing
      norm_th2_per_bin_width_per_primaries(hBeamIntensityOutTarget, total_primaries);
      hBeamIntensityOutTarget->GetXaxis()->SetMaxDigits(3);
      hBeamIntensityOutTarget->GetYaxis()->SetMaxDigits(3);
      hBeamIntensityOutTarget->GetZaxis()->SetMaxDigits(3);
       
      BeamIntensityOutTarget->Print("./Results/BeamData/BeamIntensityOutTarget.pdf");
    }

    
  
  TH2D *hRadioisotopeProduction = (TH2D*) f->Get("H22;1");
  if(hRadioisotopeProduction->GetEntries()!=0)
    {
      TCanvas *RadioisotopeProduction = new TCanvas("RadioisotopeProduction", "Radioisotope production");
      hRadioisotopeProduction->GetXaxis()->SetTitle("Z");
      hRadioisotopeProduction->GetXaxis()->SetTitleOffset(1.2);
      hRadioisotopeProduction->GetYaxis()->SetTitle("A");
      hRadioisotopeProduction->GetYaxis()->SetTitleOffset(1.3);
      hRadioisotopeProduction->GetZaxis()->SetTitle("N radioisotopes (particle^{-1}.mm^{-2})");
      hRadioisotopeProduction->GetZaxis()->SetTitleOffset(1.3);
      hRadioisotopeProduction->SetTitle("Number of radioisotopes produced in the target (particle^{-1}.mm^{-2})");
      //Normalizing
      norm_th2_per_bin_width_per_primaries(hRadioisotopeProduction, total_primaries);  
      hRadioisotopeProduction->GetXaxis()->SetMaxDigits(3);
      hRadioisotopeProduction->GetYaxis()->SetMaxDigits(3);
      hRadioisotopeProduction->GetZaxis()->SetMaxDigits(3);
      hRadioisotopeProduction->Draw("lego2");

      RadioisotopeProduction->SetLogz();
      RadioisotopeProduction->Print("./Results/IsotopesProduction/RadioisotopeProduction.pdf");
      RadioisotopeProduction->Print("./Results/IsotopesProduction/RadioisotopeProduction.jpg");
    }
  
 

  
  TH2D *hEnergyDepth = (TH2D*) f->Get("H23;1");
  if(hEnergyDepth->GetEntries()!=0)
    {
  
      TCanvas *EnergyDepth = new TCanvas("EnergyDepth", "Energy of the proton according to their depth in the target");

      hEnergyDepth->GetXaxis()->SetTitle("Depth (mm)");
      hEnergyDepth->GetYaxis()->SetTitle("Energy (MeV)");
      hEnergyDepth->SetTitle("Energy of the proton according to their depth in the target (particle^{-1}.mm^{-1}.MeV^{-1})");
      norm_th2_per_bin_width_per_primaries(hEnergyDepth, total_primaries);  
      hEnergyDepth->GetXaxis()->SetMaxDigits(3);
      hEnergyDepth->GetYaxis()->SetMaxDigits(3);
      hEnergyDepth->GetZaxis()->SetMaxDigits(3);
      hEnergyDepth->Draw("colz");
      
      EnergyDepth->SetLogz();
      EnergyDepth->Print("./Results/BeamData/EnergyDepth.pdf");
    }


    /////////////////////////////////////////
  //Stable isotope production by the beam//
  /////////////////////////////////////////

  
  //---------------------------------------------------------------//
  //                           Activity                            //
  //---------------------------------------------------------------//

  /*
  TCanvas *ActivityPrimary = new TCanvas("ActivityPerParentIsotope", "Activity in mCi per parent isotope, and total activity");
  TH1D *hActivityPrimary = (TH1D*) f->Get("H12;1");
  hActivityPrimary->GetXaxis()->Set(hActivityPrimary->GetEntries(),0.5,hActivityPrimary->GetEntries()+0.5);
  hActivityPrimary->GetXaxis()->SetTitle("Isotope");
  hActivityPrimary->GetYaxis()->SetTitle("Activity (mCi)");
  hActivityPrimary->Draw();
 
  ActivityPrimary->SetLogy();
  ActivityPrimary->Print("./Results/ActivityPerParentIsotope.pdf");

  TCanvas *ActivityDaughter = new TCanvas("ActivityPerDaughterIsotope", "Activity in mCi per daughter isotope, and total activity");

  hActivityDaughter = (TH1F*) h1DH118->Clone();
  hActivityDaughter->GetXaxis()->Set(hActivityDaughter.GetEntries(),0.5,hActivityDaughter.GetEntries()+0.5);
  hActivityDaughter->GetXaxis()->SetTitle("Isotope");
  hActivityDaughter->GetYaxis()->SetTitle("Activity (mCi)");
  hActivityDaughter->Draw();
 
  ActivityDaughter->SetLogy();
  ActivityDaughter->Print("./Results/ActivityPerDaughterIsotope.pdf");*/

  

  /*
  TCanvas *StableIsotopes = new TCanvas("StableIsotopes", "Production of stable isotopes in the target");
  
  hStableIsotopes = (TH1F*) h1DH117->Clone();
  hStableIsotopes->GetXaxis()->Set(hStableIsotopes->GetEntries(),0.5,hStableIsotopes->GetEntries()+0.5);
  hStableIsotopes->GetXaxis()->SetTitle("Stable Isotope");
  hStableIsotopes->GetYaxis()->SetTitle("Yield");
  hStableIsotopes->Draw();
  
  StableIsotopes->SetLogy();
  StableIsotopes->Print("./Results/IsotopesProduction/StableIsotopes.pdf");
  */
  /*
  /////////////////////
  //Yield per isotope//
  /////////////////////

  TCanvas *YieldParent = new TCanvas("YieldPerParentIsotope", "Yield per parent isotope");

  hYieldParent = (TH1F*) h1DH14->Clone();
  hYieldParent->GetXaxis()->Set(hYieldParent.GetEntries(),0.5,hYieldParent.GetEntries()+0.5);
  hYieldParent->GetXaxis()->SetTitle("Isotope");
  hYieldParent->GetYaxis()->SetTitle("Yield");
  hYieldParent->Draw();
 
  YieldParent->SetLogy();
  YieldParent->Print("./Results/YieldPerParentIsotope.pdf");

  TCanvas *YieldDaughter = new TCanvas("YieldPerDaughterIsotope", "Yield per daughter isotope");

  hYieldDaughter = (TH1F*) h1DH119->Clone();
  hYieldDaughter->GetXaxis()->Set(hYieldDaughter.GetEntries(),0.5,hYieldDaughter.GetEntries()+0.5);
  hYieldDaughter->GetXaxis()->SetTitle("Isotope");
  hYieldDaughter->GetYaxis()->SetTitle("Yield");
  hYieldDaughter->Draw();
 
  YieldDaughter->SetLogy();
  YieldDaughter->Print("./Results/YieldPerDaughterIsotope.pdf");

 
  TCanvas *ProductionPerSecParent = new TCanvas("ProductionPerSecParent", "Production per second per parent");

  hProdPerSec = (TH1F*) h1DH123->Clone();
  hProdPerSec->GetXaxis()->Set(hProdPerSec.GetEntries(),0.5,hProdPerSec.GetEntries()+0.5);
  hProdPerSec->GetXaxis()->SetTitle("Isotope");
  hProdPerSec->GetYaxis()->SetTitle("Production of isotope per second");
  hProdPerSec->Draw();
 
  ProductionPerSecParent->SetLogy();
  ProductionPerSecParent->Print("./Results/ProductionPerSec.pdf");

  
  TCanvas *ProductionPerSecDaughter = new TCanvas("ProductionPerSecDaughter", "Production per second per of the parent of the isotopes");
   
  hProdPerSecDaughter = (TH1F*) h1DH124->Clone();
  hProdPerSecDaughter->GetXaxis()->Set(hProdPerSecDaughter.GetEntries(),0.5,hProdPerSecDaughter.GetEntries()+0.5);
  hProdPerSecDaughter->GetXaxis()->SetTitle("Isotope");
  hProdPerSecDaughter->GetYaxis()->SetTitle("Production of the parent of the isotope per second");
  hProdPerSecDaughter->Draw();
 
  ProductionPerSecDaughter->SetLogy();
  ProductionPerSecDaughter->Print("./Results/ProductionPerSecDaughter.pdf");

  

  */  
  //////////////////
  //Decay constant//
  //////////////////
  /*
  TH1D *hYieldParent = (TH1D*) f->Get("H14;1");
  hYieldParent->GetXaxis()->Set(hYieldParent->GetEntries(),0.5,hYieldParent->GetEntries()+0.5);

  TH1D *hYieldDaughter = (TH1D*) f->Get("H119");
  hYieldDaughter->GetXaxis()->Set(hYieldDaughter->GetEntries(),0.5,hYieldDaughter->GetEntries()+0.5);

  TH1D *hProdPerSec = (TH1D*) f->Get("H123");
  hProdPerSec->GetXaxis()->Set(hProdPerSec->GetEntries(),0.5,hProdPerSec->GetEntries()+0.5);

  TH1D *hProdPerSecDaughter = (TH1D*) f->Get("H124");
  hProdPerSecDaughter->GetXaxis()->Set(hProdPerSecDaughter->GetEntries(),0.5,hProdPerSecDaughter->GetEntries()+0.5);
 
  TH1D *hConstantParent = (TH1D*) f->Get("H120;1");
  hConstantParent->GetXaxis()->Set(hConstantParent->GetEntries(),0.5,hConstantParent->GetEntries()+0.5);
 
  TH1D *hConstantDaughter = (TH1D*) f->Get("H121;1");
  hConstantDaughter->GetXaxis()->Set(hConstantDaughter->GetEntries(),0.5,hConstantDaughter->GetEntries()+0.5);

  TH1D *hConstantDaughterParent = (TH1D*) f->Get("H122;1");
  hConstantDaughterParent->GetXaxis()->Set(hConstantDaughterParent->GetEntries(),0.5,hConstantDaughterParent->GetEntries()+0.5);
  */


  /////////////////////////////
  //Plots of theorical curves//
  /////////////////////////////

  //----------------------------------------------------------------------------
  //                 CASE OF PARENT ISOTOPES
  //----------------------------------------------------------------------------

  vector<string> name;            //<---- name of the isotope
  vector<double> hDecayConstant;   //<---- in s-1
  vector<double> hHalfLifeTime;   //<---- in h
  vector<double> hProdPerSec;     //<---- nuclei per sec
  vector<double> hYieldParent;    //<---- yield at the EOB
  vector<double> hActivityParent; //<---- activity (mCi) at the EOB

  string s_tmp;
  double x_tmp;
  int i_tmp;
  bool isEoF = false;

  //------------------------ READING THE INPUTS
  G4output.open("Output_ParentIsotopes.txt");
  for(int i=0;i<4;i++)getline(G4output,endLine); //<--- read header.
  while(!isEoF)
    {
      G4output >> s_tmp; getline(G4output,endLine); //<--- name of isotope
      isEoF = G4output.eof();
      if(!isEoF)
	{
	  
	  name.push_back(s_tmp);
	  getline(G4output,endLine); //number of isotopes in the simulation
	  G4output >> x_tmp; getline(G4output,endLine); //<--- decay constant (s-1)
	  hDecayConstant.push_back(x_tmp);
	  G4output >> x_tmp; getline(G4output,endLine); //<--- half life time
	  hHalfLifeTime.push_back(x_tmp);
	  getline(G4output,endLine); //<--- process
	  G4output >> x_tmp; getline(G4output,endLine);  //<--- nuclei per sec
	  hProdPerSec.push_back(x_tmp);
	  G4output >> x_tmp; getline(G4output,endLine);  //<--- yield EOB
	  hYieldParent.push_back(x_tmp);
	  G4output >> x_tmp; getline(G4output,endLine);  //<--- activity EOB
	  hActivityParent.push_back(x_tmp);
	  getline(G4output,endLine); //<--- end of isotope case
	}
    }

  G4output.close();

  
  //------------------------ CALCULATING THE YIELDS/ACTIVITIES
  const int size_parents = hYieldParent.size(); 

  //---> Yield
  TF1* table[size_parents] ; 
  TLegend* leg = new TLegend(0.85,0.35,0.95,0.95);
  double maximum;
  
  //---> Activity
  TF1* tableActivity [size_parents] ; 
  TLegend* legActivity = new TLegend(0.85,0.35,0.95,0.95);
  double maximumActivity = 0;
  stringstream ssTotalActivity;
  
  for(int i=0;i<hYieldParent.size();i++)
    {
      string nameIsotope = name[i];
      double yield = hYieldParent[i];
      double decayConstant = hDecayConstant[i]*3600.;     //<---- s-1 to h-1
      double halfLifeTime = hHalfLifeTime[i];             //<---- h
      double nucleiPerSec = hProdPerSec[i]*3600.;         //<---- nuclei/sec to nuclei/h
      double conv = 2.70E-8;
      
      stringstream titleCanvas;
      stringstream titleHisto;
      stringstream titleLeg;
      
      titleCanvas << nameIsotope << " Production";
      titleHisto << nameIsotope << " production" ;
      titleLeg << nameIsotope ;
      
      ////CALCULATION OF SATURATION/////
      double calculationYield = nucleiPerSec/decayConstant*(1-exp(-decayConstant*tIrradiation));
      double calculationActivity = conv*nucleiPerSec*(1-exp(-decayConstant*tIrradiation)); 
      double timeSaturationCalculation = 10.*log(2)/decayConstant; //in h.
      double saturationYield = nucleiPerSec/decayConstant*(1-exp(-decayConstant*timeSaturationCalculation));
      double saturationActivity = conv*nucleiPerSec*(1-exp(-decayConstant*timeSaturationCalculation));

      //To double check: this calculation must be equal to the yield
      //cout << calculationYield << " should be equal to " << hYieldParent[i] << endl;
      //cout << calculationActivity << " should be equal to " << hActivityParent[i] << endl;
         
      stringstream ssYield,ss1Yield,ssActivity,ss1Activity;

      //STRINGSTREAM FOR NUCLEI PRODUCTION
      ssYield << "(x<="<< tIrradiation << ")*" << nucleiPerSec/decayConstant << "*(1 - exp(-" << decayConstant << "*x)) + (x>" << tIrradiation << ")*" <<  yield << "*exp(-" << decayConstant << "*(x- " << tIrradiation <<"))";
  
      //STRINGSTREAM FOR THE SATURATION
      ss1Yield << nucleiPerSec/decayConstant << "*(1 - exp(-" << decayConstant << "*x))";
     
      //STRINGSTREAM FOR ACTIVITY ACCORDING TO TIME     
      ssActivity << "(x<="<< tIrradiation << ")*" << conv*nucleiPerSec << "*(1 - exp(-" << decayConstant << "*x)) + (x>" << tIrradiation << ")*" << conv*decayConstant*yield << "*exp(-" << decayConstant << "*(x- " << tIrradiation <<"))";

      //STRINGSTREAM FOR ACTIVITY SATURATION      
      ss1Activity << conv*nucleiPerSec << "*(1 - exp(-" << decayConstant << "*x))";

      //TOTAL ACTIVITY
      if(halfLifeTime > halfLifeLimit)
	{
	  if(i == 0){
	    ssTotalActivity << "(x<="<< tIrradiation << ")*" << conv*nucleiPerSec << "*(1 - exp(-" << decayConstant << "*x)) + (x>" << tIrradiation << ")*" << conv*decayConstant*yield << "*exp(-" << decayConstant << "*(x- " << tIrradiation <<"))";
	  }
	  if(i > 0){
	    ssTotalActivity << " + (x<="<< tIrradiation << ")*" << conv*nucleiPerSec << "*(1 - exp(-" << decayConstant << "*x)) + (x>" << tIrradiation << ")*" << conv*decayConstant*yield << "*exp(-" << decayConstant << "*(x- " << tIrradiation <<"))";
	  }
	}

      double max;
      
      //PLOT OF NUCLEI PRODUCTION
      TF1 *fProd = new TF1(titleHisto.str().c_str(),ssYield.str().c_str(),tMin,tMax);
      fProd->SetTitle(titleHisto.str().c_str());
      fProd->GetXaxis()->SetTitle("Time (hour)");
      fProd->GetYaxis()->SetTitle("Number of nuclei");
     
      max = fProd->GetMaximum();
      if(max>maximum){ maximum = max;};


      TF1 *fActivity = new TF1(titleHisto.str().c_str(),ssActivity.str().c_str(),tMin,tMax);
      fActivity->SetTitle(titleHisto.str().c_str());
      fActivity->GetXaxis()->SetTitle("Time (hour)");
      fActivity->GetYaxis()->SetTitle("Activity (mCi)");
        
      max = fActivity->GetMaximum();
      if(max>maximumActivity){ maximumActivity = max;};
      
      leg->AddEntry(fProd,titleLeg.str().c_str());
      table[i]=fProd;

      legActivity->AddEntry(fActivity,titleLeg.str().c_str());
      tableActivity[i]=fActivity;


      //---->Plotting yield as a function of time
      TCanvas *canvasYield = new TCanvas(titleCanvas.str().c_str(),titleCanvas.str().c_str());
      if(halfLifeTime > halfLifeLimit) //has a life time larger than one minute to plot outputs.
	{
	  fProd->Draw();
	  stringstream saveName;
	  saveName << "./Results/IsotopesProduction/YieldOf" << nameIsotope << ".pdf";
	  canvasYield->Print(saveName.str().c_str());
	}

      //---->Plotting activity as a function of time
      TCanvas *canvasActivity = new TCanvas(titleCanvas.str().c_str(),titleCanvas.str().c_str());
      if(halfLifeTime > halfLifeLimit) //has a life time larger than one minute to plot outputs.
	{
	  fActivity->Draw();
	  stringstream saveName;
	  saveName << "./Results/IsotopesProduction/ActivityOf" << nameIsotope << ".pdf";
	  canvasActivity->Print(saveName.str().c_str());
	}


      //PLOT OF SATURATION CURVES
      stringstream titleCanvas1;
      stringstream titleHisto1;       
      titleHisto1 << nameIsotope << " saturation" ;
      titleCanvas1 << nameIsotope << " Saturation";

      TCanvas *canvasSaturationYield = new TCanvas(titleCanvas1.str().c_str(),titleCanvas1.str().c_str());
      TF1 *fSaturationYield = new TF1(titleHisto1.str().c_str(),ss1Yield.str().c_str(),tMin,timeSaturationCalculation);
      fSaturationYield->SetTitle(titleHisto1.str().c_str());
      fSaturationYield->GetXaxis()->SetTitle("Time (hour)");
      fSaturationYield->GetYaxis()->SetTitle("Number of nuclei");
      fSaturationYield->Draw();

      stringstream saveName1Yield;
      saveName1Yield << "./Results/IsotopesProduction/SaturationYieldOf" << nameIsotope << ".pdf";
      canvasSaturationYield->Print(saveName1Yield.str().c_str());

      
      TCanvas *canvasSaturationActivity = new TCanvas(titleCanvas1.str().c_str(),titleCanvas1.str().c_str());
      TF1 *fSaturationActivity = new TF1(titleHisto1.str().c_str(),ss1Activity.str().c_str(),tMin,timeSaturationCalculation);
      fSaturationActivity->SetTitle(titleHisto1.str().c_str());
      fSaturationActivity->GetXaxis()->SetTitle("Time (hour)");
      fSaturationActivity->GetYaxis()->SetTitle("Activity (mCi)");
      fSaturationActivity->Draw();

      stringstream saveName1Activity;
      saveName1Activity << "./Results/IsotopesProduction/ActivitiySaturationOf" << nameIsotope << ".pdf";
      canvasSaturationActivity->Print(saveName1Activity.str().c_str());
      
      
    }
    

  //------------------------ PLOT ALL YIELDS ON THE SAME CANVAS
  TCanvas* productionGraph = new TCanvas("ProductionOfIsotopes", "Production of isotopes");
  for(int i=0; i<size_parents; i++)
    {
      double halfLifeTime = hHalfLifeTime[i];             //<---- h
      if(halfLifeTime > halfLifeLimit)
	{
	  int color = i+2;                                    //<-- color to plot.
	  if(color == 10) color = 35;
 
	  TF1* histo = (TF1*)table[i];
	  histo->GetYaxis()->SetRangeUser(0.,maximum*1.2);
	  histo->SetTitle("Production of isotopes");
	  histo->SetLineColor(color);
	  histo->SetNpx(1000);
      
	  if(i==0)histo->Draw();
	  else histo->Draw("][sames");
	}
    }
  
  leg->Draw("C,same");
  productionGraph->SetGridy();
  productionGraph->SetTicky();
  productionGraph->SetLogy();
  productionGraph->SetTitle("Radioisotope production");
  productionGraph->Print("./Results/IsotopesProduction/Yield.pdf");
  productionGraph->Print("./Results/IsotopesProduction/Yield.jpg");


    
  //------------------------ PLOT ALL YIELDS ON THE SAME CANVAS
  TCanvas* ActivityGraph = new TCanvas("ActivityOfIsotopes", "Activity of isotopes");
  
  TF1* histoActivity = (TF1*)tableActivity[0];
  histoActivity->SetLineColor(2);
  histoActivity->GetYaxis()->SetRangeUser(0.,maximumActivity*1.2);
  histoActivity->SetTitle("Activity of isotopes");
  histoActivity->Draw();
  
  for(int i=1; i<size_parents; i++)
    {
      double halfLifeTime = hHalfLifeTime[i];             //<---- h
      if(halfLifeTime > halfLifeLimit)
	{
	  int color = i+2;                                    //<-- color to plot.
	  if(color == 10) color = 35;
 
	  TF1* histo = (TF1*)tableActivity[i];
	  histo->GetYaxis()->SetRangeUser(0.,maximumActivity*1.2);
	  histo->SetTitle("Activity of isotopes");
	  histo->SetLineColor(color);
	  histo->SetNpx(1000);
	  if(i==0) histo->Draw();
	  else histo->Draw("][sames");
	}
    }
  
  legActivity->Draw("same");
  ActivityGraph->SetGridy();
  ActivityGraph->SetTicky();
  ActivityGraph->SetLogy();
  ActivityGraph->Print("./Results/IsotopesProduction/Activity.pdf");
  ActivityGraph->Print("./Results/IsotopesProduction/Activity.jpg");
   
  
  TCanvas* TotalActivityGraph = new TCanvas("TotalActivity", "Total Activity");
  TF1 *fActivity = new TF1("TotalActivity",ssTotalActivity.str().c_str(),tMin,tMax);
  fActivity->SetTitle("Sum of the activity of each isotope");
  fActivity->GetXaxis()->SetTitle("Time (hour)");
  fActivity->GetYaxis()->SetTitle("Activity (mCi)");
  fActivity->Draw();
  TotalActivityGraph->SetGridy();
  TotalActivityGraph->SetTicky();
  TotalActivityGraph->SetLogy();
  TotalActivityGraph->Print("./Results/IsotopesProduction/TotalActivity.pdf");
   


  
  //----------------------------------------------------------------------------
  //                 CASE OF DAUGHTER ISOTOPES
  //----------------------------------------------------------------------------

  vector<string> nameParent;            //<---- name of the isotope
  vector<string> nameDaughter;            //<---- name of the isotope
  vector<double> hDecayConstantParent;   //<---- in s-1
  vector<double> hDecayConstantDaughter;   //<---- in s-1
  vector<double> hHalfLifeTimeParent;   //<---- in h
  vector<double> hHalfLifeTimeDaughter;   //<---- in h
  vector<double> hProdPerSecDecay;     //<---- nuclei per sec
  vector<double> hYieldDecay;    //<---- yield at the EOB
  vector<double> hActivityDecay; //<---- activity (mCi) at the EOB


  //------------------------ READING THE INPUTS
  G4output.open("Output_DaughterIsotopes.txt");
  isEoF=false;

  for(int i=0;i<5;i++)getline(G4output,endLine); //<--- read header.
  while(!isEoF)
    {
      G4output >> s_tmp; getline(G4output,endLine); //<--- name of daughter isotope
      isEoF = G4output.eof();
      if(!isEoF)
	{
	  
	  nameDaughter.push_back(s_tmp); //cout << s_tmp << " + " << endLine << endl;
	  G4output >> s_tmp; getline(G4output,endLine); //<--- name of parent isotope
	  nameParent.push_back(s_tmp); //cout << s_tmp << " + " << endLine << endl;
	  G4output >> x_tmp; getline(G4output,endLine); //<--- decay constant (s-1) daughter
	  hDecayConstantDaughter.push_back(x_tmp); //cout << x_tmp << " + " << endLine << endl;
	  G4output >> x_tmp; getline(G4output,endLine); //<--- decay constant (s-1) parent
	  hDecayConstantParent.push_back(x_tmp); //cout << x_tmp << " + " << endLine << endl;
	  G4output >> x_tmp; getline(G4output,endLine); //<--- half life time
	  hHalfLifeTimeDaughter.push_back(x_tmp); //cout << x_tmp << " + " << endLine << endl;
	  G4output >> x_tmp; getline(G4output,endLine); //<--- half life time
	  hHalfLifeTimeParent.push_back(x_tmp); //cout << x_tmp << " + " << endLine << endl;
	  G4output >> x_tmp; getline(G4output,endLine);  //<--- nuclei per sec
	  hProdPerSecDecay.push_back(x_tmp); //cout << x_tmp << " + " << endLine << endl;
	  G4output >> x_tmp; getline(G4output,endLine);  //<--- yield EOB
	  hYieldDecay.push_back(x_tmp); //cout << x_tmp << " + " << endLine << endl;
	  G4output >> x_tmp; getline(G4output,endLine);  //<--- activity EOB
	  hActivityDecay.push_back(x_tmp); //cout << x_tmp << " + " << endLine << endl;
	  getline(G4output,endLine); //<--- end of isotope case
	  //getchar();
	}
    }

  G4output.close();

  //------------------------ CALCULATING THE YIELDS/ACTIVITIES
  for(int i=1;i<=hDecayConstantDaughter.size();i++){
    
    string nameIsotope = nameDaughter[i];
    double yieldEOBDaughter = hYieldDecay[i];
    double decayConstantDaughter = hDecayConstantDaughter[i]*3600.;
    double decayConstantParent = hDecayConstantParent[i]*3600;
    double halfLifeDaughter = hHalfLifeTimeDaughter[i]*3600.;
    double halfLifeParent = hHalfLifeTimeParent[i]*3600;
    double nucleiPerSec = hProdPerSecDecay[i]*3600.;
    double yieldEOBParent = nucleiPerSec/decayConstantParent*(1-exp(-decayConstantParent*tIrradiation));

     
    if(halfLifeDaughter > halfLifeLimit && halfLifeParent > halfLifeLimit){
      
      stringstream titleCanvas;
      stringstream titleHisto;
      
      titleCanvas << nameIsotope << " Decay";
      titleHisto << nameIsotope << " Decay" ;
   
      double yieldEOBcalc = nucleiPerSec*((1-exp(-decayConstantDaughter*tIrradiation))/decayConstantDaughter + (exp(-decayConstantDaughter*tIrradiation)-exp(-decayConstantParent*tIrradiation))/(decayConstantDaughter - decayConstantParent));

      cout << "Isotope : " << nameIsotope << " with yield at the EOB " << yieldEOBDaughter << " calculation : " << yieldEOBcalc
    	 << " decay constant : " << decayConstantDaughter << " parent decay constant : " << decayConstantParent
   	 << " nucleiPerSec of the parent " << nucleiPerSec << " yieldEOBParent " << yieldEOBParent << endl;

   
      TCanvas *canvasYield = new TCanvas(titleCanvas.str().c_str(),titleCanvas.str().c_str());
     

      stringstream ss;
            
      ss << "(x<="<< tIrradiation << ")*" << nucleiPerSec << "*((1 - exp(-" << decayConstantDaughter << "*x))/" << decayConstantDaughter << " + (exp(-" << decayConstantDaughter << "*x)-exp(-" << decayConstantParent << "*x))/(" << decayConstantDaughter-decayConstantParent << "))+ (x>" << tIrradiation << ")*" <<  yieldEOBcalc << "*exp(-" << decayConstantDaughter << "*(x - " << tIrradiation << ")) + " << decayConstantParent*yieldEOBParent/(decayConstantDaughter-decayConstantParent) << "*(exp(- " << decayConstantParent << "*(x - " << tIrradiation << ")) - exp(- " << decayConstantDaughter << "*(x-" << tIrradiation << ")))";
      
      
      TF1 *fProd = new TF1(titleHisto.str().c_str(),ss.str().c_str(),tMin, tMax);
      fProd->SetTitle(titleHisto.str().c_str());
      fProd->GetXaxis()->SetTitle("Time (hour)");
      fProd->GetYaxis()->SetTitle("Number of nuclei");

    }
      
  }
  
  f->Close();
  results.close();


}


