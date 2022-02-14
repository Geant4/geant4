
void BinLogX(TH1* h)
{
				TAxis *axis = h->GetXaxis();
				int bins = axis->GetNbins();

				Axis_t from      = axis->GetXmin();
				Axis_t to        = axis->GetXmax();
				Axis_t width     = (to - from) / bins;
				Axis_t *new_bins = new Axis_t[bins + 1];

				for (int i = 0; i <= bins; i++) {
								new_bins[i] = TMath::Power(10, from + i * width);
				}
				axis->Set(bins, new_bins);
				delete []new_bins;
}


void plot(){

				double dens = 1000;//kg/m3

				double R =50;

				//TFile *fin           = new TFile("AuNP_Livermore.root");
				TFile *fin           = new TFile("AuNP.root");
				TH1F *h1neve         = (TH1F*)fin->Get("h1Events");
				TH1F *h1Edep         = (TH1F*)fin->Get("h1Edep");
				TH1F *h1secnp_cha    = (TH1F*)fin->Get("h1SecEnergyNP_charged");
				TH1F *h1secnp_nut    = (TH1F*)fin->Get("h1SecEnergyNP_nutral");
				TH1F *h1secnpsurf_cha= (TH1F*)fin->Get("h1SecEnergyNPSurf_charged");
				TH1F *h1secnpsurf_nut= (TH1F*)fin->Get("h1SecEnergyNPSurf_nutral");
				TH1F *h1sec_cha      = (TH1F*)fin->Get("h1Sec_charged");
				TH1F *h1sec_nut      = (TH1F*)fin->Get("h1Sec_nutral");
				TH1F *h1chem_0       = (TH1F*)fin->Get("h1Chem_0");
				TH1F *h1chem_1       = (TH1F*)fin->Get("h1Chem_1");
				TH1F *h1chem_2       = (TH1F*)fin->Get("h1Chem_2");
				TH1F *h1chem_3       = (TH1F*)fin->Get("h1Chem_3");
				TH1F *h1chem_4       = (TH1F*)fin->Get("h1Chem_4");
				TH1F *h1chem_5       = (TH1F*)fin->Get("h1Chem_5");
				TH1F *h1chem_6       = (TH1F*)fin->Get("h1Chem_6");
				TH1F *h1chem_7       = (TH1F*)fin->Get("h1Chem_7");

				TH2F *h2Edep         = (TH2F*)fin->Get("h2Edep");
				TH2F *h2sec2_cha     = (TH2F*)fin->Get("h2SecEnergyAbs_charged");
				TH2F *h2sec2_nut     = (TH2F*)fin->Get("h2SecEnergyAbs_nutral");

				double neve = h1neve->GetBinContent(1);
				h1Edep         ->Scale(1./neve);
				h1secnp_cha    ->Scale(1./neve);
				h1secnp_nut    ->Scale(1./neve);
				h1secnpsurf_cha->Scale(1./neve);
				h1secnpsurf_nut->Scale(1./neve);
				h1sec_cha      ->Scale(1./neve);
				h1sec_nut      ->Scale(1./neve);
				h1chem_0       ->Scale(1./neve);
				h1chem_1       ->Scale(1./neve);
				h1chem_2       ->Scale(1./neve);
				h1chem_3       ->Scale(1./neve);
				h1chem_4       ->Scale(1./neve);
				h1chem_5       ->Scale(1./neve);
				h1chem_6       ->Scale(1./neve);
				h1chem_7       ->Scale(1./neve);
				h2Edep         ->Scale(1./neve);
				h2sec2_cha     ->Scale(1./neve);
				h2sec2_nut     ->Scale(1./neve);

				double val_cha=0;
				double err_cha=0;
				double val_nut=0;
				double err_nut=0;

				int NR = h1sec_cha  -> GetNbinsX();
				TH1D *h1sec_tot = new TH1D("h1Sec_tot","h1Sec_tot",NR,1,5);
        BinLogX(h1sec_tot);
				for(int i=0;i<NR;i++){
								val_cha = h1sec_cha->GetBinContent(i+1);
								err_cha = h1sec_cha->GetBinError  (i+1);
								val_nut = h1sec_nut->GetBinContent(i+1);
								err_nut = h1sec_nut->GetBinError  (i+1);
								//val_cha = h1sec_cha->GetBinContent(i+1)/h1sec_cha->GetBinWidth(i+1);
								//err_cha = h1sec_cha->GetBinError  (i+1)/h1sec_cha->GetBinWidth(i+1);
								//val_nut = h1sec_nut->GetBinContent(i+1)/h1sec_nut->GetBinWidth(i+1);
								//err_nut = h1sec_nut->GetBinError  (i+1)/h1sec_nut->GetBinWidth(i+1);
								h1sec_cha  -> SetBinContent(i+1,val_cha);
								h1sec_cha  -> SetBinError  (i+1,err_cha);
								h1sec_nut  -> SetBinContent(i+1,val_nut);
								h1sec_nut  -> SetBinError  (i+1,err_nut);
								h1sec_tot  -> SetBinContent(i+1,val_cha+val_nut);
								h1sec_tot  -> SetBinError  (i+1,err_cha+err_nut);
				}

				TH1D *h1chem_tot = new TH1D("h1Chem_tot","h1Chem_tot",NR,1,5);
        BinLogX(h1chem_tot);
				for(int i=0;i<NR;i++){
								double val_0 = h1chem_0->GetBinContent(i+1)/h1chem_0->GetBinWidth(i+1);
								double val_1 = h1chem_1->GetBinContent(i+1)/h1chem_1->GetBinWidth(i+1);
								double val_2 = h1chem_2->GetBinContent(i+1)/h1chem_2->GetBinWidth(i+1);
								double val_3 = h1chem_3->GetBinContent(i+1)/h1chem_3->GetBinWidth(i+1);
								double val_4 = h1chem_4->GetBinContent(i+1)/h1chem_4->GetBinWidth(i+1);
								double val_5 = h1chem_5->GetBinContent(i+1)/h1chem_5->GetBinWidth(i+1);
								double val_6 = h1chem_6->GetBinContent(i+1)/h1chem_6->GetBinWidth(i+1);
								double val_7 = h1chem_7->GetBinContent(i+1)/h1chem_7->GetBinWidth(i+1);
								double err_0 = h1chem_0->GetBinError  (i+1)/h1chem_0->GetBinWidth(i+1);
								double err_1 = h1chem_1->GetBinError  (i+1)/h1chem_1->GetBinWidth(i+1);
								double err_2 = h1chem_2->GetBinError  (i+1)/h1chem_2->GetBinWidth(i+1);
								double err_3 = h1chem_3->GetBinError  (i+1)/h1chem_3->GetBinWidth(i+1);
								double err_4 = h1chem_4->GetBinError  (i+1)/h1chem_4->GetBinWidth(i+1);
								double err_5 = h1chem_5->GetBinError  (i+1)/h1chem_5->GetBinWidth(i+1);
								double err_6 = h1chem_6->GetBinError  (i+1)/h1chem_6->GetBinWidth(i+1);
								double err_7 = h1chem_7->GetBinError  (i+1)/h1chem_7->GetBinWidth(i+1);
								h1chem_tot  -> SetBinContent(i+1,val_0+val_1+val_2+val_3+val_4+val_5+val_6+val_7);
								h1chem_tot  -> SetBinError  (i+1,err_0+err_1+err_2+err_3+err_4+err_5+err_6+err_7);
				}


				NR = h1secnp_cha -> GetNbinsX();
				TH1D *h1secnp_tot     = new TH1D("h1SecEnergyNP"    ,"Energy Spectra in NP"        ,NR,0,6);
        BinLogX(h1secnp_tot);
				NR = h1secnpsurf_cha -> GetNbinsX();
				TH1D *h1secnpsurf_tot = new TH1D("h1SecEnergyNPSurf","Energy Spectra at NP surface",NR,0,6);
        BinLogX(h1secnpsurf_tot);
				for(int i=0;i<NR;i++){
								val_cha = h1secnp_cha->GetBinContent(i+1);
								err_cha = h1secnp_cha->GetBinError  (i+1);
								val_nut = h1secnp_nut->GetBinContent(i+1);
								err_nut = h1secnp_nut->GetBinError  (i+1);
								//val_cha = h1secnp_cha->GetBinContent(i+1)/h1secnp_cha->GetBinWidth(i+1);
								//err_cha = h1secnp_cha->GetBinError  (i+1)/h1secnp_cha->GetBinWidth(i+1);
								//val_nut = h1secnp_nut->GetBinContent(i+1)/h1secnp_nut->GetBinWidth(i+1);
								//err_nut = h1secnp_nut->GetBinError  (i+1)/h1secnp_nut->GetBinWidth(i+1);
								h1secnp_cha -> SetBinContent(i+1,val_cha);
								h1secnp_cha -> SetBinError  (i+1,err_cha);
								h1secnp_nut -> SetBinContent(i+1,val_nut);
								h1secnp_nut -> SetBinError  (i+1,err_nut);
								h1secnp_tot -> SetBinContent(i+1,val_cha+val_nut);
								h1secnp_tot -> SetBinError  (i+1,err_cha+err_nut);


								val_cha = h1secnpsurf_cha->GetBinContent(i+1);
								err_cha = h1secnpsurf_cha->GetBinError  (i+1);
								val_nut = h1secnpsurf_nut->GetBinContent(i+1);
								err_nut = h1secnpsurf_nut->GetBinError  (i+1);
								//val_cha = h1secnpsurf_cha->GetBinContent(i+1)/h1secnpsurf_cha->GetBinWidth(i+1);
								//err_cha = h1secnpsurf_cha->GetBinError  (i+1)/h1secnpsurf_cha->GetBinWidth(i+1);
								//val_nut = h1secnpsurf_nut->GetBinContent(i+1)/h1secnpsurf_nut->GetBinWidth(i+1);
								//err_nut = h1secnpsurf_nut->GetBinError  (i+1)/h1secnpsurf_nut->GetBinWidth(i+1);
								h1secnpsurf_cha -> SetBinContent(i+1,val_cha);
								h1secnpsurf_cha -> SetBinError  (i+1,err_cha);
								h1secnpsurf_nut -> SetBinContent(i+1,val_nut);
								h1secnpsurf_nut -> SetBinError  (i+1,err_nut);
								h1secnpsurf_tot -> SetBinContent(i+1,val_cha+val_nut);
								h1secnpsurf_tot -> SetBinError  (i+1,err_cha+err_nut);

				}

				NR = h1Edep  -> GetNbinsX();
				for(int i=0;i<NR;i++){
								double rmax = h1Edep->GetXaxis()->GetBinUpEdge (i+1)*1.E-9;//m
								double rmin = h1Edep->GetXaxis()->GetBinLowEdge(i+1)*1.E-9;//m
								double vol  = 4./3.*TMath::Pi()*(pow(rmax,3)-pow(rmin,3));
								double mass = dens*vol;
								h1Edep  -> SetBinContent(i+1,h1Edep->GetBinContent(i+1)/mass);
								h1Edep  -> SetBinError  (i+1,h1Edep->GetBinError  (i+1)/mass);
				}


				int Nazm =h2Edep->GetXaxis()->GetNbins()-1;
				NR   =h2Edep->GetYaxis()->GetNbins()-1;
				double    Rmax =h2Edep->GetYaxis()->GetBinUpEdge(NR);
				TH2D *h2Edep_pol = new TH2D("h2Edep_pol","h2Edep_pol",Nazm,0,360,NR,0,Rmax);
				for(int i=0;i<Nazm;i++){
								for(int j=0;j<NR;j++){
												double height = 10*1.E-9;//m
												double rmax = h2Edep_pol->GetYaxis()->GetBinUpEdge  (j+1)*1.E-9;//m
												double rmin = h2Edep_pol->GetYaxis()->GetBinLowEdge (j+1)*1.E-9;//m
												double vol  = (TMath::Pi()*pow(rmax,2)-TMath::Pi()*pow(rmin,2))*height;
								        double mass = dens*vol;
												h2Edep_pol->SetBinContent(i+1,j+1,h2Edep->GetBinContent(i+1,j+1)/mass);
								}
				}
				int NRX = h2sec2_cha->GetXaxis()->GetNbins();
				int NRY = h2sec2_cha->GetYaxis()->GetNbins();
				TH2D *h2sec2_tot = (TH2D*)h2sec2_cha->Clone("h2sec2_tot");
				for(int i=0;i<NRX;i++){
								for(int j=0;j<NRY;j++){
												val_cha = h2sec2_cha->GetBinContent(i+1,j+1);
												val_nut = h2sec2_nut->GetBinContent(i+1,j+1);
												h2sec2_tot->SetBinContent(i+1,j+1,val_cha+val_nut);
								}
				}







				TCanvas *cedep = new TCanvas("cedep","cedep",700,700);
				cedep->cd(1)->SetTheta(90);
				cedep->cd(1)->SetPhi(0);
				cedep->cd(1)->SetLogz();
				h2Edep_pol->RebinX(6);
				h2Edep_pol->RebinY(10);
				//h2Edep_pol->SetAxisRange(0,1000,"Y");
				h2Edep_pol->SetAxisRange(0,150,"Y");
				h2Edep_pol->SetTitle("");
				h2Edep_pol->SetXTitle("Angle [deg]");
				h2Edep_pol->Draw("surf1POL");

				TCanvas *cedepsec = new TCanvas("cedepsec","cedepsec",1800,700);
				cedepsec->Divide(2,1);
				cedepsec->cd(1)->SetLogy();
				cedepsec->cd(1)->SetLogx();
				h1Edep->SetYTitle("Energy Deposite [Gy/event]");
				h1Edep->SetXTitle("Distance from NP center[nm]");
				h1Edep->SetMaximum(100);
				h1Edep->SetMinimum(0.0000000001);
				h1Edep->SetAxisRange(0,1000000);
				h1Edep->Draw();
				cedepsec->cd(2)->SetLogy();
				cedepsec->cd(2)->SetLogx();
				h1sec_tot->SetLineColor(kBlack);
				h1sec_cha->SetLineColor(kRed);
				h1sec_nut->SetLineColor(kBlue);
				h1sec_tot->SetMarkerColor(kBlack);
				h1sec_cha->SetMarkerColor(kRed);
				h1sec_nut->SetMarkerColor(kBlue);
				h1sec_tot->SetAxisRange(0,1000000);
				h1sec_tot->SetMaximum(10);
				h1sec_tot->SetMinimum(0.000001);
				h1sec_tot->SetTitle ("Number of Generated Secondaries");
				h1sec_tot->SetYTitle("N/nm/event");
				h1sec_tot->SetXTitle("Distance from NP center[nm]");
				//h1sec_cha->Draw("sameL");
				//h1sec_nut->Draw("sameL");
				h1sec_tot->Draw();
				//cedepsec->cd(3)->SetLogy();
				//cedepsec->cd(3)->SetLogx();
				//h1chem_tot->SetLineColor  (kBlack);
				//h1chem_tot->SetMarkerColor(kBlack);
				//h1chem_tot->SetAxisRange(0,1000000);
				//h1chem_tot->SetMaximum(10);
				//h1chem_tot->SetMinimum(0.000001);
				//h1chem_tot->SetTitle ("Number of Generated Chemical Species");
				//h1chem_tot->SetYTitle("N/nm/event");
				//h1chem_tot->SetXTitle("Distance from NP center[nm]");
				//h1chem_tot->Draw();


				TCanvas *csec = new TCanvas("csec","csec",2400,700);
				csec->Divide(3,1);
				csec->cd(1)->SetLogx();
				csec->cd(1)->SetLogy();
				h1secnp_tot->SetLineColor(kBlack);
				h1secnp_cha->SetLineColor(kRed);
				h1secnp_nut->SetLineColor(kBlue);
				h1secnp_tot->SetMarkerColor(kBlack);
				h1secnp_cha->SetMarkerColor(kRed);
				h1secnp_nut->SetMarkerColor(kBlue);
				h1secnp_tot->RebinX(10);
				h1secnp_tot->Scale(1./10.);
				h1secnp_tot->SetMaximum(20);
				h1secnp_tot->SetMinimum(0.000000001);
				h1secnp_tot->SetAxisRange(1,1.E6);
				h1secnp_tot->SetTitle ("Secondary Energy produced in NP");
				h1secnp_tot->SetYTitle("N/event");
				h1secnp_tot->SetXTitle("Secandary Energy [eV]");
				//h1secnp_cha->DrawCopy();
				//h1secnp_nut->DrawCopy("same");
				h1secnp_tot->DrawCopy();
				csec->cd(2)->SetLogx();
				csec->cd(2)->SetLogy();
				h1secnpsurf_tot->SetLineColor(kBlack);
				h1secnpsurf_cha->SetLineColor(kRed);
				h1secnpsurf_nut->SetLineColor(kBlue);
				h1secnpsurf_tot->SetMarkerColor(kBlack);
				h1secnpsurf_cha->SetMarkerColor(kRed);
				h1secnpsurf_nut->SetMarkerColor(kBlue);
				h1secnpsurf_tot->RebinX(10);
				h1secnpsurf_tot->Scale(1./10.);
				h1secnpsurf_tot->SetMaximum(20);
				h1secnpsurf_tot->SetMinimum(0.000000001);
				h1secnpsurf_tot->SetAxisRange(1,1.E6);
				h1secnpsurf_tot->SetTitle ("Secondary Energy at NP surface");
				h1secnpsurf_tot->SetYTitle("N/event");
				h1secnpsurf_tot->SetXTitle("Secandary Energy [eV]");
				//h1secnpsurf_cha->DrawCopy();
				//h1secnpsurf_nut->DrawCopy("same");
				h1secnpsurf_tot->DrawCopy();
				csec->cd(3)->SetLogy();
				h2sec2_tot->SetTitle ("Secondary Energy vs production point");
				h2sec2_tot->SetYTitle("Secondary Energy [eV]");
				h2sec2_tot->SetXTitle("Distance from NP center [nm]");
				h2sec2_tot->GetXaxis()->SetRange(0,1000);
				h2sec2_tot->Draw("colz");

}

