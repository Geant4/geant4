#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TString.h"

#include "expdata_de.h"
#include "expdata_da.h"
#include "expdata_dd.h"

#include <map>
#include <set>

void Summarize(TTree * tDE, TTree * tDA, TTree * tDD)
{
    std::map<Int_t,std::set<Double_t>* > summary;

    expdata_de * dedata = new expdata_de;
    tDE->SetBranchAddress("DE", &dedata);
    Int_t nentradas = Int_t(tDE->GetEntries());
    Int_t theZ(-1);
    if (nentradas > 0)
      {
	tDE->GetEntry(0);
	theZ = dedata->GetData()->GetTargetZ();
	for (Int_t i = 0; i < nentradas; i++)
	  {
	    tDE->GetEntry(i);
	    Int_t TA = dedata->GetData()->GetTargetA();
	    Double_t E = dedata->GetData()->GetProjectileEnergy();
	    if (summary.find(TA) != summary.end())
	      {
		summary[TA]->insert(E);
	      }
	    else
	      {
		summary.insert(make_pair(TA,new std::set<Double_t>));
		summary[TA]->insert(E);
	      }
	  }
      }
    delete dedata;
    expdata_da * dadata = new expdata_da;
    tDA->SetBranchAddress("DA", &dadata);
    nentradas = Int_t(tDA->GetEntries());
    if (nentradas > 0)
      {
	if (theZ < 0) 
	  {
	    tDA->GetEntry(0);
	    theZ = dadata->GetData()->GetTargetZ();
	  }
	for (Int_t i = 0; i < nentradas; i++)
	  {
	    tDA->GetEntry(i);
	    Int_t TA = dadata->GetData()->GetTargetA();
	    Double_t E = dadata->GetData()->GetProjectileEnergy();
	    if (summary.find(TA) != summary.end())
	      {
		summary[TA]->insert(E);
	      }
	    else
	      {
		summary.insert(make_pair(TA,new std::set<Double_t>));
		summary[TA]->insert(E);
	      }
	  }
      }
    delete dadata;

    expdata_dd * dddata = new expdata_dd;
    tDD->SetBranchAddress("DD", &dddata);
    nentradas = Int_t(tDD->GetEntries());
    if (nentradas > 0)
      {
	if (theZ < 0)
	  {
	    tDD->GetEntry(0);
	    theZ = dddata->GetHeader()->GetTargetZ();
	  }
	for (Int_t i = 0; i < nentradas; i++)
	  {
	    tDD->GetEntry(i);
	    Int_t TA = dddata->GetHeader()->GetTargetA();
	    Double_t E = dddata->GetHeader()->GetProjectileEnergy();
	    if (summary.find(TA) != summary.end())
	      {
		summary[TA]->insert(E);
	      }
	    else
	      {
	    summary.insert(make_pair(TA,new std::set<Double_t>));
	    summary[TA]->insert(E);
	      }
	  }
      }
    delete dddata;
    
    std::map<Int_t,std::set<Double_t>* >::iterator it;   
    for (it = summary.begin(); it != summary.end(); it++)
      {
	std::cout << "Z = " << theZ 
		  << " A = " << it->first << "\n\tEnergies (MeV): ";
	for (std::set<Double_t>::iterator k = it->second->begin();
	     k != it->second->end(); ++k)
	  {
	    std::cout << *k << ", ";
	  }
	std::cout << "\b\b \n";
	it->second->clear();
	delete it->second;
    }
    summary.clear();
    return;
}

void showdata_de(TTree * tDE, TApplication * showdata_de_app, bool doplot = true)
{
    expdata_de * dedata = new expdata_de;
    tDE->SetBranchAddress("DE", &dedata);
    Int_t eleccion = -20;
    Int_t nentradas = Int_t(tDE->GetEntries());
    do 
    {
	std::cout << "\nThere are " << nentradas << " entries, show entry (";
	for (Int_t i = 0; i < nentradas; i++) std::cout << i << ",";
	std::cout << nentradas << "=all,-1=back): ";
	std::cin >> eleccion;
	if (eleccion >= nentradas) 
	{
	    for (Int_t i = 0; i < nentradas; i++) 
	    {
		std::cout << "**************\n";
		std::cout << "* Entry " << i << '\n';
		std::cout << "**************\n";
		tDE->GetEntry(i);
		dedata->ShowYourSelf();
		std::cout << "-------------------------\n";
	    }
	} 
	else if (eleccion > -1 && eleccion < nentradas) 
	{
	    std::cout << "**************\n";
	    std::cout << "* Entry " << eleccion << '\n';
	    std::cout << "**************\n";
	    tDE->GetEntry(eleccion);
	    dedata->ShowYourSelf();
	    std::cout << "-------------------------\n";
	    if (doplot) 
	    {
		TCanvas * aCanvas = new TCanvas("aCanvas","showdata DE");
		gROOT->SetStyle("Plain");
		aCanvas->SetGrid(1,1);
		TGraphErrors * gra = (TGraphErrors*)(dedata->GetGraph());
		gra->SetMarkerColor(2);
		gra->SetMarkerStyle(20);
		gra->SetMarkerSize(0.75);
		TString titulo = dedata->GetData()->GetGraphTitle() + " (" + dedata->GetData()->GetExforEntryCode() + ")";
		gra->SetTitle(titulo);
		gra->Draw("ALP");
		gra->GetXaxis()->SetTitle("T (MeV)");
		gra->GetYaxis()->SetTitle("#frac{d#sigma}{dT}  (mb/MeV)");
		gra->Draw("ALP");
		aCanvas->Update();
		showdata_de_app->Run(kTRUE);
		delete aCanvas;
	    }
	}
    } 
    while (eleccion != -1);

    delete dedata;
    return;
}

void showdata_da(TTree * tDA, TApplication * showdata_da_app, bool doplot = true)
{
    expdata_da * dadata = new expdata_da;
    tDA->SetBranchAddress("DA", &dadata);
    Int_t eleccion = -20;
    Int_t nentradas = Int_t(tDA->GetEntries());
    do 
    {
	std::cout << "\nThere are " << nentradas << " entries, show entry (";
	for (Int_t i = 0; i < nentradas; i++) std::cout << i << ",";
	std::cout << nentradas << "=all,-1=back): ";
	std::cin >> eleccion;
	if (eleccion >= nentradas) 
	{
	    for (Int_t i = 0; i < nentradas; i++) 
	    {
		std::cout << "**************\n";
		std::cout << "* Entry " << i << '\n';
		std::cout << "**************\n";
		tDA->GetEntry(i);
		dadata->ShowYourSelf();
		std::cout << "-------------------------\n";
	    }
	} 
	else if (eleccion > -1 && eleccion < nentradas) 
	{
	    std::cout << "**************\n";
	    std::cout << "* Entry " << eleccion << '\n';
	    std::cout << "**************\n";
	    tDA->GetEntry(eleccion);
	    dadata->ShowYourSelf();
	    std::cout << "-------------------------\n";
	    if (doplot) 
	    {
		TCanvas * aCanvas = new TCanvas("aCanvas","showdata DA");
		gROOT->SetStyle("Plain");
		aCanvas->SetGrid(1,1);
		TGraphErrors * gra = (TGraphErrors*)(dadata->GetGraph());
		gra->SetMarkerColor(2);
		gra->SetMarkerStyle(20);
		gra->SetMarkerSize(0.75);
		TString titulo = dadata->GetData()->GetGraphTitle() + " (" + dadata->GetData()->GetExforEntryCode() + ")";
		gra->SetTitle(titulo);
		gra->Draw("ALP");
		gra->GetXaxis()->SetTitle("#theta (deg)");
		gra->GetYaxis()->SetTitle("#frac{d#sigma}{d#theta}  (mb/rad)");
		gra->Draw("ALP");
		aCanvas->Update();
		showdata_da_app->Run(kTRUE);
		delete aCanvas;
	    }
	}
    } while (eleccion != -1);
  
    delete dadata;
    return;
}



void showdata_dd(TTree * tDD, TApplication * showdata_dd_app, bool doplot = true)
{
    expdata_dd * dddata = new expdata_dd();
    tDD->SetBranchAddress("DD", &dddata);
    Int_t eleccion = -20;
    Int_t nentradas = Int_t(tDD->GetEntries());
    do 
    {
	std::cout << "\nThere are " << nentradas << " entries, show entry (";
	for (Int_t i = 0; i < nentradas; i++) std::cout << i << ",";
	std::cout << nentradas << "=all,-1=back): ";
	std::cin >> eleccion;
	if (eleccion >= nentradas) 
	{
	    for (Int_t i = 0; i < nentradas; i++) 
	    {
		std::cout << "**************\n";
		std::cout << "* Entry " << i << '\n';
		std::cout << "**************\n";
		tDD->GetEntry(i);
		dddata->ShowYourSelf(-1);
		std::cout << "-------------------------\n";
	    }
	}
	else if (eleccion > -1 && eleccion < nentradas) 
	{
	    std::cout << "**************\n";
	    std::cout << "* Entry " << eleccion << '\n';
	    std::cout << "**************\n";
	    tDD->GetEntry(eleccion);
	    dddata->ShowYourSelf(-1);
	    std::cout << "-------------------------\n";
	    Int_t eleccion2 = -20;
	    do 
	    {
		std::cout << "\nShow Angle (";
		for (Int_t j=0; j < dddata->GetNangles(); j++) std::cout << j << ",";
		std::cout << dddata->GetNangles() << "=all,-1=back): ";
		std::cin >> eleccion2;
		if (eleccion2 >= dddata->GetNangles()) 
		{
		    TMultiGraph * mg = 0;
		    if (doplot) mg = new TMultiGraph();
		    int color = 2;
		    for (Int_t k = 0; k < dddata->GetNangles(); k++) 
		    {
			dddata->ShowYourSelf(k);
			std::cout << "-------------------------\n";
			if (doplot)
			{
			    TGraph * g = dddata->GetData(k)->GetGraph();
			    g->SetMarkerColor(color++);
			    g->SetMarkerStyle(20);
			    g->SetMarkerSize(0.75);
			    mg->Add(g);
			}
		    }
		    if (doplot)
		    {
			TCanvas * aCanvas = new TCanvas("aCanvas","showdata DD");
			gROOT->SetStyle("Plain");
			aCanvas->SetGrid(1,1);
			TString titulo = dddata->GetData(0)->GetGraphTitle();
			titulo.Remove(titulo.Index("MeV")+3);
			titulo += " #theta = ";
			char angulo[5];
			for (int idx = 0; idx < dddata->GetNangles()-1; idx++)
			{
			    sprintf(angulo,"%.1f",dddata->GetData(idx)->GetAngle());
			    titulo += TString(angulo) + ",";
			}
			sprintf(angulo,"%.1f",dddata->GetData(dddata->GetNangles()-1)->GetAngle());
			titulo += TString(angulo);
			mg->SetTitle(titulo);
			mg->Draw("ALP");
			mg->GetXaxis()->SetTitle("T (MeV)");
			mg->GetYaxis()->SetTitle("#frac{d#sigma}{dT d#theta}  (mb/MeV/rad)");
			mg->Draw("ALP");
			aCanvas->Update();
			showdata_dd_app->Run(kTRUE);
			delete aCanvas;
		    }
		    if (mg) delete mg;
		} 
		else if (eleccion2 > -1 && eleccion2 < dddata->GetNangles())
		{
		    tDD->GetEntry(eleccion);
		    dddata->ShowYourSelf(eleccion2);
		    std::cout << "-------------------------\n";
		    if (doplot) 
		    {
			TCanvas * aCanvas = new TCanvas("aCanvas","showdata DD");
			gROOT->SetStyle("Plain");
			aCanvas->SetGrid(1,1);
			TGraphErrors * gra = (TGraphErrors*)(dddata->GetData(eleccion2)->GetGraph());
			gra->SetMarkerColor(2);
			gra->SetMarkerStyle(20);
			gra->SetMarkerSize(0.75);
			TString titulo = dddata->GetData(eleccion2)->GetGraphTitle() + 
			    " (" + dddata->GetData(eleccion2)->GetExforEntryCode() + ")";
			gra->SetTitle(titulo); 
			gra->Draw("ALP");
			gra->GetXaxis()->SetTitle("T (MeV)");
			gra->GetYaxis()->SetTitle("#frac{d#sigma}{dT d#theta}  (mb/MeV/rad)");
			gra->Draw("ALP");
			aCanvas->Update();
			showdata_dd_app->Run(kTRUE);
			delete aCanvas;
		    }
		}
	    }
	    while (eleccion2 != -1);	    
	}
    } 
    while (eleccion != -1);
    
    delete dddata;
    return;
}


void PrintUsage(char * progname)
{
  
    std::cout << "Usage: " << progname << " [-g] [-s] <input filename>\n"
	      << "        -g : plot data\n"
	      << "        -s : print a summary of available data\n";
    return;
}


int main(int argc, char **argv)
{
    TROOT simple("showdata","Show data");

    TString filename;
    bool doplots = false;
    bool summary = false;

    if (argc == 1 || argc > 3) 
    {
	PrintUsage(argv[0]);
	return 1;
    } 
    else if (argc == 3)
    {
	if (TString(argv[1]) == TString("-g")) 
	{
	    doplots = true;
	    filename = argv[2];
	}
	else if (TString(argv[1]) == TString("-s"))
	{
	    summary = true;
	    filename = argv[2];
	}
	else 
	{
	    std::cout << "Syntax Error: " << argv[1] << '\n';
	    PrintUsage(argv[0]);
	    return 2;
	}
    }
    else 
    {
	filename = argv[1];
    }
  
    // Open the data file 
    TFile * file = new TFile(filename);
    if (!file) 
    {
	cerr << "Can not open " << filename << '\n';
	return 2;
    }
  
    TString basetree(filename);
    basetree.Replace(basetree.Index(".root",5,0,TString::kExact),5,"",0);
    Ssiz_t s = 0;
    do 
    {
	s = basetree.Index("/",1,0,TString::kExact);
	basetree.Replace(0,s+1,"",0);
    } 
    while (s >= 0);
  

    TString tDEname = basetree + "_DE";
    TString tDAname = basetree + "_DA";
    TString tDDname = basetree + "_DD";
  
    TTree * tDE = (TTree*)file->Get(tDEname);
    TTree * tDA = (TTree*)file->Get(tDAname);
    TTree * tDD = (TTree*)file->Get(tDDname);  

    int Argc = 1;
    char * name = "showdata";
    char * *  Argv = &name;
    TApplication showdata_app("showdata_app", &Argc, Argv);

    if (summary)
    {
	Summarize(tDE,tDA,tDD);
    }
    else 
    {
	TString eleccion("");
	do 
	{
	    std::cout << "\n\n\n\n\n\n\n";
	    std::cout << "***************************************" << '\n';
	    if (tDE)
		std::cout << "* DE has " << tDE->GetEntries() << " entries " << '\n';
	    if (tDA)
		std::cout << "* DA has " << tDA->GetEntries() << " entries " << '\n';
	    if (tDD)
		std::cout << "* DD has " << tDD->GetEntries() << " entries " << '\n';
	    std::cout << "***************************************" << '\n';
	    do 
	    {
		std::cout << '\n' << "Chose (DE,DA,DD,exit): ";
		std::cin >> eleccion;
		eleccion.ToUpper();
	    } 
	    while (eleccion != "DE" && eleccion != "DA" && eleccion != "DD" && eleccion != "EXIT");
	    if (tDE && eleccion == "DE" && tDE->GetEntries() > 0) 
	    {
		eleccion = "";
		showdata_de(tDE,&showdata_app,doplots);
	    } 
	    else if (tDA && eleccion == "DA" && tDA->GetEntries() > 0) 
	    {
		eleccion = "";
		showdata_da(tDA,&showdata_app,doplots);
	    } 
	    else if (tDD && eleccion == "DD" && tDD->GetEntries() > 0) 
	    {
		eleccion = "";
		showdata_dd(tDD,&showdata_app,doplots);
	    }
	} 
	while (eleccion != "EXIT");
    }
    
    file->Close();
    return 0;
}
  





