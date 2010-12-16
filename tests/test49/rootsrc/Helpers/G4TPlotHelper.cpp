//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// --------------------------------------------------------------------
// Class implementation
// --------------------------------------------------------------------

#include "G4TPlotHelper.h"


G4TPlotHelper *gPlotHelper = new G4TPlotHelper();


ClassImp(G4TPlotHelper)

using namespace std;
using namespace ROOT;
using namespace TMath;

//______________________________________________________________________________
TCanvas* G4TPlotHelper::PrepareCanvas()
{
	if(fCanvas != 0){
		delete fCanvas;
	}

	cout << "Preparing canvas..." << endl;
	fCanvas = new TCanvas("Canvas","HIGZ Window",800,600);
	fCanvas->Clear();
	fCanvas->Update();
	return fCanvas;
}

//______________________________________________________________________________
void G4TPlotHelper::ClearCanvas()
{
	GetCanvas(); // make sure there's a canvas

	fCanvas->Clear();
	fCanvas->SetBorderMode(0);
	fCanvas->SetBorderSize(0);
	fCanvas->SetFrameBorderSize(0);
	fCanvas->SetFillColor(fBackgroundColor);
	gStyle->SetTextFont(fTextFont);
	fCanvas->SetLeftMargin(0.15);
	fCanvas->SetBottomMargin(0.15);
	//fCanvas->SetTopMargin(0.15);
	fCanvas->Update();
}

//______________________________________________________________________________
Int_t G4TPlotHelper::DivideForNumber(Int_t number)
{
	Int_t xy = sqrt(number);
	if(xy*xy < number)
		xy += 1;

	Int_t x = xy;
	Int_t y = xy;

	if(x * (y-1) >= number)
		y -=1;

	DivideCanvas(x, y);

	return x;
}

//______________________________________________________________________________
TH1F* G4TPlotHelper::PrepareFrame(Double_t hxmin, Double_t hymin, Double_t hxmax, Double_t hymax,TString const& Title)
{
	TH1F* frame = gPad->DrawFrame(hxmin,hymin,hxmax,hymax,Title);

	// Update the title
	TPaveText* localTitle = (TPaveText*)gPad->FindObject("title");
	if(localTitle != 0)
	{
		localTitle->SetFillColor(fBackgroundColor);
		localTitle->SetTextFont(fTextFont);
		localTitle->SetTextSize(fTextSize);
	}

	// Set Titles and Size on the both axis

	TAxis* XAxis = frame->GetXaxis();
	TAxis* YAxis = frame->GetYaxis();

	YAxis->SetTitleOffset(0.8);

	XAxis->SetTitleFont(fTextFont);
	XAxis->SetLabelFont(fTextFont);
	XAxis->SetTitleSize(fAxisXTitleSize);
	XAxis->SetLabelSize(fAxisXTextSize);

	YAxis->SetTitleFont(fTextFont);
	YAxis->SetLabelFont(fTextFont);
	YAxis->SetTitleSize(fAxisYTitleSize);
	YAxis->SetLabelSize(fAxisYTextSize);

	//frame->SetNdivisions(902,"X");

	return frame;
}

//______________________________________________________________________________
void G4TPlotHelper::PreparePad(Int_t number)
{
	GetCanvas(); // make sure there's a canvas
	fCanvas->cd(number);

	// current pad options
	gPad->SetLogy();
	gPad->SetBorderSize(0);
	gPad->SetFrameFillColor(fBackgroundColor);
	gPad->SetFillColor(fBackgroundColor);
	gPad->SetFrameLineWidth(2);

}



//______________________________________________________________________________
void G4TPlotHelper::DrawAnglesLegend(G4TData* publication)
{
	vector<Double_t> cuts = publication->GetCutValues();
	Double_t y1 = 0.0308552;
	Double_t y2 = y1 + (0.0672482 * cuts.size());

	TLegend* legend = new TLegend(0.0419114, y1, 0.3383, y2);
	legend->SetBorderSize(0);
	legend->SetFillColor(fBackgroundColor);



	for(UInt_t i = 0; i < cuts.size(); ++i)
	{
		Double_t angle  = cuts[i];
		Int_t	 marker = GetMarker(i);
		TString  ltext  = TString::Format("%g", angle);

		TGraph* graph = new TGraph(1);
		graph->SetMarkerStyle((Int_t)marker);
		legend->AddEntry(graph, ltext.Data(),"P");
	}

	legend->Draw();
}

//______________________________________________________________________________
void G4TPlotHelper::DrawModelsLegend(vector<G4TData*>* models)
{
	TLegend* legend 	= new TLegend(0.189594,0.0797521,0.544206,0.435175);
	legend->SetBorderSize(0);
	legend->SetFillColor(fBackgroundColor);

	for(UInt_t i = 0; i < models->size(); ++i)
	{
		TGraph* graph = new TGraph(2);
		graph->SetLineWidth(2);
		graph->SetLineColor((*models)[i]->GetRenderColor());
		legend->AddEntry(graph, (*models)[i]->GetModelName().Data(),"L");
	}
	legend->Draw();
}



//______________________________________________________________________________
void G4TPlotHelper::DrawBigTitle(TString const& text)
{
	fCanvas->cd(0);
	TLatex title;
	title.SetTextAlign(23); // top-middle
	title.SetTextSize( 0.03);
	title.DrawLatex(0.5, 0.99, text.Data());
}

//______________________________________________________________________________
void G4TPlotHelper::DrawRightAxis(TH1F* frame, Double_t hxmax, Double_t hymin, Double_t  hymax)
{
	// Draw the axis on the right side
	TAxis*  yAxis = frame->GetYaxis();
	TGaxis* nAxis = new TGaxis(hxmax, hymin, hxmax , hymax,
			yAxis->GetXmin(), yAxis->GetXmax(), yAxis->GetNdivisions(), "+LG");
	nAxis->SetLabelSize(0); // no labels
	nAxis->Draw();
}

//______________________________________________________________________________
void G4TPlotHelper::DivideCanvas(Int_t xDiv, Int_t yDiv)
{
	GetCanvas(); // make sure there's a canvas
	fCanvas->Divide(xDiv,yDiv,0,0, fBackgroundColor);

}

//______________________________________________________________________________
TString	G4TPlotHelper::PrepareCutEnumeration(G4TData* publication)
{
	vector<Double_t> cuts = publication->GetCutValues();
	TString result = "";
	for(UInt_t i = 0; i< cuts.size(); ++i)
	{
		if(i > 0) result = result + ",";
		Double_t cut = cuts[i];
		result = result.Append(TString::Format("%g", cut).Data());
	}
	return result;
}

//______________________________________________________________________________
void G4TPlotHelper::SetTitlePosition(Int_t align , Float_t x , Float_t y )
{
	gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleAlign(align); // top-middle

	if(x != 0) gStyle->SetTitleX(x);
	if(y != 0) gStyle->SetTitleY(y);
}

//______________________________________________________________________________
TCanvas* G4TPlotHelper::GetCanvas()
{
	if(fCanvas == 0)
		return PrepareCanvas();
	return fCanvas;
}


//______________________________________________________________________________
G4TPlotHelper::~G4TPlotHelper ()
{
  if (gPlotHelper == this)
	  gPlotHelper = 0;
}


//______________________________________________________________________________
vector<Int_t>& G4TPlotHelper::GetMarkers()
{
	return fMarkers;
}

//______________________________________________________________________________
Int_t  G4TPlotHelper::GetMarker(Int_t idx)
{
	return fMarkers[idx];
}

//______________________________________________________________________________
Int_t G4TPlotHelper::GetBackgroundColor() const
{
    return fBackgroundColor;
}

//______________________________________________________________________________
void G4TPlotHelper::SetBackgroundColor(Int_t fBackgroundColor)
{
    this->fBackgroundColor = fBackgroundColor;
}

//______________________________________________________________________________
Int_t G4TPlotHelper::GetTextFont() const
{
    return fTextFont;
}

//______________________________________________________________________________
void G4TPlotHelper::SetTextFont(Int_t fTextFont)
{
    this->fTextFont = fTextFont;
}

//______________________________________________________________________________
Float_t G4TPlotHelper::GetTextSize() const
{
    return fTextSize;
}

//______________________________________________________________________________
void G4TPlotHelper::SetTextSize(Float_t fTextSize)
{
    this->fTextSize = fTextSize;
}

//______________________________________________________________________________
Float_t G4TPlotHelper::GetAxisXTextSize() const
{
    return fAxisXTextSize;
}

//______________________________________________________________________________
void G4TPlotHelper::SetAxisXTextSize(Float_t fAxisXTextSize)
{
    this->fAxisXTextSize = fAxisXTextSize;
}

//______________________________________________________________________________
Float_t G4TPlotHelper::GetAxisYTextSize() const
{
    return fAxisYTextSize;
}

//______________________________________________________________________________
void G4TPlotHelper::SetAxisYTextSize(Float_t fAxisYTextSize)
{
    this->fAxisYTextSize = fAxisYTextSize;
}


//______________________________________________________________________________
Float_t G4TPlotHelper::GetAxisYTitleSize() const
{
    return fAxisYTitleSize;
}

//______________________________________________________________________________
void G4TPlotHelper::SetAxisYTitleSize(Float_t fAxisYTitleSize)
{
    this->fAxisYTitleSize = fAxisYTitleSize;
}


//______________________________________________________________________________
Float_t G4TPlotHelper::GetAxisXTitleSize() const
{
    return fAxisXTitleSize;
}

//______________________________________________________________________________
void G4TPlotHelper::SetAxisXTitleSize(Float_t fAxisXTitleSize)
{
    this->fAxisXTitleSize = fAxisXTitleSize;
}

