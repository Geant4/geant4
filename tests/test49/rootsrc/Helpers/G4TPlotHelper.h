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
// Class G4TPlotHelper
//
// Class description:
//
// This class is a simple wrapper/helper for ROOT plotting commands. It
// also contains different features:
//   - division for n:  DivideForNumber(Int_t number) which divides the
//		  canvas automatically
//   - different 'prepare' methods in order to set the look and feel
//
// History:
// Roman Atachiants, 18/08/2009 - initial version
//
// --------------------------------------------------------------------

#ifndef G4TPlotHelper_H_
#define G4TPlotHelper_H_

#include "../CommonHeaders.h"
#include "../Database/G4TData.h"



class G4TPlotHelper : public TObject {

  private:
    Int_t 			fBackgroundColor;
    Int_t 			fTextFont;
    Float_t 		fTextSize;
    Float_t 		fAxisXTextSize;
    Float_t 		fAxisYTextSize;
    Float_t 		fAxisXTitleSize;
    Float_t 		fAxisYTitleSize;
    TCanvas*		fCanvas;
    vector<Int_t>	fMarkers;

  public:

    TCanvas*		PrepareCanvas();
    TH1F*			PrepareFrame(Double_t hxmin, Double_t hymin, Double_t hxmax, Double_t hymax,TString const& Title);
    TString			PrepareCutEnumeration(G4TData* publication);
    void			PreparePad(Int_t number);
    Int_t			DivideForNumber(Int_t number);
    void			DivideCanvas(Int_t xDiv, Int_t yDiv);
    void			ClearCanvas();
    void			DrawBigTitle(TString const& text);
    void			DrawRightAxis(TH1F* frame, Double_t hxmax,Double_t  hymin, Double_t  hymax);
    void 			DrawAnglesLegend(G4TData* publication);
    void 			DrawModelsLegend(vector<G4TData*>* models);


    // Setters/Getters
    TCanvas*		GetCanvas();
    vector<Int_t>&	GetMarkers();
	Int_t			GetMarker(Int_t idx);
    Int_t 			GetBackgroundColor() const;
    Int_t 			GetTextFont() const;
    Float_t 		GetTextSize() const;
    Float_t 		GetAxisXTextSize() const;
    Float_t 		GetAxisYTextSize() const;
    Float_t 		GetAxisXTitleSize() const;
	Float_t 		GetAxisYTitleSize() const;
    void 			SetBackgroundColor(Int_t fBackgroundColor);
    void 			SetTextFont(Int_t fTextFont);
    void 			SetAxisXTextSize(Float_t);
    void 			SetAxisYTextSize(Float_t);
    void 			SetAxisXTitleSize(Float_t);
    void 			SetAxisYTitleSize(Float_t);
    void 			SetTextSize(Float_t fTextSize);
    void			SetTitlePosition(Int_t align = 23, Float_t x = 0, Float_t y = 0);



    G4TPlotHelper()
    {
        fBackgroundColor = 10;
        fTextFont = 42;
        fTextSize = 0.07;
        fAxisXTextSize = 0.07;
        fAxisYTextSize = 0.07;
        fAxisXTitleSize = 0.07;
        fAxisYTitleSize = 0.08;


        fMarkers.push_back(20);
		fMarkers.push_back(25);
		fMarkers.push_back(22);
		fMarkers.push_back(24);
		fMarkers.push_back(21);
		fMarkers.push_back(26);
		fMarkers.push_back(29);
		fMarkers.push_back(28);
		fMarkers.push_back(23);
    }

    virtual ~G4TPlotHelper();
    ClassDef(G4TPlotHelper, 1)

};

R__EXTERN G4TPlotHelper *gPlotHelper;

#endif




