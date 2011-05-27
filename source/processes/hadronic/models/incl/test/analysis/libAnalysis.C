//
//  SATURNE experiment (n double diff cross at various angles)
//
//  A.B. 8/10/2007
//

#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TPaveLabel.h"
#include "TList.h"

#include <vector>
#include <string>
#include <algorithm>

// Wrapper class around the normal histogramming functionality of
// ROOT. This class attempts to make creating often needed histogram
// types (such, as the ones using logarithmic binning) easier and to
// hide some ugly ROOT details.
class HistoFactory {
private:
  Long_t histoID; // Unique ID number for histograms

public:
  HistoFactory() { histoID = 0; };
  ~HistoFactory() {};

  // Create a 1D histogram with desired binning.
  TH1D* createHisto1D(char *name, char *title,
			    Int_t nbins, Double_t xmin, Double_t xmax,
			    Bool_t logX = false) {
    TH1D *histo = 0;
    // First, handle some possible error conditions:
    if(xmin >= xmax) {
      cout <<"HistoFactory error: xmin is greater than xmax!" << endl;
      return 0; // Return empty pointer
    }
    if(logX && (xmin <= 0.0 || xmax <= 0.0)) {
      cout <<"HistoFactory error: Defined limits (xmin = " << xmin << ", xmax = "<< xmax << ")" << endl;
      cout <<"can't be used in a plot with logarithmic x axis!" << endl;
      cout <<"Switching to linear plot mode!" << endl;
      logX = false;
    }

    if(logX) { // Logarithmic binning
      Double_t xbins[nbins+1];
      Double_t fact=(log(xmax)-log(xmin))/nbins;
      for(Int_t i = 0; i < nbins+1; ++i) {
	xbins[i] = exp(log(xmin)+fact*i);
      }
      histo = new TH1D(name, title, nbins, xbins);
    } else { // Linear binning
      histo = new TH1D(name, title, nbins, xmin, xmax);
    }

    histoID++; // Increment the histogram ID counter

    // Common histogram style settings:
    histo->SetLineWidth(1.5);

    return histo;
  };

  // Convenience method: no need to always explicitly give a name and
  // a title to the histogram.  This method uses the automatic histoID
  // variable (a long integer that is always incremented when a
  // histogram is created).
  TH1D* createHisto1D(Int_t nbins, Double_t xmin, Double_t xmax, Bool_t logX = false) {
    return createHisto1D(Form("histo_%i", histoID), // Create name histo_number where number is 0, 1, 2, etc.
			 Form("Histogram %i", histoID), // Create title "Histogram number"
			 nbins, xmin, xmax, logX);
  };
};

TCut alpha() {
  return TCut("Avv==4&&Zvv==2");
}

TCut proton() {
  return TCut("Avv==1 && Zvv == 1");
}

TCut neutron() {
  return TCut("Avv==1 && Zvv == 0");
}

TCut thetaCut(Double_t theta, Double_t acceptance = 2.0) {
  Double_t minTheta = theta - acceptance;
  Double_t maxTheta = theta + acceptance;
  return TCut(Form("Tetlab > %f && Tetlab < %f", minTheta, maxTheta));
}

// Available default normalizations:
// - NoNormalization: not normalized at all (N)
// - CrossSectionNormalization: produce histograms in units of cross section (d\sigma/dE d\theta)
// - EventNormalization: Normalize by using the number of events (N/N_0)
enum Normalization {NoNormalization, CrossSectionNormalization, EventNormalization};

// Usage:
// HistoFactory *hf = new HistoFactory();
// CalculationAnalysis *analysis = new CalculationAnalysis(hf, 
//                                     "calculationFile.root", "ntupleName");
// analysis->fillAndNormalize("energy", true, "A==1&&Z==0&&theta>5.0&&theta<7.0")->Draw();
class CalculationAnalysis {
public:
  CalculationAnalysis(HistoFactory *aHistoFactory, char *calculationFileName,
		      char *ntupleName,
		      Double_t crossSect, Long_t events) {
    histoFactory = aHistoFactory;
    calculationFile = new TFile(calculationFileName);
    calculationTree = (TTree *) calculationFile->Get(ntupleName);
    crossSection = crossSect;
    numberOfEvents = events;
    normalizationMethod = CrossSectionNormalization;
  };
  ~CalculationAnalysis() {};

  void setNormalization(Normalization newNormalization) {
    normalizationMethod = newNormalization;
  }

  TH1D* fillHisto(TString variable, Double_t nbins,
		  Double_t xmin, Double_t xmax, Bool_t logX = false,
		  TCut cut = "", TString options = "") {
    TH1D *histo = histoFactory->createHisto1D(nbins, xmin, xmax, logX);
    if(histo != 0) {
      calculationTree->Project(histo->GetName(), // This is why histo names MUST be unique!!!
			       variable, cut, options);
      histo->GetXaxis()->SetTitle(variable);
      histo->GetYaxis()->SetTitle("Counts");
      histo->SetTitle("Variable: " + variable + " Cut: " + cut.GetTitle());
    } else {
      cout <<"CalculationAnalysis error: HistoFactory refused to create histogram!" << endl;
    }
    return histo;
  };

  // Normalize the histogram using a specified normalization
  // factor. The bin width is taken into account automatically.
  void normalizeHisto(TH1D *histo, Double_t normalizationFactor = 1.0) {
    for(Int_t bin = 0; bin <= histo->GetNbinsX(); ++bin) {
      Double_t content = histo->GetBinContent(bin);
      Double_t width = histo->GetBinWidth(bin);
      histo->SetBinContent(bin, (content/width) * normalizationFactor);
    }
  }

  // Produce neutron double-differential plot for 
  TH1D* fillAlphaDoubleDiffXS(Double_t theta, Double_t acceptance, Bool_t logE,
				Double_t emin, Double_t emax) {
    Double_t thetaRad = theta * TMath::DegToRad();
    Double_t acceptanceRad = acceptance * TMath::DegToRad();
    Double_t thetaMinRad = thetaRad - acceptanceRad;
    Double_t thetaMaxRad = thetaRad + acceptanceRad;
    Double_t angleTerm = TMath::Cos(thetaMinRad) - TMath::Cos(thetaMaxRad);
    Double_t normalization = 0.0;

    TCut theCut(alpha() + thetaCut(theta, acceptance));
    cout <<"The cut is: " << theCut << endl;
    TH1D *result = fillHisto("Enerj", 100, emin, emax, logE, theCut);

    if(normalizationMethod == CrossSectionNormalization) {
      //      normalization = (crossSection / (2.0 * numberOfEvents * TMath::Pi())) / angleTerm;
      normalization = crossSection / (2.0 * numberOfEvents * TMath::Pi() * angleTerm);
      result->SetTitle(Form("Alpha double-differential cross-section for angle %.1f", theta));
      result->GetXaxis()->SetTitle("Alpha energy [MeV]");
      result->GetYaxis()->SetTitle("#frac{d #sigma}{dE d #theta} [mb/(sr MeV)]");
    } else if(normalizationMethod == EventNormalization) {
      normalization = 1.0/(1000.0*2.0 * TMath::Pi() * numberOfEvents)/angleTerm;
      result->SetTitle(Form("Alpha double-differential cross-section for angle %.1f", theta));
      result->GetXaxis()->SetTitle("Alpha energy [MeV]");
      result->GetYaxis()->SetTitle("N/N0 [1/(msr MeV)]");
    }

    if(normalizationMethod != NoNormalization) {
      normalizeHisto(result, normalization);
    }

    return result;
  }

  TH1D* fillNeutronDoubleDiffXS(Double_t theta, Double_t acceptance, Bool_t logE,
				Double_t emin, Double_t emax) {
    Double_t thetaRad = theta * TMath::DegToRad();
    Double_t acceptanceRad = acceptance * TMath::DegToRad();
    Double_t thetaMinRad = thetaRad - acceptanceRad;
    Double_t thetaMaxRad = thetaRad + acceptanceRad;
    Double_t angleTerm = TMath::Cos(thetaMinRad) - TMath::Cos(thetaMaxRad);
    Double_t normalization = 0.0;

    TCut theCut(neutron() + thetaCut(theta, acceptance));
    cout <<"The cut is: " << theCut << endl;
    TH1D *result = fillHisto("Enerj", 100, emin, emax, logE, theCut);

    if(normalizationMethod == CrossSectionNormalization) {
      //      normalization = (crossSection / (2.0 * numberOfEvents * TMath::Pi())) / angleTerm;
      normalization = crossSection / (2.0 * numberOfEvents * TMath::Pi() * angleTerm);
      result->SetTitle(Form("Neutron double-differential cross-section for angle %.1f", theta));
      result->GetXaxis()->SetTitle("Neutron energy [MeV]");
      result->GetYaxis()->SetTitle("#frac{d #sigma}{dE d #theta} [mb/(sr MeV)]");
    } else if(normalizationMethod == EventNormalization) {
      normalization = 1.0/(1000.0*2.0 * TMath::Pi() * numberOfEvents)/angleTerm;
      result->SetTitle(Form("Neutron double-differential cross-section for angle %.1f", theta));
      result->GetXaxis()->SetTitle("Neutron energy [MeV]");
      result->GetYaxis()->SetTitle("N/N0 [1/(msr MeV)]");
    }

    if(normalizationMethod != NoNormalization) {
      normalizeHisto(result, normalization);
    }

    return result;
  }

  TH1D* fillCascadeNeutronDoubleDiffXS(Double_t theta, Double_t acceptance, Bool_t logE,
				Double_t emin, Double_t emax) {
    Double_t thetaRad = theta * TMath::DegToRad();
    Double_t acceptanceRad = acceptance * TMath::DegToRad();
    Double_t thetaMinRad = thetaRad - acceptanceRad;
    Double_t thetaMaxRad = thetaRad + acceptanceRad;
    Double_t angleTerm = TMath::Cos(thetaMinRad) - TMath::Cos(thetaMaxRad);
    Double_t normalization = 0.0;

    TCut theCut(neutron() + thetaCut(theta, acceptance) + " Ityp == 1");
    cout <<"The cut is: " << theCut << endl;
    TH1D *result = fillHisto("Enerj", 100, emin, emax, logE, theCut);

    if(normalizationMethod == CrossSectionNormalization) {
      //      normalization = (crossSection / (2.0 * numberOfEvents * TMath::Pi())) / angleTerm;
      normalization = crossSection / (2.0 * numberOfEvents * TMath::Pi() * angleTerm);
      result->SetTitle(Form("Cascade neutron double-differential cross-section for angle %.1f", theta));
      result->GetXaxis()->SetTitle("Neutron energy [MeV]");
      result->GetYaxis()->SetTitle("#frac{d #sigma}{dE d #theta} [mb/(sr MeV)]");
    } else if(normalizationMethod == EventNormalization) {
      normalization = 1.0/(1000.0*2.0 * TMath::Pi() * numberOfEvents)/angleTerm;
      result->SetTitle(Form("Cascade neutron double-differential cross-section for angle %.1f", theta));
      result->GetXaxis()->SetTitle("Neutron energy [MeV]");
      result->GetYaxis()->SetTitle("N/N0 [1/(msr MeV)]");
    }

    if(normalizationMethod != NoNormalization) {
      normalizeHisto(result, normalization);
    }

    return result;
  }

  TH1D* fillEvaporationNeutronDoubleDiffXS(Double_t theta, Double_t acceptance, Bool_t logE,
				Double_t emin, Double_t emax) {
    Double_t thetaRad = theta * TMath::DegToRad();
    Double_t acceptanceRad = acceptance * TMath::DegToRad();
    Double_t thetaMinRad = thetaRad - acceptanceRad;
    Double_t thetaMaxRad = thetaRad + acceptanceRad;
    Double_t angleTerm = TMath::Cos(thetaMinRad) - TMath::Cos(thetaMaxRad);
    Double_t normalization = 0.0;

    TCut theCut(neutron() + thetaCut(theta, acceptance) + " Ityp == 0");
    cout <<"The cut is: " << theCut << endl;
    TH1D *result = fillHisto("Enerj", 100, emin, emax, logE, theCut);

    if(normalizationMethod == CrossSectionNormalization) {
      //      normalization = (crossSection / (2.0 * numberOfEvents * TMath::Pi())) / angleTerm;
      normalization = crossSection / (2.0 * numberOfEvents * TMath::Pi() * angleTerm);
      result->SetTitle(Form("Evaporation neutron double-differential cross-section for angle %.1f", theta));
      result->GetXaxis()->SetTitle("Neutron energy [MeV]");
      result->GetYaxis()->SetTitle("#frac{d #sigma}{dE d #theta} [mb/(sr MeV)]");
    } else if(normalizationMethod == EventNormalization) {
      normalization = 1.0/(1000.0*2.0 * TMath::Pi() * numberOfEvents)/angleTerm;
      result->SetTitle(Form("Evaporation neutron double-differential cross-section for angle %.1f", theta));
      result->GetXaxis()->SetTitle("Neutron energy [MeV]");
      result->GetYaxis()->SetTitle("N/N0 [1/(msr MeV)]");
    }

    if(normalizationMethod != NoNormalization) {
      normalizeHisto(result, normalization);
    }

    return result;
  }

  TH1D* fillCascadeSpectatorNeutronDoubleDiffXS(Double_t theta, Double_t acceptance, Bool_t logE,
				Double_t emin, Double_t emax) {
    Double_t thetaRad = theta * TMath::DegToRad();
    Double_t acceptanceRad = acceptance * TMath::DegToRad();
    Double_t thetaMinRad = thetaRad - acceptanceRad;
    Double_t thetaMaxRad = thetaRad + acceptanceRad;
    Double_t angleTerm = TMath::Cos(thetaMinRad) - TMath::Cos(thetaMaxRad);
    Double_t normalization = 0.0;

    TCut theCut(neutron() + thetaCut(theta, acceptance) + " Ityp == -1");
    cout <<"The cut is: " << theCut << endl;
    TH1D *result = fillHisto("Enerj", 100, emin, emax, logE, theCut);

    if(normalizationMethod == CrossSectionNormalization) {
      //      normalization = (crossSection / (2.0 * numberOfEvents * TMath::Pi())) / angleTerm;
      normalization = crossSection / (2.0 * numberOfEvents * TMath::Pi() * angleTerm);
      result->SetTitle(Form("Spectator (cascade) neutron double-differential cross-section for angle %.1f", theta));
      result->GetXaxis()->SetTitle("Neutron energy [MeV]");
      result->GetYaxis()->SetTitle("#frac{d #sigma}{dE d #theta} [mb/(sr MeV)]");
    } else if(normalizationMethod == EventNormalization) {
      normalization = 1.0/(1000.0*2.0 * TMath::Pi() * numberOfEvents)/angleTerm;
      result->SetTitle(Form("Spectator (cascade) neutron double-differential cross-section for angle %.1f", theta));
      result->GetXaxis()->SetTitle("Neutron energy [MeV]");
      result->GetYaxis()->SetTitle("N/N0 [1/(msr MeV)]");
    }

    if(normalizationMethod != NoNormalization) {
      normalizeHisto(result, normalization);
    }

    return result;
  }

  TH1D* fillSpectatorDeExcitationNeutronDoubleDiffXS(Double_t theta, Double_t acceptance, Bool_t logE,
				Double_t emin, Double_t emax) {
    Double_t thetaRad = theta * TMath::DegToRad();
    Double_t acceptanceRad = acceptance * TMath::DegToRad();
    Double_t thetaMinRad = thetaRad - acceptanceRad;
    Double_t thetaMaxRad = thetaRad + acceptanceRad;
    Double_t angleTerm = TMath::Cos(thetaMinRad) - TMath::Cos(thetaMaxRad);
    Double_t normalization = 0.0;

    TCut theCut(neutron() + thetaCut(theta, acceptance) + " Ityp == -2");
    cout <<"The cut is: " << theCut << endl;
    TH1D *result = fillHisto("Enerj", 100, emin, emax, logE, theCut);

    if(normalizationMethod == CrossSectionNormalization) {
      //      normalization = (crossSection / (2.0 * numberOfEvents * TMath::Pi())) / angleTerm;
      normalization = crossSection / (2.0 * numberOfEvents * TMath::Pi() * angleTerm);
      result->SetTitle(Form("Spectator de-excitation neutron double-differential cross-section for angle %.1f", theta));
      result->GetXaxis()->SetTitle("Neutron energy [MeV]");
      result->GetYaxis()->SetTitle("#frac{d #sigma}{dE d #theta} [mb/(sr MeV)]");
    } else if(normalizationMethod == EventNormalization) {
      normalization = 1.0/(1000.0*2.0 * TMath::Pi() * numberOfEvents)/angleTerm;
      result->SetTitle(Form("Spectator de-excitation neutron double-differential cross-section for angle %.1f", theta));
      result->GetXaxis()->SetTitle("Neutron energy [MeV]");
      result->GetYaxis()->SetTitle("N/N0 [1/(msr MeV)]");
    }

    if(normalizationMethod != NoNormalization) {
      normalizeHisto(result, normalization);
    }

    return result;
  }

  TH1D* fillProtonDoubleDiffXS(Double_t theta, Double_t acceptance, Bool_t logE,
				Double_t emin, Double_t emax) {
    Double_t thetaRad = theta * TMath::DegToRad();
    Double_t acceptanceRad = acceptance * TMath::DegToRad();
    Double_t thetaMinRad = thetaRad - acceptanceRad;
    Double_t thetaMaxRad = thetaRad + acceptanceRad;
    Double_t angleTerm = TMath::Cos(thetaMinRad) - TMath::Cos(thetaMaxRad);
    Double_t normalization = 0.0;

    TCut theCut(proton() + thetaCut(theta, acceptance));
    cout <<"The cut is: " << theCut << endl;
    TH1D *result = fillHisto("Enerj", 100, emin, emax, logE, theCut);

    if(normalizationMethod == CrossSectionNormalization) {
      normalization = (crossSection / (2.0 * numberOfEvents * TMath::Pi())) / angleTerm;
      result->SetTitle(Form("Proton double-differential N/N0 for angle %.1f", theta));
      result->GetXaxis()->SetTitle("Proton energy [MeV]");
      result->GetYaxis()->SetTitle("#frac{d #sigma}{dE d #theta} [mb/ (sr MeV)");
    } else if(normalizationMethod == EventNormalization) {
      normalization = 1.0/(1000.0*2.0 * TMath::Pi() * numberOfEvents)/angleTerm;
      result->SetTitle(Form("Proton double-differential N/N0 for angle %.1f", theta));
      result->GetXaxis()->SetTitle("Proton energy [MeV]");
      result->GetYaxis()->SetTitle("N/N0 [1/(msr MeV)]");
    }

    if(normalizationMethod != NoNormalization) {
      normalizeHisto(result, normalization);
    }

    return result;
  }

  TH1D* fillIsotopeDistribution(Int_t Z, Int_t Zmin, Int_t Zmax) {
    Int_t numberOfBins = Zmax - Zmin + 1;
    Double_t xmin = Zmin - 0.5;
    Double_t xmax = Zmax + 0.5;
    TCut theCut(Form("Z==%i", Z));
    TH1D *result = fillHisto("A", numberOfBins, xmin, xmax, false, theCut);
    Double_t normalizationFactor = crossSection/(2.0 * numberOfEvents * TMath::Pi());
    normalizeHisto(result, normalizationFactor);
    result->SetTitle(Form("Isotopes of element %i", Z));
    result->GetXaxis()->SetTitle("A");
    result->GetYaxis()->SetTitle("#sigma (mb)");
    return result;
  }

private:
  HistoFactory *histoFactory;
  TFile *calculationFile;
  TTree *calculationTree;
  Double_t crossSection;
  Long_t numberOfEvents;
  Normalization normalizationMethod;
};

class ViewManager {
public:
  ViewManager() {
    theCanvas = new TCanvas();
    divX = 1;
    divY = 1;
    currentX = 1;
    currentY = 1;
    currentPad = 1;
    theCanvas->cd(currentPad);
  };
  ViewManager(Int_t zonesX, Int_t zonesY) {
    theCanvas = new TCanvas();
    divX = zonesX;
    divY = zonesY;
    theCanvas->Divide(divX, divY);
    currentX = 1;
    currentY = 1;
    currentPad = 1;
    theCanvas->cd(currentPad);
  };

  ~ViewManager() {};

  void gotoNextPad() {
    cout <<"currentX = " << currentX << " currentY = " << currentY << " currentPad = " << currentPad << endl;
    if(currentX < divX && currentY <= divY ) {
      currentX++;
      currentPad++;
      cout <<"Next pad: one step to the right..." << endl;
    } else if(currentX == divX && currentY <  divY) {
      currentX = 1;
      currentY++;
      currentPad++;
      cout <<"Next pad: next row of pads..." << endl;
    } else if(currentX == divX && currentY == divY) {
      currentX = 1;
      currentY = 1;
      currentPad = 1;
      cout <<"Next pad: going back to pad 1" << endl;
    } else {
      cout <<"Error! Unknown pad " << currentPad << " : x = " << currentX << " y = " << currentY << endl;
    }
    theCanvas->cd(currentPad);
  };

  void gotoPreviousPad() {
    cout <<"currentX = " << currentX << " currentY = " << currentY << " currentPad = " << currentPad << endl;
    if(currentX > 1 && currentY >= 1) {
      currentX--;
      currentPad--;
      cout <<"Next pad: one step to the right..." << endl;
    } else if(currentX == 1 && currentY > 1) {
      currentX = divX;
      currentY--;
      currentPad--;
      cout <<"Next pad: next row of pads..." << endl;
    } else if(currentX == 1 && currentY == 1) {
      currentX = divX;
      currentY = divY;
      currentPad = divX*divY; // The last pad
      cout <<"Next pad: going back to pad: " << currentPad << endl;
    } else {
      cout <<"Error! Unknown pad " << currentPad << " : x = " << currentX << " y = " << currentY << endl;
    }
    theCanvas->cd(currentPad);
  }

  // This plotting function has PAW-like behaviour:
  // If the canvas has been divided, a new histogram is plotted to the
  // next zone in order to keep the old histo visible.
  void plot(TH1D *histo, Bool_t addToSamePlot = false) {
    if(addToSamePlot) {
      gotoPreviousPad();
      histo->SetLineStyle(kDashed);
      histo->Draw("same");
    } else {
      histo->Draw();
    }
    theCanvas->Modified(true);
    theCanvas->Update();
    gotoNextPad();
  };

  void plot(TGraph *graph, Bool_t addToSamePlot = false) {
    if(addToSamePlot) {
      gotoPreviousPad();
      graph->SetLineStyle(kDashed);
      graph->Draw("p*, same");
    } else {
      graph->Draw("a, p*");
    }
    theCanvas->Modified(true);
    theCanvas->Update();
    gotoNextPad();
  };

  // Plot a list of histograms. If plottingNormalization is true, we
  // scale the first histo by 1.0, the second by 0.1, the third by
  // 0.01 and so on.
  void plot(TList *histos, Bool_t plottingNormalization = false) {
    Bool_t isFirst = true;
    TListIter *iterator = new TListIter(histos);
    TObject *histo = 0;
    Double_t scaling = 1.0;
    while(histo = iterator->Next()) {
      if(histo == 0) break;
      if(isFirst) {
	histo->DrawClone();
	isFirst = false;
      } else {
	TH1D *tmpHisto = (TH1D *) histo;
	if(plottingNormalization) {
	  scaling = scaling * 0.1;
	  tmpHisto->Scale(scaling);
	}
	tmpHisto->Draw("same");
      }
    }
    theCanvas->Modified(true);
    theCanvas->Update();
    gotoNextPad();
  }

private:
  TCanvas *theCanvas;
  Int_t divX, divY;
  Int_t currentX, currentY, currentPad;
};

class IDataReader {
public:
  IDataReader () {};
  IDataReader(char *filename) {};
  ~IDataReader() {};

  virtual TGraph* getGraph() = 0;
};

class TwoColumnBlockReader : public IDataReader {
public:
  TwoColumnBlockReader () {};
  TwoColumnBlockReader(char *filename) {
    in = new ifstream(filename, ios::in);
    theGraph =  0;
  }

  ~TwoColumnBlockReader() {}

  TGraph *getGraph() {
    char buffer[256];
    if(theGraph != 0) return theGraph;

    while(1) {
      in->getline(buffer, 256);
      if(in->eof()) break; // The proper end of file
      cout << buffer << endl;
      TString line(buffer);
      if(line.CompareTo("DATA") == 0) {
	cout <<"DATA block found!" << endl;
	parseDataBlock();
      }
    }

    Int_t numberOfPoints = x.size();
    theGraph = new TGraph(numberOfPoints);
    for(int i = 0; i < numberOfPoints; ++i) {
      theGraph->SetPoint(i, x[i], y[i]);
    }
    return theGraph;
  }

  void parseDataBlock() {
    char buffer[256];
    while(1) {
      in->getline(buffer, 256);
      if(in->eof() || in->fail() || in->bad()) {
	cout <<"Error! Error when parsing data block!" << endl;
	break;
      }
      TString data(buffer);
      if(data.CompareTo("END") == 0) {
	cout <<"End of DATA block" << endl;
	break;
      }
      const std::string realData(buffer);
      splitBufferToXAndYValues(realData);
    }
  }

  void splitBufferToXAndYValues(const std::string &data) {
    const char delim = ' ';
    vector<string> output = split(data, delim);
    double newX = ::atof((output[0]).c_str());
    double newY = ::atof((output[1]).c_str());
    cout <<"x = " << newX << " y = " << newY << endl;
    x.push_back(newX);
    y.push_back(newY);
  }

  vector<string> split( const string &str, const char delim ) {
    typedef string::const_iterator iter;
    iter beg = str.begin();
    vector<string> tokens;

    while(beg != str.end()) {
      //cout << ":" << beg._Myptr << ":" << endl;
      iter temp = find(beg, str.end(), delim);
      if(beg != str.end())
	tokens.push_back(string(beg, temp));
      beg = temp;
      while ((beg != str.end()) && (*beg == delim))
	beg++;
    }

    return tokens;
  }

private:
  ifstream *in;
  std::vector<double> x;
  std::vector<double> y;
  TList *xvalues;
  TList *yvalues;
  TGraph *theGraph;
};


class RomanoReader : public IDataReader {
public:
  RomanoReader () {};
  RomanoReader(char *filename) {
    in = new ifstream(filename, ios::in);
    theGraph =  0;
  }

  ~RomanoReader() {}

  TGraph *getGraph() {
    char buffer[256];
    if(theGraph != 0) return theGraph;

    while(1) {
      in->getline(buffer, 256);
      if(in->eof()) break; // The proper end of file
      cout << buffer << endl;
      TString line(buffer);
      if(line.CompareTo("DATA") == 0) {
	cout <<"DATA block found!" << endl;
	parseDataBlock();
      }
    }

    Int_t numberOfPoints = x.size();
    theGraph = new TGraph(numberOfPoints);
    for(int i = 0; i < numberOfPoints; ++i) {
      theGraph->SetPoint(i, x[i], y[i]);
    }
    return theGraph;
  }

  void parseDataBlock() {
    char buffer[256];
    while(1) {
      in->getline(buffer, 256);
      if(in->eof() || in->fail() || in->bad()) {
	cout <<"Error! Error when parsing data block!" << endl;
	break;
      }
      TString data(buffer);
      if(data.CompareTo("END") == 0) {
	cout <<"End of DATA block" << endl;
	break;
      }
      const std::string realData(buffer);
      splitBufferToXAndYValues(realData);
    }
  }

  void splitBufferToXAndYValues(const std::string &data) {
    const char delim = ' ';
    vector<string> output = split(data, delim);
    double newX = ::atof((output[1]).c_str());
    double newY = ::atof((output[2]).c_str());
    cout <<"x = " << newX << " y = " << newY << endl;
    x.push_back(newX);
    y.push_back(newY);
  }

  vector<string> split( const string &str, const char delim ) {
    typedef string::const_iterator iter;
    iter beg = str.begin();
    vector<string> tokens;

    while(beg != str.end()) {
      //cout << ":" << beg._Myptr << ":" << endl;
      iter temp = find(beg, str.end(), delim);
      if(beg != str.end())
	tokens.push_back(string(beg, temp));
      beg = temp;
      while ((beg != str.end()) && (*beg == delim))
	beg++;
    }

    return tokens;
  }

private:
  ifstream *in;
  std::vector<double> x;
  std::vector<double> y;
  TList *xvalues;
  TList *yvalues;
  TGraph *theGraph;
};

class ThreeColumnReader : public IDataReader {
public:
  ThreeColumnReader() {};
  ThreeColumnReader(char *filename) {
    thefile = new std::string(filename);
  };

  TGraph* getGraph() {
    ifstream in(thefile->c_str());
    while(1) {
      Double_t angle, newx, newy;
      in >> angle >> newx >> newy;
      if(!in.good()) break;
      cout <<"x = " << newx << " y = " << newy << endl;
      x.push_back(newx);
      y.push_back(newy);
    }
    TGraph *theGraph = new TGraph(x.size());
    for(int i = 0; i != x.size(); ++i) {
      theGraph->SetPoint(i, x[i], y[i]);
    }
    return theGraph;
  }

private:
  std::string *thefile;
  std::vector<double> x;
  std::vector<double> y;
};

void plotExpTheta(Char_t* file_name, Char_t* racine, Float_t fnor)
{
// Experimental points:
	Char_t toto[256];
	Int_t e;
	Float_t e_n[100];
	Float_t sig[100];
	Float_t dsig[100];
	Float_t de_n[100];
	Int_t c,iexp=0;
//Concatenation Racine+File_Name
	e=sprintf(toto,"%s%s",racine,file_name);
	cout<<toto<<endl;

	FILE* fexp160 =  fopen(toto,"r");

	iexp=0;
	do {
		c=fscanf(fexp160,"%f %f %f %f ",&e_n[iexp],&de_n[iexp],&sig[iexp],&dsig[iexp]);
		if(c==EOF) break;
		sig[iexp]=fnor*sig[iexp];
		dsig[iexp]=fnor*dsig[iexp];
//	cout<<e_n[iexp]<<"  "<<sig[iexp]<<"  "<<dsig[iexp]<<endl;
		iexp=iexp+1;
		if(iexp>100) cout<<"More experimental values than the dimension"<<endl;
	}
	while(1);

	TGraphErrors* grexp160= new TGraphErrors(iexp,e_n,sig,de_n,dsig);
	grexp160->SetMarkerColor(4);
	grexp160->SetMarkerStyle(21);
	grexp160->SetMarkerSize(0.4);

	grexp160->Draw("PZ");
	return;
}
void plotTheorTheta(TTree* ref,Char_t* titre,TCanvas* c1,Double_t emin,Double_t emax,Int_t logE,TCut neutron,TCut thet0,Double_t fnor0,Double_t fnora,Char_t* angle,Int_t first)
{
  TH1F *hist0 = 0;
  TH1F *hist_de = 0;
  Double_t de = 0.0;
	fnor0=fnor0*fnora;
	if (logE) {
// Logarithmic plot in X:
		c1->SetLogx();
		Float_t xbins[99+1];
		Float_t fact=(log(emax)-log(emin))/99;
		for(Int_t i=0;i<99+1;++i) {xbins[i]=exp(log(emin)+fact*i);}
		hist0=new TH1F("hist0","",99,xbins);
		hist_de=new TH1F("hist_de","",99,xbins);
		for(Int_t i=0;i<99;++i) {hist_de->SetBinContent(i,(Double_t)1./(xbins[i+1]-xbins[i]));}
		de=1.0;
	}
	else {
// Linear plot in X:
//		c1->SetLinx();
		TH1F* hist0=new TH1F("hist0","",99,emin,emax);
		de=(emax-emin)/99;
	}
	hist0->SetFillStyle(0);
	hist0->SetStats(kFALSE);
	hist0->SetTitle(titre);
	hist0->GetXaxis()->SetTitle(" Energy (MeV)");
   	hist0->GetXaxis()->CenterTitle(true);
   	hist0->GetXaxis()->SetLabelSize(0.03);
	hist0->GetYaxis()->SetTitle("Cross section (mb)");
   	hist0->GetYaxis()->CenterTitle(true);
   	hist0->GetYaxis()->SetLabelSize(0.03);
   	hist0->GetYaxis()->SetTitleOffset(1.3);
	hist0->SetLineColor(kRed);
	
	
	// Project in the hist histo)
	if(first){ref->Draw("energy >> hist0",neutron && thet0);}
	else {ref->Draw("energy >> hist0",neutron && thet0,"same");}
	
	fnor0=fnor0/de;
	hist0->Scale(fnor0);
	if (logE) {hist0->Multiply(hist0,hist_de,1.,1.);}
	
// Full y extension of the picture:
	hist0->SetMinimum(1.e-1);
	hist0->SetMaximum(1.e+5);
// Legend (angle):
   Float_t y1=200.;
   Float_t y2=400.;
   y1=y1*fnora;
   y2=y2*fnora;	
   TPaveLabel* p1 = new TPaveLabel(5,y1,7,y2,angle);
   p1->SetBorderSize(0);
   p1->SetFillColor(0);
   p1->SetTextSize(1.0);
   p1->Draw("same");
   
	delete hist_de;
	return;
}
