#ifndef G4InclRootLogger_hh
#define G4InclRootLogger_hh 1

#include "G4VInclLogger.hh"
#include "TROOT.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

class G4InclRootLogger : public G4VInclLogger {

public:
  G4InclRootLogger() {
    G4cout <<"Created ROOT based histogramming" << G4endl;
    numberOfHistograms = 0;
  };

  ~G4InclRootLogger() {};

  /**
   * Book 1D histogram.
   */
  void bookHistogram1D(G4String name, G4int bins, G4double xmin, G4double xmax) {
    G4cout <<"Booking histogram " << name << G4endl;
    theHistograms[numberOfHistograms] = new TH1F(name.c_str(), name.c_str(), (int) bins, (double) xmin, (double) xmax);
    numberOfHistograms++;
    G4cout <<"Done." << G4endl;
  }

  /**
   * Book 2D histogram.
   */
  void bookHistogram2D(G4String name, G4int binsx, G4double xmin, G4double xmax, G4int binsy, G4double ymin, G4double ymax) {
    G4cout <<"Booking 2D histogram " << name << G4endl;
    the2DHistograms[numberOf2DHistograms] = new TH2F(name.c_str(), name.c_str(), (int) binsx, (double) xmin, (double) xmax,
						     (int) binsy, (double) ymin, (double) ymax);
    numberOf2DHistograms++;
    G4cout <<"Done." << G4endl;
  }

  /**
   * Fill 1D histogram.
   */
  void fillHistogram1D(G4String name, G4double value) {
    for(int i = 0; i < numberOfHistograms; i++) {
      if(TString(theHistograms[i]->GetName()) == TString(name.c_str())) {
	theHistograms[i]->Fill((double) value);
      }
    }
  }

  /**
   * Fill 2D histogram.
   */
  void fillHistogram2D(G4String name, G4double xvalue, G4double yvalue) {
    for(int i = 0; i < numberOf2DHistograms; i++) {
      if(TString(the2DHistograms[i]->GetName()) == TString(name.c_str())) {
	the2DHistograms[i]->Fill((double) xvalue, (double) yvalue);
      }
    }
  }
  
  /**
   * Write the histograms to the ROOT file.
   */
  void saveHistograms() {
    for(int i = 0; i < numberOfHistograms; i++) {
      G4cout <<"Writing histogram: ";
      cout << theHistograms[i]->GetName() << endl;
      theHistograms[i]->SetLineColor(kRed); // Color 632 = red
      theHistograms[i]->Write();
    }

    for(int i = 0; i < numberOf2DHistograms; i++) {
      G4cout <<"Writing histogram: ";
      cout << the2DHistograms[i]->GetName() << endl;
      the2DHistograms[i]->Write();
    }
  }

private:
  TFile *theFile;
  TH1F *theHistograms[1000];
  TH2F *the2DHistograms[1000];
  G4int numberOfHistograms;
  G4int numberOf2DHistograms;
};

#endif
