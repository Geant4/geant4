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
// Igor Pshenichnov 01.12.2006

#include "MFTestHistoManager.hh"

#ifdef G4ANALYSIS_USE_ROOT
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h" 
#endif

MFTestHistoManager::MFTestHistoManager()
  : Z(0), A(0), lowLimitExEn(-1.), upperLimitExEn(-1.), binsExEn(1), eventsPerBin(1) 
{  
  std::cout << "######### Test the Statistical Multifragmentation model of GEANT4" << G4endl;

  while ( (Z<8) || (Z>100) )  {
  std::cout << "Please choose Z from 8 to 100: ";
  std::cin >> Z;
  }

  while ( (A<16) || (A>250) ) {
  std::cout << "Please choose A from 16 to 250: ";
  std::cin >> A;
  }

  while ( (lowLimitExEn<0.) || (lowLimitExEn>100.) ) {
  std::cout << "Please choose the low limit for the range of excitation energy (per nucleon in MeV) : from 0 to 100: ";
  std::cin >> lowLimitExEn;
  }

  while ( (upperLimitExEn<0.) || (upperLimitExEn>100.) || (upperLimitExEn<lowLimitExEn) ) {
  std::cout << "Please choose the upper limit for the range of excitation energy (per nucleon in MeV) : from 0 to 100: ";
  std::cin >> upperLimitExEn;
  }
  lowLimitExEn *= A;
  lowLimitExEn *= MeV;
  upperLimitExEn *=A;
  upperLimitExEn *=MeV;


  while ( (binsExEn<10) || (binsExEn>10000) ) {
  std::cout << "Please choose the number of bins in histograms (10 - 1000): ";
  std::cin >> binsExEn;
  }

  while ( (eventsPerBin<10) || (eventsPerBin>10000) ) {
  std::cout << "Please enter the average number of break-up events in each bin (10 - 10000): ";
  std::cin >> eventsPerBin;
  } 
  iterations = binsExEn*eventsPerBin;

  std::cout << "Please enter the file name to write histograms (.root will be supplied): ";
  std::cin >> fileName;
 
#ifdef G4ANALYSIS_USE_ROOT
  for (G4int j=0; j<20; j++) histo[j] = 0;
  for (G4int l=0; l<10; l++) histo2[l] = 0;
#endif

}


MFTestHistoManager::~MFTestHistoManager()
{}

void MFTestHistoManager::BookHisto()
{
#ifdef G4ANALYSIS_USE_ROOT

  // Open a file to keep histograms inside it ... 

  if ( fileName == "") fileName = "MFTest";
  fileType = "root";

  fileFullName = fileName+"."+fileType;
  compressionFactor = 8;
  fFile = new TFile(fileFullName, "RECREATE", fileName, compressionFactor);
  G4cout << "Histograms will be written to " << fileFullName << G4endl;

  // Book all histograms there ...

  histo[0] =  new TH1D("multip","Average multiplicity",binsExEn,lowLimitExEn/A,upperLimitExEn/A);
  histo[1] =  new TH1D("C10","Yield of C-10",binsExEn,lowLimitExEn/A,upperLimitExEn/A);
  histo[2] =  new TH1D("C11","Yield of C-11",binsExEn,lowLimitExEn/A,upperLimitExEn/A);
  histo[3] =  new TH1D("<Z>(A)","aver Z for given A", A+2, -0.5, A+1.5);
  histo[4] =  new TH1D("NofE(A)","Number of entries for each A",  A+2, -0.5, A+1.5);

  histo2[0] = new TH2D("Z", "Charge distribution of fragments",
		       binsExEn,lowLimitExEn/A,upperLimitExEn/A,
		       Z+2, -0.5, Z+1.5);
  histo2[1] = new TH2D("A", "Mass distribution of fragments",
		       binsExEn,lowLimitExEn/A,upperLimitExEn/A,
		       A+2, -0.5, A+1.5);

#endif  
}

void MFTestHistoManager::NormalizeHisto()
{
#ifdef G4ANALYSIS_USE_ROOT

  // Divide by the number of events in each energy bin:
 
  G4double factor = 1./eventsPerBin;

  histo[0]->Scale(factor);
  histo[1]->Scale(factor);
  histo[2]->Scale(factor);
  histo2[0]->Scale(factor);
  histo2[1]->Scale(factor);

  // Special procedure to calculate the average Z for given A

  G4double content = 0.;
  G4double NofEventsInTheBin = 0.;
  
  for (G4int i=0; i<=histo[3]->GetNbinsX(); i++){
   
    content = histo[3]->GetBinContent(i);
    NofEventsInTheBin = histo[4]->GetBinContent(i);

    if (NofEventsInTheBin > 0.) {
      content /= NofEventsInTheBin;
      histo[3]->SetBinContent(i,content); 
    }
  }; 
 
#endif   

}

void MFTestHistoManager::CleanHisto()
{
#ifdef G4ANALYSIS_USE_ROOT
  fFile->Write();
  G4cout << "\n----> Histograms were written into the file " << fileFullName << G4endl;
  //  delete [] histo;
  //  delete [] histo2;
  delete fFile;
#endif
}

void MFTestHistoManager::fill(G4int i, G4double x, G4double y)
{
#ifdef G4ANALYSIS_USE_ROOT
  histo[i]->Fill(x,y);
#endif
}
  
void MFTestHistoManager::fill2(G4int i, G4double x, G4double y, G4double z)
{
#ifdef G4ANALYSIS_USE_ROOT
  histo2[i]->Fill(x,y,z);
#endif
}
