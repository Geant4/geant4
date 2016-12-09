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
#define USE_CANVASINTAB

#ifdef USE_CANVASINTAB
#include "CanvasInTab.hh"
#endif

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <cstring>

#include <TApplication.h>
#include <TGApplication.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TBranch.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGFileBrowser.h>
#include <TGFileDialog.h>
#include <TChain.h>
#include <TColor.h>
using namespace std;

//------------------------------------------------------------------------------

const TGFileInfo* OpenRootFile()
{
    const char *gOpenAsTypes[] = {
    "ROOT files",   "*.root",
    "All files",    "*"
  };

  static TGFileInfo fi;
  fi.fFileTypes = gOpenAsTypes;
  //  fi.SetMultipleSelection(kTRUE);
  // User must check the box "multiple selection" in the dialog box
  //  fi.fIniDir = StrDup(".");
  new TGFileDialog(gClient->GetRoot(),
                   gClient->GetRoot(), kFDOpen, &fi);

  return &fi;
}

//------------------------------------------------------------------------------

struct SpeciesInfoAOS
{
  SpeciesInfoAOS()
  {
    fNEvent = 0;
    fNumber = 0;
    fG = 0.;
    fG2 = 0.;
  }

  SpeciesInfoAOS(const SpeciesInfoAOS& right) // Species A(B);
  {
    fNEvent = right.fNEvent;
    fNumber = right.fNumber;
    fG = right.fG;
    fG2 = right.fG2;
    fName = right.fName;
  }

  SpeciesInfoAOS& operator=(const SpeciesInfoAOS& right) // A = B
  {
    if(&right == this) return *this;
    fNEvent = right.fNEvent;
    fNumber = right.fNumber;
    fG = right.fG;
    fG2 = right.fG2;
    fName = right.fName;
    return *this;
  }

  int fNEvent;
  int fNumber;
  double fG;
  double fG2;
  string fName;
};

//------------------------------------------------------------------------------

struct SpeciesInfoSOA
{
  SpeciesInfoSOA()
  {
    fRelatErr = 0;
  }

  SpeciesInfoSOA(const SpeciesInfoSOA& right) :
  fG(right.fG),
  fGerr(right.fGerr),
  fTime(right.fTime),
  fRelatErr(right.fRelatErr),
  fName(right.fName)
  {}

  SpeciesInfoSOA& operator=(const SpeciesInfoSOA& right)
  {
    if(this == &right) return *this;
    fG = right.fG;
    fGerr = right.fGerr;
    fTime = right.fTime;
    fRelatErr = right.fRelatErr;
    fName = right.fName;
    return *this;
  }

  std::vector<double> fG;
  std::vector<double> fGerr;
  std::vector<double> fTime;
  double fRelatErr;
  string fName;
};

//------------------------------------------------------------------------------

void ProcessSingleFile(TFile* file)
{
  int speciesID;
  int number;
  int nEvent;
  char speciesName[500];
  double time;  // time
  double sumG;  // sum of G over all events
  double sumG2; // sum of G^2 over all events

  TTree* tree = (TTree*)file->Get("species");
  tree->SetBranchAddress("speciesID", &speciesID);
  tree->SetBranchAddress("number", &number);
  tree->SetBranchAddress("nEvent", &nEvent);
  tree->SetBranchAddress("speciesName", &speciesName);
  tree->SetBranchAddress("time", &time);
  tree->SetBranchAddress("sumG", &sumG);
  tree->SetBranchAddress("sumG2", &sumG2);

  Long64_t nentries = tree->GetEntries();
  // cout << nentries <<" entries" << endl;

  if(nentries == 0)
  {
    cout << "No entries found in the tree species contained in the file "
         << file->GetPath() << endl;
    exit(1);
  }

  //----------------------------------------------------------------------------
  // This first loop is used in case the processed ROOT file is issued from the
  // accumulation of several ROOT files (e.g. hadd)

  std::map<int, std::map<double, SpeciesInfoAOS>> speciesTimeInfo;

  for (int j=0; j < nentries; j++)
  {
    tree->GetEntry(j);

    SpeciesInfoAOS& infoAOS = speciesTimeInfo[speciesID][time];

    infoAOS.fNumber += number;
    infoAOS.fG += sumG;
    infoAOS.fG2 += sumG2;
    infoAOS.fNEvent += nEvent;
    infoAOS.fName = speciesName;
  }

  //----------------------------------------------------------------------------

  std::map<int, SpeciesInfoSOA> speciesInfo;

  auto it_SOA = speciesTimeInfo.begin();
  auto end_SOA = speciesTimeInfo.end();

  for (; it_SOA!=end_SOA ; ++it_SOA)
  {
    const int _speciesID = it_SOA->first;
    SpeciesInfoSOA& info = speciesInfo[_speciesID];

    auto it2 = it_SOA->second.begin();
    auto end2 = it_SOA->second.end();

    info.fName = it2->second.fName;
    const size_t size2 = it_SOA->second.size();
    info.fG.resize(size2);
    info.fGerr.resize(size2);
    info.fTime.resize(size2);

    for(int i2 = 0 ;it2!=end2;++it2, ++i2)
    {
      SpeciesInfoAOS& infoAOS = it2->second;
      
      double _SumG2 = infoAOS.fG2;
      double _MeanG = infoAOS.fG/infoAOS.fNEvent;
      double _Gerr = sqrt((_SumG2/infoAOS.fNEvent - pow(_MeanG,2))
                          /(infoAOS.fNEvent-1) );

      info.fG[i2] = _MeanG;
      info.fGerr[i2] = _Gerr;
      info.fTime[i2] = it2->first;
      info.fRelatErr += _Gerr/(_MeanG + 1e-30); // add an epsilon to prevent NAN
    }
  }

  //----------------------------------------------------------------------------

#ifdef USE_CANVASINTAB
  CanvasInTab* myFrame =
  new CanvasInTab(gClient->GetRoot(), 500, 500);
#endif

  std::map<int, SpeciesInfoSOA>::iterator it = speciesInfo.begin();
  std::map<int, SpeciesInfoSOA>::iterator end = speciesInfo.end();

  for (; it != end; ++it)
  {
    speciesID = it->first;
    SpeciesInfoSOA& info = it->second;
//    if(strstr(info.fName.c_str(), "H2O^") != 0) continue;
    
    if(info.fG.empty()) continue;

    TGraphErrors* gSpecies = new TGraphErrors(info.fG.size(),
                                              info.fTime.data(),
                                              info.fG.data(),
                                              0,
                                              info.fGerr.data());

#ifdef USE_CANVASINTAB
    int nCanvas = myFrame->AddCanvas(info.fName.c_str());
    myFrame->GetCanvas(nCanvas);
    TCanvas* cSpecies = myFrame->GetCanvas(nCanvas);
#else
    TCanvas* cSpecies = new TCanvas(info.fName.c_str(),
                                    info.fName.c_str());
#endif

    cSpecies->cd();
    int color = (2+speciesID)%TColor::GetNumberOfColors();
    if(color == 5 || color==10 || color==0) ++color;
    
    // cout << info.fName.c_str() << " " << color << endl;
    
    gSpecies->SetMarkerStyle(20+speciesID);
    gSpecies->SetMarkerColor(color);
    info.fRelatErr /= (double)info.fG.size();

    gSpecies->SetTitle((info.fName
                        + " - speciesID: "
                        + std::to_string(speciesID)+" rel. Err. "
                        + std::to_string(info.fRelatErr)).c_str() );
    gSpecies->GetXaxis()->SetTitle("Time [ns]");
    gSpecies->GetYaxis()->SetTitle("G [molecules/100 eV]");
    gSpecies->Draw("ap");
    cSpecies->SetLogx();
  }
  
#ifdef USE_CANVASINTAB
  int nCanvas = myFrame->GetNCanvas();
  for(int i = 0 ; i < nCanvas ; ++i)
  {
    myFrame->GetCanvas(i)->Update();
  }
#endif
}

//------------------------------------------------------------------------------

int ProcessSingleFile(const char* filePath)
{
  if(filePath == 0 || strlen(filePath) == 0)
  {
    perror("You must provide a valid file");
    return 1;
  }

  TFile* file = TFile::Open(filePath);

  if(file == 0)
  {
    perror ("Error opening ntuple file");
    exit(1);
  }
  
  if(!file-> IsOpen())
  {
    perror ("Error opening ntuple file");
    exit(1);
  }
  else
  {
    cout << "Opening ntple file " << filePath << endl;
  }
  ProcessSingleFile(file);
  return 0;
}

//------------------------------------------------------------------------------

#define _PROCESS_ONE_FILE_ ProcessSingleFile
//#define _PROCESS_ONE_FILE_ ProcessSingleFileTProfile

int main(int argc, char **argv)
{
  //--------------------------------
  int initialArgc = argc;
  vector<char*> initialArgv(argc);
  for(int i = 0 ; i < argc ; ++i)
  {
    initialArgv[i] = argv[i];
  }
  //--------------------------------

  TApplication* rootApp = new TApplication("PlotG",&argc, argv);

  const char* filePath = 0;

  if(initialArgc == 1) // no file provided in argument
  {
    const TGFileInfo* fileInfo = OpenRootFile();
    filePath = fileInfo->fFilename;
    if(fileInfo->fFileNamesList && fileInfo->fFileNamesList->GetSize()>1)
    {
      // several files selected
      // user has to tick "Multiple selection"
      perror("Multiple selection of files not supported, implement your own!");
     //
     // For instance, start from:
     //   TChain* tree = new TChain("species");
     //   tree->AddFileInfoList(fileInfo->fFileNamesList);
     // Or call ProcessSingleFile for each file,
     // you'll need to do some adaptation
    }
    else
    {
      if(_PROCESS_ONE_FILE_(filePath)) return 1;
    }
  }
  else // a file is provided in argument
  {
    filePath = initialArgv[1];
    if(_PROCESS_ONE_FILE_(filePath)) return 1;
  }

  rootApp->Run();
  delete rootApp;
  return 0;
}
