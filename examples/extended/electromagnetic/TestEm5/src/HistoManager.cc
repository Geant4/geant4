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
// $Id: HistoManager.cc,v 1.1 2010-09-16 16:26:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "HistoManager.hh"
#include "HistoMessenger.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
{
  fileName[0]  = "testem5";
  factoryOn = false;
  fNbHist   = 0;

  // histograms
  for (G4int k=0; k<MaxHisto; k++) {
    fHistId[k] = 0;
    fHistPt[k] = 0;    
    fExist[k] = false;
    fUnit[k]  = 1.0;
    fWidth[k] = 1.0;
    fAscii[k] = false;    
  }

  fHistoMessenger = new HistoMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete fHistoMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{
  // if no histos, do nothing
  if (fNbHist == 0) return;
  
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  G4String extension = analysisManager->GetFileType();
  fileName[1] = fileName[0] + "." + extension;
  
  // Open an output file
  //
  G4bool fileOpen = analysisManager->OpenFile(fileName[0]);
  if (!fileOpen) {
    G4cout << "\n---> HistoManager::book(): cannot open " << fileName[1] 
           << G4endl;
    return;
  }  

  // create selected histograms
  //
  analysisManager->SetFirstHistoId(1);
    
  for (G4int k=0; k<MaxHisto; k++) {
    if (fExist[k]) {
      fHistId[k] = analysisManager->CreateH1( fLabel[k], fTitle[k],
                                              fNbins[k], fVmin[k], fVmax[k]);
      fHistPt[k] = analysisManager->GetH1(fHistId[k]);
      factoryOn = true;
    }
  }

  if (factoryOn)  
    G4cout << "\n----> Histogram file is opened in " << fileName[1] << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{
  if (factoryOn) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();    
    analysisManager->Write();
    analysisManager->CloseFile();
    saveAscii();                    // Write fAscii file, if any
    G4cout << "\n----> Histograms are saved in " << fileName[1] << G4endl;
      
    delete G4AnalysisManager::Instance();
    factoryOn = false;
  }         
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double e, G4double weight)
{
  if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
           << "does not fExist; e= " << e << " w= " << weight << G4endl;
    return;
  }

  if (fHistPt[ih]) {
    G4double value = e/fUnit[ih];
    if ((ih == 4) || (ih == 5)) value = std::log10(value);
    fHistPt[ih]->fill(value, weight);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetHisto(G4int ih,
               G4int nbins, G4double vmin, G4double vmax, const G4String& unit)
{
  if (ih > MaxHisto) {
    G4cout << "---> warning from HistoManager::SetHisto() : histo " << ih
           << "does not fExist" << G4endl;
    return;
  }
  
  const G4String id[] = { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                         "10","11","12","13","14","15","16","17","18","19",
			 "20","21","22","23","24","25","26","27","28","29",
			 "30","31","32","33","34","35","36","37","38","39",
			 "40","41","42","43","44","45","46","47","48","49" 
			};
			
  const G4String title[] =
                { "dummy",						//0
                  "energy deposit in absorber",				//1
                  "energy of charged secondaries at creation",		//2
		  "energy of neutral secondaries at creation",		//3
		  "energy of charged at creation(log10(ekin))",		//4
                  "energy of neutral at creation(log10(ekin))",		//5
                  "x_vertex of charged secondaries (all)",		//6
                  "x_vertex of charged secondaries (not absorbed)",	//7
		  "dummy","dummy",					//8-9
		  "(transmit, charged) : kinetic energy at exit",	//10
		  "(transmit, charged) : ener fluence: dE(MeV)/dOmega",	//11
		  "(transmit, charged) : space angle: dN/dOmega",	//12
		  "(transmit, charged) : projected angle at exit",	//13
		  "(transmit, charged) : projected position at exit",	//14
		  "(transmit, charged) : radius at exit",		//15
		  "dummy","dummy","dummy","dummy",			//16-19
		  "(transmit, neutral) : kinetic energy at exit",	//20
		  "(transmit, neutral) : ener fluence: dE(MeV)/dOmega",	//21
		  "(transmit, neutral) : space angle: dN/dOmega",	//22
		  "(transmit, neutral) : projected angle at exit",	//23
		  "dummy","dummy","dummy","dummy","dummy","dummy",	//24-29
		  "(reflect , charged) : kinetic energy at exit",	//30
		  "(reflect , charged) : ener fluence: dE(MeV)/dOmega",	//31
		  "(reflect , charged) : space angle: dN/dOmega",	//32
		  "(reflect , charged) : projected angle at exit",	//33
		  "dummy","dummy","dummy","dummy","dummy","dummy",	//34-39
		  "(reflect , neutral) : kinetic energy at exit",	//40
		  "(reflect , neutral) : ener fluence: dE(MeV)/dOmega",	//41
		  "(reflect , neutral) : space angle: dN/dOmega",	//42
		  "(reflect , neutral) : projected angle at exit",	//43
		  "dummy","dummy","dummy","dummy","dummy","dummy"	//44-49
                 };


  G4String titl = title[ih];
  fUnit[ih] = 1.;

  if (unit != "none") {
    titl = title[ih] + " (" + unit + ")";
    fUnit[ih] = G4UnitDefinition::GetValueOf(unit);
  }
  if ((ih == 4) || (ih == 5)) { 
     vmin = std::log10(vmin); vmax = std::log10(vmax);
  }
       
  fExist[ih] = true;
  fLabel[ih] = id[ih];
  fTitle[ih] = titl;
  fNbins[ih] = nbins;
  fVmin[ih]  = vmin;
  fVmax[ih]  = vmax;
  fWidth[ih] = fUnit[ih]*(vmax-vmin)/nbins;
  
  fNbHist++;

  G4cout << "----> SetHisto " << ih << ": " << titl << ";  "
         << nbins << " bins from "
         << vmin << " " << unit << " to " << vmax << " " << unit << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Normalize(G4int ih, G4double fac)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::Normalize() : histo " << ih
           << "  fac= " << fac << G4endl;
    return;
  }

  if (fHistPt[ih]) fHistPt[ih]->scale(fac);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintHisto(G4int ih)
{
 if (ih < MaxHisto) { fAscii[ih] = true; fAscii[0] = true; }
 else
    G4cout << "---> warning from HistoManager::PrintHisto() : histo " << ih
           << "does not exist" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <fstream>

void HistoManager::saveAscii()
{
 if (!fAscii[0]) return;
 
 G4String name = fileName[0] + ".ascii";
 std::ofstream File(name, std::ios::out);
 if (!File) {
   G4cout 
       << "\n---> HistoManager::saveAscii(): cannot open " << name << G4endl;
    return;
  }  

 File.setf( std::ios::scientific, std::ios::floatfield );
  
 //write selected histograms
 for (G4int ih=0; ih<MaxHisto; ih++) {
    if (fHistPt[ih] && fAscii[ih]) {
          
      File << "\n  1D histogram " << ih << ": " << fTitle[ih] 
           << "\n \n \t     X \t\t     Y" << G4endl;
     
      for (G4int iBin=0; iBin<fNbins[ih]; iBin++) {
         File << "  " << iBin << "\t" 
	      << fHistPt[ih]->axis().bin_center(iBin) << "\t"
	      << fHistPt[ih]->bin_height(iBin) 
	      << G4endl;
      } 
    }
 } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

