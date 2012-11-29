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
// $Id$
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "HistoMessenger.hh"
#include "G4UnitsTable.hh"

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
:factoryOn(false),af(0),tree(0)
{
#ifdef G4ANALYSIS_USE
  // Creating the analysis factory
  af = AIDA_createAnalysisFactory();
  if(!af) {
    G4cout << " HistoManager::HistoManager :" 
           << " problem creating the AIDA analysis factory."
           << G4endl;
  }	   
#endif
 
  fileName[0] = "ecal";
  fileType    = "root";  
  fileOption  = "export=root";
    
  // histograms
  for (G4int k=0; k<MaxHisto; k++) { 
    histo[k] = 0;
    exist[k] = false;
    Unit[k]  = 1.0;
    Width[k] = 1.0;
    ascii[k] = false;        
  }
    
  // ntuples
  for (G4int k=0; k<MaxNtupl; k++) { 
    ntupl[k] = 0;
    existNt[k] = false;
  }
       
  histoMessenger = new HistoMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{ 
  delete histoMessenger;
  
#ifdef G4ANALYSIS_USE  
  delete af;
#endif     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::book()
{ 
#ifdef G4ANALYSIS_USE
  if(!af) return;    	    
 
 // Creating a tree mapped to an histogram file.
 fileName[1] = fileName[0] + "." + fileType;
 G4bool readOnly  = false;
 G4bool createNew = true;
 AIDA::ITreeFactory* tf  = af->createTreeFactory(); 
 tree = tf->create(fileName[1], fileType, readOnly, createNew, fileOption);
 delete tf;
 if(!tree) {
   G4cout << " HistoManager::book :" 
          << " problem creating the AIDA tree with "
          << " storeName = " << fileName[1]
          << " storeType = " << fileType
          << " readOnly = "  << readOnly
          << " createNew = " << createNew
          << " options = "   << fileOption
          << G4endl;
   return;
 }
 
 // Creating a histogram factory, whose histograms will be handled by the tree
 AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);
 
 // create selected histograms
 for (G4int k=1; k<MaxHisto; k++) {
   if (exist[k]) {
     histo[k] = hf->createHistogram1D( Label[k], Title[k],
                                                 Nbins[k], Vmin[k], Vmax[k]);
     factoryOn = true;
   }  						  					   
 }
 delete hf;
 
 // Creating a ntuple factory, handled by the tree
 AIDA::ITupleFactory* ntf = af->createTupleFactory(*tree);
 
 // create selected ntuples
 for (G4int k=1; k<MaxNtupl; k++) {
   if (existNt[k]) {
     ntupl[k] = ntf->create( LabelNt[k], TitleNt[k], ColumnNt[k]);
     factoryOn = true;
   }  						  					   
 }
 
 delete ntf;
  
 if (factoryOn) 
     G4cout << "\n----> Histogram Tree is opened in " << fileName[1]  << G4endl;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::save()
{ 
#ifdef G4ANALYSIS_USE
  if (factoryOn) {
    saveAscii();          // Write ascii file, if any     
    tree->commit();       // Writing the histograms to the file
    tree->close();        // and closing the tree (and the file)
    G4cout << "\n----> Histogram Tree is saved in " << fileName[1] << G4endl;

    delete tree;
    tree = 0;
    factoryOn = false;
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillHisto(G4int ih, G4double xbin, G4double weight)
{
  if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::FillHisto() : histo " << ih
           << " xbin= " << xbin << " weight= " << weight << G4endl;
    return;
  }
#ifdef G4ANALYSIS_USE
  if(exist[ih]) histo[ih]->fill(xbin/Unit[ih], weight);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetHisto(G4int ih, 
            G4int nbins, G4double valmin, G4double valmax, const G4String& unit)
{   
  if (ih < 1 || ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::SetHisto() : histo " << ih
           << "does not exist" << G4endl;
    return;
  }
  	 
  const G4String id[] = { "0", "1", "2", "3", "4", "5" };
  const G4String title[] = 
                { "dummy",			//0
                  "total Evis in Ecal",		//1
                  "total Edep in Ecal",		//2
                  "Evis profile",		//3
		  "Edep profile",		//4
		  "Nb of Radiation Length"	//5
                 };
		 			 
  G4String titl = title[ih];
  G4double vmin = valmin, vmax = valmax;
  Unit[ih] = 1.;

  if (unit != "none") {
    titl = title[ih] + " (" + unit + ")";
    Unit[ih] = G4UnitDefinition::GetValueOf(unit);
    vmin = valmin/Unit[ih]; vmax = valmax/Unit[ih];
  }
            
  exist[ih] = true;
  Label[ih] = id[ih];
  Title[ih] = titl;
  Nbins[ih] = nbins; 
  Vmin[ih]  = vmin; 
  Vmax[ih]  = vmax; 
  Width[ih] = (valmax-valmin)/nbins;
  
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
#ifdef G4ANALYSIS_USE
  if(exist[ih]) histo[ih]->scale(fac);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::RemoveHisto(G4int ih) 
{ 
 if (ih >= MaxHisto) {
    G4cout << "---> warning from HistoManager::RemoveHisto() : histo " << ih
           << "does not exist" << G4endl;
    return;
  }
  	  
  histo[ih] = 0; exist[ih] = false;     		 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::PrintHisto(G4int ih)
{
 if (ih < MaxHisto)  { ascii[ih] = true; ascii[0] = true; }
 else
    G4cout << "---> warning from HistoManager::PrintHisto() : histo " << ih
           << "does not exist" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <fstream>

void HistoManager::saveAscii()
{
#ifdef G4ANALYSIS_USE

 if (!ascii[0]) return;
  
 G4String name = fileName[0] + ".ascii";
 std::ofstream File(name, std::ios::out);
 File.setf( std::ios::scientific, std::ios::floatfield );
 
 //write selected histograms
 for (G4int ih=0; ih<MaxHisto; ih++) {
    if (exist[ih] && ascii[ih]) {
      File << "\n  1D histogram " << ih << ": " << Title[ih] 
           << "\n \n \t     X \t\t     Y" << G4endl;
     
      for (G4int iBin=0; iBin<Nbins[ih]; iBin++) {
         File << "  " << iBin << "\t" 
              << 0.5*(histo[ih]->axis().binLowerEdge(iBin) +
	              histo[ih]->axis().binUpperEdge(iBin)) << "\t"
	      << histo[ih]->binHeight(iBin) 
	      << G4endl;
      } 
    }
 }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::SetNtuple(G4int nt) 
{   
  if (nt < 1 || nt >= MaxNtupl) {
    G4cout << "---> warning from HistoManager::SetNtuple() : Ntuple " << nt
           << "does not exist" << G4endl;
    return;
  }
  	 
  const G4String id[] = { "100", "101" };
  const G4String title[] = 
                { "dummy",			//0
                  "Energy deposit in subModule"	//1
                 };
		 
  G4String column[MaxNtupl];
  
  column[0] = " int dum=0 ";			//0
  
  column[1] =
  "double Evis0, Evis1, Evis2, Evis3, Evis4, Evis5, Evis6, Evis7, Evis8, Evis9";
  column[1] +=
  ", Evis10, Evis11, Evis12, Evis13, Evis14, Evis15, Evis16, Evis17";
  column[1] +=
  ", Edep0, Edep1, Edep2, Edep3, Edep4, Edep5, Edep6, Edep7, Edep8, Edep9";
  column[1] +=
  ", Edep10, Edep11, Edep12, Edep13, Edep14, Edep15, Edep16, Edep17";
    			 
  G4String titl = title[nt];
            
   existNt[nt] = true;
   LabelNt[nt] = id[nt];
   TitleNt[nt] = titl;
  ColumnNt[nt] = column[nt];
  
  G4cout << "----> SetNtuple " << nt << ": " << titl << ";  " << G4endl;
   		 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::FillNtuple(G4int nt, G4int column, G4double value)
{
  if (nt >= MaxNtupl) {
    G4cout << "---> warning from HistoManager::FillNtuple() : Ntuple " << nt
           << " does not exist " << column << value << G4endl;
    return;
  }
#ifdef G4ANALYSIS_USE
  if(existNt[nt]) ntupl[nt]->fill(column, value);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::AddRowNtuple(G4int nt)
{
  if (nt >= MaxNtupl) {
    G4cout << "---> warning from HistoManager::AddRowNtuple() : Ntuple " << nt
           << " do not exist" << G4endl;
    return;
  }
#ifdef G4ANALYSIS_USE
  if(existNt[nt]) ntupl[nt]->addRow();
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


