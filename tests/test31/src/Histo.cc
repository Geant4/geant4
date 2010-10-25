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
//---------------------------------------------------------------------------
//
// ClassName:   Histo - Generic histogram/ntuple manager class
//
//
// Author:      V.Ivanchenko 30.10.03
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Histo.hh"

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#include "HistoMessenger.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo* Histo::m_instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo* Histo::GetInstance()
{
  if(m_instance == 0){
    static Histo hist; 
    m_instance = &hist;
  }
  return m_instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo::Histo()
{
  m_verbose    = 0;
  m_histName   = "test31";
  m_histType   = "root";
  m_Histo      = 0;
  m_Clouds     = 0;
  m_Tuple      = 0;
  m_defaultAct = true;
  m_messenger  = 0;
  m_tupleTitle.clear();
  m_tuplePath.clear();
  m_tupleColumns.clear();
  m_ntup.clear();
#ifdef G4ANALYSIS_USE
  m_messenger = new HistoMessenger(this);
  m_af        = 0;
  m_tree      = 0;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Histo::~Histo()
{
#ifdef G4ANALYSIS_USE
  clear();
  delete m_messenger;
  delete m_af;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::clear()
{
#ifdef G4ANALYSIS_USE
  for(G4int i=0; i<m_Histo; i++) {
    m_histo[i]  = 0;
  }
  G4int nt = m_ntup.size();
  for(G4int j=0; j<nt; j++) {
    m_ntup[j] = 0;
  }
  for(G4int k=0; k<m_Clouds; k++) {
    m_cloud[k] = 0;
  }
  if(m_tree) { delete m_tree; }
  m_tree = 0;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::book()
{
#ifdef G4ANALYSIS_USE
  G4cout << "### Histo books " << m_Histo << " histograms " << G4endl;
  // Creating the analysis factory
  if(!m_af) m_af = AIDA_createAnalysisFactory();

  // Creating the tree factory
  AIDA::ITreeFactory* tf = m_af->createTreeFactory();

  // Creating a tree mapped to a new file.
  G4String comp = "";

  G4String nam = m_histName + "." + m_histType;

  m_tree = tf->create(nam,m_histType,false,true,comp);
  if(m_tree) 
    G4cout << "Tree store  : " << m_tree->storeName() << G4endl;
  else
    G4cout << "ERROR: Tree store " << nam  << " is not created!" << G4endl;
  delete tf;

  // Creating a histogram factory, whose histograms will be handled by the tree
  AIDA::IHistogramFactory* hf = m_af->createHistogramFactory(*m_tree);

  // Creating an 1-dimensional histograms in the root directory of the tree
  G4int i;
  for(i=0; i<m_Histo; i++) {
    if(m_active[i]) {
      G4String ss = m_ids[i];
      if(m_histType == "root") ss = "h" + m_ids[i]; 
      m_histo[i] = hf->createHistogram1D(ss, m_titles[i], m_bins[i], 
					 m_xmin[i], m_xmax[i]);
      if(m_histo[i] && m_verbose>0) ListHistogram(i);
    }
  }
  for(i=0; i<m_Clouds; i++) {
    if(m_activeCl[i]) m_cloud[i] = hf->createCloud1D(m_titlesCl[i]);
  }
  delete hf;

  // Creating a tuple factory, whose tuples will be handled by the tree
  G4int nt = m_ntup.size();
  if(0 < nt) {
    G4cout << "### Histo books " << nt << " tuples " << G4endl;
    AIDA::ITupleFactory* tpf = m_af->createTupleFactory(*m_tree);
    for(G4int i=0; i<nt; i++) {
       m_ntup[i] = tpf->create(m_tuplePath[i], m_tupleTitle[i], 
			       m_tupleColumns[i]);
    }
    delete tpf;
  }
#endif
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::save()
{
#ifdef G4ANALYSIS_USE
  // Write histogram file
  m_tree->commit();
  if(m_verbose>0) G4cout << "Histo: Closing the tree " << G4endl;
  m_tree->close();
  G4cout << "Histograms and Ntuples are saved" << G4endl;
  clear();
  if(m_verbose>0) G4cout << "Tree is deleted " << G4endl;
#endif
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int Histo::add1D(const G4String& id, const G4String& name, G4int nb, 
                         G4double x1, G4double x2, G4double u)
{
  G4int ID = m_Histo;
  m_Histo++;
  G4String hid = id;
  if(hid == "0" || hid == "none") {
    char buffer [10];
    sprintf(buffer,"%d",m_Histo);
    hid = buffer;
  }
  G4double uu = 1.0;
  if(u > 0.0) uu = u;
  if(m_verbose > 0) {
    G4cout << "New histogram will be booked: id= " << ID 
           << "  " << hid << "  <" << name 
           << ">  " << nb << "  " << x1 << "  " << x2 << "  " << uu 
           << G4endl;
  }
  m_active.push_back(m_defaultAct);
  m_bins.push_back(nb);
  m_xmin.push_back(x1);
  m_xmax.push_back(x2);
  m_unit.push_back(uu);
  m_ids.push_back(hid);
  m_titles.push_back(name);
  m_histo.push_back(0);
  return ID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int Histo::addCloud1D(const G4String& name)
{
  G4int ID = m_Clouds;
  m_Clouds++;
  if(m_verbose > 0) {
    G4cout << "New cloud will be booked: id= " << ID 
           << "  <" << name 
           << ">  " 
           << G4endl;
  }
  m_activeCl.push_back(m_defaultAct);
  m_titlesCl.push_back(name);
  m_cloud.push_back(0);
  return ID;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::setHisto1D(G4int i, G4int nb, G4double x1, G4double x2, G4double u)
{
  if(i>=0 && i<m_Histo) {
    G4double uu = m_unit[i];
    if(u > 0.0) uu = u;
    if(m_verbose > 0) {
      G4cout << "Histo:Update histogram ID= " << i  
             << "  " << nb << "  " << x1 << "  " << x2 << "  " << uu 
             << G4endl;
    }
    m_bins[i] = nb;
    m_xmin[i] = x1;
    m_xmax[i] = x2;
    m_unit[i] = uu;
    m_active[i] = true;
  } else {
    G4cout << "Histo::setHisto1D: WARNING! wrong histogram index " << i << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::activate(G4int i, G4bool w)
{
  if(m_verbose > 0) {
    G4cout << "Histo: histogram: ID= " << i 
           << "  activation= " << w
           << G4endl;   
  }
  if(i>=0 && i<m_Histo) m_active[i] = w;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::activateCloud(G4int i, G4bool w)
{
  if(m_verbose > 0) {
    G4cout << "Histo: cloud: ID= " << i 
           << "  activation= " << w
           << G4endl;   
  }
  if(i>=0 && i<m_Clouds) m_activeCl[i] = w;
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::fill(G4int i, G4double x, G4double w)
{
  if(m_verbose > 1) {
    G4cout << "Histo:fill histogram ID= " << i << " at x= " << x 
           << "  weight= " << w   
           << G4endl;   
  }
#ifdef G4ANALYSIS_USE  
  if(i>=0 && i<m_Histo) {
    if(m_histo[i]) m_histo[i]->fill(x/m_unit[i], w);
  } else {
    G4cout << "Histo::fill: WARNING! wrong histogram ID " << i << G4endl;
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::fillCloud(G4int i, G4double x, G4double w)
{
  if(m_verbose > 1) {
    G4cout << "Histo:fill cloud ID= " << i << " at x= " << x 
           << "  weight= " << w   
           << G4endl;   
  }
#ifdef G4ANALYSIS_USE  
  if(i>=0 && i<m_Clouds) {
    if(m_cloud[i]) m_cloud[i]->fill(x, w);
  } else {
    G4cout << "Histo::fillCloud: WARNING! wrong histogram ID " << i << G4endl;
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::scale(G4int i, G4double x)
{
  if(m_verbose > 0) {
    G4cout << "Histo:Scale histogram ID= " << i << " by factor " << x << G4endl;   
  }
#ifdef G4ANALYSIS_USE  
  if(i>=0 && i<m_Histo) {
    if(m_histo[i]) m_histo[i]->scale(x);
  } else {
    G4cout << "Histo::scale: WARNING! wrong histogram ID " << i << G4endl;
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int Histo::addTuple(const G4String& tuplePath, const G4String& tupleTitle, 
                      const G4String& tupleColumns)
{
  m_Tuple++;
  G4String nid = tuplePath;
  if(nid == "" || nid == "0" || nid == "none") {
    char buffer [33];
    sprintf(buffer,"%d",1000+m_Tuple);
    nid = buffer;
  }
  m_tuplePath.push_back(nid);
  m_tupleTitle.push_back(tupleTitle);
  m_tupleColumns.push_back(tupleColumns);
  m_ntup.push_back(0);     
  return m_Tuple-1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::fillTuple(G4int id, const G4String& parname, G4double x)
{
  if(m_verbose > 1) {
    G4cout << "Histo:fill tuple ID= " << id 
           << " by parameter <" << parname << "> = " << x << G4endl; 
  }
#ifdef G4ANALYSIS_USE  
  if(id>=0 && id<m_Tuple) {
    if(m_ntup[id]) m_ntup[id]->fill(m_ntup[id]->findColumn(parname), (float)x);
  }
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::addRow(G4int id)
{
  if(m_verbose > 1) 
    G4cout << "Histo:addRow to ntupe ID= " << id << G4endl; 
  
#ifdef G4ANALYSIS_USE  
  if(id>=0 && id<m_Tuple) {
    if(m_ntup[id]) m_ntup[id]->addRow();
  }
#endif
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::setFileName(const G4String& nam) 
{
  m_histName = nam;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::setFileType(const G4String& nam) 
{
  if(nam == "hbook" || nam == "root") m_histType = nam;
  else m_histType = "aida";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::setVerbose(G4int val)
{
  m_verbose = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::ListHistogram(G4int val)
{
  G4int imin = val;
  G4int imax = val;
  if(val < 0 || val >= m_Histo) {
    imin = 0;
    imax = m_Histo-1;
  }
  for(G4int i=imin; i<=imax; i++) {
    G4cout << "### Histogramm #" << i << " " << m_ids[i] << " " 
           << m_titles[i] << G4endl; 
    G4cout << "                       nbins= " << m_bins[i]
           << " xmin= " << m_xmin[i]/m_unit[i]
           << " xmax= " << m_xmax[i]/m_unit[i]
           << " unit= " << m_unit[i]
           << " flag= " << m_active[i]
           << " " << m_histo[i]
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int Histo::NumberOfBins(G4int id)
{
  G4int n = 0;
  if(id>=0 && id<m_Histo) n = m_bins[id];
  return n;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double Histo::MinBin(G4int id)
{
  G4double x = 1.0e-9;
  if(id>=0 && id<m_Histo) x = m_xmin[id]*m_unit[id];
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
G4double Histo::MaxBin(G4int id)
{
  G4double x = 1.0e-9;
  if(id>=0 && id<m_Histo) x = m_xmax[id]*m_unit[id];
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool Histo::IsActive(G4int id)
{
  G4bool yes = false;
  if(id>=0 && id<m_Histo) yes = m_active[id];
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Histo::PrintHisto(G4int id)
{
  if(id>=0 && id<m_Histo) {
#ifdef G4ANALYSIS_USE  
    // if(m_active[i]) m_histo[id]->print();
#endif
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

