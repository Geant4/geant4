//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: Histo.hh,v 1.5 2004-09-21 10:41:49 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef Histo_h
#define Histo_h 1

//---------------------------------------------------------------------------
//
// ClassName:   Histo
//
// Description: Utility class to hold and manipulate histograms/nTuples
//
// Author:      V.Ivanchenko 30/10/03
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include <vector>
#include "G4DynamicParticle.hh"
#include "G4VPhysicalVolume.hh"
#include "G4DataVector.hh"
#include "G4Track.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class HistoMessenger;

namespace AIDA {
 class ITree;
 class ITuple;
 class IHistogram1D;
}

class Histo
{
public:
  static Histo* GetInstance();

  virtual ~Histo();

private:
  Histo();

public:
  void clear();
  // delete histogramms 

  void book();
  // Book predefined histogramms 

  void save();
  // Save histogramms to file

  G4int add1D(const G4String&, const G4String&, G4int nb=100, G4double x1=0., 
                                                G4double x2=1., G4double u=1.);
  // In this method histogramms are predefined

  void setHisto1D(G4int id, G4int nb=100, 
                  G4double x1=0., G4double x2=1., G4double u=1.);
  // It change bins and boundaries

  void activate(G4int, G4bool val=true);
  // Histogramms are filled

  void fill(G4int, G4double, G4double w=1.);
  // Histogramms are filled

  void scale(G4int, G4double);

  G4int addTuple(const G4String&, const G4String&, const G4String&);
  // In this method nTuple is booked

  void fillTuple(G4int, const G4String&, G4double);
  // Fill nTuple parameter

  void addRow(G4int);
  // Save tuple event 

  void setFileName(const G4String&);

  void setFileType(const G4String&);

  void PrintHisto(G4int id=0);

  void ListHistogram(G4int val);

  void setVerbose(G4int val);

  G4int NumberOfBins(G4int id);

  G4double MinBin(G4int id);

  G4double MaxBin(G4int id);

  G4bool IsActive(G4int id);

private:

  static Histo* m_instance;

  G4String      m_histName;
  G4String      m_histType;
  std::vector<G4String> m_tuplePath;
  std::vector<G4String> m_tupleTitle;
  std::vector<G4String> m_tupleColumns;

  G4int         m_Histo;
  G4int         m_Tuple;
  G4int         m_verbose;
  G4bool        m_defaultAct;

  std::vector<AIDA::IHistogram1D*> m_histo;
  std::vector<AIDA::ITuple*>       m_ntup;

  AIDA::ITree*    m_tree;
  HistoMessenger* m_messenger;

  std::vector<G4double>  m_xmin;
  std::vector<G4double>  m_xmax;
  std::vector<G4double>  m_unit;
  std::vector<G4bool>    m_active;
  std::vector<G4int>     m_bins;
  std::vector<G4String>  m_ids;
  std::vector<G4String>  m_titles;
};

#endif
