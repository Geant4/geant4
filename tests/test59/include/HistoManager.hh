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
// $Id: HistoManager.hh,v 1.1 2009-11-27 16:06:28 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Short description: A simple histogramming tool with ASCII output
//
//---------------------------------------------------------------------

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include <vector>
#include "G4UnitsTable.hh"

class HistoMessenger;

class HistoManager
{
 public:

  HistoManager();
  ~HistoManager();

  void BookHisto();
  void SaveHisto();
  void FillHisto(G4int ih, G4double e, G4double weight = 1.);

  // Temporary fake functions to compile the messenger
  void SetFileName   (const G4String& ) {;}
  void SetFileType   (const G4String& ) {;}
  void SetFileOption (const G4String& ) {;}
  void SetHisto (G4int, G4int, G4double, G4double, const G4String& ){;}
  void PrintHisto  (G4int){;}
  void RemoveHisto (G4int){;} 
   
 private:
        
  HistoMessenger*          histoMessenger;

  // vectors for filling of data

  std::vector<G4double>               minValH; // Start of the histogram
  std::vector<G4double>               stepHis; // Step of the histogram
  std::vector<G4int>                  nBinHis; // Number of bins of histograms
  std::vector<std::vector<G4double>*> hisBins; // Bin Number (*weight) of the histograms
  std::vector<std::vector<G4double>*> hisBinX; // Real bins of the histograms
};

#endif

