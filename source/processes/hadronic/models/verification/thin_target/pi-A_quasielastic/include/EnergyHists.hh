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
//
// $Id: EnergyHists.hh,v 1.1 2003-07-31 01:15:53 dwright Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#ifndef EnergyHists_h
#define EnergyHists_h 1

#include "globals.hh"
#include <vector>

class RunAction;

class EnergyHists
{
  public:
    EnergyHists();
    EnergyHists(G4double lo, G4double hi, G4double size);
    ~EnergyHists();

    void CreateHists();
    void FillHists(G4double q, G4double ke, G4int spectrum);
    void ScaleHists(G4double factor, G4int spectrum);
    
  private:

    // energy distribution histograms

    std::vector<G4double> hist19_5deg;
    std::vector<G4double> hist30deg;
    std::vector<G4double> hist42deg;

    std::vector<std::vector<G4double>* > hists;

    G4double loBinEdge;
    G4double hiBinEdge;
    G4double binSize;

};

#endif

    




