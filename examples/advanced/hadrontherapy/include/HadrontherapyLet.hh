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
// $Id: HadrontherapyLet.hh,v 1.0, May 2007;
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), M. Sallemi, A. Salvia
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#ifndef HadrontherapyLet_h
#define HadrontherapyLet_h 1

#include "globals.hh"
#include <fstream>
#include <vector>
#include <string>
using namespace std;

class HadrontherapyLetMessenger;
class HadrontherapyPrimaryGeneratorAction;

class HadrontherapyLet
{
public:
  HadrontherapyLet(HadrontherapyPrimaryGeneratorAction*);
  ~HadrontherapyLet();

public:
  void Fluence_Let(G4int,G4double);
  void LetOutput(); 
  void SetValue(G4String);

  HadrontherapyLetMessenger* letMessenger;

private:
  ofstream flet,ffluence,*fspectrum;  
  ifstream fstop; 
 
  G4int **spectrum, i, j, size, bins;
  G4double energy, *stop, n1, d1, n2, d2, lett, letd;
  G4String nome_file;

  vector<double> vetdepth;  //vector containing the points where the energy spectrum and the let will be calculated (depth in mm)

  HadrontherapyPrimaryGeneratorAction* pga; 
 
};
#endif
