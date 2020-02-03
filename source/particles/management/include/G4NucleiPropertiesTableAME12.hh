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
//
// -------------------------------------------------------------------
//
//      File name:     G4NucleiPropertiesTableAME12.cc 
//
//      Authors:       Tatsumi Koi (tkoi@slac.stanford.edu)
// 
//      Data are update to AME2012
//      "The Ame2012 atomic mass evaluation (I)"  
//      	by G.Audi, M.Wang, A.H.Wapstra, F.G.Kondev, M.MacCormick, X.Xu, and B.~Pfeiffer
//      	Chinese Physics C36 p. 1287-1602, December 2012.
//      "The Ame2012 atomic mass evaluation (II)"  
//      	by M.Wang, G.Audi, A.H.Wapstra, F.G.Kondev, M.MacCormick, X.Xu, and B.~Pfeiffer
//              Chinese Physics C36 p. 1603-2014, December 2012.
//
//      Creation date: Aug. 2016 
//                     based on G4NucleiPropertiesTableAME03
//
//      Modifications: 
// 

#ifndef G4NucleiPropertiesTableAME12_h
#define G4NucleiPropertiesTableAME12_h  1

#include <cmath>
#include "globals.hh"
//
class G4NucleiProperties;

class G4NucleiPropertiesTableAME12 
{
private:
  
  // Default constructor - this class should exist only once!
  G4NucleiPropertiesTableAME12();

public:

  // Destructor (generated)
  ~G4NucleiPropertiesTableAME12() { };

  // Following values migrate to AME12
  enum  {nEntries = 3353,MaxA = 295, ZMax = 120}; 

  // Other Operations 
  // all methods are private and can be used only by G4NucleiProperties

  friend class G4NucleiProperties;  

private:

  // Operation: GetMassExcess
  //   values imported from The Ame2003 atomic mass evaluation (II)  
  static G4double GetMassExcess(G4int Z, G4int A); 

  // Operation: GetAtomicMass .. in Geant4 Energy units!
  //      Atomic_Mass =  MassExcess + A*amu_c2
  static G4double GetAtomicMass(G4int Z, G4int A);

  // Operation: GetNuclearMass
  //      Nuclear_Mass = Atomic_Mass - electronMass
  static G4double GetNuclearMass(G4int Z, G4int A);

  // Operation: GetBindingEnergy
  static G4double GetBindingEnergy(G4int Z, G4int A);

  // Operation: GetBetaDecayEnergy
  static G4double GetBetaDecayEnergy(G4int Z, G4int A);

  // Is the nucleus (A,Z) in table?
  static G4bool IsInTable(G4int Z, G4int A);

  static G4int MaxZ(G4int A);
  static G4int MinZ(G4int A);


private:

  // Operation: GetIndex
  static G4int GetIndex(G4int Z, G4int A);
  

  // Data Members for Class Attributes
  //----------------------------------  

  // The following arrays are static for allow inicialization.
  // The inicialization is Done in G4NucleiPropertiesTableAME12.cc

  // Mass Excess
  static const G4double MassExcess[nEntries];
  
  
  // Beta Decay Energy
  static const G4double BetaEnergy[nEntries];

    
  // Table of Z (number of protons) and A (number of nucleons)
  //        indexArray[0][ ] --> Z
  //        indexArray[1][ ] --> A
  static const G4int indexArray[2][nEntries];

  // Reduced Table of A for shorter index search.
  //         The index in this table coincide with A-1
  //         For each A value shortTable[A-1] has the index of the 1st occurrence in
  //         the indexArray[][]
  static const G4int shortTable[MaxA+1];

  // electrom mass
  static G4ThreadLocal G4double electronMass[ZMax];
  static G4ThreadLocal G4bool   isIntialized;

};


#endif






