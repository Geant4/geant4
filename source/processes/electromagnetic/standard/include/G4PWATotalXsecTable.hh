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
// $Id: $
//
// -----------------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4PWATotalXsecTable
//
// Author:        Mihaly Novak
//
// Creation date: 18.05.2015
//
// Class description:
//   Class to load and handle elastic, first and second transport cross sections
//   precomputed by using ELSEPA [1] in the 100 eV - 1 GeV kinetic and Z = 1-103
//   energy range for electrons and positrons.G4PWATotalXsecZ is responsible to 
//   to handle cross sections by individual Z that are used in the current 
//   geometry and G4PWATotalXsecTable is a collection of G4PWATotalXsecZ objects. 
//
// Modifications:
//
// References:
//   [1] Francesc Salvat, Aleksander Jablonski, Cedric J Powell,
//       ELSEPA—Dirac partial-wave calculation of elastic scattering of electrons 
//       and positrons by atoms, positive ions and molecules,
//       Computer physics communications; 165, 2, (2005)
//
// -----------------------------------------------------------------------------


#ifndef G4PWATotalXsecTable_h
#define G4PWATotalXsecTable_h 1

#include "G4Types.hh"

////////////////////////////////////////////////////////////////////////////////
//  G4PWATotalXsecZ: sub-class for PWA xsec data that belong to a given Z number
//////////////////////////////////////////////////////////////////////////////// 
class G4PWATotalXsecZ 
{
 friend class G4PWATotalXsecTable;
public:
  //
  // out of energy grid cases
  G4int GetLowestEnergyBinIndex() const {return 0;}
  G4int GetHighestEnergyBinIndex()const {return fgNumTotalXsecBins-1;}
  G4double GetLowestEnergy() const {return fgPWATotalXsecEnergyGrid[0];}
  G4double GetHighestEnergy()const {return fgPWATotalXsecEnergyGrid[fgNumTotalXsecBins-1];}

  // see below what is input parameter j
  G4double GetLowestXsecValue(G4int j) const {return fPWAXsecs[j*fgNumTotalXsecBins];}
  G4double GetHighestXsecValue(G4int j)const {return fPWAXsecs[(j+1)*fgNumTotalXsecBins-1];}

  //
  // normal cases i.e. energy is within the grid
  // kinetic energy in MeV ; returns with the index of the lower energy bin edge 
  G4int GetPWATotalXsecEnergyBinIndex(G4double energy) const;

  //------------------------------------------------------------------------------//
  // The GetPWATotalXsecEnergyBinIndex(energy) will return with the lower energy  //
  // bin edge index = elowindx. Then the following formulas can be used to get the//
  // elastic,  first and second transport mean free path lower bin edge values:   //
  //  index of the lower energy bin edge = j*fgNumTotalXsecBins + elowindex       //  
  // where j is                                                                   //
  //  -elastic cross section lower bin edge index:           j = 1.5 + chrage*1.5 //
  //  -first transport cross section lower energy bin index: j = 2.5 + chrage*1.5 //
  //  -first transport cross section lower energy bin index: j = 3.5 + chrage*1.5 //
  // With this, we can avoid to use an IF over particle types (e-/e+)             //
  // Additional note: it's probably a good idea to separate the elowindex comp-   //
  // utation because it depends only on the energy of the particle while the      //
  // cross sections depends on Z and particle type as well                        //
  //------------------------------------------------------------------------------//
  G4double GetInterpXsec(G4double energy, G4int elowindex, G4int j) const ;
  G4double GetInterpXsec(G4double energy, G4int j) const ;

 private:
   // ctr and dtr can be called only by the G4PWATotalXsecTable friend
   G4PWATotalXsecZ(G4int Z);
   ~G4PWATotalXsecZ(){};

   //  hide assignment operator and cpy ctr.
   G4PWATotalXsecZ & operator=(const G4PWATotalXsecZ &right);
   G4PWATotalXsecZ(const G4PWATotalXsecZ&);

   void LoadPWATotalXsecZ(G4int Z);
 private:

   //size of the common energy grid //
   static const G4int fgNumTotalXsecBins =  106;

   // common energy grid in [1.e-4;1.e+3] MeV //
   // size is fgNumTotalXsecBins
   static const G4double fgPWATotalXsecEnergyGrid[fgNumTotalXsecBins];

   // elastic cross sections, first and second transport cross sections for e-/e+ 
   // over the common energy grid fgPWATotalXsecEnergyGrid in Geant4 internal length^2 
   G4double fPWAXsecs[fgNumTotalXsecBins*6];
   // interpolation parameters if log-log linear interpolation is used 
   G4double fInterpParamA[fgNumTotalXsecBins*6];
   G4double fInterpParamB[fgNumTotalXsecBins*6];
};


////////////////////////////////////////////////////////////////////////////////
//  G4PWATotalXsecTable
//////////////////////////////////////////////////////////////////////////////// 
class G4PWATotalXsecTable 
{
 public:
  G4PWATotalXsecTable() {};
   ~G4PWATotalXsecTable();

   void Initialise();

   const G4PWATotalXsecZ* GetPWATotalXsecForZet(G4int Z) const{
        Z = Z>fgNumZet ? fgNumZet : Z;  
        return fgPWATotalXsecTable[Z-1];
   } 

 private:
   //  hide assignment operator and cpy ctr.
   G4PWATotalXsecTable & operator=(const G4PWATotalXsecTable &right);
   G4PWATotalXsecTable(const G4PWATotalXsecTable&);

 private:
   // size of the table: Z=1-103 //
   static const G4int fgNumZet =  103;

   // G4PWATotalXsecZ pointers for Z=1-103 //
   static G4PWATotalXsecZ *fgPWATotalXsecTable[fgNumZet]; 
};

#endif

