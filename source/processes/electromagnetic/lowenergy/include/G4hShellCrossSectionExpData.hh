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
// Author: Simona Saliceti (Simona.Saliceti@ge.infn.it)
//
// History:
// -----------
// 22 Apr 2004 First committed to cvs
//
// -------------------------------------------------------------------
// Class description:
// Low Energy Electromagnetic Physics
// Parameterisation of cross sections for proton ionisation, K shell
// 1st iteration
// Further documentation available from http://www.ge.infn.it/geant4/lowE
// -------------------------------------------------------------------
// $Id: G4hShellCrossSectionExpData.hh,v 1.3 2006-06-29 19:38:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4hShellCrossSectionExpData_HH
#define G4hShellCrossSectionExpData_HH 1

#include "globals.hh"
#include <vector>
#include <map>

class G4hShellCrossSectionExpData
{
public:

  G4hShellCrossSectionExpData();

  ~G4hShellCrossSectionExpData();
 
  std::vector<G4double>* GetParam(G4int);

private:

  void FillVectorValues();
  void FillParameterMap();

  inline void InitializeVector(std::vector<G4double> &vect, G4double value1, G4double value2, G4double value3) const;

std::map<G4int,std::vector<G4double>*,std::less<G4int> > parameterMap;

std::vector<G4double> parameter6C;
std::vector<G4double> parameter7N ;
std::vector<G4double> parameter8O ;
std::vector<G4double> parameter9F ;
std::vector<G4double> parameter10Ne;
std::vector<G4double> parameter11Na;
std::vector<G4double> parameter12Mg;
std::vector<G4double> parameter13Al;
std::vector<G4double> parameter14Si;
std::vector<G4double> parameter15P ;
std::vector<G4double> parameter16S ;
std::vector<G4double> parameter17Cl;
std::vector<G4double> parameter18Ar;
std::vector<G4double> parameter19K ;
std::vector<G4double> parameter20Ca;
std::vector<G4double> parameter21Sc;
std::vector<G4double> parameter22Ti;
std::vector<G4double> parameter23V ;
std::vector<G4double> parameter24Cr;
std::vector<G4double> parameter25Mn;
std::vector<G4double> parameter26Fe;
std::vector<G4double> parameter27Co;
std::vector<G4double> parameter28Ni;
std::vector<G4double> parameter29Cu;
std::vector<G4double> parameter30Zn;
std::vector<G4double> parameter31Ga;
std::vector<G4double> parameter32Ge;
std::vector<G4double> parameter33As;
std::vector<G4double> parameter34Se;
std::vector<G4double> parameter35Br;
std::vector<G4double> parameter36Kr;
std::vector<G4double> parameter37Rb;
std::vector<G4double> parameter38Sr;
std::vector<G4double> parameter39Y ;
std::vector<G4double> parameter40Zr;
std::vector<G4double> parameter41Nb;
std::vector<G4double> parameter42Mo;
std::vector<G4double> parameter43Tc;
std::vector<G4double> parameter44Ru;
std::vector<G4double> parameter45Rh;
std::vector<G4double> parameter46Pd;
std::vector<G4double> parameter47Ag;
std::vector<G4double> parameter48Cd;
std::vector<G4double> parameter49In;
std::vector<G4double> parameter50Sn;
std::vector<G4double> parameter51Sb;
std::vector<G4double> parameter52Te;
std::vector<G4double> parameter53I ;
std::vector<G4double> parameter54Xe;
std::vector<G4double> parameter55Cs;
std::vector<G4double> parameter56Ba;
std::vector<G4double> parameter57La;
std::vector<G4double> parameter58Ce;
std::vector<G4double> parameter59Pr;
std::vector<G4double> parameter60Nd;
std::vector<G4double> parameter61Pm;
std::vector<G4double> parameter62Sm;
std::vector<G4double> parameter63Eu;
std::vector<G4double> parameter64Gd;
std::vector<G4double> parameter65Tb;
std::vector<G4double> parameter66Dy;
std::vector<G4double> parameter67Ho;
std::vector<G4double> parameter68Er;
std::vector<G4double> parameter69Tm;
std::vector<G4double> parameter70Yb;
std::vector<G4double> parameter71Lu;
std::vector<G4double> parameter72Hf;
std::vector<G4double> parameter73Ta;
std::vector<G4double> parameter74W ;
std::vector<G4double> parameter75Re;
std::vector<G4double> parameter76Os;
std::vector<G4double> parameter77Ir;
std::vector<G4double> parameter78Pt;
std::vector<G4double> parameter79Au;
std::vector<G4double> parameter80Hg;
std::vector<G4double> parameter81Tl;
std::vector<G4double> parameter82Pb;
std::vector<G4double> parameter83Bi;
std::vector<G4double> parameter84Po;
std::vector<G4double> parameter85At;
std::vector<G4double> parameter86Rn;
std::vector<G4double> parameter87Fr;
std::vector<G4double> parameter88Ra;
std::vector<G4double> parameter89Ac;
std::vector<G4double> parameter90Th;
std::vector<G4double> parameter91Pa;
std::vector<G4double> parameter92U ;
};

#endif
