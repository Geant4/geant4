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
// 2nd iteration (refined model at low energy)
// Based on Paul & Sacher database:

// Further documentation available from http://www.ge.infn.it/geant4/lowE
// -------------------------------------------------------------------
// $Id: G4hShellCrossSectionDoubleExpData.hh,v 1.3 2006-06-29 19:38:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifndef G4hShellCrossSectionDoubleExpData_HH
#define G4hShellCrossSectionDoubleExpData_HH 1

#include "globals.hh"
#include <vector>
#include <map>

class G4hShellCrossSectionDoubleExpData
{
public:

  G4hShellCrossSectionDoubleExpData();

  ~G4hShellCrossSectionDoubleExpData();
 
  std::vector<std::vector<G4double>*> GetParam(G4int);

private:
  void FillVectorValuesEnergy();
  void FillVectorValuesPar1();
  void FillVectorValuesPar2(); 
  void FillParameterMapEnergy();
  void FillParameterMapPar1(); 
  void FillParameterMapPar2();

  inline void InitializeVectorEnergy(std::vector<G4double> &vectEnergy, G4double value) const;

  inline void InitializeVectorPar1(std::vector<G4double> &vect, G4double value1, G4double value2, G4double value3) const;

  inline void InitializeVectorPar2(std::vector<G4double> &vect, G4double value1, G4double value2, G4double value3, G4double value4, G4double value5) const;

  std::map<G4int,std::vector<G4double>*,std::less<G4int> > parameterMapEnergy;

  std::map<G4int,std::vector<G4double>*,std::less<G4int> > parameterMapPar1;

  std::map<G4int,std::vector<G4double>*,std::less<G4int> > parameterMapPar2;



  // Vector conteining energy
  std::vector<G4double> energy6C;
  std::vector<G4double> energy7N ;
  std::vector<G4double> energy8O ;
  std::vector<G4double> energy9F ;
  std::vector<G4double> energy10Ne;
  std::vector<G4double> energy11Na;
  std::vector<G4double> energy12Mg;
  std::vector<G4double> energy13Al;
  std::vector<G4double> energy14Si;
  std::vector<G4double> energy15P ;
  std::vector<G4double> energy16S ;
  std::vector<G4double> energy17Cl;
  std::vector<G4double> energy18Ar;
  std::vector<G4double> energy19K ;
  std::vector<G4double> energy20Ca;
  std::vector<G4double> energy21Sc;
  std::vector<G4double> energy22Ti;
  std::vector<G4double> energy23V ;
  std::vector<G4double> energy24Cr;
  std::vector<G4double> energy25Mn;
  std::vector<G4double> energy26Fe;
  std::vector<G4double> energy27Co;
  std::vector<G4double> energy28Ni;
  std::vector<G4double> energy29Cu;
  std::vector<G4double> energy30Zn;
  std::vector<G4double> energy31Ga;
  std::vector<G4double> energy32Ge;
  std::vector<G4double> energy33As;
  std::vector<G4double> energy34Se;
  std::vector<G4double> energy35Br;
  std::vector<G4double> energy36Kr;
  std::vector<G4double> energy37Rb;
  std::vector<G4double> energy38Sr;
  std::vector<G4double> energy39Y ;
  std::vector<G4double> energy40Zr;
  std::vector<G4double> energy41Nb;
  std::vector<G4double> energy42Mo;
  std::vector<G4double> energy43Tc;
  std::vector<G4double> energy44Ru;
  std::vector<G4double> energy45Rh;
  std::vector<G4double> energy46Pd;
  std::vector<G4double> energy47Ag;
  std::vector<G4double> energy48Cd;
  std::vector<G4double> energy49In;
  std::vector<G4double> energy50Sn;
  std::vector<G4double> energy51Sb;
  std::vector<G4double> energy52Te;
  std::vector<G4double> energy53I ;
  std::vector<G4double> energy54Xe;
  std::vector<G4double> energy55Cs;
  std::vector<G4double> energy56Ba;
  std::vector<G4double> energy57La;
  std::vector<G4double> energy58Ce;
  std::vector<G4double> energy59Pr;
  std::vector<G4double> energy60Nd;
  std::vector<G4double> energy61Pm;
  std::vector<G4double> energy62Sm;
  std::vector<G4double> energy63Eu;
  std::vector<G4double> energy64Gd;
  std::vector<G4double> energy65Tb;
  std::vector<G4double> energy66Dy;
  std::vector<G4double> energy67Ho;
  std::vector<G4double> energy68Er;
  std::vector<G4double> energy69Tm;
  std::vector<G4double> energy70Yb;
  std::vector<G4double> energy71Lu;
  std::vector<G4double> energy72Hf;
  std::vector<G4double> energy73Ta;
  std::vector<G4double> energy74W ;
  std::vector<G4double> energy75Re;
  std::vector<G4double> energy76Os;
  std::vector<G4double> energy77Ir;
  std::vector<G4double> energy78Pt;
  std::vector<G4double> energy79Au;
  std::vector<G4double> energy80Hg;
  std::vector<G4double> energy81Tl;
  std::vector<G4double> energy82Pb;
  std::vector<G4double> energy83Bi;
  std::vector<G4double> energy84Po;
  std::vector<G4double> energy85At;
  std::vector<G4double> energy86Rn;
  std::vector<G4double> energy87Fr;
  std::vector<G4double> energy88Ra;
  std::vector<G4double> energy89Ac;
  std::vector<G4double> energy90Th;
  std::vector<G4double> energy91Pa;
  std::vector<G4double> energy92U ;

  // Vector conteining parameters for low energy
  std::vector<G4double> parlow6C;
  std::vector<G4double> parlow7N ;
  std::vector<G4double> parlow8O ;
  std::vector<G4double> parlow9F ;
  std::vector<G4double> parlow10Ne;
  std::vector<G4double> parlow11Na;
  std::vector<G4double> parlow12Mg;
  std::vector<G4double> parlow13Al;
  std::vector<G4double> parlow14Si;
  std::vector<G4double> parlow15P ;
  std::vector<G4double> parlow16S ;
  std::vector<G4double> parlow17Cl;
  std::vector<G4double> parlow18Ar;
  std::vector<G4double> parlow19K ;
  std::vector<G4double> parlow20Ca;
  std::vector<G4double> parlow21Sc;
  std::vector<G4double> parlow22Ti;
  std::vector<G4double> parlow23V ;
  std::vector<G4double> parlow24Cr;
  std::vector<G4double> parlow25Mn;
  std::vector<G4double> parlow26Fe;
  std::vector<G4double> parlow27Co;
  std::vector<G4double> parlow28Ni;
  std::vector<G4double> parlow29Cu;
  std::vector<G4double> parlow30Zn;
  std::vector<G4double> parlow31Ga;
  std::vector<G4double> parlow32Ge;
  std::vector<G4double> parlow33As;
  std::vector<G4double> parlow34Se;
  std::vector<G4double> parlow35Br;
  std::vector<G4double> parlow36Kr;
  std::vector<G4double> parlow37Rb;
  std::vector<G4double> parlow38Sr;
  std::vector<G4double> parlow39Y ;
  std::vector<G4double> parlow40Zr;
  std::vector<G4double> parlow41Nb;
  std::vector<G4double> parlow42Mo;
  std::vector<G4double> parlow43Tc;
  std::vector<G4double> parlow44Ru;
  std::vector<G4double> parlow45Rh;
  std::vector<G4double> parlow46Pd;
  std::vector<G4double> parlow47Ag;
  std::vector<G4double> parlow48Cd;
  std::vector<G4double> parlow49In;
  std::vector<G4double> parlow50Sn;
  std::vector<G4double> parlow51Sb;
  std::vector<G4double> parlow52Te;
  std::vector<G4double> parlow53I ;
  std::vector<G4double> parlow54Xe;
  std::vector<G4double> parlow55Cs;
  std::vector<G4double> parlow56Ba;
  std::vector<G4double> parlow57La;
  std::vector<G4double> parlow58Ce;
  std::vector<G4double> parlow59Pr;
  std::vector<G4double> parlow60Nd;
  std::vector<G4double> parlow61Pm;
  std::vector<G4double> parlow62Sm;
  std::vector<G4double> parlow63Eu;
  std::vector<G4double> parlow64Gd;
  std::vector<G4double> parlow65Tb;
  std::vector<G4double> parlow66Dy;
  std::vector<G4double> parlow67Ho;
  std::vector<G4double> parlow68Er;
  std::vector<G4double> parlow69Tm;
  std::vector<G4double> parlow70Yb;
  std::vector<G4double> parlow71Lu;
  std::vector<G4double> parlow72Hf;
  std::vector<G4double> parlow73Ta;
  std::vector<G4double> parlow74W ;
  std::vector<G4double> parlow75Re;
  std::vector<G4double> parlow76Os;
  std::vector<G4double> parlow77Ir;
  std::vector<G4double> parlow78Pt;
  std::vector<G4double> parlow79Au;
  std::vector<G4double> parlow80Hg;
  std::vector<G4double> parlow81Tl;
  std::vector<G4double> parlow82Pb;
  std::vector<G4double> parlow83Bi;
  std::vector<G4double> parlow84Po;
  std::vector<G4double> parlow85At;
  std::vector<G4double> parlow86Rn;
  std::vector<G4double> parlow87Fr;
  std::vector<G4double> parlow88Ra;
  std::vector<G4double> parlow89Ac;
  std::vector<G4double> parlow90Th;
  std::vector<G4double> parlow91Pa;
  std::vector<G4double> parlow92U ;

  // Vector conteining parameters for high energy
  std::vector<G4double> parhigh6C;
  std::vector<G4double> parhigh7N ;
  std::vector<G4double> parhigh8O ;
  std::vector<G4double> parhigh9F ;
  std::vector<G4double> parhigh10Ne;
  std::vector<G4double> parhigh11Na;
  std::vector<G4double> parhigh12Mg;
  std::vector<G4double> parhigh13Al;
  std::vector<G4double> parhigh14Si;
  std::vector<G4double> parhigh15P ;
  std::vector<G4double> parhigh16S ;
  std::vector<G4double> parhigh17Cl;
  std::vector<G4double> parhigh18Ar;
  std::vector<G4double> parhigh19K ;
  std::vector<G4double> parhigh20Ca;
  std::vector<G4double> parhigh21Sc;
  std::vector<G4double> parhigh22Ti;
  std::vector<G4double> parhigh23V ;
  std::vector<G4double> parhigh24Cr;
  std::vector<G4double> parhigh25Mn;
  std::vector<G4double> parhigh26Fe;
  std::vector<G4double> parhigh27Co;
  std::vector<G4double> parhigh28Ni;
  std::vector<G4double> parhigh29Cu;
  std::vector<G4double> parhigh30Zn;
  std::vector<G4double> parhigh31Ga;
  std::vector<G4double> parhigh32Ge;
  std::vector<G4double> parhigh33As;
  std::vector<G4double> parhigh34Se;
  std::vector<G4double> parhigh35Br;
  std::vector<G4double> parhigh36Kr;
  std::vector<G4double> parhigh37Rb;
  std::vector<G4double> parhigh38Sr;
  std::vector<G4double> parhigh39Y ;
  std::vector<G4double> parhigh40Zr;
  std::vector<G4double> parhigh41Nb;
  std::vector<G4double> parhigh42Mo;
  std::vector<G4double> parhigh43Tc;
  std::vector<G4double> parhigh44Ru;
  std::vector<G4double> parhigh45Rh;
  std::vector<G4double> parhigh46Pd;
  std::vector<G4double> parhigh47Ag;
  std::vector<G4double> parhigh48Cd;
  std::vector<G4double> parhigh49In;
  std::vector<G4double> parhigh50Sn;
  std::vector<G4double> parhigh51Sb;
  std::vector<G4double> parhigh52Te;
  std::vector<G4double> parhigh53I ;
  std::vector<G4double> parhigh54Xe;
  std::vector<G4double> parhigh55Cs;
  std::vector<G4double> parhigh56Ba;
  std::vector<G4double> parhigh57La;
  std::vector<G4double> parhigh58Ce;
  std::vector<G4double> parhigh59Pr;
  std::vector<G4double> parhigh60Nd;
  std::vector<G4double> parhigh61Pm;
  std::vector<G4double> parhigh62Sm;
  std::vector<G4double> parhigh63Eu;
  std::vector<G4double> parhigh64Gd;
  std::vector<G4double> parhigh65Tb;
  std::vector<G4double> parhigh66Dy;
  std::vector<G4double> parhigh67Ho;
  std::vector<G4double> parhigh68Er;
  std::vector<G4double> parhigh69Tm;
  std::vector<G4double> parhigh70Yb;
  std::vector<G4double> parhigh71Lu;
  std::vector<G4double> parhigh72Hf;
  std::vector<G4double> parhigh73Ta;
  std::vector<G4double> parhigh74W ;
  std::vector<G4double> parhigh75Re;
  std::vector<G4double> parhigh76Os;
  std::vector<G4double> parhigh77Ir;
  std::vector<G4double> parhigh78Pt;
  std::vector<G4double> parhigh79Au;
  std::vector<G4double> parhigh80Hg;
  std::vector<G4double> parhigh81Tl;
  std::vector<G4double> parhigh82Pb;
  std::vector<G4double> parhigh83Bi;
  std::vector<G4double> parhigh84Po;
  std::vector<G4double> parhigh85At;
  std::vector<G4double> parhigh86Rn;
  std::vector<G4double> parhigh87Fr;
  std::vector<G4double> parhigh88Ra;
  std::vector<G4double> parhigh89Ac;
  std::vector<G4double> parhigh90Th;
  std::vector<G4double> parhigh91Pa;
  std::vector<G4double> parhigh92U ;
};

#endif
