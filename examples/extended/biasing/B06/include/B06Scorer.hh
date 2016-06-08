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
// $Id: B06Scorer.hh,v 1.3 2002/04/19 10:54:31 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01 $
//

#ifndef G4Pscorer_hh
#define G4Pscorer_hh G4Pscorer_hh

#include "g4std/iostream"

#include "G4VPScorer.hh"
#include "G4PMapPtkTallys.hh"

class G4Step;
class G4PStep;
class G4VIStore;

class B06Scorer : public G4VPScorer
{
public:
  B06Scorer(const G4VIStore &IStore);
  ~B06Scorer();
  void Score(const G4Step &aStep, const G4PStep &aPStep);
  const G4PMapPtkTallys &GetMapPtkTallys() const { return fPtkTallys; }
  G4bool CorrectWeight() const {return fCorrectWeight;}
  void ResetCorrectWeight() {fCorrectWeight = true;}
private:
  G4PMapPtkTallys fPtkTallys;
  const G4VIStore &fIStore; 
  G4bool fCorrectWeight;
};

G4std::ostream& operator<<(G4std::ostream &out, const B06Scorer &ps);

#endif
