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
// $Id: Tst33VisEventAction.hh,v 1.2 2002-10-31 08:32:44 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Tst33VisEventAction_h
#define Tst33VisEventAction_h 1

#include "Tst33VEventAction.hh"
#include "globals.hh"

class Tst33VisEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Tst33VisEventAction : public Tst33VEventAction
{
public:
  Tst33VisEventAction();
  virtual ~Tst33VisEventAction();

  virtual void   BeginOfEventAction(const G4Event*);
  virtual void   EndOfEventAction(const G4Event*);
  virtual void SetCell_19_Scorer(const G4CellScorer *scorer);
  virtual void Clear();
  
  inline void SetDrawFlag   (const G4String &val);
  inline void SetPrintModulo(G4int    val);  
private:
  Tst33VisEventAction(const Tst33VisEventAction &);
  Tst33VisEventAction &operator=(const Tst33VisEventAction &);
  
  G4String                    drawFlag;
  G4int                       printModulo;                         
  Tst33VisEventActionMessenger*  eventMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void Tst33VisEventAction::SetDrawFlag   (const G4String &val){ 
  drawFlag = val;
}

inline void Tst33VisEventAction::SetPrintModulo(G4int    val) {
  printModulo = val;
}

#endif

    
