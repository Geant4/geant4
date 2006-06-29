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
// $Id: Tst33VisEventAction.hh,v 1.5 2006-06-29 22:00:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class Tst33VisEventAction
//
// Class description:
//
// Event action used for visualization.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------


#ifndef Tst33VisEventAction_h
#define Tst33VisEventAction_h 1

#include "Tst33VEventAction.hh"
#include "globals.hh"

class Tst33VisEventActionMessenger;
class G4Trajectory;

class Tst33VisEventAction : public Tst33VEventAction
{
public:
  Tst33VisEventAction();
  virtual ~Tst33VisEventAction();

  virtual void   BeginOfEventAction(const G4Event*);
  virtual void   EndOfEventAction(const G4Event*);
  virtual void SpecialCellScorer(const G4CellScorer *scorer);
  virtual void Clear();
  void DrawTrajectory(G4Trajectory &trj) const;
  
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

    
