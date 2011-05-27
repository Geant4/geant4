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
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4QDiffractionRatio -- header file
// M.V. Kossov, ITEP(Moscow), 24-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 10-Nov-2009
//
// --------------------------------------------------------------------
// Short description: Difraction excitation is a part of the incoherent
// (inelastic) interaction. This part is calculated in the class.
// --------------------------------------------------------------------

#ifndef G4QDiffractionRatio_h
#define G4QDiffractionRatio_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include <vector>
#include "G4QPDGCode.hh"
#include "G4QEnvironment.hh"
#include "G4Quasmon.hh"
#include "G4QHadronVector.hh"
#include "G4RandomDirection.hh"


class G4QDiffractionRatio
{
 protected:

  G4QDiffractionRatio()  {}                 // Constructor

 public:

  ~G4QDiffractionRatio() {}                 // Destructor

  static G4QDiffractionRatio* GetPointer(); // Gives a pointer to this singletone

  // Diffraction/Prodaction Ratio (Production=Inelastic-QuasiElastic)
  G4double GetRatio(G4double pIU, G4int prPDG, G4int tgZ, G4int tgN);

  // ==> The following will be a protected function for internal CHIPS usage
  // ProjFragment(pPDG,p4M) on a nucleus (tgZ, tgN), result: Vector of secondary hadrons
  // Whoever uses this member function is responsible for DEL/DESTROY of G4QHadronVector
  G4QHadronVector* ProjFragment(G4int pPDG, G4LorentzVector p4M, G4int tgZ, G4int tgN);

  // ==> The following will be a protected function for internal CHIPS usage
  // TargFragment(pPDG,p4M) on a nucleus (tgZ, tgN), result: Vector of secondary hadrons
  // Whoever uses this member function is responsible for DEL/DESTROY of G4QHadronVector
  G4QHadronVector* TargFragment(G4int pPDG, G4LorentzVector p4M, G4int tgZ, G4int tgN);

  // Single Diffraction Target Excitation Cross-Section (independent Units)
  G4double GetTargSingDiffXS(G4double pIU, G4int prPDG, G4int tgZ, G4int tgN);

 private:
  // These working member functions are in CHIPS units and must not be used externally
  G4double CalcDiff2Prod_Ratio(G4double s, G4int A); // R=Diff/Prod (sqrt(s)(hN) in GeV, A)

 // Body
 private:
};      
#endif
