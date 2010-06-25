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
// GEANT4 physics class: G4QNeutronCaptureRatio -- header file
// M.V. Kossov, ITEP(Moscow), 22-May-09
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 22-May-2009
//
// **********************************************************************
// This Header is a part of the CHIPS physics package (author: M. Kosov)
//=======================================================================
// Short description: (n,gamma) capture is a part of the incoherent
// (inelastic) interaction. This part is calculated in the class.
// ----------------------------------------------------------------------

#ifndef G4QNeutronCaptureRatio_h
#define G4QNeutronCaptureRatio_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include <vector>
#include "G4QPDGCode.hh"
#include "G4QEnvironment.hh"
#include "G4Quasmon.hh"
#include "G4QHadronVector.hh"

class G4QNeutronCaptureRatio
{
 protected:

  G4QNeutronCaptureRatio()  {}                 // Constructor

 public:

  ~G4QNeutronCaptureRatio() {}                 // Destructor

  static G4QNeutronCaptureRatio* GetPointer(); // Gives a pointer to this singletone

  // Capture/Inelastic Ratio (momentum is in independent units)
  G4double GetRatio(G4double pIU, G4int tgZ, G4int tgN);

 private:
  // These working member functions are in CHIPS units and must not be used externally
  G4double CalcCap2In_Ratio(G4double p, G4int Z, G4int N); // R = Capture/In (p in GeV/c)

  // Body
 private:
};      
#endif
