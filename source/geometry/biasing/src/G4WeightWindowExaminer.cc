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
// $Id: G4WeightWindowExaminer.cc,v 1.2 2006/06/29 18:18:03 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4WeightWindowExaminer.cc
//
// ----------------------------------------------------------------------

#include "G4WeightWindowExaminer.hh"
#include "G4VParallelStepper.hh"
#include "G4VWeightWindowAlgorithm.hh"
#include "G4GeometryCellStep.hh"
#include "G4VWeightWindowStore.hh"

G4WeightWindowExaminer::
G4WeightWindowExaminer(const G4VWeightWindowAlgorithm &aWWalg,
		       const G4VParallelStepper &astepper,
		       const G4VWeightWindowStore &wwstore)
 : fWWalgorithm(aWWalg),
   fPStepper(astepper),
   fWWStore(wwstore)
{}

G4WeightWindowExaminer::~G4WeightWindowExaminer()
{}

  
G4Nsplit_Weight G4WeightWindowExaminer::Examine(G4double w, 
						G4double energy) const
{
  G4GeometryCellStep pstep = fPStepper.GetPStep();

  G4Nsplit_Weight nw = fWWalgorithm.
    Calculate(w, fWWStore.GetLowerWeitgh(pstep.GetPostGeometryCell(),
					 energy));
  
  return nw;
}
