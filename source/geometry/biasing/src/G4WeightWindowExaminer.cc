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
// $Id: G4WeightWindowExaminer.cc,v 1.1 2003/08/19 15:17:26 dressel Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
