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
// $Id: G4ImportanceSampler.cc,v 1.2 2002-04-09 16:23:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ImportanceSampler.cc
//
// ----------------------------------------------------------------------

#include "G4ImportanceSampler.hh"
#include "G4ImportanceFinder.hh"
#include "G4VParallelStepper.hh"
#include "G4VImportanceAlgorithm.hh"

G4ImportanceSampler::
G4ImportanceSampler(const G4VImportanceAlgorithm &aIalg,
		    const G4VParallelStepper &astepper,
		    const G4VIStore &istore)
 : fIalgorithm(aIalg),
   fPStepper(astepper),
   fIfinder(*(new G4ImportanceFinder(istore)))
{}

G4ImportanceSampler::~G4ImportanceSampler()
{
  delete &fIfinder;
}
  
G4Nsplit_Weight G4ImportanceSampler::Sample(G4double w) const
{
  G4PStep pstep = fPStepper.GetPStep();
  return fIalgorithm.
    Calculate(fIfinder.GetIPre_over_IPost(pstep.fPreTouchableKey,
					  pstep.fPostTouchableKey), 
	      w);
}
