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
// $Id: G4VImportanceAlgorithm.hh,v 1.3 2002-04-10 13:13:07 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4VImportanceAlgorithm
//
// Class description:
//
// This is an interface used by importance sampling to get the 
// number of copies and weight a mother particle should be split 
// into when crossing a boundary of "importance cells". 
// The interface defines the input to be the ratio of the pre over
// the post importance and the weight of the mother track.
// It returns a struct containing the number of copies (including 
// the mother track) to be produced and the weight of each track.
// A user defined algorithm deriving from this interface may be used
// by the importance sampling.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VImportanceAlgorithm_hh
#define G4VImportanceAlgorithm_hh G4VImportanceAlgorithm_hh 

#include "G4Nsplit_Weight.hh"

class G4VImportanceAlgorithm
{

public:  // with description

  virtual ~G4VImportanceAlgorithm(){}
  virtual G4Nsplit_Weight Calculate(G4double ipre_over_ipost,
                                    G4double init_w) const = 0;
    // calculate the number of tracks and their weight according 
    // to the ratio of pre over post importance and the weight of
    // the mother track.
};

#endif
