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
// $Id: Tst33ScorerBuilder.hh,v 1.2 2002-10-29 16:37:09 dressel Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// Class Tst33ScorerBuilder
//
// Class description:
//
// ...

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------

#ifndef Tst33ScorerBuilder_hh
#define Tst33ScorerBuilder_hh Tst33ScorerBuilder_hh

class G4CellScorerStore;
class Tst33VGeometry;
class G4CellScorer;

class Tst33ScorerBuilder {
public:
  Tst33ScorerBuilder();
  ~Tst33ScorerBuilder();
  G4CellScorerStore *CreateScorer(Tst33VGeometry *samplegeo,
				  const G4CellScorer **specialCellScorer);  
};


#endif
