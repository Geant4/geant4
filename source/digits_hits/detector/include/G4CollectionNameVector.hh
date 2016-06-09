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
// $Id: G4CollectionNameVector.hh,v 1.1 2003/10/03 10:07:48 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//

#ifndef G4CollectionNameVector_H 
#define G4CollectionNameVector_H 1

#include "globals.hh"
#include <vector>

class G4CollectionNameVector : public std::vector<G4String>
{
  public:
    G4CollectionNameVector() {;}
    virtual ~G4CollectionNameVector() {;}

    void insert(G4String s) { push_back(s); }
};

#endif

