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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef G4SoftStringBuilder_h
#define G4SoftStringBuilder_h 1

#include "globals.hh"
#include "G4KineticTrackVector.hh"
#include "G4ExcitedStringVector.hh"
#include "G4PartonPair.hh"


class G4SoftStringBuilder 
    {
public:
    G4SoftStringBuilder();
    G4SoftStringBuilder(const G4SoftStringBuilder &right);
    ~G4SoftStringBuilder();

    int operator==(const G4SoftStringBuilder &right) const;
    int operator!=(const G4SoftStringBuilder &right) const;

    G4ExcitedString* BuildString(G4PartonPair * aPair);        

private:     
    };
     
inline int G4SoftStringBuilder::operator==(const G4SoftStringBuilder &right) const
    {
    return 1;
    }

inline int G4SoftStringBuilder::operator!=(const G4SoftStringBuilder &right) const
    {
    return  0;
    }

#endif     
