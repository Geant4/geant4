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
#ifndef G4DiffractiveStringBuilder_h
#define G4DiffractiveStringBuilder_h 1

#include "globals.hh"
#include "G4KineticTrackVector.hh"
#include "G4ExcitedStringVector.hh"
#include "G4PartonPair.hh"

class G4DiffractiveStringBuilder 
    {
public:
      G4DiffractiveStringBuilder();
      G4DiffractiveStringBuilder(const G4DiffractiveStringBuilder &right);
      ~G4DiffractiveStringBuilder();
      
      G4int operator==(const G4DiffractiveStringBuilder &right) const;
      G4int operator!=(const G4DiffractiveStringBuilder &right) const;
      
      G4ExcitedString* BuildString(G4PartonPair* aParton);
      
private:     
     };
     
inline G4int G4DiffractiveStringBuilder::operator==(const G4DiffractiveStringBuilder &) const
    {
    return 1;
    }

inline G4int G4DiffractiveStringBuilder::operator!=(const G4DiffractiveStringBuilder &) const
    {
    return  0;
    }

#endif     
