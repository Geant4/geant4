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
#ifndef G4QGSMParameters_h
#define G4QGSMParameters_h 1

#include "globals.hh"

class G4QGSMParameters 
    {
public:
      G4QGSMParameters();
      G4QGSMParameters(const G4QGSMParameters &right);
      ~G4QGSMParameters();
      
      int operator==(const G4QGSMParameters &right) const;
      int operator!=(const G4QGSMParameters &right) const;
private:     
     };
     
inline int G4QGSMParameters::operator==(const G4QGSMParameters &right) const
    {
    return 1;
    }

inline int G4QGSMParameters::operator!=(const G4QGSMParameters &right) const
    {
    return  0;
    }

#endif     
