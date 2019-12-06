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
/// \file G4String.hh
/// \brief Definition of the place-holder for the G4String class
//
/*
  This file only serves as a place-holder for the annotations of G4String class,
  which is a type alias in Geant4. 
  
  Only the marshaling code (MarshaledG4String.hh) generated from this file
   should be used in your code. Please make sure that this G4String.hh
   itself is NOT included in your code. Otherwise, it will override Geant4's
   valid G4String class and bad things would happen.

  vietha 2003.05.08
*/

#ifndef __G4String
#define __G4String

//MSH_BEGIN
class G4String : public std::string
{
    int dummy; /*MSH: manual
    { memcpy($$, param->c_str(), param->size());
    *($$+param->size()) = '\0'; 
    }
    { G4String* s = new G4String($$);
       memcpy(param, s, sizeof(G4String));
     }
     { int size = param->size()+1;
       while(size%8) size++;
       $SIZE = size; } */
};
//MSH_END
#endif
