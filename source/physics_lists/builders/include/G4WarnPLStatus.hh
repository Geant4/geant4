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
#ifndef  G4WarnPLStatus_hh
#define  G4WarnPLStatus_hh

#include "G4String.hh"

class  G4WarnPLStatus
{
public:
  G4WarnPLStatus();
  void Unsupported(const G4String aPL, const G4String Replacement ="") const ;
  void Experimental(const G4String aPL) const ;
};

inline 
G4WarnPLStatus::G4WarnPLStatus()
{}

inline
void G4WarnPLStatus::Unsupported(const G4String aPL, const G4String Replacement) const
{
    G4cout << 
"*=====================================================================" <<G4endl <<
"*                                                                     " <<G4endl <<
"*   The Physics list "<<aPL<<" is NO LONGER SUPPORTED !   " <<G4endl <<
"*   and is likely to be deleted in a future release of Geant4             " <<G4endl <<
"*                                                                     " <<G4endl;
   if (Replacement.size() > 0)
   {
   G4cout << 
"*    We recommend you try the physics lists "<<Replacement<< ","         <<G4endl <<
"*      this offers similar functionality for most use cases            " <<G4endl <<
"*                                                                      " <<G4endl;
   
   }
   G4cout << 
"*                                                                      " <<G4endl <<
"*   We invite you to report your use case for, and your experience with" <<G4endl <<
"*    this physics list on the Geant4 User Forum dedicated to physics   " <<G4endl <<
"*    lists:                                                            " <<G4endl <<
"*  http://geant4-hn.slac.stanford.edu:5090/HyperNews/public/get/phys-list.html"<<G4endl <<
"*                                                                      " <<G4endl <<
"*=====================================================================*" <<G4endl<<
G4endl;   
}
void G4WarnPLStatus::Experimental(const G4String aPL) const
{
    G4cout << 
"*=====================================================================" <<G4endl <<
"*                                                                     " <<G4endl <<
"*   The Physics list "<<aPL<<" is an experimental physics list !   " <<G4endl <<
"*                                                                      " <<G4endl <<
"*   Please  report your use case for, and your experience with this    " <<G4endl <<
"*    physics list on the Geant4 User Forum dedicated to physics lists: " <<G4endl <<
"*  http://geant4-hn.slac.stanford.edu:5090/HyperNews/public/get/phys-list.html"<<G4endl <<
"*                                                                      " <<G4endl <<
"*=====================================================================*" <<G4endl<<
G4endl;   
}

#endif
