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

#include "G4WarnPLStatus.hh"

G4WarnPLStatus::G4WarnPLStatus()
{}

G4WarnPLStatus::~G4WarnPLStatus()
{}

void G4WarnPLStatus::Replaced(const G4String aPL, const G4String Replacement) const
{
    G4cout << 
"*=====================================================================" <<G4endl <<
"*                                                                     " <<G4endl <<
"*   The Physics list "<<aPL<<" no longer exists                       " <<G4endl <<
"*   We recommend you use the physics lists "<<Replacement<< ","         <<G4endl <<
"*      this offers similar functionality for most use cases            " <<G4endl <<
"*                                                                      " <<G4endl <<
"*                                                                      " <<G4endl <<
"*   We invite you to report your use case for, and your experience with" <<G4endl <<
"*    this physics list on the Geant4 User Forum dedicated to physics   " <<G4endl <<
"*    lists:                                                            " <<G4endl <<
"*  http://hypernews.slac.stanford.edu/HyperNews/geant4/get/phys-list.html"<<G4endl <<
"*                                                                      " <<G4endl <<
"*=====================================================================*" <<G4endl<<
G4endl;   
}

void G4WarnPLStatus::OnlyFromFactory(const G4String aPL, const G4String basePL) const
{
    G4cout << 
"*=====================================================================" <<G4endl <<
"*                                                                     " <<G4endl <<
"*   The Physics list "<<aPL<<", a variation of "<< basePL<< " will be " <<G4endl <<
"*      available only via the physics list factory starting from the  " <<G4endl <<
"*      next release, Geant4 10 .                                      " <<G4endl <<
"*   We recommend you to replace code like                             " <<G4endl <<
"*                                                                     " <<G4endl <<
"       runManager->SetUserInitialization( new " << aPL << " );        " <<G4endl <<
"*                                                                     " <<G4endl <<
"*   by the following                                                  " <<G4endl <<
"*                                                                     " <<G4endl <<
"       G4PhysListFactory factory;                                     " <<G4endl <<
"       runManager->SetUserInitialization("                              <<G4endl <<
"                      factory.GetReferencePhysList(\"" << aPL << "\");" <<G4endl <<
"*                                                                      " <<G4endl <<
"*   For more information how to use G4PhysListFactory, please refer    " <<G4endl <<
"*    to the documentation available at                                 " <<G4endl <<
"*     http://cern.ch/geant4/support/physicsLists/PhysListFactory.shtml " <<G4endl <<
"*                                                                      " <<G4endl <<
"*   We invite you to report your use case for, and your experience with" <<G4endl <<
"*    this physics list on the Geant4 User Forum dedicated to physics   " <<G4endl <<
"*    lists:                                                            " <<G4endl <<
"*  http://hypernews.slac.stanford.edu/HyperNews/geant4/get/phys-list.html"<<G4endl <<
"*                                                                      " <<G4endl <<
"*=====================================================================*" <<G4endl<<
G4endl;
    //G4String txtPL;
    //txtPL=aPL + "::" + aPL;
    //G4Exception(txtPL,"PhysicsLists001", FatalException,desc);
}


void G4WarnPLStatus::Unsupported(const G4String aPL, const G4String Replacement) const
{
    G4cout << 
"*=====================================================================" <<G4endl <<
"*                                                                     " <<G4endl <<
"*   The Physics list "<<aPL<<" is NO LONGER SUPPORTED !   " <<G4endl <<
//"*   and is likely to be deleted in a future release of Geant4             " <<G4endl <<
"*   and will be deleted in the next release, Geant4 10                " <<G4endl <<
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
"*  http://hypernews.slac.stanford.edu/HyperNews/geant4/get/phys-list.html"<<G4endl <<
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
"*  http://hypernews.slac.stanford.edu/HyperNews/geant4/get/phys-list.html"<<G4endl <<
"*                                                                      " <<G4endl <<
"*=====================================================================*" <<G4endl<<
G4endl;   
}
