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
//
/// \file electromagnetic/TestEm10/include/Materials.hh
/// \brief Definition of the Materials class
//
//
// -------------------------------------------------------------
//      GEANT 4 class 
//
//      ---------- Materials-------
//    Originally Created in Test30 by Vladimir Ivanchenko, 12 March 2002 
// 
//    Modified for Test by V. Grichine, 30 Jan 2006 
//
//

#ifndef Materials_h
#define Materials_h 1

#include "globals.hh"

class G4Material;

class Materials
{
  public:
    Materials();
    ~Materials();

    static  Materials* GetInstance();
     
    G4Material* GetMaterial(const G4String&);     
                      
  private:

    void Initialise();

    static Materials* fgInstance;
};

#endif

 


