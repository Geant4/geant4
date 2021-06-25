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
/// \file eventgenerator/pythia/pythia8decayer/include/Py8Decayer.hh
/// \brief Definition of the Py8Decayer class
///
/// \author J. Yarba; FNAL
///

#ifndef Py8Decayer_H
#define Py8Decayer_H

#include "G4VExtDecayer.hh"
#include "globals.hh"

#include "Pythia8/Pythia.h"

class G4Track;
class G4DecayProducts;

class Py8Decayer : public G4VExtDecayer
{
  
   public:

      //ctor & dtor
      Py8Decayer();
      virtual ~Py8Decayer();

      virtual G4DecayProducts* ImportDecayProducts(const G4Track&);
    
   private:
   
      // data members
      Pythia8::Pythia* fDecayer;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
