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
// $Id: G4VExtDecayer.hh 105727 2017-08-16 12:47:05Z gcosmo $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//
// ------------------------------------------------------------
//  New  scheme                      23 Feb. 2001  H.Kurahige
// ------------------------------------------------------------
//
#ifndef G4VExtDecayer_h
#define G4VExtDecayer_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4DecayProducts.hh"
class G4Track;

class G4VExtDecayer 
{
 // Class Description
 //  This class is a Abstract class for external decayer
 // G4VExtDecayer has one pure virtual method of
 // ImportDecayProducts which return decay products  

  public: //With Description
    //  Constructors 
    G4VExtDecayer(const G4String& name ="");

    //  Destructor
    virtual ~G4VExtDecayer(){}

  private:
    //  copy constructor
      G4VExtDecayer(const G4VExtDecayer &){}

    //  Assignment Operation (generated)
      G4VExtDecayer & operator=(const G4VExtDecayer&){return *this;};

  public: //With Description
    virtual G4DecayProducts* ImportDecayProducts(
			         const G4Track& aTrack
                            ) = 0;

    const G4String& GetName() const;

  protected:
    G4String decayerName;
};

inline
 G4VExtDecayer::G4VExtDecayer(const G4String& name):
   decayerName(name)
{
}

inline
 const G4String& G4VExtDecayer::GetName() const
{
   return decayerName;
}

#endif










