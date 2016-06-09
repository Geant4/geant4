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
// $Id: G4QGSMFragmentation.hh,v 1.1 2003/10/07 11:25:40 hpw Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class implementation file
//
//      History: first implementation, Maxim Komogorov, 10-Jul-1998
// -----------------------------------------------------------------------------
#ifndef G4QGSMFragmentation_h
#define G4QGSMFragmentation_h 1

#include "G4VLongitudinalStringDecay.hh"

//******************************************************************************
class G4QGSMFragmentation:public G4VLongitudinalStringDecay
   {
public:
      G4QGSMFragmentation();
      ~G4QGSMFragmentation();

      const G4QGSMFragmentation & operator=(const G4QGSMFragmentation &right);
      int operator==(const G4QGSMFragmentation &right) const;
      int operator!=(const G4QGSMFragmentation &right) const;

  private:
   virtual G4double GetLightConeZ(G4double zmin, G4double zmax, G4int PartonEncoding,  G4ParticleDefinition* pHadron, G4double Px, G4double Py);      
   G4QGSMFragmentation(const G4QGSMFragmentation &right);

  private:
    // model parameters
    const G4double arho; 
    const G4double aphi;  
    const G4double an; 
    const G4double ala;  
    const G4double aksi; 
    const G4double alft;

  };

// Class G4QGSMFragmentation 
#endif


