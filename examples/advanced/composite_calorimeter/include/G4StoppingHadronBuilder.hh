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
#include "globals.hh"
#include "G4ios.hh"
#include "g4std/vector"


// At rest processes
#include "G4AntiProtonAnnihilationAtRest.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"
#include "G4PionMinusAbsorptionAtRest.hh"
#include "G4KaonMinusAbsorption.hh"

class G4StoppingHadronBuilder 
{
  public: 
    G4StoppingHadronBuilder();
    virtual ~G4StoppingHadronBuilder();

  public: 
    virtual void Build();

  private:

   // Pi -
   G4PionMinusAbsorptionAtRest thePionMinusAbsorption;

   // K -
   G4KaonMinusAbsorption theKaonMinusAbsorption;
 
   // anti-proton
   G4AntiProtonAnnihilationAtRest  theAntiProtonAnnihilation;
    
   // anti-neutron
   G4AntiNeutronAnnihilationAtRest  theAntiNeutronAnnihilation;
         
};
// 2002 by J.P. Wellisch


#endif

