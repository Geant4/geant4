// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VExtDecayer.hh,v 1.1 2001-02-22 13:29:26 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
//      For information related to this code contact:
//      CERN, CN Division, ASD group
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
    virtual ~G4VExtDecayer();

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










