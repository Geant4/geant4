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
// -------------------------------------------------------------------
//      GEANT 4 class file
//
//      CERN, Geneva, Switzerland
//
//      File name:     G4VGammaDeexcitation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
//
//      Creation date: 23 October 1998
//
//      Modifications:
//
//        18 October 2002, F. Lei
//          Added GetVaccantSN() and _vSN in order to link to ARM in low-e em
//          _vSN is updated in UpdateElectron()
//          Added SetVaccantSN(). It is need to to re-set _vSN after each 
//          IC happened.
// 
//        8 March 2002, Fan Lei (flei@space.qinetiq.com)
//          Added  SetEO () , GetEO(), UpdateElectrons() to allow the assignment 
//          and modification of electron configuration.
//          
//        21 Nov 2001, Fan Lei (flei@space.qinetiq.com)
//           Modified GenerateGamma() and UpdateUncleus() for implementation
//           of Internal Conversion processs
//
//        15 April 1999, Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)
//              Added creation time evaluation for products of evaporation
//
// -------------------------------------------------------------------

#ifndef G4VGAMMADEEXCITATION_HH
#define G4VGAMMADEEXCITATION_HH

#include "globals.hh"
#include "G4VGammaTransition.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4ElectronOccupancy.hh"

class G4VGammaDeexcitation {

public:

    G4VGammaDeexcitation();

    virtual ~G4VGammaDeexcitation();

    virtual G4VGammaTransition * CreateTransition() = 0;
    virtual G4bool CanDoTransition() const = 0;

    // Single gamma transition
    virtual G4FragmentVector * DoTransition();

    // Chain of gamma transitions
    virtual G4FragmentVector * DoChain();

    virtual G4Fragment * GenerateGamma();

    virtual const G4Fragment & GetNucleus() const;

    virtual void SetNucleus(const G4Fragment & nucleus);

    virtual void SetVerboseLevel(G4int verbose);

    void SetEO(G4ElectronOccupancy eo) { _electronO = eo; };
    void SetVaccantSN( G4int val ) { _vSN = val;};
  
    G4ElectronOccupancy GetEO() { return _electronO; };    
    G4int GetVacantSN() {return _vSN;};
protected:

    void Initialize();
    void UpdateNucleus(const G4Fragment * gamma);
    void UpdateElectrons ();
    void Update();

    G4VGammaTransition * _transition; // Owned pointer
    G4int _verbose;

private:

    G4Fragment _nucleus;
    G4ElectronOccupancy _electronO;
    G4int _vSN;

    G4VGammaDeexcitation(const G4VGammaDeexcitation & right);

    const G4VGammaDeexcitation & operator = (const G4VGammaDeexcitation & right);
    G4bool operator == (const G4VGammaDeexcitation & right) const;
    G4bool operator != (const G4VGammaDeexcitation & right) const;

};

#endif






