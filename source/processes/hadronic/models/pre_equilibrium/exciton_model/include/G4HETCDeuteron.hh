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
// $Id: G4HETCDeuteron.hh 90337 2015-05-26 08:34:27Z gcosmo $
//
// by V. Lara
//
// Modified:
// 23.08.2010 V.Ivanchenko general cleanup, move constructor and destructor 
//            the source, use G4Pow
//

#ifndef G4HETCDeuteron_h
#define G4HETCDeuteron_h 1

#include "G4HETCChargedFragment.hh"
#include "G4ReactionProduct.hh"
#include "G4DeuteronCoulombBarrier.hh"

class G4HETCDeuteron : public G4HETCChargedFragment
{
public:

  G4HETCDeuteron();

  ~G4HETCDeuteron();
  
protected:

  virtual G4double GetAlpha() const;

  virtual G4double GetBeta() const;

  virtual G4double GetSpinFactor() const;

  virtual G4double K(const G4Fragment & aFragment);

private:

  // operators  
  G4HETCDeuteron(const G4HETCDeuteron &right);
  const G4HETCDeuteron & operator=(const G4HETCDeuteron &right); 
  G4bool operator==(const G4HETCDeuteron &right) const;
  G4bool operator!=(const G4HETCDeuteron &right) const;
    
private:
  
  G4DeuteronCoulombBarrier theDeuteronCoulombBarrier;

};

#endif
 

