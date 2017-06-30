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
// $Id: G4Generator2BS.hh 104410 2017-05-30 07:17:09Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:  G4Generator2BS
//
// Author:     Andreia Trindade (andreia@lip.pt)
//             Pedro Rodrigues  (psilva@lip.pt)
//             Luis Peralta     (luis@lip.pt)
// 
// Creation date: 2 June 2003
//
// Modifications: 
// 02 Jun 2003  First implementation acording with new design
// 12 Oct 2010  V.Ivanchenko moved RejectionFunction inline
//               
//
// Class Description: 
//
// Concrete class for Bremsstrahlung Angular Distribution Generation 
// 2BS Distribution
//

// -------------------------------------------------------------------
//

#ifndef G4Generator2BS_h
#define G4Generator2BS_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4VEmAngularDistribution.hh"
#include "G4Log.hh"

class G4Pow;

class G4Generator2BS : public G4VEmAngularDistribution
{

public:

  explicit G4Generator2BS(const G4String& name="");

  virtual ~G4Generator2BS();

  virtual G4ThreeVector& SampleDirection(const G4DynamicParticle* dp,
                                         G4double out_energy,
                                         G4int Z,
                                         const G4Material* mat = nullptr);

  void PrintGeneratorInformation() const;

protected:

  inline G4double RejectionFunction(G4double value) const;

private:

  // hide assignment operator 
  G4Generator2BS & operator=(const  G4Generator2BS &right)= delete;
  G4Generator2BS(const  G4Generator2BS&) = delete;

  G4double fz;
  G4double ratio;
  G4double ratio1;
  G4double ratio2;
  G4double delta;

  G4Pow* g4pow;
  G4int  nwarn;

};

inline G4double G4Generator2BS::RejectionFunction(G4double y) const
{
  G4double y2 = (1 + y)*(1 + y);
  G4double x  = 4*y*ratio/y2;
  return 4*x - ratio1 - (ratio2 - x)*G4Log(delta + fz/y2); 
}

#endif

