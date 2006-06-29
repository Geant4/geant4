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
// $Id: G4LineSection.hh,v 1.9 2006-06-29 18:22:48 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// class G4LineSection
//
// Class description:
//
// A utility class that calculates the distance of a point from a 
// line section.

// History:
// - Created. J. Apostolakis.
// --------------------------------------------------------------------

#ifndef G4LineSection_hh
#define G4LineSection_hh

#include "G4Types.hh" 
#include "G4ThreeVector.hh"

class G4LineSection
{
  public:  // with description

     G4LineSection( const G4ThreeVector& PntA, const G4ThreeVector& PntB );

     G4double Dist( G4ThreeVector OtherPnt ) const;

     G4double GetABdistanceSq() const { return fABdistanceSq ; }

     static G4double Distline( const G4ThreeVector& OtherPnt, 
                               const G4ThreeVector& LinePntA, 
                               const G4ThreeVector& LinePntB );
  private:

     G4ThreeVector   EndpointA;
     G4ThreeVector   VecAtoB;

     G4double fABdistanceSq ;
};


#endif
