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
// $Id: G4DNAEmfietzoglouAngularDistributionPolicy.hh,v 1.2 2006/06/29 19:34:06 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//
// Nucl. Instr. Meth. 155 (1978) 145-156
// J. Phys. D: Appl. Phys. 33 (2000) 932-944
// Phys. Med. Biol. 45 (2000) 3171-3194

#ifndef   G4DNAEMFIETZOGLOUANGULARDISTRIBUTIONPOLICY_HH
 #define  G4DNAEMFIETZOGLOUANGULARDISTRIBUTIONPOLICY_HH
 
 #include "globals.hh"

 class G4DNAEmfietzoglouAngularDistributionPolicy
 {
  protected:
                                        G4DNAEmfietzoglouAngularDistributionPolicy() {}
                                       ~G4DNAEmfietzoglouAngularDistributionPolicy() {}

   G4double                             RandomizeCosTheta(G4double k, G4int z) const;
   G4bool                               KillIncomingParticle(G4double /* k */) const { return false; }
   void                                 BuildFinalStatesData(void) const {}
   
  private:
   G4double                             ScreeningFactor(G4double k, G4int z) const;

   // Hides default constructor and assignment operator as private 
                                        G4DNAEmfietzoglouAngularDistributionPolicy(const G4DNAEmfietzoglouAngularDistributionPolicy & copy); 
   G4DNAEmfietzoglouAngularDistributionPolicy & operator=(const G4DNAEmfietzoglouAngularDistributionPolicy & right);
 };
#endif /* G4DNAEMFIETZOGLOUANGULARDISTRIBUTIONPOLICY_HH */
