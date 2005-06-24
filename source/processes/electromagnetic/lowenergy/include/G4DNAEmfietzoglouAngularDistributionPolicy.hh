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
// $Id: G4DNAEmfietzoglouAngularDistributionPolicy.hh,v 1.1 2005-06-24 10:07:13 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
