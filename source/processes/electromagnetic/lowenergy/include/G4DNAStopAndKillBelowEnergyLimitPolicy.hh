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
// $Id: G4DNAStopAndKillBelowEnergyLimitPolicy.hh,v 1.2 2006/06/29 19:35:23 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $

#ifndef   G4DNASTOPANDKILLBELOWENERGYLIMITPOLICY_HH
 #define  G4DNASTOPANDKILLBELOWENERGYLIMITPOLICY_HH 1
 
 #include "globals.hh"
 
 // EnergyLimitsPolicy must provide:
 //  - [protected] const double lowEnergyLimit

 template <typename EnergyLimitsPolicy>
 class G4DNAStopAndKillBelowEnergyLimitPolicy : public EnergyLimitsPolicy
 {
  protected:
                                        G4DNAStopAndKillBelowEnergyLimitPolicy() {}
                                       ~G4DNAStopAndKillBelowEnergyLimitPolicy() {}

   G4bool                               KillIncomingParticle(G4double k) const;
   void                                 BuildFinalStatesData(void) const {}

  private:
   // Hides default constructor and assignment operator as private 
                                        G4DNAStopAndKillBelowEnergyLimitPolicy(const G4DNAStopAndKillBelowEnergyLimitPolicy & copy);
   G4DNAStopAndKillBelowEnergyLimitPolicy & operator=(const G4DNAStopAndKillBelowEnergyLimitPolicy & right);
 };
 
 #include "G4DNAStopAndKillBelowEnergyLimitPolicy.icc"
#endif /* G4DNASTOPANDKILLBELOWENERGYLIMITPOLICY_HH */

