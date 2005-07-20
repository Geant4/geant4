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
// $Id: G4DNAStopAndKillBelowEnergyLimitPolicy.hh,v 1.1 2005-07-20 10:01:54 capra Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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

