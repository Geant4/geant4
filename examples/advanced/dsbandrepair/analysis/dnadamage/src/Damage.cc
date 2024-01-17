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
/// \file Damage.cc
/// \brief Implementation of the Damage class

#include "Damage.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

Damage::Damage(DamageType pType,unsigned int pChromo,unsigned int pEvt,unsigned int pStrand,unsigned long int pCopyNb,Position pPos,DamageCause pCause,DamageChromatin pChrom):
fType(pType),
fChromo(pChromo),
fEvt(pEvt),
fStrand(pStrand),
fCopyNb(pCopyNb),
fPosition(pPos),
fCause(pCause),
fChromatin(pChrom)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool Damage::operator != (const Damage& pDmg) const
{
  return (pDmg.fType != fType)&&(pDmg.fCopyNb != fCopyNb)||(pDmg.fStrand != fStrand)||(pDmg.fEvt != fEvt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool Damage::operator == (const Damage& pDmg) const
{
  return (pDmg.fType != fType)&&(pDmg.fCopyNb == fCopyNb)&&(pDmg.fStrand == fStrand)&&(pDmg.fEvt == fEvt);
}
