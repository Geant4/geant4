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

#include "G4RadioactiveDecayMode.hh"

std::istream &operator >> (std::istream& strm, G4RadioactiveDecayMode& q)
{
  G4String a;
  strm >> a;
  if (a == "IT")
    {q = IT;}
  else if (a == "BetaMinus")
    {q = BetaMinus;}
  else if (a == "BetaPlus")
    {q = BetaPlus;}
  else if (a == "KshellEC")
    {q = KshellEC;}
  else if (a == "LshellEC")
    {q = LshellEC;}
  else if (a == "MshellEC")
    {q = MshellEC;}
  else if (a == "Alpha")
    {q = Alpha;}
  else if (a == "Proton")
    {q = Proton;}
  else if (a == "Neutron")
    {q = Neutron;}
  else if (a == "BDProton")
    {q = BDProton;}
  else if (a == "BDNeutron")
    {q = BDNeutron;}
  else if (a == "Beta2Minus")
    {q = Beta2Minus;}
  else if (a == "Beta2Plus")
    {q = Beta2Plus;}
  else if (a == "Proton2")
    {q = Proton2;}
  else if (a == "Neutron2")
    {q = Neutron2;}
  else if (a == "SpFission")
    {q = SpFission;}
  else
    {q = RDM_ERROR;}
  return strm;
}

