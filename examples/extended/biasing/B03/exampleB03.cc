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
// $Id: exampleB03.cc,v 1.12 2006/06/29 16:34:59 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// 
// --------------------------------------------------------------
//      GEANT 4 - exampleB03
//
// --------------------------------------------------------------
// Comments
// This main may not be used at all. Instead you may use python or 
// lizard an execute the script B03RunApplication in them.
//
// This main function may used to test B03AppBase::GetB03AppBase().
// 
// --------------------------------------------------------------

#include "B03AppBase.hh"
#include "globals.hh"

int main(int argc, char **argv)
{  

  B03AppBase &base = B03AppBase::GetB03AppBase();

  return 0;
}

