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
// $Id: exampleB03.cc,v 1.10 2002-11-18 13:22:48 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include "B03App.hh"
#include "globals.hh"

int main(int argc, char **argv)
{  

  B03AppBase &base = B03AppBase::GetB03AppBase();

  return 0;
}

