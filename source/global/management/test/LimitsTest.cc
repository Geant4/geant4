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
// $Id: LimitsTest.cc,v 1.2 2006-06-29 19:04:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
#include <iostream>
#include "templates.hh"

int main()
{
  std::cout << "DBL_MIN = " << DBL_MIN << std::endl
            << "DBL_DIG = " << DBL_DIG << std::endl
            << "DBL_MAX = " << DBL_MAX << std::endl
            << "DBL_EPSILON = " << DBL_EPSILON << std::endl
            << "FLT_MIN = " << FLT_MIN << std::endl
            << "FLT_DIG = " << FLT_DIG << std::endl
            << "FLT_MAX = " << FLT_MAX << std::endl
            << "FLT_EPSILON = " << FLT_EPSILON << std::endl;

  return 0;
}
