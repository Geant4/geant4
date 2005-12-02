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
// $Id: LimitsTest.cc,v 1.1 2005-12-02 13:42:16 gcosmo Exp $
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
