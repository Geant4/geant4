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
// $Id: G4VRMLNetConfig.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// G4VRMLNetConfig.hh
// Satoshi Tanaka

#ifndef WIN32

#ifdef  G4VIS_BUILD_VRML_DRIVER

#ifndef G4VRML_NET_CONFIG_HH
#define G4VRML_NET_CONFIG_HH

const int   FR_VRML_DEFAULT_PORT     = 40801              ;
const char  FR_VRML_PORT_ENV     []  = "G4VRML_PORT"      ;
const char  FR_VRML_HOST_NAME_ENV[]  = "G4VRML_HOST_NAME" ;

#endif // G4VRML_NET_CONFIG_HH
#endif // G4VIS_BUILD_VRML_DRIVER
#endif //WIN32
