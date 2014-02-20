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
/*
 * =============================================================================
 *
 *       Filename:  CexmcBasicPhysicsSettings.hh
 *
 *    Description:  basic typedefs etc. to build studied physics
 *
 *        Version:  1.0
 *        Created:  28.11.2009 15:30:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Alexey Radkov (), 
 *        Company:  PNPI
 *
 * =============================================================================
 */

#ifndef CEXMC_BASIC_PHYSICS_SETTINGS_HH
#define CEXMC_BASIC_PHYSICS_SETTINGS_HH

#ifdef CEXMC_USE_QGSP_BIC_EMY
/* this reference physics list promises higher accuracy for electrons, hadrons
 * and ions tracking */
#include <QGSP_BIC_EMY.hh>
#else
#ifdef CEXMC_USE_QGSP_BERT
#include <QGSP_BERT.hh>
#else
/* standard reference physics list */
#include <FTFP_BERT.hh>
#endif
#endif
#include <G4PionMinus.hh>
#include "CexmcProductionModelFactory.hh"
#include "CexmcHadronicPhysics.hh"
#include "CexmcChargeExchangeProductionModel.hh"


#ifdef CEXMC_USE_QGSP_BIC_EMY
typedef QGSP_BIC_EMY                  CexmcBasePhysics;
#else
#ifdef CEXMC_USE_QGSP_BERT
typedef QGSP_BERT                     CexmcBasePhysics;
#else
typedef FTFP_BERT                     CexmcBasePhysics;
#endif
#endif

typedef CexmcProductionModelFactory< CexmcBasePhysics,
                                     CexmcHadronicPhysics,
                                     CexmcChargeExchangeProductionModel >
                                      CexmcChargeExchangePMFactory;

typedef CexmcChargeExchangePMFactory  CexmcPMFactoryInstance;


#endif

