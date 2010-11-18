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
/* standard reference physics list */
#include <QGSP_BERT.hh>
#endif
#include <G4PionMinus.hh>
#include "CexmcProductionModelFactory.hh"
#include "CexmcHadronicPhysics.hh"
#include "CexmcChargeExchangeProductionModel.hh"


#ifdef CEXMC_USE_QGSP_BIC_EMY
typedef QGSP_BIC_EMY                  CexmcBasePhysics;
#else
typedef QGSP_BERT                     CexmcBasePhysics;
#endif

typedef CexmcProductionModelFactory< CexmcBasePhysics,
                                     CexmcHadronicPhysics,
                                     CexmcChargeExchangeProductionModel >
                                      CexmcChargeExchangePMFactory;

typedef CexmcChargeExchangePMFactory  CexmcPMFactoryInstance;


#endif

