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
 * File:   G4FissionFragmentGenerator.cc
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on May 11, 2011, 12:04 PM
 */
 
#include <iostream>
#include <vector>

#include "G4Ions.hh"
#include "globals.hh"
#include "G4HadFinalState.hh"
#include "G4Neutron.hh"

#include "G4FFGDebuggingMacros.hh"
#include "G4FFGDefaultValues.hh"
// Use a few select constant of CLHEP namespace 
using CLHEP::eV;
using CLHEP::keV;
using CLHEP::MeV;
using CLHEP::GeV;
#include "G4FFGEnumerations.hh"
#include "G4FPYNormalFragmentDist.hh"
#include "G4FPYBiasedLightFragmentDist.hh"
#include "G4FissionFragmentGenerator.hh"
#include "G4TableTemplate.hh"

G4FissionFragmentGenerator::
G4FissionFragmentGenerator( void )
{
    // Set the default verbosity
    Verbosity_ = G4FFGDefaultValues::Verbosity;
    
    // Initialize the class
    Initialize();
}

G4FissionFragmentGenerator::
G4FissionFragmentGenerator( G4int Verbosity )
{
    // Set the verbosity
    Verbosity_ = Verbosity;
    
    // Initialize the class
    Initialize();
}

void G4FissionFragmentGenerator::
Initialize( void )
{
G4FFG_FUNCTIONENTER__

    // Initialize the class descriptor variables to the default values. These
    // will be used unless the user redefines them.
    Isotope_ = G4FFGDefaultValues::Isotope;
    MetaState_ = G4FFGDefaultValues::MetaState;
    Cause_ = G4FFGDefaultValues::FissionCause;
    IncidentEnergy_ = G4FFGDefaultValues::ThermalNeutronEnergy;
    YieldType_ = G4FFGDefaultValues::YieldType;
    TernaryProbability_ = G4FFGDefaultValues::TernaryProbability;
    AlphaProduction_ = G4FFGDefaultValues::AlphaProduction;
    SamplingScheme_ = G4FFGDefaultValues::SamplingScheme;

    // No data class has been created yet
    YieldData_ = NULL;
    IsReconstructionNeeded_ = TRUE;

G4FFG_FUNCTIONLEAVE__
}

G4DynamicParticleVector* G4FissionFragmentGenerator::
G4GenerateFission( void )
{
G4FFG_FUNCTIONENTER__

    const G4HadProjectile Projectile(G4DynamicParticle(G4Neutron::Neutron(),
                                                       G4ThreeVector(0, 0, 0),
                                                       G4FFGDefaultValues::ThermalNeutronEnergy));

    // Call the overloaded function and generate 1 fission
    std::vector< G4DynamicParticleVector* > FissionEvent = G4GenerateFission(1, Projectile);
    G4DynamicParticleVector* Container = FissionEvent[0];

G4FFG_FUNCTIONLEAVE__
    return Container;
}

G4DynamicParticleVector* G4FissionFragmentGenerator::
G4GenerateFission( const G4HadProjectile& Projectile )
{
G4FFG_FUNCTIONENTER__

    // Call the overloaded function and generate 1 fission
    std::vector< G4DynamicParticleVector* > FissionEvent = G4GenerateFission(1, Projectile);
    G4DynamicParticleVector* const Container = FissionEvent[0];

G4FFG_FUNCTIONLEAVE__
    return Container;
}

const std::vector< G4DynamicParticleVector* > G4FissionFragmentGenerator::
G4GenerateFission( G4long NumberOfFissions,
                   const G4HadProjectile& Projectile )
{
G4FFG_FUNCTIONENTER__


    //TK Modified 131107
    //std::vector< G4DynamicParticleVector* > FissionEvents(NumberOfFissions);
    std::vector< G4DynamicParticleVector* > FissionEvents(0);

    if(Projectile.GetDefinition() == G4Neutron::Neutron())
    {
        if(IsReconstructionNeeded_ == TRUE)
        {
            // TODO Eliminate potential need for restructuring during run phase
            //InitializeFissionProductYieldClass();
        }

        for(G4long i = 0; i < NumberOfFissions; i++)
        {
            FissionEvents.push_back(YieldData_->G4GetFission());
            // FIXME Use particle momentum in balance equation
            // FissionEvents.push_back(YieldData_->G4GetFission(Projectile.Get4Momentum()));
        }
    } else
    {
        FissionEvents.push_back(NULL);
    }

G4FFG_FUNCTIONLEAVE__
    return FissionEvents;
}

G4Ions* G4FissionFragmentGenerator::
G4GenerateFissionProduct( void )
{
G4FFG_FUNCTIONENTER__

    if(IsReconstructionNeeded_ == TRUE)
    {
        // TODO Eliminate potential need for restructuring during run phase
        //InitializeFissionProductYieldClass();
    }
    
    G4Ions* Product = YieldData_->G4GetFissionProduct();

G4FFG_FUNCTIONLEAVE__
    return Product;
}

G4double G4FissionFragmentGenerator::
G4GetAlphaProduction( void )
{
G4FFG_FUNCTIONENTER__

G4FFG_FUNCTIONLEAVE__
    return AlphaProduction_;
}

G4double G4FissionFragmentGenerator::
G4GetTernaryProbability( void )
{
G4FFG_FUNCTIONENTER__

G4FFG_FUNCTIONLEAVE__
    return TernaryProbability_;
}

G4FFGEnumerations::FissionCause G4FissionFragmentGenerator::
G4GetCause( void )
{
G4FFG_FUNCTIONENTER__

G4FFG_FUNCTIONLEAVE__
    return Cause_;
}

G4double G4FissionFragmentGenerator::
G4GetIncidentEnergy( void )
{
G4FFG_FUNCTIONENTER__

G4FFG_FUNCTIONLEAVE__
    return IncidentEnergy_;
}

G4int G4FissionFragmentGenerator::
G4GetIsotope( void )
{
G4FFG_FUNCTIONENTER__

G4FFG_FUNCTIONLEAVE__
    return Isotope_;
}

G4FFGEnumerations::MetaState G4FissionFragmentGenerator::
G4GetMetaState( void )
{
G4FFG_FUNCTIONENTER__

G4FFG_FUNCTIONLEAVE__
    return MetaState_;
}

G4FFGEnumerations::FissionSamplingScheme G4FissionFragmentGenerator::
G4GetSamplingScheme( void )
{
G4FFG_FUNCTIONENTER__

G4FFG_FUNCTIONLEAVE__
    return SamplingScheme_;
}

G4FFGEnumerations::YieldType G4FissionFragmentGenerator::
G4GetYieldType()
{
G4FFG_FUNCTIONENTER__

G4FFG_FUNCTIONLEAVE__
    return YieldType_;
}

G4int G4FissionFragmentGenerator::
G4MakeIsotopeCode(G4int Z, G4int A, G4int M)
{
    // Sanity check;
    A %= 1000;
    Z %= 1000;
    M %= 10;

    return (A + Z * 1000) * 10 + M;
}

void G4FissionFragmentGenerator::
G4SetAlphaProduction( G4double WhatAlphaProduction )
{
G4FFG_FUNCTIONENTER__

    AlphaProduction_ = WhatAlphaProduction;
    if(YieldData_ != NULL)
    {
        YieldData_->G4SetAlphaProduction(AlphaProduction_);
    }

    if(Verbosity_ & G4FFGEnumerations::UPDATES)
    {
        G4FFG_SPACING__
        G4FFG_LOCATION__
        
        G4cout << " -- Alpha production set to " << AlphaProduction_ << G4endl;
    }

G4FFG_FUNCTIONLEAVE__
}

void G4FissionFragmentGenerator::
G4SetTernaryProbability( G4double WhatTernaryProbability )
{
G4FFG_FUNCTIONENTER__

    TernaryProbability_ = WhatTernaryProbability;
    if(YieldData_ != NULL)
    {
        YieldData_->G4SetTernaryProbability(TernaryProbability_);
    }

    if(Verbosity_ & G4FFGEnumerations::UPDATES)
    {
        G4FFG_SPACING__
        G4FFG_LOCATION__
        
        G4cout << " -- Ternary fission probability set to " << TernaryProbability_ << G4endl;
    }

G4FFG_FUNCTIONLEAVE__
}

void G4FissionFragmentGenerator::
G4SetCause( G4FFGEnumerations::FissionCause WhichCause )
{
G4FFG_FUNCTIONENTER__

    G4bool IsValidCause = (WhichCause == G4FFGEnumerations::SPONTANEOUS
                           || WhichCause == G4FFGEnumerations::NEUTRON_INDUCED );
    G4bool IsSameCause = (Cause_ == WhichCause);
    
    if(!IsSameCause && IsValidCause)
    {
        Cause_ = WhichCause;
        if(Cause_ == G4FFGEnumerations::SPONTANEOUS)
        {
            IncidentEnergy_ = 0;
        }
        IsReconstructionNeeded_ = TRUE;
    }
    
    if(Verbosity_ != G4FFGEnumerations::SILENT)
    {
        G4String CauseString;
        switch(WhichCause)
        {
            case G4FFGEnumerations::SPONTANEOUS:
                CauseString = "SPONTANEOUS";
                break;
            case G4FFGEnumerations::NEUTRON_INDUCED:
                CauseString = "NEUTRON_INDUCED";
                break;
            case G4FFGEnumerations::PROTON_INDUCED:
                CauseString = "PROTON_INDUCED";
                break;
            case G4FFGEnumerations::GAMMA_INDUCED:
                CauseString = "GAMMA_INDUCED";
                break;
        }
        
        if(Verbosity_ & G4FFGEnumerations::WARNING)
        {
            G4FFG_SPACING__
            G4FFG_LOCATION__
            
            if(IsValidCause)
            {
                if(IsSameCause && YieldData_ != NULL)
                {
                    G4cout << " -- Already set to use " << CauseString << " as the fission cause. Yield data class will not be reconstructed." << G4endl;
                } else if(YieldData_ == NULL)
                {
                    G4cout << " -- Yield data class not yet constructed. " << CauseString << " will be applied when it is constructed." << G4endl;
                }
            } else
            {
                G4cout << " -- Invalid cause of fission" << G4endl;
            }
        }
        
        if((Verbosity_ & G4FFGEnumerations::UPDATES)
           && IsValidCause)
        {
            G4FFG_SPACING__
            G4FFG_LOCATION__
            
            G4cout << " -- Fission cause set to " << CauseString << "." << G4endl;
        }
    }

G4FFG_FUNCTIONLEAVE__
}

void G4FissionFragmentGenerator::
G4SetIncidentEnergy( G4double WhatIncidentEnergy )
{
G4FFG_FUNCTIONENTER__

    if(Cause_ != G4FFGEnumerations::SPONTANEOUS)
    {
        IncidentEnergy_ = WhatIncidentEnergy;
        if(YieldData_ != NULL)
        {
            YieldData_->G4SetEnergy(IncidentEnergy_);
        }
    }
    
    if(Verbosity_ != G4FFGEnumerations::SILENT)
    {
        std::ostringstream EnergyString;
        if(IncidentEnergy_ / GeV > 1)
        {
            EnergyString << IncidentEnergy_ / GeV << " GeV";
        } else if(IncidentEnergy_ / MeV > 1)
        {
            EnergyString << IncidentEnergy_ / MeV << " MeV";
        } else if(IncidentEnergy_ / keV > 1)
        {
            EnergyString << IncidentEnergy_ / keV << " keV";
        } else
        {
            EnergyString << IncidentEnergy_ / eV << " eV";
        }
            
        if(Verbosity_ & G4FFGEnumerations::ENERGY_INFO
           || Verbosity_ & G4FFGEnumerations::WARNING)
        {
            if(Cause_ == G4FFGEnumerations::SPONTANEOUS && IncidentEnergy_ != 0)
            {
                G4FFG_SPACING__
                G4FFG_LOCATION__
                
                G4cout << " -- Cannot set a non-zero energy for spontaneous fission" << G4endl;
            } else if(YieldData_ == NULL)
            {
                G4FFG_SPACING__
                G4FFG_LOCATION__
                
                G4cout << " -- Yield data class not yet constructed. " << EnergyString.str() << " will be applied when it is constructed." << G4endl;
                
            }
        }
        
        if(Verbosity_ & G4FFGEnumerations::ENERGY_INFO
           || Verbosity_ & G4FFGEnumerations::UPDATES)
        {
            G4FFG_SPACING__
            G4FFG_LOCATION__
            
            G4cout << " -- Incident neutron energy set to " << EnergyString.str() << "." << G4endl;
        }
    }

G4FFG_FUNCTIONLEAVE__
}

void G4FissionFragmentGenerator::
G4SetIsotope( G4int WhichIsotope )
{
G4FFG_FUNCTIONENTER__

    G4bool IsSameIsotope = (Isotope_ == WhichIsotope);

    if(!IsSameIsotope)
    {
        Isotope_ = WhichIsotope;
        IsReconstructionNeeded_ = TRUE;
    }
    
    if(Verbosity_ != G4FFGEnumerations::SILENT)
    {
        if(Verbosity_ & G4FFGEnumerations::WARNING)
        {
            if(IsSameIsotope && YieldData_ != NULL)
            {
                G4FFG_SPACING__
                G4FFG_LOCATION__
                
                G4cout << " -- Isotope " << Isotope_ << " already in use. Yield data class will not be reconstructed." << G4endl;
            } else if(YieldData_ == NULL)
            {
                G4FFG_SPACING__
                G4FFG_LOCATION__
                
                G4cout << " -- Yield data class not yet constructed. The isotope will be set to " << Isotope_ << " when it is constructed." << G4endl;
            }
        }
        
        if(Verbosity_ & G4FFGEnumerations::UPDATES)
        {
            G4FFG_SPACING__
            G4FFG_LOCATION__
            
            G4cout << " -- Isotope set to " << Isotope_ << "." << G4endl;
        }
    }

G4FFG_FUNCTIONLEAVE__
}

void G4FissionFragmentGenerator::
G4SetMetaState( G4FFGEnumerations::MetaState WhichMetaState )
{
G4FFG_FUNCTIONENTER__

    G4bool IsValidMetaState = (WhichMetaState >= G4FFGEnumerations::MetaStateFirst
                               && WhichMetaState <= G4FFGEnumerations::MetaStateLast);
    G4bool IsSameMetaState = (MetaState_ == WhichMetaState);
    
    if(!IsSameMetaState && IsValidMetaState)
    {
        MetaState_ = WhichMetaState;
        IsReconstructionNeeded_ = TRUE;
    }
    
    if(Verbosity_ != G4FFGEnumerations::SILENT)
    {
        G4String MetaName;
        switch(MetaState_)
        {
            case G4FFGEnumerations::GROUND_STATE:
                MetaName = "GROUND_STATE";
                break;
                
            case G4FFGEnumerations::META_1:
                MetaName = "META_1";
                break;
                
            case G4FFGEnumerations::META_2:
                MetaName = "META_2";
                break;
        }
        
        if(Verbosity_ & G4FFGEnumerations::WARNING)
        {
            G4FFG_SPACING__
            G4FFG_LOCATION__

            std::ostringstream Temp;
            if(IsValidMetaState)
            {
                if(IsSameMetaState && YieldData_ != NULL)
                {
                    G4cout << " -- Already set to use " << MetaName << " as the metastable state. Yield data class will not be reconstructed" << G4endl;
                } else if(YieldData_ == NULL)
                {
                    G4cout << " -- Yield data class not yet constructed. " << MetaName << " will be applied when it is constructed." << G4endl;
                }
            } else
            {
                G4cout << " -- Invalid metastable state." << G4endl;
            }
        }
        
        if(Verbosity_ & G4FFGEnumerations::UPDATES
           && IsValidMetaState)
        {
            G4FFG_SPACING__
            G4FFG_LOCATION__
            
            G4cout << " -- Metastable state set to " << MetaName << "." << G4endl;
        }
    }

G4FFG_FUNCTIONLEAVE__
}

void G4FissionFragmentGenerator
::G4SetSamplingScheme(G4FFGEnumerations::FissionSamplingScheme NewScheme)
{
G4FFG_FUNCTIONENTER__

    G4bool IsValidScheme = (NewScheme >= G4FFGEnumerations::FissionSamplingSchemeFirst
                            && NewScheme <= G4FFGEnumerations::FissionSamplingSchemeLast);
    G4bool IsSameScheme = (NewScheme == SamplingScheme_);

    if(!IsSameScheme && IsValidScheme)
    {
        SamplingScheme_ = NewScheme;
        IsReconstructionNeeded_ = TRUE;
    }
    
    if(Verbosity_ != G4FFGEnumerations::SILENT)
    {
        G4String SchemeString;
        switch(SamplingScheme_)
        {
            case G4FFGEnumerations::NORMAL:
                SchemeString = "NORMAL";
                break;
                
            case G4FFGEnumerations::LIGHT_FRAGMENT:
                SchemeString = "LIGHT_FRAGMENT";
                break;

            default:
                SchemeString = "UNSUPPORTED";
                break;
        }
        
        if(Verbosity_ & G4FFGEnumerations::WARNING)
        {
            G4FFG_SPACING__
            G4FFG_LOCATION__
            
            if(IsValidScheme)
            {
                if(IsSameScheme && YieldData_ != NULL)
                {
                    G4cout << " -- Already set to use " << SchemeString << " as the sampling scheme. Yield data class will not be reconstructed." << G4endl;
                } else if(YieldData_ == NULL)
                {
                    G4cout << " -- Yield data class not yet constructed. " << SchemeString << " will be applied when it is constructed." << G4endl;
                }
            } else
            {
                G4cout << " -- Invalid sampling scheme." << G4endl;
            }
        }
        
        if((Verbosity_ & G4FFGEnumerations::UPDATES)
           && IsValidScheme)
        {
            G4FFG_SPACING__
            G4FFG_LOCATION__
            
            G4cout << " -- Sampling scheme set to " << SchemeString << "." << G4endl;
        }
    }

G4FFG_FUNCTIONLEAVE__
}

void G4FissionFragmentGenerator::
G4SetYieldType(G4FFGEnumerations::YieldType WhichYieldType)
{
G4FFG_FUNCTIONENTER__

    G4bool IsValidYieldType = (WhichYieldType == G4FFGEnumerations::INDEPENDENT
                               ||WhichYieldType == G4FFGEnumerations::CUMULATIVE);
    G4bool IsSameYieldType = (YieldType_ == WhichYieldType);

    if(!IsSameYieldType && IsValidYieldType)
    {
        YieldType_ = WhichYieldType;
        IsReconstructionNeeded_ = TRUE;
    }
    
    if(Verbosity_ != G4FFGEnumerations::SILENT)
    {
        G4String YieldString;
        switch((int)YieldType_)
        {
        case G4FFGEnumerations::INDEPENDENT:
            YieldString = "INDEPENDENT";
            break;
        
        case G4FFGEnumerations::SPONTANEOUS:
            YieldString = "SPONTANEOUS";
            break;

        default:
            YieldString = "UNSUPPORTED";
            break;
        }
        
        if(Verbosity_ & G4FFGEnumerations::WARNING)
        {
            G4FFG_SPACING__
            G4FFG_LOCATION__
            
            if(IsValidYieldType)
            {
                
                if(IsSameYieldType && YieldData_ != NULL)
                {
                } else if(YieldData_ == NULL)
                {
                    G4cout << " -- Yield data class not yet constructed. Yield type " << YieldString << " will be applied when it is constructed." << G4endl;
                }
            } else
            {
                G4cout << " -- Invalid yield type." << G4endl;
            }
        }
        
        if((Verbosity_ & G4FFGEnumerations::UPDATES)
           && IsValidYieldType)
        {
            G4FFG_SPACING__
            G4FFG_LOCATION__
            
            G4cout << " -- Yield type set to " << YieldString << G4endl;
        }
    }

G4FFG_FUNCTIONLEAVE__
}

void G4FissionFragmentGenerator::
G4SetVerbosity(G4int Verbosity)
{
G4FFG_FUNCTIONENTER__

    Verbosity_ = Verbosity;
    
    if(YieldData_ != NULL)
    {
        YieldData_->G4SetVerbosity(Verbosity_);
    }

G4FFG_FUNCTIONLEAVE__
}

bool G4FissionFragmentGenerator::
InitializeFissionProductYieldClass( std::istringstream& dataStream )
{
G4FFG_FUNCTIONENTER__

    if(YieldData_ != NULL)
    {
        delete YieldData_;
    
        if(Verbosity_ & G4FFGEnumerations::UPDATES)
        {
            G4FFG_SPACING__
            G4FFG_LOCATION__
            
            G4cout << " -- Old yield data class deleted." << G4endl;
        }
    }

    try
    {
        if(SamplingScheme_ == G4FFGEnumerations::NORMAL)
        {
            YieldData_ = new G4FPYNormalFragmentDist(Isotope_,
                                                     MetaState_,
                                                     Cause_,
                                                     YieldType_,
                                                     Verbosity_,
                                                     dataStream);
        } else
        {
            YieldData_ = new G4FPYBiasedLightFragmentDist(Isotope_,
                                                          MetaState_,
                                                          Cause_,
                                                          YieldType_,
                                                          Verbosity_,
                                                          dataStream);
        }

        if(AlphaProduction_ != 0 && TernaryProbability_ != 0)
        {
            YieldData_->G4SetTernaryProbability(TernaryProbability_);
            YieldData_->G4SetAlphaProduction(AlphaProduction_);
        }

        if(Verbosity_ & G4FFGEnumerations::UPDATES)
        {
            G4FFG_SPACING__
            G4FFG_LOCATION__

            G4cout << " -- Yield data class constructed with defined values." << G4endl;
        }
    } catch (std::exception& e)
    {
        YieldData_ = NULL;
    }

    IsReconstructionNeeded_ = FALSE;

G4FFG_FUNCTIONLEAVE__
    return YieldData_;
}

G4FissionFragmentGenerator::
~G4FissionFragmentGenerator()
{
G4FFG_FUNCTIONENTER__

    delete YieldData_;

G4FFG_FUNCTIONLEAVE__
}
