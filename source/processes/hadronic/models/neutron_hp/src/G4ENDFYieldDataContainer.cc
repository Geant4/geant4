/*
 * File:   G4ENDFYieldDataContainer.cc
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on September 6, 2011, 10:19 AM
 */
 
#include "globals.hh"

#include "G4ArrayOps.hh"
#include "G4ENDFYieldDataContainer.hh"
#include "G4FFGEnumerations.hh"

G4ENDFYieldDataContainer::
G4ENDFYieldDataContainer( G4int YieldSlots )
{
    YieldSlots_ = YieldSlots;

    Product_ = 0;
    MetaState_ = G4FFGEnumerations::GROUND_STATE;
    YieldProbability_ = new G4double[YieldSlots_];
    YieldError_ = new G4double[YieldSlots_];
}

G4FFGEnumerations::MetaState G4ENDFYieldDataContainer::
GetMetaState( void )
{
    return MetaState_;
}

G4int G4ENDFYieldDataContainer::
GetProduct( void )
{
    return Product_;
}

G4double* G4ENDFYieldDataContainer::
GetYieldError( void )
{
    return YieldError_;
}

G4double* G4ENDFYieldDataContainer::
GetYieldProbability( void )
{
    return YieldProbability_;
}

G4int G4ENDFYieldDataContainer::
GetYieldSlots( void )
{
    return YieldSlots_;
}

void G4ENDFYieldDataContainer::
SetMetaState( G4FFGEnumerations::MetaState MetaState )
{
    MetaState_ = MetaState;
}

void G4ENDFYieldDataContainer::
SetProduct( G4int Product )
{
    Product_ = Product;
}

void G4ENDFYieldDataContainer::
SetYieldError( G4double* YieldError )
{
    G4ArrayOps::Copy(YieldSlots_, YieldError_, YieldError);
}

void G4ENDFYieldDataContainer::
SetYieldProbability( G4double* YieldProbability )
{
    G4ArrayOps::Copy(YieldSlots_, YieldProbability_, YieldProbability);
}

void G4ENDFYieldDataContainer::
//G4ENDFYieldDataContainer::SetYieldSlots(G4int NumberOfSlots)
SetYieldSlots(G4int NumberOfSlots)
{
    YieldSlots_ = NumberOfSlots;
}

G4ENDFYieldDataContainer::
~G4ENDFYieldDataContainer( void )
{
    delete[] YieldProbability_;
    delete[] YieldError_;
}

