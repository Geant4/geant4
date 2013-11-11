/*
 * File:   G4ENDFYieldDataContainer.hh
 * Author: B. Wendt (wendbryc@isu.edu)
 *
 * Created on September 6, 2011, 10:19 AM
 */

#ifndef G4ENDFYIELDDATACONTAINER_HH
#define	G4ENDFYIELDDATACONTAINER_HH

#include "globals.hh"

#include "G4FFGEnumerations.hh"

/** G4ENDFYieldDataContainer is a simple data storage class that handles 
 *  the memory management internally. One instance stores the information
 *  for one fission product. In the event of a non-spontaneous fission, it
 *  can store the yield probabilities and errors at each of the fission-
 *  inducing particles energy levels. For ENDF data and neutron-induced
 *  fission these energies are typically 0.0253 eV, 1 MeV, and 5 MeV.
 */
class G4ENDFYieldDataContainer
{
public:
// Constructor
    /** Default constructor */
    G4ENDFYieldDataContainer ( G4int YieldSlots );

// Functions
    /** Get the meta state */
    G4FFGEnumerations::MetaState GetMetaState( void );
    /** Get the product */
    G4int GetProduct( void );
    /** Get the yield error */
    G4double* GetYieldError( void );
    /** Get the yield probability */
    G4double* GetYieldProbability( void );
    /** Get the number of yield slots */
    G4int GetYieldSlots( void );
    /** Set the meta state */
    void SetMetaState( G4FFGEnumerations::MetaState MetaState );
    /** Set the product */
    void SetProduct( G4int Product );
    /** Set the yield error */
    void SetYieldError( G4double* YieldError );
    /** Set the yield probability */
    void SetYieldProbability( G4double* YieldProbability );
    /** Set the number of yield slots */
    void SetYieldSlots( G4int NumberOfSlots );

protected:
// Data members
    /** The number of energy groups, or yield slots, that are stored */
    G4int YieldSlots_;
    /** ZZZAAA identifier of the stored isotope */
    G4int Product_;
    /** Metastable state information of the stored isotope */
    G4FFGEnumerations::MetaState MetaState_;
    /** Array of yield probabilities, one per yield slot */
    G4double* YieldProbability_;
    /** Array of the yield probability errors, one per yield slot */
    G4double* YieldError_;

// Destructor function(s)
public:
    /** Default deconstructor */
    ~G4ENDFYieldDataContainer( void );
};

#endif	/* G4ENDFYIELDDATACONTAINER_HH */

