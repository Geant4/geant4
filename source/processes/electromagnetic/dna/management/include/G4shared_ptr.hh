// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, so do not hesitate to send us your feedback!
//
// In order for Geant4-DNA to be maintained and still open-source, article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, in addition to the general paper on Geant4-DNA:
//
// The Geant4-DNA project, S. Incerti et al., Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we ask that you please cite the following papers reference papers on chemistry:
//
// Diﬀusion-controlled reactions modelling in Geant4-DNA, M. Karamitros et al., 2014 (submitted)
// Modeling Radiation Chemistry in the Geant4 Toolkit, M. Karamitros et al., Prog. Nucl. Sci. Tec. 2 (2011) 503-508

#ifndef G4SHARED_PTR_HH_
#define G4SHARED_PTR_HH_

#include "CLHEP/Utility/memory.h"

namespace G4
{
 using CLHEP::shared_ptr;
 using CLHEP::weak_ptr;
 using CLHEP::enable_shared_from_this;
 using CLHEP::static_pointer_cast;
 using CLHEP::const_pointer_cast;
 using CLHEP::dynamic_pointer_cast;
 //using CLHEP::polymorphic_cast_tag;
}


#endif /* G4SHARED_PTR_HH_ */
