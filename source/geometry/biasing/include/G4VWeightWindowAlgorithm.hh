// ----------------------------------------------------------------------
// Class G4VWeightWindowAlgorithm
// This is a base class for an algorithm used for 
// weight window importance sampling. It calculates
// the new weight of particles and the number of copies
// it should be split into according to the initial weight 
// of the track and the importance of the cell. The window
// is defined by the values of Upper and Lower.
// 
// Class description:
//

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4VWeightWindowAlgorithm_hh
#define G4VWeightWindowAlgorithm_hh G4VWeightWindowAlgorithm_hh 

#include "G4Nsplit_Weight.hh"

class G4VWeightWindowAlgorithm
{

public:  // with description

  virtual ~G4VWeightWindowAlgorithm(){}

  virtual void SetUpperLimit(G4double Upper) = 0;
    // set upper limiting factor for window
  
  virtual void SetLowerLimit(G4double Lower) = 0;
    // set lower limiting facotr for window

  virtual G4Nsplit_Weight Calculate(G4double init_w, 
				    G4double importance) const = 0;
    // calculate the number of tracks and their weight according 
    // to the upper and lower limmiting factors of the window
    // and the initial weight
  
};

#endif
