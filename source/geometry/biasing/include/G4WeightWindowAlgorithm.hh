// ----------------------------------------------------------------------
// Class G4WeightWindowAlgorithm
// see also G4VWeightWindowAlgorithm.
// 
// Class description:
//

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4WeightWindowAlgorithm_hh
#define G4WeightWindowAlgorithm_hh G4WeightWindowAlgorithm_hh 

#include "G4Nsplit_Weight.hh"
#include "G4VWeightWindowAlgorithm.hh"

class G4WeightWindowAlgorithm : public G4VWeightWindowAlgorithm
{

public:  // with description
  
  G4WeightWindowAlgorithm();
  
  ~G4WeightWindowAlgorithm(){}

  void SetUpperLimit(G4double Upper);
    // set upper limiting factor for window
    // - Upper is the maximum factor by which the weight of a particle may
    //   be higher than it should be according to the importance
  
  void SetLowerLimit(G4double Lower);
    // set lower limiting facotr for window
    // - Lower is the minimal factor by which the wight may be
    //   lower than it should be according to the importance



  G4Nsplit_Weight Calculate(G4double init_w, 
			    G4double importance) const;
    // calculate the number of tracks and their weight according 
    // to the upper and lower limmiting factors of the window
    // and the initial weight
    // - init_w is the initial weight of the particle
    // - importance is the importance of the cell

private:
  G4double fUpper;
  G4double fLower;
};

#endif
