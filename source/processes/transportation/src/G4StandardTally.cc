#include "G4StandardTally.hh"
#include "G4Sigma.hh"

G4StandardTally::G4StandardTally(const G4String &tallyname,
				 const G4String &rawtallyname,
				 const G4String &SigmaSpec) : 
  fName(tallyname),
  fRawTallyName(rawtallyname),
  fValue(0.),
  fHasTallied(false),
  fSigmaSpec(SigmaSpec)
{}

void G4StandardTally::Reset(){
  fValue= 0.;
  fHasTallied = false;
}

void G4StandardTally::Tally(const G4String &rawtallyname, 
			    G4Sigma &sigma){
  if (rawtallyname==fRawTallyName) {
    if (!fHasTallied) {
      fValue = ReadSigma(sigma, fSigmaSpec);
    }
    else {
      Error("Tally: already tallied!");
    }
    fHasTallied = true;
  }
}

G4String G4StandardTally::GetName(){
  return fName;
}
G4double G4StandardTally::GetValue(){
  return fValue;
}
G4bool G4StandardTally::HasTallied(){
  return fHasTallied;
}

G4double G4StandardTally::
ReadSigma(G4Sigma &sigma, const G4String &sigspec){
  if (sigspec=="Mean") {
    return sigma.GetMean();
  }
  else if (sigspec=="Sigma") {
    return sigma.GetSigma();
  }
  else if (sigspec=="Entries") {
    return sigma.GetEntries();
  }
  else if (sigspec=="Xsum") {
    return sigma.GetXsum();
  }
  else if (sigspec=="XXsum") {
    return sigma.GetXXsum();
  }
  else if (sigspec=="SumOfWeights") {
    return sigma.GetSumOfWeights();
  }
  else if (sigspec=="WeightedXsum") {
    return sigma.GetWeightedXsum();
  }
  else if (sigspec=="WeightedXXsum") {
    return sigma.GetWeightedXXsum();
  }
  else {
    Error("ReadSigma: can't read sigmas: Get" + sigspec + "() function");
    return -1;
  }
};
