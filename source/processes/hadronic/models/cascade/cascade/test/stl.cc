#include "../utils/globals.hh"  // iterator defined here

// STL example

int main(){
  const int SIZE = 10;
  vector <G4double> crossSection; 

  for(int i = 0; i < SIZE; i++) crossSection.insert(crossSection.begin(), static_cast<G4double>(UniformRand()));

  sort(crossSection.begin(), crossSection.end());

  cout << SIZE << " uniform random numbers sorted:" << endl;

  for (iterator i = crossSection.begin(); i < crossSection.end(); i++) cout << *i << " ";

  cout << endl;
  return 1;
}



