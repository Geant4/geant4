#include "../cascade/src/G4NucleusModel.cc"
#include "../cascade/src/G4VRegionModel.cc"
#include "../cascade/src/G4RegionModel.cc"
#include "../cascade/src/G4BertiniData.cc"

//test method for G4NucleusModel::ApplyRecoil()
//One proton and one neutron are expelled from nucleus at each round.

int main(){
  vector<G4KineticTrack*> trackVector;
  G4KineticTrack *track1, *track2;
  G4NucleusModel testNucleus;

  ifstream in;
  G4int A, Z, pix, piy, piz, p1x, p1y, p1z, p2x, p2y, p2z, n;
  char choice;

  in.open("recoil.in", ios::in);
  in >> A >> Z >> pix >> piy >> piz >> p1x >> p1y >> p1z >> p2x >> p2y >> p2z >> n >> choice;

   testNucleus.CreateModel(A,Z); 
   track1 = new G4KineticTrack();
   track2 = new G4KineticTrack();

   G4LorentzVector p1(p1x, p1y, p1z, 0.0);
   G4LorentzVector p2(p2x, p2y, p2z, 0.0);
   G4LorentzVector pInit(pix, piy, piz, 0.0);

   //cout << p1 << endl;
   //cout << p2 << endl;

   track1->Set4Momentum(p1);
   track2->Set4Momentum(p2);

   trackVector.push_back(track1);
   trackVector.push_back(track2);

   testNucleus.SetMomentum(pInit);

   if(choice=='e'){
     cout <<"0\t" << testNucleus.GetKineticEnergy() << endl;
    for(int i=0; i < n; i++){ 
      if( i >= (Z - 1)) break;
      testNucleus.ApplyRecoil(trackVector);
      cout << i+1 <<"\t"<< testNucleus.GetKineticEnergy() << endl;
      testNucleus.ChangeParameters(-2, -1); //reduce nucleus mass
    }
  }
  if(choice=='m'){
    G4ThreeVector momentum = testNucleus.GetMomentum();
    cout << momentum[0] <<"\t"<< momentum[1] << endl; 

    for(int i=0; i < n; i++){
      if(i>= (Z - 1)) break;
      testNucleus.ApplyRecoil(trackVector);
      momentum = testNucleus.GetMomentum();
      cout << momentum[0] <<"\t"<< momentum[1] << endl; 
      if(i <= (n-1)) testNucleus.ChangeParameters(-2, -1); //reduce nucleus mass

    }
  }

  if(choice=='t'){
    
    G4ThreeVector momentum;
    G4double time = 0.0001; //time step

    G4double x=0, y=0;
    cout << x << "\t" << y << endl;

    G4double xVelocity = pInit[0]/testNucleus.GetAtomicMass()*3E+08;
    x += xVelocity*time; 

    G4double yVelocity = pInit[1]/testNucleus.GetAtomicMass()*3E+08;
    y += yVelocity*time;     

    cout << x << "\t" << y << endl;
    
    for(int i=0; i < n; i++){
      if(i>= (Z - 1)) break;
      testNucleus.ApplyRecoil(trackVector);
      momentum = testNucleus.GetMomentum();

      xVelocity = momentum[0]/testNucleus.GetAtomicMass()*3E+08;
      x += xVelocity*time; 

      yVelocity = momentum[1]/testNucleus.GetAtomicMass()*3E+08;
      y += yVelocity*time;     

      cout << x << "\t" << y << endl;
      
      if(i <= (n-1)) testNucleus.ChangeParameters(-2, -1); //reduce nucleus mass

    }

  }
    
  return 0;
}
