#ifndef DumpFrame_h
#define DumpFrame_h

class DumpFrames
{
public:
static void DumpFrame(G4FastVector<G4ReactionProduct,128> &vec, G4int vecLen)
{
//  cout << vecLen<<endl;
//  for(G4int i=0; i<vecLen; i++)
//  {
//    cout << vec[i]->GetDefinition()->GetPDGEncoding()<<" ";
//    cout << vec[i]->GetPositionInNucleus()<<" ";
//    G4double x,y,z;
//    do
//    {
//      x = G4UniformRand()*10-5;
//      y = G4UniformRand()*10-5;
//      z = G4UniformRand()*10-5;
//    }
//    while(sqrt(x*x+y*y+z*z)>5);
//    cout << x<<" "<<y<<" "<<z;
//    cout << vec[i]->GetMomentum()<<" ";
//    cout << vec[i]->GetTotalEnergy();
//    cout << endl;
//  }
}
};

#endif
