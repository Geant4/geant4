#include <fstream>
#include <vector>
#include <math.h>

main()
{
  vector<double> th;
  th.push_back(7.5);
  th.push_back(30.);
  th.push_back(60.);
  th.push_back(120.);
  
 for(int i=0; i<th.size(); i++)
 {
  cout << "Theta = "<<th[i];
  ifstream data("../data/p.pb.256mev.dsig_n.dt.domega_1");
  double e0, en, theta, sig, dsig;
  while(data)
  {
    data >> e0 >> en >> theta >> sig >> dsig;
    if(fabs(theta-th[i])<.1) 
    {
      cout << e0 <<" "<<theta<<" "<<en << " "<<sig<<" "<<dsig<<endl;
    }
  }
  cout << endl;
 }
}

