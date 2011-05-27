#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>

class MEvent
{
  public:
  bool Init(ifstream &data)
  {
    if(! data) return false;
    int dummy;
    data >> dummy>> ip >> energy >> px >> py >> pz >> dummy;
    return true;
  }
  double getEkin() {return energy-939;}
  double getTheta() {return atan2(sqrt(px*px+py*py),pz)*180./3.14159265;}
  int getIP() {return ip;}
  private:
  
  int ip;
  double energy;
  double px;
  double py;
  double pz;
};

class analyse
{
  public:
    analyse(double max)
    {
      pair<double, int> aBin;
      int i;
      aBin.second = 0;
      aBin.first = 0;
      for(i=0; i<25; i++)
      {
        aBin.first += max/25.;
	bins.push_back(aBin);
      }
      aBin.second = 0;
      aBin.first = 0;
      for(i=0; i<25; i++)
      {
        aBin.first += max/25.;
	th1.push_back(aBin);
	th2.push_back(aBin);
	th3.push_back(aBin);
	th4.push_back(aBin);
	th5.push_back(aBin);
      }
    }
    
    operator()(MEvent * it) 
    {
      double weight = 1437./10000./1.; // sgima/10000./binwidth
      int ip = it->getIP();
      if(ip!=2112) return;
      double energy = it->getEkin();
      for(int j=0; j<bins.size(); j++)
      {
	if(energy<bins[j].first)
	{
	  bins[j].second+=weight;
	  break;
	}
      }
      double accept = 5.;
      for(int j=0; j<th1.size(); j++)
      {
	if(energy<th1[j].first)
	{
	  if( fabs(it->getTheta()-0.) <accept )  th1[j].second+=weight;
	  if( fabs(it->getTheta()-20.)<accept  ) th2[j].second+=weight;
	  if( fabs(it->getTheta()-60.)<accept  ) th3[j].second+=weight;
	  if( fabs(it->getTheta()-120.)<accept ) th4[j].second+=weight;
	  if( fabs(it->getTheta()-160.)<accept ) th5[j].second+=weight;
	  break;
	}
      }
    }
    
    struct print{ operator()(pair<double, double>& it){cout << it.first<<" "<<it.second<<endl;}};
    void dump()
    {
      cout << "energy"<<endl;
      for_each(bins.begin(),  bins.end(),  print());
      cout << "th=0"<<endl;
      for_each(th1.begin(),   th1.end(),   print());
      cout << "th=20"<<endl;
      for_each(th2.begin(),   th2.end(),  print());
      cout << "th=60"<<endl;
      for_each(th3.begin(),   th3.end(),  print());
      cout << "th=120"<<endl;
      for_each(th4.begin(),   th4.end(), print());
      cout << "th=160"<<endl;
      for_each(th5.begin(),   th5.end(), print());
    }
  private:
  
    vector<pair<double, double> > bins;
    
    vector<pair<double, double> > th1;
    vector<pair<double, double> > th2;
    vector<pair<double, double> > th3;
    vector<pair<double, double> > th4;
    vector<pair<double, double> > th5;
};

main()
{
  ifstream theData("current");
  vector<MEvent *> data;
  MEvent * current = new MEvent;
  while(current->Init(theData))
  {
    data.push_back(current);
    current = new MEvent;
    if(10000*(data.size()/10000) == data.size()) cerr << data.size() << endl;
  }
  cout << "Event count "<<data.size()<<endl;
  analyse theA = for_each(data.begin(), data.end(), analyse(25.));
  theA.dump();
}
