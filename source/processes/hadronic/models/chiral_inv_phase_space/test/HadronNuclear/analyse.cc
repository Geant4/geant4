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
      for(i=0; i<100; i++)
      {
        aBin.first += max/100.;
	bins.push_back(aBin);
      }
      aBin.second = 0;
      aBin.first = 0;
      for(i=0; i<25; i++)
      {
        aBin.first += max/25.;
	th7.push_back(aBin);
	th30.push_back(aBin);
	th60.push_back(aBin);
	th120.push_back(aBin);
	th150.push_back(aBin);
      }
    }
    
    operator()(MEvent * it) 
    {
      int ip = it->getIP();
      if(ip!=2112) return;
      double energy = it->getEkin();
      for(int j=0; j<bins.size(); j++)
      {
	if(energy<bins[j].first)
	{
	  bins[j].second++;
	  break;
	}
      }
      for(int j=0; j<th7.size(); j++)
      {
	if(energy<th7[j].first)
	{
	  if( fabs(it->getTheta()-7.5)<10. ) th7[j].second++;
	  if( fabs(it->getTheta()-30)<10.  ) th30[j].second++;
	  if( fabs(it->getTheta()-60)<10.  ) th60[j].second++;
	  if( fabs(it->getTheta()-120)<10. ) th120[j].second++;
	  if( fabs(it->getTheta()-150)<10. ) th150[j].second++;
	  break;
	}
      }
    }
    
    struct print{ operator()(pair<double, int>& it){cout << it.first<<" "<<it.second<<endl;}};
    void dump()
    {
      cout << "energy"<<endl;
      for_each(bins.begin(),  bins.end(),  print());
      cout << "th=7.5"<<endl;
      for_each(th7.begin(),   th7.end(),   print());
      cout << "th=30"<<endl;
      for_each(th30.begin(),  th30.end(),  print());
      cout << "th=60"<<endl;
      for_each(th60.begin(),  th60.end(),  print());
      cout << "th=120"<<endl;
      for_each(th120.begin(), th120.end(), print());
      cout << "th=150"<<endl;
      for_each(th150.begin(), th150.end(), print());
    }
  private:
  
    vector<pair<double, int> > bins;
    
    vector<pair<double, int> > th7;
    vector<pair<double, int> > th30;
    vector<pair<double, int> > th60;
    vector<pair<double, int> > th120;
    vector<pair<double, int> > th150;
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
  analyse theA = for_each(data.begin(), data.end(), analyse(800.));
  theA.dump();
}
