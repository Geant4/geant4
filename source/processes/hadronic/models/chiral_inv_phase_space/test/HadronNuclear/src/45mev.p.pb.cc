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
      for(i=0; i<15; i++)
      {
        aBin.first += max/15.;
	bins.push_back(aBin);
      }
      aBin.second = 0;
      aBin.first = 0;
      for(i=0; i<15; i++)
      {
        aBin.first += max/15.;
	th1.push_back(aBin);
	th2.push_back(aBin);
	th3.push_back(aBin);
	th4.push_back(aBin);
	th5.push_back(aBin);
	th6.push_back(aBin);
	th7.push_back(aBin);
	th8.push_back(aBin);
	th9.push_back(aBin);
	th10.push_back(aBin);
	th11.push_back(aBin);
	th12.push_back(aBin);
      }
    }
    
    operator()(MEvent * it) 
    {
      double weight = 1803./10000./3.; // sgima/10000./binwidth
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
	  if( fabs(it->getTheta()-10.)<accept  ) th2[j].second+=weight;
	  if( fabs(it->getTheta()-20.)<accept  ) th3[j].second+=weight;
	  if( fabs(it->getTheta()-30.)<accept )  th4[j].second+=weight;
	  if( fabs(it->getTheta()-40.)<accept )  th5[j].second+=weight;
	  if( fabs(it->getTheta()-60.)<accept )  th6[j].second+=weight;
	  if( fabs(it->getTheta()-80.)<accept )  th7[j].second+=weight;
	  if( fabs(it->getTheta()-100.)<accept ) th8[j].second+=weight;
	  if( fabs(it->getTheta()-120.)<accept ) th9[j].second+=weight;
	  if( fabs(it->getTheta()-140.)<accept ) th10[j].second+=weight;
	  if( fabs(it->getTheta()-150.)<accept ) th11[j].second+=weight;
	  if( fabs(it->getTheta()-160.)<accept ) th12[j].second+=weight;
	  break;
	}
      }
    }
    
    struct print{ operator()(pair<double, double>& it){cout << it.first<<" "<<it.second<<endl;}};
    void dump()
    {
      cout << "energy"<<endl; for_each(bins.begin(),  bins.end(),  print());
      cout << "th=0"<<endl;   for_each(th1.begin(),   th1.end(),   print());
      cout << "th=10"<<endl;  for_each(th2.begin(),   th2.end(),  print());
      cout << "th=20"<<endl;  for_each(th3.begin(),   th3.end(),  print());
      cout << "th=30"<<endl;  for_each(th4.begin(),   th4.end(), print());
      cout << "th=40"<<endl;  for_each(th5.begin(),   th5.end(), print());
      cout << "th=60"<<endl;  for_each(th6.begin(),   th6.end(), print());
      cout << "th=80"<<endl;  for_each(th7.begin(),   th7.end(), print());
      cout << "th=100"<<endl; for_each(th8.begin(),   th8.end(), print());
      cout << "th=120"<<endl; for_each(th9.begin(),   th9.end(), print());
      cout << "th=140"<<endl; for_each(th10.begin(),  th10.end(), print());
      cout << "th=150"<<endl; for_each(th11.begin(),  th11.end(), print());
      cout << "th=160"<<endl; for_each(th12.begin(),  th12.end(), print());
    }
  private:
  
    vector<pair<double, double> > bins;
    
    vector<pair<double, double> > th1;
    vector<pair<double, double> > th2;
    vector<pair<double, double> > th3;
    vector<pair<double, double> > th4;
    vector<pair<double, double> > th5;
    vector<pair<double, double> > th6;
    vector<pair<double, double> > th7;
    vector<pair<double, double> > th8;
    vector<pair<double, double> > th9;
    vector<pair<double, double> > th10;
    vector<pair<double, double> > th11;
    vector<pair<double, double> > th12;
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
  analyse theA = for_each(data.begin(), data.end(), analyse(45.));
  theA.dump();
}
