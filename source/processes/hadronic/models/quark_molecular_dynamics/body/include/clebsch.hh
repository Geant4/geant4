class ClebschGordan
{
  double logfak[101];
  double w3j(double,double,double,double,double,double);
public:
  ClebschGordan();
  double operator()(double,double,double,double,double,double);
  static ClebschGordan Clebsch;
};

