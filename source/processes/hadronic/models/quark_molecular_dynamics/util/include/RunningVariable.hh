#ifndef __RUNNING__
#define __RUNNING__

#pragma interface

#define RunningVariable(x) for ( x.rewind(); x.more(); x.advance() )

template<class t>
class LoopVariable
{
public:
  LoopVariable(t x) : minValue(x),maxValue(x),dValue(x),cur(x),isRange(0) {}
  LoopVariable(t x1,t x2,t dx) 
    : minValue(x1),maxValue(x2),dValue(dx),cur(x1),isRange(1) {}
  
  operator t() { return cur; }
  t current() { return cur; }
  t rewind() { leftRange=0; return cur = minValue; }
  int more() { return cur <= maxValue && !leftRange; }
  t advance() { if (!isRange) ++leftRange; return cur += dValue; }
  t& operator=(const t& y) { return cur = y; }
  t& operator=(const LoopVariable<t>& y) { return cur = y.cur; }
  G4std::ostream& Print(G4std::ostream& o) { return o << minValue << " .. " << maxValue << ", " << dValue<< " "; }
private:
  int leftRange,isRange;
  t minValue;
  t dValue,maxValue,cur;
};

#endif
