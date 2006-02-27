// $Id: Track.hh,v 1.1 2006-02-27 10:05:24 kmura Exp $
// ====================================================================
//   Track.hh
//
//                                         2005 Q
// ====================================================================
#ifndef TRACK_H
#define TRACK_H


// ====================================================================
//
// class definition
//
// ====================================================================
class Step;

class Track {
private:
  Step* step;

public:
  Track();
  ~Track();

  const Step* GetStep() const;
  const Step* GetStep1() const;
  const Step* GetStep2() const;
  const Step* GetStep3() const;
  const Step* GetStep4() const;

};

#endif
