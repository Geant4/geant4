#ifndef G4Pair_h
#define G4Pair_h

template <class t1, class t2>
class G4Pair
{
  typedef t1 first;
  typedef t2 rest;
};

class G4Terminator
{
};

#define GROUP1(t1) G4Pair<t1, G4Terminator>
#define GROUP2(t1,t2) G4Pair<t1, GROUP1(t2)>
#define GROUP3(t1,t2,t3) G4Pair<t1, GROUP2(t2,t3)>
#define GROUP4(t1,t2,t3,t4) G4Pair<t1, GROUP3(t2,t3,t4)>
#define GROUP5(t1,t2,t3,t4,t5) G4Pair<t1, GROUP4(t2,t3,t4,t5)>
#define GROUP6(t1,t2,t3,t4,t5,t6) G4Pair<t1, GROUP5(t2,t3,t4,t5,t6)>
#define GROUP7(t1,t2,t3,t4,t5,t6,t7) G4Pair<t1, GROUP6(t2,t3,t4,t5,t6,t7)>
#define GROUP8(t1,t2,t3,t4,t5,t6,t7,t8) G4Pair<t1, GROUP7(t2,t3,t4,t5,t6,t7,t8)>
#define GROUP9(t1,t2,t3,t4,t5,t6,t7,t8,t9) G4Pair<t1, GROUP8(t2,t3,t4,t5,t6,t7,t8,t9)>
#define GROUP10(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10) G4Pair<t1, GROUP9(t2,t3,t4,t5,t6,t7,t8,t9,t10)>
#define GROUP11(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11) G4Pair<t1, GROUP10(t2,t3,t4,t5,t6,t7,t8,t9,t10,t11)>
#define GROUP12(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12) G4Pair<t1, GROUP11(t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12)>
#define GROUP13(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13) G4Pair<t1, GROUP12(t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13)>
#define GROUP14(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14) G4Pair<t1, GROUP13(t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14)>
#define GROUP15(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15) G4Pair<t1, GROUP14(t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15)>

template <class group, class functor>
class G4ForEach
{
  private:
  template <class T1> static void ActOn(T1* aT1)
  {
    typename group::first aL;
    functor aF;
    aF(&aL, aT1);
  }
  public:
  template <class T1> static void Apply(T1* aT1=0)
  {
    ActOn(aT1);
    G4ForEach<typename group::rest, functor>::Apply(aT1);
  }
  static void Apply()
  {
    ActOn<int>(0);
    G4ForEach<typename group::rest, functor>::Apply();
  }
};

template <class functor>
struct G4ForEach<G4Terminator, functor>
{
  template <class T1> static void Apply(T1* ){};
  static void Apply(){};
};

template <class t1, int i2, class t3>
class G4INT
{
  typedef t1 it;
  enum {i = i2};
  typedef t3 rest;
};

#define INT1(t1,i2) G4INT<t1, i2, G4Terminator>
#define INT2(t1,i2,i3) G4INT<t1, i3, INT1(t1,i2)>
#define INT3(t1,i2,i3,i4) G4INT<t1, i4, INT2(t1,i2,i3)>
#define INT4(t1,i2,i3,i4,i5) G4INT<t1, i5, INT3(t1,i2,i3,i4)>
#define INT5(t1,i2,i3,i4,i5,i6) G4INT<t1, i6, INT4(t1,i2,i3,i4,i5)>
#define INT6(t1,i2,i3,i4,i5,i6,i7) G4INT<t1, i7, INT5(t1,i2,i3,i4,i5,i6)>
#define INT7(t1,i2,i3,i4,i5,i6,i7,i8) G4INT<t1, i8, INT6(t1,i2,i3,i4,i5,i6,i7)>
#define INT8(t1,i2,i3,i4,i5,i6,i7,i8,i9) G4INT<t1, i9, INT7(t1,i2,i3,i4,i5,i6,i7,i8)>
#define INT9(t1,i2,i3,i4,i5,i6,i7,i8,i9,i10) G4INT<t1 i10, INT8(t1,i2,i3,i4,i5,i6,i7,i8,i9)>
#define INT10(t1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11) G4INT<t1, i11, INT9(t1,i2,i3,i4,i5,i6,i7,i8,i9,i10)>
#define INT11(t1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12) G4INT<t1, i12, INT10(t1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11)>
#define INT12(t1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13) G4INT<t1, i13, INT11(t1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12)>
#define INT13(t1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14) G4INT<t1, i14, INT12(t1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13)>
#define INT14(t1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15) G4INT<t1, i15, INT13(t1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14)>
#define INT15(t1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16) G4INT<t1, i16, INT14(t1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15)>
#endif

