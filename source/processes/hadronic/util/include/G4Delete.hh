#ifndef Delete_h
#define Delete_h

template <class T> struct Delete{void operator()(T * aT){delete aT;}};

#endif
