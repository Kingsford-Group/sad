#ifndef __MASSERT_H__
#define __MASSERT_H__

#include <iostream>

#ifdef NDEBUG
#define massert(condition, message)
#else
#define massert(condition, message)                              \
  (!(condition)) ?                                              \
  (std::cerr << "Assertion failed: (" << #condition << "), "    \
   << "function " << __FUNCTION__                               \
   << ", file " << __FILE__                                     \
   << ", line " << __LINE__ << "."                              \
   << '\n' << message << std::endl, abort(), 0) : 1
#endif


#endif /* __MASSERT_H__ */
