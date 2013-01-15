/**************************************************************/
/* C++ memory allocation                                      */
/**************************************************************/

#ifndef _CPPMEM_TEMPLATE_H
#define _CPPMEM_TEMPLATE_H

#include <cstddef>
#include <assert.h>
//#include <stdio.h> //debug

#if !defined(NDEBUG)
#define NDEBUG  //disables the assert calls below. Comment this line for debugging double deletion issues
#endif

// Use zap instead of delete (avoids issues with double deletion)
template <class T>
inline void zap(T & x) {
  {assert(x != NULL);}
  delete x;
  x = NULL;
}

// Use zaparray instead of delete [] (avoid issues with double deletion)
template <class T>
inline void zaparray(T & x) {
  {assert(x != NULL);}
  delete [] x;
  x = NULL;
}

//Free dynamically allocated objects in a container
//Example
// vector <Chromosome *> vc;
// list <Chromosome *> lc;
// (populate & use)
// std::for_each(vc.begin(), vc.end(), DeleteFromVector());
// vc.clear();
// std::for_each(lc.begin(), lc.end(), DeleteFromVector());
// lc.clear();

struct DeleteFromVector {
  template <class T> void operator() ( T* ptr) const {
    delete ptr;
  }
};


// Free any container of dynamically allocated objects
// Example:
//   vector <Chromosome *> vc;
//   list <Chromosome *> lc;
//   (populate & use)
//   FreeClear( &lc );
//   FreeClear( &vc );

//template <class C> void FreeClear( C * cntr ) {
//  for ( typename C::iterator it = cntr.begin(); it != cntr.end(); ++it ) {
//    delete *it;
//  }
//  cntr.clear();
//}


#endif
