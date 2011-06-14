#include "Indexing.h"
int _Indexing_base::init=0;
int _Indexing_base::offset=1;
template<> const int _Indexing<int>::flag_index=-1;
template<> const string _Indexing<string>::flag_index("*INVALID*");
