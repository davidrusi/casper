#include <map>
#include <string>
#include "cstat.h"
using namespace std;


class dropVariant {

public:

  dropVariant(int nvars);

  ~dropVariant();  //frees memory in submodels->second

  //Slots

  int nvars;  //number of variants

  map <string, int*> submodels;  //Each elem indicates a submodel, e.g. <"011",[0,1,1]> indicates variant 1 is dropped, vars 2-3 kept

  //Methods

  int size(); //returns this->submodels->size()

  void add(int *varsin); //Add varsin to this->submodels, if already there then free_ivector(varsin, 0, nvars-1)

  void erase(string s);  //Erase submodel from this->submodels (and free memory from submodels->second)

  dropVariant* combinations();  //Create new set of submodels by combining dropped variables in this->submodels

};
