#include "dropVariant.h"
#include "cstat.h"
#include <stdlib.h>
using namespace std;

dropVariant::dropVariant(int nvars) {

  this->nvars= nvars;

}


dropVariant::~dropVariant() {  

  map<string, int*>::const_iterator mi;

  for (mi = this->submodels.begin(); mi != this->submodels.end(); mi++) {

    free_ivector(mi->second, 0, nvars-1);

  }

  submodels.clear();

}


int dropVariant::size() {

  return (this->submodels).size();

}

//Add submodel to this->submodels
void dropVariant::add(int *varsin) { 

  int i;
  char *zerochar;

  //Turn varsin into string

  zerochar = (char *) calloc(nvars + 1, sizeof(char));

  for (i=0; i < this->nvars; i++) if (varsin[i]==1) zerochar[i]= '1'; else zerochar[i]= '0';
  //printf("%d %d %s\n", nvars, this->nvars, zerochar);
  std::string s (zerochar);

  //Add varsin to submodels. If already there, free memory

  if (submodels.count(s) == 0) {

    submodels[s]= varsin;

  } else {

    free_ivector(varsin, 0, nvars-1);

  }

  free((char  *) zerochar);

}

void dropVariant::erase(string s) {
//Erase submodel from this->submodels

  free_ivector(submodels.find(s)->second, 0, nvars-1);

  submodels.erase(s);

}

dropVariant* dropVariant::combinations() {
  //Create new set of submodels by combining dropped variables in this->submodels

  dropVariant* ans= new dropVariant(this->nvars);

  //Fill ans->submodels from all pairwise multiplications in this->submodels
  if (this->size() > 1) {

    int i, *varsin;
    std::map<string, int*>::const_iterator it1, it2, it1end;
    it1end= this->submodels.end();
    it1end--;

    for (it1 = this->submodels.begin(); it1 != it1end; ++it1) {

      it2= it1;
      it2++;

      while(it2 != (this->submodels.end())) {

	varsin= ivector(0,nvars-1);

	for (i=0; i<nvars; i++) varsin[i]= (it1->second[i]) * (it2->second[i]);

	ans->add(varsin);

	it2++;

      }
     
    }

  }

  return ans;

}
