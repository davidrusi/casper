#include "discretedf.h"

#include "cppmemory.h"

#include <stdio.h>



DiscreteDF::DiscreteDF(double* data, int* vals, int size)

{

  //this->size = size;

        this->size = vals[size-1]-vals[0]+1;

        this->values= new int[this->size];

	this->prob = new double[this->size];

	this->cumu = new double[this->size];

	double sum = 0;

	int pos=0;

	for (int i = 0; i < this->size; i++)

	{

	  if(vals[pos] == i+vals[0]){

	    values[i] = vals[pos];

	    sum = sum + data[pos];

	    prob[i] = data[pos];

	    cumu[i] = sum;

	    pos++;

	  } else {

	    values[i] = i+vals[0];

	    prob[i] = 0;

	    cumu[i] = sum;

	  }

	  

	}

}



DiscreteDF::~DiscreteDF() {

  zaparray(values);   //delete [] values;

  zaparray(prob);   //delete [] prob;

  zaparray(cumu);   //delete [] cumu;

}



int DiscreteDF::value(int i) { return values[i]; }

double DiscreteDF::probability(int i) { return prob[i]; }

double DiscreteDF::cumulativeProbability(int i) { return cumu[i]; }

