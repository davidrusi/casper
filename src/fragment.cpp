#include "fragment.h"

#include "cppmemory.h"



Fragment::Fragment(int leftc, int rightc, int count)

{

	this->left = new int[leftc];

	this->right = new int[rightc];

	this->leftc = leftc;

	this->rightc = rightc;

	this->count = count;

}



Fragment::~Fragment() {

  zaparray(left); //delete [] left;

  zaparray(right); //delete [] right;

}

