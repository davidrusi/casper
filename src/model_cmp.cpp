#include "model_cmp.h"



size_t ModelCmp::operator()(const Model* a) const

{

	Model* va = (Model*)a;

	return va->hashcode;

}

bool ModelCmp::operator()(const Model* a, const Model* b) const

{

	Model* va = (Model*)a;

	Model* vb = (Model*)b;



	return va->compare(vb) < 0;

}

