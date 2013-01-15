#include "seppel.h"

#include "cppmemory.h"



const double SmartModelDist::exon_weight = 1.0;

const double SmartModelDist::create_prob = 0.5;



SmartModelDist::SmartModelDist(Seppel* seppel, DataFrame* frame, Model* center, double exp_exons, set<Model*> *models)

{

	this->seppel = seppel;

	this->center = center;

	this->exp_exons = exp_exons;

	this->frame = frame;



	updatepks();

	buildrmtable(models);



	pnull = 0;

	vector<Variant*>::const_iterator vi;

	for (vi = center->items.begin(); vi != center->items.end(); vi++)

	{

		Variant* v = *vi;

		pnull += prob(v);

	}

	pnull += prob(0);



	if (center->items.size() == 1)

	{

		pcreate = 1;

	}

	else if (frame->exons.size() == (unsigned int)(log((center->count() + 1.0) / log(2.0)) + 0.5))

	{

		pcreate = 0;

	}

	else if (removeprobs.size() == 0)

	{

		pcreate = 1;

	}

	else

	{

		pcreate = create_prob;

	}

}



SmartModelDist::~SmartModelDist() {

  zaparray(exon_used);

  zaparray(exon_prob);

}





void SmartModelDist::updatepks()

{

	int maxexused = 0;

	int sumexused = 0;

	exon_used = new int[frame->exons.size()];

	for (unsigned int u = 0; u < frame->exons.size(); u++)

	{

		exon_used[u] = 0;

	}

	

	vector<Variant*>::const_iterator vi;



	for (vi = center->items.begin(); vi != center->items.end(); vi++)

	{

		Variant* v = *vi;



		for (int ei = 0; ei < v->exonCount; ei++)

		{

			Exon* e = frame->exons[ei];

			int i = e->num;

			exon_used[i]++;



			if (exon_used[i] > maxexused)

			{

				maxexused = exon_used[i];

			}

		}

		sumexused += v->exonCount;

	}



	double gex = frame->exons.size();

	double explen = exp_exons;

	if (explen > gex * 0.9)

	{

		explen = gex * 0.9;

	}

	double esti = (explen * maxexused - sumexused) / (gex - explen);

	double fact = max(esti + 0.001, exon_weight);



	exon_prob = new double[frame->exons.size()];

	for (unsigned int i = 0; i < frame->exons.size(); i++)

	{

		double k = exon_used[i];

		exon_prob[i] = explen * ((double)k + fact) / (sumexused + fact * gex);

	}

}



void SmartModelDist::buildrmtable(set<Model*> *models)

//void SmartModelDist::buildrmtable(vector<Model*> *models)

{

  Model** possible = new Model*[center->items.size()];  //vector of model pointers

  double* integrals = new double[center->items.size()];

  int n = 0;



	list<Variant*>* copy = new list<Variant*>(center->items.begin(), center->items.end());

	for (unsigned int i = 0; i < center->items.size(); i++)

	{

		Variant* v = copy->front();

		copy->pop_front();

		Model* m = new Model(copy);

		copy->push_back(v);



		models->insert(m);

		//models->push_back(m);

		//varisSet->insert(v);  //not needed, as these variants have already been added previously



		double like = seppel->calcIntegral(m,center); 
		//double like = seppel->calcIntegral(m); //debug

		if (like != 1)

		{

			possible[n] = m;

			integrals[n] = like;

			n++;

		}

	}



	if (n == 0) { return; }



	double *probs = new double[n];

	Seppel::normalizeIntegrals(probs, integrals, n);

	for (int i = 0; i < n; i++)

	{

		Model* m = possible[i];

		removeprobs[m] = probs[i];

	}



	zaparray(possible);

	zaparray(integrals);

	zaparray(probs);

}



Variant* SmartModelDist::makevar()

{



	vector<Exon*> *nex = new vector<Exon*>();

	for (unsigned int i = 0; i < frame->exons.size(); i++)

	{

		double pk = exon_prob[i];



		double r = runif();

		if (r < pk)

		{

			nex->push_back(frame->exons[i]);

		}

	}



	Variant* nva = new Variant(nex);



	delete nex;

	return nva;

}

double SmartModelDist::prob(Variant* v)

{

	double p = 1;

	for (unsigned int i = 0; i < frame->exons.size(); i++)

	{

		double pk = exon_prob[i];

		if (v != NULL && v->contains(frame->exons[i]))

		{

			p *= pk;

		}

		else

		{

			p *= 1.0 - pk;

		}

	}

	return p;

}



Model* SmartModelDist::sample(set<Variant*> *varisSet)

{

	vector<Variant*>* newm = new vector<Variant*>();

	vector<Variant*>::const_iterator vi;



	double x = runif();

	if (x < pcreate)

	{

		Variant* v;

		do

		{

			v = makevar();

			varisSet->insert(v);

		}

		while (v->exonCount == 0 || center->contains(v));

		

		for (vi = center->items.begin(); vi != center->items.end(); vi++)

		{

			newm->push_back(*vi);

		}

		newm->push_back(v);

		Model *ans= new Model(newm);

		delete newm;

		return ans;

	}

	else

	{

		double r = runif();

		double psum = 0;

		map<Model*, double, ModelCmp>::const_iterator mi;

		Model* last = NULL;

		for (mi = removeprobs.begin(); mi != removeprobs.end(); mi++)

		{

			psum = mi->second + psum;

			if (r <= psum)

			{

				return mi->first;

			}

			last = mi->first;

		}

		return last;



		/*int ri = runifdisc(0, varis.size() - 1);

		Variant* r = varis[ri];



		for (vi = center->items.begin(); vi != center->items.end(); vi++)

		{

			Variant* v = *vi;

			if (r->compare(v) != 0)

			{

				newm->push_back(v);

			}

		}

		return new Model(newm);*/

	}

}

double SmartModelDist::densityLn(Model* model)

{

	if (model->count() < center->count())

	{

		//return log(1.0 - pcreate) - log((double)center->count());

		return log(1.0 - pcreate) + log(removeprobs[model]);

	}

	/*else if (center->count() == model->count())

	{

		return log(pcreate) + log(

	}*/

	else

	{

		Variant* var = NULL;

		vector<Variant*>::const_iterator vi;

		for (vi = model->items.begin(); vi != model->items.end(); vi++)

		{

			Variant* v = *vi;

			if (center->contains(v) == false)

			{

				var = v;

				break;

			}

		}



		return log(pcreate) + log(prob(var)) - log(1.0 - pnull);

	}

	return DBL_MIN;

}

