#include "model_cmp.h"

#include "discretedf.h"

#include <vector>

using namespace std;



class DataFrame

{

public:

	DataFrame(DiscreteDF* fraglen_dist, double (*fragsta_cumu)(double));

	~DataFrame();



	// all exons mapped by their id

	vector<Exon*> exons;

	map<int, Exon*> id2exon;

	// all fragments

	list<Fragment*> data;



	void addData(Fragment* f);

	void addExon(Exon* e);



	// probabilities of all the fragments given a variant

	map<Fragment*, double> probabilities(Variant* v);  //uses cache if available, otherwise fills cache and returns prob

	double probability(Variant* v, Fragment* f);  //does not use cache



	// create variant from fragment

	Variant* path2Variant(Fragment* f); 

	int fixUnexplFrags(set<Variant*, VariantCmp>* initvars, int denovo);



	// returns a list of all possible models that could explain this data

	void allModels(vector<Variant*> *varis, vector<Model*> *models);

        void allVariants(vector<Variant*> *varis);

	int frag_readlen;



	//	void debugprint();

private:



	int fraglen_minx;

        int fraglen_maxx;

	DiscreteDF* fraglen_dist;

	double (*fragsta_cumu)(double x);

        map<Variant*, map<Fragment*, double>, VariantCmp > cache;



	double prob(int fs, int fe, int bs, int be, int* pos, double T);

	void allVariantsRec(vector<Exon*>* stack, unsigned int level, vector<Variant*>* vars);

	void allModelsRec(vector<Variant*>* stack, unsigned int level, vector<Variant*>* vars, vector<Model*>* models);

};

