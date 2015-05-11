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
	
	list<Fragment*> dataM;



	void addData(Fragment* f);

	void addDataM(Fragment* f);

	void addExon(Exon* e);



	// probabilities of all the fragments given a variant

	map<Fragment*, double> probabilities(Variant* v);  //uses cache if available, otherwise fills cache and returns prob

	double probability(Variant* v, Fragment *f); //doesn't use cache. Sets checkFragSense to false
	double probability(Variant* v, Fragment* f, bool checkFragSense); //if checkFragSense==true it checks that direction of v and f align



	Variant* path2Variant(Fragment* f); 	// create single variant from fragment

	//create >1 variants from initvaris/fragment and add to newvaris. Return best in bestvar. Proposed variants already in allvarnames are not added.
	void path2Variants(set <Variant*, VariantCmp> *newvaris, set <Variant*, VariantCmp> *bestvar, set <string> *allvarnames, bool *explained, set <Variant*, VariantCmp> *initvaris, Fragment* f); 



	int fixUnexplFrags(set<Variant*, VariantCmp>* initvars, std::map<Variant*,std::string>* varshortnames, int* geneid, int denovo);



	// returns a list of all possible models

	void allModels(vector<Variant*> *varis, vector<Model*> *models, vector<Variant*> *initvaris);
	void allModels(vector<Variant*> *varis, vector<Model*> *models, set<Variant*, VariantCmp> *initvaris);

        void allVariants(vector<Variant*> *varis, set<Variant*, VariantCmp> *initvaris);
        void allVariants(vector<Variant*> *varis, vector<Variant*> *initvaris);

	int frag_readlen;

	// returns total number of counts in data and dataM objects

	int totCounts();


	//	void debugprint();

private:



	int fraglen_minx;

        int fraglen_maxx;

	DiscreteDF* fraglen_dist;

	double (*fragsta_cumu)(double x);

        map<Variant*, map<Fragment*, double>, VariantCmp > cache;



	double prob(int fs, int fe, int bs, int be, int* pos, double T);

	void allVariantsRec(vector<Exon*>* stack, unsigned int level, vector<Variant*>* vars, set<string>* inithash);

	void allModelsRec(vector<Variant*>* stack, unsigned int level, vector<Variant*>* vars, vector<Model*>* models);

};

