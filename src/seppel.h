#include "casper.h"



class Seppel

{

public:

	Seppel(DataFrame* frame);

	Seppel(DataFrame* frame, double* nvarPrior, double* nexonPrior);

	~Seppel();



	double calcIntegral(Model* model);  //compute integral (requires computing mode)

	double calcIntegral(Model* model, Model* similarModel); //same but initializes mode using mode from a similar model (saves time)

	double* initMode(Model* model, Model* similarModel);



	void exploreExact(); // exhaustive enumeration of all possible models

	void exploreUnif(int runs); //Metropolis-Hastings MCMC with independent proposals (uniform)

	void exploreSmart(Model* startmodel, int runs); // Metropolis-Hastings MCMC with random walk (uses SeppelSmartDist as proposal)



	map<Model*, double*, ModelCmp> resultModes();

	map<Model*, double, ModelCmp> resultPPIntegral();  //compute post prob using marginal likelihoods. Sets integralSum & integralMax

	map<Model*, double, ModelCmp> resultPPMCMC();



	static void normalizeIntegrals(double *probs, double *values, int n);

	double integralSum; //sum integrals/exp(integralMax)

	double integralMax; //maximum log(integrals), i.e. log(marginal likelihood) + log(prior)



	int modelUnifPrior; //set to 1 to assign uniform prior on model space

	double calculatePrior(Model* model); //compute log-prior probability



private:

	DataFrame* frame;



	vector<Variant*> *varis;  //stores all considered variants linearly (good for explicit enumeration)

	set<Variant*> *varisSet; //same but stores in a set (good for random search where 1 variant/model can be visited twice)



	vector<Model*> *models; //stores all considered models linearly (good for explicit enumeration)

	set<Model*> *modelsSet; //same but stores in a set (good for random search where 1 model can be visited twice)



	vector<double> priorpNbExons; //prior prob for a single variant to have 1,2... exons

	vector<double> priorpNbVars;  //log prior prob for 1,2... variants

	vector<int> nvarsPoibin; //partial results needed to compute poibinProbs

	vector<double> Tvector;  //partial results needed to compute poibinProbs

	vector<double> poibinProbs; //Poisson-Binomial probs, needed by calculatePrior



	map<Model*, double, ModelCmp> counts;

	map<Model*, double, ModelCmp> integrals; //stores log(marginal likelihood) + log(model prior prob)

        map<Model*, double, ModelCmp> priorprobs; //stores log(model prior prob)

	map<Model*, double*, ModelCmp> modes;

};



// moved header of SmartModelDist to here because Seppel and SmartModelDist now use eachother

class SmartModelDist

{

public:

  SmartModelDist(Seppel* seppel, DataFrame* frame, Model* center, double exp_exons, set<Model*> *models);

  ~SmartModelDist();

	

  // sample a proposal

  Model* sample(set<Variant*> *varisSet);

  // densitiy of the proposal

  double densityLn(Model* model);



private:

  // extra weight of an exon if its already used

  static const double exon_weight;

  // probability to create, delete is 1-pcreate

  static const double create_prob;

	

  DataFrame* frame;

  Seppel* seppel;

  Model* center;



    // expected number of exons for a gene of this length

    double exp_exons;



    // number of times each exon of the gene was used in a variant

    int* exon_used;



    // probability of an exon with a given used count to be in a variant

    double* exon_prob;



    double pnull;

    double pcreate;

	

	map<Model*, double, ModelCmp> removeprobs;



	void updatepks();

	void buildrmtable(set<Model*> *models);

	//void buildrmtable(vector<Model*> *models);

	Variant* makevar();

	double prob(Variant* v);

};

