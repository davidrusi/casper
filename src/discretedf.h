

// Probability mass function for a discrete random variable X

// size= number of possible values for X; values= vector enumerating the possible values; prob[i]=P(X=values[i]); cumu[i]=P(X<=values[i]);  

class DiscreteDF

{

public:

        DiscreteDF(double* data, int* values, int size);

	~DiscreteDF();



	int size;

        int value(int i);

	double probability(int i);

	double cumulativeProbability(int i);



private:

        int* values;

	double* prob;

	double* cumu;

};

