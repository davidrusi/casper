#include "fragment.h"

#include "exon.h"

#include <map>

#include <string>

#include <vector>

#include <sstream>



using namespace std;



class Variant

{

public:



	int id;

	std::string name;



	// ordered array of exons of this variant

	Exon** exons;

	int exonCount;

	// positions of the exons

	int* positions;

	// bp length of this variant

	int length;

	// forward or backward strand

	//bool strand;



	// used for hashing

	/*int codelen;

	int* codes;*/

	int hashcode;



	Variant(vector<Exon*>* exons);

	~Variant ();



	// index of the exon in the exon list of this variant

	int indexOf(int exonid);

	// checks whether this variant could explain the fragment

	bool contains(Fragment* frag);

	// checks whether the exon is used in this variant

	bool contains(Exon* v);

	

	void toString(char *str);



	// compares two variants, 0 if equal. -1 and +1 used for sorting in maps

	int compare(const Variant* other);

	// hash for this variant

	int gethash();

	

private:

	// mapping of exonid to position in variant

	map<int, int> idmap;

};

