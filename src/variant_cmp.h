#include "variant.h"



class VariantCmp 

{

public:

	static size_t const bucket_size = 8;

	size_t operator()(const Variant* a) const;

    bool operator()(const Variant* a, const Variant* b) const;

};

