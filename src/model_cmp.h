#include "model.h"



struct ModelCmp

{

public:

	static size_t const bucket_size = 8;

	size_t operator()(const Model* a) const;

    bool operator()(const Model* a, const Model* b) const;

};

