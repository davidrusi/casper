void build_genes(gene_t *genes, double *ve, int *vn, int *vl, int *en, int *es, int *ee, int *ei, int *txstr, int ngenes, SEXP chr);
int choose_gene(double *exp, int ngenes);
int choose_var(gene_t gene);
int choose_len(int *ldv, double *ldd, int ldlen);
int choose_st(int fraglen, int varlen, double *sdv, double *sdd, int sdlen, int strand);
unsigned NextPow2( unsigned x );
int *build_cigar(var_t var, int len, int st, int rl, char **cigars, int strand);
int *build_path(var_t var, int len, int st, int rl, hash_t *path, int strand, int *starts);
