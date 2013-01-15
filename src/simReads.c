#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "header.h"
#include "hash.h"
#include "simReadsfunc.h"


SEXP casperSimC(SEXP gene_exp, SEXP var_exp, SEXP var_num, SEXP var_len, SEXP exon_num, SEXP exon_st, SEXP exon_end, SEXP exon_id, SEXP len_distrV, SEXP len_distrD, SEXP st_distrV, SEXP st_distrD, SEXP read_len, SEXP nn, SEXP gene_strand, SEXP lr_fileR, SEXP rr_fileR, SEXP chr, SEXP rseed, SEXP rbam, SEXP rinsideBam){

  int i=0, *ge, *vn, *vl, *en, *es, *ee, *ei, *gs, ngenes, *ldv, rl, ldlen, sdlen, n, bam, insideBam;
  double *ve, *ldd, *sdv, *sdd;
  FILE *LRFILE=NULL, *RRFILE=NULL;
  SEXP startsTmp;

  PROTECT(gene_exp = coerceVector(gene_exp, INTSXP));
  PROTECT(var_exp = coerceVector(var_exp, REALSXP));
  PROTECT(var_num = coerceVector(var_num, INTSXP));
  PROTECT(var_len = coerceVector(var_len, INTSXP));
  PROTECT(exon_num = coerceVector(exon_num, INTSXP));
  PROTECT(exon_st = coerceVector(exon_st, INTSXP));
  PROTECT(exon_end = coerceVector(exon_end, INTSXP));
  PROTECT(exon_id = coerceVector(exon_id, INTSXP));
  PROTECT(gene_strand = coerceVector(gene_strand, INTSXP));
  PROTECT(len_distrV = coerceVector(len_distrV, INTSXP));
  PROTECT(len_distrD = coerceVector(len_distrD, REALSXP));
  PROTECT(st_distrV = coerceVector(st_distrV, REALSXP));
  PROTECT(st_distrD = coerceVector(st_distrD, REALSXP));
  PROTECT(read_len = coerceVector(read_len, INTSXP));
  PROTECT(nn = coerceVector(nn, INTSXP));
  PROTECT(lr_fileR = AS_CHARACTER(lr_fileR));
  PROTECT(rr_fileR = AS_CHARACTER(rr_fileR));
  PROTECT(chr = coerceVector(chr, STRSXP));
  PROTECT(rbam = coerceVector(rbam, INTSXP));
  PROTECT(rinsideBam = coerceVector(rinsideBam, INTSXP));

  ge = INTEGER(gene_exp);
  ve = REAL(var_exp);
  vn = INTEGER(var_num);
  vl = INTEGER(var_len);
  en = INTEGER(exon_num);
  es = INTEGER(exon_st);
  ee = INTEGER(exon_end);
  ei = INTEGER(exon_id);
  gs = INTEGER(gene_strand);
  ldv = INTEGER(len_distrV);
  ldd = REAL(len_distrD);
  sdv = REAL(st_distrV);
  sdd = REAL(st_distrD);
  rl = INTEGER(read_len)[0];
  n = INTEGER(nn)[0];
  bam = INTEGER(rbam)[0];
  insideBam = INTEGER(rinsideBam)[0];

  ngenes = length(var_num);
  ldlen = length(len_distrD);
  sdlen = length(st_distrD);
  gene_t *genes;
  genes = malloc((ngenes+1) * sizeof(gene_t));
  build_genes(genes, ve, vn, vl, en, es, ee, ei, gs, ngenes, chr);

  int gene, var, len, *gansS, *vansS, *lansS, *strS, seed=INTEGER(rseed)[0], *pos, st;
  double *sansS;// st;
  srand(seed);

  SEXP gans, vans, lans, sans, ans, strs, qname, rname, strand, posr, cigar;
  PROTECT(gans = allocVector(INTSXP, n));
  PROTECT(vans = allocVector(INTSXP, n));
  PROTECT(lans = allocVector(INTSXP, n));
  PROTECT(sans = allocVector(REALSXP, n));
  PROTECT(strs = allocVector(INTSXP, n));
  PROTECT(ans = allocVector(VECSXP, 12));
  if(insideBam==1){
    PROTECT(qname = allocVector(STRSXP, n*2)); 
    PROTECT(rname = allocVector(STRSXP, n*2));
    PROTECT(strand = allocVector(STRSXP, n*2));
    PROTECT(posr = allocVector(INTSXP, n*2));
    PROTECT(cigar = allocVector(STRSXP, n*2));
  } else {
    PROTECT(qname = allocVector(STRSXP, 1));
    PROTECT(rname = allocVector(STRSXP, 1));
    PROTECT(strand = allocVector(STRSXP, 1));
    PROTECT(posr = allocVector(INTSXP, 1));
    PROTECT(cigar = allocVector(STRSXP, 1));
  }
  gansS = INTEGER(gans);
  vansS = INTEGER(vans);
  lansS = INTEGER(lans);
  sansS = REAL(sans);
  strS = INTEGER(strs);
  pos = INTEGER(posr);

   if(bam==1){
    LRFILE = fopen(CHAR(STRING_ELT(lr_fileR, 0)), "a");
    if(strcmp("", CHAR(STRING_ELT(rr_fileR, 0)))!=0) {
      RRFILE = fopen(CHAR(STRING_ELT(rr_fileR, 0)), "a");
    }
  }

  int j, *starts, gap, totp=0;
  char ***cigars=NULL, seqstr[rl+1], seqnuc[2]="C", tmpchar[100];
  hash_t *paths, paths_pted;
  PROTECT(startsTmp=allocVector(INTSXP, 3));
  starts = INTEGER(startsTmp);

  if(bam==1){
    cigars = malloc((n+1) * sizeof(char **));
    for(i=0; i<n; i++) cigars[i] = malloc(3 * sizeof(char *));
    for(j=0; j<n; j++)  for(i=0; i<3; i++) cigars[j][i] = malloc(100 * sizeof(char));
    strcpy(seqstr, seqnuc);
    for (i=0; i<rl-1; i++) strcat(seqstr, seqnuc);
  }
  //starts = malloc(3 * sizeof(int));
  paths = &paths_pted;
  hash_init(paths, NextPow2(n)); 

  int l=0, cnt=0;
  char geStr[2];

  i=0;
  while(i<n) {
    j=0;
    l=0;
    cnt=0;
    //gene = choose_gene(ge, ngenes);
    gene = ge[i];
        //var = (int) ldd[i];
    var = choose_var(genes[gene]);
    //printf("Generating for Gene: %d and tx: %d. var exp: %f nreads: %d strand: %d\n", gene, var, genes[gene].vars[var].exp, n, genes[gene].strand);
    st=-1;
    //printf("Choosing len (varlen=%d)\n", genes[gene].vars[var].len);
    len = ldv[i];
    if(len>genes[gene].vars[var].len-1) len=genes[gene].vars[var].len;
    while(j==0){
      cnt++;
      st = choose_st(len, genes[gene].vars[var].len, sdv, sdd, sdlen, genes[gene].strand);
      if(st>=0) {
	if(bam==1){
	  starts=build_cigar(genes[gene].vars[var], len, st, rl, cigars[i], genes[gene].strand);
	  if(genes[gene].strand==1) vansS[i] = starts[0];
	  else vansS[i] = starts[1];
	  if(genes[gene].strand == 1) gap = starts[1]-(starts[0]+rl);
	  else gap = starts[0]-(starts[1]+rl);
	  if(genes[gene].strand==1) strcpy(geStr, "+");
	  else strcpy(geStr, "-");
	  if((strcmp(cigars[i][0], "\0")!=0)&&(strcmp(cigars[i][1], "\0")!=0)) fprintf(LRFILE, "%d.%d\t147\t%s\t%d\t64\t%s\t=\t%d\t%d\t%s\t%s\tXS:A:%s\tRG:Z:%d\n", i, var+1, genes[gene].chr, starts[0], cigars[i][0], starts[1], gap, seqstr, seqstr, geStr, var+1);
	  if(strcmp("", CHAR(STRING_ELT(rr_fileR, 0)))!=0) {
	    if((strcmp(cigars[i][0], "\0")!=0)&&(strcmp(cigars[i][1], "\0")!=0)) fprintf(LRFILE, "%d.%d\t83\t%s\t%d\t64\t%s\t=\t%d\t%d\t%s\t%s\tXS:A:%s\tRG:Z:%d\n", i, var+1, genes[gene].chr, starts[1], cigars[i][1], starts[0], gap, seqstr, seqstr, geStr, var+1);
	  }
	} 	
	build_path(genes[gene].vars[var], len, st, rl, paths, genes[gene].strand, starts);
	if(genes[gene].strand==1) vansS[i] = starts[0];
        else vansS[i] = starts[1];
	if(insideBam==1){
	  sprintf(tmpchar, "%d.%d", i, var+1); SET_STRING_ELT(qname, i*2, mkChar(tmpchar)); SET_STRING_ELT(qname, i*2+1, mkChar(tmpchar));
	  SET_STRING_ELT(rname, i*2, mkChar(genes[gene].chr)); SET_STRING_ELT(rname, i*2+1, mkChar(genes[gene].chr));
	  if(genes[gene].strand==0) { SET_STRING_ELT(strand, i*2, mkChar("-")); SET_STRING_ELT(strand, i*2+1, mkChar("-")); }
	  else { SET_STRING_ELT(strand, i*2, mkChar("+")); SET_STRING_ELT(strand, i*2+1, mkChar("+")); }
	  pos[i*2] = starts[3]; pos[i*2+1] = starts[4];
	  SET_STRING_ELT(cigar, i*2, mkChar(cigars[i][0])); SET_STRING_ELT(cigar, i*2+1, mkChar(cigars[i][1]));
	}
	totp+=starts[2];
	strS[i] = genes[gene].strand;
	gansS[i] = genes[gene].vars[var].len;
	lansS[i] = len;
	sansS[i] = st;
	i++;
	j=1;
	if(i % 1000000 ==0) Rprintf("%d fragments simulated\n", i);
      }
      if(cnt==100) {
	//printf("%f %d %d %d %d %d %f\n", st, len, genes[gene].vars[var].len, gene, var, i, genes[gene].exp); 
	break;
      }
    }
  }

  hash_node_t *bucket;
  int count=0, *p_pathc;
  SEXP key, pathc;
  PROTECT(key = allocVector(STRSXP, totp));  
  PROTECT(pathc = allocVector(INTSXP, totp));
  p_pathc = INTEGER(pathc);

  
  for(i=0; i<paths_pted.size; i++) {
    if(paths_pted.bucket[i]!=NULL) {
      bucket=paths_pted.bucket[i];
      while(bucket) {
	SET_STRING_ELT(key, count, mkChar(bucket->key));
	p_pathc[count] = bucket->data;
	bucket=bucket->next;
	count++;
      }
    }
  }
  
  
  SET_VECTOR_ELT(ans, 0, gans);
  SET_VECTOR_ELT(ans, 1, vans);
  SET_VECTOR_ELT(ans, 2, lans);
  SET_VECTOR_ELT(ans, 3, sans);
  SET_VECTOR_ELT(ans, 4, strs);
  SET_VECTOR_ELT(ans, 5, key);
  SET_VECTOR_ELT(ans, 6, pathc);
  SET_VECTOR_ELT(ans, 7, qname);
  SET_VECTOR_ELT(ans, 8, rname);
  SET_VECTOR_ELT(ans, 9, strand);
  SET_VECTOR_ELT(ans, 10, posr);
  SET_VECTOR_ELT(ans, 11, cigar);

  UNPROTECT(34);

  //free(starts);   
  if(bam==1){
    fclose(LRFILE);
    fclose(RRFILE);
  }
  return(ans);

}
