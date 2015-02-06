#include <Rcpp.h>
#include "dada.h"
// [[Rcpp::interfaces(cpp)]]

/*
 methods for "B" objects.
 The "B" object is the partition of a DNA data set into
 clusters that is updated until convergence during DADA's
 inner two loops. It does not include the current error model,
 which is updated after the covergence of "B".
 */

#define RAWBUF 50
#define FAMBUF 50
#define CLUSTBUF 50
#define TAIL_APPROX_CUTOFF 1e-7 // Should test to find optimal

/* private function declarations */
Raw *raw_new(char *seq, int reads);
Fam *fam_new();
Bi *bi_new(int totraw);
void raw_free(Raw *raw);
void fam_free(Fam *fam);
void bi_free(Bi *bi);

int fam_add_raw(Fam *fam, Raw *raw);
int bi_add_fam(Bi *bi, Fam *fam);
int b_add_bi(B *b, Bi *bi);
void bi_shove_raw(Bi *bi, Raw *raw);
void bi_add_raw(Bi *bi, Raw *raw);   // NOT WORKING, MAY BE DELETED
Raw *fam_pop_raw(Fam *fam, int r);
Raw *bi_pop_raw(Bi *bi, int f, int r);
Fam *bi_pop_fam(Bi *bi, int f);

void bi_census(Bi *bi);
void bi_consensus_update(Bi *bi, double err[4][4]);
void fam_consensus_update(Fam *fam);
void bi_fam_update(Bi *bi, double err[4][4], double score[4][4], double gap_pen);
double get_self(char *seq, double err[4][4]);
double compute_lambda(Sub *sub, double self, double t[4][4]);
Sub *al2subs(char **al);

/*
 raw_new:
 The constructor for the Raw object.
 */
Raw *raw_new(char *seq, int reads) {
  Raw *raw = (Raw *) malloc(sizeof(Raw));
  raw->seq = (char *) malloc(strlen(seq)+1);
  strcpy(raw->seq, seq);
  raw->kmer = get_kmer(seq, KMER_SIZE);
  raw->reads = reads;
  return raw;
}

void raw_free(Raw *raw) {
  free(raw->seq);
  free(raw->kmer);
  free(raw);  
}
/*
 fam_new:
 The constructor for the Fam object.
 */
Fam *fam_new() {
  Fam *fam = (Fam *) malloc(sizeof(Fam));
  fam->seq = (char *) malloc(SEQLEN);
  fam->raw = (Raw **) malloc(RAWBUF * sizeof(Raw *));
  fam->maxraw = RAWBUF;
  fam->nraw = 0;
  fam->reads = 0;
  return fam;
}

/*
 fam_free:
 The destructor for the Fam object.
 */
void fam_free(Fam *fam) {
  free(fam->seq);
  free(fam->raw);
  free(fam);
}

/*
 fam_add_raw:
 Add a raw to a fam object.
 Update reads/nraw.
 Return index to raw in fam.
 */
int fam_add_raw(Fam *fam, Raw *raw) {
//  if(VERBOSE) { printf("\tfam_add_raw - enter (%i, %i)\n", fam->nraw, fam->maxraw); }
  if(fam->nraw >= fam->maxraw) {    // Extend Raw* buffer
    fam->raw = (Raw **) realloc(fam->raw, (fam->maxraw+RAWBUF) * sizeof(Raw *));
    fam->maxraw+=RAWBUF;
  }

  fam->raw[fam->nraw] = raw;
  fam->reads += raw->reads;
  return(fam->nraw++);
}

/*
 fam_pop_raw:
 Removes a raw from a fam object.
 Update fam reads/nraw.
 Returns pointer to that raw.
 THIS REORDERS THE FAMS LIST OF RAWS. CAN BREAK INDICES.
 THIS DOES NOT UPDATE READS OR NRAW IN THE PARENT BI.
 */
Raw *fam_pop_raw(Fam *fam, int r) {
  Raw *pop;
  
  if(r<fam->nraw) {
    pop = fam->raw[r];
    fam->raw[r] = fam->raw[fam->nraw-1];
    fam->raw[fam->nraw-1] = NULL;
    fam->nraw--;
    fam->reads -= pop->reads;
  } else {  // Trying to pop an out of range raw
    printf("fam_pop_raw: Not enough raws %i (%i)\n", r, fam->nraw);
    pop = NULL;
  }
  return pop;
}

/*
 bi_pop_raw:
 Removes a raw from the child object.
 Update bi reads/nraw.
 Returns pointer to that raw.
 SEE NOTES ON FAM_POP_RAW
 */
Raw *bi_pop_raw(Bi *bi, int f, int r) {
  Raw *pop;
  
  if(f<bi->nfam) {
    pop = fam_pop_raw(bi->fam[f], r);
    if(pop != NULL) {
      bi->nraw--;
      bi->reads -= pop->reads;
    }
  } else {  // Trying to pop an out of range raw
    printf("bi_pop_raw: not enough fams %i (%i)\n", f, bi->nfam);
    pop = NULL;
  }
  return pop;
}

/*
 bi_pop_fam:
 Removes a fam from cluster.
 Update bi reads/nraw.
 Returns pointer to the removed fam.
 */
Fam *bi_pop_fam(Bi *bi, int f) {
  Fam *pop;
  
  if(f<bi->nfam) {
    pop = bi->fam[f];
    bi->fam[f] = bi->fam[bi->nfam-1];
    bi->fam[bi->nfam-1] = NULL;
    bi->nfam--;
    bi->nraw -= pop->nraw;
    for(int r=0;r<pop->nraw;r++) {  // Doesn't rely on fam reads
      bi->reads -= pop->raw[r]->reads;
    }
  } else {  // Trying to pop an out of range fam
    printf("bi_pop_fam: not enough fams %i (%i)\n", f, bi->nfam);
    pop = NULL;
  }
  return pop;
}

/*
 bi_shove_raw:
 Hackily just adds raw to fam[0].
 consensus_update and fam_update are needed on this cluster afterwards.
 */
void bi_shove_raw(Bi *bi, Raw *raw) {
  if(bi->nfam == 0) {
    bi_add_fam(bi, fam_new());
  }
  fam_add_raw(bi->fam[0], raw);
  bi->nraw++;
  bi->reads += raw->reads;
  bi->update_fam = TRUE;
}

/*
 bi_add_raw:
 Adds raw to appropriate fam. Makes new fam if necessary.
 ISSUE WITH ALIGN, NEED ERR. UNUSED FOR NOW.
 UNUSED UNUSED UNUSED UNUSED UNUSED UNUSED
 */
void bi_add_raw(Bi *bi, Raw *raw) {
//  printf("bi_add_raw\n");
  Sub *sub;
  int foo, result;
  char buf[10];  // holds sprintf("%i", fam_index)
  
  sub = bi->sub[raw->index];
  result = sm_exists(bi->sm, sub->key);
  
  if (result == 0) {                  // Handle value not found
    foo = bi_add_fam(bi, fam_new());
    fam_add_raw(bi->fam[foo], raw);
    bi->fam[foo]->sub = sub;                   // Sub set on new fam formation.
    sprintf(buf, "%i", foo);
    sm_put(bi->sm, sub->key, buf);
  } else {                            // Joining existing family
    sm_get(bi->sm, sub->key, buf, sizeof(buf));
    foo = atoi(buf);
    fam_add_raw(bi->fam[foo], raw);
  }
  
  bi->nraw++;
  bi->reads += raw->reads;
}

/* bi_new:
 The constructor for the Bi object.
 */
Bi *bi_new(int totraw) {
  Bi *bi = (Bi *) malloc(sizeof(Bi));
  bi->seq = (char *) malloc(SEQLEN);
  bi->fam = (Fam **) malloc(FAMBUF * sizeof(Fam *));
  bi->maxfam = FAMBUF;
  bi->sub = (Sub **) malloc(totraw * sizeof(Sub *));
  bi->lambda = (double *) malloc(totraw * sizeof(double));
  bi->e = (double *) malloc(totraw * sizeof(double));
  bi->sm = sm_new(HASHOCC);
  bi->update_lambda = TRUE;
  bi->update_fam = TRUE;
  
  bi->nfam = 0;
  bi->reads = 0;
  bi->nraw = 0;
  return bi;
}

/* bi_free:
 Destructs bi objects.
*/
void bi_free(Bi *bi) {
  for(int f=0;f<bi->nfam;f++) { fam_free(bi->fam[f]); }
  free(bi->fam);
  free(bi->seq);
  free(bi->sub);
  free(bi->lambda);
  free(bi->e);
  sm_delete(bi->sm);
  free(bi);
}

/* bi_add_fam:
 Adds a fam to an existing Bi. Returns its index.
 */
int bi_add_fam(Bi *bi, Fam *fam) {
  if(bi->nfam >= bi->maxfam) {    // Extend Fam* buffer
    bi->fam = (Fam **) realloc(bi->fam, (bi->maxfam+FAMBUF) * sizeof(Fam *));
    bi->maxfam+=FAMBUF;
  }

  bi->fam[bi->nfam] = fam;
  bi->reads += fam->reads;
  return(bi->nfam++);
}

void bi_census(Bi *bi) {
  int f, r;
  int reads=0, nraw=0;
  for(f=0;f<bi->nfam;f++) {
    for(r=0;r<bi->fam[f]->nraw;r++) {
      reads += bi->fam[f]->raw[r]->reads;
      nraw++;
    }
  }
  if(VERBOSE && (bi->reads != reads)) { printf("CENSUS: Census changed reads: %i -> %i\n", bi->reads, reads); }
  if(VERBOSE && (bi->nraw != nraw)) { printf("CENSUS: Census changed nraw: %i -> %i\n", bi->nraw, nraw); }
  bi->reads = reads;
  bi->nraw = nraw;
}

/* B_new:
 The constructor for the B object. Takes in a Uniques object.
 Places all sequences into the same family within the same cluster.
 New objects have flags set in such a way to force an update.
*/
B *b_new(Uniques *uniques, double err[4][4], double score[4][4], double gap_pen) {
  int i, j, index;
  char *seq;
  
  B *b = (B *) malloc(sizeof(B));
  b->bi = (Bi **) malloc(CLUSTBUF * sizeof(Bi *));
  b->maxclust = CLUSTBUF;
  b->nclust = 0;
  b->reads = 0;
  b->nraw = uniques_nseqs(uniques);
  b->gap_pen = gap_pen;
  
  // Allocate the list of raws and then create them all.
  seq = (char *) malloc(SEQLEN);
  b->raw = (Raw **) malloc(b->nraw * sizeof(Raw *));
  for (index = 0; index < b->nraw; index++) {
    uniques_sequence(uniques, index, (char *) seq);
    b->raw[index] = raw_new(seq, uniques_reads(uniques, index));
    b->raw[index]->index = index;
    b->reads += b->raw[index]->reads;
  }
  free(seq);

  // Copy the error matrix
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      b->err[i][j] = err[i][j];
    }
  }

  // Copy the score matrix
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      b->score[i][j] = score[i][j];
    }
  }
  
  // Initialize with one cluster/one-family containing all the raws.
  b_init(b);
  return b;
}

/* b_init:
 Initialized b with all raws in one cluster, with consensus, flagged to update lambda and fam.
*/
void b_init(B *b) {
  // Destruct existing clusters/fams
  for(int i=0; i<b->nclust; i++) {
    bi_free(b->bi[i]);
  }
  b->nclust=0;
  
  // Add the one cluster
  b_add_bi(b, bi_new(b->nraw));

  // Add all raws to that cluster
  for (int index=0; index<b->nraw; index++) {
    bi_shove_raw(b->bi[0], b->raw[index]);
  }

  bi_census(b->bi[0]);
  b_consensus_update(b); // Makes cluster consensus sequence
  b_lambda_update(b, FALSE, 1.);
  b_e_update(b);
  if(VERBOSE) { printf("b_init - exit\n"); }
}

/* b_free:
  Destruct the B object.
*/
void b_free(B *b) {
  for(int i=0;i<b->nclust;i++) { bi_free(b->bi[i]); }
  free(b->bi);

  for (int index = 0; index < b->nraw; index++) {
    raw_free(b->raw[index]);
  }
  free(b->raw);

  free(b);
}

/* b_add_bi:
 Adds a new cluster Bi to the clustering B. Returns its index.
 */
int b_add_bi(B *b, Bi *bi) {
  if(b->nclust >= b->maxclust) {    // Extend Bi* buffer
    b->bi = (Bi **) realloc(b->bi, (b->maxclust+CLUSTBUF) * sizeof(Bi *));
    b->maxclust+=CLUSTBUF;
  }
  b->bi[b->nclust] = bi;
  return(b->nclust++);
}

/*
 lambda_update:
 updates the alignments and lambda of all raws to Bi with
 updated consensus sequences. 
*/
void b_lambda_update(B *b, bool use_kmers, double kdist_cutoff) {
  int i, index;
  double lambda;
  char **al; // stores alignments
  Sub *sub; // stores Sub structs
  for (i = 0; i < b->nclust; i++) {
    if(b->bi[i]->update_lambda) {   // consensus sequence for Bi[i] has changed
      // update alignments and lambda of all raws to this sequence
      if(tVERBOSE) printf("C%iLU:", i);
      for(index=0; index<b->nraw; index++) {
        /* perform alignment */
        al = raw_align(b->bi[i]->center, b->raw[index], b->score, b->gap_pen, use_kmers, kdist_cutoff);
        sub = al2subs(al);
        
        /* Store sub and lambda in the cluster object Bi */
        b->bi[i]->sub[index] = sub;
        lambda = compute_lambda(sub, b->bi[i]->self, b->err);  // WHERE IS SELF SET!??  -- IN BI_CONSENSUS_UPDATE
        b->bi[i]->lambda[index] = lambda;
//        if(index == TARGET_RAW) printf("lam(TARG)=%.2e; ", b->bi[i]->lambda[index]);
      }
      b->bi[i]->update_lambda = FALSE;
    } // if(b->bi[i]->update_lambda)
  }
  b_e_update(b);
}

/* b_fam_update(B *b):
   Creates and allocates new fams based on the sub between each raw and the
   cluster consensus.
  Currently completely destructive of old fams.
   */
void bi_fam_update(Bi *bi, double err[4][4], double score[4][4], double gap_pen) {
  int foo, f, r, result, r_c;
  Sub *sub;
  char buf[10];
  
  sm_delete(bi->sm);   // May want to do better here at some point for performance
  bi->sm = sm_new(HASHOCC);
  bi_census(bi);
  
  // Make list of pointers to the raws
  Raw **raws = (Raw **) malloc(bi->nraw * sizeof(Raw *));
  r_c=0;
  for(f=0;f<bi->nfam;f++) {
    for(r=0;r<bi->fam[f]->nraw;r++) {
      raws[r_c] = bi->fam[f]->raw[r];
      r_c++;
    }
  }
  
  if(r_c != bi->nraw) {
    printf("Warning: bi_fam_update --- nraw inconsistent (%i, %i).\n", r_c, bi->nraw);
  }
  
  // Destruct old fams
  for(f=0;f<bi->nfam;f++) { fam_free(bi->fam[f]); }
  bi->nfam = 0;
  
  // Construct hash from raw->sub->key to new fam index.
  // Make fams, that contain all raws with the same substitution pattern.
  for(r_c=0;r_c<bi->nraw;r_c++) {
    // Make hash key from substitutions
    sub = bi->sub[raws[r_c]->index];
    // Place raw in fams.
    if(!sub) { // Protect from null subs, but this should never arise...
      printf("Warning: bi_fam_update hit a null sub. THIS SHOULDNT HAPPEN.\n");
      sub = al2subs(raw_align(bi->center, raws[r_c], score, gap_pen, FALSE, 1.));
    }
    result = sm_exists(bi->sm, sub->key);

    if (result == 0) {                  // Handle value not found
      foo = bi_add_fam(bi, fam_new());
      fam_add_raw(bi->fam[foo], raws[r_c]);
      bi->fam[foo]->sub = sub;                   // Sub set on new fam formation.
      sprintf(buf, "%i", foo); // strmap only takes strings as values
      sm_put(bi->sm, sub->key, buf);
    } else {                            // Joining existing family
      sm_get(bi->sm, sub->key, buf, sizeof(buf));
      foo = atoi(buf);
      fam_add_raw(bi->fam[foo], raws[r_c]);
    }
  }
  
  // consensus_update/align the fams
  for(f=0;f<bi->nfam;f++) {
    fam_consensus_update(bi->fam[f]);
//    al = raw_align(bi->center, bi->fam[f]->center, score, gap_pen, FALSE, 1.);
//    sub = al2subs(al);

    sub = bi->sub[bi->fam[f]->center->index];
    if(!sub) { // Protect from null subs, but this should never arise...
      printf("WARNING: bi_fam_update hit a null sub. THIS SHOULDNT HAPPEN (2).\n");
      sub = al2subs(raw_align(bi->center, bi->fam[f]->center, score, gap_pen, FALSE, 1.));
    }

    bi->fam[f]->sub = sub;
    bi->fam[f]->lambda = compute_lambda(sub, bi->self, err);
  }
  if(tVERBOSE) printf("(nraw=%d,nfam=%d), ", bi->nraw, bi->nfam);
  bi->update_fam = FALSE;
}

void b_fam_update(B *b) {
  for (int i=0; i<b->nclust; i++) {
    if(b->bi[i]->update_fam) {  // Consensus has changed OR??? a raw has been shoved EITHER DIRECTION breaking the fam structure
      if(tVERBOSE) printf("C%iFU:", i);
      bi_fam_update(b->bi[i], b->err, b->score, b->gap_pen);
    }
  }
  
}

void b_e_update(B *b) {
  for(int i=0; i < b->nclust; i++) {
    for(int index=0; index < b->nraw; index++) {
      b->bi[i]->e[index] = b->bi[i]->lambda[index]*b->bi[i]->reads;
//      if(index == TARGET_RAW) {
//        printf("E_%i(TARG|r=%i) = %.3e; ", i, b->bi[i]->reads, b->bi[i]->e[index]);
//      }
    }
  }
}

/* b_shuffle:
 move each sequence to the bi that produces the highest expected
 number of that sequence. The center of a Bi cannot leave.
*/
void b_shuffle(B *b) {
  int ibest, index;
  Raw *raw;
  double e, maxe;
  // Iterate over raws via clusters/fams
  for(int i=0; i<b->nclust; i++) {
    for(int f=0; f<b->bi[i]->nfam; f++) {
      // IMPORTANT TO ITERATE BACKWARDS DUE TO FAM_POP_RAW!!!!!!
      for(int r=b->bi[i]->fam[f]->nraw-1; r>=0; r--) {
        // Find cluster with best e for this raw
        maxe = 0.0; ibest=-99;
        index = b->bi[i]->fam[f]->raw[r]->index;
        for(int j=0;j<b->nclust; j++) {
          e = b->bi[j]->e[index];
          if(e > maxe) {
            maxe = e;
            ibest = j;
          }
        }
        
        // Check to see if a cluster was assigned, complain if not.
        if(ibest == -99) {  // Shouldn't see this
          printf("shuffle: Failed to assign raw %i to cluster. Defaulting to no move.\n", index);
          printf("\t (e_i=%.4e, lam_i=%.4e).\n", b->bi[i]->e[index], b->bi[i]->lambda[index]);
          ibest=i;
        }
        
        // If different, move the raw to the new bi
        if(ibest != i) {
          if(index == b->bi[i]->center->index) {  // Check if center
            printf("Warning: Shuffle blocked the center of a Bi from leaving.\n");
            if(tVERBOSE) {
              printf("Attempted: Raw %i from C%i (center=%i, self=%.2e) to C%i (%.4e (lam=%.2e,n=%i) -> %.4e (lam=%.2e,n=%i))\n", \
                  index, i, b->bi[i]->center->index, get_self(b->raw[index]->seq, b->err), ibest, \
                  b->bi[i]->e[index], b->bi[i]->lambda[index], b->bi[i]->reads, \
                  b->bi[ibest]->e[index], b->bi[ibest]->lambda[index], b->bi[ibest]->reads);
              printf("%s\n", ntstr(b->raw[index]->seq));
            }
          } else {
            raw = bi_pop_raw(b->bi[i], f, r);
            bi_shove_raw(b->bi[ibest], raw);
            b->bi[i]->update_fam = TRUE;
            b->bi[ibest]->update_fam = TRUE;  // DUPLICATIVE FLAGGING FROM SHOVE_RAW FUNCTION
            if(VERBOSE) { printf("shuffle: Raw %i from C%i to C%i (%.4e -> %.4e)\n", index, i, ibest, b->bi[i]->e[index], b->bi[ibest]->e[index]); }
          }  
        }
      } //End loop(s) over raws (r).
    } // End loop over fams (f).
  } // End loop over clusters (i).

}

/* b_p_update:
 Calculates the abundance p-value for each family in the clustering.
 Depends on the lambda between the fam and its cluster, and the reads of each.
*/
void b_p_update(B *b) {
  int i, f, reads;
  double self, mu, lambda, pval, norm;
  Fam *fam;
  for(i=0;i<b->nclust;i++) {
    for(f=0;f<b->bi[i]->nfam;f++) {
      fam = b->bi[i]->fam[f];

      // Calculate abundance pval
      if(fam->reads < 1) {
        printf("Warning: No or negative reads (%i) in fam %i.\n", fam->reads, f);
        pval=1.;
      } 
      else if(fam->reads == 1) {   // Singleton. No abundance pval.
        pval=1.;
      } 
      else if(!(fam->sub)) { // Outside kmer threshhold
        pval=0.;
      } 
      else if(fam->sub->nsubs == 0) { // Cluster center
        pval=1.;
      }
      else if(fam->lambda == 0) { // Zero expected reads of this fam
        pval = 0.;
      } else { // Calculate abundance pval.
        // mu is the expected number of reads for this fam
        mu = fam->lambda * b->bi[i]->reads;
        
        // Calculate norm (since conditioning on sequence being present).
        norm = (1.0 - exp(-mu));
        if(norm < TAIL_APPROX_CUTOFF) {
          norm = mu - 0.5*pow(mu,2.);    // Assumption: TAIL_APPROX_CUTOFF is small enough to terminate taylor expansion here
        }
        
        // Calculate pval from poisson cdf.
        if(IMPLEMENTATION == 'R') {
          Rcpp::IntegerVector n_repeats(1);
          n_repeats(0) = fam->reads-1;
          Rcpp::NumericVector res = Rcpp::ppois(n_repeats, mu, false);  // lower.tail = false
          pval = Rcpp::as<double>(res);
          //*(res.begin());
          
/*          double gslval = 1 - gsl_cdf_poisson_P(reads-1, mu);
          double minval = (pval < gslval) ? pval : gslval;
          if( fabs(pval-gslval)/minval > 1e-5 && fabs(pval-gslval) > 1e-25 ) {
            Rcpp::Rcout << "Pval disagreement (gsl/R) for mu=" << mu << " and n=" << reads-1 << ": " << gslval << ", " << pval << "\n";
          }
        } else {
          pval = 1 - gsl_cdf_poisson_P(reads-1, mu);  // THIS MUST EQUAL ZERO WHEN MU=DBL_MIN ... DOES IT?
          // THIS NEEDS TO BE CHANGED. DOES NOT LIMIT APPROPRIATELY!!!! */
        }
        pval = pval/norm;
      }
      
      // Assign (abundance) pval to fam->pval
      b->bi[i]->fam[f]->p = pval;
    }
  }
}

/* b_bud:
 Finds the minimum p-value. If significant, creates a new cluster and moves the
 raws from the fam with the minimum p-value to the new cluster.
 Returns index of new cluster, or 0 if no new cluster added.
*/

int b_bud(B *b, double omegaA) {
  int rval=0;
  int i, f, r;
  int mini=0, minf=0, totfams=0;
  double minp = 1.;
  int minreads = 0;
  Fam *fam;

  // Find i, f indices and value of minimum pval.
  for(i=0;i<b->nclust;i++) {
    for(f=0; f<b->bi[i]->nfam; f++) {
      totfams++;
      if(b->bi[i]->fam[f]->p < minp) { // Most significant
        mini = i; minf = f;
        minp = b->bi[i]->fam[f]->p;
        minreads = b->bi[i]->fam[f]->reads;
      } 
      else if((b->bi[i]->fam[f]->p == minp) && (b->bi[i]->fam[f]->reads > minreads)) {
        // Ties occur at p=0 (underflow). In that case choose the fam with most reads.
        mini = i; minf = f;
        minp = b->bi[i]->fam[f]->p;
        minreads = b->bi[i]->fam[f]->reads;
      }
    }
  }
  
  // Bonferroni correct the abundance pval by the number of fams
  // (quite conservative, although probably unimportant given the abundance model issues)
  if(minp*totfams >= omegaA) {  // Not significant, return 0
    rval = 0;
  }
  else {  // A significant abundance pval
    fam = bi_pop_fam(b->bi[mini], minf);
    i = b_add_bi(b, bi_new(b->nraw));
    
    // Move raws into new cluster, could be more elegant.
    for(r=0;r<fam->nraw;r++) {
      bi_shove_raw(b->bi[i], fam->raw[r]);
    }

    if(tVERBOSE) { 
      printf("\nNew cluster from C%iF%i: p*=%.3e\n", mini, minf, minp*totfams);
      printf("Contains raws: ");
      for(r=0;r<fam->nraw;r++) { printf("%i,", fam->raw[r]->index); }
      printf("reads = %i\n", fam->reads);
      if(fam->sub) {
        printf("Subs: %s\n", ntstr(fam->sub->key));
      } else { // Fam has a NULL sub -- Should not happen... Error?
        printf("Warning: Budded Fam had a NULL sub. THIS SHOULD NEVER BE SEEN!\n");
      }
    }
    
    bi_consensus_update(b->bi[i], b->err);
    fam_free(fam);
    rval = i;
  }
  return rval;
}

/* fam_consensus_update:
 Takes a fam object, and calculates and assigns its consensus sequence.
 Currently this is done trivially by choosing the most abundant sequence.
 NOTE: IT IS NOT CLEAR THAT THIS IS THE SAME AS THE MATLAB FUNCTIONALITY
 */
void fam_consensus_update(Fam *fam) {
  int max_reads = 0;
  int r, maxr= -99;
  
  for(r=0;r<fam->nraw;r++) {
    if(fam->raw[r]->reads > max_reads) {
      maxr=r;
      max_reads = fam->raw[r]->reads;
    }
  }
  
  if(maxr>=0) {       // Make sure there is at least one raw before strcpy
    fam->center = fam->raw[maxr];
    strcpy(fam->seq, fam->raw[maxr]->seq);
  } else {            // No raw, make empty string.
    fam->center = NULL;
    fam->seq[0] = '\0';
  }
}

/* Bi_consensus_update:
 Takes a Bi object, and calculates and assigns its consensus sequence.
 Currently this is done trivially by choosing the most abundant sequence.
 If consensus changes...
    Updates .center and .self
    Flags update_fam and update_lambda
*/
void bi_consensus_update(Bi *bi, double err[4][4]) {
  int max_reads = 0;
  int f, r;
  int maxf = 0, maxr= 0;
  
  for(f=0;f<bi->nfam;f++) {
    for(r=0;r<bi->fam[f]->nraw;r++) {
      if(bi->fam[f]->raw[r]->reads > max_reads) { // Most abundant
         maxf=f; maxr=r;
         max_reads = bi->fam[f]->raw[r]->reads;
      }
    }
  }
  
  if(strcmp(bi->seq, bi->fam[maxf]->raw[maxr]->seq) != 0) {  // strings differ
    bi->update_lambda = TRUE;
    bi->update_fam = TRUE;
    if(VERBOSE) {
      printf("bi_c_u: New Consensus from raw %i (%i)\n", bi->fam[maxf]->raw[maxr]->index, bi->fam[maxf]->raw[maxr]->reads);
      printf("\tNew(%i): %s\n", (int) strlen(bi->fam[maxf]->raw[maxr]->seq), ntstr(bi->fam[maxf]->raw[maxr]->seq));
      printf("\tOld(%i): %s\n", (int) strlen(bi->seq), ntstr(bi->seq));
    }
    bi->center = bi->fam[maxf]->raw[maxr];
    strcpy(bi->seq,bi->fam[maxf]->raw[maxr]->seq);
    bi->self = get_self(bi->seq, err);
  }
}

/* B_consensus_update:
 Updates all Bi consensi.
 */
void b_consensus_update(B *b) {
  for (int i=0; i<b->nclust; i++) {
    bi_consensus_update(b->bi[i], b->err);
  }
}

void b_update_err(B *b, double err[4][4]) {
  int nti0, nti1, i, j, f, s;
  Fam *fam;
  int32_t counts[4] = {4, 4, 4, 4};  // PSEUDOCOUNTS
  int32_t obs[4][4] = {{1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}, {1, 1, 1, 1}};
  
  // Count up all observed transitions
  for(i=0;i<b->nclust;i++) {
    // Initially add all counts to the no-error slots
    for(j=0;j<strlen(b->bi[i]->seq);j++) {
      nti0 = (int) (b->bi[i]->seq[j] - 1);
      counts[nti0] += b->bi[i]->reads;
      obs[nti0][nti0] += b->bi[i]->reads;
    }
    
    // Move counts corresponding to each substitution with the fams
    for(f=0;f<b->bi[i]->nfam;f++) {
      fam = b->bi[i]->fam[f];
      for(s=0;s<fam->sub->nsubs;s++) { // ASSUMING SUB IS NOT NULL!!
        nti0 = fam->sub->nt0[s]-1;
        nti1 = fam->sub->nt1[s]-1;
        obs[nti0][nti0] -= fam->reads;
        obs[nti0][nti1] += fam->reads;
      }
    }
  } // for(i=0;i<b->nclust;i++)
  
  // Calculate observed error rates and place in err
  for(nti0=0;nti0<4;nti0++) {
    for(nti1=0;nti1<4;nti1++) {
      err[nti0][nti1] = ((double) obs[nti0][nti1]) / ((double) counts[nti0]);
    }
  }
}

// Function to get the number of transitions observed for each possible nt combo
void b_get_trans_matrix(B *b, int32_t obs[4][4]) {
  int nti0, nti1, i, j, f, s;
  int32_t total = 0; // Will be used to check for overflows
  int32_t prev;
  Fam *fam;
  
  // Initialize obs - no pseudocounts
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      obs[i][j] = 0;
    }
  }
  
  // Count up all observed transitions
  for(i=0;i<b->nclust;i++) {
    // Initially add all counts to the no-error slots
    for(j=0;j<strlen(b->bi[i]->seq);j++) {
      nti0 = (int) (b->bi[i]->seq[j] - 1);
      obs[nti0][nti0] += b->bi[i]->reads;
      
      prev = total;
      total += b->bi[i]->reads;
      if(total < prev) { // OVERFLOW
        printf("OVERFLOW IN b_get_trans_matrix!!!!!!\n");
      }
    }
    
    // Move counts corresponding to each substitution with the fams
    for(f=0;f<b->bi[i]->nfam;f++) {
      fam = b->bi[i]->fam[f];
      for(s=0;s<fam->sub->nsubs;s++) { // ASSUMING SUB IS NOT NULL!
        nti0 = (int) (fam->sub->nt0[s]-1);
        nti1 = (int) (fam->sub->nt1[s]-1);
        obs[nti0][nti0] -= fam->reads;
        obs[nti0][nti1] += fam->reads;
      }
    }
  } // for(i=0;i<b->nclust;i++)
}


char **b_get_seqs(B *b) {
  int i;
  char **seqs = (char **) malloc(b->nclust * sizeof(char *));
  for(i=0;i<b->nclust;i++) {
    seqs[i] = (char *) malloc((strlen(b->bi[i]->seq)+1) * sizeof(char));
    ntcpy(seqs[i], b->bi[i]->seq);
  }
  return seqs;
}

int *b_get_abunds(B *b) {
  int i;
  int *abunds = (int *) malloc(b->nclust * sizeof(int));
  for(i=0;i<b->nclust;i++) {
    abunds[i] = b->bi[i]->reads;
  }
  return abunds;
}












