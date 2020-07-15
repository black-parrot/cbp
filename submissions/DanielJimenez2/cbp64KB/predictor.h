/*
 * Multiperspective Perceptron with TAGE Predictor
 * Daniel A. Jimenez
 *
 * This file simulates the 64.25KB "multiperspective perceptron with TAGE predictor" for CBP2015.
 * 
 * About 3/4 of this code is directly copied and lightly modified from
 * code provided with the MICRO 2015 paper from Seznec, San Miguel, and
 * Albericio: The Inner Most Loop Iteration counter: a new dimension in
 * branch history. The rest is Jimenez's augmentation to the hashed perceptron
 * predictor Seznec et al. refer to as the "statistical corrector"
 *
 * The IMLI-SIC and IMLI-OH techniques implemented in Seznec's code are disabled.
 * IMLI-SIC is re-implemented slightly differently in Jimenez's code. IMLI-OH
 * showed no benefit.
 */
#ifndef _PREDICTOR_H_
#define _PREDICTOR_H_

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include "utils.h"
#include "bt9.h"
#include "bt9_reader.h"

// this class contains information that is kept between prediction and update
// the only thing it has here is the prediction.  it should be subclassed to hold
// other data e.g. indices in tables, branch PC, etc.


class branch_info {
	bool _prediction;
public:
	unsigned int address;

	branch_info (void) { }

	// set the prediction for this branch

	void prediction (bool p) {
		_prediction = p;
	}

	// return the prediction for this branch

	bool prediction (void) {
		return _prediction;
	}
};

// this class represents a branch predictor

class branch_predictor {
public:
	branch_predictor (void) { }

	// the first parameter is the branch address
	// the second parameter is true if the branch was ever taken
	// the third parameter is true if the branch was ever not taken
	// the return value is a pointer to a branch_info object that gives
	// the prediction and presumably contains information needed to
	// update the predictor

	virtual branch_info *lookup (unsigned int pc) = 0;

	// the first parameter is the branch info returned by the corresponding predict
	// the second parameter is the outcome of the branch

	virtual void update (branch_info *p, unsigned int target, bool taken, int type = 0) = 0;

	// if this returns false, then we will do lookup even for always/never taken branches

	virtual bool filter_always_never (void) { return true; }

	virtual void nonconditional_branch (unsigned int pc, unsigned int target, int type) { }

	virtual const char *name (bool) { return ""; }

	virtual ~branch_predictor (void) { }
};

class multiperspective_perceptron_withTAGEbig;

class PREDICTOR {
        multiperspective_perceptron_withTAGEbig *ex;
        branch_info *u;
public:

        PREDICTOR(void);

        bool GetPrediction (UINT64 PC);
        void UpdatePredictor (UINT64 PC, OpType opType, bool resolveDir, bool predDir, UINT64 branchTarget);
        void TrackOtherInst (UINT64 PC, OpType opType, bool branchDir, UINT64 branchTarget);
};


#include "utils.h"
#include <inttypes.h>
#include <math.h>
#include <vector>
#include <string.h>

//#define PRINTSIZE

// hash for PC

#define HASHPC ((PC ^ (PC >> 4)))

// Seznec's code is full of global variables, so we wrap it and our code 
// in a namespace to avoid collisions when including it in our infrastructure

namespace multiperspective_perceptron_withTAGEbig_namespace {

// djimenez's additions here

// feature types

enum history_type {
        ACYCLIC = 1,            // "acylic history" - put the most recent history of pc into array[pc%parameter] and use this array as a feature
        MODHIST = 2,            // "modulo history" - shift history bit in if PC%modulus==0
        RECENCY = 3,            // hash of a recency stack of PCs
        MODPATH = 4,            // like modhist but with path history
        GHISTMODPATH = 5,      // (mod path history << 1) | global history
        BLURRYPATH = 6,        // "page" history sort of
        RECENCYPOS = 7,        // position of this PC in the recency stack
	IMLI = 8,		// innermost loop iteration counter
	BIAS = 9,		// plain old branch bias (mostly for testing)
        MAXTYPE = 10
};

// this struct represents one feature that is an input to the hashed
// perceptron predictor

struct history_spec {
        // the type of the feature
        enum history_type type;

        // up to three parameters for this feature

        int p1, p2, p3;

        // a coefficient multiplied by the weight value

        double coeff;

        // the size of the array of weights for this feature (not used)

        int
                size,

        // bit width of weights for this feature (5 or 6) (not used)
                width;
};

// this is a reasonable set of features derived from the multiperspective perceptron
// predictor features

history_spec tuned_specs[] = {
{ BLURRYPATH, 5, 15, -1, 2.25, 0, 6 },
{ BLURRYPATH, 8, 10, -1, 2.25, 0, 6 },
{ RECENCYPOS, 31, -1, -1, 3.5, 0, 6 },
{ GHISTMODPATH, 3, 7, 1, 2.24, 0, 6 },
{ MODPATH, 3, 20, 3, 2.24, 0, 6 },
{ IMLI, 1, -1, -1, 2.23, 0, 6 },
{ IMLI, 4, -1, -1, 1.98, 0, 6 },
{ RECENCY, 9, 3, -1, 2.51, 0, 6},
{ ACYCLIC, 12, -1, -1, 2.0, 0, 6 },
};

// number of features

int num_tables = 0;

// the history specifications, provided as a parameter to a constructor later

history_spec *specs = NULL;

// some defines

#define MAX_TABLES	64
#define MAX_GHIST	4096
#define MAX_PATHHIST	4096
#define WEIGHT_BITS               6
#define MAX_WEIGHT	((1<<(WEIGHT_BITS-1))-1)
#define MIN_WEIGHT	(-(1<<(WEIGHT_BITS-1)))
#define MAX_LG_TABLE_SIZE       13
#define MAX_TABLE_SIZE  (1<<MAX_LG_TABLE_SIZE)
#define MAX_ACYCLIC     20
#define MAX_MOD         10
#define MAX_BLURRY      16
#define MAX_BLURRY2     16
#define MAX_ASSOC       256

unsigned int 
	imli_mask1,	// which tables should have their indices hashed with the first IMLI counter
	imli_mask4;	// which tables should have their indices hashed with the fourth IMLI counter (turns out none)

// associativity of the recency stack

int assoc = 0;

// the recency stack of lower-order bits of PCs

unsigned short int recency_stack[MAX_ASSOC];

// tables of weights

int8_t tables[MAX_TABLES][MAX_TABLE_SIZE];


int 
	table_size = 1024,	// size of table (dummy value for testing)
	block_size = 21,	// number of bits in a "block"; this is the width of an initial hash of histories
	path_length = 0,	// initial maximum known path length
	modghist_length = 0;	// initial maximum modulo history length

// "blurry" path histories

unsigned int blurrypath_histories[MAX_BLURRY][MAX_BLURRY2];

// acyclic histories

bool acyclic_histories[MAX_ACYCLIC][32];

// modulo pattern history

bool mod_histories[MAX_MOD][MAX_GHIST];

// modulo path history

unsigned short int modpath_histories[MAX_MOD][MAX_PATHHIST];

int

// many different moduli could be used for building
// modulo history. these arrays keep track of which moduli
// are actually used so only the relevant tables are updated
// for performance reasons. in a real implementation, only
// those tables would be implemented.

modhist_indices[MAX_MOD],
modpath_indices[MAX_MOD],
modpath_lengths[MAX_MOD],
modhist_lengths[MAX_MOD];

// number of modulo histories and modulo path histories

int nmodhist_histories, nmodpath_histories;

// insert an item into an array without duplication. used for
// building the _indices arrays

int insert (int *v, int *n, int x) {
	for (int i=0; i<*n; i++) if (v[i] == x) return i;
	v[(*n)] = x;
	return (*n)++;
}

bool punting = false; // for testing; set to true if not enough space for perceptron tables

// analyze the specification of the features to see what the various maximum
// history lengths are, how much space is taken by various history structures,
// and how much space is left over in the hardware budget for tables. performs
// other initializations as well.

int analyze_spec (int totalbits) {
#ifdef PRINTSIZE
	fprintf (stderr, "total size so far is %d bits\n", totalbits);
#endif

	int blurrypath_bits[MAX_BLURRY][MAX_BLURRY2];
	bool acyclic_bits[MAX_ACYCLIC][32];
	bool doing_recency = false;

	// initialize predictor state

	memset (recency_stack, 0, sizeof (recency_stack));
	memset (tables, 0, sizeof (tables));
	memset (mod_histories, 0, sizeof (mod_histories));
	memset (modpath_histories, 0, sizeof (modpath_histories));
	memset (acyclic_histories, 0, sizeof (acyclic_histories));
	memset (blurrypath_histories, 0, sizeof (blurrypath_histories));

	// initialize accounting variables

	memset (blurrypath_bits, 0, sizeof (blurrypath_bits));
	memset (acyclic_bits, 0, sizeof (acyclic_bits));
	memset (modhist_lengths, 0, sizeof (modhist_lengths));

	nmodpath_histories = 0;
	for (int i=0; i<num_tables; i++) {
		if (specs[i].type == RECENCY || specs[i].type == RECENCYPOS) {
			doing_recency = true;
			if (assoc < specs[i].p1) assoc = specs[i].p1;
		}
		if (specs[i].type == ACYCLIC) {
			for (int j=0; j<specs[i].p1+2; j++) {
				acyclic_bits[specs[i].p1][j] = true;
			}
		}

		// count blurry path bits (assuming full 32-bit addresses less shifted bits)

		if (specs[i].type == BLURRYPATH) for (int j=0; j<specs[i].p2; j++) blurrypath_bits[specs[i].p1][j] = 32 - specs[i].p1;

		// count modulo history bits

		if (specs[i].type == MODHIST || specs[i].type == GHISTMODPATH) {
			int j = insert (modhist_indices, &nmodhist_histories, specs[i].p1);
			if (modhist_lengths[j] < specs[i].p2+1) modhist_lengths[j] = specs[i].p2 + 1;
			if (specs[i].p2 >= modghist_length) modghist_length = specs[i].p2+1;
		}

		for (int i=0; i<num_tables; i++) {
			if (specs[i].type == MODPATH || specs[i].type == GHISTMODPATH) {
				int j = insert (modpath_indices, &nmodpath_histories, specs[i].p1);
				if (modpath_lengths[j] < specs[i].p2+1) modpath_lengths[j] = specs[i].p2 + 1;
				if (path_length <= specs[i].p2) path_length = specs[i].p2 + 1;
			}
		}
	}

	assert (path_length < MAX_PATHHIST);
	assert (modghist_length < MAX_GHIST);

	// acount for IMLI counters

	totalbits += 64;

	// account for blurry path bits

	for (int i=0; i<MAX_BLURRY; i++) for (int j=0; j<MAX_BLURRY2; j++) totalbits += blurrypath_bits[i][j];

	// account for acyclic histories

	for (int i=0; i<MAX_ACYCLIC; i++) for (int j=0; j<32; j++) totalbits += acyclic_bits[i][j];

	// account modhist bits

	for (int i=0; i<nmodhist_histories; i++) totalbits += modhist_lengths[i];

	// account for recency stack entries

	if (doing_recency)
		totalbits += assoc * 16; // recency stack of short ints

	// figure out how much storage we have left for tables

	int contest_bits = 65536 * 8 + 2048;
	int table_bits = contest_bits - totalbits;
#ifdef PRINTSIZE
	fprintf (stderr, "%d bits left for tables\n", table_bits);
#endif
	table_size = table_bits / (num_tables * WEIGHT_BITS);
	//if ((table_size * num_tables) < (1024 * 6)) { // I want 6 kiloentries of 6-bit entries to divide among the tables.
		// now I want to see how big or small this thing needs to be
	if (table_size < 0) {
		fprintf (stderr, "too big. punting.\n");
		punting = true;
		table_size = 1024;
		fflush (stderr);
	} else {
#ifdef PRINTSIZE
		fprintf (stderr, "table size is %d entries\n", table_size);
		fflush (stderr);
#endif
	}
	totalbits += table_size * num_tables * WEIGHT_BITS;
#ifdef PRINTSIZE
	fprintf (stderr, "ending up with %d total bits\n", totalbits);
#endif
	return totalbits;
}

// insert a (shifted) PC into the recency stack with LRU replacement

void insert_recency (unsigned short int pc) {
	int i = 0;
	for (i=0; i<assoc; i++) {
		if (recency_stack[i] == pc) break;
	}
	if (i == assoc) {
		i = assoc-1;
		recency_stack[i] = pc;
	}
	int j;
	unsigned int b = recency_stack[i];
	for (j=i; j>=1; j--) recency_stack[j] = recency_stack[j-1];
	recency_stack[0] = b;
}

// hash the recency position where we find this PC

unsigned int hash_recencypos (unsigned short int pc, int l) {
	unsigned int r = 0;
	int i;

	// search for the PC

	for (i=0; i<l; i++) {
		if (recency_stack[i] == pc) {
			r = i * table_size / l;
			break;
		}
	}

	// return last index in table on a miss

	if (i == l) r = table_size-1;
	return r;
}

// hash the items in the recency stack to a given depth, shifting by a
// certain amount each iteration, with two different ways of mixing the bits

unsigned int hash_recency (int depth, int shift, int style) {
	if (style == -1) {
		unsigned int x = 0;
		for (int i=0; i<depth; i++) {
			x <<= shift;
			x += recency_stack[i];
		}
		return x;
	} else {
		unsigned int x = 0, k = 0;
		for (int i=0; i<depth; i++) {
			x ^= (!!(recency_stack[i] & (1<<shift))) << k;
			k++;
			k %= block_size;
		}
		return x;
	}
}

// hash modulo history

unsigned int hash_modhist (int a, int b, int n) {
	unsigned int x = 0, k = 0;
	for (int i=0; i<b; i++) {
		x ^= mod_histories[a][i] << k;
		k++;
		k %= n;
	}
	return x;
}

// hash acyclic history

unsigned int hash_acyclic (int a, int shift, int style) {
	unsigned int x = 0;
	assert (style == -1);
	unsigned int k = 0;
	for (int i=0; i<a+2; i++) {
		x ^= acyclic_histories[a][i] << k;
		k++;
		k %= block_size;
	}
	return x;
}

// hash "blurry" history

unsigned int hash_blurry (int scale, int depth, int shiftdelta) {
	if (shiftdelta == -1) shiftdelta = 0;
	int sdint = shiftdelta >> 2;
	int sdfrac = shiftdelta & 3;
	unsigned int x = 0;
	int shift = 0;
	int count = 0;
	for (int i=0; i<depth; i++) {
		x += blurrypath_histories[scale][i] >> shift;
		count++;
		if (count == sdfrac) {
			shift += sdint;
			count = 0;
		}
	}
	return x;
}

// hash modulo history together with modulo path history

unsigned int hash_ghistmodpath (int a, int depth, int shift) {
	unsigned int x = 0;
	for (int i=0; i<depth; i++) {
		x <<= shift;
		x += (modpath_histories[a][i] << 1) | mod_histories[a][i];
	}
	return x;
}

// hash modulo path history

unsigned int hash_modpath (int a, int depth, int shift) {
	unsigned int x = 0;
	for (int i=0; i<depth; i++) {
		x <<= shift;
		x += modpath_histories[a][i];
	}
	return x;
}

// Jimenez's own IMLI counters, separate from Seznec's

unsigned int 
	imli_counter1 = 0, 
	imli_counter4 = 0;

// use a history specification to call the corresponding history hash function

unsigned int get_hash (history_spec *s, unsigned int pc) {
	unsigned int x = 0;
	switch (s->type) {
	case BIAS:
		break;
	case IMLI:
		if (s->p1 == 1) 
			x = imli_counter1;
		else if (s->p1 == 4)
			x = imli_counter4;
		else 
			assert (0);
		break;
	case ACYCLIC:
		x = hash_acyclic (s->p1, s->p2, s->p3);
		break;
	case MODHIST:
		x = hash_modhist (s->p1, s->p2, block_size);
		break;
	case GHISTMODPATH:
		x = hash_ghistmodpath (s->p1, s->p2, s->p3);
		break;
	case MODPATH:
		x = hash_modpath (s->p1, s->p2, s->p3);
		break;
	case RECENCY:
		x = hash_recency (s->p1, s->p2, s->p3);
		break;
	case BLURRYPATH:
		x = hash_blurry (s->p1, s->p2, s->p3);
		break;
	case RECENCYPOS:
		x = hash_recencypos (pc >> 2, s->p1);
		break;
	default: assert (0);
	}
	return x;
}

// the perceptron sum

int yout;

// hash the PC

unsigned int hash_pc (unsigned int pc) {
	return pc ^ (pc >> 2);
}

// get the index into a given table using this pc and hashed pc

int get_index (unsigned int pc, unsigned int hpc, int t) {
	unsigned long long int h;

	// get the hash for the feature

	unsigned int g = get_hash (&specs[t], pc);

	// shift it and xor it with the hashed PC

	h = g;
	h <<= 20;
	h ^= hpc;

	// maybe xor in an IMLI counter

	if ((1ull<<t) & imli_mask1) h += imli_counter1;
	if ((1ull<<t) & imli_mask4) h += imli_counter4;

	// return it modulo the table size

	return h % table_size;
}

// compute yout, a partial sum to be used by Seznec's perceptron predictor

int compute_partial_sum (unsigned int pc) {
	yout = 0;
	unsigned int hpc = hash_pc (pc);
	for (int i=0; i<num_tables; i++) {
		yout += specs[i].coeff * tables[i][get_index(pc, hpc, i)];
	}
	return yout;
}

// update the tables, called from Seznec's update to the perceptron predictor

void update_partial (unsigned int pc, bool taken) {
	unsigned int hpc = hash_pc (pc);

	// update tables

	for (int i=0; i<num_tables; i++) {
		int8_t *c = &tables[i][get_index(pc,hpc,i)];
		if (taken) {
			if (*c < MAX_WEIGHT) (*c)++;
		} else {
			if (*c > MIN_WEIGHT) (*c)--;
		}
	}
}

// update histories, called after updates to the tables have been done

void update_histories (unsigned int pc, bool taken) {
	unsigned int hpc = hash_pc (pc);

	// update recency stack

	insert_recency (pc >> 2);

	// update acyclic history

	for (int i=0; i<MAX_ACYCLIC; i++) {
		acyclic_histories[i][hpc%(i+2)] = taken;
	}

	// update modpath histories

	for (int ii=0; ii<nmodpath_histories; ii++) {
		int i = modpath_indices[ii];
		if (hpc % (i+2) == 0) {
			memmove (&modpath_histories[i][1], &modpath_histories[i][0], sizeof (unsigned short int) * (modpath_lengths[ii]-1));
			modpath_histories[i][0] = hpc;
		}
	}

	// update modulo histories

	for (int ii=0; ii<nmodhist_histories; ii++) {
		int i = modhist_indices[ii];
		if (hpc % (i+2) == 0) {
			memmove (&mod_histories[i][1], &mod_histories[i][0], modhist_lengths[ii]-1);
			mod_histories[i][0] = taken;
		}
	}

	// update blurry history

	for (int i=0; i<MAX_BLURRY; i++) {
		unsigned int z = pc >> i;
		if (blurrypath_histories[i][0] != z) {
			memmove (&blurrypath_histories[i][1], &blurrypath_histories[i][0], sizeof (unsigned int) * (MAX_BLURRY2-1));
			blurrypath_histories[i][0] = z;
		}
	}
}

// Seznec's (modified) code begins here. Jimenez has changed some constants
// and added some hooks into the perceptron code above, but declines to change
// it too much lest he introduce bugs. So, it has a lot of #ifdefs guarding code
// that is never executed.
 
/*Copyright (c) <2006>, INRIA : Institut National de Recherche en Informatique et en Automatique (French National Research Institute for Computer Science and Applied Mathematics)
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


/* author A. Seznec (2015) */
// this code is derived from the code of TAGE-SC-L by A. Seznec at CBP4
// it allows to reproduce the results in the IMLI  Micro 2015 paper on the CBP4 traces
// for more realistic design (and cleaner code), use predictorTAGE-GSC-IMLI.h

// total misprediction numbers
// TAGE-SC-L (CBP4 predictor): 2.365 MPKI
// TAGE-SC-L  + IMLI:  2.226 MPKI
// TAGE-GSC + IMLI: 2.313 MPKI
// TAGE-GSC : 2.473 MPKI
// TAGE alone: 2.563 MPKI



#define SC //use the statistical corrector: comment to get TAGE alone
#define LOCALH			// use local histories

#define LOOPPREDICTOR		//  use loop  predictor
#define IMLI			// using IMLI component for hashing into things
#undef IMLISIC			// we are doing our own IMLISIC
#undef IMLIOH			// do not use IMLI-OH; not worth the space for the CBP2016 benchmarks

//#define STRICTSIZE
//uncomment to get the 256 Kbits record predictor mentioned in the paper achieves 2.228 MPKI

//#define REALISTIC
// uncomment to get a realistic predictor around 256 Kbits , with 12 1024 entries tagged tables in the TAGE predictor, and a global history and single local history GEHL statistical corrector
// total misprediction numbers
// TAGE-SC-L : 2.435 MPKI
// TAGE-SC-L  + IMLI:  2.294 MPKI
// TAGE-GSC + IMLI: 2.370 MPKI 
// TAGE-GSC : 2.531 MPKI
// TAGE alone: 2.602 MPKI

//To get the predictor storage budget on stderr  uncomment the next line

#ifdef REALISTIC 
#define LOGG 10/* logsize of the  tagged TAGE tables*/
#define TBITS 8 /* minimum tag width*/
#define  POWER
//use geometric history length
#define MINHIST 7
#define MAXHIST 1000
//probably not the best history length, but nice
#endif

//IMLI related data declaration
long long IMLIcount;
#define MAXIMLIcount 1023
#ifdef IMLI
#ifdef IMLISIC
//IMLI-SIC related data declaration
//#define LOGINB 10 // (LOG of IMLI-SIC table size +1)
#define LOGINB 11 // (LOG of IMLI-SIC table size +1) // djimenez
#define INB 1
int Im[INB] = { 10 }; // the IMLIcounter is limited to 10 bits
int8_t IGEHLA[INB][(1 << LOGINB)];
int8_t *IGEHL[INB];
#endif
//IMLI-OH related data declaration
#ifdef IMLIOH
long long localoh; //intermediate data to recover the two bits needed from the past outer iteration
#define SHIFTFUTURE 6 // (PC<<6) +IMLIcount to index the Outer History table
#define PASTSIZE 16
int8_t PIPE[PASTSIZE]; // the PIPE vector
#define OHHISTTABLESIZE 1024 //
int8_t ohhisttable[OHHISTTABLESIZE];
#ifdef STRICTSIZE
#define LOGFNB 7 // 64 entries
#define LOGFNB 7 // 64 entries
#else
#define LOGFNB 9 //256 entries
#endif
#define FNB 1
int Fm[FNB] = { 2};
int8_t FGEHLA[FNB][(1 << LOGFNB)];
int8_t *FGEHL[FNB];
#endif
#endif



#define NNN 1			// number of entries allocated on a TAGE misprediction

// number of tagged tables
#ifndef REALISTIC
#define NHIST 15
#else
#define NHIST 12
#endif

#define HYSTSHIFT 2		// bimodal hysteresis shared by 4 entries
#define LOGB 14			// log of number of entries in bimodal predictor


#if 0
#ifndef STRICTSIZE
#define PERCWIDTH 6		//Statistical corrector maximum counter width
#else
#define PERCWIDTH 7 // appears as a reasonably efficient way to use the last available bits
#endif
#endif
#define PERCWIDTH 6
//The statistical corrector components from CBP4

//global branch GEHL
#ifdef REALISTIC
#define LOGGNB 10
#else
//#define LOGGNB 9
#define LOGGNB 10
#endif

#if 0 // djimenez
#ifdef IMLI
#ifdef STRICTSIZE
#define GNB 2
int Gm[GNB] = { 17, 14 };
#else
#ifndef REALISTIC
#define GNB 4
int Gm[GNB] = { 27,22, 17, 14 };
#else
#define GNB 2
int Gm[GNB] = {  17, 14 };
#endif
#endif
#else
#ifndef REALISTIC
#define GNB 4
int Gm[GNB] = { 27,22, 17, 14 };
#else
#define GNB 2
int Gm[GNB] = {  17, 14 };
#endif
#endif

#endif
// djimenez: let's just stick with this configuration regardless of IMLI/no IMLI
#define GNB 4
int Gm[GNB] = { 27,22, 17, 14 };

/*effective length is  -11,  we use (GHIST<<11)+IMLIcount; we force the IMLIcount zero when IMLI is not used*/


int8_t GGEHLA[GNB][(1 << LOGGNB)];
int8_t *GGEHL[GNB];

//large local history
#define  LOGLOCAL 8
#define NLOCAL (1<<LOGLOCAL)
//#define INDLOCAL (PC & (NLOCAL-1)) djimenez
#define INDLOCAL (HASHPC & (NLOCAL-1))
#ifdef REALISTIC
//only one local history
#define LOGLNB 10
#define LNB 4
int Lm[LNB] = { 16, 11, 6, 3 };
int8_t LGEHLA[LNB][(1 << LOGLNB)];
int8_t *LGEHL[LNB];
#else
//three different local histories (just completely crazy :-)

#define LOGLNB 10
#define LNB 3
int Lm[LNB] = { 11, 6, 3 };
int8_t LGEHLA[LNB][(1 << LOGLNB)];
int8_t *LGEHL[LNB];
#endif

// small local history
#define LOGSECLOCAL 4
#define NSECLOCAL (1<<LOGSECLOCAL)	//Number of second local histories
#define INDSLOCAL  (((PC ^ (PC >>5))) & (NSECLOCAL-1))
#define LOGSNB 9
#define SNB 4
int Sm[SNB] = { 16, 11, 6, 3 };
int8_t SGEHLA[SNB][(1 << LOGSNB)];
int8_t *SGEHL[SNB];


//third local history
#define LOGTNB 9
#ifdef STRICTSIZE
#define TNB 2
int Tm[TNB] = { 17, 14 };
#else
#define TNB 3
int Tm[TNB] = { 22, 17, 14 };
#endif
//effective local history size +11: we use IMLIcount + (LH) << 11
int8_t TGEHLA[TNB][(1 << LOGTNB)];
int8_t *TGEHL[TNB];
#define INDTLOCAL  (((PC ^ (PC >>3))) & (NSECLOCAL-1))	// different hash for the 3rd history


long long L_shist[NLOCAL];
long long S_slhist[NSECLOCAL];
long long T_slhist[NSECLOCAL];
long long HSTACK[16];
int pthstack;
#ifndef REALISTIC
//return-stack associated history component
#ifdef STRICTSIZE
#define LOGPNB 8
#else
#define LOGPNB 9
#endif
#define PNB 4
int Pm[PNB] = { 16, 11, 6, 3 };
int8_t PGEHLA[PNB][(1 << LOGPNB)];
int8_t *PGEHL[PNB];
#else
//in this case we don't use the call stack
#define PNB 2
#define LOGPNB 11
int Pm[PNB] = { 16, 11 };

int8_t PGEHLA[PNB][(1 << LOGPNB)];
int8_t *PGEHL[PNB];
#endif


//parameters of the loop predictor
#if 0
#define LOGL 5
#define WIDTHNBITERLOOP 10	// we predict only loops with less than 1K iterations
#define LOOPTAG 10		//tag width in the loop predictor
#endif

// djimenez

#define LOGL	6
#define WIDTHNBITERLOOP 12
#define LOOPTAG 10

//update threshold for the statistical corrector
#ifdef REALISTIC
#define LOGSIZEUP 0
#else
#define LOGSIZEUP 5
#endif
int Pupdatethreshold[(1 << LOGSIZEUP)];	//size is fixed by LOGSIZEUP
//#define INDUPD (PC & ((1 << LOGSIZEUP) - 1))
#define INDUPD (HASHPC & ((1 << LOGSIZEUP) - 1))

// The three counters used to choose between TAGE ang SC on High Conf TAGE/Low Conf SC
int8_t FirstH, SecondH, ThirdH;
#define CONFWIDTH 7		//for the counters in the choser
#define PHISTWIDTH 27		// width of the path history used in TAGE




#define UWIDTH 2		// u counter width on TAGE
#define CWIDTH 3		// predictor counter width on the TAGE tagged tables

#define HISTBUFFERLENGTH 4096	// we use a 4K entries history buffer to store the branch history


//the counter(s) to chose between longest match and alternate prediction on TAGE when weak counters
#ifdef REALISTIC
#define LOGSIZEUSEALT 0
#else
#define LOGSIZEUSEALT 8
#endif
#define SIZEUSEALT  (1<<(LOGSIZEUSEALT))
//#define INDUSEALT (PC & (SIZEUSEALT -1))
#define INDUSEALT (HASHPC & (SIZEUSEALT -1))
int8_t use_alt_on_na[SIZEUSEALT][2];



long long GHIST;


//The two BIAS tables in the SC component
#define LOGBIAS 7
int8_t Bias[(1 << (LOGBIAS + 1))];
#define INDBIAS (((PC<<1) + pred_inter) & ((1<<(LOGBIAS+1)) -1))
int8_t BiasSK[(1 << (LOGBIAS + 1))];
#define INDBIASSK ((((PC^(PC>>LOGBIAS))<<1) + pred_inter) & ((1<<(LOGBIAS+1)) -1))

bool HighConf;
int LSUM;
int8_t BIM;

// utility class for index computation
// this is the cyclic shift register for folding 
// a long global history into a smaller number of bits; see P. Michaud's PPM-like predictor at CBP-1
class folded_history
{
public:


  unsigned comp;
  int CLENGTH;
  int OLENGTH;
  int OUTPOINT;

    folded_history ()
  {
  }


  void init (int original_length, int compressed_length, int N)
  {
    comp = 0;
    OLENGTH = original_length;
    CLENGTH = compressed_length;
    OUTPOINT = OLENGTH % CLENGTH;

  }

  void update (uint8_t * h, int PT)
  {
    comp = (comp << 1) ^ h[PT & (HISTBUFFERLENGTH - 1)];
    comp ^= h[(PT + OLENGTH) & (HISTBUFFERLENGTH - 1)] << OUTPOINT;
    comp ^= (comp >> CLENGTH);
    comp = (comp) & ((1 << CLENGTH) - 1);
  }

};

#ifdef LOOPPREDICTOR
class lentry			//loop predictor entry
{
public:
  uint16_t NbIter;		//10 bits
  uint8_t confid;		// 4bits
  uint16_t CurrentIter;		// 10 bits

  uint16_t TAG;			// 10 bits
  uint8_t age;			// 4 bits
  bool dir;			// 1 bit

  //39 bits per entry    
    lentry ()
  {
    confid = 0;
    CurrentIter = 0;
    NbIter = 0;
    TAG = 0;
    age = 0;
    dir = false;



  }

};
#endif


class bentry			// TAGE bimodal table entry  
{
public:
  int8_t hyst;
  int8_t pred;


    bentry ()
  {
    pred = 0;

    hyst = 1;
  }

};
class gentry			// TAGE global table entry
{
public:
  int8_t ctr;
  uint tag;
  int8_t u;

    gentry ()
  {
    ctr = 0;
    tag = 0;
    u = 0;
  }
};




int TICK;			// for the reset of the u counter


uint8_t ghist[HISTBUFFERLENGTH];
int ptghist;
long long phist;		//path history
folded_history ch_i[NHIST + 1];	//utility for computing TAGE indices
folded_history ch_t[2][NHIST + 1];	//utility for computing TAGE tags

//For the TAGE predictor
bentry *btable;			//bimodal TAGE table
gentry *gtable[NHIST + 1];	// tagged TAGE tables
int logg[NHIST+1] = {  0, 10, 11, 11, 11, 11, 11, 12, 12, 10, 11, 11, 9, 7, 7, 8,};
int m[NHIST + 1] = { 0, 5, 12, 15, 21, 31, 43, 64, 93, 137, 200, 292, 424, 612, 877, 1241, };	// history lengths tuned by djimenez 5/3/2016
int TB[NHIST + 1] = { 0, 7, 9, 9, 9, 10, 11, 11, 12, 12, 12, 13, 14, 15, 15, 15 };	// tag width for the different tagged tables
int GI[NHIST + 1];		// indexes to the different tables are computed only once  
uint GTAG[NHIST + 1];		// tags for the different tables are computed only once  
int BI;				// index of the bimodal table
bool pred_taken;		// prediction
bool alttaken;			// alternate  TAGEprediction
bool tage_pred;			// TAGE prediction
bool LongestMatchPred;
int HitBank;			// longest matching bank
int AltBank;			// alternate matching bank
int Seed;			// for the pseudo-random number generator

bool pred_inter;

#ifdef LOOPPREDICTOR
lentry *ltable;			//loop predictor table
//variables for the loop predictor
bool predloop;			// loop predictor prediction
int LIB;
int LI;
int LHIT;			//hitting way in the loop predictor
int LTAG;			//tag on the loop predictor
bool LVALID;			// validity of the loop predictor prediction
int8_t WITHLOOP;		// counter to monitor whether or not loop prediction is beneficial

#endif

int
predictorsize ()
{
  int STORAGESIZE = 0;

  int inter = 0;

  for (int i = 1; i <= NHIST; i += 1)
    {
      STORAGESIZE += (1 << (logg[i])) * (CWIDTH + UWIDTH + TB[i]);

    }
  STORAGESIZE += 2 * (SIZEUSEALT) * 4;
  STORAGESIZE += (1 << LOGB) + (1 << (LOGB - HYSTSHIFT));
  STORAGESIZE += m[NHIST];
  STORAGESIZE += PHISTWIDTH;
  STORAGESIZE += 10;		//the TICK counter

#ifdef PRINTSIZE
  fprintf (stderr, " (TAGE %d) ", STORAGESIZE);
#endif

#ifdef LOOPPREDICTOR

  inter = (1 << LOGL) * (2 * WIDTHNBITERLOOP + LOOPTAG + 4 + 4 + 1);
#ifdef PRINTSIZE
  fprintf (stderr, " (LOOP %d) ", inter);
#endif
  STORAGESIZE += inter;

#endif
#ifdef SC
  inter =0;
  
  inter += 16;			//global histories for SC
  inter = 8 * (1 << LOGSIZEUP);	//the update threshold counters
  inter += (PERCWIDTH) * 4 * (1 << (LOGBIAS));
  inter +=
    (GNB - 2) * (1 << (LOGGNB)) * (PERCWIDTH - 1) +
    (1 << (LOGGNB - 1)) * (2 * PERCWIDTH - 1);

  inter +=
    (PNB - 2) * (1 << (LOGPNB)) * (PERCWIDTH - 1) +
    (1 << (LOGPNB - 1)) * (2 * PERCWIDTH - 1);


#ifdef LOCALH
  inter +=
    (LNB - 2) * (1 << (LOGLNB)) * (PERCWIDTH - 1) +
    (1 << (LOGLNB - 1)) * (2 * PERCWIDTH - 1);
  inter += NLOCAL * Lm[0];

#ifndef REALISTIC
  inter +=
    (SNB - 2) * (1 << (LOGSNB)) * (PERCWIDTH - 1) +
    (1 << (LOGSNB - 1)) * (2 * PERCWIDTH - 1);
  inter +=
    (TNB - 2) * (1 << (LOGTNB)) * (PERCWIDTH - 1) +
    (1 << (LOGTNB - 1)) * (2 * PERCWIDTH - 1);
  inter += 16 * 16;		// the history stack
  inter += 4;			// the history stack pointer

  inter += NSECLOCAL * Sm[0];
  inter += NSECLOCAL * (Tm[0] - 11);
/* Tm[0] is artificially increased by 11 to accomodate IMLI*/
#endif // REALISTIC
#endif // LOCALH

#ifdef IMLI
#ifdef IMLIOH
  inter += OHHISTTABLESIZE;
  inter += PASTSIZE;
/*the PIPE table*/
// in cases you add extra tables to IMLI OH, the formula is correct
  switch (FNB)
    {
    case 1:
      inter += (1 << (LOGFNB - 1)) * PERCWIDTH;
      break;
    default: inter +=
	(FNB - 2) * (1 << (LOGFNB)) * (PERCWIDTH - 1) +
	(1 << (LOGFNB - 1)) * (2 * PERCWIDTH-1);
    }
#endif // IMLIOH
#ifdef IMLISIC  // in cases you add extra tables to IMLI SIC, the formula is correct
  switch (INB)
    {
    case 1:
      inter += (1 << (LOGINB - 1)) * PERCWIDTH;
      break;

    default:
      inter +=
	(INB - 2) * (1 << (LOGINB)) * (PERCWIDTH - 1) +
	(1 << (LOGINB - 1)) * (2 * PERCWIDTH-1);
    }
#endif // IMLISIC

#endif // IMLI

  inter += 3 * CONFWIDTH;	//the 3 counters in the choser
  STORAGESIZE += inter;

#ifdef PRINTSIZE
  fprintf (stderr, " (SC %d) ", inter);
#endif
#endif // SC

  int oldsize = STORAGESIZE;
  STORAGESIZE = analyze_spec (STORAGESIZE);
#ifdef PRINTSIZE
  fprintf (stderr, " (djimenez %d) ", STORAGESIZE-oldsize);
#endif
#ifdef PRINTSIZE
  fprintf (stderr, " (TOTAL %d) (%0.2fKB)\n ", STORAGESIZE, STORAGESIZE / 8192.0);
#endif

  return (STORAGESIZE);
}

class PREDICTOR
{
public:
  PREDICTOR (void)
  {
    reinit ();
    predictorsize (); // this must be called to figure out bits left for perceptron
  }


  void reinit ()
  {
#ifdef POWER

    m[1] = MINHIST;
    m[NHIST] = MAXHIST;
    for (int i = 2; i <= NHIST; i++)
      {
	m[i] =
	  (int) (((double) MINHIST *
		  pow ((double) (MAXHIST) / (double) MINHIST,
		       (double) (i - 1) / (double) ((NHIST - 1)))) + 0.5);
      }
#endif
#ifdef REALISTIC
    for (int i = 1; i <= NHIST; i++)
    {
         TB[i]= TBITS + (i/2);
         logg[i]= LOGG;
         
    }
#endif    
#ifdef LOOPPREDICTOR
    ltable = new lentry[1 << (LOGL)];
#endif


    for (int i = 1; i <= NHIST; i++)
      {
	gtable[i] = new gentry[1 << (logg[i])];
      }

    btable = new bentry[1 << LOGB];

    for (int i = 1; i <= NHIST; i++)
      {
	ch_i[i].init (m[i], (logg[i]), i - 1);
	ch_t[0][i].init (ch_i[i].OLENGTH, TB[i], i);
	ch_t[1][i].init (ch_i[i].OLENGTH, TB[i] - 1, i + 2);
      }
#ifdef LOOPPREDICTOR
    LVALID = false;
    WITHLOOP = -1;
#endif
    Seed = 0;

    TICK = 0;
    phist = 0;
    Seed = 0;

    for (int i = 0; i < HISTBUFFERLENGTH; i++)
      ghist[0] = 0;
    ptghist = 0;

    for (int i = 0; i < (1 << LOGSIZEUP); i++)
      Pupdatethreshold[i] = 35;

    for (int i = 0; i < GNB; i++)
      GGEHL[i] = &GGEHLA[i][0];
    for (int i = 0; i < LNB; i++)
      LGEHL[i] = &LGEHLA[i][0];

    for (int i = 0; i < GNB; i++)
      for (int j = 0; j < ((1 << LOGGNB) - 1); j++)
	{
	  if (!(j & 1))
	    {
	      GGEHL[i][j] = -1;

	    }
	}
    for (int i = 0; i < LNB; i++)
      for (int j = 0; j < ((1 << LOGLNB) - 1); j++)
	{
	  if (!(j & 1))
	    {
	      LGEHL[i][j] = -1;

	    }
	}
    for (int i = 0; i < PNB; i++)
      PGEHL[i] = &PGEHLA[i][0];
    for (int i = 0; i < PNB; i++)
      for (int j = 0; j < ((1 << LOGPNB) - 1); j++)
	{
	  if (!(j & 1))
	    {
	      PGEHL[i][j] = -1;

	    }
	}

#ifdef IMLI
#ifdef IMLIOH
    for (int i = 0; i < FNB; i++)
      FGEHL[i] = &FGEHLA[i][0];
    for (int i = 0; i < FNB; i++)
      for (int j = 0; j < ((1 << LOGFNB) - 1); j++)
	{
	  if (!(j & 1))
	    {
	      FGEHL[i][j] = -1;

	    }
	}
#endif
#ifdef IMLISIC
    for (int i = 0; i < INB; i++) {
      IGEHL[i] = &IGEHLA[i][0];
    }
    for (int i = 0; i < INB; i++)
      for (int j = 0; j < ((1 << LOGINB) - 1); j++)
	{
	  if (!(j & 1))
	    {
	      IGEHL[i][j] = -1;
	    }

	}
#endif    
#endif
#ifndef REALISTIC
    for (int i = 0; i < SNB; i++)
      SGEHL[i] = &SGEHLA[i][0];
    for (int i = 0; i < TNB; i++)
      TGEHL[i] = &TGEHLA[i][0];

    for (int i = 0; i < SNB; i++)
      for (int j = 0; j < ((1 << LOGSNB) - 1); j++)
	{
	  if (!(j & 1))
	    {
	      SGEHL[i][j] = -1;

	    }
	}
    for (int i = 0; i < TNB; i++)
      for (int j = 0; j < ((1 << LOGTNB) - 1); j++)
	{
	  if (!(j & 1))
	    {
	      TGEHL[i][j] = -1;

	    }
	}
#endif

    for (int i = 0; i < (1 << LOGB); i++)
      {
	btable[i].pred = 0;
	btable[i].hyst = 1;
      }


    for (int j = 0; j < (1 << (LOGBIAS + 1)); j++)
      Bias[j] = (j & 1) ? 15 : -16;
    for (int j = 0; j < (1 << (LOGBIAS + 1)); j++)
      BiasSK[j] = (j & 1) ? 15 : -16;


    for (int i = 0; i < NLOCAL; i++)
      {
	L_shist[i] = 0;
      }


    for (int i = 0; i < NSECLOCAL; i++)
      {
	S_slhist[i] = 0;

      }
    GHIST = 0;

    for (int i = 0; i < SIZEUSEALT; i++)
      {
	use_alt_on_na[i][0] = 0;
	use_alt_on_na[i][1] = 0;
      }

    TICK = 0;
    ptghist = 0;
    phist = 0;

  }




  // index function for the bimodal table

  int bindex (uint32_t PC)
  {
    return ((PC) & ((1 << (LOGB)) - 1));
  }


// the index functions for the tagged tables uses path history as in the OGEHL predictor
//F serves to mix path history: not very important impact

  int F (long long A, int size, int bank)
  {
    int A1, A2;
    A = A & ((1 << size) - 1);
    A1 = (A & ((1 << logg[bank]) - 1));
    A2 = (A >> logg[bank]);
    A2 =
      ((A2 << bank) & ((1 << logg[bank]) - 1)) + (A2 >> (logg[bank] - bank));
    A = A1 ^ A2;
    A = ((A << bank) & ((1 << logg[bank]) - 1)) + (A >> (logg[bank] - bank));
    return (A);
  }

// gindex computes a full hash of PC, ghist and phist
  int gindex (unsigned int PC, int bank, long long hist,
	      folded_history * ch_i)
  {
    int index;
    int M = (m[bank] > PHISTWIDTH) ? PHISTWIDTH : m[bank];
    index =
      PC ^ (PC >> (abs (logg[bank] - bank) + 1))
      ^ ch_i[bank].comp ^ F (hist, M, bank);
    return (index & ((1 << (logg[bank])) - 1));
  }

  //  tag computation
  uint16_t gtag (unsigned int PC, int bank, folded_history * ch0,
		 folded_history * ch1)
  {
    int tag = PC ^ ch0[bank].comp ^ (ch1[bank].comp << 1);
    return (tag & ((1 << TB[bank]) - 1));
  }

  // up-down saturating counter
  void ctrupdate (int8_t & ctr, bool taken, int nbits)
  {
    if (taken)
      {
	if (ctr < ((1 << (nbits - 1)) - 1))
	  ctr++;
      }
    else
      {
	if (ctr > -(1 << (nbits - 1)))
	  ctr--;
      }
  }


#ifdef LOOPPREDICTOR
  int lindex (uint32_t PC)
  {
//    return ((PC & ((1 << (LOGL - 2)) - 1)) << 2);
    return ((HASHPC & ((1 << (LOGL - 2)) - 1)) << 2);
  }


//loop prediction: only used if high confidence
//skewed associative 4-way
//At fetch time: speculative
#define CONFLOOP 15
  bool getloop (uint32_t PC)
  {
    LHIT = -1;

    LI = lindex (PC);
    LIB = ((PC >> (LOGL - 2)) & ((1 << (LOGL - 2)) - 1));
    LTAG = (PC >> (LOGL - 2)) & ((1 << 2 * LOOPTAG) - 1);
    LTAG ^= (LTAG >> LOOPTAG);
    LTAG = (LTAG & ((1 << LOOPTAG) - 1));

    for (int i = 0; i < 4; i++)
      {
	int index = (LI ^ ((LIB >> i) << 2)) + i;

	if (ltable[index].TAG == LTAG)
	  {
	    LHIT = i;
	    LVALID = ((ltable[index].confid == CONFLOOP)
		      || (ltable[index].confid * ltable[index].NbIter > 128));
	    {

	    }
	    if (ltable[index].CurrentIter + 1 == ltable[index].NbIter)
	      return (!(ltable[index].dir));
	    else
	      {
		
		return ((ltable[index].dir));
	      }

	  }
      }

    LVALID = false;
    return (false);

  }



  void loopupdate (uint32_t PC, bool Taken, bool ALLOC)
  {
    if (LHIT >= 0)
      {
	int index = (LI ^ ((LIB >> LHIT) << 2)) + LHIT;
//already a hit 
	if (LVALID)
	  {
	    if (Taken != predloop)
	      {
// free the entry
		ltable[index].NbIter = 0;
		ltable[index].age = 0;
		ltable[index].confid = 0;
		ltable[index].CurrentIter = 0;
		return;

	      }
	    else if ((predloop != tage_pred) || ((MYRANDOM () & 7) == 0))
	      if (ltable[index].age < CONFLOOP)
		ltable[index].age++;
	  }

	ltable[index].CurrentIter++;
	ltable[index].CurrentIter &= ((1 << WIDTHNBITERLOOP) - 1);
	//loop with more than 2** WIDTHNBITERLOOP iterations are not treated correctly; but who cares :-)
	if (ltable[index].CurrentIter > ltable[index].NbIter)
	  {
	    ltable[index].confid = 0;
	    ltable[index].NbIter = 0;
//treat like the 1st encounter of the loop 
	  }
	if (Taken != ltable[index].dir)
	  {
	    if (ltable[index].CurrentIter == ltable[index].NbIter)
	      {
		if (ltable[index].confid < CONFLOOP)
		  ltable[index].confid++;
		if (ltable[index].NbIter < 3)
		  //just do not predict when the loop count is 1 or 2     
		  {
// free the entry
		    ltable[index].dir = Taken;
		    ltable[index].NbIter = 0;
		    ltable[index].age = 0;
		    ltable[index].confid = 0;
		  }
	      }
	    else
	      {
		if (ltable[index].NbIter == 0)
		  {
// first complete nest;
		    ltable[index].confid = 0;
		    ltable[index].NbIter = ltable[index].CurrentIter;
		  }
		else
		  {
//not the same number of iterations as last time: free the entry
		    ltable[index].NbIter = 0;
		    ltable[index].confid = 0;
		  }
	      }
	    ltable[index].CurrentIter = 0;
	  }

      }
    else if (ALLOC)

      {
	uint32_t X = MYRANDOM () & 3;

	if ((MYRANDOM () & 3) == 0)
	  for (int i = 0; i < 4; i++)
	    {
	      int LHIT = (X + i) & 3;
	      int index = (LI ^ ((LIB >> LHIT) << 2)) + LHIT;
	      if (ltable[index].age == 0)
		{
		  ltable[index].dir = !Taken;
// most of mispredictions are on last iterations
		  ltable[index].TAG = LTAG;
		  ltable[index].NbIter = 0;
		  ltable[index].age = 7;
		  ltable[index].confid = 0;
		  ltable[index].CurrentIter = 0;
		  break;

		}
	      else
		ltable[index].age--;
	      break;
	    }
      }
  }
#endif

  bool getbim ()
  {
    BIM = (btable[BI].pred << 1) + (btable[BI >> HYSTSHIFT].hyst);
    HighConf = (BIM == 0) || (BIM == 3);
    return (btable[BI].pred > 0);
  }

  void baseupdate (bool Taken)
  {
    int inter = BIM;
    if (Taken)
      {
	if (inter < 3)
	  inter += 1;
      }
    else if (inter > 0)
      inter--;
    btable[BI].pred = inter >> 1;
    btable[BI >> HYSTSHIFT].hyst = (inter & 1);
  };

//just a simple pseudo random number generator: use available information
// to allocate entries  in the loop predictor
  int MYRANDOM ()
  {
    Seed++;
    Seed ^= phist;
    Seed = (Seed >> 21) + (Seed << 11);
    return (Seed);
  };


  //  TAGE PREDICTION: same code at fetch or retire time but the index and tags must recomputed
  void Tagepred (UINT32 PC)
  {
    HitBank = 0;
    AltBank = 0;
    for (int i = 1; i <= NHIST; i++)
      {
	GI[i] = gindex (PC, i, phist, ch_i);
	GTAG[i] = gtag (PC, i, ch_t[0], ch_t[1]);
      }

//    BI = PC & ((1 << LOGB) - 1);
    BI = HASHPC & ((1 << LOGB) - 1);


//Look for the bank with longest matching history
    for (int i = NHIST; i > 0; i--)
      {
	if (gtable[i][GI[i]].tag == GTAG[i])
	  {
	    HitBank = i;
	    LongestMatchPred = (gtable[HitBank][GI[HitBank]].ctr >= 0);
	    break;
	  }
      }

//Look for the alternate bank
    for (int i = HitBank - 1; i > 0; i--)
      {
	if (gtable[i][GI[i]].tag == GTAG[i])
	  {

	    AltBank = i;
	    break;
	  }
      }
//computes the prediction and the alternate prediction

    if (HitBank > 0)
      {
	if (AltBank > 0)
	  alttaken = (gtable[AltBank][GI[AltBank]].ctr >= 0);
	else
	  alttaken = getbim ();

//if the entry is recognized as a newly allocated entry and 
//USE_ALT_ON_NA is positive  use the alternate prediction
	int index = INDUSEALT ^ LongestMatchPred;
	bool Huse_alt_on_na =
	  (use_alt_on_na[index][HitBank > (NHIST / 3)] >= 0);

	if ((!Huse_alt_on_na)
	    || (abs (2 * gtable[HitBank][GI[HitBank]].ctr + 1) > 1))
	  tage_pred = LongestMatchPred;
	else
	  tage_pred = alttaken;

	HighConf = (abs (2 * gtable[HitBank][GI[HitBank]].ctr + 1) >= (1 << CWIDTH) - 1);
      }
    else
      {
	alttaken = getbim ();
	tage_pred = alttaken;
	LongestMatchPred = alttaken;
      }



  }
//compute the prediction

  bool GetPrediction (UINT32 PC)
  {
 if (punting) return false; 

// computes the TAGE table addresses and the partial tags

    Tagepred (PC);
    pred_taken = tage_pred;

#ifdef LOOPPREDICTOR
    predloop = getloop (PC);	// loop prediction
    pred_taken = ((WITHLOOP >= 0) && (LVALID)) ? predloop : pred_taken;  
#endif

    pred_inter = pred_taken;
#ifndef SC
    //just the TAGE predictor
    return(pred_taken);

#endif    
//Compute the SC prediction

// begin to bias the sum towards TAGE predicted direction
    LSUM = 22; // djimenez
    if (!pred_inter) LSUM = -LSUM;

// djimenez
    LSUM += compute_partial_sum (PC);

////////////////////////////////////
//integrate BIAS prediction   
#if 0
    int8_t ctr = Bias[INDBIAS];
    LSUM += (2 * ctr + 1);
    ctr = BiasSK[INDBIASSK];
    LSUM += (2 * ctr + 1);
#endif
	// djimenez
    int8_t ctr = Bias[INDBIAS];
    LSUM += 2.09 * ctr;
    ctr = BiasSK[INDBIASSK];
    LSUM += 2.08 * ctr;

//integrate the GEHL predictions
#ifdef IMLI
#ifdef IMLIOH
    localoh = 0;
     localoh = PIPE[(PC ^ (PC >> 4)) & (PASTSIZE - 1)] + (localoh << 1);
    for (int i = 0; i >= 0; i--)
      {
	localoh =
	  ohhisttable[(((PC ^ (PC >> 4)) << SHIFTFUTURE) + IMLIcount +
		     i) & (OHHISTTABLESIZE - 1)] + (localoh << 1);
      }


    if (IMLIcount >= 2)
      LSUM += 2 * Gpredict ((PC << 2), localoh, Fm, FGEHL, FNB, LOGFNB);
#endif

#ifdef IMLISIC
    LSUM += Gpredict (PC, IMLIcount, Im, IGEHL, INB, LOGINB);
#else
    long long interIMLIcount= IMLIcount;
/* just a trick to disable IMLIcount*/
    IMLIcount=0;
#endif
#endif

#ifndef IMLI
    IMLIcount = 0;
    
//just  a trick to simplify the programming
#endif

    LSUM +=
      Gpredict ((PC << 1) + pred_inter /*PC*/, (GHIST << 11) + IMLIcount, Gm, GGEHL, GNB, LOGGNB);

#ifdef LOCALH
    LSUM += 2.02 * Gpredict (PC, L_shist[INDLOCAL], Lm, LGEHL, LNB, LOGLNB); // djimenez
    if (L_shist[INDLOCAL] == 2047) LSUM += 4;  // djimenez
    if (L_shist[INDLOCAL] == 0) LSUM += -4;  // djimenez
#ifndef REALISTIC
    LSUM +=
      Gpredict (PC, (T_slhist[INDTLOCAL] << 11) + IMLIcount, Tm, TGEHL, TNB,
		LOGTNB);
    LSUM += Gpredict (PC, S_slhist[INDSLOCAL], Sm, SGEHL, SNB, LOGSNB);

#endif
#endif
#ifdef IMLI    
#ifndef IMLISIC
IMLIcount= interIMLIcount;
#endif
#endif

#ifdef REALISTIC
    LSUM += Gpredict (PC, GHIST, Pm, PGEHL, PNB, LOGPNB);
#else
    LSUM += Gpredict (PC, HSTACK[pthstack], Pm, PGEHL, PNB, LOGPNB);
#endif

    bool SCPRED = (LSUM >= 0);

// chose between the SC output and the TAGE + loop  output

    if (pred_inter != SCPRED)
      {
//Choser uses TAGE confidence and |LSUM|
	pred_taken = SCPRED;
	if (HighConf)
	  {
	    if ((abs (LSUM) < Pupdatethreshold[INDUPD] / 3))
	      pred_taken = (FirstH < 0) ? SCPRED : pred_inter;

	    else if ((abs (LSUM) < 2 * Pupdatethreshold[INDUPD] / 3))
	      pred_taken = (SecondH < 0) ? SCPRED : pred_inter;
	    else if ((abs (LSUM) < Pupdatethreshold[INDUPD]))
	      pred_taken = (ThirdH < 0) ? SCPRED : pred_inter;
	  }

      }

    return pred_taken;
  }

  void HistoryUpdate (uint32_t PC, uint16_t brtype, bool taken,
		      uint32_t target, long long &X, int &Y,
		      folded_history * H, folded_history * G,
		      folded_history * J, long long &LH,
		      long long &SH, long long &TH, long long &PH,
		      long long &GBRHIST)
  {
//special treatment for unconditional branchs;
    int maxt;
    if (brtype == OPTYPE_JMP_DIRECT_COND||brtype == OPTYPE_JMP_INDIRECT_COND||brtype==OPTYPE_RET_COND||brtype==OPTYPE_CALL_DIRECT_COND||brtype==OPTYPE_CALL_INDIRECT_COND)
      maxt = 1;
    else
      maxt = 4;

    // djimenez
	if (brtype == OPTYPE_JMP_DIRECT_COND) {
		if (target < PC) {
			if (!taken) {
				imli_counter1 = 0;
			} else {
				imli_counter1++;
			}
		} else {
			if (taken) {
				imli_counter4 = 0;
			} else {
				imli_counter4++;
			}
		}
	}

    // the return stack associated history

#ifdef IMLI
    // this is the only kind of branch that is likely to be useful with IMLI
    if (brtype == OPTYPE_JMP_DIRECT_COND) {
	if (target < PC) {
		//This branch is a branch "loop" 
		if (!taken) {
			//exit of the "loop"
			IMLIcount = 0;
		}
		if (taken) {
			if (IMLIcount < (MAXIMLIcount))
			IMLIcount++;
		}
	}
    }

#ifdef IMLIOH
    if (IMLIcount >= 1)
      // if (brtype == OPTYPE_JMP_DIRECT_COND||brtype == OPTYPE_JMP_INDIRECT_COND||brtype==OPTYPE_RET_COND||brtype==OPTYPE_CALL_DIRECT_COND||brtype==OPTYPE_CALL_INDIRECT_COND)
      if (brtype == OPTYPE_JMP_DIRECT_COND)
	{
	  if (target >= PC)
	    {
	      PIPE[(PC ^ (PC >> 4)) & (PASTSIZE - 1)] =
		ohhisttable[(((PC ^ (PC >> 4)) << SHIFTFUTURE) +
			   IMLIcount) & (OHHISTTABLESIZE - 1)];
	      ohhisttable[(((PC ^ (PC >> 4)) << SHIFTFUTURE) +
			 IMLIcount) & (OHHISTTABLESIZE - 1)] = taken;
	    }

	}
#endif    
#endif

    if (brtype == OPTYPE_JMP_DIRECT_COND||brtype == OPTYPE_JMP_INDIRECT_COND||brtype==OPTYPE_RET_COND||brtype==OPTYPE_CALL_DIRECT_COND||brtype==OPTYPE_CALL_INDIRECT_COND)
      {
	GBRHIST = (GBRHIST << 1) + taken;
	LH = (LH << 1) + (taken);
#ifndef REALISTIC
	SH = (SH << 1) + (taken);
//	SH ^= (PC & 15);
	SH ^= (HASHPC & 15);
	TH = (TH << 1) + (taken);
#endif
      }
#ifndef REALISTIC
    PH = (PH << 1) ^ (target ^ (target >> 5) ^ taken);
    if (brtype == OPTYPE_RET_COND || brtype == OPTYPE_RET_UNCOND)
      {
	pthstack = (pthstack - 1) & 15;
      }

    if (brtype == OPTYPE_CALL_DIRECT_COND||brtype == OPTYPE_CALL_INDIRECT_COND||brtype == OPTYPE_CALL_DIRECT_UNCOND||brtype == OPTYPE_CALL_DIRECT_UNCOND)
      {
	int index = (pthstack + 1) & 15;
	HSTACK[index] = HSTACK[pthstack];
	pthstack = index;
      }
#endif

    int T = ((PC) << 1) + taken;
    int PATH = PC;


    for (int t = 0; t < maxt; t++)
      {
	bool DIR = (T & 1);
	T >>= 1;
	int PATHBIT = (PATH & 127);
	PATH >>= 1;
//update  history
	Y--;
	ghist[Y & (HISTBUFFERLENGTH - 1)] = DIR;
	X = (X << 1) ^ PATHBIT;
	for (int i = 1; i <= NHIST; i++)
	  {

	    H[i].update (ghist, Y);
	    G[i].update (ghist, Y);
	    J[i].update (ghist, Y);


	  }
      }


//END UPDATE  HISTORIES
  }

// PREDICTOR UPDATE

  void UpdatePredictor (UINT32 PC, bool resolveDir, bool predDir,
			UINT32 branchTarget)
  {
if (punting) return;

#ifdef LOOPPREDICTOR
    if (LVALID)
      {

	if (pred_taken != predloop)
	  ctrupdate (WITHLOOP, (predloop == resolveDir), 7);
      }

    loopupdate (PC, resolveDir, (pred_taken != resolveDir));
#endif
#ifdef SC
    bool SCPRED = (LSUM >= 0);
    if (pred_inter != SCPRED)
      {
	if ((abs (LSUM) < Pupdatethreshold[INDUPD]))
	  if ((HighConf))
	    {

	      if ((abs (LSUM) < Pupdatethreshold[INDUPD] / 3))
		ctrupdate (FirstH, (pred_inter == resolveDir), CONFWIDTH);
	      else if ((abs (LSUM) < 2 * Pupdatethreshold[INDUPD] / 3))
		ctrupdate (SecondH, (pred_inter == resolveDir), CONFWIDTH);
	      else if ((abs (LSUM) < Pupdatethreshold[INDUPD]))
		ctrupdate (ThirdH, (pred_inter == resolveDir), CONFWIDTH);

	    }
      }

    if ((SCPRED != resolveDir) || ((abs (LSUM) < Pupdatethreshold[INDUPD])))
      {
	update_partial (PC, resolveDir);
	{
	  if (SCPRED != resolveDir)
	    Pupdatethreshold[INDUPD] += 1;
	  else
	    Pupdatethreshold[INDUPD] -= 1;

	  if (Pupdatethreshold[INDUPD] >= 1024) // djimenez; this was 256
	    Pupdatethreshold[INDUPD] = 1023;
          
	  if (Pupdatethreshold[INDUPD] < 0)
	    Pupdatethreshold[INDUPD] = 0;
	}

	ctrupdate (Bias[INDBIAS], resolveDir, PERCWIDTH);
	ctrupdate (BiasSK[INDBIASSK], resolveDir, PERCWIDTH);
#ifdef IMLI
#ifndef IMLISIC
    long long interIMLIcount= IMLIcount;
/* just a trick to disable IMLIcount*/
    IMLIcount=0;
#endif
#endif
        Gupdate ((PC << 1) + pred_inter /*PC*/, resolveDir, (GHIST << 11) + IMLIcount,
		 Gm, GGEHL, GNB, LOGGNB);
	Gupdate (PC, resolveDir, L_shist[INDLOCAL], Lm, LGEHL, LNB, LOGLNB);
#ifdef REALISTIC
	Gupdate (PC, resolveDir, GHIST, Pm, PGEHL, PNB, LOGPNB);
#endif

#ifndef REALISTIC
	Gupdate (PC, resolveDir, S_slhist[INDSLOCAL], Sm, SGEHL, SNB, LOGSNB);
	Gupdate (PC, resolveDir, (T_slhist[INDTLOCAL] << 11) + IMLIcount, Tm,
		 TGEHL, TNB, LOGTNB);
	Gupdate (PC, resolveDir, HSTACK[pthstack], Pm, PGEHL, PNB, LOGPNB);
#endif

#ifdef IMLI
#ifdef IMLISIC
        Gupdate (PC, resolveDir, IMLIcount, Im, IGEHL, INB, LOGINB);
#else
       IMLIcount=interIMLIcount;
#endif
#ifdef IMLIOH

        if (IMLIcount >= 2)
	  Gupdate ((PC << 2), resolveDir, localoh, Fm, FGEHL, FNB, LOGFNB);
#endif

#endif
      }
//ends update of the SC states
#endif

//TAGE UPDATE
    bool ALLOC = ((tage_pred != resolveDir) & (HitBank < NHIST));
    if (pred_taken == resolveDir)
        if ((MYRANDOM () & 31) != 0)
	ALLOC = false;
    //do not allocate too often if the overall prediction is correct 

    if (HitBank > 0)
      {
// Manage the selection between longest matching and alternate matching
// for "pseudo"-newly allocated longest matching entry
	bool PseudoNewAlloc =
	  (abs (2 * gtable[HitBank][GI[HitBank]].ctr + 1) <= 1);
// an entry is considered as newly allocated if its prediction counter is weak
	if (PseudoNewAlloc)
	  {
	    if (LongestMatchPred == resolveDir)
	      ALLOC = false;
// if it was delivering the correct prediction, no need to allocate a new entry
//even if the overall prediction was false
	    if (LongestMatchPred != alttaken)
	      {
		int index = (INDUSEALT) ^ LongestMatchPred;
		ctrupdate (use_alt_on_na[index]
			   [HitBank > (NHIST / 3)],
			   (alttaken == resolveDir), 4);

	      }

	  }
      }


    if (ALLOC)
      {

	int T = NNN;

	int A = 1;
	if ((MYRANDOM () & 127) < 32)
	  A = 2;
	int Penalty = 0;
	int NA = 0;
	for (int i = HitBank + A; i <= NHIST; i += 1)
	  {
	    if (gtable[i][GI[i]].u == 0)
	      {
		gtable[i][GI[i]].tag = GTAG[i];
		gtable[i][GI[i]].ctr = (resolveDir) ? 0 : -1;
		NA++;
		if (T <= 0)
		  {
		    break;
		  }
		i += 1;
		T -= 1;
	      }
	    else
	      {
		Penalty++;
	      }
	  }
	TICK += (Penalty - NA);
//just the best formula for the Championship
	if (TICK < 0)
	  TICK = 0;
	if (TICK > 1023)
	  {
	    for (int i = 1; i <= NHIST; i++)
	      for (int j = 0; j <= (1 << logg[i]) - 1; j++)
// substracting 1 to a whole array is not that realistic
#ifdef REALISTIC
		gtable[i][j].u >>= 1;
#else
		if (gtable[i][j].u)
		  gtable[i][j].u--;
#endif
	    TICK = 0;

	  }


      }
//update predictions
    if (HitBank > 0)
      {
	if (abs (2 * gtable[HitBank][GI[HitBank]].ctr + 1) == 1)
	  if (LongestMatchPred != resolveDir)

	    {			// acts as a protection 
	      if (AltBank > 0)
		{
		  ctrupdate (gtable[AltBank][GI[AltBank]].ctr,
			     resolveDir, CWIDTH);

		}
	      if (AltBank == 0)
		baseupdate (resolveDir);
	    }
	ctrupdate (gtable[HitBank][GI[HitBank]].ctr, resolveDir, CWIDTH);
//sign changes: no way it can have been useful
	if (abs (2 * gtable[HitBank][GI[HitBank]].ctr + 1) == 1)
	  gtable[HitBank][GI[HitBank]].u = 0;

      }
    else
      baseupdate (resolveDir);

    if (LongestMatchPred != alttaken)
      if (LongestMatchPred == resolveDir)
	{
	  if (gtable[HitBank][GI[HitBank]].u < (1 << UWIDTH) - 1)
	    gtable[HitBank][GI[HitBank]].u++;
	}
//END TAGE UPDATE

    // djimenez
    update_histories (PC, resolveDir);
    HistoryUpdate (PC, OPTYPE_JMP_DIRECT_COND, resolveDir, branchTarget, phist,
		   ptghist, ch_i, ch_t[0],
		   ch_t[1], L_shist[INDLOCAL],
		   S_slhist[INDSLOCAL], T_slhist[INDTLOCAL], HSTACK[pthstack],
		   GHIST);
//END PREDICTOR UPDATE

  }
#define GINDEX (((long long) PC) ^ bhist ^ (bhist >> (8 - i)) ^ (bhist >> (16 - 2 * i)) ^ (bhist >> (24 - 3 * i)) ^ (bhist >> (32 - 3 * i)) ^ (bhist >> (40 - 4 * i))) & ((1 << (logs - (i >= (NBR - 2)))) - 1)

  int Gpredict (UINT32 PC, long long BHIST, int *length, int8_t ** tab,
		int NBR, int logs)
  {

    int PERCSUM = 0;
    for (int i = 0; i < NBR; i++)
      {
	long long bhist = BHIST & ((long long) ((1 << length[i]) - 1));
	int16_t ctr = tab[i][GINDEX];
	PERCSUM += (2 * ctr + 1);
      }
    return ((PERCSUM));

  }
  void Gupdate (UINT32 PC, bool taken, long long BHIST, int *length,
		int8_t ** tab, int NBR, int logs)
  {


    for (int i = 0; i < NBR; i++)
      {
	long long bhist = BHIST & ((long long) ((1 << length[i]) - 1));
	ctrupdate (tab[i][GINDEX], taken, PERCWIDTH- (i < (NBR-1)));
      }



  }


  void TrackOtherInst (UINT32 PC, OpType opType, UINT32 branchTarget)
  {

    bool taken = true;

    switch (opType)
      {
      case OPTYPE_CALL_DIRECT_UNCOND:
      case OPTYPE_CALL_INDIRECT_UNCOND:
      case OPTYPE_JMP_DIRECT_UNCOND:
      case OPTYPE_JMP_INDIRECT_UNCOND:
      case OPTYPE_RET_UNCOND:
	HistoryUpdate (PC, opType, taken, branchTarget, phist,
		       ptghist, ch_i,
		       ch_t[0], ch_t[1],
		       L_shist[INDLOCAL],
		       S_slhist[INDSLOCAL], T_slhist[INDTLOCAL],
		       HSTACK[pthstack], GHIST);
	break;


      default:;
      }


  }

};

}; // namespace

class multiperspective_perceptron_withTAGEbig : public branch_predictor {
        branch_info u;
        multiperspective_perceptron_withTAGEbig_namespace::PREDICTOR *p;
        unsigned int oldpc;

public:
        multiperspective_perceptron_withTAGEbig (multiperspective_perceptron_withTAGEbig_namespace::history_spec *_specs = NULL, int _num_tables = 0, unsigned int _imli_mask1 = 0x70 , unsigned int _imli_mask4 = 0) {
		multiperspective_perceptron_withTAGEbig_namespace::num_tables = _num_tables;
		multiperspective_perceptron_withTAGEbig_namespace::specs = _specs;
		multiperspective_perceptron_withTAGEbig_namespace::imli_mask1 = _imli_mask1;
		multiperspective_perceptron_withTAGEbig_namespace::imli_mask4 = _imli_mask4;
	
		if (!_specs) {
			multiperspective_perceptron_withTAGEbig_namespace::specs = multiperspective_perceptron_withTAGEbig_namespace::tuned_specs;
			multiperspective_perceptron_withTAGEbig_namespace::num_tables = sizeof (multiperspective_perceptron_withTAGEbig_namespace::tuned_specs) / sizeof (multiperspective_perceptron_withTAGEbig_namespace::history_spec);
			fflush (stdout);
		}
                p = new multiperspective_perceptron_withTAGEbig_namespace::PREDICTOR ();
        }

        branch_info *lookup (unsigned int pc) {
		oldpc = pc;
                u.prediction (p->GetPrediction (pc));
                return &u;
        }

        void update (branch_info *pp, unsigned int target, bool taken, int type) {
                p->UpdatePredictor (oldpc, taken, pp->prediction(), target);
        }

        void nonconditional_branch (unsigned int pc, unsigned int target, int type) {
                p->TrackOtherInst (pc, (OpType) type, target);
        }
};

PREDICTOR::PREDICTOR(void){
	ex = new multiperspective_perceptron_withTAGEbig ();
}

bool PREDICTOR::GetPrediction(UINT64 PC){
	u = ex->lookup (PC);
	return u->prediction ();
}

void PREDICTOR::UpdatePredictor(UINT64 PC, OpType opType, bool resolveDir, bool predDir, UINT64 branchTarget){
	ex->update (u, branchTarget, resolveDir, opType);
}

void PREDICTOR::TrackOtherInst(UINT64 PC, OpType opType, bool branchDir, UINT64 branchTarget){
	ex->nonconditional_branch (PC, branchTarget, opType);
}

#endif
