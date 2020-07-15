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

	virtual void nonconditional_branch (unsigned int pc) { }

	virtual const char *name (bool) { return ""; }

	virtual ~branch_predictor (void) { }
};

class multiperspective_perceptron;

class PREDICTOR {
	multiperspective_perceptron *ex;
	branch_info *u;
public:

	PREDICTOR(void);

	bool GetPrediction (UINT64 PC);  
	void UpdatePredictor (UINT64 PC, OpType opType, bool resolveDir, bool predDir, UINT64 branchTarget);
	void TrackOtherInst (UINT64 PC, OpType opType, bool branchDir, UINT64 branchTarget);
};

//#define VERBOSE

#define MAX_TABLES	37

typedef char weight;

class multiperspective_perceptron_update : public branch_info {
public:
	unsigned int pc;
	unsigned short int hpc, pc2;
	int yout; 
	int indices[MAX_TABLES];
	bool filtered;
	multiperspective_perceptron_update (void) { }
};

#if 1
struct bestpair {
	int index;
	int mpreds;
};

// feature types
enum history_type {
        GHIST = 1,              // global pattern history
        ACYCLIC = 2,            // "acylic history" - put the most recent history of pc into array[pc%parameter] and use this array as a feature
        MODHIST = 3,            // "modulo history" - shift history bit in if PC%modulus==0
        BIAS = 4,               // bias of this branch
        RECENCY = 5,            // hash of a recency stack of PCs
        IMLI = 6,               // inner most loop iteration counter(s)
        PATH = 7,               // path history
        LOCAL = 8,              // local history (pshare)
        MODPATH = 9,            // like modhist but with path history
        GHISTPATH = 10,         // (path history << 1) | global history
        GHISTMODPATH = 11,	// (mod path history << 1) | global history
        BLURRYPATH = 12,        // "page" history sort of
        RECENCYPOS = 13,        // position of this PC in the recency stack
        SGHISTPATH = 14,        // alternate ghistpath inspired by strided sampling
        MAXTYPE = 16
};

// stupid C++
#ifndef _MAIN_CC
int paircmp (const bestpair *a, const bestpair *b) {
	return a->mpreds - b->mpreds;
}
#else
int paircmp (const bestpair *a, const bestpair *b);
#endif

#ifndef _MAIN_CC

const char *type_names[] = {
"nothing", "GHIST", "ACYCLIC", "MODHIST", "BIAS", "RECENCY", "IMLI", "PATH", "LOCAL", "MODPATH", "GHISTPATH", "GHISTMODPATH", "BLURRYPATH", "RECENCYPOS", "SGHISTPATH",
};
#endif

// this struct represents one feature that is an input to the hashed
// perceptron predictor

struct history_spec {
	// the type of the feature

	enum history_type type;

	// up to three parameters for this feature

	int p1, p2, p3;

	// a coefficient multiplied by the weight value

	double coeff;

	// the size of the array of weights for this feature
	int 
		size, 

	// bit width of weights for this feature (5 or 6)
		width;
};

#ifdef _BRANCHY_CC
#include "spec.h"
#endif

// an entry in the "branch filter", a structure that keeps track of whether
// a branch has been take or not taken so far

struct filter_entry {
	bool	seen_taken, seen_untaken;

	filter_entry (void) {
		seen_taken = false;
		seen_untaken = false;
	}
};
#endif

// main class for the prediction code, contained in Jimenez's little branch
// prediction infrastructure

class multiperspective_perceptron : public branch_predictor {
public:

	// some values carried from prediction to update. there's no good
	// reason that some stuff is put in here and other stuff isn't. it
	// just happened that way as a historical artifact of supporting
	// speculative update for previous predictors.

	multiperspective_perceptron_update u;

	// for simplicity we keep stuff in large arrays and then only use
	// the part of them needed for the predictor. these defines give the
	// maximum sizes of the arrays. they are very large to allow for flexible
	// search of the design space, but the actual tuned values are represented
	// in variables passed as parameters to the constructor and fit within the
	// given hardware budgets

#define MAX_PATHHIST	4096
#define MAX_GHIST	4096
#define MAX_LG_TABLE_SIZE	21
#define MAX_TABLE_SIZE	(1<<MAX_LG_TABLE_SIZE)
#define MAX_FILTER	262144
#define MAX_LOCAL_HISTORIES	262144
#define MAX_ACYCLIC	20
#define MAX_MOD		10
#define MAX_BLURRY	16
#define MAX_BLURRY2	16
#define MAX_ASSOC	256

	// STATE THAT COUNTS AGAINST THE HARDWARE BUDGET. these variables
	// are the part of the predictor that contain mutable state that
	// persists from one prediction to the next. the parts of these variables
	// that are actually used by the predictor count against the hardware budget.

	// seed for random number generator

	unsigned int my_seed;

	// branch filter

	filter_entry filter_table[MAX_FILTER];

	// table of per-branch (local) histories

	unsigned int local_histories[MAX_LOCAL_HISTORIES];

	// counter for adaptive theta setting

	int tc;

	int occupancy;

	// tables of weight magnitudes

	weight tables[MAX_TABLES][MAX_TABLE_SIZE];

	// tables of weight signs

	bool sign_bits[MAX_TABLES][MAX_TABLE_SIZE][2];

	// global history bits divided into block_size bits per array element

	unsigned int 
		ghist_words[MAX_GHIST/MAX_LG_TABLE_SIZE+1];

	// count of mispredictions per table

	int mpreds[MAX_TABLES];

	// acyclic histories
	
	bool acyclic_histories[MAX_ACYCLIC][32];

	// alternate acyclic histories that use path instead of pattern history

	unsigned int acyclic2_histories[MAX_ACYCLIC][32];

	// modulo pattern history
	bool mod_histories[MAX_MOD][MAX_GHIST];

	// modulo path history

	unsigned short int modpath_histories[MAX_MOD][MAX_PATHHIST];

	// "page" history

	unsigned int blurrypath_histories[MAX_BLURRY][MAX_BLURRY2];

	// recency stack of recenctly visited PCs (hashed)

	unsigned int short recency_stack[MAX_ASSOC];

	// the last global history outcome

	bool last_ghist_bit;

	// history of recent PCs (hashed)

	unsigned int short path_history[MAX_PATHHIST];

	// innermost loop iteration counters; 4 versions

	unsigned int imli_counter1, imli_counter2, imli_counter3, imli_counter4;

	// RUN-TIME CONSTANTS THAT DO NOT COUNT AGAINST HARDWARE BUDGET

	int

		// many different moduli could be used for building
		// modulo history. these arrays keep track of which moduli
		// are actually used in the features for this predictor
		// so only the relevant tables are updated for performance
		// reasons. in a real implementation, only those tables would
		// be implemented.

		modhist_indices[MAX_MOD],
		modpath_indices[MAX_MOD],	
		modpath_lengths[MAX_MOD], 
		modhist_lengths[MAX_MOD];

	int 
		// maximum global history required by any feature

		ghist_length, 

		// maximum modulo history required by any feature

		modghist_length, 

		// maximum path length required by any feature

		path_length,

		// associativity of the recency stack, derived from the
		// maximum depth any RECENCY feature hashes

		assoc;

	long long int

		// total bits used by the predictor (for sanity check)

		totalbits;

	// number of modulo histories and modulo path histories

	int nmodhist_histories, nmodpath_histories;

	// insert an item into an array without duplication. used for
	// building the _indices arrays

	int insert (int *v, int *n, int x) {
		for (int i=0; i<*n; i++) if (v[i] == x) return i;
		v[(*n)] = x;
		return (*n)++;
	}

	// analyze the specification of the features to see what the
	// various maximum history lengths are, how much space is taken by
	// various history structures, and how much space is left over in
	// the hardware budget for tables

	void analyze_spec (void) {
		bool 
			// true if at least one feature requires the recency stack

			doing_recency = false, 

			// true if at least one feature uses local history

			doing_local = false;

	
		int 
			// how many bits are allocated to the IMLI counters

			imli_counter_bits[4];

		// initially assume a recency stack of depth 0

		assoc = 0;

		// accounting for bits for the blurry path history

		int blurrypath_bits[MAX_BLURRY][MAX_BLURRY2];

		// accouting for the bits for the acyclic path history

		bool acyclic_bits[MAX_ACYCLIC][32][2];

		// set these accounting arrays to 0 initially

		memset (blurrypath_bits, 0, sizeof (blurrypath_bits));
		memset (imli_counter_bits, 0, sizeof (imli_counter_bits));
		memset (acyclic_bits, 0, sizeof (acyclic_bits));

		// initially assume very short histories, later find out
		// what the maximum values are from the features

		ghist_length = 1;
		modghist_length = 1;
		nmodhist_histories = 0;
		path_length = 1;

		// go through each feature in the specification finding the requirements

		for (int i=0; i<num_tables; i++) {
			// find the maximum associativity of the recency stack required

			if (specs[i].type == RECENCY || specs[i].type == RECENCYPOS) {
				if (assoc < specs[i].p1) assoc = specs[i].p1;
			}

			// find out how much and what kind of history is needed for acyclic feature

			if (specs[i].type == ACYCLIC) {
				for (int j=0; j<specs[i].p1+2; j++) {
					acyclic_bits[specs[i].p1][j][!specs[i].p3] = true;
				}
			}

			// do we need local history?

			if (specs[i].type == LOCAL) doing_local = true;

			// how many IMLI counter bits do we need for different versions of IMLI

			if (specs[i].type == IMLI) imli_counter_bits[specs[i].p1-1] = 32;

			// do we require a recency stack?

			if (specs[i].type == RECENCY || specs[i].type == RECENCYPOS) doing_recency = true;

			// count blurry path bits (assuming full 32-bit addresses less shifted bits)

			if (specs[i].type == BLURRYPATH) for (int j=0; j<specs[i].p2; j++) blurrypath_bits[specs[i].p1][j] = 32 - specs[i].p1;

			// if we are doing modulo history, figure out which and how much history we need

			if (specs[i].type == MODHIST || specs[i].type == GHISTMODPATH) {
				int j = insert (modhist_indices, &nmodhist_histories, specs[i].p1);
				if (modhist_lengths[j] < specs[i].p2+1) modhist_lengths[j] = specs[i].p2 + 1;
				if (specs[i].p2 >= modghist_length) modghist_length = specs[i].p2+1;
			}
		}

		// figure out how much history we need for modulo path, modulo+global path, and regular path

		nmodpath_histories = 0;
		for (int i=0; i<num_tables; i++) {
			if (specs[i].type == MODPATH || specs[i].type == GHISTMODPATH) {
				int j = insert (modpath_indices, &nmodpath_histories, specs[i].p1);
				if (modpath_lengths[j] < specs[i].p2+1) modpath_lengths[j] = specs[i].p2 + 1;
				if (path_length <= specs[i].p2) path_length = specs[i].p2 + 1;
			}
			if (specs[i].type == PATH) {
				if (path_length <= specs[i].p1) path_length = specs[i].p1 + 1;
			}
		}

		// how much global history and global path history do we need

		for (int i=0; i<num_tables; i++) {
			switch (specs[i].type) {
			case GHIST:
			if (ghist_length <= specs[i].p2) ghist_length = specs[i].p2 + 1;
			break;
			case GHISTPATH:
			if (ghist_length <= specs[i].p1) ghist_length = specs[i].p1 + 1;
			if (path_length <= specs[i].p1) path_length = specs[i].p1 + 1;
			break;
			default: ;
			}
		}

		// sanity check

		assert (ghist_length <= MAX_GHIST);
		assert (modghist_length <= MAX_GHIST);

		// account for IMLI counters

		for (int i=0; i<4; i++) totalbits += imli_counter_bits[i];

		// account for ghist bits

		totalbits += ghist_length;

		// account for global path bits (represented as an array of 16 bit integers)
	
		totalbits += path_length * 16;

		// account for misprediction monitoring counters if the "threshold" idea is being used

		if (threshold >= 0) totalbits += tunebits * num_tables;

		// plus one bit to record each prediction (can avoid by recomputing)
		// totalbits += num_tables;

		// account for modulo history bits

		for (int i=0; i<nmodhist_histories; i++) totalbits += modhist_lengths[i];

		// account for modulo path bits

		for (int i=0; i<nmodpath_histories; i++) totalbits += 16 * modpath_lengths[i];

		// account for local histories

		if (doing_local) totalbits += local_history_length * nlocal_histories;

		// account for recency stack

		if (doing_recency) totalbits += assoc * 16;

		// account for blurry path bits

		for (int i=0; i<MAX_BLURRY; i++) for (int j=0; j<MAX_BLURRY2; j++) totalbits += blurrypath_bits[i][j];

		// account for filter bits

		totalbits += num_filter * 2;

		// account for acyclic bits

		for (int i=0; i<MAX_ACYCLIC; i++) for (int j=0; j<32; j++) for (int k=0; k<2; k++) totalbits += acyclic_bits[i][j][k];

		// how many bits are left for the tables?

		long long int remaining = budgetbits - totalbits;

		// count the tables that have already been assigned sizes

		int num_sized = 0;
		for (int i=0; i<num_tables; i++) {
			if (table_sizes[i] != 0) {
				int sz = table_sizes[i] * (specs[i].width + (n_sign_bits-1));
				totalbits += sz;
				remaining -= sz;
				num_sized++;
			}
		}

		// whatever is left, we divide among the rest of the tables

		long long int table_size_bits = (remaining / (num_tables-num_sized));
		for (int i=0; i<num_tables; i++) {

			// if a table doesn't have a size yet, give it one and count those bits

			if (!table_sizes[i]) {
				int my_table_size = table_size_bits / (specs[i].width + (n_sign_bits-1)); // extra sign bits
				table_sizes[i] = my_table_size;
				totalbits += my_table_size * (specs[i].width + (n_sign_bits-1));
			}
		}
#ifdef VERBOSE
		printf ("%lld bits of metadata so far, %lld left out of %lld total budget\n", totalbits, remaining, budgetbits);
		printf ("table size is %lld bits, ", table_size_bits);
		int my_table_size = table_size_bits / (5 + (n_sign_bits-1));
		printf ("%d entries for 5 bit, ", my_table_size);
		my_table_size = table_size_bits / (6 + (n_sign_bits-1));
		printf ("%d entries for 6 bit\n", my_table_size);
		printf ("%lld total bits (%0.2fKB)\n", totalbits, totalbits / 8192.0);
#endif
	}


	// these variables set from parameters

	history_spec *specs;
	int 
		num_tables; 		// number of features (each one gets its own table)

	long long int
		budgetbits; 		// hardware budget in bits

	int
		nlocal_histories, 	// number of local histories
		local_history_length, 	// local history length
		theta;			// initial threshold for adaptive theta adjusting

	double 
		fudge;			// fudge factor to multiply by perceptron output
	int 
		*xlat,			// transfer function for 6-bit weights (5-bit magnitudes)
		*xlat4,			// transfer function for 5-bit weights (4-bit magnitudes)
		pcbit,			// bit from the PC to use for hashing global history
		pcshift,		// shift for hashing PC
		block_size,		// number of ghist bits in a "block"; this is the width of an initial hash of ghist
		num_filter;		// number of entries in the filter; 0 means don't use the filter

	bool 
		hash_taken;		// should we hash the taken/not taken value with a PC bit?
	int 
		n_sign_bits,		// number of sign bits per magnitude
		extra_rounds;		// number of extra rounds of training a single weight on a low-confidence prediction

	long long int
		table_sizes[MAX_TABLES];// size of each table

	unsigned int 
		record_mask;		// which histories are updated with filtered branch outcomes

	int 
		speed,			// speed for adaptive theta training
		threshold,		// threshold for deciding low/high confidence
		tunebits;		// number of bits in misprediction counters

	bool 
		tuneonly;		// if true, only count mispredictions of low-confidence branches
	int 
		nbest, 			// use this many of the top performing tables on a low-confidence branch
		bias0, 			// bias perceptron output this much on all-bits-zero local history
		bias1, 			// bias perceptron output this much on all-bits-one local history
		biasmostly0, 		// bias perceptron output this much on almost-all-bits-zero local history
		biasmostly1, 		// bias perceptron output this much on almost-all-bits-one local history
		hshift;			// how much to shift initial feauture hash before XORing with PC bits

	int
		decay;			// whether and how often to decay a random weight

	unsigned long long int 
		imli_mask1, 		// which tables should have their indices hashed with the first IMLI counter
		imli_mask4, 		// which tables should have their indices hashed with the fourth IMLI counter
		recencypos_mask;	// which tables should have their indices hashed with the recency position

	// constructor; parameters described above 

	multiperspective_perceptron (
		history_spec *_specs, int _num_tables = 38, long long int _budgetbits = 65536 * 8 + 2048, int _nlocal_histories = 512, 
		int _local_history_length = 11, int _theta = 10, double _fudge = 0.26, int *_xlat = NULL, int *_xlat4 = NULL, 
		int _pcbit = 2, int _pcshift = -10, int _block_size = 21, int _num_filter = 18000, bool _hash_taken = false,
		int _n_sign_bits = 2, int _extra_rounds = 1, unsigned int _record_mask= 191, int _speed = 9, int _threshold = 0,
		int _tunebits = 30, bool _tuneonly = false, int _nbest = 1, int _bias0 = 0, int _bias1 = 0, int _biasmostly0 = 0, 
		int _biasmostly1 = 0, int _hshift = 0, int _decay = 0, unsigned long long int _imli_mask1 = 0, 
		unsigned long long int _imli_mask4 = 0, unsigned long long int _recencypos_mask = 0) :
			specs(_specs), num_tables(_num_tables), budgetbits(_budgetbits), nlocal_histories(_nlocal_histories),
			local_history_length(_local_history_length), theta(_theta), fudge(_fudge), xlat(_xlat), xlat4(_xlat4), pcbit(_pcbit), pcshift(_pcshift),
			block_size(_block_size), num_filter(_num_filter), hash_taken(_hash_taken), n_sign_bits(_n_sign_bits), extra_rounds(_extra_rounds),
			record_mask(_record_mask), speed(_speed), threshold( _threshold), tunebits(_tunebits), tuneonly(_tuneonly), nbest(_nbest), bias0(_bias0),
			bias1(_bias1), biasmostly0(_biasmostly0), biasmostly1(_biasmostly1), hshift(_hshift), decay(_decay), imli_mask1(_imli_mask1),
			imli_mask4(_imli_mask4), recencypos_mask (_recencypos_mask) {

		// initially nothing in the filter

		occupancy = 0;

		// copy sizes of tables from feature specifications into a convenient array

		for (int i=0; i<num_tables; i++) table_sizes[i] = specs[i].size;

		// no bits used so far

		totalbits = 0;

		// seed the random number generator

		my_seed = 0xdeadbeef;

		// initialize various structures to 0

		imli_counter1 = 0;
		imli_counter2 = 0;
		imli_counter3 = 0;
		imli_counter4 = 0;
		memset (modhist_lengths, 0, sizeof (modhist_lengths));
		memset (modpath_lengths, 0, sizeof (modpath_lengths));

		// analyze the specification to figure out how many bits are allocated to what

		analyze_spec ();

		// set more structures to 0

		memset (ghist_words, 0, sizeof (ghist_words));
		memset (tables, 0, sizeof (tables));
		memset (sign_bits, 0, sizeof (sign_bits));
		memset (acyclic_histories, 0, sizeof (acyclic_histories));
		memset (acyclic2_histories, 0, sizeof (acyclic2_histories));
		memset (local_histories, 0, sizeof (local_histories));
		memset (mod_histories, 0, sizeof (mod_histories));
		memset (modpath_histories, 0, sizeof (modpath_histories));
		memset (recency_stack, 0, sizeof (recency_stack));
		memset (path_history, 0, sizeof (path_history));
		memset (blurrypath_histories, 0, sizeof (blurrypath_histories));

		// threshold counter for adaptive theta training

		tc = 0;
		memset (mpreds, 0, sizeof (mpreds));

		// initialize tables and signs

		for (int i=0; i<num_tables; i++) {
			for (int j=0; j<table_sizes[i]; j++) {
				for (int k=0; k<n_sign_bits; k++) 
					sign_bits[i][j][k] = (i & 1) | (k & 1);
				tables[i][j] = 0;
			}
		}
	}

	// insert a (shifted) PC into the recency stack with LRU replacement

	void insert_recency (unsigned int pc) {
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

	// hash a PC

	unsigned int hash_pc (unsigned int pc) {
		if (pcshift < 0) {
			return hash (pc, -pcshift);
		} else if (pcshift < 11) {
			unsigned int x = pc;
			x ^= (pc >> pcshift);
			return x;
		} else {
			return pc >> (pcshift-11);
		}
	}

	// hash global history from position a to position b

	unsigned int hash_ghist (int a, int b) {
		unsigned int x = 0;
		// am is the next multiple of block_size after a
		int am = (((a/block_size)*block_size)+block_size);
		// bm is the previous multiple of block_size before b
		int bm = (b/block_size)*block_size;

		// the 0th bit of ghist_words[a/block_size] is the most recent bit.
		// so the number of bits between a and am is the number to shift right?

		// start out x as remainder bits from the beginning:
		// x = [ . . . . . b b b b b ]
		x += ghist_words[a/block_size] >> (a-am);
		// add in bits from the middle
		for (int i=am; i<bm; i+=block_size)
			x += ghist_words[i/block_size];
		// add in remainder bits from end:
		// x += [ b b b b b . . . . . ]
		unsigned int y = ghist_words[bm/block_size] & ((1<<(b-bm))-1);
		x += y << (block_size-(b-bm));
		return x;
	}

	// hash the items in the recency stack to a given depth, shifting
	// by a certain amount each iteration, with two different ways of
	// mixing the bits

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

	// hash the "blurry" path history of a given scale to a given
	// depth and given shifting parameter

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

	// hahs path history of a given depth, with a given shift, and in
	// two different mixing styles

	unsigned int hash_path (int depth, int shift, int style) {
		if (style == -1) {
			unsigned int x = 0;
			for (int i=0; i<depth; i++) {
				x <<= shift;
				x += path_history[i];
			}
			return x;
		} else {
			unsigned int x = 0;
			int bm = (depth/block_size)*block_size;
			for (int i=0; i<bm; i+=block_size) {
				for (int j=0; j<block_size; j++)
					x ^= (!!(path_history[i+j] & (1<<shift))) << j;
			}
			int k = 0;
			for (int i=bm; i<depth; i++)
				x ^= (!!(path_history[i] & (1<<shift))) << k++;
			return x;
		}
	}

	// hash ghist together with path history, with a given shift
	// and with two different mixing styles

	unsigned int hash_ghistpath (int depth, int shift, int style) {
		if (style == -1) {
			unsigned int x = 0;
			int bm = (depth/block_size)*block_size;
			unsigned int w;
			for (int i=0; i<bm; i+=block_size) {
				w = ghist_words[i/block_size];
				for (int j=0; j<block_size; j++) {
					x <<= shift;
					x += (path_history[i+j] << 1) | (w & 1);
					w >>= 1;
				}
			}
			w = ghist_words[bm/block_size];
			for (int i=bm; i<depth; i++) {
				x <<= shift;
				x += (path_history[i] << 1) | (w & 1);
				w >>= 1;
			}
			return x;
		} else {
			unsigned int x = 0;
			int bm = (depth/block_size)*block_size;
			unsigned int w = 0;
			for (int i=0; i<bm; i+=block_size) {
				w = ghist_words[i/block_size];
				for (int j=0; j<block_size; j++) {
					x ^= (!!(path_history[i+j] & (1<<shift))) << j;
					x ^= (w & 1) << j;
					w >>= 1;
				}
			}
			w = ghist_words[bm/block_size];
			int k = 0;
			for (int i=bm; i<depth; i++) {
				x ^= (!!(path_history[i] & (1<<shift))) << k;
				x ^= (w & 1) << k;
				w >>= 1;
				k++;
			}
			return x;
		}
	}

	// an alternate formulation of ghistpath: using ranges of path
	// and ghist

	unsigned int hash_sghistpath (int a, int b, int shift) {
		unsigned int x = 0;
		int bm = (b/block_size)*block_size;
		unsigned int w;
		for (int i=a; i<bm; i+=block_size) {
			w = ghist_words[i/block_size];
			for (int j=0; j<block_size; j++) {
				x <<= shift;
				x += (path_history[i+j] << 1) | (w & 1);
				w >>= 1;
			}
		}
		w = ghist_words[bm/block_size];
		for (int i=bm; i<b; i++) {
			x <<= shift;
			x += (path_history[i] << 1) | (w & 1);
			w >>= 1;
		}
		return x;
	}

	// hash acyclic histories with a given modulus, shift, and style

	unsigned int hash_acyclic (int a, int shift, int style) {
		unsigned int x = 0;
		if (style == -1) {
			unsigned int k = 0;
			for (int i=0; i<a+2; i++) {
				x ^= acyclic_histories[a][i] << k;
				k++;
				k %= block_size;
			}
		} else {
			for (int i=0; i<a+2; i++) {
				x <<= shift;
				x += acyclic2_histories[a][i];
			}
		}
		return x;
	}

	// hash modulo history with a given modulus, length, and style

	unsigned int hash_modhist (int a, int b, int n) {
		unsigned int x = 0, k = 0;
		for (int i=0; i<b; i++) {
			x ^= mod_histories[a][i] << k;
			k++;
			k %= n;
		}
		return x;
	}

	// hash modulo path history with a given modulus, depth, and shift

	unsigned int hash_modpath (int a, int depth, int shift) {
		unsigned int x = 0;
		for (int i=0; i<depth; i++) {
			x <<= shift;
			x += modpath_histories[a][i];
		}
		return x;
	}

	// hash modulo path history together with modulo (outcome) history

	unsigned int hash_ghistmodpath (int a, int depth, int shift) {
		unsigned int x = 0;
		for (int i=0; i<depth; i++) {
			x <<= shift;
			x += (modpath_histories[a][i] << 1) | mod_histories[a][i];
		}
		return x;
	}

	// hash the recency position where we find this PC

	unsigned int hash_recencypos (unsigned int pc, int l, int t) {
	
		// search for the PC

		for (int i=0; i<l; i++) {
			if (recency_stack[i] == pc) return i * table_sizes[t] / l;
		}

		// return last index in table on a miss

		return table_sizes[t] - 1;
	}

	// use a history specification to call the corresponding history hash function

	unsigned int get_hash (history_spec *s, unsigned int pc, unsigned int pc2, int t) {
		unsigned int x = 0;
		switch (s->type) {
		case SGHISTPATH:
			x = hash_sghistpath (s->p1, s->p2, s->p3);
			break;
		case GHIST:
			x = hash_ghist (s->p1, s->p2);
			break;
		case GHISTPATH:
			x = hash_ghistpath (s->p1, s->p2, s->p3);
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
		case BIAS:
			x = 0;
			break;
		case RECENCY:
			x = hash_recency (s->p1, s->p2, s->p3);
			break;
		case IMLI:
			if (s->p1 == 1) x = imli_counter1;
			else if (s->p1 == 2) x = imli_counter2;
			else if (s->p1 == 3) x = imli_counter3;
			else if (s->p1 == 4) x = imli_counter4;
			else assert (0);
			break;
		case PATH:
			x = hash_path (s->p1, s->p2, s->p3);
			break;
		case LOCAL:
			x = local_histories[hashlocal() % nlocal_histories];
			if (s->p1 != -1) x &= ((1<<s->p1)-1);
			break;
		case BLURRYPATH:
			x = hash_blurry (s->p1, s->p2, s->p3);
			break;
		case RECENCYPOS:
			x = hash_recencypos (pc2, s->p1, t);
			break;
		default: assert (0);
		}
		return x;
	}

	// for quicksort

	typedef int (*compar_t)(const void *, const void *);

	// find the best performing subset of features (estimated by
	// sorting by mispredictionr rate)

	void findbest (int *bestpreds) {
		if (threshold < 0) return;
		bestpair pairs[num_tables];
		for (int i=0; i<num_tables; i++) {
			pairs[i].index = i;
			pairs[i].mpreds = mpreds[i];
		}
		qsort (pairs, num_tables, sizeof (bestpair), (compar_t) paircmp);
		for (int i=0; i<nbest; i++) bestpreds[i] = pairs[i].index;
	}

	// compute the perceptron output as u.yout

	void compute_output (int *bestval) {

		// list of best predictors

		int bestpreds[MAX_TABLES];

		// initialize sum

		u.yout = 0;

		// bias the prediction by whether the local history is
		// one of four distinctive patterns

		int x = local_histories[hashlocal() % nlocal_histories];
		if (x == 0) 
			u.yout = bias0;
		else if (x == ((1<<local_history_length)-1)) 
			u.yout = bias1;
		else if (x == (1 << (local_history_length-1)))
			u.yout = biasmostly0;
		else if (x == ((1<<(local_history_length-1))-1))
			u.yout = biasmostly1;

		// finds the best subset of features to use in case of a low-confidence branch

		findbest (bestpreds);

		// begin computation of the sum for low-confidence branch

		*bestval = 0;

		// for each feature...

		for (int i=0; i<num_tables; i++) {
			// get the hash to index the table

			unsigned int g = get_hash (&specs[i], u.pc, u.pc2, i);
			unsigned long long int h;

			// shift the hash from the feature to xor with the hashed PC

			if (hshift < 0) {
				h = g;
				h <<= -hshift;
				h ^= u.pc2;
			} else {
				h = g;
				h <<= hshift;
				h ^= u.hpc;
			}

			// xor in the imli counter(s) and/or recency position based on the masks

			if ((1ull<<i) & imli_mask1) h ^= imli_counter1;
			if ((1ull<<i) & imli_mask4) h ^= imli_counter4;
			if ((1ull<<i) & recencypos_mask) h ^= hash_recencypos (u.pc2, 31, i);
			h %= table_sizes[i];
			u.indices[i] = h;

			// add the weight; first get the weight's magnitude

			int c = tables[i][h];

			// get the sign

			bool sign = sign_bits[i][h][u.hpc % n_sign_bits];

			// apply the transfer function and multiply by a coefficient

			int w = specs[i].coeff * ((specs[i].width == 5) ? xlat4[c] : xlat[c]);

			// apply the sign

			int val = sign ? -w : w;

			// add the value

			u.yout += val;

			// if this is one of those good features, add the value to bestval
	
			if (threshold >= 0) for (int j=0; j<nbest; j++) if (bestpreds[j] == i) {
				*bestval += val;
				break;
			}
		}

		// apply a fudge factor to affect when training is triggered

		u.yout *= fudge;
	}

	// hash functions by Thomas Wang, best live link is http://burtleburtle.net/bob/hash/integer.html

        unsigned int hash1 (unsigned int a) {
		a = (a ^ 0xdeadbeef) + (a<<4);
		a = a ^ (a>>10);
		a = a + (a<<7);
		a = a ^ (a>>13);
		return a;
        }

        unsigned int hash2 (unsigned int key) {
		int c2=0x27d4eb2d; // a prime or an odd constant
		key = (key ^ 61) ^ (key >> 16);
		key = key + (key << 3);
		key = key ^ (key >> 4);
		key = key * c2;
		key = key ^ (key >> 15);
		return key;
        }

        // hash a key with the i'th hash function using the common Bloom filter trick
        // of linearly combining two hash functions with i as the slope

        unsigned int hash (unsigned int key, unsigned int i) {
		return hash2 (key) * i + hash1 (key);
        }

	// hash for indexing the table of local histories

	unsigned int hashlocal (void) {
		return u.pc >> 2;
	}

	// hash for indexing the filter

	unsigned int hashf (void) {
		unsigned int x;
		x = last_ghist_bit ^ u.hpc;
		return x;
	}

	// make a prediction

	branch_info *lookup (unsigned int pc) {

		// get different PC hashes

		u.pc = pc;
		u.pc2 = pc >> 2;
		u.hpc = hash_pc (pc);

		// initially assume we have not filtered this branch

		u.filtered = false;
		bool use_static = false;

		// if we have a filter...

		if (num_filter) {

			// find the right filter entry for this branch

			unsigned int x = hashf ();
			filter_entry *f = &filter_table[x % num_filter];

			// always not taken so far?

			bool ansf = f->seen_untaken & !(f->seen_taken);

			// always taken so far?

			bool atsf = f->seen_taken & !(f->seen_untaken);

			// if not taken yet, predict it will not be taken and be done

			if (ansf) {
				u.prediction (false);
				u.filtered = true;
				return &u;
			} else if (atsf) {

			// if not untaken yet, predict it will be taken and be done

				u.prediction (true);
				u.filtered = true;
				return &u;
			}

			// if this is the first time we have encountered
			// this branch, use the static prediction and train

			if (!(f->seen_taken) && !(f->seen_untaken)) use_static = true;
				
		}

		// compute yout and bestval

		int bestval;
		compute_output (&bestval);

		// use the static prediction?
		if (use_static) {
			u.prediction (false); // static prediction of not taken
		} else {

			// if we have a low confidence branch, use the
			// prediction given by the best subset

			if (abs (u.yout) <= threshold) {
				u.prediction (bestval >= 1);
			}
			else

			// otherwise use the full perceptron prediction

				u.prediction (u.yout >= 1);
		}
		return &u;
	}


	// adaptive theta training, adapted from O-GEHL

	void theta_setting (bool correct, int a) {
		if (!correct) {
			tc++;
			if (tc >= speed) {
				theta++;
				tc = 0;
			}
		}
		if (correct && a < theta) {
			tc--;
			if (tc <= -speed) {
				theta--;
				tc = 0;
			}
		}
	}

	// train the perceptron predictor

	void train (bool taken) {

		// was the prediction correct?

		bool correct = (u.yout >= 1) == taken;

		// what is the magnitude of yout?

		int a = abs (u.yout);

		// keep track of mispredictions per table

		if (threshold >= 0) if (!tuneonly || (a <= threshold)) {
			bool halve = false;

			// for each table, figure out if there was a misprediction

			for (int i=0; i<num_tables; i++) {
				bool sign = sign_bits[i][u.indices[i]][u.hpc % n_sign_bits];
				int c = tables[i][u.indices[i]];
				int w = specs[i].coeff * ((specs[i].width == 5) ? xlat4[c] : xlat[c]);
				if (sign) w = -w;
				bool pred = w >= 1;
				if (pred != taken) {
					mpreds[i]++;
					if (mpreds[i] == (1<<tunebits)-1) halve = true;
				}
			}

			// if we reach the maximum counter value, halve all the counters

			if (halve)
				for (int i=0; i<num_tables; i++) mpreds[i] /= 2;
		}

		// if the branch was predicted incorrectly or the correct
		// prediction was weak, update the weights

		bool do_train = !correct || (a <= theta);
		if (!do_train) return;

		// adaptive theta tuning

		theta_setting (correct, a);

		// train the weights, computing what the value of yout
		// would have been if these updates had been applied before

		int newyout = 0;
		for (int i=0; i<num_tables; i++) {

			// get the magnitude

			int c = tables[i][u.indices[i]];

			// get the sign

			bool sign = sign_bits[i][u.indices[i]][u.hpc % n_sign_bits];

			// increment/decrement if taken/not taken

			satincdec (taken, &sign, &c, i);

			// update the magnitude and sign

			tables[i][u.indices[i]] = c;
			sign_bits[i][u.indices[i]][u.hpc % n_sign_bits] = sign;

			// update the new version of yout

			if (sign) newyout -= xlat[c]; else newyout += xlat[c];
		}

		// if the prediction still would have been incorrect even
		// with the updated weights, update some more weights to
		// try to fix the problem

		if ((newyout >= 1) != taken) {
			if (extra_rounds != -1) {
				int z = 0;
over:
				// udpate a random weight

				int besti = -1;
				int r = rand_r (&my_seed) % num_tables;
				int pout;
				for (int j=0; j<num_tables; j++) {
					int i = (r + j) % num_tables;
					int c = tables[i][u.indices[i]];
					bool sign = sign_bits[i][u.indices[i]][u.hpc % n_sign_bits];
					int d = sign ? -xlat[c] : xlat[c];
					pout = newyout - d;
					if ((pout >= 1) == taken) {
						// we have found a weight that if we blow it away will help!
						besti = i;
						break;
					}
				}
				if (besti != -1) {
					int i = besti;
					int c = tables[i][u.indices[i]];
					bool sign = sign_bits[i][u.indices[i]][u.hpc % n_sign_bits];
					if (c > 1) {
						c--;
						tables[i][u.indices[i]] = c;
					}
					int d = sign ? -xlat[c] : xlat[c];
					int out = pout + d;
					z++;
					if ((out >= 1) != taken) {
						if (z < extra_rounds) goto over;
					}
				}
			}
		}
	}

	// saturating increment/decrement for sign/magnitude

	void satincdec (bool taken, bool *sign, int *c, int t) {
		int max_weight = (1<<(specs[t].width-1))-1;
		if (taken) {
			// increment sign/magnitude
			if (*sign) {
				// go toward 0 away from negative max weight
				if (*c == 0)
					*sign = false; // moved to positive 0
				else
					(*c)--;
			} else {
				// go toward max weight away from 0
				if (*c < max_weight) (*c)++;
			}
		} else {
			// decrement sign/magnitude
			if (*sign) {
				// go toward negative max weight down from 0
				if (*c < max_weight) (*c)++;
			} else {
				// go toward 0 away from max weight
				if (*c == 0) 
					*sign = true; // negative 0
				else 
					(*c)--;
			}
		}
	}

	// update the predictor

	void update (branch_info *p, unsigned int target, bool taken, int type) {
		unsigned int hpc = u.hpc;
		bool do_train = true;

		// update the filter

		if (num_filter) {
			unsigned int x = hashf ();
			filter_entry *f = &filter_table[x % num_filter];
			bool ansf;
			bool atsf;

			// compute this first, so we don't not train on the
			// first time a branch is seen.

			ansf = f->seen_untaken & !(f->seen_taken);
			atsf = f->seen_taken & !(f->seen_untaken);
			bool transition = false;
			if (taken) {
				if (f->seen_taken == false) if (f->seen_untaken) transition = true;
				f->seen_taken = true;
			} else {
				if (f->seen_untaken == false) if (f->seen_taken) transition = true;
				f->seen_untaken = true;
			}
			if (ansf) do_train = false;
			if (atsf) do_train = false;

			// is this the first time time the branch has gone both ways?

			if (transition) occupancy++;

			// for every new dynamic branch, when there are
			// more than 'decay' number of branches in the
			// filter, blow a random filter entry away

			if (decay && transition && ((occupancy > decay) || (decay == 1))) {
				int r = rand_r (&my_seed) % num_filter;
				if (filter_table[r].seen_taken && filter_table[r].seen_untaken) occupancy--;
				filter_table[r].seen_taken = false;
				filter_table[r].seen_untaken = false;
			}
		}
		if (do_train) train (taken);

		// mask values for whether or not to record a filtered branch into a history register 

#define RECORD_FILTERED_IMLI	1
#define RECORD_FILTERED_GHIST	2
#define RECORD_FILTERED_PATH	4
#define RECORD_FILTERED_ACYCLIC	8
#define RECORD_FILTERED_MOD	16
#define RECORD_FILTERED_BLURRY	32
#define RECORD_FILTERED_LOCAL	64	// should never record a filtered local branch - duh!
#define RECORD_FILTERED_RECENCY	128

		// four different styles of IMLI

		if (!u.filtered || (record_mask & RECORD_FILTERED_IMLI)) {
			if (target < u.pc) {
				if (taken)
					imli_counter1++;
				else
					imli_counter1 = 0;
				if (!taken)
					imli_counter2++;
				else
					imli_counter2 = 0;
			} else {
				if (taken)
					imli_counter3++;
				else
					imli_counter3 = 0;
				if (!taken)
					imli_counter4++;
				else
					imli_counter4 = 0;
			}
		}

		// we can hash the branch outcome with a PC bit. doesn't really help.

		bool hashed_taken = hash_taken ? (taken ^ !!(u.pc & (1<<pcbit))) : taken;

		// record into ghist

		if (!u.filtered || (record_mask & RECORD_FILTERED_GHIST)) {
			bool ab = hashed_taken;
			for (int i=0; i<ghist_length/block_size+1; i++) {
				unsigned int a = ghist_words[i];
				bool ab_new = (a >> (block_size-1)) & 1;
				a <<= 1;
				a |= ab;
				ab = ab_new;
				a &= (1<<block_size)-1;
				ghist_words[i] = a;
			}
		}

		// record into path history

		if (!u.filtered || (record_mask & RECORD_FILTERED_PATH)) {
			memmove (&path_history[1], &path_history[0], sizeof (unsigned short int) * (path_length-1));
			path_history[0] = u.pc2;
		}

		// record into acyclic history

		if (!u.filtered || (record_mask & RECORD_FILTERED_ACYCLIC)) {
			for (int i=0; i<MAX_ACYCLIC; i++) {
				acyclic_histories[i][hpc%(i+2)] = hashed_taken;
				acyclic2_histories[i][hpc%(i+2)] = hpc;
			}
		}

		// record into modulo path history

		if (!u.filtered || (record_mask & RECORD_FILTERED_MOD)) {
			for (int ii=0; ii<nmodpath_histories; ii++) {
				int i = modpath_indices[ii];
				if (hpc % (i+2) == 0) {
					memmove (&modpath_histories[i][1], &modpath_histories[i][0], sizeof (unsigned short int) * (modpath_lengths[ii]-1));
					modpath_histories[i][0] = u.pc2;
				}
			}
		}

		// update blurry history

		if (!u.filtered || (record_mask & RECORD_FILTERED_BLURRY)) {
			for (int i=0; i<MAX_BLURRY; i++) {
				unsigned int z = u.pc >> i;
				if (blurrypath_histories[i][0] != z) {
					memmove (&blurrypath_histories[i][1], &blurrypath_histories[i][0], sizeof (unsigned int) * (MAX_BLURRY2-1));
					blurrypath_histories[i][0] = z;
				}
			}
		}

		// record into modulo pattern history

		if (!u.filtered || (record_mask & RECORD_FILTERED_MOD)) {
			for (int ii=0; ii<nmodhist_histories; ii++) {
				int i = modhist_indices[ii];
				if (hpc % (i+2) == 0) {
					memmove (&mod_histories[i][1], &mod_histories[i][0], modhist_lengths[ii]-1);
					mod_histories[i][0] = hashed_taken;
				}
			}
		}

		// insert this PC into the recency stack

		if (!u.filtered || (record_mask & RECORD_FILTERED_RECENCY)) {
			insert_recency (u.pc2);
		}

		// record into a local history

		if (!u.filtered || (record_mask & RECORD_FILTERED_LOCAL)) {
			unsigned int *lo = &local_histories[hashlocal() % nlocal_histories];
			*lo <<= 1;
			*lo |= hashed_taken;
			*lo &= ((1<<local_history_length)-1);
		}

		// update last ghist bit, used to index filter

		last_ghist_bit = taken;
	}

	// update ghist and path history on branches that aren't conditional

	void nonconditional_branch (unsigned int pc) {
		unsigned short int pc2 = pc >> 2;

		bool ab = !(pc & (1<<pcbit));
		for (int i=0; i<ghist_length/block_size+1; i++) {
			bool ab_new = (ghist_words[i] >> (block_size-1)) & 1;
			ghist_words[i] <<= 1;
			ghist_words[i] |= ab;
			ghist_words[i] &= (1<<block_size)-1;
			ab = ab_new;
		}
		memmove (&path_history[1], &path_history[0], sizeof (unsigned short int) * (path_length-1));
		path_history[0] = pc2;
	}
};

// tuned features for 64KB predictor

history_spec big[37] = {
{ ACYCLIC, 10, -1, -1, 1.0, 0, 6 },
{ BLURRYPATH, 10, 7, -1, 1.0, 0, 6 },
{ GHIST, 0, 19, 1, 1.3125, 0, 6 },
{ GHIST, 0, 65, 1, 0.850000, 0, 6 },
{ GHIST, 111, 222, 1, 1.0, 0, 6 },
{ GHIST, 115, 206, 1, 1.0, 0, 6 },
{ GHIST, 190, 336, 1, 1.0, 0, 6 },
{ GHIST, 21, 64, 1, 1.0, 0, 6 },
{ GHIST, 48, 119, 1, 1.0, 0, 6 },
{ GHIST, 67, 203, 1, 0.75, 0, 6 },
{ GHIST, 75, 150, 1, 1.0, 0, 6 },
{ GHISTMODPATH, 1, 5, 4, 1.450000, 0, 6 },
{ GHISTMODPATH, 3, 13, 1, 1.0, 0, 6 },
{ GHISTPATH, 105, 4, 0, 1.0, 0, 6 },
{ GHISTPATH, 11, 2, -1, 1.25, 0, 6 },
{ GHISTPATH, 15, 4, -1, 1.125, 0, 6 },
{ GHISTPATH, 23, 4, -1, 1.0, 0, 6 },
{ GHISTPATH, 31, 1, -1, 1.0, 0, 6 },
{ GHISTPATH, 32, 2, -1, 1.0, 0, 6 },
{ GHISTPATH, 36, 4, -1, 1.0, 0, 6 },
{ GHISTPATH, 51, 1, 0, 1.0, 0, 6 },
{ GHISTPATH, 7, 1, -1, 1.5, 0, 6 },
{ GHISTPATH, 72, 1, -1, 1.0, 0, 6 },
{ GHISTPATH, 86, 4, 0, 1.0, 0, 6 },
{ IMLI, 1, -1, -1, 1.8125, 0, 6 },
{ IMLI, 4, -1, -1, 1.78125, 0, 6 },
{ LOCAL, -1, -1, -1, 1.0, 0, 6 },
{ LOCAL, -1, -1, -1, 2.0, 0, 6 },
{ MODHIST, 1, 16, -1, 0.9375, 0, 6 },
{ MODPATH, 3, 20, 1, 1.0, 0, 6 },
{ PATH, 13, 2, -1, 1.4375, 0, 6 },
{ PATH, 27, 5, -1, 1.00000, 0, 6 },
{ RECENCY, 10, 3, -1, 0.550000, 0, 6 },
{ RECENCY, 14, 4, -1, 1.3125, 0, 6 },
{ RECENCYPOS, 31, -1, -1, 1.5, 0, 6 },
{ SGHISTPATH, 1, 2, 5, 1.25, 0, 6 },
{ SGHISTPATH, 1, 5, 2, 1.3125, 0, 6 },
};

PREDICTOR::PREDICTOR(void){
	static int branchyxlat[] = { 1,3,4,5,7,8,9,11,12,14,15,17,19,21,23,25,27,29,32,34,37,41,45,49,53,58,63,69,76,85,94,106,};
	static int xlat4[] = { 0,4,5,7,9,11,12,14,16,17,19,22,28,33,39,45, };
	ex = new multiperspective_perceptron (
		big,
		37,
		450000000,	// large budget
		262144,         // number of local histories
		11,             // local history length
		10,             // initial theta
		0.245,          // fudge factor 
		branchyxlat,    // xlat vector
		xlat4,          // xlat vector
		2,              // pcbit
		-10,            // pcshift
		21,             // block size
		262144,         // number of filter entries
		false,          // hash taken
		2,              // number of sign bits
		1,              // number of rounds correcting after mispred training
		191,            // mask for recording filtered branches
		9,              // adaptive theta learning speed
		1,              // threshold
		24,             // tunebits
		1,              // tune only
		20,             // nbest
		-5,             // bias for lhist=all bits zero
		5,              // bias for lhist=all bits one
		-1,             // bias for lhist=most recent bits 0
		1,              // bias for lhist=most recent bits 1
		-6,             // hash shift
		0, 
		0xc1000, 	// imli count 1 mask
		0x80008000, 	// imli count 4 mask
		0x100000090);	// recency position mask
}

bool PREDICTOR::GetPrediction(UINT64 PC){
	u = ex->lookup (PC);
	return u->prediction ();
}

void PREDICTOR::UpdatePredictor(UINT64 PC, OpType opType, bool resolveDir, bool predDir, UINT64 branchTarget){
	ex->update (u, branchTarget, resolveDir, opType);
}

void PREDICTOR::TrackOtherInst(UINT64 PC, OpType opType, bool branchDir, UINT64 branchTarget){
	ex->nonconditional_branch (PC);
}

#endif
