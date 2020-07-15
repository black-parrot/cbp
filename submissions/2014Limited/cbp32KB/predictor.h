///////////////////////////////////////////////////////////////////////
////  Copyright 2015 Samsung Austin Semiconductor, LLC.                //
/////////////////////////////////////////////////////////////////////////
//


//#define REALISTIC // uncomment to get a realistic predictor within the 256 Kbits limit , with only 12 1024 entries tagged tables in the TAGE predictor, and a global history and single local history GEHL statistical corrector:  total misprediction numbers 2.430 MPKI



#ifndef _PREDICTOR_H_
#define _PREDICTOR_H_

//To get the predictor storage budget on stderr  uncomment the next line
//#define PRINTSIZE

#include <cstdlib>
#include <cstring>
#include <cassert>
#include <inttypes.h>
#include <cmath>
#include "utils.h"
#include "bt9.h"
#include "bt9_reader.h"

#define LOOPPREDICTOR	//loop predictor enable
#define NNN 1			// number of entries allocated on a TAGE misprediction

// number of tagged tables
#ifndef REALISTIC
#define NHIST 15
#else 
#define NHIST 12
#endif

#define HYSTSHIFT 2 // bimodal hysteresis shared by 4 entries
#define LOGB 14 // log of number of entries in bimodal predictor 



#define PERCWIDTH 6 //Statistical coorector maximum counter width 

//The statistical corrector components

//global branch GEHL
#ifdef REALISTIC
#define LOGGNB 10
#else
#define LOGGNB 9
#endif
#define GNB 4
int Gm[GNB] = {16,11, 6,3};
int8_t GGEHLA[GNB][(1 << LOGGNB)];
int8_t *GGEHL[GNB];

//large local history
#ifdef REALISTIC
#define LOGLNB 10
#define LNB 4
int Lm[LNB] = {16,11,6,3};
int8_t LGEHLA[LNB][(1 << LOGLNB)];
int8_t *LGEHL[LNB];
#else
#define LOGLNB 10
#define LNB 3
int Lm[LNB] = {11,6,3};
int8_t LGEHLA[LNB][(1 << LOGLNB)];
int8_t *LGEHL[LNB];
#endif
#define  LOGLOCAL 8
#define NLOCAL (1<<LOGLOCAL)
#define INDLOCAL (PC & (NLOCAL-1))
long long L_shist[NLOCAL];


// small local history
#define LOGSNB 9
#define SNB 4
int Sm[SNB] =  {16,11, 6, 3};
int8_t SGEHLA[SNB][(1 << LOGSNB)];
int8_t *SGEHL[SNB];
#define LOGSECLOCAL 4
#define NSECLOCAL (1<<LOGSECLOCAL)	//Number of second local histories
#define INDSLOCAL  (((PC ^ (PC >>5))) & (NSECLOCAL-1))	
long long S_slhist[NSECLOCAL];


#define LOGTNB 9
#define TNB 3
int Tm[TNB] =  {11, 6, 3};
int8_t TGEHLA[TNB][(1 << LOGTNB)];
int8_t *TGEHL[TNB];
#define INDTLOCAL  (((PC ^ (PC >>3))) & (NSECLOCAL-1))	// differen hash for thehistory
long long T_slhist[NSECLOCAL];

//return-stack associated history component
#define PNB 4
#define LOGPNB 9
int Pm[PNB] ={16,11,6,3};
int8_t PGEHLA[PNB][(1 << LOGPNB)];
int8_t *PGEHL[PNB];
long long HSTACK[16];
int pthstack;


//parameters of the loop predictor
#define LOGL 5
#define WIDTHNBITERLOOP 10	// we predict only loops with less than 1K iterations
#define LOOPTAG 10		//tag width in the loop predictor


//update threshold for the statistical corrector
#ifdef REALISTIC
#define LOGSIZEUP 0
#else
#define LOGSIZEUP 5
#endif
int Pupdatethreshold[(1 << LOGSIZEUP)];	//size is fixed by LOGSIZEUP
#define INDUPD (PC & ((1 << LOGSIZEUP) - 1))

// The three counters used to choose between TAGE ang SC on High Conf TAGE/Low Conf SC
int8_t  FirstH, SecondH, ThirdH;
#define CONFWIDTH 7	//for the counters in the choser


#define PHISTWIDTH 27		// width of the path history used in TAGE




#define UWIDTH 2 // u counter width on TAGE		
#define CWIDTH 3		// predictor counter width on the TAGE tagged tables

#define HISTBUFFERLENGTH 4096	// we use a 4K entries history buffer to store the branch history


//the counter(s) to chose between longest match and alternate prediction on TAGE when weak counters
#ifdef REALISTIC
#define LOGSIZEUSEALT 0
#else
#define LOGSIZEUSEALT 8 
#endif
#define SIZEUSEALT  (1<<(LOGSIZEUSEALT))
#define INDUSEALT (PC & (SIZEUSEALT -1))
int8_t use_alt_on_na[SIZEUSEALT][2];



long long GHIST;


//The two BIAS tables in the SC component
#define LOGBIAS 7
int8_t Bias[(1<<(LOGBIAS+1))];
#define INDBIAS (((PC<<1) + pred_inter) & ((1<<(LOGBIAS+1)) -1))
int8_t BiasSK[(1<<(LOGBIAS+1))];
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




int TICK;// for the reset of the u counter


uint8_t ghist[HISTBUFFERLENGTH];
int ptghist;
long long phist;		//path history
folded_history ch_i[NHIST + 1];	//utility for computing TAGE indices
folded_history ch_t[2][NHIST + 1];	//utility for computing TAGE tags

//For the TAGE predictor
bentry *btable;			//bimodal TAGE table
gentry *gtable[NHIST + 1];	// tagged TAGE tables
#ifdef REALISTIC
int m[NHIST+1]={0,8,12,18,27,40,60,90,135,203,305,459,690}; // history lengths
int TB[NHIST + 1]  ={0,8, 9, 9,10,10,11,11,12,12,13,13,14};		// tag width for the different tagged tables
int logg[NHIST + 1]={0,10,10,10,10,10,10,10,10,10,10,10,10};		// log of number entries of the different tagged tables
#else
int m[NHIST+1]={0,6,10,18,25,35,55,69,105,155,230,354,479,642,1012,1347}; // history lengths
int TB[NHIST + 1]={0,7,9,9,9,10,11,11,12,12,12,13,14,15,15,15};		// tag width for the different tagged tables
int logg[NHIST + 1]={0,10,10,10,11,10,10,10,10,10,9,9,9,8,7,7};		// log of number entries of the different tagged tables
#endif 
int GI[NHIST + 1];		// indexes to the different tables are computed only once  
uint GTAG[NHIST + 1];	// tags for the different tables are computed only once  
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
  int inter=0;
  
   for (int i = 1; i <= NHIST; i += 1)
    {
      STORAGESIZE += (1 << (logg[i])) * (CWIDTH + UWIDTH + TB[i]);

    }
 STORAGESIZE += 2 * (SIZEUSEALT) * 4;
STORAGESIZE += (1 << LOGB) + (1 << (LOGB - HYSTSHIFT));
STORAGESIZE+= m[NHIST];
STORAGESIZE += PHISTWIDTH;
STORAGESIZE += 10 ; //the TICK counter

fprintf (stderr, " (TAGE %d) ", STORAGESIZE);  

#ifdef LOOPPREDICTOR

  inter= (1 << LOGL) * (2 * WIDTHNBITERLOOP + LOOPTAG + 4 + 4 + 1);fprintf (stderr, " (LOOP %d) ", inter); 
  STORAGESIZE+= inter;
  
#endif


  inter = 8 * (1 << LOGSIZEUP) ; //the update threshold counters
inter += (PERCWIDTH) * 4 * (1 << (LOGBIAS));
inter += (GNB-2) * (1 << (LOGGNB)) * (PERCWIDTH - 1) + (1 << (LOGGNB-1))*(2*PERCWIDTH-1);
inter += (LNB-2) * (1 << (LOGLNB)) * (PERCWIDTH - 1) + (1 << (LOGLNB-1))*(2*PERCWIDTH-1);

#ifndef REALISTIC
inter += (SNB-2) * (1 << (LOGSNB)) * (PERCWIDTH - 1) + (1 << (LOGSNB-1))*(2*PERCWIDTH-1);
inter += (TNB-2) * (1 << (LOGTNB)) * (PERCWIDTH - 1) + (1 << (LOGTNB-1))*(2*PERCWIDTH-1);
inter += (PNB-2) * (1 << (LOGPNB)) * (PERCWIDTH - 1) + (1 << (LOGPNB-1))*(2*PERCWIDTH-1);



inter += 16*16; // the history stack
inter += 4; // the history stack pointer
     inter += 16; //global histories for SC
  
  inter += NSECLOCAL * (Sm[0]+Tm[0]);
#endif
  inter += NLOCAL * Lm[0];
inter += 3*CONFWIDTH; //the 3 counters in the choser
STORAGESIZE+= inter;
fprintf (stderr, " (SC %d) ", inter);



#ifdef PRINTSIZE
  fprintf (stderr, " (TOTAL %d) ", STORAGESIZE);
#endif


  return (STORAGESIZE);


}






class PREDICTOR
{
public:
  PREDICTOR (void)
  {

    reinit ();
#ifdef PRINTSIZE    
    predictorsize ();
#endif
  }


  void reinit ()
  {

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
#ifndef REALISTIC
for (int i = 0; i < SNB; i++)
      SGEHL[i] = &SGEHLA[i][0];
    for (int i = 0; i < TNB; i++)
      TGEHL[i] = &TGEHLA[i][0];
    for (int i = 0; i < PNB; i++)
      PGEHL[i] = &PGEHLA[i][0];
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
    for (int i = 0; i < PNB; i++)
      for (int j = 0; j < ((1 << LOGPNB) - 1); j++)
	{
	  if (!(j & 1))
	    {
	      PGEHL[i][j] = -1;

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
    return ((PC & ((1 << (LOGL - 2)) - 1)) << 2);
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
	    if (ltable[index].CurrentIter + 1 == ltable[index].NbIter)
	      return (!(ltable[index].dir));
	    else
	      return ((ltable[index].dir));
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
  void Tagepred (UINT64 PC)
  {
    HitBank = 0;
    AltBank = 0;
for (int i = 1; i <= NHIST; i++)
      {
	GI[i] = gindex (PC, i, phist, ch_i);
	GTAG[i] = gtag (PC, i, ch_t[0], ch_t[1]);
      }

    BI = PC & ((1 << LOGB) - 1);


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

	HighConf =
	  (abs (2 * gtable[HitBank][GI[HitBank]].ctr + 1) >=
	   (1 << CWIDTH) - 1);
      }
    else
      {
	alttaken = getbim ();
	tage_pred = alttaken;
	LongestMatchPred = alttaken;
      }



  }
//compute the prediction

  bool GetPrediction (UINT64 PC)
  {
// computes the TAGE table addresses and the partial tags
    
    Tagepred (PC);
    pred_taken = tage_pred;

#ifdef LOOPPREDICTOR
    predloop = getloop (PC);	// loop prediction
    pred_taken = ((WITHLOOP >= 0) && (LVALID)) ? predloop : pred_taken;
#endif

    pred_inter = pred_taken;

//Compute the SC prediction


// begin to bias the sum towards TAGE predicted direction
    LSUM = 1;
    LSUM += 2 * (GNB + PNB);
#ifndef REALISTIC
    LSUM += 2*(SNB+ LNB+TNB);
#endif
    if (!pred_inter)
      LSUM = -LSUM;

//integrate BIAS prediction   
 int8_t ctr = Bias[INDBIAS];
    LSUM += (2 * ctr + 1);
    ctr = BiasSK[INDBIASSK];
    LSUM += (2 * ctr + 1);

//integrate the GEHL predictions
LSUM += Gpredict ((PC<<1)+pred_inter, GHIST, Gm, GGEHL, GNB, LOGGNB);
LSUM += Gpredict (PC, L_shist[INDLOCAL], Lm, LGEHL, LNB, LOGLNB);
#ifndef REALISTIC
LSUM += Gpredict (PC, S_slhist[INDSLOCAL], Sm, SGEHL, SNB, LOGSNB);
LSUM += Gpredict (PC, T_slhist[INDTLOCAL], Tm, TGEHL, TNB, LOGTNB);
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
	    if ((abs (LSUM) <
		 Pupdatethreshold[INDUPD] / 3))
	      pred_taken = (FirstH < 0) ? SCPRED : pred_inter;

	    else
	      if ((abs (LSUM) <
		   2 * Pupdatethreshold[INDUPD] / 3))
	      pred_taken = (SecondH < 0) ? SCPRED : pred_inter;
	    else
	      if ((abs (LSUM) <
		   Pupdatethreshold[INDUPD]))
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
if (brtype == OPTYPE_RET_COND || brtype == OPTYPE_JMP_DIRECT_COND || brtype == OPTYPE_JMP_INDIRECT_COND || brtype == OPTYPE_CALL_DIRECT_COND || brtype == OPTYPE_CALL_INDIRECT_COND)
      maxt = 1;
    else
      maxt = 4;

    // the return stack associated history

PH = (PH << 1) ^ (target ^(target >> 5) ^ taken);

    if (brtype == OPTYPE_RET_COND || brtype == OPTYPE_JMP_DIRECT_COND || brtype == OPTYPE_JMP_INDIRECT_COND || brtype == OPTYPE_CALL_DIRECT_COND || brtype == OPTYPE_CALL_INDIRECT_COND)
      {
            GBRHIST = (GBRHIST << 1) +   taken;       
	LH = (LH << 1) + (taken);
	SH = (SH << 1) + (taken);
	SH ^= (PC & 15);
        TH=  (TH << 1) + (taken);
     }
    if (brtype ==  OPTYPE_RET_COND) {pthstack= (pthstack-1) & 15; 
    }
    
    if (brtype ==  OPTYPE_CALL_DIRECT_COND){
         int index= (pthstack+1) & 15; HSTACK[index]= HSTACK[pthstack];
         pthstack= index;
    }
    

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

  void UpdatePredictor (UINT64 PC, bool branchDir, bool resolveDir, bool predDir,
			UINT64 branchTarget)
  {

    
   
#ifdef LOOPPREDICTOR
    if (LVALID)
      {


	if (pred_taken != predloop)
	  ctrupdate (WITHLOOP, (predloop == resolveDir), 7);

      }

    loopupdate (PC, resolveDir, (pred_taken != resolveDir));


#endif

bool  SCPRED = (LSUM >= 0);    
if (pred_inter != SCPRED)
	{
     	if ((abs (LSUM) <
		     Pupdatethreshold[INDUPD]))        if ((HighConf))
	    {

	      if ((abs (LSUM) <
		   Pupdatethreshold[INDUPD] / 3))
		ctrupdate (FirstH, (pred_inter == resolveDir), CONFWIDTH);
	      else
		if ((abs (LSUM) <
		     2 * Pupdatethreshold[INDUPD] / 3))
		ctrupdate (SecondH, (pred_inter == resolveDir), CONFWIDTH);
	      else
		if ((abs (LSUM) <
		     Pupdatethreshold[INDUPD]))
		ctrupdate (ThirdH, (pred_inter == resolveDir), CONFWIDTH);

	    }
	}

      if ((SCPRED != resolveDir)
	  || ((abs (LSUM) < Pupdatethreshold[INDUPD])))
	{
	  {
	    if (SCPRED != resolveDir)
	      Pupdatethreshold[INDUPD] += 1;
	    else
	      Pupdatethreshold[INDUPD] -= 1;

	    if (Pupdatethreshold[INDUPD] >= 256)
	      Pupdatethreshold[INDUPD] = 255;
	    if (Pupdatethreshold[INDUPD] < 0)
	      Pupdatethreshold[INDUPD] = 0;
	  }

	  ctrupdate (Bias
		     [INDBIAS], resolveDir, PERCWIDTH);
ctrupdate (BiasSK[INDBIASSK], resolveDir, PERCWIDTH); 
	  Gupdate ((PC<<1)+pred_inter, resolveDir, GHIST, Gm, GGEHL, GNB, LOGGNB);
	  Gupdate (PC, resolveDir, L_shist[INDLOCAL], Lm, LGEHL,
		   LNB, LOGLNB);
	  
#ifndef REALISTIC
Gupdate (PC, resolveDir, S_slhist[INDSLOCAL], Sm, SGEHL, SNB, LOGSNB);
Gupdate (PC, resolveDir,T_slhist[INDTLOCAL], Tm, TGEHL, TNB, LOGTNB);
	  Gupdate (PC, resolveDir, HSTACK[pthstack], Pm, PGEHL, PNB, LOGPNB);
#endif

	}


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
                  i +=  1;
		  T -= 1;
		}
	      else
		{
		  Penalty++;
		}
	    }
	  TICK += (Penalty -NA );
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
                     if (gtable[i][j].u)  gtable[i][j].u--;
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
        {if (gtable[HitBank][GI[HitBank]].u < (1 << UWIDTH) - 1)
                  gtable[HitBank][GI[HitBank]].u++;
        }
//END TAGE UPDATE
      
      HistoryUpdate (PC, OPTYPE_JMP_DIRECT_COND, resolveDir, branchTarget,  phist,
		       ptghist, ch_i, ch_t[0],
		       ch_t[1], L_shist[INDLOCAL], 
		       S_slhist[INDSLOCAL], T_slhist[INDTLOCAL], HSTACK[pthstack], GHIST);

//END PREDICTOR UPDATE


  }

  int Gpredict (UINT64 PC, long long BHIST, int *length, int8_t ** tab,
		int NBR, int logs)
  {
//calcul de la somme, commence par  le biais du branchement
    int PERCSUM = 0;
    for (int i = 0; i < NBR; i++)
      {
	long long bhist = BHIST & ((long long) ((1 << length[i]) - 1));

	long long index =
	  (((long long) PC) ^ bhist ^ (bhist >> (8 - i)) ^
	   (bhist >> (16 - 2 * i)) ^ (bhist >> (24 - 3 * i)) ^ (bhist >>
								(32 -
								 3 *
								 i)) ^ (bhist
									>> (40
									    -
									    4
									    *
									    i)))
             & ((1 << (logs-(i >=(NBR-2)))) - 1);
                
	int16_t ctr = tab[i][index];
	PERCSUM += (2 * ctr + 1);
      }
return ((PERCSUM));

  }
  void Gupdate (UINT64 PC, bool taken, long long BHIST, int *length,
		int8_t ** tab, int NBR, int logs)
  {


    for (int i = 0; i < NBR; i++)
      {
	long long bhist = BHIST & ((long long) ((1 << length[i]) - 1));
	long long index =
	  (((long long) PC) ^ bhist ^ (bhist >> (8 - i)) ^
	   (bhist >> (16 - 2 * i)) ^ (bhist >> (24 - 3 * i)) ^ (bhist >>
								(32 -
								 3 *
								 i)) ^ (bhist
									>> (40
									    -
									    4
									    *
									    i)))
  & ((1 << (logs-(i >=(NBR-2)))) - 1);

        ctrupdate (tab[i][index], taken, PERCWIDTH -(i<(NBR-1)));
         
 
      }

	

  }


  void TrackOtherInst (UINT64 PC, OpType opType, bool branchDir, UINT64 branchTarget)
  {

    bool taken = true;
  
    switch (opType)
      {

   case    OPTYPE_RET_UNCOND:
   case    OPTYPE_JMP_DIRECT_UNCOND:
   case    OPTYPE_JMP_INDIRECT_UNCOND:
   case    OPTYPE_CALL_DIRECT_UNCOND:
   case    OPTYPE_CALL_INDIRECT_UNCOND:

	HistoryUpdate (PC, opType, taken, branchTarget, phist,
		       ptghist, ch_i,
		       ch_t[0], ch_t[1],
		       L_shist[INDLOCAL], 
                       S_slhist[INDSLOCAL], T_slhist[INDTLOCAL],HSTACK[pthstack], GHIST);
	break;


      default:;
      }
  

  }

};





#endif
