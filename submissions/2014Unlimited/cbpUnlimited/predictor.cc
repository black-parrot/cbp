///////////////////////////////////////////////////////////////////////
////  Copyright 2015 Samsung Austin Semiconductor, LLC.                //
/////////////////////////////////////////////////////////////////////////
#include "predictor.h"

/////////////// STORAGE BUDGET JUSTIFICATION ////////////////

// unlimited size

/////////////////////////////////////////////////////////////

#define GPSTEP 6
#define LPSTEP 6
#define BPSTEP 6
#define PPSTEP 6
#define SPSTEP 6
#define YPSTEP 6
#define TPSTEP 6

#define GWIDTH 60
#define LWIDTH 60
#define BWIDTH 42
#define PWIDTH 60
#define SWIDTH 60
#define YWIDTH 60
#define TWIDTH 60

#define LOGTAB 19
#define TABSIZE (1<<LOGTAB)
#define LOGB LOGTAB
#define LOGG LOGTAB
#define LOGSIZE (10)
#define LOGSIZEG (LOGSIZE)
#define LOGSIZEL (LOGSIZE)
#define LOGSIZEB (LOGSIZE)
#define LOGSIZES (LOGSIZE)
#define LOGSIZEP (LOGSIZE)
#define LOGSIZEY (LOGSIZE)
#define LOGSIZET (LOGSIZE)  


// The perceptron-inspired components
int8_t PERC[(1 << LOGSIZEP)][10 * (1 << GPSTEP)] = {{0}};
int8_t PERCLOC[(1 << LOGSIZEL)][10 * (1 << LPSTEP)] = {{0}};
int8_t PERCBACK[(1 << LOGSIZEB)][10 * (1 << BPSTEP)] = {{0}};
int8_t PERCYHA[(1 << LOGSIZEY)][10 * (1 << YPSTEP)] = {{0}};
int8_t PERCPATH[(1 << LOGSIZEP)][10 * (1 << PPSTEP)] = {{0}};
int8_t PERCSLOC[(1 << LOGSIZES)][10 * (1 << SPSTEP)] = {{0}};
int8_t PERCTLOC[(1 << LOGSIZET)][10 * (1 << TPSTEP)] = {{0}};

// the three 16 local histories components

#define LOGLOCAL 4
#define NLOCAL (1<<LOGLOCAL)
#define INDLOCAL (PC & (NLOCAL-1))
long long L_shist[NLOCAL] = {0};
#define LNB 15
int Lm[LNB] = { 2,4,6,9,12,16,20,24,29,34,39,44,50,56,63};
int8_t LGEHLA[LNB][TABSIZE] = {{0}};
int8_t *LGEHL[LNB] = {0};

#define  LOGSECLOCAL 4
#define NSECLOCAL (1<<LOGSECLOCAL)	
#define NB 3
#define INDSLOCAL  (((PC ^ (PC >>5)) >> NB) & (NSECLOCAL-1))
long long S_slhist[NSECLOCAL] = {0};
#define SNB 15
int Sm[SNB] = { 2,4,6,9,12,16,20,24,29,34,39,44,50,56,63};
int8_t SGEHLA[SNB][TABSIZE] = {{0}};
int8_t *SGEHL[SNB] = {0};

#define  LOGTLOCAL 4
#define NTLOCAL (1<<LOGTLOCAL)
#define INDTLOCAL  (((PC ^ (PC >>3)  ^(PC >> 6))) & (NTLOCAL-1)) 
long long T_slhist[NTLOCAL] = {0};
#define TNB 15
int Tm[TNB] = { 2,4,6,9,12,16,20,24,29,34,39,44,50,56,63};
int8_t TGEHLA[TNB][TABSIZE] = {{0}};
int8_t *TGEHL[TNB] = {0};

// the 32768 local histories component

#define QPSTEP 6
#define QWIDTH 60
#define LOGSIZEQ (LOGSIZE)
int8_t PERCQLOC[(1 << LOGSIZEQ)][10 * (1 << QPSTEP)] = {{0}};
#define  LOGQLOCAL 15
#define NQLOCAL (1<<LOGQLOCAL)	//Number of second local histories
#define INDQLOCAL  (((PC ^ (PC >>2) ^(PC>> 4) ^ (PC >> 8))) & (NQLOCAL-1))	
long long Q_slhist[NQLOCAL] = {0};
#define QNB 15
int Qm[QNB] = { 2,4,6,9,12,16,20,24,29,34,39,44,50,56,63};
int8_t QGEHLA[QNB][TABSIZE] = {{0}};
int8_t *QGEHL[QNB] = {0};


//about the skeleton histories
#define YNB 15
int Ym[SNB] = { 2,4,6,9,12,16,20,24,29,34,39,44,50,56,63};

int8_t YGEHLA[SNB][TABSIZE] = {{0}};
int8_t *YGEHL[SNB] = {0};
long long YHA = 0;
int LastBR[8] = {0};

#define BNB 10
int Bm[BNB] = { 2,4,6,9,12,16,20,24,29,34};

int8_t BGEHLA[BNB][TABSIZE] = {{0}};
int8_t *BGEHL[BNB] = {0};


//for the choser
#define CGNB 8
int CGm[CGNB] = { 2,4,6,9,12,16,21,26 };

int8_t CGGEHLA[CGNB][TABSIZE] = {{0}};
int8_t *CGGEHL[CGNB] = {0};
#define CLNB 8
int CLm[CLNB] = { 2,4,6,9,12,16,21,26 };

int8_t CLGEHLA[CLNB][TABSIZE] = {{0}};
int8_t *CLGEHL[CLNB] = {0};
int indexchose;
int LCHOSE;
int Cupdatethreshold;
#define LOGSIZECHOSE LOGTAB
#define MAXSIZECHOSE  TABSIZE
#define SIZECHOSE  (1<<(LOGSIZECHOSE))
#define INDCHOSE (PC & (SIZECHOSE -1))
int8_t GCHOSE[MAXSIZECHOSE * 4] = {0};
int8_t GCHOSECOLT[TABSIZE] = {0};
int indexchosecolt;
// end choser


#define LOGSIZEUSEALT 15
#define MAXSIZEUSEALT  (1<<15)
int8_t use_alt_on_na[MAXSIZEUSEALT][2] = {{0}};
#define SIZEUSEALT  (1<<(LOGSIZEUSEALT))
#define INDUSEALT (PC & (SIZEUSEALT -1))



int PERCSUM;
int Pupdatethreshold[(1 << LOGSIZE)] = {0};
#define INDUPD (PC & ((1 << LOGSIZE) - 1))
int updatethreshold;


//the GEHL predictor 
#define MAXNHISTGEHL 256 
#ifndef LOGGEHL
// base 2 logarithm of number of entries  on each GEHL  component
#define LOGGEHL (LOGTAB+1)
#endif
#define MINSTEP 2
#define MINHISTGEHL 1
static int8_t GEHL[1 << LOGGEHL][MAXNHISTGEHL + 1] = {{0}};	//GEHL tables
int mgehl[MAXNHISTGEHL + 1] = {0};	//GEHL history lengths
int GEHLINDEX[MAXNHISTGEHL + 1] = {0};


//The MACRHSP inspired predictor
#define MAXNRHSP 128
#define LOGRHSP LOGGEHL
static int8_t RHSP[1 << LOGRHSP][MAXNRHSP + 1] = {{0}};	//RHSP tables
int mrhsp[MAXNRHSP + 1] = {0};	//RHSP history lengths
int RHSPINDEX[MAXNRHSP + 1] = {0};
int SUMRHSP;


#define PERCWIDTH 8	
#define CHOSEWIDTH 6	

// local history management

long long BHIST = 0;
long long lastaddr = 0;
long long P_phist = 0;
long long GHIST = 0;


int SUMGEHL;


#define LOGBIAS (LOGB-4)
int8_t Bias[(1 << (LOGBIAS + 1))] = {0};
#define INDBIAS (((PC<<1) ^ PRED) & ((1<<(LOGBIAS+1))-1))
#define LOGBIASCOLT (LOGTAB)
int8_t BiasColt[TABSIZE] = {0};
#define INDBIASCOLT (((PC<<5) ^ PRED ^ ((predtaken[0] ^ (predtaken[1] <<1) ^ (predtaken[2] <<2) ^ (predtaken[3] <<3)) <<1)) & ((1<<(LOGBIASCOLT))-1))


int LSUM;
bool predSC;
bool predchose;
bool pred_inter;

// global history
#define HISTBUFFERLENGTH (1<<18)
uint8_t ghist[HISTBUFFERLENGTH];
int ptghist;

//utilities for computing GEHL indices
folded_history chgehl_i[MAXNHISTGEHL + 1]; 
folded_history chrhsp_i[MAXNRHSP + 1];


int NGEHL;
int NRHSP;
int MAXHISTGEHL;



#define ASSERT(cond) if (!(cond)) {fprintf(stderr,"file %s assert line %d\n",__FILE__,__LINE__); abort();}


#define DECPTR(ptr,size)                        \
{                                               \
  ptr--;                                        \
  if (ptr==(-1)) {                              \
    ptr = size-1;                               \
  }                                             \
}

#define INCSAT(ctr,max) {if (ctr < (max)) ctr++;}

#define DECSAT(ctr,min) {if (ctr > (min)) ctr--;}


// for updating up-down saturating counters
bool
ctrupdate(int8_t & ctr, bool inc, int nbits)
{
  ASSERT(nbits<=8);
  int ctrmin = -(1 << (nbits-1));
  int ctrmax = -ctrmin-1;
  bool issat = (ctr==ctrmax) || (ctr==ctrmin);
  if (inc) {
    INCSAT(ctr,ctrmax);
  } else {
    DECSAT(ctr,ctrmin);
  }
  return issat && ((ctr==ctrmax) || (ctr==ctrmin));
}


void
path_history::init(int hlen)
{
  hlength = hlen;
  h = new unsigned [hlen];
  for (int i=0; i<hlength; i++) {
    h[i] = 0;
  }
  ptr = 0;
}


void 
path_history::insert(unsigned val)
{
  DECPTR(ptr,hlength);
  h[ptr] = val;
}


unsigned & 
path_history::operator [] (int n)
{
  ASSERT((n>=0) && (n<hlength));
  int k = ptr + n;
  if (k >= hlength) {
    k -= hlength;
  }
  ASSERT((k>=0) && (k<hlength));
  return h[k];
}


compressed_history::compressed_history()
{
  reset();
}


void
compressed_history::reset()
{
  comp = 0; // must be consistent with path_history::reset()
}


void 
compressed_history::init(int original_length, int compressed_length, int injected_bits)
{
  olength = original_length;
  clength = compressed_length;
  nbits = injected_bits;
  outpoint = olength % clength;
  ASSERT(clength < 32);
  ASSERT(nbits <= clength);
  mask1 = (1<<clength)-1;
  mask2 = (1<<nbits)-1;
  reset();
}


void 
compressed_history::rotateleft(unsigned & x, int m)
{
  ASSERT(m < clength);
  ASSERT((x>>clength) == 0);
  unsigned y = x >> (clength-m);
  x = (x << m) | y;
  x &= mask1;
}


void 
compressed_history::update(path_history & ph)
{
  rotateleft(comp,1);
  unsigned inbits = ph[0] & mask2;
  unsigned outbits = ph[olength] & mask2;
  rotateleft(outbits,outpoint);
  comp ^= inbits ^ outbits;
}


coltentry::coltentry()
{
  for (int i=0; i<(1<<NPRED); i++) {
    c[i] = ((i>>(NPRED-1)) & 1)? 1:-2;
  }
}


int8_t & 
coltentry::ctr(bool predtaken[NPRED])
{
  int v = 0;
  for (int i=0; i<NPRED; i++) {
    v = (v << 1) | ((predtaken[i])? 1:0);
  }
  return c[v];
}


int8_t & 
colt::ctr(UINT64 pc, bool predtaken[NPRED])
{
  int i = pc & ((1<<LOGCOLT)-1);
  return c[i].ctr(predtaken);
}


bool
colt::predict(UINT64 pc, bool predtaken[NPRED])
{
  return (ctr(pc,predtaken) >= 0);
}


void 
colt::update(UINT64 pc, bool predtaken[NPRED], bool taken)
{
  ctrupdate(ctr(pc,predtaken),taken,COLTBITS);
}


bftable::bftable()
{
  for (int i=0; i<BFTSIZE; i++) {
    freq[i] = 0;
  }
}


int & 
bftable::getfreq(UINT64 pc)
{
  int i = pc % BFTSIZE;
  ASSERT((i>=0) && (i<BFTSIZE));
  return freq[i];
}



void 
subpath::init(int ng, int hist[], int logg, int tagbits, int pathbits, int hp)
{
  ASSERT(ng>0);
  numg = ng;
  ph.init(hist[numg-1]+1);
  chg = new compressed_history [numg];
  chgg = new compressed_history [numg];
  cht = new compressed_history [numg];
  chtt = new compressed_history [numg];
  int ghlen = 0;
  for (int i=numg-1; i>=0; i--) {
    ghlen = (ghlen < hist[numg-1-i]) ? hist[numg-1-i] : ghlen+1;
    chg[i].init(ghlen,logg,pathbits);
    chgg[i].init(ghlen,logg-hp,pathbits);
    cht[i].init(ghlen,tagbits,pathbits);
    chtt[i].init(ghlen,tagbits-1,pathbits);
  }
}


void 
subpath::init(int ng, int minhist, int maxhist, int logg, int tagbits, int pathbits, int hp)
{
  int * h = new int [ng];
  for (int i=0; i<ng; i++) {
    h[i] = minhist * pow((double)maxhist/minhist,(double)i/(ng-1));
  }
  init(ng,h,logg,tagbits,pathbits,hp);
}



void 
subpath::update(UINT64 targetpc, bool taken)
{
  ph.insert((targetpc<<1)|taken);
  for (int i=0; i<numg; i++) {
    chg[i].update(ph);
    chgg[i].update(ph);
    cht[i].update(ph);
    chtt[i].update(ph);
  }
}


unsigned 
subpath::cg(int bank)
{
  ASSERT((bank>=0) && (bank<numg));
  return chg[bank].comp;
}


unsigned 
subpath::cgg(int bank)
{
  ASSERT((bank>=0) && (bank<numg));
  return chgg[bank].comp << (chg[bank].clength-chgg[bank].clength);
}


unsigned 
subpath::ct(int bank)
{
  ASSERT((bank>=0) && (bank<numg));
  return cht[bank].comp;
}


unsigned 
subpath::ctt(int bank)
{
  ASSERT((bank>=0) && (bank<numg));
  return chtt[bank].comp << (cht[bank].clength-chtt[bank].clength);
}


spectrum::spectrum()
{
  size = 0;
  p = NULL;
}


void
spectrum::init(int sz, int ng, int minhist, int maxhist, int logg, int tagbits, int pathbits, int hp)
{
  size = sz;
  p = new subpath [size];
  for (int i=0; i<size; i++) {
    p[i].init(ng,minhist,maxhist,logg,tagbits,pathbits,hp);
  }
}


void
freqbins::init(int nb)
{
  nbins = nb;
  maxfreq = 0;
}


int
freqbins::find(int bfreq)
{
  // find in which frequency bin the input branch frequency falls
  ASSERT(bfreq>=0);
  int b = -1;
  int f = maxfreq;
  for (int i=0; i<nbins; i++) {
    f = f >> FRATIOBITS;
    if (bfreq >= f) {
      b = i;
      break; 
    }
  }
  if (b < 0) {
    b = nbins-1;
  }
  return b;
}


void 
freqbins::update(int bfreq)
{
  if (bfreq > maxfreq) {
    ASSERT(bfreq==(maxfreq+1));
    maxfreq = bfreq;
  }
}


gentry::gentry()
{
  ctr = 0;
  tag = 0;
  u = 0;
}


potage::potage()
{
  b = NULL;
  g = NULL;
  gi = NULL;
  postp = NULL;
  nmisp = 0;
}



void 
potage::init(const char * nm, int ng, int logb, int logg, int tagb, int ctrb, int ppb, int ru, int caph)
{
  ASSERT(ng>1);
  ASSERT(logb<30);
  ASSERT(logg<30);
  name = nm;
  numg = ng;
  bsize = 1 << logb;
  gsize = 1 << logg;
  tagbits = tagb;
  ctrbits = ctrb;
  postpbits = ppb;
  postpsize = 1 << ((1+POSTPEXTRA)*ctrbits+1);
  b = new int8_t [bsize];
  for (int i=0; i<bsize; i++) {
    b[i] = 0;
  }
  g = new gentry * [numg];
  for (int i=0; i<numg; i++) {
    g[i] = new gentry [gsize];
  }
  gi = new int [numg];
  postp = new int8_t [postpsize];
  for (int i=0; i<postpsize; i++) {
    postp[i] = -(((i>>1) >> (ctrbits-1)) & 1);
  }
  allocfail = 0;
  rampup = ru;
  caphist = caph;
}


int 
potage::bindex(UINT64 pc)
{
  return pc & (bsize-1);
}


int 
potage::gindex(UINT64 pc, subpath & p, int bank)
{
  return (pc ^ p.cg(bank) ^ p.cgg(bank)) & (gsize-1);
}


int 
potage::gtag(UINT64 pc, subpath & p, int bank)
{
  return (pc ^ p.ct(bank) ^ p.ctt(bank)) & ((1<<tagbits)-1);
}


int
potage::postp_index()
{
  // post predictor index function
  int ctr[POSTPEXTRA+1];
  for (int i=0; i<=POSTPEXTRA; i++) {
    ctr[i] = (i < (int)hit.size())? getg(hit[i]).ctr : b[bi];
  }
  int v = 0;
  for (int i=POSTPEXTRA; i>=0; i--) {
    v = (v << ctrbits) | (ctr[i] & (((1<<ctrbits)-1)));
  }
  int u0 = (hit.size()>0)? getg(hit[0]).u : 1;
  v = (v << 1) | u0;
  v &= postpsize-1;
  return v;
}


gentry &
potage::getg(int i)
{
  ASSERT((i>=0) && (i<numg));
  return g[i][gi[i]];
}


bool 
potage::condbr_predict(UINT64 pc, subpath & p)
{
  hit.clear();
  bi = bindex(pc);
  for (int i=0; i<numg; i++) {
    gi[i] = gindex(pc,p,i);
    if (g[i][gi[i]].tag == gtag(pc,p,i)) {
      hit.push_back(i);
    }
  }
  predtaken = (hit.size()>0)? (getg(hit[0]).ctr>=0) : (b[bi]>=0); 
  altpredtaken = (hit.size()>1)? (getg(hit[1]).ctr>=0) : (b[bi]>=0);
  ppi = postp_index();
  ASSERT(ppi<postpsize);
  postpredtaken = (postp[ppi] >= 0);
  return postpredtaken;
}


void
potage::uclear()
{
  for (int i=0; i<numg; i++) {
    for (int j=0; j<gsize; j++) {
      g[i][j].u = 0;
    }
  }
}


void
potage::galloc(int i, UINT64 pc, bool taken, subpath & p)
{
  getg(i).tag = gtag(pc,p,i);
  getg(i).ctr = (taken)? 0 : -1;
  getg(i).u = 0;
}


void
potage::aggressive_update(UINT64 pc, bool taken, subpath & p)
{
  // update policy used during ramp up
  bool allsat = true;
  for (int i=0; i<(int)hit.size(); i++) {
    allsat &= ctrupdate(getg(hit[i]).ctr,taken,ctrbits);
  }
  if (hit.size()==0) {
    allsat = ctrupdate(b[bi],taken,ctrbits);
  }
  int i = (hit.size()>0)? hit[0] : numg;
  while (--i >= 0) {
    if (getg(i).u != 0) continue;
    if (! allsat || (p.chg[i].olength <= caphist)) {
      galloc(i,pc,taken,p);
    }
  }
}


void
potage::careful_update(UINT64 pc, bool taken, subpath & p)
{
  // update policy devised by Andre Seznec for the ISL-TAGE predictor (MICRO 2011)
  if (hit.size()>0) {
    ctrupdate(getg(hit[0]).ctr,taken,ctrbits);
    if (getg(hit[0]).u==0) {
      if (hit.size()>1) {
  	ctrupdate(getg(hit[1]).ctr,taken,ctrbits);
      } else {
  	ctrupdate(b[bi],taken,ctrbits);
      }
    }
  } else {
    ctrupdate(b[bi],taken,ctrbits);
  }

  if (mispred) {
    int nalloc = 0;
    int i = (hit.size()>0)? hit[0] : numg;
    while (--i >= 0) {
      if (getg(i).u == 0) {
	galloc(i,pc,taken,p);
	DECSAT(allocfail,0);
	i--;
	nalloc++;
	if (nalloc==MAXALLOC) break;
      } else {
	INCSAT(allocfail,ALLOCFAILMAX);
	if (allocfail==ALLOCFAILMAX) {
	  uclear();
	}
      }
    }
  }

}


bool
potage::condbr_update(UINT64 pc, bool taken, subpath & p)
{
  mispred = (postpredtaken != taken);

  if (mispred) {
    nmisp++;
  }

  if (nmisp < rampup) {
    aggressive_update(pc,taken,p);
  } else {
    careful_update(pc,taken,p);
  }

  // update u bit (see TAGE, JILP 2006)
  if (predtaken != altpredtaken) {
    ASSERT(hit.size()>0);
    if (predtaken == taken) {
      getg(hit[0]).u = 1;
    } else {
      getg(hit[0]).u = 0;
    }
  }

  // update post pred
  ctrupdate(postp[ppi],taken,postpbits);

  return mispred;
}


void 
potage::printconfig(subpath & p)
{
  printf("%s path lengths: ",name.c_str());
  for (int i=numg-1; i>=0; i--) {
    printf("%d ",p.chg[i].olength);
  }
  printf("\n");
}



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


void folded_history::init(int original_length, int compressed_length, int N)
{
  comp = 0;
  OLENGTH = original_length;
  CLENGTH = compressed_length;
  OUTPOINT = OLENGTH % CLENGTH;
}

void folded_history::update (uint8_t * h, int PT)
{
  comp = (comp << 1) ^ h[PT & (HISTBUFFERLENGTH - 1)];
  comp ^= h[(PT + OLENGTH) & (HISTBUFFERLENGTH - 1)] << OUTPOINT;
  comp ^= (comp >> CLENGTH);
  comp = (comp) & ((1 << CLENGTH) - 1);
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////



PREDICTOR::PREDICTOR(void)
{
  sp[0].init(P0_SPSIZE,P0_NUMG,P0_MINHIST,P0_MAXHIST,P0_LOGG,TAGBITS,PATHBITS,P0_HASHPARAM);
  sp[1].init(P1_SPSIZE,P1_NUMG,P1_MINHIST,P1_MAXHIST,P1_LOGG,TAGBITS,PATHBITS,P1_HASHPARAM);
  sp[2].init(P2_SPSIZE,P2_NUMG,P2_MINHIST,P2_MAXHIST,P2_LOGG,TAGBITS,PATHBITS,P2_HASHPARAM);
  sp[3].init(P3_SPSIZE,P3_NUMG,P3_MINHIST,P3_MAXHIST,P3_LOGG,TAGBITS,PATHBITS,P3_HASHPARAM);
  sp[4].init(P4_SPSIZE,P4_NUMG,P4_MINHIST,P4_MAXHIST,P4_LOGG,TAGBITS,PATHBITS,P4_HASHPARAM);

  pred[0].init("G",P0_NUMG,P0_LOGB,P0_LOGG,TAGBITS,CTRBITS,POSTPBITS,P0_RAMPUP,CAPHIST);
  pred[1].init("A",P1_NUMG,P1_LOGB,P1_LOGG,TAGBITS,CTRBITS,POSTPBITS,P1_RAMPUP,CAPHIST);
  pred[2].init("S",P2_NUMG,P2_LOGB,P2_LOGG,TAGBITS,CTRBITS,POSTPBITS,P2_RAMPUP,CAPHIST);
  pred[3].init("s",P3_NUMG,P3_LOGB,P3_LOGG,TAGBITS,CTRBITS,POSTPBITS,P3_RAMPUP,CAPHIST);
  pred[4].init("F",P4_NUMG,P4_LOGB,P4_LOGG,TAGBITS,CTRBITS,POSTPBITS,P4_RAMPUP,CAPHIST);

  bfreq.init(P4_SPSIZE); // number of frequency bins = P4 spectrum size

  initSC();
}



bool   
PREDICTOR::GetPrediction(UINT64 PC)
{
  subp[0] = & sp[0].p[0]; // global path
  subp[1] = & sp[1].p[PC % P1_SPSIZE]; // per-address subpath
  subp[2] = & sp[2].p[(PC>>P2_PARAM) % P2_SPSIZE]; // per-set subpath
  subp[3] = & sp[3].p[(PC>>P3_PARAM) % P3_SPSIZE]; // another per-set subpath
  int f = bfreq.find(bft.getfreq(PC));
  ASSERT((f>=0) && (f<P4_SPSIZE));
  subp[4] = & sp[4].p[f]; // frequency subpath 

  for (int i=0; i<NPRED; i++) {
    predtaken[i] = pred[i].condbr_predict(PC,*subp[i]);
  }

  pred_inter = co.predict(PC,predtaken);
  
  predSC = SCpredict(PC,pred_inter);

  bool finalpred = predSC;

  if (predSC != pred_inter){
    predchose = ChoseSCpredict(PC,pred_inter);
    if (predchose) finalpred = pred_inter;
  }

  return finalpred;
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

void  
PREDICTOR::UpdatePredictor(UINT64 PC, OpType opType, bool resolveDir, bool predDir, UINT64 branchTarget)
{
  UpdateChoseSC(PC, resolveDir,pred_inter);
  UpdateSC(PC,resolveDir,pred_inter);
     
  for (int i=0; i<NPRED; i++) {
    pred[i].condbr_update(PC,resolveDir,*subp[i]);
    subp[i]->update(branchTarget,resolveDir);
  }

  co.update(PC,predtaken,resolveDir);

  bfreq.update(bft.getfreq(PC));
  bft.getfreq(PC)++;

  HistoryUpdate (PC, opType, resolveDir, branchTarget, ptghist, chgehl_i, chrhsp_i, 
                 L_shist[INDLOCAL], BHIST, lastaddr, S_slhist[INDSLOCAL], T_slhist[INDTLOCAL],
                 Q_slhist[INDQLOCAL], P_phist,YHA, GHIST);
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


void    
PREDICTOR::TrackOtherInst(UINT64 PC, OpType opType, bool branchDir, UINT64 branchTarget)
{
  switch (opType) {
  case OPTYPE_RET_UNCOND:
  case OPTYPE_CALL_INDIRECT_UNCOND:
  case OPTYPE_JMP_INDIRECT_UNCOND:
  case OPTYPE_CALL_DIRECT_UNCOND:
  case OPTYPE_JMP_DIRECT_UNCOND:    
    // also update the global path with unconditional branches
    sp[0].p[0].update(branchTarget,true);
    HistoryUpdate (PC, 0, true, branchTarget, ptghist,  chgehl_i, chrhsp_i,
		   L_shist[INDLOCAL], BHIST, lastaddr, S_slhist[INDSLOCAL], T_slhist[INDTLOCAL],
                   Q_slhist[INDQLOCAL], P_phist,YHA, GHIST);
    break;
  default: break;
  }
}




void 
PREDICTOR::initSC()
{

  NRHSP = 80;
  NGEHL = 209;
  MAXHISTGEHL = 1393;
    
  for (int i = 0; i < HISTBUFFERLENGTH; i++)
    ghist[0] = 0;

  ptghist = 0;

  //GEHL initialization
  mgehl[0] = 0;
  mgehl[1] = MINHISTGEHL;
  mgehl[NGEHL] = MAXHISTGEHL;

  for (int i = 2; i <= NGEHL; i++)
    mgehl[i] = (int) (((double) MINHISTGEHL * pow ((double) MAXHISTGEHL / (double) MINHISTGEHL, (double) (i-1) / (double) (NGEHL-1))) + 0.5);

  // just guarantee that all history lengths are distinct

  for (int i = 1; i <= NGEHL; i++)
    if (mgehl[i] <= mgehl[i - 1] + MINSTEP)
      mgehl[i] = mgehl[i - 1] + MINSTEP;

  for (int i = 1; i <= NGEHL; i++)
    chgehl_i[i].init (mgehl[i], LOGGEHL, ((i & 1)) ? i : 1);

  // initialization of GEHL tables

  for (int j = 0; j < (1 << LOGGEHL); j++)
    for (int i = 0; i <= NGEHL; i++)
      GEHL[j][i] = (i & 1) ? -4 : 3;

  // RHSP initialization

  for (int i = 1; i <= NRHSP; i++)
    mrhsp[i] = 6 * i;

  for (int i = 1; i <= NRHSP; i++)
    chrhsp_i[i].init (mrhsp[i], LOGRHSP, ((i & 1)) ? i : 1);

  // initialization of RHSP tables

  for (int j = 0; j < (1 << LOGRHSP); j++)
    for (int i = 0; i <= NRHSP; i++)
      RHSP[j][i] = (i & 1) ? -4 : 3;

  updatethreshold = 100;
  Cupdatethreshold = 11;
    
  for (int i = 0; i < (1 << LOGSIZE); i++)
    Pupdatethreshold[i] = 0;

  for (int i = 0; i < LNB; i++)
    LGEHL[i] = &LGEHLA[i][0];

  for (int i = 0; i < SNB; i++)
    SGEHL[i] = &SGEHLA[i][0];

  for (int i = 0; i < QNB; i++)
    QGEHL[i] = &QGEHLA[i][0];

  for (int i = 0; i < TNB; i++)
    TGEHL[i] = &TGEHLA[i][0];


  for (int i = 0; i < BNB; i++)
    BGEHL[i] = &BGEHLA[i][0];

  for (int i = 0; i < YNB; i++)
    YGEHL[i] = &YGEHLA[i][0];

  for (int i = 0; i < LNB; i++)
    for (int j = 0; j < TABSIZE; j++)
      {
	if (j & 1)
	  {
	    LGEHL[i][j] = -1;
	  }
      }

  for (int i = 0; i < SNB; i++)
    for (int j = 0; j < TABSIZE; j++)
      {
	if (j & 1)
	  {
	    SGEHL[i][j] = -1;
	  }
      }

  for (int i = 0; i < QNB; i++)
    for (int j = 0; j < TABSIZE; j++)
      {
	if (j & 1)
	  {
	    QGEHL[i][j] = -1;
	  }
      }

  for (int i = 0; i < TNB; i++)
    for (int j = 0; j < TABSIZE; j++)
      {
	if (j & 1)
	  {
	    TGEHL[i][j] = -1;
	  }
      }

  for (int i = 0; i < BNB; i++)
    for (int j = 0; j < TABSIZE; j++)
      {
	if (j & 1)
	  {
	    BGEHL[i][j] = -1;
	  }
      }

  for (int i = 0; i < YNB; i++)
    for (int j = 0; j < TABSIZE; j++)
      {
	if (j & 1)
	  {
	    YGEHL[i][j] = -1;
	  }
      }


  //for the choser
  for (int i = 0; i < CGNB; i++)
    CGGEHL[i] = &CGGEHLA[i][0];

  for (int i = 0; i < CLNB; i++) 
    CLGEHL[i] = &CLGEHLA[i][0];

  for (int i = 0; i < CGNB; i++)
    for (int j = 0; j < TABSIZE; j++)
      {
	if (j & 1)
	  {
	    CGGEHL[i][j] = -1;
	  }
      }

  for (int i = 0; i < CLNB; i++)
    for (int j = 0; j < TABSIZE; j++)
      {
	if (j & 1)
	  {
	    CLGEHL[i][j] = -1;
	  }
      }

  for (int j = 0; j < (1 << (LOGBIAS + 1)); j++)
    Bias[j] = (j & 1) ? 15 : -16;

  for (int j = 0; j < (1 << LOGBIASCOLT); j++)
    BiasColt[j] = (j & 1) ? 0: -1;

  for (int i = 0; i < (1 << LOGSIZES); i++)
    for (int j = 0; j < (SWIDTH / SPSTEP) * (1 << SPSTEP); j++)
      {
	if (j & 1)
	  {
	    PERCSLOC[i][j] = -1;
	  }
      }

  for (int i = 0; i < (1 << LOGSIZEQ); i++)
    for (int j = 0; j < (QWIDTH / QPSTEP) * (1 << QPSTEP); j++)
      {
	if (j & 1)
	  {
	    PERCQLOC[i][j] = -1;
	  }
      }

  for (int i = 0; i < (1 << LOGSIZEP); i++)
    for (int j = 0; j < (GWIDTH / GPSTEP) * (1 << GPSTEP); j++)
      {
	if (j & 1)
	  {
	    PERC[i][j] = -1;
	  }
      }

  for (int i = 0; i < (1 << LOGSIZEL); i++)
    for (int j = 0; j < (LWIDTH / LPSTEP) * (1 << LPSTEP); j++)
      {
	if (j & 1)
	  {
	    PERCLOC[i][j] = -1;
	  }
      }

  for (int i = 0; i < (1 << LOGSIZEB); i++)
    for (int j = 0; j < ((BWIDTH / BPSTEP)) * (1 << BPSTEP); j++)
      {
	if (j & 1)
	  {
	    PERCBACK[i][j] = -1;
	  }
      }

  for (int i = 0; i < (1 << LOGSIZEY); i++)
    for (int j = 0; j < (YWIDTH / YPSTEP) * (1 << YPSTEP); j++)
      {
	if (j & 1)
	  {
	    PERCYHA[i][j] = -1;
	  }
      }

  for (int i = 0; i < (1 << LOGSIZEP); i++)
    for (int j = 0; j < (PWIDTH / PPSTEP) * (1 << PPSTEP); j++)
      {
	if (j & 1)
	  {
	    PERCPATH[i][j] = -1;
	  }
      }

}



void 
PREDICTOR::HistoryUpdate (uint64_t PC, uint8_t brtype, bool taken, uint64_t target, int &Y,
			  folded_history * K, folded_history * L, long long &LH, long long &BH,
			  long long &lastaddr, long long &XH, long long &NH,
			  long long &QH, long long &TH,long long &YH, long long &BHIST)
{

  // History skeleton
  bool V = false;

  for (int i=0; i<=7; i++)
    if (LastBR[i]==(int)PC) V = true;

  for (int i=7; i>=1; i--)
    LastBR[i]= LastBR[i-1];

  LastBR[0]=PC;

  if (!V)   YH = (YH << 1) ^ (taken ^ ((PC >> 5) & 1));

  if ((PC < lastaddr - 16) || (abs (PC - lastaddr) >= 128))
    {
      BH = (BH << 1) ^ (PC & 15);
    }

  lastaddr = PC;

  //Path history

  TH = (TH << 1) ^ (taken ^ ((PC >> 5) & 1));


  if (brtype == OPTYPE_RET_COND || brtype == OPTYPE_CALL_INDIRECT_COND || brtype == OPTYPE_JMP_INDIRECT_COND || brtype == OPTYPE_CALL_DIRECT_COND || brtype == OPTYPE_JMP_DIRECT_COND )
    {
      // local history 
      LH = (LH << 1) + (taken);
      QH = (QH << 1) + (taken);
      // second local history + a little bit of path 
      XH = (XH << 1) + (taken);
      XH ^= ((PC >> LOGSECLOCAL) & 15);
      NH = (NH << 1) + (taken);
      NH ^= ((PC >> LOGTLOCAL) & 15);
      // branch history
      BHIST = (BHIST << 1) + taken;
    }


  int T = ((target ^ (target >> 3) ^ PC) << 1) + taken;
  int PATH = PC;
  int8_t DIR = (T & 127);
  T >>= 1;
  PATH >>= 1;
  //update  history
  Y--;
  int inter = (T ^ PC);
  inter = (inter ^ (inter >> 4)) & 1;
  ghist[Y & (HISTBUFFERLENGTH - 1)] = DIR;

 //prepare next index and tag computations 
  for (int i = 1; i <= NGEHL; i++)
    {
      K[i].update (ghist, Y);
    }

  for (int i = 1; i <= NRHSP; i++)
    {
      L[i].update (ghist, Y);
    }

}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

void 
PREDICTOR::UpdateChoseSC(UINT64 PC, bool taken, bool PRED)
{
  bool CRES = (PRED == taken);

  if (PRED != predSC) {
    if ((predchose != CRES) || (std::abs(LCHOSE) < Cupdatethreshold))
      {
	if (predchose != CRES) {
	  Cupdatethreshold += 1;
	} else {
	  Cupdatethreshold -= 1;
	}
	ctrupdate(GCHOSE[indexchose], CRES, CHOSEWIDTH);
	ctrupdate(GCHOSECOLT[indexchosecolt], CRES, CHOSEWIDTH);
	Gupdate (indexchose, CRES, GHIST, CGm, CGGEHL, CGNB,CHOSEWIDTH);
	Gupdate (indexchose, CRES,L_shist[INDLOCAL], CLm,CLGEHL, CLNB,CHOSEWIDTH);
      }
  }
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

bool 
PREDICTOR::ChoseSCpredict(UINT64 PC, bool Tpred)
{
  int CLASS = Tpred +  ((std::abs(LSUM) <  updatethreshold + Pupdatethreshold[INDUPD]) <<1);
  int CLASSCOLT = Tpred + ((predtaken[0] ^ (predtaken[1] <<1) ^ (predtaken[2] <<2) ^ (predtaken[3] <<3)) <<1);

  indexchose = ((INDCHOSE << 2) + CLASS) & (TABSIZE-1);
  indexchosecolt = ((PC<<5) + CLASSCOLT) & (TABSIZE-1);
  LCHOSE = 2*(GCHOSE[indexchose]+1); 
  LCHOSE += 2*(GCHOSECOLT[indexchosecolt]+1);
  LCHOSE += Gpredict (indexchose, GHIST, CGm, CGGEHL, CGNB);
  LCHOSE += Gpredict (indexchose, L_shist[INDLOCAL], CLm, CLGEHL, CLNB);
  return (LCHOSE>=0);
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

void 
PREDICTOR::UpdateSC(UINT64 PC, bool taken, bool PRED)
{
  if ((predSC != taken) || ((abs (LSUM) < updatethreshold + Pupdatethreshold[INDUPD])))
    {
      if (predSC!= taken) {
	updatethreshold += 1;
      } else {
	updatethreshold -= 1;
      }

      if (predSC != taken) {
	Pupdatethreshold[INDUPD] += 1;
      } else {
	Pupdatethreshold[INDUPD] -= 1;
      }

      gehlupdate (PC, taken);
      rhspupdate (PC, taken);
      ctrupdate (Bias[INDBIAS], taken, PERCWIDTH);
      ctrupdate (BiasColt[INDBIASCOLT], taken, PERCWIDTH);
      updateperc (taken, PERC[PC & ((1 << LOGSIZEG) - 1)], GHIST, GPSTEP, GWIDTH);
      updateperc (taken, PERCLOC[PC & ((1 << LOGSIZEL) - 1)], L_shist[INDLOCAL], LPSTEP, LWIDTH);
      updateperc (taken, PERCBACK[PC & ((1 << LOGSIZEB) - 1)], BHIST, BPSTEP, BWIDTH);
      updateperc (taken, PERCYHA[PC & ((1 << LOGSIZEB) - 1)], YHA, YPSTEP, YWIDTH);
      updateperc (taken, PERCPATH[PC & ((1 << LOGSIZEP) - 1)], P_phist, PPSTEP, PWIDTH);
      updateperc (taken, PERCSLOC[PC & ((1 << LOGSIZES) - 1)], S_slhist[INDSLOCAL], SPSTEP, SWIDTH);
      updateperc (taken, PERCTLOC[PC & ((1 << LOGSIZES) - 1)], T_slhist[INDTLOCAL], SPSTEP, SWIDTH);
      updateperc (taken, PERCQLOC[PC & ((1 << LOGSIZEQ) - 1)], Q_slhist[INDQLOCAL], QPSTEP, QWIDTH);

      Gupdate (PC, taken, L_shist[INDLOCAL], Lm, LGEHL, LNB, PERCWIDTH);
      Gupdate (PC, taken, S_slhist[INDSLOCAL], Sm, SGEHL, SNB, PERCWIDTH);
      Gupdate (PC, taken, T_slhist[INDTLOCAL], Tm, TGEHL, TNB, PERCWIDTH);
      Gupdate (PC, taken, Q_slhist[INDQLOCAL], Qm, QGEHL, QNB, PERCWIDTH);
      Gupdate (PC, taken, BHIST, Bm, BGEHL, BNB, PERCWIDTH);
      Gupdate (PC, taken, YHA, Ym, YGEHL, YNB, PERCWIDTH);
    }
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


bool 
PREDICTOR::SCpredict(UINT64 PC, bool PRED)
{
  LSUM = 0;
  predict_gehl (PC);
  predict_rhsp (PC);
  LSUM += SUMGEHL;
  LSUM += SUMRHSP;

  int16_t ctr = Bias[INDBIAS];
  LSUM += 2* (2 * ctr + 1);

  ctr = BiasColt[INDBIASCOLT];
  LSUM += 2 * (2 * ctr + 1);

  LSUM += percpredict (PC, GHIST, PERC[PC & ((1 << LOGSIZEG) - 1)], GPSTEP, GWIDTH);
  LSUM += percpredict (PC, L_shist[INDLOCAL], PERCLOC[PC & ((1 << LOGSIZEL) - 1)], LPSTEP, LWIDTH);
  LSUM += percpredict (PC, BHIST, PERCBACK[PC & ((1 << LOGSIZEB) - 1)], BPSTEP, BWIDTH);
  LSUM += percpredict (PC, YHA, PERCYHA[PC & ((1 << LOGSIZEY) - 1)], YPSTEP, YWIDTH);
  LSUM += percpredict (PC, P_phist, PERCPATH[PC & ((1 << LOGSIZEP) - 1)], PPSTEP, PWIDTH);
  LSUM += percpredict (PC, S_slhist[INDSLOCAL], PERCSLOC[PC & ((1 << LOGSIZES) - 1)], SPSTEP, SWIDTH);
  LSUM += percpredict (PC, T_slhist[INDTLOCAL], PERCTLOC[PC & ((1 << LOGSIZET) - 1)], TPSTEP, TWIDTH);
  LSUM += percpredict (PC, Q_slhist[INDQLOCAL], PERCQLOC[PC & ((1 << LOGSIZEQ) - 1)], QPSTEP, QWIDTH);

  LSUM += Gpredict (PC, L_shist[INDLOCAL], Lm, LGEHL, LNB);
  LSUM += Gpredict (PC, S_slhist[INDSLOCAL], Sm, SGEHL, SNB);
  LSUM += Gpredict (PC, T_slhist[INDTLOCAL], Tm, TGEHL, TNB);
  LSUM += Gpredict (PC, Q_slhist[INDQLOCAL], Qm, QGEHL, QNB);
  LSUM += Gpredict (PC, BHIST, Bm, BGEHL, BNB); 
  LSUM += Gpredict (PC, YHA, Ym, YGEHL, BNB); 
        
  return (LSUM>=0);
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
//Functions  for the statiscal corrector

void 
PREDICTOR::predict_gehl (uint64_t PC)
{
  //index computation     
  for (int i = 1; i <= NGEHL; i++) {
    GEHLINDEX[i] = gehlindex (PC, i);
  }
  GEHLINDEX[0] = PC & ((1 << LOGGEHL) - 1);

  // SUMGEHL is centered
  SUMGEHL = 0;
  for (int i = 0; i <= NGEHL; i++) {
    SUMGEHL += 2 * GEHL[GEHLINDEX[i]][i] + 1;
  }
}


void 
PREDICTOR::gehlupdate(uint64_t PC, bool taken)
{
  //update the GEHL  predictor tables
  for (int i = NGEHL; i >= 0; i--)
    ctrupdate (GEHL[GEHLINDEX[i]][i], taken, PERCWIDTH);
}


void 
PREDICTOR::predict_rhsp (uint64_t PC)
{
  //index computation     
  for (int i = 1; i <= NRHSP; i++) {
    RHSPINDEX[i] = rhspindex (PC, i);
  }
  RHSPINDEX[0] = PC & ((1 << LOGRHSP) - 1);

  // SUMRHSP is centered
  SUMRHSP = 0;
  for (int i = 1; i <= NRHSP; i++)
    SUMRHSP += 2 * RHSP[RHSPINDEX[i]][i] + 1;
}


void 
PREDICTOR::rhspupdate (uint64_t PC, bool taken)
{
  for (int i = NRHSP; i >= 1; i--)
    ctrupdate (RHSP[RHSPINDEX[i]][i], taken, PERCWIDTH);
}


int 
PREDICTOR::percpredict (int PC, long long BHIST, int8_t * line, int PSTEP, int WIDTH)
{
  PERCSUM = 0;
  long long bhist = BHIST;
  int PT = 0;

  for (int i = 0; i < WIDTH; i += PSTEP)
    {
      int index = bhist & ((1 << PSTEP) - 1);
      int8_t ctr = line[PT + index];

      PERCSUM += 2 * ctr + 1;

      bhist >>= PSTEP;
      PT += (1 << PSTEP);
    }

  return PERCSUM;
}


void 
PREDICTOR::updateperc (bool taken, int8_t * line, long long BHIST, int PSTEP, int WIDTH)
{
  int PT = 0;
  long long bhist = BHIST;

  for (int i = 0; i < WIDTH; i += PSTEP)
    {
      int index = bhist & ((1 << PSTEP) - 1);
      ctrupdate (line[PT + index], taken, PERCWIDTH);
      bhist >>= PSTEP;
      PT += (1 << PSTEP);
    }
}


int 
PREDICTOR::Gpredict (UINT64 PC, long long BHIST, int *length, int8_t ** tab, int NBR)
{
  PERCSUM = 0;

  for (int i = 0; i < NBR; i++)
    {
      long long bhist = BHIST & ((long long) ((1 << length[i]) - 1));

      int index = (((long long) PC) ^ bhist ^ (bhist >> (LOGTAB - i)) ^ (bhist >> (40 - 2 * i)) ^ (bhist >> (60 - 3 * i))) & (TABSIZE - 1);

      int8_t ctr = tab[i][index];
      PERCSUM += 2 * ctr + 1;
    }
  return PERCSUM;
}


void 
PREDICTOR::Gupdate (UINT64 PC, bool taken, long long BHIST, int *length, int8_t ** tab, int NBR, int WIDTH)
{
  for (int i = 0; i < NBR; i++)
    {
      long long bhist = BHIST & ((long long) ((1 << length[i]) - 1));

      int index = (((long long) PC) ^ bhist ^ (bhist >> (LOGTAB - i)) ^ (bhist >> (40 - 2 * i)) ^ (bhist >> (60 - 3 * i))) & (TABSIZE - 1);

      ctrupdate (tab[i][index], taken, WIDTH);
    }
}


// index function for the GEHL tables
//FGEHL serves to mix path history

int 
PREDICTOR::FGEHL (int A, int size, int bank)
{
  int SH = (bank % LOGGEHL);
  A = A & ((1 << size) - 1);
  int A1 = (A & ((1 << LOGGEHL) - 1));
  int A2 = (A >> LOGGEHL);
  A2 = ((A2 << SH) & ((1 << LOGGEHL) - 1)) + (A2 >> (LOGGEHL - SH));
  A = A1 ^ A2;
  A = ((A << SH) & ((1 << LOGGEHL) - 1)) + (A >> (LOGGEHL - SH));
  return A;
}


int
PREDICTOR::gehlindex (uint64_t PC, int bank)
{
  int index = PC ^ (PC >> ((mgehl[bank] % LOGGEHL) + 1)) ^ chgehl_i[bank].comp;
  return index & ((1 << LOGGEHL) - 1);
}


// index function for the RHSP tables
//FRHSP serves to mix path history
int 
PREDICTOR::FRHSP (int A, int size, int bank)
{
  int SH = (bank % LOGRHSP);
  A = A & ((1 << size) - 1);
  int A1 = (A & ((1 << LOGRHSP) - 1));
  int A2 = (A >> LOGRHSP);
  A2 = ((A2 << SH) & ((1 << LOGRHSP) - 1)) + (A2 >> (LOGRHSP - SH));
  A = A1 ^ A2;
  A = ((A << SH) & ((1 << LOGRHSP) - 1)) + (A >> (LOGRHSP - SH));
  return A;
}


int 
PREDICTOR::rhspindex (uint64_t PC, int bank)
{
  int index = PC ^ (PC >> ((mrhsp[bank] % LOGRHSP) + 1)) ^ chrhsp_i[bank].comp;
  if (bank > 1) {
    index ^= chrhsp_i[bank - 1].comp;
  }
  if (bank > 3) {
    index ^= chrhsp_i[bank / 3].comp;
  }
  return index & ((1 << LOGRHSP) - 1);
}

