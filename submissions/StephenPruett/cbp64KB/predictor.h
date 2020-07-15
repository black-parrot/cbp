#include "utils.h"

#include "bt9.h"
#include "bt9_reader.h"

#include <inttypes.h>
#include <inttypes.h>
#include <math.h>
#include <vector>
#include <map>

#ifndef _PREDICTOR_H_
#define _PREDICTOR_H_

/*SC PREDICTOR DEFINES*/
#define SC //Statistical Corrector, SC
#define LOCALH //Local History component

//large local history elements
#define LOGLOCAL 9
#define NLOCAL (1<<LOGLOCAL)
#define LOGLNB 8
#define LNB 3

//small local history elements
#define LOGSECLOCAL 2
#define NSECLOCAL (1<<LOGSECLOCAL)	//Number of second local histories
#define LOGSNB 7
#define SNB 4

//third local history elements
#define LOGTNB 7
#define TNB 3

//return-stack associated history component
#define LOGPNB 7
#define PNB 4
#define HSTACKSIZE 4
#define logHSTACKSIZE 2

#define GNB 4
#define LOGGNB 7  //global branch GEHL
#define CONFWIDTH 7		//for the counters in the choser
#define LOGSIZEUP 3

#define LOGBIAS 5
#define PERCWIDTH 6 //Statistical corrector maximum counter width 

/*IMLI COMPONENT DEFINES*/
//#define IMLI
#define IMLISIC
#define IMLIOH

#define MAXIMLIcount 1023

#define LOGINB 10 // (LOG of IMLI-SIC table size +1)
#define INB 1

#define SHIFTFUTURE 6 // (PC<<6) +IMLIcount to index the Outer History table
#define PASTSIZE 16
#define OHHISTTABLESIZE 1024 //

#define LOGFNB 9 //256 entries
#define FNB 1

/*LOOP PREDICTOR DEFINES*/
#define LOOPPREDICTOR 
// log of number of table entries for loop predictor
#define LOGL 7
// we predict only loops with less than 1K iterations
#define WIDTHNBITERLOOP 10
//tag width in the loop predictor
#define LOOPTAG 10	
//confidence treshold for the confidence of the loop predictor
#define CONFLOOP 15   




//number of entries to allocate on mispredict = NNN+1
#define NNN 1

//counter and useful bit width in tagged tables
#define CWIDTH 3
#define UWIDTH 2

//Size of Bimodal table
#define LOGB 14
#define ALLOCFILTER_EN false
#define VCACHE_EN true
#define VCACHE_FP (1.0/16)

//number of global history registers
#define NHIST 16
//number of tagged tables
#define NUM_TABLES 64
//log of the max number of tables that can be combined into 1
#define LOG_MAX_TABLES 3

//History Sizes
uint32_t m[NHIST+1]={0,4,6,8,12,18,25,37,53,77,111,161,233,338,489,708,1024}; // history lengths


//tag and table sizes for HPS
#define LOG_TABLE_SIZE 9
#define TAG_SIZE 10
uint32_t bestConfig[NHIST+1] =     {0,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};

#define PRINT_CONFIG false

//reconfiguration parameters for HPS
#define USEFUL_THRESHOLD 1.0
#define NUM_RECONFIGS 6
#define CONFIG_THRESH_MISPRED 1
#define BIMOD_ALLOC false
#define CONTINUOUS_RECONFIG true
#define SNAPSHOT_PERIOD 50000
#define ALLOC_POLICY (true || HYBRID_POLICY)
#define CONF_POLICY (true || HYBRID_POLICY)
#define UPDATE_OVERWRITE false
#define OVERWRITE_KNOB true
#define HYBRID_POLICY true
#define HYBRID_MISSES 8000

//To get the predictor storage budget on stderr  uncomment the next line
//#define PRINTSIZE


/* TAGE PARTS: */
class SaturatingSignedCounter {
  private:
    int32_t width;
    int32_t value;
    int32_t max_value;
    int32_t min_value;
    SaturatingSignedCounter() {}
    inline int32_t wrap () {
      return (value > max_value)? max_value : (value < min_value)? min_value : value;
    }
  public:
    SaturatingSignedCounter (int32_t w) {
      width = w;
      value = 0;
      max_value = (1 << (width-1))-1;
      min_value = -(1 << (width-1));
    }
    int32_t inc (bool p) {
      value = p? value+1 : value-1;
      value = wrap();
      return value;
    }
    bool getPrediction() {
      return (value >= 0);
    }
    bool isUnbiased() {
      return (value == 0) || (value == -1);
    }
    int32_t getSize () {
      return width;
    }
    int32_t get () {
      return value;
    }
    void reset (bool p) {
      value = p? 0 : -1;
    }
    int32_t getValue() {
      return value;
    }
    void setValue(int32_t v) {
      value = v;
    }
};

class SaturatingUnsignedCounter {
  private:
    int32_t width;
    int32_t value;
    int32_t max_value;
    SaturatingUnsignedCounter() {}
    inline int32_t wrap () {
      return (value > max_value)? max_value : (value < 0)? 0 : value;
    }
  public:
    SaturatingUnsignedCounter (int32_t w) {
      width = w;
      value = 0;
      max_value = (1 << width) - 1;
    }
    int32_t inc (bool p) {
      value = p? value+1 : value-1;
      value = wrap();
      return value;
    }
    int32_t get () {
      return value;
    }
    int32_t getSize () {
      return width;
    }
    bool isValid () {
      return value != 0;
    }
    void reset () {
      value = 0;
    }
    void setMax () {
      value = max_value;
    }
    bool isMax () {
      return value == max_value;
    }
};

class BimodalTable {
  private:
    int32_t tableSize;
    int32_t hysteresisRatio;
    int32_t predictionSize;
    int32_t hysteresisSize;
    int32_t maxValue;
    int32_t minValue;
    int8_t *prediction;
    int8_t *hysteresis;
    bool highConfidence;
  public:
    BimodalTable () {
      tableSize = (1<<LOGB);
      hysteresisRatio = 2; //2^x:1 prediction:hysteresis
      predictionSize = tableSize;
      hysteresisSize = predictionSize >> hysteresisRatio;
      prediction = new int8_t [predictionSize];
      hysteresis = new int8_t [hysteresisSize];
      // counter states: 00 -> 01 -> 10 -> 11
      maxValue = 0x03;
      minValue = 0x00;
      // 01 is init state
      for (int i = 0; i < predictionSize; ++i)
        prediction[i] = 0x00;
      for (int i = 0; i < hysteresisSize; ++i)
        hysteresis[i] = 0x01;
      highConfidence = false;
    }
    ~BimodalTable () {
      delete [] prediction;
      delete [] hysteresis;
    }
    int32_t getSize () {
      return predictionSize + hysteresisSize;
    }
    int32_t hash (int32_t PC) {
      return PC & (tableSize - 1);
    }

    bool getPrediction (int32_t PC)
    {
      int32_t index = hash( PC );
      int32_t value = prediction[index] << 1;
      assert((index >> hysteresisRatio) < hysteresisSize);
      value |= hysteresis[index >> hysteresisRatio];
      highConfidence = ((value == 3) || (value == 0));
      return (prediction[index] == 0x01);
    }

    bool isConfident(void) {
      return highConfidence;
    }

    int32_t update (int32_t pc, bool p) {
      int32_t index = hash(pc);
      int32_t value = prediction[index] << 1;
      assert((index >> hysteresisRatio) < hysteresisSize);
      value |= hysteresis[index >> hysteresisRatio];
      value = p? value+1 : value-1;
      value = (value > maxValue)? maxValue : (value < minValue)? minValue : value;
      hysteresis[index >> hysteresisRatio] = value & 0x01;
      prediction[index] = (value >> 1) & 0x01;
      return value;
    }
};

struct TaggedTableEntry {
  SaturatingSignedCounter counter;
  SaturatingUnsignedCounter useful;
  uint32_t tag;
  uint32_t lru;
  TaggedTableEntry() : counter(CWIDTH), useful(UWIDTH) {
    tag = 0;
    lru = 0;
  }
  uint32_t getSize() {
    return counter.getSize() + useful.getSize();
  }
  bool isUnused() {
#if OVERWRITE_KNOB
    return useful.isMax() && counter.isUnbiased();
#else
    return false;
#endif
  }
  bool isUseful() {
    return useful.isValid() && !isUnused();
  }
  void set (uint32_t tag_, bool resolveDir) {
    tag = tag_;
    counter.reset(resolveDir);
#if OVERWRITE_KNOB
    useful.setMax();
#endif
  }
};

class TaggedTable {
  private:
    int32_t index;
    int32_t logTableSize;
    int32_t tableSize;
    int32_t tagWidth;
    TaggedTableEntry *table;
    TaggedTable() {}

  public:
    TaggedTable (int32_t i, int32_t tagWidth_, int32_t logTableSize_) {
      index = i;
      logTableSize = logTableSize_;
      tableSize = 1<<logTableSize;
      tagWidth = tagWidth_;
      table = new TaggedTableEntry [tableSize];
    }
    virtual ~TaggedTable () {
      delete [] table;
    }
    int32_t getSize() {
      return tableSize*(table[0].getSize() + tagWidth);
    }
    inline TaggedTableEntry* getEntry(int32_t idx) {
      assert(idx < tableSize);
      return &table[idx];
    }
    uint32_t decAllUseful () {
      uint32_t numUsefulBeforeDec = 0;
      for (int i = 0; i < tableSize; ++i) {
        if (table[i].useful.get() > 0) numUsefulBeforeDec++;
        if (!table[i].isUnused())
          table[i].useful.inc(false);
      }
      return numUsefulBeforeDec;
    }
};


//==========================================================
// History Code:

// [PPM, page 4] Discusses how to build a low latency FoldedHistory register
struct HistoryRegister {
  public:
  uint32_t size;
  uint32_t head;
  std::vector<bool> history;
  long long history_l;

  void init(uint32_t s) {
    size = s;
    history.resize(size);

    for (uint32_t i = 0; i < size; ++i) {
      history[i] = false;
    }
    history_l = 0;

    head = 0;
  }

  HistoryRegister () {}

  HistoryRegister (uint32_t s) {
    init (s);
  }

  void push (bool p) {
    head = (head+1)%size;
    history[head] = p;

    history_l <<= 1;
    history_l += (p & 0x1);
  }

  bool operator[] (const uint32_t i) {
    uint32_t index = (head+size-i)%size;
    assert(index < size);
    return history[index];
  }

  void print () {
    printf("History");
    for (uint32_t i = 0; i < size; ++i) {
      printf("%d, ", (bool) history[(head-i)%size]);
    }
    printf("\n");
  }

  long long getHistory() {
    return history_l;
  }

  uint32_t getSize() {
    return size;
  }
};

class FoldedHistory {
  private:
    uint32_t inputWidth; // size of emulated history register
    uint32_t outputWidth; // size of folded register
    uint32_t maxOutputWidth; //first width register is set to. Used to calc size.
    int32_t remainder;
    int32_t value;
    HistoryRegister *ghr;

    FoldedHistory () {}

  public:
    FoldedHistory (HistoryRegister *g, uint32_t iw, uint32_t ow) {
      inputWidth = iw;
      outputWidth = ow;
      maxOutputWidth = outputWidth;
      ghr = g;

      //using a 32-bit integer as register
      // -need an extra bit, so max is 31 bits...
      assert (outputWidth < 32);
      assert (outputWidth != 0);
      remainder = inputWidth % outputWidth;
      value = 0;
    }

    //Expectation is that FoldedHistory push is called
    //  after HistoryRegister push
    void update () {
      //input bit most recent shifted into ghr
      bool inBit = (*ghr)[0];

      //Shift in bit
      value = (value << 1) | (inBit? 0x01 : 0x00);

      //Fold shifted-out bit in
      value = value ^ (value >> outputWidth);
      value = value & ((1 << outputWidth) - 1);

      //Get bit to shift out from Global History
      bool outputBit = (*ghr)[inputWidth];
      int32_t outputValue = (outputBit)? (0x01 << (remainder)) : 0x0;

      //Shift out bit
      value = value ^ outputValue;
    }

    inline int32_t get () {
      return value;
    }

    void reset () {
      value = 0;
    }

    uint32_t getSize() {
      return maxOutputWidth;
    }
};

class History {
  public:
    uint32_t numHistories;
  private: // maybe should be protected
    uint32_t *historySizes, *hashOutputSizes, *tagWidths;

    //local -- always private
    const uint32_t pathHistoryWidth;
    const uint32_t HISTBUFFERLENGTH;
    const bool useTargetForPath;

    HistoryRegister *ghist;
    FoldedHistory **tableFHist;
    FoldedHistory **tag1FHist;
    FoldedHistory **tag2FHist;
    uint64_t phist;

    //Default constructor is not intended to be used
    History() : pathHistoryWidth(27), HISTBUFFERLENGTH(4096), useTargetForPath(false) {}

    // the index functions for the tagged tables uses path history as in the OGEHL predictor
    //F serves to mix path history: not very important impact
    uint32_t pathHash (uint64_t pathHistory, uint32_t minWidth, uint32_t bank, uint32_t idxSize)
    {
      int32_t A1, A2; //temp variables

      //truncate path history to index size
      pathHistory = (pathHistory & ((1 << minWidth) - 1));

      A1 = (pathHistory & ((1 << idxSize) - 1));

      //Take high part of path history and left rotate it by "bank" ammount
      // this is just to generate a unique hash for each bank
      A2 = (pathHistory >> idxSize);
      A2 = ((A2 << bank) & ((1 << idxSize) - 1)) + (A2 >> (idxSize - bank));

      //fold rotated high part of path into low part of path
      pathHistory = A1 ^ A2;

      //left rotate that chunk by "bank"
      pathHistory = ((pathHistory << bank) & ((1 << idxSize) - 1)) + (pathHistory >> (idxSize - bank));
      return pathHistory;
    }

  public:
    // constructor for HPS tage, equal table and tag sizes
    History (uint32_t numHistories_, uint32_t *historySizes_, uint32_t *tagWidths_, uint32_t *hashOutputSizes_)
                                       : pathHistoryWidth(27), HISTBUFFERLENGTH(historySizes_[numHistories_]+1), useTargetForPath(false) {
      numHistories = numHistories_;
      historySizes = historySizes_;
      hashOutputSizes = hashOutputSizes_;
      tagWidths = tagWidths_;
      phist = 0;
      ghist = new HistoryRegister (HISTBUFFERLENGTH);
      tableFHist = new FoldedHistory* [numHistories];
      tag1FHist = new FoldedHistory* [numHistories];
      tag2FHist = new FoldedHistory* [numHistories];
      for (uint32_t i = 1; i <= numHistories; ++i) {
        tableFHist[i] = new FoldedHistory (ghist, historySizes[i], hashOutputSizes[i]);
        tag1FHist[i] = new FoldedHistory (ghist, historySizes[i], tagWidths[i]);
        tag2FHist[i] = new FoldedHistory (ghist, historySizes[i], tagWidths[i]-1);
      }
    }

    ~History () {
      delete ghist;
    }

    void print () {
      ghist->print();
    }

    long long getGlobalHistory() {
      return ghist->getHistory();
    }

    uint64_t getPHist () {
      return phist;
    }

    // getIndex computes a full hash of PC, global and path histories
    uint32_t getIndex (UINT32 PC, uint32_t bank) {
      uint32_t index;
      uint32_t idxSize = hashOutputSizes[bank];
      uint32_t historySize = historySizes[bank];

      // assumes that pathHistoryWidth is always greater than idxSize
      uint32_t minWidth = (historySize > pathHistoryWidth) ? pathHistoryWidth : idxSize;

      //index is a mix of three components: PC, history, and path
      uint32_t pcPart =  PC ^ (PC >> (std::abs((int32_t) idxSize - (int32_t) bank) + 1));
      uint32_t histPart =  tableFHist[bank]->get();
      uint32_t pathPart =  pathHash(phist, minWidth, bank, idxSize);
      index = pcPart ^ histPart ^ pathPart;
      return (index & ((1 << idxSize) - 1));
    }

    //  tag computation
    uint32_t getTag (UINT32 PC, uint32_t bank) {
      uint32_t tagWidth = tagWidths[bank];
      //uint32_t historySize = historySizes[bank];

      uint32_t tag = PC //PC component
        ^ tag1FHist[bank]->get() ^ (tag2FHist[bank]->get() << 1); // history component

      return (tag & ((1 << tagWidth) - 1));
    }

    uint32_t getSize () {
      uint32_t size = 0;
      for (uint32_t i = 1; i <= numHistories; ++i) {
        size += tableFHist[i]->getSize();
        size += tag1FHist[i]->getSize();
        size += tag2FHist[i]->getSize();
      }
      return ghist->getSize() + pathHistoryWidth + size;
    }

    void updateHistory (UINT32 PC, uint16_t brtype, bool taken, uint32_t target) {
      //special treatment for unconditional branchs;
      int maxt;

      bool isConditionalBranch = false;
      bool isUnConditionalBranch = false;
      bool isReturn = false;
      bool isCallDirect = false;
      switch (brtype) {
        case OPTYPE_RET_COND:
          isReturn = true;
          isConditionalBranch = true;
          break;
        case OPTYPE_CALL_DIRECT_COND:
          isCallDirect = true;
        case OPTYPE_JMP_DIRECT_COND:
        case OPTYPE_JMP_INDIRECT_COND:
        case OPTYPE_CALL_INDIRECT_COND:
          isConditionalBranch = true;
          break;

        case OPTYPE_RET_UNCOND:
          isReturn = true;
          isUnConditionalBranch = true;
          break;
        case OPTYPE_CALL_DIRECT_UNCOND:
          isCallDirect = true;
        case OPTYPE_JMP_DIRECT_UNCOND:
        case OPTYPE_JMP_INDIRECT_UNCOND:
        case OPTYPE_CALL_INDIRECT_UNCOND:
          isUnConditionalBranch = true;
          break;

        default: break;
      }
      (void)isUnConditionalBranch ; //surpress unused paramter warning for inUnConditionalBranch. Delete the line if inUnConditionalBranch is used
      (void)isReturn ; //surpress unused paramter warning for inUnConditionalBranch. Delete the line if inUnConditionalBranch is used
      (void)isCallDirect ; //surpress unused paramter warning for inUnConditionalBranch. Delete the line if inUnConditionalBranch is used

      if (isConditionalBranch) {
        //1) Push branch direction into history register
        //2) Push low bit of PC (target?) into path history register
        maxt = 1;
      } else {
       // 1) Push branch direction and low bits of PC into history register
       // 2) Push several bits of the PC (target?) into the path history register
        maxt = 4;
      }

      int historyBits = ((PC) << 1) + taken; // combination of path and history bits
      int PATH = (useTargetForPath)? target : PC;

      for (int t = 0; t < maxt; t++) {
        //update  history
        bool brDir = (historyBits & 1);
        ghist->push(brDir); // global history register
        // update folded history
        for (uint32_t i = 1; i <= numHistories; ++i) {
          tableFHist[i]->update();
          tag1FHist[i]->update();
          tag2FHist[i]->update();
        }

        //update path
        int PATHBIT = (PATH & 127);
        phist = (phist << 1) ^ PATHBIT;

        //Get next bit (if unconditional branch we want to use more of path)
        historyBits >>= 1;
        PATH >>= 1;
      }
    }
};

class AltPredTable {
  private:
    uint32_t counterSize;
    uint32_t logTableSize;
    SaturatingSignedCounter **table[2];

    AltPredTable () {}

  public:
  AltPredTable (uint32_t logTableSize_) {
    counterSize = 4;
    logTableSize = logTableSize_;
    for (int i = 0; i < 2; ++i) {
      table[i] = new SaturatingSignedCounter* [1 << logTableSize];
      for (int j = 0; j < (1 << logTableSize); ++j) {
        table[i][j] = new SaturatingSignedCounter(counterSize);
      }
    }
  }
  ~AltPredTable () {
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < (1 << logTableSize); ++j) {
        delete table[i][j];
      }
      delete table[i];
    }
  }

  bool index (UINT32 PC, bool pred, uint32_t bank, uint32_t numHistories) {
    int index = (PC & ((1 << logTableSize) - 1)) ^ (pred? 0x01 : 0x00);
    return table[bank > (numHistories/3)][index]->getPrediction();
  }

  void update (UINT32 PC, bool pred, uint32_t bank, uint32_t numHistories, bool p) {
    int index = (PC & ((1 << logTableSize) - 1)) ^ (pred? 0x01 : 0x00);
    table[bank > (numHistories/3)][index]->inc(p);
  }

  uint32_t getSize() {
    return (1 << logTableSize) * 2 * table[0][0]->getSize();
  }
};

struct LinearRandomRegister {
  private :
    uint32_t size;
    LinearRandomRegister () {}

  public:
  static int32_t seed;
  LinearRandomRegister (uint32_t size_) {
    size = size_;
  }

  uint32_t getSize() {
    return size;
  }

  int32_t get (uint64_t phist) {
    seed++;
    seed ^= phist;
    seed = (seed >> 21) + (seed << 11);
    return (seed) & ((1<<size)-1);
  }
};

int32_t LinearRandomRegister::seed = 0;

class IdxRegister {
  public:
    uint32_t numSlots;
    uint32_t slotSize;
    uint32_t *value;

    IdxRegister (uint32_t numSlots_, uint32_t slotSize_) {
      numSlots = numSlots_;
      slotSize = slotSize_;
      //extra entry just for allignment, 0'th entry is not used
      //  and thus does not count towards size;
      value = new uint32_t [numSlots+1];
      for (uint32_t i = 0; i <= numSlots; ++i) {
        value[i] = i;
      }
    }
    ~IdxRegister () {
      delete [] value;
    }
    uint32_t getSize () {
      return slotSize * numSlots;
    }
    IdxRegister& operator= (const uint32_t *a) {
      for (uint32_t i = 1; i <= numSlots; ++i) {
        value[i] = a[i];
      }
      return *this;
    }
    uint32_t& operator[] (const uint32_t& i) {
      assert (i != 0);
      assert (i <= numSlots);
      return value[i];
    }
    //same as operator [], but more intuitive use for pointers
    uint32_t& get (const uint32_t& i) {
      assert (i != 0);
      assert (i <= numSlots);
      return value[i];
    }
};

class Stats {
  public:
    uint32_t predictions;
    uint32_t mispredictions;
    uint32_t allocations; 
    uint32_t conflicts; //used
    uint32_t attempted; // used
    uint32_t overwriteUnused;

    void reset() {
      predictions = 0;
      mispredictions = 0;
      allocations = 0;
      conflicts = 0;
      attempted = 0;
      overwriteUnused = 0;
    }
    Stats () {
      reset();
    }
    ~Stats () {}
    uint32_t getSize() {
      //stats cannot exceed the shapshot period
      return 2*log2(SNAPSHOT_PERIOD);
      //Only two of the stats are acutally used by predictor. Others are just nice to look at. 
    }

};

//http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5066970
class VictimCache {
  private:
    bool *active_v[2];
    bool *active_p[2];
    int32_t *hashingOutputs;
    uint32_t size;
    int32_t active;
    int32_t dataWidth;
    int32_t numInserted[2];
    int32_t activeMaxInserts;
    int32_t numHashingFunctions;
    double falsePositiveRatio;
    double activeFalsePositiveRatio;
    const double ln2 = 0.69314718; //natural log of 2
    VictimCache () {}
    inline int32_t next () {
      return (active+1) & 0x01;
    }
    void flush (uint32_t a) {
      for (uint32_t i = 0; i < size; ++i) {
        active_v[a][i] = false;
      }
      numInserted[a] = 0;
    }
    bool isFull () {
      return numInserted[active] > activeMaxInserts;
    }
    void swap () {
      active = next();
      flush(active);
    }
    int32_t hash (uint32_t index, int32_t hashNum) {
      //1) right rotate index by hashNum
      index = index & ((1 << dataWidth) - 1); // dataWidth can at most be 31 to avoid sign extension in the next phase
      index = (index >> hashNum) | (index << (dataWidth - hashNum));
      index = index & ((1 << dataWidth) - 1);

      //2) Fold dataWidth down to size
      uint32_t output = 0;
      while (index != 0) {
        output ^= (index & (size - 1));
        index >>= ((int32_t) log2(size));
      }
      return output;
    }
    int32_t contains (uint32_t idx, int32_t a) {
      int32_t numFalse = 0;
      for (int32_t i = 0; i < numHashingFunctions; ++i) {
        hashingOutputs[i] = hash(idx, i);
        
        if (active_v[a][hashingOutputs[i]] == false) return -1;

        if (active_p[a][hashingOutputs[i]] == false) numFalse++;
      }
      //majority function
      if (numFalse > numHashingFunctions - numFalse) return 0;
      else return 1;
    }
  public:
    void init (uint32_t size_) {
      //in this implementation, size = m/2 (total capacity/2)
      size = size_; // size must be a power of 2
      falsePositiveRatio = VCACHE_FP;
      dataWidth = TAG_SIZE+LOG_TABLE_SIZE+LOG_MAX_TABLES;
      activeFalsePositiveRatio = 1 - sqrt(1-falsePositiveRatio);
      numHashingFunctions = -log2(activeFalsePositiveRatio);
      activeMaxInserts = (size*ln2/numHashingFunctions);
      printf("VictimCache: m=%d, f=%f, k=%d, n=%d\n", size*2, activeFalsePositiveRatio, numHashingFunctions, activeMaxInserts);
      hashingOutputs = new int32_t [numHashingFunctions];
      active_v[0] = new bool [size];
      active_p[0] = new bool [size];
      active_v[1] = new bool [size];
      active_p[1] = new bool [size];
      active = 0;
    }
    VictimCache (uint32_t size_) {
      init(size_);
    }
    ~VictimCache() {
      delete active_v[0];
      delete active_p[0];
      delete active_v[1];
      delete active_p[1];
      delete hashingOutputs;
    }
    void insert (uint32_t idx, bool p) {
      numInserted[active]++;
      if (isFull()) {
        swap();
        numInserted[active]++;
      }

      for (int32_t i = 0; i < numHashingFunctions; ++i) {
        hashingOutputs[i] = hash(idx, i);

        active_v[active][hashingOutputs[i]] = true;
        active_p[active][hashingOutputs[i]] = p;
      }
    }
    int32_t search(uint32_t idx) {
      idx = idx & (size - 1);

      //data is in active 1
      int32_t return_value = contains(idx, active);
      if (return_value >= 0) 
        return return_value; 

      //data is in active 2
      return_value = contains(idx, next());
      if (return_value >= 0) {
        bool p = (return_value == 1)? true : false;
        //copy data to active 1 (undo aging)
        insert(idx, p);
        return return_value;
      }
      return -1;
    }
    void flush() {
      flush(0);
      flush(1);
    }
    uint32_t getSize() {
      return 4*size; // 2 bits per entry (valid bit and prediction), 2 tables (active1 and active2)
    }
};

class Tage {
  private:
    BimodalTable *bim;
    int32_t predBank, altpredBank, usedBank;
    TaggedTableEntry *pred, *altpred;
    bool predDir, altpredDir, tagePred;
    AltPredTable *altPredTable;
    LinearRandomRegister *randomSeven;

    bool highConfidence;
  protected:
    uint32_t numTables;
    TaggedTable **table;
    History *hist;
    uint32_t *logTableSizes, *tagWidths;
    Stats *stats;
    VictimCache *vcache;
    bool vcacheEnable;
    int32_t victimCacheBank;

    /*Used for getting entries, for prediction or for updating*/
    virtual TaggedTableEntry* getEntry (UINT32 PC, uint32_t historyIdx) {
      return table[historyIdx]->getEntry(hist->getIndex(PC, historyIdx));
    }

    /*Used for getting entries, for prediction or for updating*/
    virtual TaggedTableEntry* getMatch (UINT32 PC, uint32_t historyIdx) {
      TaggedTableEntry *tte = getEntry(PC, historyIdx);
      if (tte->tag == hist->getTag(PC, historyIdx)) return tte;
      return nullptr;
    }

    /*Used for getting newly allocated entries*/
    virtual TaggedTableEntry* getAllocEntry (UINT32 PC, uint32_t historyIdx) {
        return getEntry(PC, historyIdx);
    }

  public:
    Tage(History *h, uint32_t numTables_, uint32_t *logTableSizes_, uint32_t *tagWidths_) {
      //Sizes
      numTables = numTables_;
      logTableSizes = logTableSizes_;
      tagWidths = tagWidths_;

      //allocate structures
      table = new TaggedTable* [numTables+1];
      bim = new BimodalTable();
      hist = h;//new History(numTables, historySizes, tagWidths, logTableSizes);
      altPredTable = new AltPredTable(8);
      randomSeven = new LinearRandomRegister(7);
      for (uint32_t i = 1; i <= numTables; ++i) {
        table[i] = new TaggedTable(i, tagWidths[i], logTableSizes[i]);
      }
      stats = new Stats [hist->numHistories+1];
      highConfidence = false;
      vcache = new VictimCache (256);
      victimCacheBank = 0;
      vcacheEnable = false;
    }
    ~Tage () {
      delete bim;
      delete hist;
      delete altPredTable;
      delete randomSeven;
      for (uint32_t i = 1; i <= numTables; ++i) {
        delete table[i];
      }
      delete [] table;
      delete [] stats;
      delete vcache;
    }
    int32_t getSize() {
      int32_t size = 0;
      size += bim->getSize();
      size += 10; // TICK register
      size += altPredTable->getSize();
      size += randomSeven->getSize();
      for (uint32_t i = 1; i <= numTables; ++i) {
        size += table[i]->getSize();
      }
      size += (stats->getSize() * (hist->numHistories+1));
      size += vcache->getSize();
      return size;
    }

    bool isConfident () {
      return highConfidence;
    }

    bool getPrediction (UINT32 PC) {
      TaggedTableEntry *temp=nullptr;
      predBank = 0; altpredBank = 0;
      pred=nullptr; altpred=nullptr;
      predDir = false; altpredDir = false;
      tagePred = false;
      highConfidence = false;

      //look for pred -- match in longest history table
      for (int i = hist->numHistories; i > 0; --i) {
        temp = getMatch(PC, i);
        if (temp) {
          pred = temp;
          predBank = i;
          break;
        }
      }

      //look for altpred -- match in second longest history table
      for (int i = predBank-1; i > 0; --i) {
        temp = getMatch(PC, i);
        if (temp) {
          altpred = temp;
          altpredBank = i;
          break;
        }
      }

      //set predDir and altPred dir based on matches
      bool predIsUnbiased = false;
      if (!pred) {
        predDir =  bim->getPrediction(PC);
        altpredDir = predDir;
        highConfidence = bim->isConfident();
      } else {
        predDir = pred->counter.getPrediction();
        predIsUnbiased = pred->counter.isUnbiased();

        if (!altpred) {
          altpredDir = bim->getPrediction(PC);
        } else {
          altpredDir = altpred->counter.getPrediction();
        }

        highConfidence = (abs (2 * pred->counter.getValue() + 1) >= (1 << CWIDTH) - 1);
      }

      //decide between pred and altpred
      bool useAltPred = altPredTable->index(PC, predDir, predBank, hist->numHistories);

      /*return prediction
        putting prediction inside tage pred variable so tage can track what IT predicted.. in case
        other predictors are added along side tage.
        */
      tagePred = (!useAltPred || !predIsUnbiased)? predDir : altpredDir;

      // Update prediction counters (only used during reconfiguration)
      usedBank = (!useAltPred || !predIsUnbiased)? predBank : altpredBank;
      if (usedBank != 0) stats[usedBank].predictions += 1;

      return tagePred;
    }

    void update (uint32_t PC, bool resolveDir, bool predictedDir, UINT32 branchTarget) {
      (void)predictedDir; //surpress unused paramter warning
      static int32_t TICK = 0;

      //We want to allocate a new entry if the prediction was wrong.
      bool ALLOC = ((tagePred != resolveDir) && (predBank < (int32_t) hist->numHistories));

      // Update missprediction counters (only used during reconfiguration)
      if (tagePred != resolveDir) stats[usedBank].mispredictions += 1;

      // 1) Update the AltPred table
      if (predBank > 0) {
        /*Only update AltPred table if the counter is weakly biased, i.e. newly allocated.
           The alternate prediction is meant to act as certainty before switching to a newly allocated
           entry. That entry needs to become strongly biased first.
        */
        if (pred->counter.isUnbiased()) {
          /*ALLOC is only true if a misprediction occurred, if the highest matching table was correct,
            that means the AltPred was used and was wrong. Therefore a new correct entry has already been
            allocated and we do not need to allocate a another new one.
          */
          if (predDir == resolveDir)
            ALLOC = false;

          /*If the predictions are different.. update AltPred table to point towards the one that
            was correct */
          if (predDir != altpredDir) {
            //Update the altpred table
            altPredTable->update(PC, predDir, predBank, hist->numHistories, (altpredDir == resolveDir));
          }
        }
      }

      // 2) Allocate new entries, periodic decrement of useful
      if (ALLOC) {
        TaggedTableEntry *temp;
        /* We may want to be more aggressive with how many entries we allocate. This was discused in original
           TAGE paper [Tage, page 7] as one of the differences between the PPM predictor [PPM]. TAGE is much more conservative,
           only allocating one entry on a mispredict, while PPM was more aggressive allocating many entires.
           The insight behind TAGE is that we don't want to cause a big disruption if the mispredict was just a
           random accident.*/
        int numEntriesToAllocate = (BIMOD_ALLOC && predBank == 0)? hist->numHistories : NNN;

        /*These variables are used to avoid the ping-pong phenomenon, explained in Rule B in the paragraph
         below. hopDistance is set to 2 with a low probability to release some contention for a possibly hot
         entry*/
        int hopDistance = 1;
        int rand;
        if ( ( rand = randomSeven->get(hist->getPHist()) ) < 32)
          hopDistance = 2;

        /*The variables are used to estimate when the number of conflicts are too high and thus the
          useful bits need to be decremented
        */
        int numConflict = 0;
        int numAlloc = 0;

        //Check the VictimCache to see if we would have gotten the prediction correct
#if VCACHE_EN
        if (vcacheEnable) {
          int32_t vcache_out = -1;
          if (victimCacheBank != 0) vcache_out = vcache->search(hist->getTag(PC, victimCacheBank)  | (hist->getIndex(PC, victimCacheBank) << TAG_SIZE));
          if (vcache_out >= 0) {
            bool vcachePredDir = (vcache_out == 1)? true : false;
            if (vcachePredDir == resolveDir) {
              //restore the entry if the prediction would have been correct
              temp = getAllocEntry(PC, victimCacheBank);
              if (temp) {
                if (!temp->isUseful()) {
                  //Entry is available, so initialize new entry here.
                  //printf("want to allocate\n");
                  if (temp->isUnused()) {
                    //stats[victimCacheBank].overwriteUnused++;
                    vcache->insert(hist->getTag(PC, victimCacheBank)  | (hist->getIndex(PC, victimCacheBank) << TAG_SIZE), temp->counter.getPrediction());
                  }
                  //printf("new entry\n");
                  temp->set( hist->getTag(PC, victimCacheBank), resolveDir);
                  temp->useful.inc(true);

                  //Bookkeeping
                  //stats[victimCacheBank].allocations++;
                  numAlloc++;
                  numEntriesToAllocate--;
                } else {
                  numConflict++;
                }
              }
            }
          }
        }
#endif

        /*Search through higher ranking tables for: [Tage, page 7]
           A.1) If there exists some entry in a higher table s.t. u=0
           A.2) else the u counters are decremented (specific policy is different here)
           B) Avoiding ping-ponging: if two-tables are available, then the higher order
              one should be selected with a lower probability than the lower order table.
              If say, we were to always pick the lower order table, then we could end up
              in a cycle where two branches are fighting for a single entry (repeatedly
              kicking the other out before it warms up to be useful). Having the probability
              of mapping to a higher table limits the effect of such a scenario. This is
              implemented with the variable "hopDistance".
        */
        for (int i = predBank + hopDistance; i <= NHIST; i += 1) {
          stats[i].attempted++;
          temp = getAllocEntry(PC, i);
          if (!temp) continue;
          if (!temp->isUseful()) {
            //Entry is available, so initialize new entry here.
            //printf("want to allocate\n");
            if (temp->isUnused()) {
              stats[i].overwriteUnused++;
#if VCACHE_EN
              if (vcacheEnable && (i == victimCacheBank)) vcache->insert(hist->getTag(PC, victimCacheBank)  | (hist->getIndex(PC, victimCacheBank) << TAG_SIZE), temp->counter.getPrediction());
#endif
            }
            //printf("new entry\n");
            temp->set( hist->getTag(PC, i), resolveDir);

            //Bookkeeping
            stats[i].allocations++;
            numAlloc++;
            if (numEntriesToAllocate <= 0)
              break;
            i +=  1; // skip a table
            numEntriesToAllocate -= 1;
          } else {
            //Entry is not available, so increase conflict
              numConflict++;
              stats[i].conflicts++;
          }
        }

        /*Formula for determening when to decrease useful bits:
            Intuition: if the number of conflicts is greatly exceeding then number of
            allocations, then the tables are getting too full. Useful entries stay useful
            indefinitely by default, so periodically decrimenting them will force them
            to continually prove their usefulness as time goes on.*/
        TICK += (numConflict - numAlloc);
        if (TICK < 0)
          TICK = 0;
        if (TICK > 1023) {
          for (uint32_t i = 1; i <= numTables; i++) {
            table[i]->decAllUseful();
          }
          TICK = 0;
        }
      }

      // 3) Update the prediction counters
      if (predBank > 0) {
        //Sometimes update AltPred prediction counter
        if (pred->counter.isUnbiased()) {
          //the pred is unbiased, so we are unsure about the prediction. Maybe
          // the case that altpred is better suited. Therefore update altpred.
          if (predDir != resolveDir) {
            if (altpredBank > 0) {
              if (altpred->isUnused()) altpred->useful.reset();
              altpred->counter.inc(resolveDir);
            }
            if (altpredBank == 0)
              bim->update(PC, resolveDir);
          }
        }

        //Always update pred
        if (pred->isUnused()) pred->useful.reset();
        pred->counter.inc(resolveDir);

        //sign changes: no way it can have been useful
        if (pred->counter.isUnbiased())
          pred->useful.reset();
        //Not originally in TAGE, but unbiased altpred entries are also not useful
#if OVERWRITE_KNOB
        if (altpred && altpred->counter.isUnbiased())
          altpred->useful.reset();
#endif

      } else {
        //bimodal table was used for the prediciton
        bim->update(PC, resolveDir);
      }

      // 4) Increment useful entries
      /*     Prediction is only useful if it is different that the altpred.
             If it isnt different, then we might as well use the altpred and
             have this branch take up one less slot in the predictor.
       */
      if (predDir != altpredDir) { //different from alt pred
        if (predDir == resolveDir) { // the correct prediction
          pred->useful.inc(true);
        }
      }

      // 5) Update the history/path information
      hist->updateHistory (PC, OPTYPE_JMP_DIRECT_COND, resolveDir, branchTarget);

    }
    void trackOtherInst (UINT32 PC, OpType opType, bool resolveDir, UINT32 branchTarget) {
      (void)resolveDir; //surpress unused paramter warning
      switch (opType) {
        case OPTYPE_RET_UNCOND:
        case OPTYPE_JMP_DIRECT_UNCOND:
        case OPTYPE_JMP_INDIRECT_UNCOND:
        case OPTYPE_CALL_DIRECT_UNCOND:
        case OPTYPE_CALL_INDIRECT_UNCOND:
          hist->updateHistory (PC, opType, true, branchTarget);
          break;

        default:;
      }
    }
};

/*References
 [PPM]   A PPM-like, tag-based branch predictor. http://www.jilp.org/vol7/v7paper10.pdf
 [Tage]  A case for (partially) TAgged GEometric history length branch prediction.
          http://www.irisa.fr/caps/people/seznec/JILP-COTTAGE.pdf
 [Mudge] Analysis of Branch Prediction via Data Compression.
          http://web.eecs.umich.edu/~tnm/trev_test/papersPDF/1996.10.Analysis%20of%20Branch%20Prediction%20via%20Data%20Compression.pdf
*/

//==========================================================
// LOOP PREDICTOR CODE:
struct LoopTableEntry {     //loop predictor entry
    uint16_t numIterations; //10 bits =====> length = WIDTHNBITERLOOP
    uint16_t currIteration; // 10 bits =====> length = WIDTHNBITERLOOP
    int iterationCounterWidth;

    uint16_t tag; // 10 bits =====> length = LOOPTAG
    int tagWidth;

    uint8_t confidence;  // 4bits
    uint8_t age; // 4 bits
    bool dir;  // 1 bit =====> +sum = 39

    LoopTableEntry() : numIterations(0), currIteration(0), iterationCounterWidth(WIDTHNBITERLOOP),
                       tag(0), tagWidth(LOOPTAG), confidence(0), age(0), dir(false) {}

    static int32_t getSize() {
      return (2 * WIDTHNBITERLOOP + LOOPTAG + 4 + 4 + 1);
    }
};

class LoopPredictor {
  private:
    SaturatingSignedCounter withLoop; // counter to monitor whether or not loop prediction is beneficial
    History* hist;
    LoopTableEntry* loopTable; // loop predictor table
    uint32_t tableSize;
    LinearRandomRegister *randomSeven;

    //variables for the loop predictor
    bool pred; // loop predictor prediction
    int loopIndex;
    int LIB;
    int hitWay; // hitting way in the loop predictor
    int loopTag; // tag on the loop predictor
    bool loopValid; // validity of the loop predictor prediction

    /*loop with more than 2** WIDTHNBITERLOOP iterations are not treated correctly; but who cares :-)*/
    void incIterationCounter( int index ) {
      loopTable[index].currIteration++;
      loopTable[index].currIteration &= ((1 << WIDTHNBITERLOOP) - 1);

      if(loopTable[index].currIteration > loopTable[index].numIterations) {
        //treat like the 1st encounter of the loop
        loopTable[index].confidence = 0;
        loopTable[index].numIterations = 0;
      }
    }

    /*Getting the index to the Loop table: Loop table is a 4-way associative table (so LOGL-2)*/
    int computeLoopIndex(uint32_t PC) {
      return ((PC & ((1 << (LOGL - 2)) - 1)) << 2);
    }
    int computeTag(uint32_t PC) {
      int temp_tag = (PC >> (LOGL - 2)) & ((1 << 2 * LOOPTAG) - 1);
      temp_tag ^= (temp_tag >> LOOPTAG);
      temp_tag = (temp_tag & ((1 << LOOPTAG) - 1));
      return temp_tag;
    }

    void freeEntry(int index) {
      loopTable[index].numIterations = 0;
      loopTable[index].age = 0;
      loopTable[index].confidence = 0;
      loopTable[index].currIteration = 0;
    }

  public:
    LoopPredictor(History* h) : withLoop( SaturatingSignedCounter(7) ) {
      hist = h;
      tableSize =  1 << LOGL;
      loopTable = new LoopTableEntry[ tableSize ];
      randomSeven = new LinearRandomRegister(7);

      withLoop = false;
      loopValid = -1;
    }


    int8_t getWithLoop() {
      return withLoop.getPrediction();
    }
    void updateWithLoop(bool pos_neg) {
      withLoop.inc( pos_neg );
    }
    bool getValid() {
      return loopValid;
    }

    bool getPrediction(uint32_t PC) {
      hitWay = -1;

      loopIndex = computeLoopIndex(PC);
      LIB = ((PC >> (LOGL - 2)) & ((1 << (LOGL - 2)) - 1));
      loopTag = computeTag(PC);

      for(int i = 0; i < 4; i++) {   // searching the 4 ways
        int index = (loopIndex ^ ((LIB >> i) << 2)) + i; //??

        if(loopTable[index].tag == loopTag) {
          hitWay = i;
          loopValid = ((loopTable[index].confidence == CONFLOOP) ||
                       (loopTable[index].confidence * loopTable[index].numIterations > 128));

          if (loopTable[index].currIteration + 1 == loopTable[index].numIterations)
            return ( pred = !(loopTable[index].dir) );
          else
            return ( pred = (loopTable[index].dir) );
        }
      }

      loopValid = false;
      return (pred = false);
    }

    void update (uint32_t PC, bool resolveDir, bool tagePrediction, bool allocate) {
      (void) PC; //for warning, remove if you want to use.

      if(hitWay >= 0) {    // if any entry found for the loop branch in the loop predictor.
        int index = (loopIndex ^ ((LIB >> hitWay) << 2)) + hitWay;

        /*here if the prediction taken from  the loop predictro is wrong, predictor discards the
         *entry. Otherwise, in case entry is giving a different prediction that the Tage, loop
         *predictor tries to keep it in the table, by "age++". Note that in replacement, entry
         *with age == 0 can be allocated.*/
        if(loopValid) {  // Hit in the loop predictor
          if(resolveDir != pred) {
            freeEntry( index );
            /*loopTable[index].numIterations = 0;
            loopTable[index].age = 0;
            loopTable[index].confidence = 0;
            loopTable[index].currIteration = 0;*/
            return;
          }
          else if ((pred != tagePrediction) || ((randomSeven->get(hist->getPHist()) & 7) == 0)) {   
            if (loopTable[index].age < CONFLOOP)
              loopTable[index].age++;
          }
        }

        //warming up the entry in the loop predictor.
        incIterationCounter(index);

        //Note: pred and the direction of the entry in the loopTable are not necessariliy the same.
        if(resolveDir != loopTable[index].dir) {
          if(loopTable[index].currIteration == loopTable[index].numIterations) {
            if(loopTable[index].confidence < CONFLOOP) loopTable[index].confidence++;

            //just do not predict when the loop count is 1 or 2
            if(loopTable[index].numIterations < 3) {
              /*loopTable[index].dir = resolveDir;
              loopTable[index].numIterations = 0;
              loopTable[index].age = 0;
              loopTable[index].confidence = 0;*/
              freeEntry( index ); loopTable[index].dir = resolveDir;
            }
          } else {
            if (loopTable[index].numIterations == 0) {   // first complete nest;
              loopTable[index].confidence = 0;
              loopTable[index].numIterations = loopTable[index].currIteration;
            } else {  //not the same number of iterations as last time: free the entry
              /*loopTable[index].numIterations = 0;
              loopTable[index].confidence = 0;*/
              freeEntry( index );
            }
          }

          loopTable[index].currIteration = 0;
        }
      } else if(allocate) {
        uint32_t random = randomSeven->get(hist->getPHist()) & 3;

        if((randomSeven->get(hist->getPHist()) & 3) == 0) {
          for(int i = 0; i < 4; i++) {
            int way = (random + i) & 3;
            int index = (loopIndex ^ ((LIB >> way) << 2)) + way;
            if(loopTable[index].age == 0) {  //Least recently used
              loopTable[index].dir = !resolveDir; // most of mispredictions are on last iterations
              loopTable[index].tag = loopTag;
              loopTable[index].numIterations = 0;
              loopTable[index].age = 7;
              loopTable[index].confidence = 0;
              loopTable[index].currIteration = 0;
              break;
            } else {
              loopTable[index].age--;
            }
            break;
          }
        }
      }
    }

    static int32_t getSize() {
      return ((1 << LOGL) * LoopTableEntry::getSize());
    }

    ~LoopPredictor() {
      delete [] loopTable;
      delete randomSeven;
    }
};
//END OF LOOP PREDICTOR CODE
//==========================================================

//==========================================================
// STATISTICAL CORRECTOR CODE:
class Statistical_Corrector {
  private:
    //GEHL predictor components
    SaturatingSignedCounter confidenceCounter1, confidenceCounter2, confidenceCounter3; // width =====> CONFWIDTH
    History* hist;

    int LSUM;

    std::vector <SaturatingSignedCounter> Bias;
    std::vector <SaturatingSignedCounter> BiasSK;
    std::vector <int> Pupdatethreshold;

    std::vector <int> Gm = { 27,22, 17, 14 };
    std::vector <std::vector <SaturatingSignedCounter>> GGEHL;

    std::vector <int> Lm = { 11, 6, 3 };
    std::vector <std::vector <SaturatingSignedCounter>> LGEHL;

    std::vector <int> Pm = { 16, 11, 6, 3 };
    std::vector <std::vector <SaturatingSignedCounter>>  PGEHL;

    std::vector <int> Sm = { 16, 11, 6, 3 };
    std::vector <std::vector <SaturatingSignedCounter>> SGEHL;

    //effective local history size + 11: we use IMLIcount + (LH) << 11
    std::vector <int> Tm = { 22, 17, 14 };
    std::vector <std::vector <SaturatingSignedCounter>> TGEHL;

    long long* L_shist;
    long long* S_slhist;
    long long* T_slhist;

    long long HSTACK[ HSTACKSIZE ];
    int pthstack;

    //IMLI related variables
    long long IMLIcount;

#ifdef IMLI
    std::vector <int> Im = { 10 }; // the IMLIcounter is limited to 10 bits
    std::vector <std::vector <SaturatingSignedCounter>> IGEHL;

    long long localoh; //intermediate data to recover the two bits needed from the past outer iteration
    int8_t PIPE[PASTSIZE]; // the PIPE vector
    int8_t ohhisttable[OHHISTTABLESIZE];

    std::vector <int> Fm = { 2 };
    std::vector <std::vector <SaturatingSignedCounter>> FGEHL;
#endif

    int INDBIAS(uint32_t PC, bool pred) {
      return (((PC<<1) + pred) & ((1<<(LOGBIAS+1)) -1));
    }
    int INDBIASSK(uint32_t PC, bool pred) {
      return ((((PC^(PC>>LOGBIAS))<<1) + pred) & ((1<<(LOGBIAS+1)) -1));
    }
    int INDUPD(uint32_t PC) {
      return (PC & ((1 << LOGSIZEUP) - 1));
    }
    int INDTLOCAL(uint32_t PC) {
      return  (((PC ^ (PC >>3))) & (NSECLOCAL-1));  // different hash for the 3rd history
    }
    int INDSLOCAL(uint32_t PC) {
      return (((PC ^ (PC >>5))) & (NSECLOCAL-1));
    }
    int INDLOCAL(uint32_t PC) {
      return (PC & (NLOCAL-1));
    }

    #define GINDEX (((long long) PC) ^ bhist ^ (bhist >> (8 - i)) ^ (bhist >> (16 - 2 * i)) ^ (bhist >> (24 - 3 * i)) ^ (bhist >> (32 - 3 * i)) ^ (bhist >> (40 - 4 * i))) & (size / (1 << (i >= (NBR - 2))) - 1)
    int GEHL_predict(uint32_t PC, long long history, std::vector <int> history_lengths, std::vector <std::vector <SaturatingSignedCounter>>& GEHL_table ) {
      int PERCSUM = 0;
      int NBR = history_lengths.size();
      int size = GEHL_table[ 0 ].size();
      for (int i = 0; i < NBR; i++) {
        long long bhist = history & ((long long) ((1 << history_lengths[ i ]) - 1));
        int16_t ctr = GEHL_table[ i ][ GINDEX ].getValue();
        PERCSUM += (2 * ctr + 1);
      }
      return (PERCSUM);
    }

    void GEHL_update(uint32_t PC, bool taken, long long history, std::vector <int> history_lengths, std::vector <std::vector <SaturatingSignedCounter>>& GEHL_table ) {
      int NBR = history_lengths.size();
      int size = GEHL_table[ 0 ].size();
      for (int i = 0; i < NBR; i++) {
        long long bhist = history & ((long long) ((1 << history_lengths[ i ]) - 1));
        GEHL_table[ i ][ GINDEX ].inc( taken );
      }
    }

	void updateHistories (uint32_t PC, uint16_t brtype, bool resolveDir, uint32_t branchTarget) {
      bool isConditionalBranch = false;
      bool isReturn = false;
      bool isCallDirect = false;
      switch (brtype) {
#ifndef CBP_FOURTEEN
        case OPTYPE_RET_COND:
          isReturn = true;
          isConditionalBranch = true;
          break;
        case OPTYPE_CALL_DIRECT_COND:
          isCallDirect = true;
        case OPTYPE_JMP_DIRECT_COND:
        case OPTYPE_JMP_INDIRECT_COND:
        case OPTYPE_CALL_INDIRECT_COND:
          isConditionalBranch = true;
          break;
        case OPTYPE_RET_UNCOND:
          isReturn = true;
          break;
        case OPTYPE_CALL_DIRECT_UNCOND:
          isCallDirect = true;
        case OPTYPE_JMP_DIRECT_UNCOND:
        case OPTYPE_JMP_INDIRECT_UNCOND:
        case OPTYPE_CALL_INDIRECT_UNCOND:
          break;
#else
        case OPTYPE_CALL_DIRECT:
        case OPTYPE_INDIRECT_BR_CALL:
        case OPTYPE_RET:
        case OPTYPE_BRANCH_UNCOND:
#endif
        default: break;
      }

#ifdef IMLI
      if(isConditionalBranch) {
        if (branchTarget < PC) {
          if (!resolveDir) {
            IMLIcount = 0;
          }
          else { if (IMLIcount < (MAXIMLIcount)) IMLIcount++; }
        }
      }
#ifdef IMLIOH
      if (IMLIcount >= 1) {
        if(isConditionalBranch) {
          if (branchTarget >= PC) {
            PIPE[(PC ^ (PC >> 4)) & (PASTSIZE - 1)] = ohhisttable[(((PC ^ (PC >> 4)) << SHIFTFUTURE) + IMLIcount) & (OHHISTTABLESIZE - 1)];
            ohhisttable[(((PC ^ (PC >> 4)) << SHIFTFUTURE) + IMLIcount) & (OHHISTTABLESIZE - 1)] = resolveDir;
          }
        }
      }
#endif    
#endif

      if(isConditionalBranch) {
        L_shist[INDLOCAL(PC)] = (L_shist[INDLOCAL(PC)] << 1) + resolveDir;
        S_slhist[INDSLOCAL(PC)] = (S_slhist[INDSLOCAL(PC)] << 1) + resolveDir;
        S_slhist[INDSLOCAL(PC)] ^= (PC & 15);
        T_slhist[INDTLOCAL(PC)] = (T_slhist[INDTLOCAL(PC)] << 1) + resolveDir;
      }

      HSTACK[pthstack] = (HSTACK[pthstack] << 1) ^ (branchTarget ^ (branchTarget >> 5) ^ resolveDir);

      if(isReturn)
      pthstack = (pthstack - 1) & (HSTACKSIZE-1);

    if(isCallDirect) {
      int index = (pthstack + 1) & (HSTACKSIZE-1);
      HSTACK[index] = HSTACK[pthstack];
      pthstack = index;
    }
	}

  public:
    Statistical_Corrector(History* h) : confidenceCounter1( SaturatingSignedCounter (CONFWIDTH) ),
                                        confidenceCounter2( SaturatingSignedCounter (CONFWIDTH) ),
                                        confidenceCounter3( SaturatingSignedCounter (CONFWIDTH) ),
                                        hist( h ) {
    Bias = std::vector <SaturatingSignedCounter> ((1 << (LOGBIAS + 1)), SaturatingSignedCounter(PERCWIDTH));
    for (int i = 0; i < (1 << (LOGBIAS + 1)); ++i)
      Bias[ i ] = (i & 1) ? 15 : -16;

    BiasSK = std::vector <SaturatingSignedCounter> ((1 << (LOGBIAS + 1)), SaturatingSignedCounter(PERCWIDTH));
    for (int i = 0; i < (1 << (LOGBIAS + 1)); ++i)
      BiasSK[ i ] = (i & 1) ? 15 : -16;

    Pupdatethreshold = std::vector <int> ((1 << LOGSIZEUP), 35);

    GGEHL = std::vector <std::vector <SaturatingSignedCounter>> (GNB);
    for(int i = 0; i < GNB; ++i)
      GGEHL[ i ] = std::vector <SaturatingSignedCounter> ((1 << LOGGNB), PERCWIDTH - (i < (GNB-1)) );

    LGEHL = std::vector <std::vector <SaturatingSignedCounter>> (LNB);
    for(int i = 0; i < LNB; ++i)
      LGEHL[ i ] = std::vector <SaturatingSignedCounter> ((1 << LOGLNB), PERCWIDTH - (i < (LNB-1)));

    PGEHL = std::vector <std::vector <SaturatingSignedCounter>> (PNB);
    for(int i = 0; i < PNB; ++i)
      PGEHL[ i ] = std::vector <SaturatingSignedCounter> ((1 << LOGPNB), PERCWIDTH - (i < (PNB-1)));

    SGEHL = std::vector <std::vector <SaturatingSignedCounter>> (SNB);
    for(int i = 0; i < SNB; ++i)
      SGEHL[ i ] = std::vector <SaturatingSignedCounter> ((1 << LOGSNB), PERCWIDTH - (i < (SNB-1)));

    TGEHL = std::vector <std::vector <SaturatingSignedCounter>> (TNB);
    for(int i = 0; i < TNB; ++i)
      TGEHL[ i ] = std::vector <SaturatingSignedCounter> ((1 << LOGTNB), PERCWIDTH - (i < (TNB-1)));

    for(int i = 0; i < GNB; ++i) {
      for(int j = 0; j < ((1 << LOGGNB) - 1); ++j) {
        if(!(j & 1)) GGEHL[ i ][ j ].setValue( -1 );
      }
    }

    for (int i = 0; i < LNB; ++i) {
      for (int j = 0; j < ((1 << LOGLNB) - 1); ++j) {
        if (!(j & 1)) LGEHL[ i ][ j ].setValue( -1 );
      }
    }

    for (int i = 0; i < PNB; ++i) {
      for (int j = 0; j < ((1 << LOGPNB) - 1); ++j) {
        if (!(j & 1)) PGEHL[ i ][ j ].setValue( -1 );
      }
    }

    for (int i = 0; i < SNB; i++) {
      for (int j = 0; j < ((1 << LOGSNB) - 1); j++) {
        if (!(j & 1)) SGEHL[ i ][ j ].setValue( -1 );
      }
    }

    for (int i = 0; i < TNB; i++) {
      for (int j = 0; j < ((1 << LOGTNB) - 1); j++) {
        if (!(j & 1)) TGEHL[ i ][ j ].setValue( -1 );
      }
    }

#ifdef IMLI
#ifdef IMLIOH
			FGEHL = std::vector <std::vector <SaturatingSignedCounter>> (FNB);
      for (int i = 0; i < FNB; i++)
        FGEHL[ i ] = std::vector <SaturatingSignedCounter> ((1 << LOGFNB), PERCWIDTH - (i < (FNB-1)));

      for (int i = 0; i < FNB; i++) {
        for (int j = 0; j < ((1 << LOGFNB) - 1); j++) {
          if (!(j & 1)) FGEHL[ i ][ j ].setValue( -1 );
        }
      }
#endif
#ifdef IMLISIC
			IGEHL = std::vector <std::vector <SaturatingSignedCounter>> (INB);
      for (int i = 0; i < INB; i++)
        IGEHL[ i ] = std::vector <SaturatingSignedCounter> ((1 << LOGINB), PERCWIDTH - (i < (INB-1)));

      for (int i = 0; i < INB; i++) {
        for (int j = 0; j < ((1 << LOGINB) - 1); j++) {
          if (!(j & 1)) IGEHL[ i ][ j ].setValue( -1 );
        }
      }
#endif
#endif

    L_shist = new long long [ NLOCAL ];
    for (int i = 0; i < NLOCAL; i++) {
      L_shist[i] = 0;
    }

    S_slhist = new long long [ NSECLOCAL ];
    for (int i = 0; i < NSECLOCAL; i++) {
      S_slhist[i] = 0;
    }

    T_slhist = new long long [ NSECLOCAL ]; 
  }

  bool getPrediction(uint32_t PC, bool tageLoopPrediction, bool tageConfidence) {
    //Compute the SC prediction
    LSUM = 0;

    // Very marginal effect
    // begin to bias the sum towards TAGE predicted direction
    LSUM = 1;
    LSUM += (2 * PNB);

#ifdef LOCALH
    LSUM += (2 * LNB);
#endif

    LSUM += 2 * (SNB + TNB);

#ifdef IMLI
		LSUM += 8;
#endif    

    if( !tageLoopPrediction )
      LSUM = -LSUM;

    //integrate BIAS prediction
    int8_t ctr = (int8_t) Bias[ INDBIAS(PC, tageLoopPrediction) ].getValue();
    LSUM += (2 * ctr + 1);
    ctr = (int8_t) BiasSK[ INDBIASSK(PC, tageLoopPrediction) ].getValue();
    LSUM += (2 * ctr + 1);

#ifdef IMLI
#ifdef IMLIOH
		localoh = 0;
		localoh = PIPE[(PC ^ (PC >> 4)) & (PASTSIZE - 1)] + (localoh << 1);

		for (int i = 0; i >= 0; i--)
			localoh = ohhisttable[(((PC ^ (PC >> 4)) << SHIFTFUTURE) + IMLIcount + i) & (OHHISTTABLESIZE - 1)] + (localoh << 1);

		if (IMLIcount >= 2)
			LSUM += 2 * GEHL_predict((PC << 2), localoh, Fm, FGEHL);
#endif

#ifdef IMLISIC
		LSUM += 2 * GEHL_predict(PC, IMLIcount, Im, IGEHL);
#else
		long long interIMLIcount= IMLIcount;
		IMLIcount = 0;/* just a trick to disable IMLIcount*/
#endif
#endif

#ifndef IMLI
    IMLIcount = 0;
#endif

    LSUM += GEHL_predict((PC << 1) + tageLoopPrediction /*PC*/, (hist->getGlobalHistory() << 11) + IMLIcount, Gm, GGEHL);

#ifdef LOCALH
    LSUM += GEHL_predict(PC, L_shist[ INDLOCAL(PC) ], Lm, LGEHL);
    LSUM += GEHL_predict(PC, (T_slhist[ INDTLOCAL(PC) ] << 11) + IMLIcount, Tm, TGEHL);
    LSUM += GEHL_predict(PC, S_slhist[ INDSLOCAL(PC) ], Sm, SGEHL);
#endif

#ifdef IMLI    
#ifndef IMLISIC
		IMLIcount = interIMLIcount;
#endif
#endif

    LSUM += GEHL_predict(PC, HSTACK[pthstack], Pm, PGEHL);

    bool SCPRED = (LSUM >= 0);
    bool temp_prediction = tageLoopPrediction;

    // choose between the SC output and the TAGE + loop  output
    if(tageLoopPrediction != SCPRED) {
      //Chooser uses TAGE confidence and |LSUM|
      temp_prediction = SCPRED;

      if( tageConfidence ) {
        if( (abs(LSUM) < Pupdatethreshold[INDUPD(PC)] / 3) ) {
          //temp_prediction = confidenceCounter1.getPrediction() ? SCPRED : tageLoopPrediction;
          temp_prediction = confidenceCounter1.getPrediction() ? tageLoopPrediction : SCPRED;
        } else if( (abs(LSUM) < 2 * Pupdatethreshold[INDUPD(PC)] / 3) ) {
          //temp_prediction = confidenceCounter2.getPrediction() ? SCPRED : tageLoopPrediction;
          temp_prediction = confidenceCounter2.getPrediction() ? tageLoopPrediction : SCPRED;
        } else if( (abs(LSUM) < Pupdatethreshold[INDUPD(PC)]) ) {
          //temp_prediction = confidenceCounter3.getPrediction() ? SCPRED : tageLoopPrediction;
          temp_prediction = confidenceCounter3.getPrediction() ? tageLoopPrediction : SCPRED;
        }
      }
    }
    return temp_prediction;
  }

  void update (uint32_t PC, uint16_t brtype, uint32_t branchTarget, bool resolveDir, bool tagePrediction, bool tageConfidence) {
    bool SCPRED = (LSUM >= 0);
    if( tagePrediction != SCPRED ) {
      if ((abs(LSUM) < Pupdatethreshold[INDUPD(PC)])) {
        if (tageConfidence) {
          if ((abs(LSUM) < Pupdatethreshold[INDUPD(PC)] / 3))
            confidenceCounter1.inc( resolveDir == tagePrediction );
          else if ((abs (LSUM) < 2 * Pupdatethreshold[INDUPD(PC)] / 3))
            confidenceCounter2.inc( resolveDir == tagePrediction );
          else if ((abs (LSUM) < Pupdatethreshold[INDUPD(PC)]))
            confidenceCounter3.inc( resolveDir == tagePrediction );
        }
      }
    }

    if((SCPRED != resolveDir) || ((abs(LSUM) < Pupdatethreshold[INDUPD(PC)]))) {
      if(SCPRED != resolveDir)
        Pupdatethreshold[INDUPD(PC)] += 1;
      else
        Pupdatethreshold[INDUPD(PC)] -= 1;

      if(Pupdatethreshold[INDUPD(PC)] >= 256)
        Pupdatethreshold[INDUPD(PC)] = 255;

      if(Pupdatethreshold[INDUPD(PC)] < 0)
        Pupdatethreshold[INDUPD(PC)] = 0;


      Bias[INDBIAS(PC, tagePrediction)].inc( resolveDir );
      BiasSK[INDBIASSK(PC, tagePrediction)].inc( resolveDir );

#ifdef IMLI
#ifndef IMLISIC
			long long interIMLIcount= IMLIcount;
			IMLIcount = 0; /* just a trick to disable IMLIcount*/
#endif
#endif

      GEHL_update((PC << 1) + tagePrediction /*PC*/, resolveDir, (hist->getGlobalHistory() << 11) + IMLIcount, Gm, GGEHL);
      GEHL_update (PC, resolveDir, L_shist[INDLOCAL(PC)], Lm, LGEHL);

      GEHL_update(PC, resolveDir, S_slhist[INDSLOCAL(PC)], Sm, SGEHL);
      GEHL_update(PC, resolveDir, (T_slhist[INDTLOCAL(PC)] << 11) + IMLIcount, Tm, TGEHL);
      GEHL_update(PC, resolveDir, HSTACK[pthstack], Pm, PGEHL);

#ifdef IMLI
#ifdef IMLISIC
			GEHL_update(PC, resolveDir, IMLIcount, Im, IGEHL);
#else
			IMLIcount = interIMLIcount;
#endif

#ifdef IMLIOH
			if(IMLIcount >= 2) GEHL_update((PC << 2), resolveDir, localoh, Fm, FGEHL);
#endif
#endif
    }

		//History Update of the SC + IMLI
		updateHistories(PC, brtype, resolveDir, branchTarget);
  }

    void trackOtherInst (UINT32 PC, OpType opType, bool resolveDir, UINT32 branchTarget) {
      (void)resolveDir; //surpress unused paramter warning
      switch (opType) {
#ifndef CBP_FOURTEEN
        case OPTYPE_RET_UNCOND:
        case OPTYPE_JMP_DIRECT_UNCOND:
        case OPTYPE_JMP_INDIRECT_UNCOND:
        case OPTYPE_CALL_DIRECT_UNCOND:
        case OPTYPE_CALL_INDIRECT_UNCOND:
#else
        case OPTYPE_CALL_DIRECT:
        case OPTYPE_INDIRECT_BR_CALL:
				case OPTYPE_RET:
				case OPTYPE_BRANCH_UNCOND:
#endif
					updateHistories(PC, opType, true, branchTarget);
					break;

				default:;
			}
    }

  int32_t getSize() {
    int inter = 0;

    inter += 16;      //global histories for SC
    inter = 8 * (1 << LOGSIZEUP);  //the update threshold counters
    inter += (PERCWIDTH) * 4 * (1 << (LOGBIAS));
    inter +=
      (GNB - 2) * (1 << (LOGGNB)) * (PERCWIDTH - 1) +
      (1 << (LOGGNB - 1)) * (2 * PERCWIDTH - 1);

    inter +=
      (PNB - 2) * (1 << (LOGPNB)) * (PERCWIDTH - 1) +
      (1 << (LOGPNB - 1)) * (2 * PERCWIDTH - 1);

    inter +=
      (LNB - 2) * (1 << (LOGLNB)) * (PERCWIDTH - 1) +
      (1 << (LOGLNB - 1)) * (2 * PERCWIDTH - 1);
    inter += NLOCAL * Lm[0];

    inter +=
      (SNB - 2) * (1 << (LOGSNB)) * (PERCWIDTH - 1) +
      (1 << (LOGSNB - 1)) * (2 * PERCWIDTH - 1);
    inter +=
      (TNB - 2) * (1 << (LOGTNB)) * (PERCWIDTH - 1) +
      (1 << (LOGTNB - 1)) * (2 * PERCWIDTH - 1);
    inter += HSTACKSIZE * 16;    // the history stack
    inter += logHSTACKSIZE;      // the history stack pointer

    inter += NSECLOCAL * Sm[0];
    inter += NSECLOCAL * (Tm[0] - 11);

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
#endif
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
#endif
#endif // IMLI

    inter += 3 * CONFWIDTH;  //the 3 counters in the choser
    return inter;
  }


  ~Statistical_Corrector() {
    delete [] L_shist;
    delete [] S_slhist;
    delete [] T_slhist;
  }
};
//END OF STATISTICAL CORRECTOR CODE
//==========================================================


class HPSPredictor : public Tage {
  public:
  protected:
    //Maps the historyIdx to the corresponding history index.
    IdxRegister *idxToHist;
    IdxRegister *numTablesInGroup;
    IdxRegister *numUseful;

    IdxRegister *numTablesInGroupConfig[NUM_RECONFIGS+1];
    uint32_t mispredictionsConfig[NUM_RECONFIGS+1];
    IdxRegister *historyScore;
    uint32_t configNum;
    uint32_t currentConfig;
    uint32_t configTick;
    uint32_t avgAttempted = 0; // this maybe should count toward storage
    uint32_t avgConflicts = 0;
    uint32_t avgOverwriteUnused = 0;
    uint32_t hybridConflictEn = true;
    SaturatingUnsignedCounter switchConfig = SaturatingUnsignedCounter(2);

    TaggedTableEntry* getMatch(UINT32 PC, uint32_t historyIdx) {
      //return if there is no storage allocated to this history
      if ( (*numTablesInGroup)[historyIdx] == 0) return nullptr;

      //get the full index from the hashing function... some may be used to select tables
      uint32_t idx = hist->getIndex(PC, historyIdx);
      //index of table in table group (note: All tables are the same size)
      uint32_t tableNum = (idx >> logTableSizes[1]) & ((*numTablesInGroup)[historyIdx]-1);
      idx = idx & ((1 << logTableSizes[1]) - 1);

      //Search for matching entries
      for (uint32_t i = 1,  j = 0; i <= numTables; i++) {
        if ((*idxToHist)[i] == historyIdx) { //search for tables matching in history length
          //looking for the table "j" that our hash functions told us
          if ( j == tableNum) {
            TaggedTableEntry *tte = table[i]->getEntry(idx);
            if (tte->tag == hist->getTag(PC, historyIdx))
              return tte;
            else
              return nullptr;
          }
          j += 1;
        }
      }

      printf("historyIdx: %d, TableNum: %d\n", historyIdx, tableNum);
      assert (0); //should never happen, an entry always should be found
      return nullptr;
    }

    TaggedTableEntry* getAllocEntry(UINT32 PC, uint32_t historyIdx) {
      //return if there is no storage allocated to this history
      if ( (*numTablesInGroup)[historyIdx] == 0) return nullptr;

      //get the full index from the hashing function... some may be used to select tables
      uint32_t idx = hist->getIndex(PC, historyIdx);
      //index of table in table group (note: All tables are the same size)
      uint32_t tableNum = (idx >> logTableSizes[1]) & ((*numTablesInGroup)[historyIdx]-1);
      idx = idx & ((1 << logTableSizes[1]) - 1);

      for (uint32_t i = 1, j = 0; i <= numTables; i++) {
        if ((*idxToHist)[i] == historyIdx) { //search for tables matching in history length
          if (j == tableNum) {
            TaggedTableEntry *tte = table[i]->getEntry(idx);
            return tte;
          }
          j += 1;
        }
      }
      assert (0); //should never happen, an entry always should be found
      return nullptr;
    }

    uint32_t nextPow2 (uint32_t a) {
      while ( ((a-1)&a) != 0 ) {
        a += 1;
      }
      return a; //(a==0)?1:a;
    }

    void updateIdxToHist () {
      uint32_t totalTables = 0;
      uint32_t ntig = 0;
      for (uint32_t i = 1, j = 1; i <= hist->numHistories; ++i) {
        ntig = numTablesInGroup->get(i);
        assert (ntig == nextPow2(ntig)); // all tables must be a power of 2
        totalTables += ntig;
        for (uint32_t k = 0; k < ntig; ++k) {
          (*idxToHist)[j] = i;
          j += 1;
        }
      }
      if (totalTables != NUM_TABLES) {
        printf("totalTables: %d\n", totalTables);
      }
      assert (totalTables == NUM_TABLES);
    }

    void countUseful () {
      for (uint32_t i = 1; i <= hist->numHistories; ++i) {
        (*numUseful)[i] = 0;
      }
      for (uint32_t i = 1; i <= numTables; ++i) {
        for (uint32_t j = 0; j < (1<<LOG_TABLE_SIZE); ++j) {
          if (table[i]->getEntry(j)->isUseful()) {
            (*numUseful)[ (*idxToHist)[i] ] += 1;
          }
        }
      }
    }

    void updateScore () {
      // count the number of useful entries at each history
      countUseful();

      uint32_t predRatio = 0;
      uint32_t maxAttempted = 0;
      int32_t maxAttemptedBank = 0;
      for (uint32_t i = 1; i <= hist->numHistories; ++i) {
        (*historyScore)[i] = 1;
        predRatio = (1000*stats[i].predictions)/SNAPSHOT_PERIOD;
        if (predRatio < 5 && (*numUseful)[i] < 10) {
          //(*historyScore)[i] = 0;
        }
#if ALLOC_POLICY
        if (!HYBRID_POLICY || hybridConflictEn == false) {
          if (stats[i].attempted > avgAttempted)
            (*historyScore)[i]++;
          if (stats[i].attempted > maxAttempted) {
            maxAttempted = stats[i].attempted;
            maxAttemptedBank = i;
          }
        }
#endif
#if CONF_POLICY
        if (!HYBRID_POLICY || hybridConflictEn == true) {
          if (stats[i].conflicts > avgConflicts)
            (*historyScore)[i]++;
        }
#endif
#if UPDATE_OVERWRITE
        if (stats[i].overwriteUnused > avgOverwriteUnused)
          (*historyScore)[i]++;
#endif
      }
#if ALLOC_POLICY
      if (maxAttemptedBank != 0) {
        for (uint32_t i = maxAttemptedBank; i <= hist->numHistories; ++i)
          (*historyScore)[i]++;
      }
#endif
    }

    /*Use this historyScore produced by the scoring algorithm (updateScore) to reconfigure the table.
     *
     * Rules for reconfiguration: There are 5 possible scores.
     *   Score   :  Meaning
     *     0     :  Remove all storage from this history.
     *
     * Function works in two phases: 1) Reclaim storage. 2) Distribute storage.
     * Note that each history must have a power of two storage associated with it when
     * in direct mapped mode. Also note that becase the number of tiles is a power of 2
     * This will always be possible, no matter what order we remove or add tiles in (assuming
     * our initial configuration is valid).
     */
    void reconfigure (uint32_t configNum) {
      uint32_t reallocSize = 0;

      /*Phase 1:
       *  Attempt to reclaim storage from all tables that have a
       *  score of 0 or 1. All storage will be taken from tables with
       *  a score of 0. Tables with a score of 1 will be reduced to
       *  the number of useful entries they house times USEFUL_THRESHOLD.
       */
      for (uint32_t i = 1; i <= hist->numHistories; ++i) {
        switch ((*historyScore)[i]) {
          case 0: { //Remove all storage
            reallocSize += (*numTablesInGroup)[i];
            (*numTablesInGroupConfig[configNum])[i] = 0;
            break;
          }
          /*case 1:*/default: { //Minimize storage.
            uint32_t usefulTemp = (*numUseful)[i];
            usefulTemp = (USEFUL_THRESHOLD * usefulTemp);
            uint32_t usefulSize = nextPow2( (uint32_t) ceil ( (double) usefulTemp / (1<<LOG_TABLE_SIZE) ) );
            usefulSize = (usefulSize > (*numTablesInGroup)[i])? (*numTablesInGroup)[i] : usefulSize;

            (*numTablesInGroupConfig[configNum])[i] = usefulSize;
            reallocSize += ( (*numTablesInGroup)[i] - usefulSize );
            break;
          }
        }
      }

      /*Phase 2:
       * Attempt to double-up tables as much as possible staring from highest priority.
       */
      for (int32_t i = hist->numHistories+4; i >= 0; --i) {
        while (true) {
          // find the table with minimum size, with priority of i
          uint32_t minIdx = 0, minSize = LOG_MAX_TABLES+10; // <- should be +1 when LOG_MAX_TABLES is fixed
          for (uint32_t j = 1; j <= hist->numHistories; ++j) {
            if ((*historyScore)[j] == (uint32_t) i) {
              if ( (*numTablesInGroupConfig[configNum])[j] < minSize ) {
                minIdx = j;
                minSize = (*numTablesInGroupConfig[configNum])[j];
              }
            }
          }
          // if there are no table of priority i, go to next priority
          if (minIdx == 0) break;

          uint32_t incSize = (*numTablesInGroupConfig[configNum])[minIdx];
          incSize = (incSize==0)? 1 : incSize; //if current size is zero, add 1, otherwise double-up
          // if table is already as large as maximum, go to next priority
          if ( incSize == (1<<LOG_MAX_TABLES) ) break;

          //double-up the table if possible
          if ( incSize <= reallocSize ) {
            (*numTablesInGroupConfig[configNum])[minIdx] += incSize;
            reallocSize -= incSize;
          } else {
            break; //if minimum-size table cannot be doubled-up, no other table of this priority can be doubled up
          }
        }
      }
    }

    void printCurrentConfig () {
      printf("\n");
      printf("%14s:\t%d\n", "Config", currentConfig);
      printf("%14s:\t%d\n", "vcacheBank", victimCacheBank);
      printf("%14s:\t%d\n", "Mispredicts", mispredictionsConfig[currentConfig]);
      printf("%14s:\t\t%d\n", "avgAttempted", avgAttempted);
      printf("%14s:\t\t%d\n", "avgConflicts", avgConflicts);
      printf("%14s:\t\t%d\n", "avgUnused", avgOverwriteUnused);
      printf("%14s:\t", "tablesInGroup");
      for (uint32_t i = 1; i <= hist->numHistories; ++i) printf("%6d\t", (*numTablesInGroupConfig[currentConfig])[i]);
      printf("\n");

      printf("%14s:\t", "Score");
      for (uint32_t i = 1; i <= hist->numHistories; ++i) printf("%6d\t", (*historyScore)[i]);
      printf("\n");
      printf("%14s:\t", "NumUseful");
      for (uint32_t i = 1; i <= hist->numHistories; ++i) printf("%6d\t", (*numUseful)[i]);
      printf("\n");
      printf("%14s:\t", "UsefulRatio");
      for (uint32_t i = 1; i <= hist->numHistories; ++i) printf("%6d\t",  100*
          (*numUseful)[i]/(1 + (1 << LOG_TABLE_SIZE) * (*numTablesInGroupConfig[currentConfig])[i]));
      printf("\n");
      printf("%14s:\t", "Mispredicts");
      for (uint32_t i = 1; i <= hist->numHistories; ++i) printf("%6d\t", stats[i].mispredictions);
      printf("\n");
      printf("%14s:\t", "MispredRatio");
      for (uint32_t i = 1; i <= hist->numHistories; ++i) printf("%6d\t",
          100*stats[i].mispredictions/(mispredictionsConfig[currentConfig]+1));
      printf("\n");
      printf("%14s:\t", "Predicts");
      for (uint32_t i = 1; i <= hist->numHistories; ++i) printf("%6d\t", stats[i].predictions);
      printf("\n");
      printf("%14s:\t", "PredRatio");
      for (uint32_t i = 1; i <= hist->numHistories; ++i) printf("%6d\t",
          100*stats[i].predictions/(SNAPSHOT_PERIOD));
      printf("\n");
      printf("%14s:\t", "Attempted");
      for (uint32_t i = 1; i <= hist->numHistories; ++i) printf("%6d\t", stats[i].attempted);
      printf("\n");
      printf("%14s:\t", "Allocations");
      for (uint32_t i = 1; i <= hist->numHistories; ++i) printf("%6d\t", stats[i].allocations);
      printf("\n");
      printf("%14s:\t", "Conflicts");
      for (uint32_t i = 1; i <= hist->numHistories; ++i) printf("%6d\t", stats[i].conflicts);
      printf("\n");
      printf("%14s:\t", "Overwrites");
      for (uint32_t i = 1; i <= hist->numHistories; ++i) printf("%6d\t", stats[i].overwriteUnused);
      printf("\n");

      printf("\n");
      printf("\n");
    }

  public:
    HPSPredictor(History *h, uint32_t numTables_, uint32_t *logTableSizes_, uint32_t *tagWidths_, uint32_t logMaxTableGroup)
      : Tage (h, numTables_, logTableSizes_, tagWidths_) {

      for (uint32_t i = 0; i <= NUM_RECONFIGS; ++i) {
        numTablesInGroupConfig[i] = new IdxRegister(hist->numHistories, logMaxTableGroup);
      }
      historyScore           = new IdxRegister(hist->numHistories, hist->numHistories+3);
      numUseful              = new IdxRegister(hist->numHistories, 1<<(logTableSizes_[0]+logMaxTableGroup));
      idxToHist              = new IdxRegister(numTables, log2(numTables));

      configTick = 0;
      configNum = 0;
      currentConfig = configNum;
      (*numTablesInGroupConfig[0]) = bestConfig;
      numTablesInGroup = numTablesInGroupConfig[currentConfig];
      updateIdxToHist();
    }
    virtual ~HPSPredictor() {
      delete idxToHist;
      delete historyScore;
      delete numUseful;
      for (uint32_t i = 0; i <= NUM_RECONFIGS; ++i) {
        delete numTablesInGroupConfig[i];
      }
    }
    uint32_t getSize() {
      uint32_t size = Tage::getSize();
      size += numUseful->getSize();
      size += idxToHist->getSize();
      size += historyScore->getSize();
      size += log2(SNAPSHOT_PERIOD)*(NUM_RECONFIGS+1); //mispredictionsConfig
      for (uint32_t i = 0; i <= NUM_RECONFIGS; ++i) {
        size += numTablesInGroupConfig[i]->getSize();
      }
      size += log2(SNAPSHOT_PERIOD); //configTick counter
      size += log2(NUM_RECONFIGS+1); //configNum counter
      size += log2(NUM_RECONFIGS+1); //currentConfig counter
      return size;
    }

    bool getPrediction (UINT32 PC) {
      static uint32_t lastConfig = 0;
      //Use the TAGE alogrithm to get a prediction.
      bool result = Tage::getPrediction (PC);
      bool reconfigFlag = CONTINUOUS_RECONFIG || (configNum <= NUM_RECONFIGS);

      /*The configTick counter counts the number of predictions
       * that have elapsed since the last reconfiguration. Once
       * the counter reaches the SNAPSHOT_PERIOD, we will collect
       * some stats and try a new configuration
       */
      configTick += 1;
      configTick = configTick % SNAPSHOT_PERIOD; //wrap counter

      /*Check to see if it is the end of a reconfig period and
       * that there are still more reconfigs to go. We are only doing
       * NUM_RECONFIGS+1 recongis. The last reconfig is reserved to
       * revert back to the most successful reconfig.
       */
      if ((configTick == 0) && reconfigFlag) {
        uint32_t numMisses = 0;
        uint32_t maxAlloc = 0;
        uint32_t maxAllocBank = 0;
        uint32_t maxOverwrite = 0;
        uint32_t maxOverwriteBank = 0;
        uint32_t minConflict = 0;
        uint32_t minConflictBank = 0;
        avgAttempted = 0;
        avgConflicts = 0;
        avgOverwriteUnused = 0;
        /* 1) Accumulate stats from the previous run and clear
         *    the counters for the next run (include the bimodal
         *    mispredicts)*/
        for (uint32_t i = 0; i <= hist->numHistories; ++i) {
          numMisses += stats[i].mispredictions;
          avgAttempted += stats[i].attempted;
          avgConflicts += stats[i].conflicts;
          avgOverwriteUnused += stats[i].overwriteUnused;
        }
        avgAttempted /= hist->numHistories;
        avgConflicts /= hist->numHistories; // may want to account for tables that dont exist
        avgOverwriteUnused /= hist->numHistories;
        for (uint32_t i = 0; i <= hist->numHistories; ++i) {
          if (maxAlloc < stats[i].attempted && stats[i].conflicts < avgConflicts) {
            maxAlloc = stats[i].attempted;
            maxAllocBank = i;
          }
          if (maxOverwrite < stats[i].overwriteUnused && stats[i].conflicts < avgConflicts) {
            maxOverwrite = stats[i].overwriteUnused;
            maxOverwriteBank = i;
          }
          if (minConflict < stats[i].conflicts && stats[i].attempted > avgAttempted) {
            minConflict = stats[i].conflicts;
            minConflictBank = i;
          }
        }
        if (maxAllocBank != 0 && victimCacheBank != maxAllocBank) {
          vcache->flush();
          victimCacheBank = maxAllocBank;
        }
        if (HYBRID_POLICY) {
          if (numMisses < HYBRID_MISSES) hybridConflictEn = false;
          else hybridConflictEn = true; 
        }

        /* 2) Save the total number of mispredictions that happened in
         *    the current configuration. We will need this information
         *    to revert back to the best config.*/
        mispredictionsConfig[ currentConfig ] = numMisses;

        //debug
        countUseful();
#if PRINT_CONFIG
        printf("Finishing Config:\n");
        printCurrentConfig();
        printf("==============================================================\n\n");
#endif
        //end debug

        /* 3) Point to the next configuration*/
        if ( configNum <= NUM_RECONFIGS ) configNum += 1;
        currentConfig = configNum;
 

        /* 4) Try a new reconfiguration based on the stats.
         *    If this is the last reconfiguration, revert back
         *    to the best performatin config.
         *
         *    Note: Config 0 has already been set*/
        if (configNum <= NUM_RECONFIGS) {
          if ( (1000*mispredictionsConfig[currentConfig - 1])/SNAPSHOT_PERIOD >= CONFIG_THRESH_MISPRED) {
            //if (avgAttempted < (1 << LOG_TABLE_SIZE)) return;
            updateScore();
            reconfigure(configNum);
          } else {
            for (uint32_t i = 1; i <= hist->numHistories; ++i) {
              (*numTablesInGroupConfig[currentConfig])[i] = (*numTablesInGroupConfig[currentConfig - 1])[i];
            }
          }	
        } else {
          vcacheEnable = true;
          /*This is the last reconfig, so choose the previous config
           * that had the least number of mispredictions. */
          uint32_t minMisses = SNAPSHOT_PERIOD+1, minMissConfig = 0;
          for (uint32_t i = 0; i <= NUM_RECONFIGS; ++i) {
            if ( mispredictionsConfig[i] < minMisses) {
              minMisses = mispredictionsConfig[i];
              minMissConfig = i;
            }
          }
#if PRINT_CONFIG
          printf("cc:%d, mc:%d\n", lastConfig, minMissConfig);
#endif
          if (lastConfig == minMissConfig)
            switchConfig.inc(true);// higher means switch
          else
            switchConfig.inc(false); // lower means dont switch

          if (!switchConfig.isValid())
            currentConfig = minMissConfig;
          else
            currentConfig = lastConfig;

#if PRINT_CONFIG
          //Print for debug purposes
          printf("\n\nFinal Config: %d\n", currentConfig);
#endif
        }

        //reset counters
        for (uint32_t i = 0; i <= hist->numHistories; ++i) {
          stats[i].reset();
        }

        /* 5) Set configuration vectors*/
        numTablesInGroup = numTablesInGroupConfig[currentConfig];
        updateIdxToHist();

#if PRINT_CONFIG
        printf("Starting Config:\n");
        printCurrentConfig();
#endif
        lastConfig = currentConfig;
      }

      //return prediction from TAGE
      return result;
    }
};

class PREDICTOR {
  private:
    //parameters
    uint32_t numTables;
    uint32_t tableSize;
    uint32_t tagWidth;
    uint32_t logMaxTableGroup;

    //other variables
    uint32_t *hashOutSizes;
    uint32_t *logTableSizes;
    uint32_t *tagWidths;
    History *h;

    HPSPredictor *tage;
    bool tage_prediction;

    #ifdef LOOPPREDICTOR
      LoopPredictor* loop;
      bool loop_prediction;
    #endif

    #ifdef SC
      Statistical_Corrector* sc;
      bool sc_prediction;
    #endif

    bool prediction;

  public:
    PREDICTOR () {
      logMaxTableGroup = LOG_MAX_TABLES;
      numTables = NUM_TABLES;
      tableSize = LOG_TABLE_SIZE;
      tagWidth = TAG_SIZE;

      hashOutSizes = new uint32_t [NHIST+1];
      logTableSizes = new uint32_t [numTables+1];
      tagWidths = new uint32_t [numTables+1];
      for (int i = 1; i <= NHIST; ++i) {
          hashOutSizes[i] = tableSize + logMaxTableGroup;
      }
      for (uint32_t i = 1; i <= numTables; ++i) {
        logTableSizes[i] = tableSize;
        tagWidths[i] = tagWidth;
      }

      h = new History(NHIST, m, tagWidths, hashOutSizes);
      //tage = new HPSPredictor(h, numTables, logTableSizes, TB, logMaxTableGroup);
      tage = new HPSPredictor(h, numTables, logTableSizes, tagWidths, logMaxTableGroup);
      tage_prediction = false;

      int predictor_size = tage->getSize();
      predictor_size += h->getSize();

      #ifdef LOOPPREDICTOR
        loop = new LoopPredictor(h);
        loop_prediction = false;
        predictor_size += loop->getSize();
      #endif

      #ifdef SC
        sc = new Statistical_Corrector(h);
        sc_prediction = false;
        predictor_size += sc->getSize();
      #endif

      prediction = false;

      printf("\nPredictor Size: %d bits, %.2f KB\n", predictor_size, (double)predictor_size/8/1024);
      printf("Breakdowns:\n");
      printf("\tGlobal History: %d bits, %.2f KB\n", h->getSize(), (double)h->getSize()/8/1024);
      printf("\tTage: %d bits, %.2f KB\n", tage->getSize(), (double)tage->getSize()/8/1024);

      #ifdef SC
      printf("\tSC: %d bits, %.2f KB\n", sc->getSize(), (double)sc->getSize()/8/1024);
      #endif

      #ifdef LOOPPREDICTOR
      printf("\tLoop: %d bits, %.2f KB\n", loop->getSize(), (double)loop->getSize()/8/1024);
      #endif
    }
    ~PREDICTOR () {
      delete [] hashOutSizes;
      delete [] logTableSizes;
      delete [] tagWidths;
      delete h;
      delete tage;
      #ifdef LOOPPREDICTOR
        delete loop;
      #endif
      #ifdef SC
        delete sc;
      #endif
    }

    bool GetPrediction (uint32_t PC) {
      prediction = tage_prediction = tage->getPrediction(PC);

      #ifdef LOOPPREDICTOR
        loop_prediction = loop->getPrediction(PC);
        if(loop->getWithLoop() && loop->getValid()) prediction = loop_prediction;
      #endif

      #ifdef SC
        /*Based on the confidence of SC and Tage/Loop, SC may return it's prediction or Tage/Loop*/
        prediction = sc_prediction = sc->getPrediction(PC, prediction, tage->isConfident());
      #endif
      return prediction;
    }

    void UpdatePredictor (uint32_t PC, OpType opType, bool resolveDir, bool predDir, uint32_t branchTarget) {
      (void) opType ; //surpress unused paramter warning for inUnConditionalBranch. Delete the line if inUnConditionalBranch is used

      #ifdef LOOPPREDICTOR
        if(loop->getValid()) {
          if( prediction != loop_prediction ) {
            loop->updateWithLoop( loop_prediction == resolveDir );
          }
        }
        loop->update(PC, resolveDir, tage_prediction, (prediction != resolveDir) );
      #endif

      #ifdef SC
        sc->update(PC, opType, branchTarget, resolveDir, tage_prediction, tage->isConfident());
      #endif

      tage->update(PC, resolveDir, predDir, branchTarget);
    }

    void TrackOtherInst (UINT32 PC, OpType opType, bool branchTaken, UINT32 branchTarget) {
      tage->trackOtherInst (PC, opType, branchTaken, branchTarget);

      #ifdef SC
			sc->trackOtherInst(PC, opType, branchTaken, branchTarget);
      #endif
    }
};
#endif

