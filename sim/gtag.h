
///////////////////////////////////////////////////////////////////////
////  Copyright 2015 Samsung Austin Semiconductor, LLC.                //
/////////////////////////////////////////////////////////////////////////
//
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

#define PHT_CTR_MAX  3
#define PHT_CTR_INIT 2

#define HIST_LEN   11
#define TAG_LEN  6

//NOTE competitors are allowed to change anything in this file include the following two defines
//ver2 #define FILTER_UPDATES_USING_BTB     0     //if 1 then the BTB updates are filtered for a given branch if its marker is not btbDYN
//ver2 #define FILTER_PREDICTIONS_USING_BTB 0     //if 1 then static predictions are made (btbANSF=1 -> not-taken, btbATSF=1->taken, btbDYN->use conditional predictor)

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

class PREDICTOR{

  // The state is defined for Gshare, change for your design

 private:
  UINT32  ghr;           // global history register
  UINT32  *pht;          // pattern history table
  UINT32  *pht2;         //
  UINT32  tagLength;     // tag length
  UINT32  historyLength; // history length
  UINT32  numPhtEntries; // entries in pht
  UINT32  *tags;
  
 public:


  PREDICTOR(void);
  //NOTE contestants are NOT allowed to use these versions of the functions
  //ver2   bool    GetPrediction(UINT64 PC, bool btbANSF, bool btbATSF, bool btbDYN);  
  //ver2   void    UpdatePredictor(UINT64 PC, OpType opType, bool resolveDir, bool predDir, UINT64 branchTarget, bool btbANSF, bool btbATSF, bool btbDYN);
  //ver2   void    TrackOtherInst(UINT64 PC, OpType opType, bool branchDir, UINT64 branchTarget);

  // The interface to the functions below CAN NOT be changed
  bool    GetPrediction(UINT64 PC);  
  void    UpdatePredictor(UINT64 PC, OpType opType, bool resolveDir, bool predDir, UINT64 branchTarget);
  void    TrackOtherInst(UINT64 PC, OpType opType, bool branchDir, UINT64 branchTarget);

  // Contestants can define their own functions below
};

/////////////// STORAGE BUDGET JUSTIFICATION ////////////////
// Total storage budget: 32KB + 17 bits
// Total PHT counters: 2^17 
// Total PHT size = 2^17 * 2 bits/counter = 2^18 bits = 32KB
// GHR size: 17 bits
// Total Size = PHT size + GHR size
/////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PREDICTOR::PREDICTOR(void){

  historyLength    = HIST_LEN;
  tagLength        = (1 << TAG_LEN);
  ghr              = 0;
  numPhtEntries    = (1<< HIST_LEN);

  pht = new UINT32[numPhtEntries]; // Secondary Predictor 2-bit counters
  pht2 = new UINT32[numPhtEntries]; // Default Predictor 2-bit counters
  tags = new UINT32[numPhtEntries]; // Tags

  for(UINT32 ii=0; ii< numPhtEntries; ii++){
    pht2[ii]=PHT_CTR_INIT; // Initialize default predictor 2-bit counters 
  }
  
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

//ver2 version of infrastructure. Contestants are NOT allowed to use this version of the function
//ver2 bool   PREDICTOR::GetPrediction(UINT64 PC, bool btbANSF, bool btbATSF, bool btbDYN){
//ver2 
//ver2   UINT32 phtIndex   = (PC^ghr) % (numPhtEntries);
//ver2   UINT32 phtCounter = pht[phtIndex];
//ver2 
//ver2 //  printf(" ghr: %x index: %x counter: %d prediction: %d\n", ghr, phtIndex, phtCounter, phtCounter > PHT_CTR_MAX/2);
//ver2 
//ver2   //NOTE taking advantage of btbANFS and btbATSF to make static predictions based upon previously observed behaviour. You don't have to take advantage of this if you don't want to.
//ver2   if(btbANSF && FILTER_PREDICTIONS_USING_BTB) {
//ver2     return NOT_TAKEN; //statically predict N if it has always been N since discovery if BTB prediction filtering is enabled
//ver2   }
//ver2   else if(btbATSF && FILTER_PREDICTIONS_USING_BTB) {
//ver2     return TAKEN; //statically predict T if it has always been T since discovery if BTB prediction filtering is enabled
//ver2   }
//ver2   else if (btbDYN || !FILTER_PREDICTIONS_USING_BTB) { //use PHT if it is marked as dynamic in the btb
//ver2     if(phtCounter > (PHT_CTR_MAX/2)){ 
//ver2       return TAKEN; 
//ver2     }
//ver2     else{
//ver2       return NOT_TAKEN; 
//ver2     }
//ver2   }
//ver2 
//ver2   //statically predict N if is not marked at all in the btb
//ver2   return NOT_TAKEN; 
//ver2 }

//NOTE contestants are not allowed to use the btb* information from ver2 of the simulation infrastructure. Interface to this function CAN NOT be changed.
bool   PREDICTOR::GetPrediction(UINT64 PC){

  UINT32 phtIndex   = (PC^ghr) % (numPhtEntries); // Secondary Predictor Index - Gshare1
  UINT32 phtIndex2 = (PC) % (numPhtEntries); // Default Predictor - Bimodal
  UINT32 phtCounter = pht[phtIndex]; // Secondary Predictor Counter
  UINT32 phtCounter2 = pht2[phtIndex2]; // Default Predictor Counter
  UINT32 tag = (PC) % (tagLength);
  
//  printf(" ghr: %x index: %x counter: %d prediction: %d\n", ghr, phtIndex, phtCounter, phtCounter > PHT_CTR_MAX/2);
  if(tag == tags[phtIndex]) { // If the Gshare3 Tag matches the tag in the Secondary Predictor at Gshare1 index
    if(phtCounter > (PHT_CTR_MAX/2)){ // Produce appropriate prediction from 2-bc in secondary predictor
      return TAKEN;
    } else {
      return NOT_TAKEN;
    }
  } else { // No tag hit
    if(phtCounter2 > (PHT_CTR_MAX/2)){ // Use appropriate default prediction from appropriate counter
      return TAKEN;
    } else {
      return NOT_TAKEN; 
    }
  }
}


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

//ver2 version of infrastructure. Contestants are NOT allowed to use this version of the function
//ver2 void  PREDICTOR::UpdatePredictor(UINT64 PC, OpType opType, bool resolveDir, bool predDir, UINT64 branchTarget, bool btbANSF, bool btbATSF, bool btbDYN){
//ver2 
//ver2   UINT32 phtIndex   = (PC^ghr) % (numPhtEntries);
//ver2   UINT32 phtCounter = pht[phtIndex];
//ver2 
//ver2   //NOTE only updating PHT and GHR if branch has exhibited dynamic behaviour in the past OR has always been taken until the current not-taken. This is optional, just showing what you could do with this info
//ver2 
//ver2   if(btbDYN || (resolveDir != predDir) || !FILTER_UPDATES_USING_BTB) { 
//ver2   // update the PHT if btb indicates dynamic behaviour OR the outcome does not match the prediction
//ver2     if(resolveDir == TAKEN){
//ver2       pht[phtIndex] = SatIncrement(phtCounter, PHT_CTR_MAX);
//ver2     }else{
//ver2       pht[phtIndex] = SatDecrement(phtCounter);
//ver2     }
//ver2 
//ver2   // update the GHR
//ver2     ghr = (ghr << 1);
//ver2 
//ver2     if(resolveDir == TAKEN){
//ver2       ghr++; 
//ver2     }
//ver2   }
//ver2 
//ver2 }

//NOTE contestants are not allowed to use the btb* information from ver2 of the simulation infrastructure. Interface to this function CAN NOT be changed.
void  PREDICTOR::UpdatePredictor(UINT64 PC, OpType opType, bool resolveDir, bool predDir, UINT64 branchTarget){

  UINT32 phtIndex   = (PC^ghr) % (numPhtEntries); // Secondary Predictor Index - Gshare1
  UINT32 phtIndex2 = (PC) % (numPhtEntries); // Default Predictor - Bimodal
  UINT32 phtCounter = pht[phtIndex]; // Secondary Predictor Counter
  UINT32 phtCounter2 = pht2[phtIndex2]; // Default Predictor Counter
  UINT32 tag = (PC) % (tagLength); // Tag - Gshare3
  
  if(tag == tags[phtIndex]) { // If Secondary Predictor was used
    if(resolveDir == TAKEN) { // If Branch Outcome was taken
      pht[phtIndex] = SatIncrement(phtCounter, PHT_CTR_MAX); // Increment appropriate counter in secondary predictor
      if((phtCounter <= (PHT_CTR_MAX/2)) & (phtCounter2 <= (PHT_CTR_MAX/2))) { // If both secondary and default incorrectly predicted NOT TAKEN
	pht2[phtIndex2] = SatIncrement(phtCounter2, PHT_CTR_MAX); // Increment appropriate counter in default predictor
      } // else dont update default predictor
    }else { // If Branch Outcome was not taken
      pht[phtIndex] = SatDecrement(phtCounter); // Decrement appropriate counter in secondary predictor
      if((phtCounter > (PHT_CTR_MAX/2)) & (phtCounter2 > (PHT_CTR_MAX/2))) { // If both secondary and default incorrectly predicted TAKEN
	pht2[phtIndex2] = SatDecrement(phtCounter2); // Decrement appropriate counter in default predictor
      } // else dont update default predictor
    }
  } else { // If Default Predictor was used
    if(resolveDir == TAKEN) { // If Branch Outcome was taken
      if(phtCounter2 > (PHT_CTR_MAX/2)) { // Update appropriate counter if correctly predicted TAKEN
	pht2[phtIndex2] = SatIncrement(phtCounter2, PHT_CTR_MAX);
      } else { // Otherwise add tag to secondary predictor and initialize two bit counter to WEAKLY TAKEN
        tags[phtIndex] = tag;
	pht[phtIndex] = 2;
      }
    }else { // If Branch outcome was NOT TAKEN
      if(phtCounter2 > (PHT_CTR_MAX/2)) { // Update appropraite default predicor counter if correctly predicted NOT TAKEN
	pht2[phtIndex2] = SatDecrement(phtCounter2);
      } else { // Otherwise add tag to secondary predictor and initialize two bit counter to WEAKLY NOT TAKEN
	tags[phtIndex] = tag;
	pht[phtIndex] = 1;
      }
    }
  }
  
  // update the GHR
  ghr = (ghr << 1);
  
  if(resolveDir == TAKEN){
    ghr++;
  }
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

void    PREDICTOR::TrackOtherInst(UINT64 PC, OpType opType, bool branchDir, UINT64 branchTarget){

  // This function is called for instructions which are not
  // conditional branches, just in case someone decides to design
  // a predictor that uses information from such instructions.
  // We expect most contestants to leave this function untouched.

  return;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


/***********************************************************/
#endif

