#include <Rcpp.h>
// [[Rcpp::depends(dqrng, BH, sitmo)]]
#include <pcg_random.hpp>
#include <dqrng_distribution.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <Rmath.h>
#include<math.h>
#include<map>
using namespace std;
using namespace Rcpp;
using namespace R;

struct RandomFill : public RcppParallel::Worker {
  RcppParallel::RMatrix<double> output;
  uint64_t seed;
  dqrng::normal_distribution dist{0.0, 1.0};

  RandomFill(NumericMatrix output, const uint64_t seed) : output(output), seed(seed) {};

  void operator()(std::size_t begin, std::size_t end) {
    pcg64 rng(seed, end);
    auto gen = std::bind(dist, std::ref(rng));
    for (std::size_t col = begin; col < end; ++col) {
      RcppParallel::RMatrix<double>::Column column = output.column(col);
      std::generate(column.begin(), column.end(), std::ref(gen));
    }
  }
};

// [[Rcpp::export]]
NumericMatrix parallel_random_matrix(const int n, const int m, const uint64_t seed, const int ncores) {
  NumericMatrix res(n, m);
  RandomFill randomFill(res, seed);
  RcppParallel::parallelFor(0, m, randomFill, m/ncores + 1);
  return res;
    
}

//

class SimulantTypeKey {
	/** Key identifying a SIMULANT TYPE.
	*/
public:
	const char *_gender,*_region,*_bc; 
	int _imd,_age;

	SimulantTypeKey(const char *gender,const char *region, const char *bc,int imd,int age):
		_gender(gender),_region(region),_bc(bc),_imd(imd),_age(age) {}
	
	
	SimulantTypeKey(StringVector vGender,StringVector vRegion,StringVector vBc,IntegerVector vImd,IntegerVector vAge,int iRow) {
		_gender= ((String)vGender[iRow]).get_cstring();
		_region= ((String)vRegion[iRow]).get_cstring();
		_bc= ((String)vBc[iRow]).get_cstring();
		_imd= vImd[iRow];
		_age= vAge[iRow];
	}
	
		
	friend bool operator<(const SimulantTypeKey& a,const SimulantTypeKey& b) {
		int iCmpGender= strcmp(a._gender,b._gender);
		if(iCmpGender<0)return true;
		else 
		{
			int iCmpRegion= strcmp(a._region,b._region);
		  int iCmpBc= strcmp(a._bc,b._bc);
		  return iCmpGender==0 && (iCmpRegion<0 || iCmpRegion==0 && (
		    iCmpBc<0 || iCmpBc==0 && (
				a._imd<b._imd || a._imd==b._imd && (
				a._age<b._age ) ) ) );
		}
	}
	friend bool operator==(const SimulantTypeKey& a,const SimulantTypeKey& b) {
		return strcmp(a._gender,b._gender)==0 && strcmp(a._region,b._region)==0&&
		  strcmp(a._bc,b._bc)==0&& a._imd==b._imd && a._age==b._age;
	}
	
	void Print() {
		//Rprintf("gender:%s; region:%s; imd:%d; yob:%d; ageenter:%d\n",_gender,_region,_imd,_yob,_age);
		
	}
	
	
};


class SimulantTypeWithQuantileKey:public SimulantTypeKey {
	/** Key identifying a SIMULANT TYPE and quantile.
	*/
public:
	float _quantile;
	
	SimulantTypeWithQuantileKey(SimulantTypeKey &other,double quantile):
		SimulantTypeKey(other),_quantile(quantile) {}

	SimulantTypeWithQuantileKey(const StringVector vGender,const StringVector vRegion,const StringVector vBc,const IntegerVector vImd,
		const IntegerVector vAge,const NumericVector vQuantile,int iRow):
				SimulantTypeKey(vGender,vRegion,vBc,vImd,vAge,iRow) {
			_quantile= vQuantile[iRow];
		}
		
	friend bool operator<(const SimulantTypeWithQuantileKey& a,const SimulantTypeWithQuantileKey& b) {
		return ((SimulantTypeKey)a)<((SimulantTypeKey)b) || ((SimulantTypeKey)a)==((SimulantTypeKey)b) && a._quantile<b._quantile;
	  //Rprintf("SimulantTypeWithQuantileKey=%c",SimulantTypeKey);
	  
	}
};

/* Linear Interpolation : returns y-position of point to determine */
double lin_interp(
  const double qp, // x-position of point to determine.
  const double q0, const double q1, // set points x-coords.
  const double y0, const double y1 // set points y-coords.
)
{
  return y0 + ((y1-y0)/(q1-q0)) * (qp - q0);
}

List Debug_CreateRetVal(const int n) {
	/** Just create a return value which allows us to examine variables subsequently in R.
	@param n Number of rows in [pop] table.
	@return Test data.frame with columns {"healthy", "incCond", "bmm", "cmm"} giving zero time in each state, and [last_state] column reporting last state.
	*/
	NumericVector healthy(n), incCond(n), bmm(n), cmm(n);
	IntegerVector death(n);
	
	List out = List::create(Named("healthy") = healthy , _["incCond"] = incCond,
	Named("bmm") = bmm , _["cmm"] = cmm, _["last_state"] = death);

	return out;
}

bool FindQuantileBySimulantTypeAndQuantileRow(int simulantIndex,SimulantTypeKey &simulantKey,int rowOffsetFromFirstSimulantTypeToQuantile,map<SimulantTypeWithQuantileKey,
                                              int> &mapSimulantTypeQuantileToRow,int& lowerSimulantTypeToQuantileRow,bool debugMode=false)
	/** Find a 0-based row index in the dfTransitionTimeQuantilesBySimulantType (former lookup) table for the given [SimulantTypeKey] and [rowOffsetFromFirstSimulantTypeToQuantile].
	@param simulantIndex The simulant's 0-based row index.
	@param simulantKey The SimulantTypeKey identifying this simulant.
	@param rowOffsetFromFirstSimulantTypeToQuantile The offset value 0<=rowOffsetFromFirstSimulantTypeToQuantile<=10 from the first row in [dfTransitionTimeQuantilesBySimulantType] which then together identifies the row giving appropriate quantile and transition times.
	@param mapSimulantTypeQuantileToRow A map taking [SimulantTypeWithQuantileKey]s to the appropriate row in the dfTransitionTimeQuantilesBySimulantType (former lookup) table holding their quantity transition times.
	@param lowerSimulantTypeToQuantileRow [out] Gives the 0-based row in dfTransitionTimeQuantilesBySimulantType (former lookup) table for the required [quantile] - so not always the [quantile]=0 row.
  @param debugMode Print some debug information.
	 @return 
	*/
{
	// lowerSimulantTypeToQuantileRow= mapSimulantTypeQuantileToRow[SimulantTypeWithQuantileKey(simulantKey,rowOffsetFromFirstSimulantTypeToQuantile*0.1)];
	map<SimulantTypeWithQuantileKey,int>::iterator iterSimulantTypeWithQuantileKey;
	SimulantTypeWithQuantileKey simulantTypeWithQuantileKey(simulantKey,rowOffsetFromFirstSimulantTypeToQuantile*0.1);
	iterSimulantTypeWithQuantileKey= mapSimulantTypeQuantileToRow.find(simulantTypeWithQuantileKey);
	if(iterSimulantTypeWithQuantileKey==mapSimulantTypeQuantileToRow.end()) {
		Rprintf("Failed finding key for simulant index %d:\n",simulantIndex);
		simulantTypeWithQuantileKey.Print();
		return false;
		
	}
	lowerSimulantTypeToQuantileRow= iterSimulantTypeWithQuantileKey->second; // value (iterSimulantTypeWithQuantileKey->first is key).
	//if(simulantIndex!=-1)Rprintf("Debug - simulant %d, rowOffsetFromFirstSimulantTypeToQuantile %d, lu row %d\n",simulantIndex,rowOffsetFromFirstSimulantTypeToQuantile,lowerSimulantTypeToQuantileRow);

	//Printing some things to try and understand... 
	if(debugMode){
		if(simulantIndex%200000==0)Rprintf("simulantIndex=%d\n", simulantIndex);
		if(simulantIndex%200000==0)Rprintf("simulantKey=%d\n", simulantKey);
	  if(simulantIndex%200000==0)Rprintf("rowOffsetFromFirstSimulantTypeToQuantile=%d\n", rowOffsetFromFirstSimulantTypeToQuantile);
	  if(simulantIndex%200000==0)Rprintf("lowerSimulantTypeToQuantileKey=%d\n", iterSimulantTypeWithQuantileKey->first);
	  if(simulantIndex%200000==0)Rprintf("lowerSimulantTypeToQuantileRow=%d\n", lowerSimulantTypeToQuantileRow);
	  if(simulantIndex%200000==0)Rprintf("simulantTypeWithQuantileKey=%d\n", simulantTypeWithQuantileKey);
	}
	return true;

}



bool TransitionToNextHealthState(int simulantIndex,int currentHealthState,const NumericMatrix& rank,int gender,
	int transitionAliveColIndex,int transitionDeathColIndex,const IntegerVector& rownumStart,
	const NumericVector& quant,double Mcalib,double Fcalib,NumericVector& timeInCurrentState,
	IntegerVector& stateBeforeDeath,const NumericVector& transitionTimesToAlive,
	const NumericVector& transitionTimesToDeath,map<SimulantTypeWithQuantileKey,int> &mapSimulantTypeQuantileToRow,
	SimulantTypeKey &simulantKey,bool finalAliveState,double Mnudge=1.0,bool debugMode=false)
	/** Perform transition to next health-state for this SIMULANT.
	@param simulantIndex This simulant's 0-based row index. 
	@param currentHealthState The simulant's initial health state.
	@param rank
	@param gender
	@param transitionAliveColIndex The relevant column giving the transition-time to the simulant's next (alive) state.
	@param transitionDeathColIndex The relevant column giving the transition-time to the simulant's death state.
	@param rownumStart 
	@param quant
	@param Mcalib
	@param Fcalib
	@param timeInCurrentState
	@param stateBeforeDeath
	@param transitionTimesToAlive A vector listing appropriate transition-times to alive states for all simulants.
	@param transitionTimesToDeath A vector listing appropriate transition-times to death states for all simulants.
	@param mapSimulantTypeQuantileToRow A map taking [SimulantTypeWithQuantileKey]s to the appropriate row in the dfTransitionTimeQuantilesBySimulantType (former lookup) table holding their quantity transition times.
	@param simulantKey The SimulantTypeKey which identifies this simulant.
	@param finalAliveState
	@param Mnudge
	@param debugMode Print some debug information.
	@return bool whether transition has taken simulant to death state.
	*/
{
	const int GenderMale=1,GenderFemale=2,granularity=1/quant[1]; // quant is 0-based; [lu]'s key has no duplicates and sorts quantile lastmost; quant[0]==0, quant[1]=delta, so 1/delta=numSteps
	int lowerSimulantTypeToQuantileAliveRow,lowerSimulantTypeToQuantileDeathRow;
	double rndTransitionAlive, rndTransitionDeath, timeTransitionAlive, timeTransitionDeath;
	
	

	
	// get random numbers for transitions to alive and death states
	if (gender == GenderMale) {
		if(!finalAliveState)rndTransitionAlive= rank(simulantIndex, transitionAliveColIndex) * Mcalib*Mnudge; //random number for tAlive (calibrated)
		rndTransitionDeath= rank(simulantIndex, transitionDeathColIndex) * Mcalib; //random number for tDeath (calibrated)
	} else {
		if(!finalAliveState)rndTransitionAlive= rank(simulantIndex, transitionAliveColIndex) * Fcalib; //random number for tAlive (calibrated)
		rndTransitionDeath= rank(simulantIndex, transitionDeathColIndex) * Fcalib; //random number for tDeath (calibrated)
	}
	if(!finalAliveState && rndTransitionAlive > 1.0) rndTransitionAlive = 0.999999;
	if (rndTransitionDeath > 1.0) rndTransitionDeath = 0.999999;
	//if(debugMode)Rprintf("debug: Simulant:%d, rndTransitionAlive:%f, rndTnsitionDeath:%f\n", simulantIndex, rndTransitionAlive,rndTransitionDeath);
	
	// get floor row numbers for these transitions
	if(!finalAliveState) {
		// rndTransitionAlive/Death is random number 0..1
		// granularity, is num quantiles
		int rowOffsetFromFirstSimulantTypeToQuantile= rndTransitionAlive*granularity; // {min,max}= {0,granularity}.
		// TODO: error, require [result]<granularity, as with rowOffsetFromFirstSimulantTypeToQuantile==granularity, use incorrect row during lowerSimulantTypeToQuantileAliveRow+1 in lin_interp() below.
		//if(rowOffsetFromFirstSimulantTypeToQuantile>9||rowOffsetFromFirstSimulantTypeToQuantile<0) {
			//Rprintf("Error rowOffsetFromFirstSimulantTypeToQuantile ALIVE %d", rowOffsetFromFirstSimulantTypeToQuantile);
			//return 0;
			//}
			if(rowOffsetFromFirstSimulantTypeToQuantile>9) {
			  //Rprintf("Error rowOffsetFromFirstSimulantTypeToQuantile ALIVE %d", rowOffsetFromFirstSimulantTypeToQuantile);
			  return 9;
		}
			if(rowOffsetFromFirstSimulantTypeToQuantile<0) {
			  //Rprintf("Error rowOffsetFromFirstSimulantTypeToQuantile ALIVE %d", rowOffsetFromFirstSimulantTypeToQuantile);
			  return 0;
			}
		if(!FindQuantileBySimulantTypeAndQuantileRow(simulantIndex,simulantKey,rowOffsetFromFirstSimulantTypeToQuantile,mapSimulantTypeQuantileToRow,lowerSimulantTypeToQuantileAliveRow))return 0;

	}
	

	
	int rowOffsetFromFirstSimulantTypeToQuantile= rndTransitionDeath*granularity; // {min,max}= {0,granularity}.
	// TODO: error, require [result]<granularity. See above.
	//if(rowOffsetFromFirstSimulantTypeToQuantile>9||rowOffsetFromFirstSimulantTypeToQuantile<0) {
		//Rprintf("Error rowOffsetFromFirstSimulantTypeToQuantile DEATH %d", rowOffsetFromFirstSimulantTypeToQuantile);
		//return 0;
	//}
	if(rowOffsetFromFirstSimulantTypeToQuantile>9) {
	  //Rprintf("Error rowOffsetFromFirstSimulantTypeToQuantile DEATH %d", rowOffsetFromFirstSimulantTypeToQuantile);
	  return 9;
	}
	if(rowOffsetFromFirstSimulantTypeToQuantile<0) {
	  //Rprintf("Error rowOffsetFromFirstSimulantTypeToQuantile DEATH %d", rowOffsetFromFirstSimulantTypeToQuantile);
	  return 0;
	}
	if(!FindQuantileBySimulantTypeAndQuantileRow(simulantIndex,simulantKey,rowOffsetFromFirstSimulantTypeToQuantile,mapSimulantTypeQuantileToRow,lowerSimulantTypeToQuantileDeathRow)){
	  return 0;
	  Rprintf("some error");
	}
	// interpolate times for these transitions (TODO SE)
	
	
  if(!finalAliveState)timeTransitionAlive= lin_interp(rndTransitionAlive, quant[lowerSimulantTypeToQuantileAliveRow], 
		quant[lowerSimulantTypeToQuantileAliveRow + 1], transitionTimesToAlive[lowerSimulantTypeToQuantileAliveRow], 
		transitionTimesToAlive[lowerSimulantTypeToQuantileAliveRow + 1]);
	timeTransitionDeath= lin_interp(rndTransitionDeath, quant[lowerSimulantTypeToQuantileDeathRow], 
		quant[lowerSimulantTypeToQuantileDeathRow + 1], transitionTimesToDeath[lowerSimulantTypeToQuantileDeathRow], 
		transitionTimesToDeath[lowerSimulantTypeToQuantileDeathRow + 1]);
	
	bool transitionToDeath= finalAliveState || timeTransitionDeath<=timeTransitionAlive;
	if(transitionToDeath)
	{
		//Rprintf("simulant %d, age %d, transit DEATH in %fy, deceased age %d (lu row %d+1; rnd %f)\n",simulantIndex,simulantKey._age,timeTransitionDeath,(int)(simulantKey._age+lround(timeTransitionDeath)),lowerSimulantTypeToQuantileDeathRow,rndTransitionDeath);
		timeInCurrentState[simulantIndex]= timeTransitionDeath;
		stateBeforeDeath[simulantIndex]= currentHealthState;
		//if(debugMode)Rprintf("debug: timeTransitionDeath:%f\n",(float)timeTransitionDeath);
		
	
	}
	else
	{
		long iTimeTransitionAlive= lround(timeTransitionAlive); // lround rounds +-0.5s away from zero).
		//Rprintf("simulant %d, age %d, transit ALIVE in %fy, new age %d (lu row %d+1; rnd %f)\n",simulantIndex,simulantKey._age,timeTransitionAlive,(int)(simulantKey._age+iTimeTransitionAlive),lowerSimulantTypeToQuantileAliveRow,rndTransitionAlive);
		timeInCurrentState[simulantIndex]= timeTransitionAlive;
		simulantKey._age+= iTimeTransitionAlive;
		//if(debugMode)Rprintf("debug: timeTransitionAlive:%f, simulantKey._age:%d\n",(float)timeTransitionAlive,simulantKey._age);
	}
	return transitionToDeath;
	
	
	if(debugMode){
	  if(!finalAliveState)Rprintf("DEBUG ALIVE: Simulant:%d, rnd %f, lowerSimulantTypeToQuantileAliveRow %d, quantLower %f (tranTime %f),  quantUpper %f (tranTime %f), linInter %f, timeTransitionAlive:%f, simulantKey._age:%d\n",
      simulantIndex,
      rndTransitionAlive,
      lowerSimulantTypeToQuantileAliveRow,
      quant[lowerSimulantTypeToQuantileAliveRow],
           transitionTimesToAlive[lowerSimulantTypeToQuantileAliveRow],
                                 quant[lowerSimulantTypeToQuantileAliveRow + 1],
                                      transitionTimesToAlive[lowerSimulantTypeToQuantileAliveRow + 1],
                                                            lin_interp(rndTransitionAlive, quant[lowerSimulantTypeToQuantileAliveRow], 
                                                                       quant[lowerSimulantTypeToQuantileAliveRow + 1], transitionTimesToAlive[lowerSimulantTypeToQuantileAliveRow], 
                                                                                                                                             transitionTimesToAlive[lowerSimulantTypeToQuantileAliveRow + 1]),
                                                                                                                                             (float)timeTransitionAlive, simulantKey._age);
	  if(!finalAliveState)Rprintf("DEBUG DEATH: Simulant:%d, rnd %f, lowerSimulantTypeToQuantileDeathRow %d, quantLower %f (tranTime %f),  quantUpper %f (tranTime %f), linInter %f, timeTransitionDeath:%f, simulantKey._age:%d\n",
      simulantIndex,
      rndTransitionDeath,
      lowerSimulantTypeToQuantileDeathRow,
      quant[lowerSimulantTypeToQuantileDeathRow],
           transitionTimesToDeath[lowerSimulantTypeToQuantileDeathRow],
                                 quant[lowerSimulantTypeToQuantileDeathRow + 1],
                                      transitionTimesToDeath[lowerSimulantTypeToQuantileDeathRow + 1],
                                                            lin_interp(rndTransitionDeath, quant[lowerSimulantTypeToQuantileDeathRow], 
                                                                       quant[lowerSimulantTypeToQuantileDeathRow + 1], transitionTimesToDeath[lowerSimulantTypeToQuantileDeathRow], 
                                                                                                                                             transitionTimesToDeath[lowerSimulantTypeToQuantileDeathRow + 1]),
                                                                                                                                             (float)timeTransitionDeath, simulantKey._age);
	  if(finalAliveState)Rprintf("DEBUG DEATH: Simulant:%d, rnd %f, lowerSimulantTypeToQuantileDeathRow %d, quantLower %f (tranTime %f),  quantUpper %f (tranTime %f), linInter %f, timeTransitionDeath:%f\n",
      simulantIndex,
      rndTransitionDeath,
      lowerSimulantTypeToQuantileDeathRow,
      quant[lowerSimulantTypeToQuantileDeathRow],
           transitionTimesToDeath[lowerSimulantTypeToQuantileDeathRow],
                                 quant[lowerSimulantTypeToQuantileDeathRow + 1],
                                      transitionTimesToDeath[lowerSimulantTypeToQuantileDeathRow + 1],
                                                            lin_interp(rndTransitionDeath, quant[lowerSimulantTypeToQuantileDeathRow], 
                                                                       quant[lowerSimulantTypeToQuantileDeathRow + 1], transitionTimesToDeath[lowerSimulantTypeToQuantileDeathRow], 
                                                                                                                                             transitionTimesToDeath[lowerSimulantTypeToQuantileDeathRow + 1]),
                                                                                                                                             (float)timeTransitionDeath);
	}
	
	
	
}

bool SimulantHasInvalidAge(int simulantIndex,SimulantTypeKey &simulantKey,bool &bSomeSimulantsOverMaxAge,NumericVector& timeInNextState,IntegerVector stateBeforeDeath,int lastAliveState)
	/** Identify invalid simulant ages.
	Ages may be calculated which are greater than those given in the [dfTransitionTimeQuantilesBySimulantType] (formerly lookup) table. In such cases these are marked in the output results.
	@param simulantIndex This simulant's 0-based row index. 
	@param simulantKey The SimulantTypeKey which identifies this simulant.
	@param bSomeSimulantsOverMaxAge [out] boolean indicating that some simulants have invalid keys.
	@param stateBeforeDeath
	@param timeInNextState list of all simulants times in the next ALIVE health state; used to identify an age error on this simulant.
	@return Boolean this simulant had an age error.
	*/
{
	if(simulantKey._age>100) {
		bSomeSimulantsOverMaxAge=true;
		timeInNextState[simulantIndex]=-1;
		stateBeforeDeath[simulantIndex]= lastAliveState;
		return true;
	}
	return false;
}

// [[Rcpp::export]]
List mmsim(const DataFrame& pop, const DataFrame& dfTransitionTimeQuantilesBySimulantType,const NumericMatrix& rank,const double& Mcalib_t12=1,const double& Mcalib_t34=1,const double& Mcalib_t56=1,const double& Mcalib_t7=1,const double& Fcalib_t12=1,const double& Fcalib_t34=1,const double& Fcalib_t56=1,const double& Fcalib_t7=1,const double& Mnudge_t5=1)
	/** Run multi-state model for each SIMULANT using age-specific health-state transition times.
	@param pop Table giving SIMULANT details. Expected columns: ageenter, init_state, gender, imd, region, bc, rownum. See the documentation for further details.
	@param dfTransitionTimeQuantilesBySimulantType Table giving transition time quantiles for each SIMULANT TYPE. Expected columns: gender, imd, region, yob, quantile, t1, t2, t3, t4, t5, t6, t7, ageenter. See the documentation for further details.
	@param rank Matrix giving probabilities for each SIMULANT.
	@param Mcalib_t12 ... Mnudge_t5 calibration values for each transition.
	@return Additional data.frame with columns {"healthy", "incCond", "bmm", "cmm"} giving time in each state, and [last_state] column reporting last state.
	*/
{
	// NumericMatrix out(5 , n); // filled with 0
	const int n= pop.nrows();
	NumericVector timeInHealthyState(n), timeIn1ConditionState(n), timeInMultiConditionBasicState(n), timeInMultiConditionComplexState(n);
	IntegerVector stateBeforeDeath(n);
	const IntegerVector rownumStart= pop["rownum"], init_state= pop["init_state"];
	enum {HealthState_Healthy=1,HealthState_1Condition,HealthState_MultiConditionBasic,HealthState_MultiConditionComplex};
	const IntegerVector vSimulantAgeEnter=pop["ageenter"], vSimulantImd=pop["imd"],  iGender=pop["gender"];
	const StringVector vSimulantRegion=pop["region"], vSimulantBc=pop["bc"], vSimulantGender=pop["gender"];

	// map transition-time quantile rows by simulant-type 
	map<SimulantTypeWithQuantileKey,int> mapSimulantTypeQuantileToRow;
	const int numRows= dfTransitionTimeQuantilesBySimulantType.nrows();
	const StringVector vSimulantTypeGender=dfTransitionTimeQuantilesBySimulantType["gender"], vSimulantTypeRegion=dfTransitionTimeQuantilesBySimulantType["region"], vSimulantTypeBc=dfTransitionTimeQuantilesBySimulantType["bc"];
	const IntegerVector vSimulantTypeImd=dfTransitionTimeQuantilesBySimulantType["imd"],
		vSimulantTypeAgeEnter=dfTransitionTimeQuantilesBySimulantType["ageenter"];
	const NumericVector t1=dfTransitionTimeQuantilesBySimulantType["t1"], t2=dfTransitionTimeQuantilesBySimulantType["t2"],
		t3=dfTransitionTimeQuantilesBySimulantType["t3"], t4=dfTransitionTimeQuantilesBySimulantType["t4"],
		t5=dfTransitionTimeQuantilesBySimulantType["t5"], t6=dfTransitionTimeQuantilesBySimulantType["t6"],
		t7=dfTransitionTimeQuantilesBySimulantType["t7"], vSimulantTypeQuantile=dfTransitionTimeQuantilesBySimulantType["quantile"];
	const NumericVector quant= dfTransitionTimeQuantilesBySimulantType["quantile"];
	for(int iRow=0;iRow<numRows;++iRow) {
		SimulantTypeWithQuantileKey stqk= SimulantTypeWithQuantileKey(vSimulantTypeGender,vSimulantTypeRegion,vSimulantTypeBc,vSimulantTypeImd,vSimulantTypeAgeEnter,vSimulantTypeQuantile,iRow);
		mapSimulantTypeQuantileToRow[stqk]= iRow;
		
		// testing: can we recover the SIMULANT added on the previous line?
		int iGotRow;
		int iQuantile= vSimulantTypeQuantile[iRow]*10;
		if(!FindQuantileBySimulantTypeAndQuantileRow(-1,stqk,iQuantile,mapSimulantTypeQuantileToRow,iGotRow)) {
			Rprintf("Failed finding key just added on row %d\n",iRow);
			return Debug_CreateRetVal(n);
		}
		if(iGotRow!=iRow) {
			Rprintf("Failed matching key value just added on row %d : got row %d\n",iRow,iGotRow);
			return Debug_CreateRetVal(n);
		}
		//Rprintf("simulantTypeWithQuantileKey=%d");
	}
	
	bool bSomeSimulantsOverMaxAge=false;
	for (int simulantIndex=0; simulantIndex<n; simulantIndex++) // consider each simulant
	{
	  int simulantInitialAge= (int)vSimulantAgeEnter[simulantIndex];
		//if(simulantIndex%200000==0)Rprintf("Simulant:%d\n",simulantIndex);
			
		// TODO: no longer need [rownum], so resolve cases were fail to find [lu] row for this [pop] row (and after this, remove below condition).
		// handle [pop].[rownum]==NA issue case. Arises where [pop] RIGHT JOIN [directory] (essentially to [lu]) fails.
		if(IntegerVector::is_na(rownumStart[simulantIndex])) { 
			timeInHealthyState[simulantIndex]= -2; // error marker
			stateBeforeDeath[simulantIndex]= init_state[simulantIndex];
			continue;
		}
		SimulantTypeKey simulantKey= SimulantTypeKey( ((String)vSimulantGender[simulantIndex]).get_cstring(),((String)vSimulantRegion[simulantIndex]).get_cstring(),((String)vSimulantBc[simulantIndex]).get_cstring(),(int)vSimulantImd[simulantIndex],simulantInitialAge);
		
		// consider simulant's life path.
		bool bSimulantDead=false;
		if (init_state[simulantIndex] == HealthState_Healthy)
			bSimulantDead= TransitionToNextHealthState(simulantIndex,HealthState_Healthy,rank,iGender[simulantIndex],
				0,1,rownumStart,quant,Mcalib_t12,Fcalib_t12,timeInHealthyState,stateBeforeDeath,t1,t2,mapSimulantTypeQuantileToRow,simulantKey,false,1.0,true);
				
		// TODO: below horrid SimulantHasInvalidAge() code avoids such issues - but find better solution.
		if(SimulantHasInvalidAge(simulantIndex,simulantKey,bSomeSimulantsOverMaxAge,timeIn1ConditionState,stateBeforeDeath,HealthState_Healthy))continue;
		if(!bSimulantDead) {
			if(init_state[simulantIndex] <= HealthState_1Condition)
				bSimulantDead= TransitionToNextHealthState(simulantIndex,HealthState_1Condition,rank,iGender[simulantIndex],
					2,3,rownumStart,quant,Mcalib_t34,Fcalib_t34,timeIn1ConditionState,stateBeforeDeath,t3,t4,mapSimulantTypeQuantileToRow,simulantKey,false,1.0,true);
			if(SimulantHasInvalidAge(simulantIndex,simulantKey,bSomeSimulantsOverMaxAge,timeInMultiConditionBasicState,stateBeforeDeath,HealthState_1Condition))continue;
		}
		if(!bSimulantDead) {
			if(init_state[simulantIndex] <= HealthState_MultiConditionBasic)
				bSimulantDead= TransitionToNextHealthState(simulantIndex,HealthState_MultiConditionBasic,rank,iGender[simulantIndex],
					4,5,rownumStart,quant,Mcalib_t56,Fcalib_t56,timeInMultiConditionBasicState,stateBeforeDeath,t5,t6,mapSimulantTypeQuantileToRow,simulantKey,false,Mnudge_t5, true);
			if(SimulantHasInvalidAge(simulantIndex,simulantKey,bSomeSimulantsOverMaxAge,timeInMultiConditionBasicState,stateBeforeDeath,HealthState_MultiConditionBasic))continue;
		}
		if(!bSimulantDead && init_state[simulantIndex] <= HealthState_MultiConditionComplex)
			TransitionToNextHealthState(simulantIndex,HealthState_MultiConditionComplex,rank,iGender[simulantIndex],
				-1,6,rownumStart,quant,Mcalib_t7,Fcalib_t7,timeInMultiConditionComplexState,stateBeforeDeath,t7,t7,mapSimulantTypeQuantileToRow,simulantKey,true,1.0,true);
	}
	//if(bSomeSimulantsOverMaxAge)Rprintf("WARNING! some simulants had life paths whose transitions gave ages exceeding those provided in the [lookup] (dfTransitionTimeQuantilesBySimulantType) table. Ages in relevant {IncCond, BMM, or CMM} column will show -1 for these.\n");
	
	// convert death to a factor
	CharacterVector ch = {"Healthy", "IncCond", "BMM", "CMM"};
	stateBeforeDeath.attr("class") = "factor";
	stateBeforeDeath.attr("levels") = ch;
	List out = List::create(Named("healthy") = timeInHealthyState , _["incCond"] = timeIn1ConditionState,
		Named("bmm") = timeInMultiConditionBasicState , _["cmm"] = timeInMultiConditionComplexState, _["last_state"] = stateBeforeDeath);

	return out;

	//Rprintf("simulantTypeWithQuantileKey=%d");
	
}


