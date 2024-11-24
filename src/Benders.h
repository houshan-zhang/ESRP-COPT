#ifndef __BENDERS_HEADER__
#define __BENDERS_HEADER__


#include "Variables.h"
#include "Functions.h"

using namespace std;


/* Build the Benders master problem formulation of the NRP */
void BuildMasterModel(InstStruct *Inst);

/* Solve the NRP under BD framework */
void SolveMasterModel(InstStruct *Inst);

/* Get the right hand side and coefficients of Benders cut */
void GetBendersCut(InstStruct* Inst, int t, double& mincost);

/* Build the Benders subproblem formulation of the NRP */
void BuildSubModel(InstStruct *Inst);

/* Free and close COPT environment */
void CleanBendersModel(InstStruct *Inst);

/* COPT CALLBACK of adding Benders optimality cuts for binary relaxed solution */
int COPT_CALL BinCallback(copt_prob* Lp, void* cb_data, int cbctx, void* cb_handle);

/* COPT CALLBACK of adding Benders optimality cuts for fractional relaxed solution */
int COPT_CALL FracCallback(copt_prob* Lp, void* cb_data, int cbctx, void* cb_handle);


#endif
