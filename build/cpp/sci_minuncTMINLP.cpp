// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// Author: Harpreet Singh, Pranav Deshpande and Akshay Miterani
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

#define __USE_DEPRECATED_STACK_FUNCTIONS__
#include "minuncTMINLP.hpp"
//#include "IpIpoptData.hpp"
#include "sci_iofunc.hpp"

extern "C"
{
#include "call_scilab.h"
#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <localization.h>
#include <sciprint.h>
#include <string.h>
#include <assert.h>
#include <iostream>
}

using namespace std;
using namespace Ipopt;
using namespace Bonmin;

minuncTMINLP::~minuncTMINLP()
{
	free(finalX_);
	free(finalGradient_);
	free(finalHessian_);
}

// Set the type of every variable - CONTINUOUS or INTEGER
bool minuncTMINLP::get_variables_types(Index n, VariableType* var_types)
{
  n = numVars_;
  for(int i=0; i < n; i++)
    var_types[i] = CONTINUOUS;
  for(int i=0 ; i < intconSize_ ; ++i)
  	var_types[(int)(intcon_[i]-1)] = INTEGER;
  return true;
}

// The linearity of the variables - LINEAR or NON_LINEAR
bool minuncTMINLP::get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types)
{
	/*
	n = numVars_;
	for(int i = 0; i < n; i++)
		var_types[i] = Ipopt::TNLP::LINEAR;
  */
  return true;
}

// The linearity of the constraints - LINEAR or NON_LINEAR
bool minuncTMINLP::get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types)
{
 /* m = numConstr_;
  for(int i = 0; i < m; i++)
    const_types[i] = Ipopt::TNLP::LINEAR;
  */
  return true;
}

//get NLP info such as number of variables,constraints,no.of elements in jacobian and hessian to allocate memory
bool minuncTMINLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
{
	finalGradient_ = (double*)malloc(sizeof(double) * numVars_ * 1);
	finalHessian_ = (double*)malloc(sizeof(double) * numVars_ * numVars_);
	n=numVars_; // Number of variables
	m=numConstr_; // Number of constraints
	nnz_jac_g = 0; // No. of elements in Jacobian of constraints 
	nnz_h_lag = n*(n+1)/2; // No. of elements in lower traingle of Hessian of the Lagrangian.
	index_style=TNLP::C_STYLE; // Index style of matrices
	return true;
}

//get variable and constraint bound info
bool minuncTMINLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
	unsigned int i;
	for(i=0;i<n;i++)
	{
		x_l[i]=-1.0e19;
		x_u[i]=1.0e19;
	}

        g_l=NULL;
        g_u=NULL;
	return true;
}

// return the value of the constraints: g(x)
bool minuncTMINLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  	// return the value of the constraints: g(x)
  	g=NULL;
  	return true;
}

// return the structure or values of the jacobian
bool minuncTMINLP::eval_jac_g(Index n, const Number* x, bool new_x,Index m, Index nele_jac, Index* iRow, Index *jCol,Number* values)
{
 	if (values == NULL) 
 	{
    		// return the structure of the jacobian of the constraints
    		iRow=NULL; 
    		jCol=NULL;
  	}
  	else 
	{
		values=NULL;
  	}

  	return true;
}

//get value of objective function at vector x
bool minuncTMINLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  	double check;
  	int* funptr=NULL;  
  	if(getFunctionFromScilab(1,&funptr))
  	{
		return 1;
 	}
  	char name[20]="f";
  	double obj=0;
  	double *xNew=x;
  	createMatrixOfDouble(pvApiCtx, 7, 1, numVars_, xNew);
  	int positionFirstElementOnStackForScilabFunction = 7;
  	int numberOfRhsOnScilabFunction = 1;
  	int numberOfLhsOnScilabFunction = 2;
  	int pointerOnScilabFunction     = *funptr;
  
  	C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
                                                               &numberOfLhsOnScilabFunction,
                                                               &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));
                               
  	if(getDoubleFromScilab(8,&check))
  	{
		return true;
	}
	if (check==1)
	{
		return true;
	}	
	else
	{ 
		if(getDoubleFromScilab(7,&obj))
  		{
			sciprint("No obj value");
			return 1;
  		}
  		obj_value=obj;  
	}
  	return true;
}

//get value of gradient of objective function at vector x.
bool minuncTMINLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  	double check;
  	if (flag1_==0)
  	{	
  		int* gradhessptr=NULL;
  		if(getFunctionFromScilab(2,&gradhessptr))
  		{
			return 1;
  		}  
  		double *xNew=x;
  		double t=1;
  		createMatrixOfDouble(pvApiCtx, 7, 1, numVars_, xNew);
  		createScalarDouble(pvApiCtx, 8,t);
  		int positionFirstElementOnStackForScilabFunction = 7;
  		int numberOfRhsOnScilabFunction = 2;
  		int numberOfLhsOnScilabFunction = 2;
  		int pointerOnScilabFunction     = *gradhessptr;
		char name[20]="gradhess";
 
  		C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
                                                               &numberOfLhsOnScilabFunction,
                                                               &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));
  	}

  	else if (flag1_==1)
  	{
  		int* gradptr=NULL;
  		if(getFunctionFromScilab(4,&gradptr))
  		{
			return 1;
  		}  
  		double *xNew=x;
	  	createMatrixOfDouble(pvApiCtx, 7, 1, numVars_, xNew);
  		int positionFirstElementOnStackForScilabFunction = 7;
  		int numberOfRhsOnScilabFunction = 1;
  		int numberOfLhsOnScilabFunction = 2;
  		int pointerOnScilabFunction     = *gradptr;
		char name[20]="fGrad1";
 	
  		C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
        	                                                       &numberOfLhsOnScilabFunction,
        	                                                       &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));
   	}

	if(getDoubleFromScilab(8,&check))
  	{
		return true;
	}
	if (check==1)
	{
		return true;
	}	
	else
	{ 
		double* resg;  
  		int x0_rows,x0_cols;                           
  		if(getDoubleMatrixFromScilab(7, &x0_rows, &x0_cols, &resg))
  		{
			sciprint("No results");
			return 1;
			
  		}
	
  		Index i;
  		for(i=0;i<numVars_;i++)
  		{
			grad_f[i]=resg[i];
    	    	finalGradient_[i]=resg[i];
  		}
	}
  	return true;
}

// This method sets initial values for required vectors . For now we are assuming 0 to all values. 
bool minuncTMINLP::get_starting_point(Index n, bool init_x, Number* x,bool init_z, Number* z_L, Number* z_U,Index m, bool init_lambda,Number* lambda)
{
 	assert(init_x == true);
  	assert(init_z == false);
  	assert(init_lambda == false);
	if (init_x == true)
	{ //we need to set initial values for vector x
		for (Index var=0;var<n;var++)
			x[var]=varGuess_[var];//initialize with 0 or we can change.
	}

	return true;
}

/*
 * Return either the sparsity structure of the Hessian of the Lagrangian, 
 * or the values of the Hessian of the Lagrangian  for the given values for
 * x,lambda,obj_factor.
*/

bool minuncTMINLP::eval_h(Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,bool new_lambda, Index nele_hess, Index* iRow,Index* jCol, Number* values)
{
	double check;
	if (values==NULL)
	{
		Index idx=0;
		for (Index row = 0; row < numVars_; row++) 
		{
			for (Index col = 0; col <= row; col++)
			{
				iRow[idx] = row;
				jCol[idx] = col;
				idx++;
		  	}
		}
	}

	else 
	{	
		if(flag2_==0)
	  	{
			int* gradhessptr=NULL;
			if(getFunctionFromScilab(2,&gradhessptr))
			{
				return 1;
			}          	
			double *xNew=x;
  			double t=2;
			createMatrixOfDouble(pvApiCtx, 7, 1, numVars_, xNew);
  			createScalarDouble(pvApiCtx, 8,t);
  			int positionFirstElementOnStackForScilabFunction = 7;
  			int numberOfRhsOnScilabFunction = 2;
  			int numberOfLhsOnScilabFunction = 2;
  			int pointerOnScilabFunction     = *gradhessptr;
			char name[20]="gradhess";
  	
  			C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
                	                                               &numberOfLhsOnScilabFunction,
                	                                               &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));
 	    }	

 	    else if (flag2_==1)
 	    {		
			int* hessptr=NULL;
			if(getFunctionFromScilab(6,&hessptr))
			{
				return 1;
			}          	
			double *xNew=x;	
  			createMatrixOfDouble(pvApiCtx, 7, 1, numVars_, xNew);
  			int positionFirstElementOnStackForScilabFunction = 7;
  			int numberOfRhsOnScilabFunction = 1;
  			int numberOfLhsOnScilabFunction = 2;
  			int pointerOnScilabFunction     = *hessptr;
			char name[20]="fHess1";
  	
  			C2F(scistring)(&positionFirstElementOnStackForScilabFunction,name,
                	                                               &numberOfLhsOnScilabFunction,
                	                                               &numberOfRhsOnScilabFunction,(unsigned long)strlen(name));	
 	    }	

		if(getDoubleFromScilab(8,&check))
  		{
			return true;
		}
		if (check==1)
		{
			return true;
		}	
		else
		{ 
	        double* resh;  
  			int x0_rows,x0_cols;                           
  			if(getDoubleMatrixFromScilab(7, &x0_rows, &x0_cols, &resh))
			{
				sciprint("No results");
				return 1;
			}
			
			Index index=0;
			for (Index row=0;row < numVars_ ;++row)
			{
				for (Index col=0; col <= row; ++col)
				{
					values[index++]=obj_factor*(resh[numVars_*row+col]);
				}
			}
	
			Index i;
			for(i=0;i<numVars_*numVars_;i++)
			{
       				finalHessian_[i]=resh[i];
			}
		}	
	}	

       	return true;
}


void minuncTMINLP::finalize_solution(SolverReturn status,Index n, const Number* x, Number obj_value)
{
	finalX_ = (double*)malloc(sizeof(double) * numVars_ * 1);
	for (Index i=0; i<numVars_; i++) 
	{
    		 finalX_[i] = x[i];
	}

	finalObjVal_ = obj_value;
	status_ = status;
}

const double * minuncTMINLP::getX()
{	
	return finalX_;
}

double minuncTMINLP::getObjVal()
{	
	return finalObjVal_;
}

int minuncTMINLP::returnStatus()
{	
	return status_;
}
