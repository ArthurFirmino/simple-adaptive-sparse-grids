#include "mex.h"

#include <class_handle.hpp>
#include <sasg.h>

#include <array>

#ifndef NDIMENSIONS
#error NDIMENSIONS not defined
#endif
using Grid = SASG::Grid<double, NDIMENSIONS>;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Get the command string
    char cmd[64];
	if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
		mexErrMsgTxt("First input should be a command string less than 64 characters long.");

    // New
    if (!strcmp("new", cmd)) {
        // Check parameters
        if (nlhs != 1)
            mexErrMsgTxt("New: One output expected.");
        if (nrhs != 2)
            mexErrMsgTxt("New: Starting grid level expected.");
        // Return a handle to a new C++ instance
        int n = (int) mxGetScalar(prhs[1]);
        plhs[0] = convertPtr2Mat<Grid>(new Grid(n));
        return;
    }

    // Check there is a second input, which should be the class instance handle
    if (nrhs < 2)
		mexErrMsgTxt("Second input should be a class instance handle.");
    
    // Delete
    if (!strcmp("delete", cmd)) {
        // Destroy the C++ object
        destroyObject<Grid>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs != 0 || nrhs != 2)
            mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
        return;
    }

    // Get the class instance pointer from the second input
    Grid* grid_instance = convertMat2Ptr<Grid>(prhs[1]);

    // Call the various class methods
    // Hierarchize
    if(!strcmp("hierarchize", cmd)) {
        if( !mxIsClass( prhs[2] , "function_handle")) {
            mexErrMsgTxt("Argument is not a function handle.");
        }
        mxArray* rhs[NDIMENSIONS + 1];
        rhs[0] = const_cast<mxArray*>(prhs[2]); 
        auto f = [&rhs](const std::array<double,NDIMENSIONS>& x) {
            for(int i = 0; i < NDIMENSIONS; ++i) {
                rhs[i+1] = mxCreateDoubleScalar(x[i]); 
            }
            mxArray* lhs;
            mexCallMATLAB(1,&lhs,NDIMENSIONS+1,rhs,"feval");
            return *mxGetPr(lhs);
        };
        grid_instance->hierarchize(f);
        return;
    }
    // Evaluate
    if(!strcmp("eval", cmd)) {
        std::array<double, NDIMENSIONS> x;
        for(int i = 0; i < NDIMENSIONS; ++i) {
            x[i] = mxGetScalar(prhs[i+2]);
        }
        double f = grid_instance->evaluate(x);
        plhs[0] = mxCreateDoubleScalar(f);
        return;
    }
    // List Nodes
    if(!strcmp("listnodes", cmd)) {
        int N = grid_instance->getsurpluses().size();
        std::vector<double> v;
        v.reserve(NDIMENSIONS*N);
        for(const auto& p : grid_instance->getsurpluses()) {
            auto x = p.first.x();
            for(int i = 0; i < NDIMENSIONS; ++i) {
                v.push_back(x[i]);
            }
        }
        plhs[0] = mxCreateDoubleMatrix(NDIMENSIONS,N,mxREAL);
        double* out = mxGetPr(plhs[0]);
        memcpy(out, &v[0], v.size()*sizeof(double));
        return;
    } 
    // Unrerefine
    if(!strcmp("unrerefine", cmd)) {
        if( !mxIsClass( prhs[2] , "function_handle")) {
            mexErrMsgTxt("Argument is not a function handle.");
        }
        mxArray* rhs[NDIMENSIONS+1];
        rhs[0] = const_cast<mxArray*>(prhs[2]); 
        auto f = [&rhs](const std::array<double,NDIMENSIONS>& x) {
            for(int i = 0; i < NDIMENSIONS; ++i) {
                rhs[i+1] = mxCreateDoubleScalar(x[i]); 
            }
            mxArray* lhs;
            mexCallMATLAB(1,&lhs,NDIMENSIONS+1,rhs,"feval");
            return *mxGetPr(lhs);
        };
        double eps = mxGetScalar(prhs[3]);
        grid_instance->unrerefine(f, eps);
        return;
    }

    // Got here, so command not recognized
    mexErrMsgTxt("Command not recognized.");
}