#ifndef _svm_H
#define _svm_H

#include <stdlib.h> /* abs */
#include <vector>
#include <cstdio>
#include <cstdlib>
#include "simple_sparse_vec_hash.h"
#include "util.h"

using namespace std;

class svm
{
private:
    /* data */
    vector<double> *alpha;
    vector<simple_sparse_vector> data;
    vector<int> target;
    vector<int> *I_active;

    
public:
    svm(vector<double>* const alpha, vector<int>* const I_active) : alpha(alpha), I_active(I_active) {};
    svm(){};
    ~svm();
    double _decision_func(simple_sparse_vector s, vector<simple_sparse_vector> *data, vector<int> *label);
    vector<double> predict(vector<double> test_set);
    void _fit(vector<double> s, vector<int> t);
    void eval(vector<simple_sparse_vector> *data, vector<int> *target);

};

//   def _decision_function(self, x, I_active):
//         n_samples = self.X.shape[0]
//         ret = 0
//         # for i in range(n_samples):
//         for i in I_active:
//             ret = ret + self.y[i, 0] * self.alpha[i] * (gaussian_kernel_function(self.X[i, :], x) + 1)
//         return ret

double svm::_decision_func(simple_sparse_vector s, vector<simple_sparse_vector> *data, vector<int> *label)
{
    double ret = 0;
    int size = I_active->size();
    int dim = alpha->size();
    vector<double> x(dim, 0.0);
    // this->I_active = I_active;
    s.put(&x);
    for (int i =0;i < size;i++)
    {
        int idx = I_active->at(i);
        vector<double> tt(dim, 0.0);
        simple_sparse_vector v = data->at(idx);
        v.put(&tt);
        ret = ret + (*label)[idx] * (*alpha)[idx] * (_guassian_feature_map(tt, x, 10) + 1);

    }   
    if (ret > 0)
    {
        return 1;
    }
    return -1;
}

svm::~svm()
{
    delete alpha;
    delete I_active;
}

void svm::eval(vector<simple_sparse_vector> *data, vector<int> *target)
{
    double ret = 0;
    int dim = target->size();
    double correct = 0; 
    for (int i =0;i < dim;i++)
    {
        sparse_vector x = data->at(i);
        double pred = _decision_func(x, data, target);
        if ((int)pred == target->at(i))
        {
            correct = correct + 1;
        }
    }   
    double precision = correct / dim;
    printf("Precision score: ------------- %f --------------\n", precision);
}

#endif