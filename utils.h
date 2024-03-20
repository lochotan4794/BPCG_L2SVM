
#ifndef _utils_H
#define _utils_H

#include <stdlib.h> /* abs */
#include <vector>
#include <string>
#include <bits/stdc++.h>
#include <cstdio>
#include <cstdlib>
#include <errno.h>
#include <iostream>
#include <fstream>

using namespace std;

std::ostream &writemap(std::ostream &os, double *map[], int size)
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            os << map[i][j] << " ";
        }
        os << "\n";
    }
    return os;
}
template <typename Iter_T>

double vectorNorm(Iter_T first, Iter_T last)
{
    return sqrt(inner_product(first, last, first, 0.0L));
}

vector<double> operator-(vector<double> x, vector<double> y)
{
    vector<double> ret = vector<double>();
    for (int i = 0; i < x.size(); i++)
    {
        ret.push_back(x[i] - y[i]);
    }
    return ret;
}

double norm(vector<double> x)
{
    double ret = 0;
    for (int i = 0; i < x.size(); i++)
    {
        ret = ret + x[i] * x[i];
    }
    return abs(sqrt(ret));
}

double _guassian_feature_map(vector<double> x, vector<double> y, double sigma)
{

    vector<double> diff = x - y;

    double n2 = norm(diff);

    // float result = exp(-(n2 * n2) / (2*sigma*sigma));

    double r, s = 2.0 * sigma * sigma; 

    double result = exp((-(n2 * n2) / s)) / (M_PI * s); 

    // if (errno == UNDERFLOW) {
    //     printf("exp(%f) overflows\n", n2);
    //     result = n2;
    // }
    return result;
};

double _MEB_kernel(vector<double> x1, int y1, vector<double> x2, int y2, int i, int j, double C, double sigma)
{
    int xi = 0;
    if (i == j)
        xi = 1;
    return (double)y1 * (double)y2 * (_guassian_feature_map(x1, x2, 10) + 1) + (xi / C);
};

void _save_K_file(string file_name, vector<vector<double> > *K)
{
    std::fstream of("Map.txt", std::ios::out | std::ios::app);
    int size = K->size();
    if (of.is_open())
    {
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                of << K->at(i)[j] << " ";
            }
            of << "\n";
        }
    }
    of.close();
}

void _load_K_file(string file_name, vector<vector<double> > *K)
{
    // vec4 *F = new vec4[size * size];
    int size = K->size();
    ifstream fp("Map.txt");
    if (!fp)
    {
        cout << "Error, file couldn't be opened" << endl;
        return;
    }
    for (int row = 0; row < size; row++)
    { // stop loops if nothing to read
        for (int column = 0; column < size; column++)
        {
            fp >> K->at(row)[column];
            if (!fp)
            {
                cout << "Error reading file for element " << row << "," << column << endl;
                return;
            }
        }
    }
}

void _cacl_K_matrix(vector<vector<double> > *s, vector<int> *t, vector<vector<double> > *K)
{
    int size = K->size();
    for (int i = 0; i < size; i++)
    {
        vector<double> q;
        for (int j = 0; j < size; j++)
        {
            q.push_back(_MEB_kernel(s->at(i), t->at(i), s->at(j), t->at(j), i, j, 1.0, 0.2));
        }
        K->push_back(q);
    }
};

void _cacl_z_vector(vector<vector<double> > *K, vector<double> *z)
{
    int size = K->size();
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (i == j)
            {
                z->push_back(K->at(i)[i]);
            }
        }
    }
};

void _precision_score(vector<double> *s, vector<int> *t, vector<int> *pred){

};
#endif