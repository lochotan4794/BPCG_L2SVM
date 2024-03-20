#ifndef _BPCG_H
#define _BPCG_H

#include <vector>
#include <stdlib.h> /* abs */
#include <numeric>
#include "simple_sparse_vec_hash.h"
#include <algorithm>

using namespace std;
using sparse_vector = simple_sparse_vector;

inline void logg(string s)
{
    cout << s << endl;
}

struct vertex_info
{
    /* data */
    sparse_vector vertex;
    int index;
};

struct bpcg_result
{
    /* data */
    sparse_vector x;
    vector<double> alpha;
    vector<int> I_active;
    vector<double> f_vls;
};

class bpcg
{
private:
    /* data */
    vector<double> alpha;
    vector<double> f_vls;
    vector<int> I_active;
    vector<sparse_vector> s_list;
    sparse_vector init_vertex;
    vector<vector<double> > K;
    vector<double> z;
    int dim;

public:
    bpcg(int dim, vector<vector<double> > K, vector<double> z);
    ~bpcg();
    sparse_vector p_list;
    sparse_vector dual_list;
    double primal(sparse_vector x);
    sparse_vector vector_2_spare(vector<double> x);
    vector<double> spare_2_vector(sparse_vector x);
    sparse_vector gradient_func(sparse_vector x);
    bpcg_result optimize(int max_iter);
    vertex_info LMO_solver(sparse_vector grad);
    vertex_info away_vertex(sparse_vector grad);
    vertex_info fw_vertex(sparse_vector grad);
    sparse_vector update_x(int it);
    vector<sparse_vector> remove(vector<sparse_vector> ls, int pos);
    vector<int> find_pos();
    double line_search(sparse_vector dt, sparse_vector x_t, double g_t, double gama_max, double tau, double L);
};

bpcg::bpcg(int dim, vector<vector<double> > K, vector<double> z)
{
    this->dim = dim;
    this->K = K;
    this->z = z;
    vector<double> x_0(this->dim, (1/dim));
    this->init_vertex = vector_2_spare(x_0);
    this->alpha = vector<double>(this->dim, 0);
}

bpcg::~bpcg()
{
}

double bpcg::line_search(sparse_vector dt, sparse_vector x_t, double g_t, double gama_max, double tau, double L)
{
    double mu = 0.9;
    double M = mu * 44;

    double n_dt = dt.snorm();

    double gama = min((g_t / (M * n_dt)), gama_max);

    double Q_t = primal(x_t) - gama * g_t + 0.5 * M * n_dt * gama * gama;

    double f_new = primal(x_t - gama * dt);

    while (f_new > Q_t) // && M < 100000)
    {

        M = tau * M;

        gama = min((g_t / (M * n_dt)), gama_max);

        Q_t = primal(x_t) - gama * g_t + 0.5 * M * n_dt * gama * gama;

        f_new = primal(x_t - gama * dt);
    }
    // printf("exit line search....\n");
    // printf("step size %f \n", gama);
    return gama;
};

vector<double> bpcg::spare_2_vector(sparse_vector x)
{
    vector<double> v(this->dim, 0.0);
    int ix = 0;
    if (x.my_vec.size() == 0)
    {
        return v;
    }
    if (x.my_vec.begin() == x.my_vec.end())
    {
        v[(*x.my_vec.begin()).first] = (*x.my_vec.begin()).second;
        return v;
    }
    for (simple_sparse_vector_iterator it = x.my_vec.begin();
         it != x.my_vec.end(); it++)
    {
        v[(*it).first] = (*it).second;
    }

    return v;
};

sparse_vector bpcg::vector_2_spare(vector<double> x)
{
    sparse_vector sp = sparse_vector();
    for (int i = 0; i < this->dim; i++)
    {
        if (x[i] != 0)
        {
            sp.my_vec.push_back(IndexValuePair(i, x[i]));
        }
    }
    return sp;
}

vector<double> operator*(vector<double> lhs, double rhs)
{
    vector<double> x = vector<double>();
    for (int i = 0; i < lhs.size(); i++)
    {
        x.push_back(lhs[i] * rhs);
    }
    return x;
}

vector<int> bpcg::find_pos()
{
    vector<int> idx;
    int ac_size = 0;
    this->s_list = vector<sparse_vector>();
    for (int i = 0; i < this->alpha.size(); i++)
    {
        if (this->alpha[i] > 0)
        {
            idx.push_back(i);
            ac_size += 1;
            vector<IndexValuePair> myv;
            myv.push_back(IndexValuePair(i, 1.0));
            sparse_vector v = sparse_vector(myv);
            this->s_list.push_back(v);
        }
    }
    printf("active set size %d \n", ac_size);
    return idx;
};

sparse_vector bpcg::update_x(int it)
{
    sparse_vector x;
    vector<double> c(this->dim, 0.0);
    for (int i = 0; i < this->I_active.size(); i++)
    {
        int alpha_idx = this->I_active[i];
        sparse_vector tmp = this->alpha[alpha_idx] * this->s_list[i];
        // print_v(this->s_list[i]);
        vector<double> b = spare_2_vector(tmp);
        for (int i = 0; i < this->dim; i++)
        {
            c[i] = c[i] + b[i];
        }
    }
    x = vector_2_spare(c);
    return x;
};

vertex_info bpcg::LMO_solver(sparse_vector grad)
{
    vertex_info info = vertex_info();
    sparse_vector x;
    int min_idx = 0;
    min_idx = grad.arg_min_index();
    x.my_vec.push_back(IndexValuePair(min_idx, 1));
    info.vertex = x;
    info.index = min_idx;
    return info;
};

vertex_info bpcg::away_vertex(sparse_vector grad)
{
    vertex_info info = vertex_info();
    sparse_vector s_0 = this->s_list[0];
    double dot = s_0.dot(grad);
    int idx = 0;
    for (int i = 0; i < this->s_list.size(); i++)
    {
        double d = this->s_list[i].dot(grad);
        if (d > dot)
        {
            dot = d;
            s_0 = this->s_list[i];
            idx = i;
        }
    };
    info.vertex = s_0;
    info.index = this->I_active[idx];
    return info;
};

vertex_info bpcg::fw_vertex(sparse_vector grad)
{
    vertex_info info = vertex_info();
    sparse_vector s_0 = this->s_list[0];
    double dot = s_0.dot(grad);

    int idx = 0;
    for (int i = 0; i < this->s_list.size(); i++)
    {
        double d = this->s_list[i].dot(grad);
        if (d < dot)
        {
            dot = d;
            s_0 = this->s_list[i];
            idx = i;
        }
    };
    info.vertex = s_0;
    info.index = this->I_active[idx];
    return info;
};

double bpcg::primal(sparse_vector x)
{
    // double ret = x.T @ K @ x;
    vector<double> gg = spare_2_vector(x);
    vector<double> ret = vector<double>(dim, 0.0);
    for (int i = 0; i < this->dim; i++)
    {
        double tmp;
        for (int j = 0; j < this->dim; j++)
        {
            tmp = tmp + K[i][j] * gg[j];
        }
        ret[i] = 2 * tmp - z[i];
    }
    sparse_vector sret = vector_2_spare(ret);
    return sret.dot(x);
};

sparse_vector bpcg::gradient_func(sparse_vector x)
{
    // double ret = 2 * K @ x - z;
    vector<double> gg = spare_2_vector(x);
    vector<double> ret = vector<double>(this->dim, 0);
    for (int i = 0; i < this->dim; i++)
    {
        double tmp = 0;
        for (int j = 0; j < this->dim; j++)
        {
            tmp = tmp + K[i][j] * gg[j];
        }
        ret[i] = 2 * tmp - z[i];
    }
    return vector_2_spare(ret);
};

vector<sparse_vector> bpcg::remove(vector<sparse_vector> ls, int pos)
{
    vector<sparse_vector> ret;
    for (int i = 0; i < ls.size(); i++)
    {
        if (i != pos)
        {
            ret.push_back(ls[i]);
        };
    }
    return ret;
}
bpcg_result bpcg::optimize(int max_iter)
{
    int it = 0;
    double eps = 9.87654e-0222;
    sparse_vector grad = gradient_func(init_vertex);
    vertex_info start = LMO_solver(grad);
    sparse_vector x = start.vertex;
    this->alpha[0] = 1.0;
    this->s_list.push_back(x);
    this->I_active.push_back(0);
    bpcg_result res = {x, this->alpha, this->I_active, this->f_vls};

    while (it < max_iter)
    {
        printf("..........................start iteration %d ............................. \n", it);

        sparse_vector grad = gradient_func(x);

        vertex_info w = LMO_solver(grad);

        vertex_info a = away_vertex(grad);

        vertex_info s = fw_vertex(grad);

        /* code */
        sparse_vector d_fw = x - w.vertex;

        sparse_vector d_aw = a.vertex - s.vertex;
        // return bpcg_result();

        double dual_gap = grad.dot(d_fw);

        if (dual_gap < eps)
        {
            printf("Dual gap reach small value %f \n", dual_gap);
            return res;
        }

        double away_gap = grad.dot(d_aw);

        // printf("dual_gap %f \n", dual_gap);
        // printf("away_gap %f \n", away_gap);

        if (away_gap >= dual_gap)
        {
            double d = this->alpha[a.index];
            double step = line_search(d_aw, x, away_gap, d, 1.5, 1.0);
            // printf("step %f \n", step);
            printf("Descent step ... \n");

            if (step < d)
            {
                this->alpha[a.index] = this->alpha[a.index] - step;
                this->alpha[s.index] = this->alpha[s.index] + step;
            }
            else
            {
                printf("Drop step ... \n");
                int remove_id = -1;
                std::vector<int>::iterator position = std::find(this->I_active.begin(), this->I_active.end(), a.index);
                if (position != this->I_active.end())
                { // == myVector.end() means the element was not found
                    this->I_active.erase(position);
                    remove_id += 1;
                }
                this->s_list = remove(this->s_list, remove_id);
                this->alpha[s.index] = this->alpha[s.index] + d;
                this->alpha[a.index] = 0;
            };
        }
        else
        {
            printf("LMO step ... \n");

            double step = line_search(d_fw, x, dual_gap, 1, 1.5, 1);

            this->alpha = this->alpha * (1 - step);

            this->alpha[w.index] = this->alpha[w.index] + step;

            // printf("LMO step %f... \n", this->alpha[w.index]);

            if (step > 1 - eps)
            {
                printf("only fw step...\n");
                this->I_active = vector<int>();
                this->I_active.push_back(w.index);
                this->alpha = this->alpha * 0;
                this->alpha[w.index] = 1;
                this->s_list = vector<sparse_vector>();
                this->s_list.push_back(w.vertex);
            }
        }
        this->I_active = find_pos();
        x = update_x(it);
        double f_v = primal(x);
        printf("primal value %f \n", f_v);
        this->f_vls.push_back(f_v);
        it = it + 1;
    };
    res.alpha = this->alpha;
    res.I_active = this->I_active;
    return res;
};

#endif
