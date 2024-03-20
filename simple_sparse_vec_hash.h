// Distributed under GNU General Public License (see license.txt for details).
//
//  Copyright (c) 2007 Shai Shalev-Shwartz.
//  All Rights Reserved.
//==============================================================================
// File Name: simple_sparse_vec_hash.h
// Written by: Shai Shalev-Shwartz (28.01.07)
// implement a very simple hash table and sparse vector based on stl vector
// and the modulu function
//==============================================================================

#ifndef _SHAI_SIMPLE_HASH_TABLE
#define _SHAI_SIMPLE_HASH_TABLE

//*****************************************************************************
// Included Files
//*****************************************************************************
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <sstream>
// #include <WeightVector.h>

using namespace std;

class IndexValuePair
{
public:
  IndexValuePair() {}
  IndexValuePair(uint i, double v) : first(i), second(v) {}
  uint first;
  double second;
};

typedef std::vector<IndexValuePair>::iterator simple_sparse_vector_iterator;

/** Implements a simple sparse vector based on stl vector.
    @author Shai Shalev-Shwartz (shais@cs.huji.ac.il)
*/
class simple_sparse_vector
{
public:
  /** Default Constructor. Allocates an all zeros sparse vector
   */
  simple_sparse_vector() {}

  /** Constructor. Read a line from a file of pairs of string and value.
      Construct a vector from this line.
      @param is The input stream to read from.
  */
  simple_sparse_vector(std::istream &is);

  /** Constructor. Construct a vector from an istringstream.
      @param is The input stringstream to read from.
      @param n The number of elements to read
  */
  simple_sparse_vector(std::istringstream &is, int n);

  simple_sparse_vector(std::vector<IndexValuePair> my_vec);

  // Multiply a sparse vector by a scalar s
  void scale(double s);

  simple_sparse_vector multiply(double s);

  /** Returns the squared l_2 norm of the vector
      @return the squared l_2 norm of the vector
  */
  double snorm();

  void put(vector<double> *x);

  double dot(std::vector<double> x);

  double dot(simple_sparse_vector x);
  void add(simple_sparse_vector x);

  /** Convert the vector to a binary vector that just indicates
      which elements are non-zero
  */
  void make_binary();

  double get_second(int i);

  // return the maximal non-zero index
  uint max_index();

  uint min_index();

  uint arg_max_index();

  uint arg_min_index();

  // Print the content of a sparse instance to an output stream
  void print(std::ostream &os);

  // Zero all indexed elements of a sparse vector
  void zero();

  std::vector<IndexValuePair> my_vec;
};

void simple_sparse_vector::put(vector<double> *x)
{
  for (simple_sparse_vector_iterator it = my_vec.begin();
       it != my_vec.end(); it++)
  {
    int idx = (*it).first;
    double value = (*it).second;
    (*x)[idx] = value;
  }
}

simple_sparse_vector operator*(double s, simple_sparse_vector t)
{
  simple_sparse_vector ret = simple_sparse_vector();
  for (simple_sparse_vector_iterator it = t.my_vec.begin();
       it != t.my_vec.end(); it++)
  {
    int idx = (*it).first;
    double value = (*it).second * s;
    ret.my_vec.push_back(IndexValuePair(idx, value));
  }
  return ret;
};

simple_sparse_vector operator+(simple_sparse_vector t, double s)
{
  simple_sparse_vector ret = t;
  for (simple_sparse_vector_iterator it = ret.my_vec.begin();
       it != ret.my_vec.end(); it++)
  {
    (*it).second += s;
  }
  return ret;
};


void print_v(simple_sparse_vector v)
{

  if (v.my_vec.size() == 0)
  {
    return;
  }
  for (simple_sparse_vector_iterator it = v.my_vec.begin();
       it != v.my_vec.end(); it++)
  {

    cout << " " << (*it).first << ":" << (*it).second << " ";
  }

  cout << endl;
};
/*---------------------------------------------------------------------------*/
simple_sparse_vector operator+(simple_sparse_vector lhs, simple_sparse_vector rhs)
{
  if (lhs.my_vec.empty())
  {
    std::vector<IndexValuePair> rhs_vec = rhs.my_vec;
    return simple_sparse_vector(rhs_vec);
  }
  if (rhs.my_vec.empty())
  {
    return lhs;
  }
  std::vector<IndexValuePair> lhs_vec = lhs.my_vec;
  std::vector<IndexValuePair> rhs_vec = rhs.my_vec;
  uint mx1 = lhs.max_index();

  uint mx2 = rhs.max_index();
  // return simple_sparse_vector();
  uint mx;
  if (mx1 > mx2)
  {
    mx = mx1;
  }
  else
  {
    mx = mx2;
  }
  mx = mx + 1;

  std::vector<double> w(mx, 0.0);
  std::vector<IndexValuePair> my_vec;

  if (lhs_vec.begin() == lhs_vec.end())
  {
    w[(*lhs_vec.begin()).first] = (*lhs_vec.begin()).second;
  }
  else
  {
    for (simple_sparse_vector_iterator it = lhs_vec.begin();
         it != lhs_vec.end(); it++)
    {
      uint id = (*it).first;
      w[id] = w[id] + (*it).second;
    }
  }

  if (rhs_vec.begin() == rhs_vec.end())
  {
    w[(*rhs_vec.begin()).first] += (*rhs_vec.begin()).second;
  }
  else
  {
    for (simple_sparse_vector_iterator it = rhs_vec.begin();
         it != rhs_vec.end(); it++)
    {
      uint id = (*it).first;
      w[id] = w[id] + (*it).second;
    }
  }

  for (int i = 0; i < mx; ++i)
  {

    // read the pair (key,val)
    if (w[i] != 0.0)
    {
      my_vec.push_back(IndexValuePair(i, w[i]));
    }
  };

  return simple_sparse_vector(my_vec);
};

simple_sparse_vector operator-(simple_sparse_vector lhs, simple_sparse_vector rhs)
{
  if (lhs.my_vec.empty())
  {
    std::vector<IndexValuePair> rhs_vec = rhs.my_vec;
    simple_sparse_vector ret = simple_sparse_vector(rhs_vec);
    ret.scale(-1.0);
    return ret;
  }
  if (rhs.my_vec.empty())
  {
    std::vector<IndexValuePair> lhs_vec = lhs.my_vec;
    simple_sparse_vector ret = simple_sparse_vector(lhs_vec);
    return ret;
  }
  std::vector<IndexValuePair> lhs_vec = lhs.my_vec;
  std::vector<IndexValuePair> rhs_vec = rhs.my_vec;
  uint mx1 = lhs.max_index();

  uint mx2 = rhs.max_index();

  // return simple_sparse_vector();
  uint mx;
  if (mx1 > mx2)
  {
    mx = mx1;
  }
  else
  {
    mx = mx2;
  }
  mx = mx + 1;

  std::vector<double> w(mx, 0.0);
  std::vector<IndexValuePair> my_vec;

  if (lhs_vec.begin() == lhs_vec.end())
  {
    w[(*lhs_vec.begin()).first] = (*lhs_vec.begin()).second;
  }
  else
  {
    for (simple_sparse_vector_iterator it = lhs_vec.begin();
         it != lhs_vec.end(); it++)
    {
      uint id = (*it).first;
      w[id] = w[id] + (*it).second;
    }
  }

  if (rhs_vec.begin() == rhs_vec.end())
  {
    w[(*rhs_vec.begin()).first] -= (*rhs_vec.begin()).second;
  }
  else
  {
    for (simple_sparse_vector_iterator it = rhs_vec.begin();
         it != rhs_vec.end(); it++)
    {
      uint id = (*it).first;
      w[id] = w[id] - (*it).second;
    }
  }

  for (int i = 0; i < mx; ++i)
  {

    // read the pair (key,val)
    if (w[i] != 0.0)
    {
      my_vec.push_back(IndexValuePair(i, w[i]));
    }
  };

  return simple_sparse_vector(my_vec);
};

double simple_sparse_vector::dot(simple_sparse_vector x)
{
  int i = 0;
  int j = 0;
  double res = 0.0;
  while (i < my_vec.size() && j < x.my_vec.size())
  {
    /* code */
    if (my_vec[i].first < x.my_vec[j].first)
    {
      i++;
    }
    else if (my_vec[i].first > x.my_vec[j].first)
    {
      j++;
    }
    else
    {
      res += (my_vec[i].second * x.my_vec[j].second);
      i++;
      j++;
    }
  }
  return res;
};

void simple_sparse_vector::add(simple_sparse_vector rhs)
{
  if (this->my_vec.size() == 0)
  {
    if (rhs.my_vec.size() == 0)
    {
      return;
    }
    for (int i = 0; i < rhs.my_vec.size(); i++)
    {
      this->my_vec.push_back(rhs.my_vec[i]);
    }
    return;
  }
  if (rhs.my_vec.size() == 0)
  {
    return;
  }
  std::vector<IndexValuePair> rhs_vec = rhs.my_vec;
  uint mx1 = max_index();
  uint mx2 = rhs.max_index();
  // return simple_sparse_vector();
  uint mx;
  if (mx1 > mx2)
  {
    mx = mx1;
  }
  else
  {
    mx = mx2;
  }
  mx = mx + 1;
  vector<IndexValuePair> new_vec;
  std::vector<double> w(mx, 0.0);
  if (my_vec.begin() == my_vec.end())
  {
    w[(*my_vec.begin()).first] = (*my_vec.begin()).second;
  }
  else
  {
    for (simple_sparse_vector_iterator it = my_vec.begin();
         it != my_vec.end(); it++)
    {
      uint id = (*it).first;
      w[id] = w[id] + (*it).second;
    }
  }

  if (rhs_vec.begin() == rhs_vec.end())
  {
    w[(*rhs_vec.begin()).first] += (*rhs_vec.begin()).second;
  }
  else
  {
    for (simple_sparse_vector_iterator it = rhs_vec.begin();
         it != rhs_vec.end(); it++)
    {
      uint id = (*it).first;
      w[id] = w[id] + (*it).second;
    }
  }

  for (int i = 0; i < mx; ++i)
  {

    // read the pair (key,val)
    if (w[i] != 0.0)
    {
      new_vec.push_back(IndexValuePair(i, w[i]));
    }
  };
  my_vec = new_vec;
};
/*---------------------------------------------------------------------------*/
simple_sparse_vector::simple_sparse_vector(std::istream &is) : my_vec()
{

  // read the number of elements
  int n = 0;
  is >> n;

  // read the elements
  for (int i = 0; i < n; ++i)
  {

    // read the pair (key,val)
    uint key;
    is >> key;
    double val;
    is >> val;
    // insert to the map
    my_vec.push_back(IndexValuePair(key, val));
  }
}
/*---------------------------------------------------------------------------*/
simple_sparse_vector::simple_sparse_vector(std::vector<IndexValuePair> new_vec) : my_vec()
{
  my_vec = new_vec;
}
/*---------------------------------------------------------------------------*/
simple_sparse_vector::simple_sparse_vector(std::istringstream &is, int n) : my_vec()
{

  // read the elements
  for (int i = 0; i < n; ++i)
  {

    // read the pair (key,val)
    uint key;
    is >> key;
    double val;
    is >> val;

    // insert to the map
    my_vec.push_back(IndexValuePair(key, val));
  }
}

/*---------------------------------------------------------------------------*/
double simple_sparse_vector::dot(std::vector<double> x)
{
  double res = 0.0;
  for (simple_sparse_vector_iterator it = my_vec.begin();
       it != my_vec.end(); it++)
  {
    res += (*it).second * x[(*it).first];
  }
  return res;
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
void simple_sparse_vector::scale(double s)
{

  for (simple_sparse_vector_iterator it = my_vec.begin();
       it != my_vec.end(); it++)
  {
    (*it).second *= s;
  }
}

simple_sparse_vector simple_sparse_vector::multiply(double s)
{

  for (simple_sparse_vector_iterator it = my_vec.begin();
       it != my_vec.end(); it++)
  {
    (*it).second *= s;
  }

  return simple_sparse_vector(my_vec);
}

double simple_sparse_vector::get_second(int i)
{

  for (simple_sparse_vector_iterator it = my_vec.begin();
       it != my_vec.end(); it++)
  {
    if ((*it).first == i)
    {
      return (*it).second;
    }
  }
  return 0.0;
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

uint simple_sparse_vector::arg_max_index()
{
  double se;
  uint fst = 0;
  int i = 0;
  for (const IndexValuePair &in : my_vec)
  {
    if (i == 0)
    {
      se = in.second;
      i = 1;
    };
    // std::cout << in.second;
    if (in.second > se)
    {
      fst = in.first;
      se = in.second;
    };
    i = i + 1;
  };
  return fst;
}
/*---------------------------------------------------------------------------*/

uint simple_sparse_vector::arg_min_index()
{

  double se;
  uint fst = 0;
  int i = 0;
  for (const IndexValuePair &in : my_vec)
  {
    if (i == 0)
    {
      se = in.second;
      i = 1;
    };
    // std::cout << in.second;
    if (in.second < se)
    {
      fst = in.first;
      se = in.second;
    };
    i = i + 1;
  };
  return fst;
}

uint simple_sparse_vector::max_index()
{
  int d = 0;
  if (my_vec.begin() == my_vec.end())
  {
    return (*my_vec.begin()).first;
  }
  for (simple_sparse_vector_iterator it = my_vec.begin();
       it != my_vec.end(); it++)
  {
    int inx = (*it).first;
    if (inx > d)
    {
      d = inx;
    };
  };
  return d;
}
/*---------------------------------------------------------------------------*/

uint simple_sparse_vector::min_index()
{
  if ((int)my_vec.size() == 0)
  {
    return 0;
  }
  return (*my_vec.begin()).first;
}

/*---------------------------------------------------------------------------*/
double simple_sparse_vector::snorm()
{

  double output = 0.0;
  for (simple_sparse_vector_iterator it = my_vec.begin();
       it != my_vec.end(); it++)
  {
    double tmp = (*it).second;
    output += tmp * tmp;
  }

  return output;
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
void simple_sparse_vector::make_binary()
{

  for (simple_sparse_vector_iterator it = my_vec.begin();
       it != my_vec.end(); it++)
  {
    (*it).second = 1.0;
  }
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
void simple_sparse_vector::print(std::ostream &os)
{

  for (simple_sparse_vector_iterator it = my_vec.begin();
       it != my_vec.end(); it++)
  {
    os << "(" << (*it).first << "," << (*it).second << ") ";
  }
  os << std::endl;
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
void simple_sparse_vector::zero()
{
  my_vec.clear();
}
/*---------------------------------------------------------------------------*/

#ifdef nodef

typedef std::vector<std::vector<IndexValuePair>>::iterator simple_hash_table_iterator;

/** Implements a simple hash table based on stl vector.
    The hash function is the modulu %
    @author Shai Shalev-Shwartz (shais@cs.huji.ac.il)
*/
class simple_hash_table
{
public:
  /** Default Constructor. Allocates an all zeros hash table
   */
  simple_hash_table() : my_vec(1087) {}

  /** Constructor. Allocates an all zeros hash table with m entries
   */
  simple_hash_table(uint m) : my_vec(m) {}

  /** retreive the value of the i'th element of the hash
      @param i the index to retrieve
  */
  double get(uint i);

  /** reference operator. Retrieves the i'th element of the hash and
      create it if it doesn't exist.
      @param i the index to retrieve
  */
  double &get_ref(uint i);

  // Multiply a sparse vector by a scalar s
  void scale(double s);

  // this = a * this + b * other;
  void scale_and_add(simple_sparse_vector &other, double a, double b);

  // this = a * this + b * other;
  void scale_and_add(simple_hash_table &other, double a, double b);

  // this += b*other;
  void add(simple_hash_table &other, double b);

  // this += b*other;
  void add(simple_sparse_vector &other, double b);

  /** Returns the squared l_2 norm of the vector
      @return the squared l_2 norm of the vector
  */
  double snorm();

  // Print the content of a sparse instance to an output stream
  void print(std::ostream &os);

  // Zero all indexed elements of a sparse vector
  void zero();

  std::vector<std::vector<IndexValuePair>> my_vec;
};

//-----------------------------------------------------------------------------
/** Operator * for vector-vector multiplication
    @param u A reference to a simple_sparse_vector
    @param v A reference to a simple_hash_table
    @return The product (double)
*/
double operator*(simple_sparse_vector &u, simple_hash_table &v);

//-----------------------------------------------------------------------------
/** Operator * for vector-vector multiplication
    @param u A reference to a simple_hash_table
    @param v A reference to a simple_sparse_vector
    @return The product (double)
*/
double operator*(simple_hash_table &u, simple_sparse_vector &v);

#endif

#endif
