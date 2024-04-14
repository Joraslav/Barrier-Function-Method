#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <functional>

using namespace std;
using type = long double;
template <class T>
using vector_tmpl = vector<T>;

using vector_type = vector_tmpl<type>;

template <class T>
vector_tmpl<T> operator+(vector_tmpl<T> const& l, vector_tmpl<T> const& r)
{
  if (l.size() != r.size()){throw std::invalid_argument("wrong sizes");}  //Проверка размерности

  auto Rez{l};
  auto const nRows{Rez.size()};

  for (auto i{0u}; i < nRows; ++i)
  {
    Rez[i] += r[i];
  }
  
  return Rez;
}

template <class T>
vector_tmpl<T> operator-(vector_tmpl<T> const& l, vector_tmpl<T> const& r)
{
  if (l.size() != r.size()){throw std::invalid_argument("wrong sizes");}  //Проверка размерности

  auto Rez{l};
  auto const nRows{Rez.size()};

  for (auto i{0u}; i < nRows; ++i)
  {
    Rez[i] -= r[i];
  }
  
  return Rez;
}

template<class T>
type operator*(vector_tmpl<T> const& l, vector_tmpl<T> const& r)
{
  if (l.size() != r.size()){throw std::invalid_argument("wrong sizes");}  //Проверка размерности

  type Rez = 0.0;
  auto const nRows{l.size()};

  for (auto i{0u}; i < nRows; ++i)
  {
    Rez += l[i] * r[i];
  }
  
  return Rez;
}

template<class T>
vector_tmpl<T> operator*(type const& c, vector_tmpl<T> const& r)
{
  vector_tmpl<T> Rez{r};
  auto const nRows{Rez.size()};
  for (auto i{0u}; i < nRows; ++i)
  {
    Rez[i] = c*Rez[i];
  }
  return Rez;
}

template<class T>
ostream& operator<<(ostream& out, vector_tmpl<T> const& v)
{
  for (auto const &i : v)
  {
    out << i << '\t';
  }
  return out;
}

type Norm(vector_type const& v)
{
  auto Rez = 0.0;
  auto ScalarProd = v*v;

  Rez = sqrt(ScalarProd);

  return Rez;
}

type f(vector_type const& x)
{
  return pow(x[0]-2,4)+pow(x[0]-2*x[1],2);
}

type g1(vector_type const& x)
{
  return 1.-2*x[0]-x[1];
}

type g2(vector_type const& x)
{
  return 1.-2*x[0]-pow(x[1],3);
}

type P(vector_type const& x, type const& r)
{
  return -r/(g1(x)+g2(x));
}

type BF(vector_type const& x, type const& r)
{
  return f(x) + P(x,r);
}

template <class T>
vector_tmpl<T> Grad_P(vector_tmpl<T> const& x, type const& r)
{
  vector_tmpl<T> Rez;
  Rez.push_back(-r*(2./pow(g1(x),2)+2./pow(g2(x),2)));
  Rez.push_back(-r*(1./pow(g1(x),2)+(3*pow(x[1],2))/pow(g2(x),2)));
  return Rez;
}

template <class T>
vector_tmpl<T> Grad(vector_tmpl<T> const& x, type const& r)
{
  vector_tmpl<T> Rez;
  Rez.push_back(4.*pow(x[0]-2,3)+2.*(x[0]-2.*x[1])+Grad_P(x,r)[0]);
  Rez.push_back(-4.*(x[0]-2.*x[1])+Grad_P(x,r)[1]);
  return Rez;
}

int main()
{
  type r = 100;
  type c = 10;

  return 0;
}
