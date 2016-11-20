#include <iostream>
// test of 2 point BVP solver

#include "gal2ptBVP.h"
#include "utility.h"
#include <cmath>
#include <fstream>

double zero( double x);
double one( double x );
double two( double x );
double five( double x );
double step( double x );
double xminus1( double x );
std::vector<double> bisect( const std::vector<double> x );

int main()
{

  // Example 1
  std::vector<double> x1(5);
  x1[0] = 0.0; x1[1] = 0.2; x1[2] = 0.7;
  x1[3] = 0.8; x1[4] = 1.0;
  Matrix Y01(galerkinBVP(one,zero,zero,zero,x1,0.0,Dirichlet,1.0,Dirichlet));
  std::ofstream example01("example01");
  writexy(x1,Y01,example01);
  example01.close();

  // Example 2
  Matrix Y02(galerkinBVP(one,zero,zero,two,x1,0.0,Dirichlet,1.0,Dirichlet));
  std::ofstream example02("example02");
  writexy(x1,Y02,example02);
  example02.close();

  // Example 3
  Matrix Y03(galerkinBVP(step,zero,zero,zero,x1,0.0,Dirichlet,1.0,Dirichlet));
  std::ofstream example03("example03");
  writexy(x1,Y03,example03);
  example03.close();

  // Example 4
  Matrix Y04(galerkinBVP(one,five,zero,zero,x1,0.0,Dirichlet,1.0,Dirichlet));
  std::ofstream example04("example04");
  writexy(x1,Y04,example04);
  example04.close();

  // Example 5
  std::vector<double> t(bisect(x1));
  std::vector<double> x5(bisect(t));
  Matrix Y05(galerkinBVP(one,five,zero,zero,x5,0.0,Dirichlet,1.0,Dirichlet));
  std::ofstream example05("example05");
  writexy(x5,Y05,example05);
  example05.close();

  // Example 6
  Matrix Y06(galerkinBVP(one,five,five,zero,x5,0.0,Dirichlet,1.0,Dirichlet));
  std::ofstream example06("example06");
  writexy(x5,Y06,example06);
  example06.close();

  // Example 7
  Matrix Y07(galerkinBVP(one,zero,zero,two,x1,0.0,Dirichlet,0.0,Neumann));
  std::ofstream example07("example07");
  writexy(x1,Y07,example07);
  example07.close();

  // Example 8
  Matrix Y08(galerkinBVP(one,zero,zero,two,x1,0.0,Neumann,0.0,Dirichlet));
  std::ofstream example08("example08");
  writexy(x1,Y08,example08);
  example08.close();

  // Example 9
  Matrix Y09(galerkinBVP(one,zero,zero,zero,x1,0.0,Dirichlet,1.0,Neumann));
  std::ofstream example09("example09");
  writexy(x1,Y09,example09);
  example09.close();

  // Example 10
  std::vector<double> x10(21);
  unsigned int i;
  for( i=0; i<21; i++) x10[i] = 0.1*i;
  Matrix Y10(galerkinBVP(one,zero,one,xminus1,x10,0.0,Neumann,0.0,Neumann));
  std::ofstream example10("example10");
  writexy(x10,Y10,example10);
  example10.close();

  // Example 11
  std::vector<double> x11(81);
  for( i=0; i<81; i++) x11[i] = 0.025*i;
  Matrix Y11(galerkinBVP(one,zero,one,xminus1,x11,0.0,Neumann,0.0,Neumann));
  std::ofstream example11("example11");
  writexy(x11,Y11,example11);
  example10.close();

  return 0;
}

double zero( double x)
{
  return 0.0;
}

double one( double x )
{
  return 1.0;
}

double two( double x )
{
  return 2.0;
}

double five( double x )
{
  return 5.0;
}

double step( double x )
{
  if( x<=0.7 ) return 4.0;
  else return 1.0;
}

double xminus1( double x )
{
  return x-1.0;
}

std::vector<double> bisect( const std::vector<double> x )
{
  unsigned long len = 0;
  if( x.size() > 0 )
    len = 2*x.size() -1;
  std::vector<double> t(len);
  unsigned int i;
  for(i=0; i<x.size(); i++ ){
    t[2*i] = x[i];
    if( i> 0 )
      t[2*i-1] = 0.5*(x[i-1]+x[i]);
  }
  return t;
}



