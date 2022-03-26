#include <Rcpp.h>
#include<stdio.h>
using namespace Rcpp;

// [[Rcpp::export]]
bool is_odd_cpp(int n=10){
  bool v = (n %2 == 1); //body
  return v; //return value
}

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x){
  return x*2;
}


// [[Rcpp::export]]
void han(int n,char a,char b, char c)
{
  if(n==1)
  {
    int move(char x,char y);
    move(a,b);
  }else
  {
    int move(char x,char y);
    han(n-1,a,c,b);
    move(a,c);
    han(n-1,b,a,c);
  }
}

int move(char x,char y)
{
  printf("%c-->%c\n",x,y);
  return 0;
}
