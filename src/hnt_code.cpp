#include <Rcpp.h>
#include<stdio.h>
using namespace Rcpp;

int move(char x,char y);

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
    move(a,c);
  }
  else
  {
    int move(char x,char y);
    han(n-1,a,c,b); //将A上的n-1个盘借助c移动到b
    move(a,c); //将A最后一个移动到C
    han(n-1,b,a,c); //将b上的n-1个借助a移动到c上
  }
}

int move(char x,char y)
{
  printf("%c-->%c\n",x,y);
  return 0;
}
