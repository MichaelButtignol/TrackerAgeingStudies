#include <math.h>

Double_t fitfunction(Double_t *x, Double_t *par){
  Double_t value, value_th, value_0;
  if (x[0]<par[2]) {
    value = par[0] + par[1]*sqrt(x[0]) + par[3]*x[0]*sqrt(x[0]); // x^1/2 + x^3/2
  }
  else
  {
    value_th = par[0] + par[1]*sqrt(par[2]) + par[3]*par[2]*sqrt(par[2]);
    value_0 = value_th - par[4]*par[2];
    value = value_0 + par[4]*x[0]; // pol1
  }
  //else { value=0; TF1::RejectPoint();}
  
  return value;
}

Double_t fitfunction2(Double_t *x, Double_t *par){
  Double_t value, value_th, value_0;
  if (x[0]<par[2]) {
    value = par[0]/(1+exp(-1*par[1]*(x[0]+par[3]))); // sigmoid
  }
  else
  {
    value_th = par[0]/(1+exp(-1*par[1]*(par[2]+par[3])));
    value_0 = value_th - par[4]*par[2];
    value = value_0 + par[4]*x[0]; // pol1
  }
  //else { value=0; TF1::RejectPoint();}
  
  return value;
}

Double_t fitfunctionderiv(Double_t *x, Double_t *par){
  Double_t value, value_th, value_0;
  if (x[0]<par[2]) {
    value = par[0] + par[1]*exp(-1*par[3]*x[0]);
  }
  else
  {
    value_th = par[0] + par[1]*exp(-1*par[3]*par[2]);
    value_0 = value_th - par[4]*par[2];
    value = value_0 + par[4]*x[0]; // pol1
  }
  //else { value=0; TF1::RejectPoint();}
  
  return value;
}

Double_t fitfunctionderiv2(Double_t *x, Double_t *par){
  Double_t value, value_th, value_0;
  if (x[0]<par[2]) {
    value = par[0] + par[1]*x[0];
  }
  else
  {
    value_th = par[0] + par[1]*par[2];
    value_0 = value_th - par[3]*par[2];
    value = value_0 + par[3]*x[0]; // pol1
  }
  
  return value;
}
