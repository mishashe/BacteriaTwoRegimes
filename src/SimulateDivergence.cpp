#include <Rcpp.h>
#include <iostream>
#include <random>
#include <chrono>
#include <math.h>
#include <stdint.h>

// [[Rcpp::depends(RcppArmadillo)]]
// #define ARMA_WARN_LEVEL 0
// #include <RcppArmadillo.h>
using namespace Rcpp;
// using namespace arma;
using namespace std;


string Mutate(string letter, int ran)
{
  string newletter;
  if (letter=="a")
  {
    if (ran==1) newletter = "c";
    else if (ran==2) newletter = "g";
    else newletter = "t";
  }
  else if (letter=="c")
  {
    if (ran==1) newletter = "a";
    else if (ran==2) newletter = "t";
    else  newletter = "g";
  }
  else if (letter=="t")
  {
    if (ran==1) newletter = "c";
    else if (ran==2) newletter = "g";
    else newletter = "a";
  }
  else if (letter=="g")
  {
    if (ran==1) newletter = "a";
    else if (ran==2) newletter = "t";
    else newletter = "c";
  }
  return(newletter);
}


// [[Rcpp::export]]
CharacterVector SimulateDivergence(CharacterVector Genome, long int a, double tau)
{
  double muc = 1.0e-2/100.0e6; //https://www.pnas.org/doi/epdf/10.1073/pnas.96.22.12638 1.5e-10
  double mus = 8.9e-11*200.0*941000.0/4.6e6;
  double rho = 1.0e6/tau/tau;
  long int block = 1e4;
  double alpha = 1.0;
  long int L = Genome.size();

  std::mt19937 gen_lognor(a);
  std::lognormal_distribution<double> dis_lognor((log(mus)+log(muc))/2.0,(log(mus)-log(muc))/2.0/3.0);

  // std::mt19937 gen_gamma(a);
  // std::gamma_distribution<> dis_gamma(1.5,mus);

  //std::random_device rd_int;
  std::mt19937 gen_int(a);
  std::uniform_int_distribution<> dis_int(1, 3);

  //std::random_device rd;
  std::mt19937 gen(a);
  std::uniform_real_distribution<> dis(0.0, 1.0);

  //std::random_device rd_exp;
  std::mt19937 gen_exp(a);
  std::exponential_distribution<> dis_exp(rho);

  //std::random_device rd_poi;
  std::mt19937 gen_poi(a);

  NumericVector mu(L);
  double tauHGT = dis_exp(gen_exp);

   mu(0) = pow(dis(gen),1.0/(1.0+alpha))*(mus-muc)+muc; //mu(0) = mus;
   for (long int i=1;i<L;i++)
   {
     if (dis(gen)*block < 1.0)
     {
       mu(i) = pow(dis(gen),1.0/(1.0+alpha))*(mus-muc)+muc;  //mu(i) = mus;
     }
     else
     {
       mu(i) = mu(i-1);
     }
   }
   //
   // for (long int i=0;i<(long int)((1.0*L)/100.0);i++)
   // {
   //   long int pos = (int)(dis(gen)*(1.0*L-12.0));
   //   for (long int j=pos;j<pos+20+1;j++)
   //   {
   //     mu(j) = muc;
   //   }
   // }


  for (long int i=0;i<L;i++)
  {

    if (dis(gen)*block < 1.0)
    {
      tauHGT = dis_exp(gen_exp);
    }

    std::poisson_distribution<> dis_poi(mu(i)*min(tau,tauHGT));
    long int nmuts = min(dis_poi(gen_poi),10);
    if (nmuts>0)
    {
      for (long int n=1;n <= nmuts; n++) Genome(i) = Mutate(as<std::string>(Genome(i)), dis_int(gen_int));
    }
  }
  return(Genome);
}
