#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
double SCR2(NumericVector Q0, NumericVector chd, NumericVector NN){
  
  NumericVector Q = pnorm(Q0);
  
  double Vt = Q[0];
  double Gt = Q[1];
  double Vr = Q[2];
  double Gr = Q[3];
  double a  = Q[4];
  double b  = Q[5];
  double Vtx = Q[6];
  double Gtx = Q[7];
  double Vrx = Q[8];
  double Grx = Q[9];
  
  Rcpp::NumericVector ee(15);
  
  ee[0] = Vt + (1-Vt)*Gt*a + (1-Vt)*(1-Gt)*b*a;
  ee[1] = (1-Vt)*Gt*(1-a) + (1-Vt)*(1-Gt)*b*(1-a);
  ee[2] = (1-Vt)*(1-Gt)*(1-b);
  
  ee[3] = (1-Vr)*Gr*a + (1-Vr)*(1-Gr)*b*a;
  ee[4] = Vr + (1-Vr)*Gr*(1-a) + (1-Vr)*(1-Gr)*b*(1-a);
  ee[5] = (1-Vr)*(1-Gr)*(1-b);
  
  ee[6] = Vtx + (1-Vtx)*Gtx*a + (1-Vtx)*(1-Gtx)*b*a;
  ee[7] = (1-Vtx)*Gtx*(1-a) + (1-Vtx)*(1-Gtx)*b*(1-a);
  ee[8] = (1-Vtx)*(1-Gtx)*(1-b);
  
  ee[9] = (1-Vrx)*Grx*a + (1-Vrx)*(1-Grx)*b*a;
  ee[10] = Vrx + (1-Vrx)*Grx*(1-a) + (1-Vrx)*(1-Grx)*b*(1-a);
  ee[11] = (1-Vrx)*(1-Grx)*(1-b);
  
  ee[12] = b*a;
  ee[13] = b*(1-a);
  ee[14] = (1-b);
  
  ee = ee*NN;
  
  double LL = 0.0;
  
  for(int ii=0; ii<15; ii++){
    if(chd[ii] > 0.0){
      LL = LL + chd[ii]*(log(chd[ii])-log(ee[ii]));
    }
  }
  return 2.0*LL;
}