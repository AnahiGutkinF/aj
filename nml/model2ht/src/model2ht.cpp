#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double M2HT(NumericVector Q0, NumericVector chd, NumericVector NN){
  
  NumericVector Q = pnorm(Q0);
  
  double d_o = Q[0];
  double dn = Q[1];
  double g3 = Q[2];
  double s_g2 = Q[3];
  double s_g1 = Q[4];
  double Ldo = Q[5];
  double Ldn = Q[6];
  double s_Lgo = Q[7];
  double s_Lgn = Q[8];
  double Hdof = Q[9];
  double Hdos = Q[10];
  double Hdnf = Q[11];
  double Hdns = Q[12];
  double s_Hgof = Q[13];
  double s_Hgos = Q[14];
  double s_Hgnf = Q[15];
  double s_Hgns = Q[16];
  
  Rcpp::NumericVector ee(48);  
  
  // Condition 1
  ee[0] = d_o*Ldo*Hdof + (1-d_o)*(g3*s_g2*s_g1)*(Ldo*s_Lgo)*(Hdof*s_Hgof);                    // Hit | high0
  ee[1] = d_o*(1-Ldo)*Hdos + (1-d_o)*(g3*s_g2*s_g1)*(1-(Ldo*s_Lgo))*(Hdos*s_Hgos);            // Hit | low0
  ee[2] = d_o*Ldo*(1-Hdof) + (1-d_o)*(g3*s_g2*s_g1)*(Ldo*s_Lgo)*(1-(Hdof*s_Hgof));            // Hit | high1
  ee[3] = d_o*(1-Ldo)*(1-Hdos) + (1-d_o)*(g3*s_g2*s_g1)*(1-(Ldo*s_Lgo))*(1-(Hdos*s_Hgos));    // Hit | low1
  ee[4] = (1-d_o)*(1-(g3*s_g2*s_g1))*(Ldn*s_Lgn)*(Hdnf*s_Hgnf);                             // Miss | high0
  ee[5] = (1-d_o)*(1-(g3*s_g2*s_g1))*(1-(Ldn*s_Lgn))*(Hdns*s_Hgns);                         // Miss | low0
  ee[6] = (1-d_o)*(1-(g3*s_g2*s_g1))*(Ldn*s_Lgn)*(1-(Hdnf*s_Hgnf));                         // Miss | high1
  ee[7] = (1-d_o)*(1-(g3*s_g2*s_g1))*(1-(Ldn*s_Lgn))*(1-(Hdns*s_Hgns));                     // Miss | low1
  
  ee[8] = (1-dn)*(g3*s_g2*s_g1)*(Ldo*s_Lgo)*(Hdof*s_Hgof);                                 // FA  | high0
  ee[9] = (1-dn)*(g3*s_g2*s_g1)*(1-(Ldo*s_Lgo))*(Hdos*s_Hgos);                             // FA  | low0
  ee[10] = (1-dn)*(g3*s_g2*s_g1)*(Ldo*s_Lgo)*(1-(Hdof*s_Hgof));                            // FA  | high1
  ee[11] = (1-dn)*(g3*s_g2*s_g1)*(1-(Ldo*s_Lgo))*(1-(Hdos*s_Hgos));                        // FA  | low1
  ee[12] = dn*Ldn*Hdnf + (1-dn)*(1-(g3*s_g2*s_g1))*(Ldn*s_Lgn)*(Hdnf*s_Hgnf);              // CR  | high0
  ee[13] = dn*(1-Ldn)*Hdns + (1-dn)*(1-(g3*s_g2*s_g1))*(1-(Ldn*s_Lgn))*(Hdns*s_Hgns);      // CR  | low0
  ee[14] = dn*Ldn*(1-Hdnf) + (1-dn)*(1-(g3*s_g2*s_g1))*(Ldn*s_Lgn)*(1-(Hdnf*s_Hgnf));      // CR  | high1
  ee[15] = dn*(1-Ldn)*(1-Hdns) + (1-dn)*(1-(g3*s_g2*s_g1))*(1-(Ldn*s_Lgn))*(1-(Hdns*s_Hgns));  // CR  | low1
  
  // Condition 2
  ee[16] = d_o*Ldo*Hdof + (1-d_o)*g3*s_g2*(Ldo*s_Lgo)*(Hdof*s_Hgof);                         // Hit | high0
  ee[17] = d_o*(1-Ldo)*Hdos + (1-d_o)*g3*s_g2*(1-(Ldo*s_Lgo))*(Hdos*s_Hgos);                 // Hit | low0
  ee[18] = d_o*Ldo*(1-Hdof) + (1-d_o)*g3*s_g2*(Ldo*s_Lgo)*(1-(Hdof*s_Hgof));                 // Hit | high1
  ee[19] = d_o*(1-Ldo)*(1-Hdos) + (1-d_o)*g3*s_g2*(1-(Ldo*s_Lgo))*(1-(Hdos*s_Hgos));         // Hit | low1
  ee[20] = (1-d_o)*(1-(g3*s_g2))*(Ldn*s_Lgn)*(Hdnf*s_Hgnf);                                // Miss | high0
  ee[21] = (1-d_o)*(1-(g3*s_g2))*(1-(Ldn*s_Lgn))*(Hdns*s_Hgns);                            // Miss | low0
  ee[22] = (1-d_o)*(1-(g3*s_g2))*(Ldn*s_Lgn)*(1-(Hdnf*s_Hgnf));                            // Miss | high1
  ee[23] = (1-d_o)*(1-(g3*s_g2))*(1-(Ldn*s_Lgn))*(1-(Hdns*s_Hgns));                        // Miss | low1
  
  ee[24] = (1-dn)*(g3*s_g2)*(Ldo*s_Lgo)*(Hdof*s_Hgof);                                    // FA  | high0
  ee[25] = (1-dn)*(g3*s_g2)*(1-(Ldo*s_Lgo))*(Hdos*s_Hgos);                                // FA  | low0
  ee[26] = (1-dn)*(g3*s_g2)*(Ldo*s_Lgo)*(1-(Hdof*s_Hgof));                                // FA  | high1
  ee[27] = (1-dn)*(g3*s_g2)*(1-(Ldo*s_Lgo))*(1-(Hdos*s_Hgos));                            // FA  | low1
  ee[28] = dn*Ldn*Hdnf + (1-dn)*(1-(g3*s_g2))*(Ldn*s_Lgn)*(Hdnf*s_Hgnf);                  // CR  | high0
  ee[29] = dn*(1-Ldn)*Hdns + (1-dn)*(1-(g3*s_g2))*(1-(Ldn*s_Lgn))*(Hdns*s_Hgns);          // CR  | low0
  ee[30] = dn*Ldn*(1-Hdnf) + (1-dn)*(1-(g3*s_g2))*(Ldn*s_Lgn)*(1-(Hdnf*s_Hgnf));          // CR  | high1
  ee[31] = dn*(1-Ldn)*(1-Hdns) + (1-dn)*(1-(g3*s_g2))*(1-(Ldn*s_Lgn))*(1-(Hdns*s_Hgns));  // CR  | low1
  
  // Condition 3
  ee[32] = d_o*Ldo*Hdof + (1-d_o)*g3*(Ldo*s_Lgo)*(Hdof*s_Hgof);                             // Hit | high0
  ee[33] = d_o*(1-Ldo)*Hdos + (1-d_o)*g3*(1-(Ldo*s_Lgo))*(Hdos*s_Hgos);                     // Hit | low0
  ee[34] = d_o*Ldo*(1-Hdof) + (1-d_o)*g3*(Ldo*s_Lgo)*(1-(Hdof*s_Hgof));                     // Hit | high1
  ee[35] = d_o*(1-Ldo)*(1-Hdos) + (1-d_o)*g3*(1-(Ldo*s_Lgo))*(1-(Hdos*s_Hgos));             // Hit | low1
  ee[36] = (1-d_o)*(1-g3)*(Ldn*s_Lgn)*(Hdnf*s_Hgnf);                                       // Miss | high0
  ee[37] = (1-d_o)*(1-g3)*(1-(Ldn*s_Lgn))*(Hdns*s_Hgns);                                   // Miss | low0
  ee[38] = (1-d_o)*(1-g3)*(Ldn*s_Lgn)*(1-(Hdnf*s_Hgnf));                                   // Miss | high1
  ee[39] = (1-d_o)*(1-g3)*(1-(Ldn*s_Lgn))*(1-(Hdns*s_Hgns));                               // Miss | low1
  
  ee[40] = (1-dn)*g3*(Ldo*s_Lgo)*(Hdof*s_Hgof);                                           // FA  | high0
  ee[41] = (1-dn)*g3*(1-(Ldo*s_Lgo))*(Hdos*s_Hgos);                                       // FA  | low0
  ee[42] = (1-dn)*g3*(Ldo*s_Lgo)*(1-(Hdof*s_Hgof));                                       // FA  | high1
  ee[43] = (1-dn)*g3*(1-(Ldo*s_Lgo))*(1-(Hdos*s_Hgos));                                   // FA  | low1
  ee[44] = dn*Ldn*Hdnf + (1-dn)*(1-g3)*(Ldn*s_Lgn)*(Hdnf*s_Hgnf);                         // CR  | high0
  ee[45] = dn*(1-Ldn)*Hdns + (1-dn)*(1-g3)*(1-(Ldn*s_Lgn))*(Hdns*s_Hgns);                 // CR  | low0
  ee[46] = dn*Ldn*(1-Hdnf) + (1-dn)*(1-g3)*(Ldn*s_Lgn)*(1-(Hdnf*s_Hgnf));                 // CR  | high1
  ee[47] = dn*(1-Ldn)*(1-Hdns) + (1-dn)*(1-g3)*(1-(Ldn*s_Lgn))*(1-(Hdns*s_Hgns));         // CR  | low1
  
  // Multiplicamos por NN
  ee = ee * NN;
  
  // Cálculo de la log-verosimilitud acumulada
  double LL = 0.0;
  for (int ii = 0; ii < 48; ii++) {
    if (chd[ii] > 0.0) {
      LL += chd[ii] * (log(chd[ii]) - log(ee[ii]));
    }
  }
  
  return 2.0 * LL;
}
