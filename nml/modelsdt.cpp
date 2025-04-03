#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double MSDT(NumericVector Q0, NumericVector chd, NumericVector NN) {
  
  // Aplicamos pnorm a Q0 para obtener el vector Q
  NumericVector Q = Rcpp::pnorm(Q0);
  
  // Definimos las variables basadas en Q
  double cr11 = Q[0];  // Umbral cr11 para la condición 1
  double cr12 = Q[1];  // Umbral cr12 para la condición 2
  double cr13 = Q[2];  // Umbral cr13 para la condición 3
  double cr2 = Q[3];   // Umbral cr2 común para todas las condiciones
  double cr3 = Q[4];   // Umbral cr3 común para todas las condiciones
  double d = Q[5];     // Desplazamiento d
  double ss = Q[6];    // Desviación estándar ss
  double Lh = Q[7];    // Parámetro Lh
  double Lr = Q[8];    // Parámetro Lr
  double s_Lm = Q[9];  // Parámetro s_Lm
  double s_Lf = Q[10]; // Parámetro s_Lf
  
  Rcpp::NumericVector ee(48);  // Para almacenar los valores de Hit, Miss, FA, RC para todas las condiciones
  
  // ---------------------- CONDITION 1 ----------------------
  ee[0] = (1 - R::pnorm((cr3 + cr2 + cr11 - d) / ss, 0.0, 1.0, true, false)) * Lh;
  ee[1] = (R::pnorm((cr3 + cr2 + cr11 - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr2 + cr11 - d) / ss, 0.0, 1.0, true, false)) * Lh;
  ee[2] = (1 - R::pnorm((cr3 + cr2 + cr11 - d) / ss, 0.0, 1.0, true, false)) * (1 - Lh);
  ee[3] = (R::pnorm((cr3 + cr2 + cr11 - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr2 + cr11 - d) / ss, 0.0, 1.0, true, false)) * (1 - Lh);
  ee[4] = R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false) * (Lh * s_Lm);
  ee[5] = (R::pnorm((cr11 + cr2 - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false)) * (Lh * s_Lm);
  ee[6] = R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false) * (1 - Lh * s_Lm);
  ee[7] = (R::pnorm((cr11 + cr2 - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false)) * (1 - Lh * s_Lm);
  
  ee[8] = (1 - R::pnorm(cr3 + cr2 + cr11, 0.0, 1.0, true, false)) * (Lr * s_Lf);
  ee[9] = (R::pnorm(cr3 + cr2 + cr11, 0.0, 1.0, true, false) - R::pnorm(cr2 + cr11, 0.0, 1.0, true, false)) * (Lr * s_Lf);
  ee[10] = (1 - R::pnorm(cr3 + cr2 + cr11, 0.0, 1.0, true, false)) * (1 - Lr * s_Lf);
  ee[11] = (R::pnorm(cr3 + cr2 + cr11, 0.0, 1.0, true, false) - R::pnorm(cr2 + cr11, 0.0, 1.0, true, false)) * (1 - Lr * s_Lf);
  ee[12] = R::pnorm(cr11, 0.0, 1.0, true, false) * Lr;
  ee[13] = (R::pnorm(cr2 + cr11, 0.0, 1.0, true, false) - R::pnorm(cr11, 0.0, 1.0, true, false)) * Lr;
  ee[14] = R::pnorm(cr11, 0.0, 1.0, true, false) * (1 - Lr);
  ee[15] = (R::pnorm(cr2 + cr11, 0.0, 1.0, true, false) - R::pnorm(cr11, 0.0, 1.0, true, false)) * (1 - Lr);
  
  
  // ---------------------- CONDITION 2 ----------------------
  ee[16] = (1 - R::pnorm((cr3 + cr2 + cr12 - d) / ss, 0.0, 1.0, true, false)) * Lh;
  ee[17] = (R::pnorm((cr3 + cr2 + cr12 - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr2 + cr12 - d) / ss, 0.0, 1.0, true, false)) * Lh;
  ee[18] = (1 - R::pnorm((cr3 + cr2 + cr12 - d) / ss, 0.0, 1.0, true, false)) * (1 - Lh);
  ee[19] = (R::pnorm((cr3 + cr2 + cr12 - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr2 + cr12 - d) / ss, 0.0, 1.0, true, false)) * (1 - Lh);
  ee[20] = R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false) * (Lh * s_Lm);
  ee[21] = (R::pnorm((cr12 + cr2 - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false)) * (Lh * s_Lm);
  ee[22] = R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false) * (1 - Lh * s_Lm);
  ee[23] = (R::pnorm((cr12 + cr2 - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false)) * (1 - Lh * s_Lm);
  
  ee[24] = (1 - R::pnorm(cr3 + cr2 + cr12, 0.0, 1.0, true, false)) * (Lr * s_Lf);
  ee[25] = (R::pnorm(cr3 + cr2 + cr12, 0.0, 1.0, true, false) - R::pnorm(cr2 + cr12, 0.0, 1.0, true, false)) * (Lr * s_Lf);
  ee[26] = (1 - R::pnorm(cr3 + cr2 + cr12, 0.0, 1.0, true, false)) * (1 - Lr * s_Lf);
  ee[27] = (R::pnorm(cr3 + cr2 + cr12, 0.0, 1.0, true, false) - R::pnorm(cr2 + cr12, 0.0, 1.0, true, false)) * (1 - Lr * s_Lf);
  ee[28] = R::pnorm(cr12, 0.0, 1.0, true, false) * Lr;
  ee[29] = (R::pnorm(cr2 + cr12, 0.0, 1.0, true, false) - R::pnorm(cr12, 0.0, 1.0, true, false)) * Lr;
  ee[30] = R::pnorm(cr12, 0.0, 1.0, true, false) * (1 - Lr);
  ee[31] = (R::pnorm(cr2 + cr12, 0.0, 1.0, true, false) - R::pnorm(cr12, 0.0, 1.0, true, false)) * (1 - Lr);
  
  
  // ---------------------- CONDITION 3 ----------------------
  ee[32] = (1 - R::pnorm((cr3 + cr2 + cr13 - d) / ss, 0.0, 1.0, true, false)) * Lh;
  ee[33] = (R::pnorm((cr3 + cr2 + cr13 - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr2 + cr13 - d) / ss, 0.0, 1.0, true, false)) * Lh;
  ee[34] = (1 - R::pnorm((cr3 + cr2 + cr13 - d) / ss, 0.0, 1.0, true, false)) * (1 - Lh);
  ee[35] = (R::pnorm((cr3 + cr2 + cr13 - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr2 + cr13 - d) / ss, 0.0, 1.0, true, false)) * (1 - Lh);
  ee[36] = R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false) * (Lh * s_Lm);
  ee[37] = (R::pnorm((cr13 + cr2 - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false)) * (Lh * s_Lm);
  ee[38] = R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false) * (1 - Lh * s_Lm);
  ee[39] = (R::pnorm((cr13 + cr2 - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false)) * (1 - Lh * s_Lm);
  
  ee[40] = (1 - R::pnorm(cr3 + cr2 + cr13, 0.0, 1.0, true, false)) * (Lr * s_Lf);
  ee[41] = (R::pnorm(cr3 + cr2 + cr13, 0.0, 1.0, true, false) - R::pnorm(cr2 + cr13, 0.0, 1.0, true, false)) * (Lr * s_Lf);
  ee[42] = (1 - R::pnorm(cr3 + cr2 + cr13, 0.0, 1.0, true, false)) * (1 - Lr * s_Lf);
  ee[43] = (R::pnorm(cr3 + cr2 + cr13, 0.0, 1.0, true, false) - R::pnorm(cr2 + cr13, 0.0, 1.0, true, false)) * (1 - Lr * s_Lf);
  ee[44] = R::pnorm(cr13, 0.0, 1.0, true, false) * Lr;
  ee[45] = (R::pnorm(cr2 + cr13, 0.0, 1.0, true, false) - R::pnorm(cr13, 0.0, 1.0, true, false)) * Lr;
  ee[46] = R::pnorm(cr13, 0.0, 1.0, true, false) * (1 - Lr);
  ee[47] = (R::pnorm(cr2 + cr13, 0.0, 1.0, true, false) - R::pnorm(cr13, 0.0, 1.0, true, false)) * (1 - Lr);

  ee = ee * NN;
  
  // Cálculo de la log-verosimilitud acumulada
  double LL = 0.0;
  for (int ii = 0; ii < 48; ii++) {
    if (chd[ii] > 0.0) {
      LL += chd[ii] * (log(chd[ii]) - log(ee[ii]));
    }
  }
  
  // Retornamos el valor de la log-verosimilitud multiplicado por 2
  return 2.0 * LL;
}
