#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double MAJ(NumericVector Q0, NumericVector chd, NumericVector NN) {
  
  // Aplicamos pnorm a Q0 para obtener el vector Q
  NumericVector Q = Rcpp::pnorm(Q0);
  
  // Definimos las variables basadas en Q
  double cr11 = Q[0];  // Umbral cr11 para la condición 1
  double cr12 = Q[1];  // Umbral cr12 para la condición 2
  double cr13 = Q[2];  // Umbral cr13 para la condición 3
  double d = Q[3];     // Desplazamiento d
  double ss = Q[4];    // Desviación estándar ss
  double c = Q[5];     // Parámetro c
  double r = Q[6];     // Parámetro r
  double Lf = Q[7];    // Parámetro Lf
  double Hsf = Q[8];   // Parámetro Hsf
  double Hss = Q[9];   // Parámetro Hss
  double g3 = Q[10];   // Parámetro g3
  double s_g2 = Q[11]; // Parámetro s_g2
  double s_g1 = Q[12]; // Parámetro s_g1
  double s_Ls = Q[13]; // Parámetro s_Ls
  double s_Lgo = Q[14]; // Parámetro s_Lgo
  double s_Hff = Q[15]; // Parámetro s_Hff
  double s_Hgf = Q[16]; // Parámetro s_Hgf
  double s_Hgs = Q[17]; // Parámetro s_Hgs
  double s_Hfs = Q[18]; // Parámetro s_Hfs
  double s_Lgn = Q[19]; // Parámetro s_Lgn
  
  Rcpp::NumericVector ee(48);  // Para almacenar 16 valores (Hit, Miss, FA, CR) para Fast/Slow y High/Low en Condition 1
  
  // Condición 1 - Hit Fast High
  ee[0] = ((1 - R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false)) * Lf * (Hsf * s_Hff)) +
    ((R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false)) * r * ((Lf * s_Ls) * Hsf)) +
    ((R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (g3 * s_g2 * s_g1) * ((Lf * s_Ls * s_Lgo)) * (Hsf * s_Hgf));
  ee[1] = ((1 - R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false)) * Lf * (1 - (Hsf * s_Hff))) +
    ((R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false)) * r * ((Lf * s_Ls) * (1 - Hsf))) +
    ((R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (g3 * s_g2 * s_g1) * ((Lf * s_Ls * s_Lgo)) * (1 - (Hsf * s_Hgf)));
  ee[2] = ((1 - R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false)) * (1 - Lf) * (Hss * s_Hfs)) +
    ((R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false)) * r * (1 - (Lf * s_Ls)) * Hss) +
    ((R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (g3 * s_g2 * s_g1) * (1 - (Lf * s_Ls * s_Lgo)) * (Hss * s_Hgs));
  ee[3] = ((1 - R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false)) * (1 - Lf) * (1 - (Hss * s_Hfs))) +
    ((R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false)) * r * (1 - (Lf * s_Ls)) * (1 - Hss)) +
    ((R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (g3 * s_g2 * s_g1) * (1 - (Lf * s_Ls * s_Lgo)) * (1 - (Hss * s_Hgs)));
  ee[4] = ((R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (1 - (g3 * s_g2 * s_g1)) * ((Lf * s_Ls * s_Lgn)) * (Hsf * s_Hgf)) +
    (R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false) * Lf * (Hsf * s_Hff));
  ee[5] = ((R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (1 - (g3 * s_g2 * s_g1)) * ((Lf * s_Ls * s_Lgn)) * (1 - (Hsf * s_Hgf))) +
    (R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false) * Lf * (1 - (Hsf * s_Hff)));
  ee[6] = ((R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (1 - (g3 * s_g2 * s_g1)) * (1 - (Lf * s_Ls * s_Lgn)) * (Hss * s_Hgs)) +
    (R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false) * (1 - Lf) * (Hss * s_Hfs));
  ee[7] = ((R::pnorm((cr11 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (1 - (g3 * s_g2 * s_g1)) * (1 - (Lf * s_Ls * s_Lgn)) * (1 - (Hss * s_Hgs))) +
    (R::pnorm((cr11 - d) / ss, 0.0, 1.0, true, false) * (1 - Lf) * (1 - (Hss * s_Hfs)));
  
  ee[8] = (R::pnorm(cr11 + c, 0.0, 1.0, true, false) - R::pnorm(cr11, 0.0, 1.0, true, false)) * (g3 * s_g2 * s_g1) * (Lf * s_Ls * s_Lgo) * (Hsf * s_Hgf) +
    ((1 - R::pnorm(cr11 + c, 0.0, 1.0, true, false)) * Lf * (Hsf * s_Hff));
  ee[9] = (R::pnorm(cr11 + c, 0.0, 1.0, true, false) - R::pnorm(cr11, 0.0, 1.0, true, false)) * (g3 * s_g2 * s_g1) * (Lf * s_Ls * s_Lgo) * (1 - (Hsf * s_Hgf)) +
    ((1 - R::pnorm(cr11 + c, 0.0, 1.0, true, false)) * Lf * (1 - (Hsf * s_Hff)));
  ee[10] = (R::pnorm(cr11 + c, 0.0, 1.0, true, false) - R::pnorm(cr11, 0.0, 1.0, true, false)) * (g3 * s_g2 * s_g1) * (1 - (Lf * s_Ls * s_Lgo)) * (Hss * s_Hgs) +
    ((1 - R::pnorm(cr11 + c, 0.0, 1.0, true, false)) * (1 - Lf) * (Hss * s_Hfs));
  ee[11] = (R::pnorm(cr11 + c, 0.0, 1.0, true, false) - R::pnorm(cr11, 0.0, 1.0, true, false)) * (g3 * s_g2 * s_g1) * (1 - (Lf * s_Ls * s_Lgo)) * (1 - (Hss * s_Hgs)) +
    ((1 - R::pnorm(cr11 + c, 0.0, 1.0, true, false)) * (1 - Lf) * (1 - (Hss * s_Hfs)));
  ee[12] = R::pnorm(cr11, 0.0, 1.0, true, false) * Lf * (Hsf * s_Hff) +
    ((R::pnorm(cr11 + c, 0.0, 1.0, true, false) - R::pnorm(cr11, 0.0, 1.0, true, false)) * (1 - (g3 * s_g2 * s_g1)) * (Lf * s_Ls * s_Lgn) * (Hsf * s_Hgf));
  ee[13] = R::pnorm(cr11, 0.0, 1.0, true, false) * Lf * (1 - (Hsf * s_Hff)) +
    ((R::pnorm(cr11 + c, 0.0, 1.0, true, false) - R::pnorm(cr11, 0.0, 1.0, true, false)) * (1 - (g3 * s_g2 * s_g1)) * (Lf * s_Ls * s_Lgn) * (1 - (Hsf * s_Hgf)));
  ee[14] = R::pnorm(cr11, 0.0, 1.0, true, false) * (1 - Lf) * (Hss * s_Hfs) +
    ((R::pnorm(cr11 + c, 0.0, 1.0, true, false) - R::pnorm(cr11, 0.0, 1.0, true, false)) * (1 - (g3 * s_g2 * s_g1)) * (1 - (Lf * s_Ls * s_Lgn)) * (Hss * s_Hgs));
  ee[15] = R::pnorm(cr11, 0.0, 1.0, true, false) * (1 - Lf) * (1 - (Hss * s_Hfs)) +
    ((R::pnorm(cr11 + c, 0.0, 1.0, true, false) - R::pnorm(cr11, 0.0, 1.0, true, false)) * (1 - (g3 * s_g2 * s_g1)) * (1 - (Lf * s_Ls * s_Lgn)) * (1 - (Hss * s_Hgs)));
  
  // Condición 2
  ee[16] = ((1 - R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false)) * Lf * (Hsf * s_Hff)) +
    ((R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false)) * r * ((Lf * s_Ls) * Hsf)) +
    ((R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (g3 * s_g2 * s_g1) * ((Lf * s_Ls * s_Lgo)) * (Hsf * s_Hgf));
  ee[17] = ((1 - R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false)) * Lf * (1 - (Hsf * s_Hff))) +
    ((R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false)) * r * ((Lf * s_Ls) * (1 - Hsf))) +
    ((R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (g3 * s_g2 * s_g1) * ((Lf * s_Ls * s_Lgo)) * (1 - (Hsf * s_Hgf)));
  ee[18] = ((1 - R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false)) * (1 - Lf) * (Hss * s_Hfs)) +
    ((R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false)) * r * (1 - (Lf * s_Ls)) * Hss) +
    ((R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (g3 * s_g2 * s_g1) * (1 - (Lf * s_Ls * s_Lgo)) * (Hss * s_Hgs));
  ee[19] = ((1 - R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false)) * (1 - Lf) * (1 - (Hss * s_Hfs))) +
    ((R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false)) * r * (1 - (Lf * s_Ls)) * (1 - Hss)) +
    ((R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (g3 * s_g2 * s_g1) * (1 - (Lf * s_Ls * s_Lgo)) * (1 - (Hss * s_Hgs)));
  ee[20] = ((R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (1 - (g3 * s_g2 * s_g1)) * ((Lf * s_Ls * s_Lgn)) * (Hsf * s_Hgf)) +
    (R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false) * Lf * (Hsf * s_Hff));
  ee[21] = ((R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (1 - (g3 * s_g2 * s_g1)) * ((Lf * s_Ls * s_Lgn)) * (1 - (Hsf * s_Hgf))) +
    (R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false) * Lf * (1 - (Hsf * s_Hff)));
  ee[22] = ((R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (1 - (g3 * s_g2 * s_g1)) * (1 - (Lf * s_Ls * s_Lgn)) * (Hss * s_Hgs)) +
    (R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false) * (1 - Lf) * (Hss * s_Hfs));
  ee[23] = ((R::pnorm((cr12 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (1 - (g3 * s_g2 * s_g1)) * (1 - (Lf * s_Ls * s_Lgn)) * (1 - (Hss * s_Hgs))) +
    (R::pnorm((cr12 - d) / ss, 0.0, 1.0, true, false) * (1 - Lf) * (1 - (Hss * s_Hfs)));
  
  ee[24] = (R::pnorm(cr12 + c, 0.0, 1.0, true, false) - R::pnorm(cr12, 0.0, 1.0, true, false)) * (g3 * s_g2 * s_g1) * (Lf * s_Ls * s_Lgo) * (Hsf * s_Hgf) +
    ((1 - R::pnorm(cr12 + c, 0.0, 1.0, true, false)) * Lf * (Hsf * s_Hff));
  ee[25] = (R::pnorm(cr12 + c, 0.0, 1.0, true, false) - R::pnorm(cr12, 0.0, 1.0, true, false)) * (g3 * s_g2 * s_g1) * (Lf * s_Ls * s_Lgo) * (1 - (Hsf * s_Hgf)) +
    ((1 - R::pnorm(cr12 + c, 0.0, 1.0, true, false)) * Lf * (1 - (Hsf * s_Hff)));
  ee[26] = (R::pnorm(cr12 + c, 0.0, 1.0, true, false) - R::pnorm(cr12, 0.0, 1.0, true, false)) * (g3 * s_g2 * s_g1) * (1 - (Lf * s_Ls * s_Lgo)) * (Hss * s_Hgs) +
    ((1 - R::pnorm(cr12 + c, 0.0, 1.0, true, false)) * (1 - Lf) * (Hss * s_Hfs));
  ee[27] = (R::pnorm(cr12 + c, 0.0, 1.0, true, false) - R::pnorm(cr12, 0.0, 1.0, true, false)) * (g3 * s_g2 * s_g1) * (1 - (Lf * s_Ls * s_Lgo)) * (1 - (Hss * s_Hgs)) +
    ((1 - R::pnorm(cr12 + c, 0.0, 1.0, true, false)) * (1 - Lf) * (1 - (Hss * s_Hfs)));
  ee[28] = R::pnorm(cr12, 0.0, 1.0, true, false) * Lf * (Hsf * s_Hff) +
    ((R::pnorm(cr12 + c, 0.0, 1.0, true, false) - R::pnorm(cr12, 0.0, 1.0, true, false)) * (1 - (g3 * s_g2 * s_g1)) * (Lf * s_Ls * s_Lgn) * (Hsf * s_Hgf));
  ee[29] = R::pnorm(cr12, 0.0, 1.0, true, false) * Lf * (1 - (Hsf * s_Hff)) +
    ((R::pnorm(cr12 + c, 0.0, 1.0, true, false) - R::pnorm(cr12, 0.0, 1.0, true, false)) * (1 - (g3 * s_g2 * s_g1)) * (Lf * s_Ls * s_Lgn) * (1 - (Hsf * s_Hgf)));
  ee[30] = R::pnorm(cr12, 0.0, 1.0, true, false) * (1 - Lf) * (Hss * s_Hfs) +
    ((R::pnorm(cr12 + c, 0.0, 1.0, true, false) - R::pnorm(cr12, 0.0, 1.0, true, false)) * (1 - (g3 * s_g2 * s_g1)) * (1 - (Lf * s_Ls * s_Lgn)) * (Hss * s_Hgs));
  ee[31] = R::pnorm(cr12, 0.0, 1.0, true, false) * (1 - Lf) * (1 - (Hss * s_Hfs)) +
    ((R::pnorm(cr12 + c, 0.0, 1.0, true, false) - R::pnorm(cr12, 0.0, 1.0, true, false)) * (1 - (g3 * s_g2 * s_g1)) * (1 - (Lf * s_Ls * s_Lgn)) * (1 - (Hss * s_Hgs)));
  
  // Condición 3
  ee[32] = ((1 - R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false)) * Lf * (Hsf * s_Hff)) +
    ((R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false)) * r * ((Lf * s_Ls) * Hsf)) +
    ((R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (g3 * s_g2 * s_g1) * ((Lf * s_Ls * s_Lgo)) * (Hsf * s_Hgf));
  ee[33] = ((1 - R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false)) * Lf * (1 - (Hsf * s_Hff))) +
    ((R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false)) * r * ((Lf * s_Ls) * (1 - Hsf))) +
    ((R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (g3 * s_g2 * s_g1) * ((Lf * s_Ls * s_Lgo)) * (1 - (Hsf * s_Hgf)));
  ee[34] = ((1 - R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false)) * (1 - Lf) * (Hss * s_Hfs)) +
    ((R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false)) * r * (1 - (Lf * s_Ls)) * Hss) +
    ((R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (g3 * s_g2 * s_g1) * (1 - (Lf * s_Ls * s_Lgo)) * (Hss * s_Hgs));
  ee[35] = ((1 - R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false)) * (1 - Lf) * (1 - (Hss * s_Hfs))) +
    ((R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false)) * r * (1 - (Lf * s_Ls)) * (1 - Hss)) +
    ((R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (g3 * s_g2 * s_g1) * (1 - (Lf * s_Ls * s_Lgo)) * (1 - (Hss * s_Hgs)));
  ee[36] = ((R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (1 - (g3 * s_g2 * s_g1)) * ((Lf * s_Ls * s_Lgn)) * (Hsf * s_Hgf)) +
    (R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false) * Lf * (Hsf * s_Hff));
  ee[37] = ((R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (1 - (g3 * s_g2 * s_g1)) * ((Lf * s_Ls * s_Lgn)) * (1 - (Hsf * s_Hgf))) +
    (R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false) * Lf * (1 - (Hsf * s_Hff)));
  ee[38] = ((R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (1 - (g3 * s_g2 * s_g1)) * (1 - (Lf * s_Ls * s_Lgn)) * (Hss * s_Hgs)) +
    (R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false) * (1 - Lf) * (Hss * s_Hfs));
  ee[39] = ((R::pnorm((cr13 + c - d) / ss, 0.0, 1.0, true, false) - R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false)) * (1 - r) * (1 - (g3 * s_g2 * s_g1)) * (1 - (Lf * s_Ls * s_Lgn)) * (1 - (Hss * s_Hgs))) +
    (R::pnorm((cr13 - d) / ss, 0.0, 1.0, true, false) * (1 - Lf) * (1 - (Hss * s_Hfs)));
  
  ee[40] = (R::pnorm(cr13 + c, 0.0, 1.0, true, false) - R::pnorm(cr13, 0.0, 1.0, true, false)) * (g3 * s_g2 * s_g1) * (Lf * s_Ls * s_Lgo) * (Hsf * s_Hgf) +
    ((1 - R::pnorm(cr13 + c, 0.0, 1.0, true, false)) * Lf * (Hsf * s_Hff));
  ee[41] = (R::pnorm(cr13 + c, 0.0, 1.0, true, false) - R::pnorm(cr13, 0.0, 1.0, true, false)) * (g3 * s_g2 * s_g1) * (Lf * s_Ls * s_Lgo) * (1 - (Hsf * s_Hgf)) +
    ((1 - R::pnorm(cr13 + c, 0.0, 1.0, true, false)) * Lf * (1 - (Hsf * s_Hff)));
  ee[42] = (R::pnorm(cr13 + c, 0.0, 1.0, true, false) - R::pnorm(cr13, 0.0, 1.0, true, false)) * (g3 * s_g2 * s_g1) * (1 - (Lf * s_Ls * s_Lgo)) * (Hss * s_Hgs) +
    ((1 - R::pnorm(cr13 + c, 0.0, 1.0, true, false)) * (1 - Lf) * (Hss * s_Hfs));
  ee[43] = (R::pnorm(cr13 + c, 0.0, 1.0, true, false) - R::pnorm(cr13, 0.0, 1.0, true, false)) * (g3 * s_g2 * s_g1) * (1 - (Lf * s_Ls * s_Lgo)) * (1 - (Hss * s_Hgs)) +
    ((1 - R::pnorm(cr13 + c, 0.0, 1.0, true, false)) * (1 - Lf) * (1 - (Hss * s_Hfs)));
  ee[44] = R::pnorm(cr13, 0.0, 1.0, true, false) * Lf * (Hsf * s_Hff) +
    ((R::pnorm(cr13 + c, 0.0, 1.0, true, false) - R::pnorm(cr13, 0.0, 1.0, true, false)) * (1 - (g3 * s_g2 * s_g1)) * (Lf * s_Ls * s_Lgn) * (Hsf * s_Hgf));
  ee[45] = R::pnorm(cr13, 0.0, 1.0, true, false) * Lf * (1 - (Hsf * s_Hff)) +
    ((R::pnorm(cr13 + c, 0.0, 1.0, true, false) - R::pnorm(cr13, 0.0, 1.0, true, false)) * (1 - (g3 * s_g2 * s_g1)) * (Lf * s_Ls * s_Lgn) * (1 - (Hsf * s_Hgf)));
  ee[46] = R::pnorm(cr13, 0.0, 1.0, true, false) * (1 - Lf) * (Hss * s_Hfs) +
    ((R::pnorm(cr13 + c, 0.0, 1.0, true, false) - R::pnorm(cr13, 0.0, 1.0, true, false)) * (1 - (g3 * s_g2 * s_g1)) * (1 - (Lf * s_Ls * s_Lgn)) * (Hss * s_Hgs));
  ee[47] = R::pnorm(cr13, 0.0, 1.0, true, false) * (1 - Lf) * (1 - (Hss * s_Hfs)) +
    ((R::pnorm(cr13 + c, 0.0, 1.0, true, false) - R::pnorm(cr13, 0.0, 1.0, true, false)) * (1 - (g3 * s_g2 * s_g1)) * (1 - (Lf * s_Ls * s_Lgn)) * (1 - (Hss * s_Hgs)));
  
  
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
