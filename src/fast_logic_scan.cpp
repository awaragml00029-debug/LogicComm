// src/fast_logic_scan.cpp
// High-performance C++ functions for LogicComm REO matrix computation
#include <Rcpp.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

//' Compute per-cell expression quantile anchors
//'
//' For each cell (column), computes the quantile of expressed (non-zero)
//' gene values. This serves as the dynamic background anchor for REO conversion.
//'
//' @param mat A numeric matrix (genes x cells).
//' @param prob Quantile probability, e.g. 0.5 for median. Default: 0.5.
//' @return Numeric vector of length ncol(mat), one anchor per cell.
// [[Rcpp::export]]
NumericVector cell_quantiles(NumericMatrix mat, double prob = 0.5) {
    int n_genes = mat.nrow();
    int n_cells = mat.ncol();
    NumericVector result(n_cells, 0.0);

    for (int j = 0; j < n_cells; ++j) {
        // Collect expressed values (> 0) for this cell
        std::vector<double> expressed;
        expressed.reserve(n_genes);
        for (int i = 0; i < n_genes; ++i) {
            double v = mat(i, j);
            if (v > 0.0) expressed.push_back(v);
        }
        if (expressed.empty()) {
            result[j] = 0.0;
            continue;
        }
        std::sort(expressed.begin(), expressed.end());
        double idx_f = prob * (expressed.size() - 1);
        int idx_lo  = (int)idx_f;
        int idx_hi  = std::min(idx_lo + 1, (int)expressed.size() - 1);
        double frac = idx_f - idx_lo;
        result[j] = expressed[idx_lo] * (1.0 - frac) + expressed[idx_hi] * frac;
    }
    return result;
}

//' Fast binary logic scan (REO conversion)
//'
//' Converts a gene expression matrix to a logical 0/1 matrix by comparing
//' each cell's gene values against that cell's dynamic anchor.
//'
//' @param mat A numeric matrix (genes x cells).
//' @param anchors A numeric vector of per-cell thresholds (length = ncols).
//' @return A logical matrix (genes x cells); TRUE where expression > anchor.
// [[Rcpp::export]]
LogicalMatrix fast_logic_scan(NumericMatrix mat, NumericVector anchors) {
    int n_genes = mat.nrow();
    int n_cells = mat.ncol();
    if (anchors.size() != n_cells) {
        stop("Length of anchors must equal ncol(mat).");
    }
    LogicalMatrix res(n_genes, n_cells);
    for (int j = 0; j < n_cells; ++j) {
        double thr = anchors[j];
        for (int i = 0; i < n_genes; ++i) {
            res(i, j) = (mat(i, j) > thr);
        }
    }
    return res;
}

//' Compute KNN-based Logic Consensus Score (LCS)
//'
//' For each ligand-receptor pair, counts directed KNN edges where the sender
//' cell expresses the ligand above anchor AND the receiver cell expresses the
//' receptor above anchor.
//'
//' @param ligand_logic Logical vector (length = n_cells): TRUE if cell expresses ligand.
//' @param receptor_logic Logical vector (length = n_cells): TRUE if cell expresses receptor.
//' @param knn_i Integer vector of sender cell indices (1-based), from sparse KNN matrix.
//' @param knn_j Integer vector of receiver cell indices (1-based).
//' @return Numeric scalar: fraction of active KNN directed edges.
// [[Rcpp::export]]
double compute_lcs_knn(LogicalVector ligand_logic,
                       LogicalVector receptor_logic,
                       IntegerVector knn_i,
                       IntegerVector knn_j) {
    int n_edges = knn_i.size();
    if (knn_j.size() != n_edges) {
        stop("knn_i and knn_j must have the same length.");
    }
    int n_cells = ligand_logic.size();
    if (receptor_logic.size() != n_cells) {
        stop("ligand_logic and receptor_logic must have the same length.");
    }
    if (n_edges == 0) return 0.0;

    int active = 0;
    for (int e = 0; e < n_edges; ++e) {
        int si = knn_i[e] - 1;  // 1-based -> 0-based
        int rj = knn_j[e] - 1;
        if (si < 0 || si >= n_cells || rj < 0 || rj >= n_cells) {
            stop("KNN edge index out of bounds.");
        }
        if (ligand_logic[si] == TRUE && receptor_logic[rj] == TRUE) {
            active++;
        }
    }
    return (double)active / n_edges;
}

//' Compute global co-expression Logic Consensus Score (no KNN)
//'
//' Fraction of cells that co-express both ligand and receptor above anchor.
//' This is the "global mode" LCS (no spatial/neighborhood constraint).
//'
//' @param ligand_logic Logical vector: TRUE if cell expresses ligand.
//' @param receptor_logic Logical vector: TRUE if cell expresses receptor.
//' @return Numeric scalar: fraction of cells co-expressing both.
// [[Rcpp::export]]
double compute_lcs_global(LogicalVector ligand_logic,
                          LogicalVector receptor_logic) {
    int n = ligand_logic.size();
    if (receptor_logic.size() != n) {
        stop("ligand_logic and receptor_logic must have the same length.");
    }
    if (n == 0) return 0.0;
    int co_active = 0;
    for (int i = 0; i < n; ++i) {
        if (ligand_logic[i] == TRUE && receptor_logic[i] == TRUE) co_active++;
    }
    return (double)co_active / n;
}
