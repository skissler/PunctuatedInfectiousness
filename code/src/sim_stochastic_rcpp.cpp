#include <Rcpp.h>
#include <queue>
#include <vector>
#include <functional>

// Use R's sample.int for RNG compatibility with the pure-R implementation
static Rcpp::Function r_sample_int("sample.int");

// [[Rcpp::export]]
Rcpp::NumericVector sim_stochastic_rcpp(int n, Rcpp::Function gen_inf_attempts) {

  // Initialize infection times to Inf
  Rcpp::NumericVector tinf_vec(n, R_PosInf);

  // Min-heap priority queue for event times
  std::priority_queue<double, std::vector<double>, std::greater<double>> pq;

  // Seed infection: pick random index case using R's sample.int (1-based)
  int indexcase = Rcpp::as<int>(r_sample_int(n, 1)) - 1; // 0-based
  tinf_vec[indexcase] = 0.0;

  // Generate initial infection attempts from index case
  Rcpp::NumericVector init_attempts = gen_inf_attempts(0.0);
  for (int i = 0; i < init_attempts.size(); i++) {
    pq.push(init_attempts[i]);
  }

  // Process events in chronological order (min-heap gives smallest first)
  while (!pq.empty()) {
    double t_attempt = pq.top();
    pq.pop();

    // Pick a random target using R's sample.int (1-based -> 0-based)
    int target = Rcpp::as<int>(r_sample_int(n, 1)) - 1;

    if (tinf_vec[target] == R_PosInf) {
      // Successful infection
      tinf_vec[target] = t_attempt;

      // Generate new infection attempts from this newly infected person
      Rcpp::NumericVector new_attempts = gen_inf_attempts(t_attempt);
      for (int i = 0; i < new_attempts.size(); i++) {
        pq.push(new_attempts[i]);
      }
    }
  }

  return tinf_vec;
}
