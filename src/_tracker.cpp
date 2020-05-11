#include "_tracker.h"
#include <Rcpp.h>

Tracker::Tracker(int n_, int frequency_) {
  i = 0;
  n = n_;
  frequency = frequency_;
}

Tracker::~Tracker() {
  std::cerr << std::endl;
}

void Tracker::it() {
  if (i % 100)
    Rcpp::checkUserInterrupt();

  // if ((frequency * i) % n == 0) {
  //   std::cerr << '.';
  // }
  i++;
}