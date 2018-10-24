#ifndef TRACKER_H
#define TRACKER_H

class Tracker {
private:
  int i, n, frequency;
public:
  Tracker(int n_, int frequency_ = 80);
  ~Tracker();
  void it();
};

#endif