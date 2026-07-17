#ifndef ELECTRON_H
#define ELECTRON_H

#include "TObject.h"

class fCluster : public TObject {
public:
  double x;
  double y;
  double z;
  double t;
  double energy;
  double weight;
  int size;
  ClassDef(fCluster, 1);
};

class fElectron : public TObject {
public:
  double x;
  double y;
  double z;
  double t;
  double energy;

  ClassDef(fElectron, 1);
};

#endif
