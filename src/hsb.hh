#ifndef SIGNIF_HSB
#define SIGNIF_HSB

#include <vector>
#include <string>
#include <TH1.h>

struct hsb {
  TH1D *sig, *bkg, *tmp;
  hsb(const std::string& name, unsigned nbins, double xmin, double xmax);
  hsb(TH1D* sig, TH1D* bkg);
  ~hsb() { delete tmp; }
  static std::vector<hsb*> all;
};

#endif
