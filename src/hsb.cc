#include "hsb.hh"

std::vector<hsb*> hsb::all;

hsb::hsb(const std::string& name, unsigned nbins, double xmin, double xmax)
: sig(new TH1D((name+"_sig").c_str(),"",nbins,xmin,xmax)),
  bkg(new TH1D((name+"_bkg").c_str(),"",nbins,xmin,xmax)),
  tmp(new TH1D("","",nbins,xmin,xmax))
{
  all.push_back(this);
}

hsb::hsb(TH1D* sig, TH1D* bkg): sig(sig), bkg(bkg) {
  all.push_back(this);
}
