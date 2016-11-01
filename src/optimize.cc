#include <iostream>
#include <algorithm>
#include <cmath>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TVectorT.h>

#include "runtime_exception.hh"
#include "timed_counter.hh"

using std::cout;
using std::cerr;
using std::endl;
namespace po = boost::program_options;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

template <typename H>
inline unsigned nbins(const H* h) noexcept {
  return h->GetNbinsX();
}

template <typename H1, typename H2, typename... Hs>
inline unsigned nbins(const H1* h1, const H2* h2, const Hs*... hs) {
  const unsigned n1 = nbins(h1);
  const unsigned n2 = nbins(h2,hs...);
  if (n1!=n2) throw rte(
    "histogram ",h2->GetName()," has ",n2," bins instead of expected ",n1);
  return n1;
}

struct var_sb {
  struct sb {
    double sig, bkg, up;
    sb(): sig(0.), bkg(0.), up(0.) { }
    sb(double sig, double bkg, double up): sig(sig), bkg(bkg), up(up) { }
    inline void operator+=(const sb& x) noexcept {
      sig += x.sig;
      bkg += x.bkg;
      up = x.up;
    }
    inline double signif() const noexcept {
      const double spb = sig + bkg;
      return (__builtin_expect(spb>0.,1) ? sig/sqrt(spb) : 0.);
    }
  };
  std::string name;
  std::vector<sb> bins;
  template <typename Name>
  var_sb(Name&& name, TH1D* sig, TH1D* bkg)
  : name(std::forward<Name>(name)) {
    const unsigned n = nbins(sig,bkg);
    bins.reserve(n);
    auto * sig_arr = sig->GetArray();
    auto * bkg_arr = bkg->GetArray();
    for (unsigned i=1; i<=n; ++i) {
      bins.emplace_back(sig_arr[i],bkg_arr[i],sig->GetBinLowEdge(i+1));
    }
  }
  inline sb& operator[](unsigned i) noexcept { return bins[i]; }
};

int main(int argc, char* argv[])
{
  std::string ifname, ofname, cfname;
  double scaled_lumi;

  // options ---------------------------------------------------
  try {
    po::options_description desc("Options");
    desc.add_options()
      ("input,i", po::value(&ifname)->required(),
       "input root data file")
      ("output,o", po::value(&ofname)->required(),
       "output root file")
      ("conf,c", po::value(&cfname),
       "configuration file")
      ("lumi,l", po::value(&scaled_lumi)->default_value(0.),
       "scale to luminosity of interest")
    ;

    po::positional_options_description pos;
    pos.add("conf",1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
      .options(desc).positional(pos).run(), vm);
    if (argc == 1) {
      cout << desc << endl;
      return 0;
    }
    if (vm.count("conf")) {
      po::store( po::parse_config_file<char>(
        vm["conf"].as<std::string>().c_str(), desc), vm);
    }
    po::notify(vm);
  } catch (std::exception& e) {
    cerr << "\033[31m" << argv[0]
         << " options: " <<  e.what() <<"\033[0m"<< endl;
    return 1;
  }
  // end options ---------------------------------------------------

  auto fin = std::make_unique<TFile>(ifname.c_str(),"read");
  if (fin->IsZombie()) return 1;
  cout << "Input file: " << fin->GetName() << endl;

  auto& mass_range  = *dynamic_cast<TVectorD*>(fin->Get("mass_range"));
  auto& mass_window = *dynamic_cast<TVectorD*>(fin->Get("mass_window"));
  auto& input_lumi  = *dynamic_cast<TVectorD*>(fin->Get("lumi"));

  if (scaled_lumi==0.) scaled_lumi = input_lumi[0];
  const double fl = std::sqrt(scaled_lumi/input_lumi[0]);

  cout << "mass_range  = ("<<mass_range[0]/1e3 <<','<<mass_range[1]/1e3<<") GeV"<<endl;
  cout << "mass_window = ("<<mass_window[0]/1e3 <<','<<mass_window[1]/1e3<<") GeV"<<endl;
  cout << "input  lumi = " << input_lumi[0] << " pb" << endl;
  cout << "scaled lumi = " << scaled_lumi   << " pb" << endl;
  cout << "lumi factor = " << fl << endl;

  std::vector<var_sb> vars;

  {
    TIter next(fin->GetListOfKeys());
    TKey *key;
    std::vector<std::pair<std::string,std::array<TKey*,2>>> hbuff;
    while ((key = static_cast<TKey*>(next()))) {
      std::string name(key->GetName());
      auto sep = name.rfind('_');
      std::array<std::string,2> var {
        name.substr(0,sep), name.substr(sep+1)
      };

      if (!(var[1]=="sig" || var[1]=="bkg")) continue;

      auto rit = std::find_if(hbuff.rbegin(),hbuff.rend(),
        [&name = var[0]](decltype(hbuff)::const_reference x){
          return (x.first == name);
        }
      );
      if (rit==hbuff.rend()) {
        using second_t = decltype(hbuff)::value_type::second_type;
        hbuff.emplace_back(var[0],second_t{nullptr,nullptr});
        rit = hbuff.rbegin();
      }
      rit->second[var[1]=="sig" ? 0 : 1] = key;
    }
    cout << "Variables:";
    for (auto& buff : hbuff) {
      if (buff.second[0]==nullptr || buff.second[1]==nullptr)
        throw rte("missing histograms for variable \'",buff.first,"\'");
      cout << ' ' << buff.first;
      auto hist = [&buff](unsigned i){ // get histogram from key
        return dynamic_cast<TH1D*>(buff.second[i]->ReadObj());
      };
      vars.emplace_back(buff.first,hist(0),hist(1));
    }
  }
  cout << '\n' << endl;

  for (auto& var : vars) {
    cout << "\033[32m" << var.name << "\033[0m" << endl;

    // for (auto& sb : var.bins) cout << sb.signif() << endl;

    unsigned a = 0, b = 1, n = var.bins.size();
    while (b<n) {
      while (fl*var[a].signif()<2. || var[b].sig+var[b].bkg==0.) {
        var[a] += var[b++];
        if (b==n) break;
      }
      cout << fl*var[a].signif() << ' ' << var[a].up << endl;
      a = b;
    }
  }

  return 0;
}
