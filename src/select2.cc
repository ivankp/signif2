#include <iostream>
#include <iomanip>
#include <array>
#include <cmath>
#include <memory>
#include <initializer_list>
#include <experimental/optional>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TVectorT.h>
#include <TLorentzVector.h>

#include "runtime_exception.hh"
#define IVANP_ARRAY_BOOST_PO
#include "array_ops.hh"
#include "timed_counter.hh"

using std::cout;
using std::cerr;
using std::endl;
namespace po = boost::program_options;

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

inline TH1D* make_hist(const std::string& name,
  unsigned nbins, double min, double max
) {
  return new TH1D(name.c_str(),"",nbins,min,max);
}

double weight;

struct hsb {
  static constexpr unsigned nbins = 10000;

  TH1D *sig, *bkg, *tmp, *signif, *purity;
  hsb(const std::string& name, double min, double max)
  : sig(make_hist(name+"_sig",nbins,min,max)),
    bkg(make_hist(name+"_bkg",nbins,min,max)),
    tmp(make_hist({},nbins,min,max))
  { all.push_back(this); }
  ~hsb() {
    delete tmp;
  }
  static std::vector<hsb*> all;

  inline void operator()(double x) { tmp->Fill(x,weight); }

  void set_dir(TDirectory *dir) {
    sig->SetDirectory(dir);
    bkg->SetDirectory(dir);
  }
};
std::vector<hsb*> hsb::all;

std::ostream& operator<<(std::ostream& o, const hsb& h) {
  const auto prec = o.precision();
  const std::ios::fmtflags f( o.flags() );
  std::string name(h.signif->GetName());
  o << "\033[32m" << name.substr(0,name.rfind('_'))  << "\033[0m\n";

  const int n = h.sig->GetNbinsX();
  const auto *_sig = dynamic_cast<TArrayD*>(h.sig)->GetArray();
  const auto *_bkg = dynamic_cast<TArrayD*>(h.bkg)->GetArray();

  int edge = 1;
  double s = 0., b = 0.;

  for (int i=1; i<=n; ++i) {
    s += _sig[i];
    b += _bkg[i];
    const double signif = s/sqrt(s+b);
    if (signif > 2. && i == n) {
      o << "\033[35m[" << h.sig->GetBinLowEdge(edge) << ','
        << h.sig->GetBinLowEdge(i+1) << ")\033[0m "
        << round(s) << ' ' << round(b) << ' '
        << std::setprecision(2) << std::fixed
        << signif
        << std::setprecision(prec) << endl;
      o.flags( f );

      edge = i+1;
      s = b = 0.;
    }
  }
  return o;
}

struct ifname_t { const std::string *name; bool is_mc; };

int main(int argc, char* argv[])
{
  std::vector<std::string> ifname_data, ifname_mc;
  std::string ofname, cfname;
  double lumi_in, lumi_need;
  std::array<double,2> mass_range, mass_window;

  // options ---------------------------------------------------
  try {
    po::options_description desc("Options");
    desc.add_options()
      ("data", po::value(&ifname_data)->multitoken()->required(),
       "input root data files")
      ("mc", po::value(&ifname_mc)->multitoken()->required(),
       "input root Monte Carlo files")
      ("output,o", po::value(&ofname)->required(),
       "output root file")
      ("conf,c", po::value(&cfname),
       "configuration file")
      ("lumi.in,l", po::value(&lumi_in)->required(),
       "luminosity of input data")
      ("lumi.need,l", po::value(&lumi_need)->required(),
       "target luminosity")
      ("mass-range,r", po::value(&mass_range)
       ->default_value({105.,160.},"(105,160) GeV"),
       "select mass range (data)")
      ("mass-window,w", po::value(&mass_window)
       ->default_value({121.,129.},"(121,129) GeV"),
       "select mass window (mc)")
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
  mass_range  *= 1e3;
  mass_window *= 1e3;

  const double fw = len(mass_window)/(len(mass_range)-len(mass_window));
  const double fl = lumi_need/lumi_in;

  std::vector<ifname_t> ifname;
  for (const auto& f : ifname_data) ifname.emplace_back(ifname_t{&f,false});
  for (const auto& f : ifname_mc  ) ifname.emplace_back(ifname_t{&f,true });

  TH1::AddDirectory(kFALSE);

  hsb
    h_pT_yy("pT_yy",0.,400.),
    h_yAbs_yy("yAbs_yy",0.,2.4),
    h_cosTS_yy("cosTS_yy",0.,1.),
    // h_pTt_yy("pTt_yy",),
    h_Dy_y_y("Dy_y_y",0.,8.8),

    h_pT_j1("pT_j1",0.,400.), h_yAbs_j1("yAbs_j1",0.,2.4),
    h_pT_j2("pT_j2",0.,400.), h_yAbs_j2("yAbs_j2",0.,2.4),
    h_pT_j3("pT_j3",0.,400.), h_yAbs_j3("yAbs_j3",0.,2.4),

    h_m_jj("m_jj",0.,1e3),
    h_HT("HT",0.,1e3),
    h_Dy_j_j("Dy_j_j",0.,8.8),
    h_Dphi_j_j("Dphi_j_j",0.,M_PI),

    h_pT_yyjj("pT_yyjj",0.,1e3),
    h_Dphi_yy_jj("Dphi_yy_jj",0.,M_PI);

    // h_tau("tau",),
    // h_tau1("tau1",);

  std::array<TLorentzVector,3> jets;

  for (const auto& f : ifname) { // loop over input files
    auto file = std::make_unique<TFile>(f.name->c_str(),"read");
    if (file->IsZombie()) return 1;

    cout << ( f.is_mc ? "MC:" : "Data:" ) << ' '
         << file->GetName() << endl;

    double n_all_inv = 1;
    if (f.is_mc) {
      TIter next(file->GetListOfKeys());
      TKey *key;
      while ((key = static_cast<TKey*>(next()))) {
        std::string name(key->GetName());
        if (name.substr(0,8)!="CutFlow_" ||
            name.substr(name.size()-18)!="_noDalitz_weighted") continue;
        TH1 *h = static_cast<TH1*>(key->ReadObj());
        cout << h->GetName() << endl;
        n_all_inv = h->GetBinContent(3);
        cout << h->GetXaxis()->GetBinLabel(3) << " = " << n_all_inv << endl;
        n_all_inv = 1./n_all_inv;
        break;
      }
    }

    TTreeReader reader("CollectionTree",file.get());
    std::experimental::optional<TTreeReaderValue<Float_t>> cs_br_fe, _weight;
    if (f.is_mc) {
      cs_br_fe.emplace(reader,"HGamEventInfoAuxDyn.crossSectionBRfilterEff");
      _weight .emplace(reader,"HGamEventInfoAuxDyn.weight");
    }
    TTreeReaderValue<Char_t>  isPassed (reader,"HGamEventInfoAuxDyn.isPassed");
    TTreeReaderValue<Int_t>   N_j_30   (reader,"HGamEventInfoAuxDyn.N_j_30");
    TTreeReaderValue<Float_t> _m_yy    (reader,"HGamEventInfoAuxDyn.m_yy");
    TTreeReaderValue<Float_t> _pT_yy   (reader,"HGamEventInfoAuxDyn.pT_yy");
    TTreeReaderValue<Float_t> _yAbs_yy (reader,"HGamEventInfoAuxDyn.yAbs_yy");
    TTreeReaderValue<Float_t> _cosTS_yy(reader,"HGamEventInfoAuxDyn.cosTS_yy");
    TTreeReaderValue<Float_t> _Dy_y_y  (reader,"HGamEventInfoAuxDyn.Dy_y_y");
    TTreeReaderValue<Float_t> _pT_j1   (reader,"HGamEventInfoAuxDyn.pT_j1");
    TTreeReaderValue<Float_t> _Dphi_j_j(reader,"HGamEventInfoAuxDyn.Dphi_j_j");
    TTreeReaderValue<Float_t> _Dy_j_j  (reader,"HGamEventInfoAuxDyn.Dy_j_j");
    TTreeReaderValue<Float_t> _m_jj    (reader,"HGamEventInfoAuxDyn.m_jj");

    TTreeReaderArray<Float_t>
      _pt_j (reader,"HGamAntiKt4EMTopoJetsAuxDyn.pt"),
      _eta_j(reader,"HGamAntiKt4EMTopoJetsAuxDyn.eta"),
      _phi_j(reader,"HGamAntiKt4EMTopoJetsAuxDyn.phi"),
      _m_j  (reader,"HGamAntiKt4EMTopoJetsAuxDyn.m");

    // weight is a global variable
    weight = 1.; // Assuming weight in data files is always 1

    using tc = ivanp::timed_counter<Long64_t>;
    for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
      if (!*isPassed) continue;

      const double m_yy = *_m_yy;
      if (!in(m_yy,mass_range)) continue;

      if (f.is_mc) { // signal from MC
        if (!in(m_yy,mass_window)) continue;
        weight = (**_weight)*(**cs_br_fe);
      } else { // background from data
        if (in(m_yy,mass_window)) continue;
      }

      // FILL HISTOGRAMS =====================================================
    // h_pT_yy("pTt_yy",),

    // h_pT_j1("pT_j1",400.), h_yAbs_j1("yAbs_j1",2.4),
    // h_pT_j2("pT_j2",400.), h_yAbs_j2("yAbs_j2",2.4),
    // h_pT_j3("pT_j3",400.), h_yAbs_j3("yAbs_j3",2.4),
    //
    // h_m_jj("m_jj",1e3),
    // h_HT("HT",1e3),
    // h_Dy_j_j("Dy_j_j",8.8),
    // h_Dphi_j_j("Dphi_j_j",M_PI),
    //
    // h_pT_yyjj("pT_yyjj",1e3),
    // h_pT_yyjj("Dphi_yy_jj",M_PI);

      const auto pT_yy = *_pT_yy/1e3;
      const auto yAbs_yy = *_yAbs_yy;
      const auto cosTS_yy = abs(*_cosTS_yy);
      const auto Dy_y_y = abs(*_Dy_y_y);

      h_pT_yy( pT_yy ); h_yAbs_yy( yAbs_yy ); h_cosTS_yy( cosTS_yy );
      h_Dy_y_y( Dy_y_y );

      const unsigned nj = *N_j_30;

      if (nj < 1) continue;

      for (unsigned i=0; i<std::min(nj,3u); ++i) {
        jets[i].SetPtEtaPhiM(_pt_j[i],_eta_j[i],_phi_j[i],_m_j[i]);
      }

      const auto pT_j1 = *_pT_j1/1e3;
      // const auto pT_j2 = *_pT_j2/1e3;

      if (jets[0].Pt() != pT_j1) test( (jets[0].Pt()-pT_j1) )

      // h_pT_j1( pT_j1 ); h_yAbs_j1( yAbs_j1 );


    }

    for (auto* h : hsb::all) { // merge histograms
      (f.is_mc ? h->sig : h->bkg)->Add(h->tmp,n_all_inv);
      h->tmp->Reset();
    }
  }

  for (auto* h : hsb::all) {
    h->bkg->Scale(fw*fl);      // scale data by window factor
    h->sig->Scale(lumi_in*fl); // scale MC to lumi
    cout << *h << endl;
  }

  auto fout = std::make_unique<TFile>(ofname.c_str(),"recreate");
  if (fout->IsZombie()) return 1;
  cout << "\nOutput file: " << fout->GetName() << endl;
  fout->cd();

  std::stringstream ss;
  for (const auto& f : ifname_data) ss << f << '\n';
  for (const auto& f : ifname_mc  ) ss << f << '\n';
  TNamed("input_files", ss.str().c_str()).Write();

  fout->Write();

  return 0;
}
