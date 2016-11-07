#include <iostream>
#include <array>
#include <memory>
#include <experimental/optional>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
#include <TKey.h>
#include <TVectorT.h>

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

struct hsb {
  TH1D *sig, *bkg, *tmp, *signif;
  hsb(const std::string& name, const std::vector<double>& bins)
  : sig(new TH1D((name+"_sig").c_str(),"",bins.size()-1,bins.data())),
    bkg(new TH1D((name+"_bkg").c_str(),"",bins.size()-1,bins.data())),
    tmp(new TH1D("","",bins.size()-1,bins.data())),
    signif(new TH1D(name.c_str(),"",bins.size()-1,bins.data()))
  { all.push_back(this); }
  ~hsb() {
    delete tmp;
    delete bkg;
    delete sig;
  }
  static std::vector<hsb*> all;
};
std::vector<hsb*> hsb::all;

struct ifname_t { const std::string *name; bool is_mc; };

void last_bin_incl(TH1* h) noexcept {
  const auto b = h->GetNbinsX(); // last bin
  h->SetBinContent(b,
    h->GetBinContent(b) +
    h->GetBinContent(b+1)
  );
}

void label_jet_bins(TH1* h) {
  TAxis *xa = h->GetXaxis();
  const unsigned n = h->GetNbinsX();

  for (unsigned i=1; ; ++i) {
    std::stringstream ss;
    if (n-i) {
      ss << " = " << i-1;
      xa->SetBinLabel(i,ss.str().c_str());
    } else {
      ss << " #geq " << i-1;
      xa->SetBinLabel(i,ss.str().c_str());
      break;
    }
  }
  xa->SetLabelSize(0.05);
}

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
  const double fl = std::sqrt(lumi_need/lumi_in);

  std::vector<ifname_t> ifname;
  for (const auto& f : ifname_data) ifname.emplace_back(ifname_t{&f,false});
  for (const auto& f : ifname_mc  ) ifname.emplace_back(ifname_t{&f,true });

  TH1::AddDirectory(kFALSE);

  hsb h_total("total",{0.,1.}),
      h_pT_yy("pT_yy",{0,10,15,20,30,45,60,80,115,200,400}),
      h_yAbs_yy("yAbs_yy",{0.0,0.3,0.6,0.9,1.2,1.6,2.4}),
      h_cosTS_yy("cosTS_yy",{0.,0.125,0.250,0.375,0.500,0.625,0.750,0.875,1.0}),
      h_N_j_30("N_j_30",{0, 1, 2, 3, 4}),
      h_pT_j1("pT_j1",{30,40,55,75,120,400}),
      h_Dphi_j_j("Dphi_j_j",{0., 1.0472, 2.0944, 2.61799, 3.14159}),
      h_Dy_j_j("Dy_j_j",{0., 2., 4., 5.5, 8.8}),
      h_m_jj("m_jj",{0,200,400,600,1e3});

      // h_pT_yy("pT_yy",{0,10,15,20,30,45,60,80,115,150,200,260,400}),
      // h_yAbs_yy("yAbs_yy",{0.0,0.15,0.3,0.45,0.6,0.75,0.9,1.2,1.6,2.4}),
      // h_cosTS_yy("cosTS_yy",{0.,0.0625,0.125,0.1875,0.250,0.3125,0.375,0.4375,0.500,0.625,1.0}),
      // h_N_j_30("N_j_30",{0, 1, 2, 3, 4}),
      // h_pT_j1("pT_j1",{30,40,55,75,93,120,167,400}),
      // h_Dphi_j_j("Dphi_j_j",{0., 1.0472, 2.0944, 3.14159}),
      // h_Dy_j_j("Dy_j_j",{0., 2., 4., 5.5, 8.8}),
      // h_m_jj("m_jj",{0,200,400,600,1e3});

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
    TTreeReaderValue<Float_t> _pT_j1   (reader,"HGamEventInfoAuxDyn.pT_j1");
    TTreeReaderValue<Float_t> _Dphi_j_j(reader,"HGamEventInfoAuxDyn.Dphi_j_j");
    TTreeReaderValue<Float_t> _Dy_j_j  (reader,"HGamEventInfoAuxDyn.Dy_j_j");
    TTreeReaderValue<Float_t> _m_jj    (reader,"HGamEventInfoAuxDyn.m_jj");

    double weight = 1.; // Assuming weight in data files is always 1

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

      const auto nj = *N_j_30;

      h_total   .tmp->Fill( 0.5             ,weight);
      h_pT_yy   .tmp->Fill( *_pT_yy/1e3     ,weight);
      h_yAbs_yy .tmp->Fill( *_yAbs_yy       ,weight);
      h_cosTS_yy.tmp->Fill( abs(*_cosTS_yy) ,weight);
      h_N_j_30  .tmp->Fill( nj              ,weight);

      if (nj < 1) continue;
      h_pT_j1   .tmp->Fill( *_pT_j1/1e3     ,weight);

      if (nj < 2) continue;
      h_Dphi_j_j.tmp->Fill( abs(*_Dphi_j_j) ,weight);
      h_Dy_j_j  .tmp->Fill( abs(*_Dy_j_j)   ,weight);
      h_m_jj    .tmp->Fill( *_m_jj/1e3      ,weight);
    }

    for (auto* h : hsb::all) { // merge histograms
      (f.is_mc ? h->sig : h->bkg)->Add(h->tmp,n_all_inv);
      h->tmp->Reset();
    }
  }

  last_bin_incl(h_N_j_30.sig);
  last_bin_incl(h_N_j_30.bkg);
  label_jet_bins(h_N_j_30.signif);

  for (auto* h : hsb::all) {
    h->bkg->Scale(fw);      // scale data by window factor
    h->sig->Scale(lumi_in); // scale MC to lumi

    // calculate significances
    for (unsigned i=0, n=h->sig->GetNbinsX()+2; i<n; ++i) {
      const double sig = h->sig->GetBinContent(i),
                   bkg = h->bkg->GetBinContent(i);
      h->signif->SetBinContent(i, fl*sig/sqrt(sig+bkg) );
    }
  }

  // print expected numbers of events
  for (auto* h : hsb::all) {
    cout << "\033[32m" << h->signif->GetName()  << "\033[0m\n";
    for (int i=1, n=h->sig->GetNbinsX()+2; i<n; ++i) {
      cout << "[" << h->sig->GetBinLowEdge(i) << ',';
      if (i==n-1) cout << "\\infty";
      else cout << h->sig->GetBinLowEdge(i+1);
      cout << ") "
           << h->sig->GetBinContent(i) << ' '
           << h->bkg->GetBinContent(i) << ' '
           << h->signif->GetBinContent(i) << endl;
    }
    cout << endl;
  }

  h_N_j_30.signif->SetBinContent(5,0);

  auto fout = std::make_unique<TFile>(ofname.c_str(),"recreate");
  if (fout->IsZombie()) return 1;
  cout << "\nOutput file: " << fout->GetName() << endl;
  fout->cd();

  std::stringstream ss;
  ss << "Background is estimated from data, "
    "using events inside mass range, but outside mass window.\n"
    "Signal is estimated from MC, using events inside mass window.\n"
    "signif = (fl×s)/sqrt(s+fw×b)\n"
    "fl = sqrt(wanted lumi/data lumi)\n"
    "fw = |mass window|/(|mass range|-|mass window|)\n"
    "The signal is written for fl = 1.\n"
    "The background is written multiplied by fw.";
  TNamed("note", ss.str().c_str()).Write();

  for (auto* h : hsb::all) h->signif->SetDirectory(fout.get());

  auto write_d = [](const char* name, const std::vector<double>& d){
    TVectorD v(d.size(),d.data());
    v.Write(name);
  };

  write_d("lumi_in",{lumi_in});
  write_d("lumi_need",{lumi_need});
  write_d("mass_range",{mass_range[0],mass_range[1]});
  write_d("mass_window",{mass_window[0],mass_window[1]});

  ss.str({});
  for (const auto& f : ifname_data) ss << f << '\n';
  for (const auto& f : ifname_mc  ) ss << f << '\n';
  TNamed("input_files", ss.str().c_str()).Write();

  fout->Write();

  cout << "\033[32m Total SIG:\033[0m "
       << h_total.sig->GetBinContent(1) << endl;
  cout << "\033[32m Total BKG:\033[0m "
       << h_total.bkg->GetBinContent(1) << endl;

  return 0;
}
