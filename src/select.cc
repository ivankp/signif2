#include <iostream>
#include <array>
#include <experimental/optional>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1.h>
#include <TKey.h>

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

struct ifname_t { const std::string *name; bool is_mc; };

struct hsb {
  TH1D *sig, *bkg, *tmp;
  hsb(const std::string& name, unsigned nbins, double xmin, double xmax)
  : sig(new TH1D((name+"_sig").c_str(),"",nbins,xmin,xmax)),
    bkg(new TH1D((name+"_bkg").c_str(),"",nbins,xmin,xmax)),
    tmp(new TH1D("","",nbins,xmin,xmax))
  {
    all.push_back(this);
  }
  //~hsb() { delete tmp; }
  static std::vector<hsb*> all;
};
std::vector<hsb*> hsb::all;

int main(int argc, char* argv[])
{
  std::vector<std::string> ifname_data, ifname_mc;
  std::string ofname, cfname;
  double lumi;
  unsigned nbins;
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
      ("nbins,n", po::value(&nbins)->default_value(10000),
       "number of bins in intermediate histograms")
      ("lumi,l", po::value(&lumi)->required(),
       "luminosity of input data")
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

  std::vector<ifname_t> ifname;
  for (const auto& f : ifname_data) ifname.emplace_back(ifname_t{&f,false});
  for (const auto& f : ifname_mc  ) ifname.emplace_back(ifname_t{&f,true });

  TFile *fout = new TFile(ofname.c_str(),"recreate");
  if (fout->IsZombie()) return 1;
  cout << "Output file: " << fout->GetName() << endl;

  hsb h_Dphi_j_j("Dphi_j_j",nbins,0,M_PI),
      h_Dy_j_j("Dy_j_j",nbins,0,8.8);

  for (const auto& f : ifname) {
    TFile* file = new TFile(f.name->c_str(),"read");
    if (file->IsZombie()) return 1;

    cout << ( f.is_mc ? "MC:" : "Data:" ) << ' '
         << file->GetName() << endl;

    double n_all = 1;
    if (f.is_mc) {
      TIter next(file->GetListOfKeys());
      TKey *key;
      while ((key = static_cast<TKey*>(next()))) {
        std::string name(key->GetName());
        if (name.substr(0,8)!="CutFlow_" ||
            name.substr(name.size()-18)!="_noDalitz_weighted") continue;
        TH1 *h = static_cast<TH1*>(key->ReadObj());
        cout << h->GetName() << endl;
        n_all = h->GetBinContent(3);
        cout << h->GetXaxis()->GetBinLabel(3) << " = " << n_all << endl;
        break;
      }
    }

    TTreeReader reader("CollectionTree", file);
    TTreeReaderValue<Char_t>  isPassed (reader, "HGamEventInfoAuxDyn.isPassed");
    TTreeReaderValue<Float_t> weight   (reader, "HGamEventInfoAuxDyn.weight");
    std::experimental::optional<TTreeReaderValue<Float_t>> cs_br_fe;
    if (f.is_mc) cs_br_fe.emplace(reader,
      "HGamEventInfoAuxDyn.crossSectionBRfilterEff");
    TTreeReaderValue<Float_t> m_yy     (reader, "HGamEventInfoAuxDyn.m_yy");
    TTreeReaderValue<Int_t>   N_j_30   (reader, "HGamEventInfoAuxDyn.N_j_30");
    TTreeReaderValue<Float_t> Dphi_j_j (reader, "HGamEventInfoAuxDyn.Dphi_j_j");
    TTreeReaderValue<Float_t> Dy_j_j   (reader, "HGamEventInfoAuxDyn.Dy_j_j");

    using tc = ivanp::timed_counter<Long64_t>;
    for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
      if (!*isPassed) continue;
      if (*N_j_30 < 2) continue;

      if (!f.is_mc) { // background from data
        if (!in(*m_yy,mass_range)) continue;
        // Assuming weight in data files is always 1
        h_Dphi_j_j.tmp->Fill(*Dphi_j_j);
        h_Dy_j_j  .tmp->Fill(*Dy_j_j  );
      } else { // signal from MC
        if (!in(*m_yy,mass_window)) continue;
        h_Dphi_j_j.tmp->Fill(*Dphi_j_j,*weight***cs_br_fe);
        h_Dy_j_j  .tmp->Fill(*Dy_j_j  ,*weight***cs_br_fe);
      }
    }

    for (auto* h : hsb::all)
      (f.is_mc ? h->sig : h->bkg)->Add(h->tmp,1./n_all);

    for (auto* h : hsb::all) h->tmp->Reset();

    file->Close();
    delete file;
  }
  for (auto* h : hsb::all) h->sig->Scale(lumi);

  // TODO: store info in the root file

  fout->Write();
  fout->Close();
  delete fout;

  return 0;
}
