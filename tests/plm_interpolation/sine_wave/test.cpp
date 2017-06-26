#include <iostream>
#include <vector>
#include <cmath>
#include "utilities.hpp"
#include "ideal_gas.hpp"
#include "spatial_distribution.hpp"
#include "spatial_reconstruction.hpp"
#include "srhydro.hpp"
#include "pcm.hpp"
#include "van_leer.hpp"
#include "diagnostics.hpp"
#include "advanced_hydrodynamic_variables.hpp"

using namespace std;

namespace {
  vector<double> interpolated_values(const HydroSnapshot& hs, 
				     const SpatialReconstruction& sr)
  {
    const vector<pair<Primitive,Primitive> > temp =
      sr.interpolateAll(hs,0);

    vector<double> res;
    for(size_t i=0;i<temp.size();++i){
      res.push_back(temp[i].first.Density);
      res.push_back(temp[i].second.Density);
    }

    return res;
  }

  vector<double> calculated_values(vector<double> const& edges,
				   SpatialDistribution const& velocity)
  {
    vector<double> res;
    for(size_t i=1;i<edges.size()-1;++i){
      res.push_back(velocity(edges[i]));
      res.push_back(velocity(edges[i]));
    }
    return res;
  }

  double l1_error_norm(vector<double> const& v1, 
		       vector<double> const& v2)
  {
    if(v1.size()!=v2.size())
      throw "Vectors have different length";

    double res = 0;
    for(size_t i=0;i<v1.size();++i)
      res += abs(v1[i]-v2[i]);
    return res/static_cast<double>(v1.size());
  }

  double test_interpolation(SpatialReconstruction& sr,
			    size_t np)
  {
    vector<double> edges = linspace(0,1,np);
    SineWave density(0.5,1,0,1);
    Uniform pressure(1);
    Uniform velocity(0);
    IdealGas eos(4./3.);
    vector<Primitive> cells = InitCells(edges,
					density,
					pressure,
					velocity);
    HydroSnapshot hs(edges, cells);
    vector<double> numeric = interpolated_values(hs, sr);
    vector<double> analytic = calculated_values(edges,density);
    return l1_error_norm(numeric, analytic);
  }

  vector<double> test_series(SpatialReconstruction& sr,
			     vector<size_t> const& np_range)
  {
    vector<double> res(np_range.size(),0);
    for(size_t i=0;i<np_range.size();++i)
      res[i] = test_interpolation(sr,np_range[i]);
    return res;
  }

  vector<size_t> arange(size_t nstart, size_t nstop, size_t nskip)
  {
    vector<size_t> res((nstop-nstart)/nskip,0);
    for(size_t i=0;i<res.size();++i)
      res[i] = nstart + i*nskip;
    return res;
  }

  vector<double> order_edges(vector<double> const& edges)
  {
    vector<double> res;
    for(size_t i=1;i<edges.size()-1;++i){
      res.push_back(edges[i]);
      res.push_back(edges[i]);
    }
    return res;
  }

  void interpolation_snapshot(size_t np)
  {
    PCM pcm;
    VanLeer plm;
    vector<double> edges = linspace(0,1,np);
    SineWave density(0.5,1,0,1);
    Uniform pressure(1);
    Uniform velocity(0);
    IdealGas eos(4./3.);
    vector<Primitive> cells = InitCells(edges,
					density,
					pressure,
					velocity);
    const HydroSnapshot hs(edges, cells);
    vector<double> pcm_data = interpolated_values(hs, pcm);
    vector<double> plm_data = interpolated_values(hs, plm);
    vector<double> analytic = calculated_values(edges,velocity);

    vector<double> oedges = order_edges(edges);

    write_to_file(oedges,"edges.txt");
    write_to_file(pcm_data,"pcm_data.txt");
    write_to_file(plm_data,"plm_data.txt");
    write_to_file(analytic,"analytic_data.txt");
  }

  void interpolation_comparison(void)
  {
    PCM pcm;
    VanLeer plm;
    vector<size_t> np_range = arange(100,1001,100);
    vector<double> l1_pcm = test_series(pcm,np_range);
    vector<double> l1_plm = test_series(plm,np_range);
    write_to_file(np_range,"points.txt");
    write_to_file(l1_pcm,"l1_pcm.txt");
    write_to_file(l1_plm,"l1_plm.txt");  
  }
}

int main(void)
{
  interpolation_snapshot(10);

  interpolation_comparison();

  ofstream("test_terminated_normally.res").close();
  return 0;
}
