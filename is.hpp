#ifndef _ITERSCAN_HPP_
#define _ITERSCAN_HPP_
/************************************************************
Iterative Scan
************************************************************/
#include <iostream>
#include <fstream>
#include <vector>
#include <set>

#include <boost/random.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map.hpp>
#include "clusterdens.hpp"

using std::vector;
using std::set;

template <typename Graph, typename ClusterDensity, typename EdgeWeightMap>
void iterative_scan_core(
  const Graph & g, 
  const set<typename graph_traits<Graph>::vertex_descriptor> & seed,
  set<typename graph_traits<Graph>::vertex_descriptor> & final,
  ClusterDensity density,
  EdgeWeightMap edge_weight)
{
  final = seed;
  int e_in, e_out;
  double dens = density(g, edge_weight, final, &e_in, &e_out);
  double new_dens = 0;
  bool increased = true;
  while (increased) {
    increased = false;
    typename graph_traits<Graph>::vertex_iterator vi, vi_end;
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi) {
      new_dens = density(g, edge_weight, final, *vi, dens, &e_in, &e_out);
      if (new_dens - dens > 0) {
        increased = true;
        if (final.contains(*vi)) {
          final.insert(*vi);
        }
        else {
          final.remove(*vi);
        }
      }
    }
  }
}

template <typename Graph, typename ClusterDensity>
void iterative_scan_core(
  const Graph & g, 
  const set<typename graph_traits<Graph>::vertex_descriptor> & seed,
  set<typename graph_traits<Graph>::vertex_descriptor> & final)
{
  unit_weight_map<Graph> u;
  iterative_scan_core<Graph, ClusterDensity, unit_weight_map<Graph> >
    (g, seed, final, u);
}

template <typename Graph, typename ClusterDensity, typename EdgeWeightMap>
void iterative_scan_clustering(
  Graph & g, 
  vector<set<typename graph_traits<Graph>::vertex_descriptor> > & clusters,
  EdgeWeightMap edge_weight,
  int failures = 5)
{
  ecuyer1988 rand_gen;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::edge_descriptor Edge;

  vector<vector<Vertex> > processed;
  set<Vertex>* final;
  set<Vertex> cur;
  int fail = 0;
  while (fail < failures) {
    // Choose a random edge.
    cur.clear();
    Edge e = random_edge<Graph, ecuyer1988>(g, rand_gen);
    cur.insert(source(e, g));
    cur.insert(target(e, g));

    // Improve that edge.
    final = new set<Vertex>;
    iterative_scan_core<Graph, ClusterDensity, EdgeWeightMap>
      (g, cur, *final, edge_weight);

    // If this cluster has not been found, add it to the list.
    // Otherwise, increase the failures.
    bool found = false;
    for (int c = 0; c < clusters.size(); ++c) {
      if (clusters[c] == *final) {
        found = true;
      }
    }
    if (found) {
      ++fail;
    } else {
      fail = 0;
      clusters.push_back(*final);
    }
  }
}

template <typename Graph, typename ClusterDensity>
void iterative_scan_clustering(
  Graph & g, 
  vector<set<typename graph_traits<Graph>::vertex_descriptor> > & clusters,
  int failures = 5)
{
  unit_weight_map<Graph> u;
  iterative_scan_clustering<Graph, ClusterDensity, unit_weight_map<Graph> >
    (g, clusters, u, failures);
}

template <typename Graph, typename ClusterDensity, typename EdgeWeightMap>
void iterative_scan_seed(
  Graph & g, 
  const vector<set<typename graph_traits<Graph>::vertex_descriptor> > & seeds,
  vector<set<typename graph_traits<Graph>::vertex_descriptor> > & clusters,
  EdgeWeightMap edge_weight)
{
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typename vector<set<Vertex> >::iterator ci, ci_end;
  ci_end = seeds.end();
  set<Vertex>* final;
  for (ci = seeds.begin(); ci != ci_end; ++ci) {
    final = new set<Vertex>;
    iterative_scan_core<Graph, ClusterDensity, EdgeWeightMap>
      (g, edge_weight, *ci, *final);
    clusters.push_back(*final);
  }
}

template <typename Graph, typename ClusterDensity>
void iterative_scan_seed(
  Graph & g, 
  const vector<set<typename graph_traits<Graph>::vertex_descriptor> > & seeds,
  vector<set<typename graph_traits<Graph>::vertex_descriptor> > & clusters)
{
  unit_weight_map<Graph> u;
  iterative_scan_seed<Graph, ClusterDensity, unit_weight_map<Graph> >
    (g, seeds, clusters, u);
}

#endif // _ITERSCAN_HPP_
