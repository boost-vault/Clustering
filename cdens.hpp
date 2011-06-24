#ifndef _CLUSTERDENS_HPP_
#define _CLUSTERDENS_HPP_

#include <set>
#include <boost/graph/graph_traits.hpp>

using namespace boost;
using std::set;

template<typename Graph>
struct unit_weight_map {
  float operator[](graph_traits<Graph>::edge_descriptor e) {
    return 1.0f;
  }
};

template<typename Graph>
float get(const unit_weight_map<Graph> & m, 
    const graph_traits<Graph>::edge_descriptor & e) {
  return 1.0f;
}

abstract class density {
  // Initialize function.
  // g <in> is the graph the cluster belongs to.
  virtual void operator()(
    const Graph & g,
    const EdgeWeightMap & edge_weight,
    const set<typename graph_traits<Graph>::vertex_descriptor> & members,
    float & dens,
    float & w_in,
    float & w_out
  ) = 0;

  // Update function.
  virtual void operator()(
    const Graph & g,
    const EdgeWeightMap & edge_weight,
    const set<typename graph_traits<Graph>::vertex_descriptor> & members,
    typename graph_traits<Graph>::vertex_descriptor v,
    float & dens,
    float & w_in,
    float & w_out
  ) = 0;
}

template <typename Graph, typename EdgeWeightMap = unit_weight_map<Graph> >
class basic_density : public cluster_density<Graph, EdgeWeightMap> {
private:
  float m_weight_in;
  float m_weight_out;

public:
  virtual void init(
    const Graph & g, 
    const EdgeWeightMap & edge_weight, 
    const set<typename graph_traits<Graph>::vertex_descriptor> & members)
  {
    typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename graph_traits<Graph>::out_edge_iterator OutEdgeIt;
    m_members = members;
    m_weight_in = 0;
    m_weight_out = 0;
    typename set<Vertex>::iterator vi;
    for (vi = m_members.begin(); vi != m_members.end(); ++vi) {
      OutEdgeIt ei, ei_end;
      for (tie(ei, ei_end) = out_edges(*vi, g); 
        ei != ei_end; ++ei) 
      {
        if (m_members.count(target(*ei, g)) > 0) {
          m_weight_in += get(edge_weight, *ei);
        } else {
          m_weight_out += get(edge_weight, *ei);
        }
      }
    }
    // Every interior edge has been counted twice.
    // TODO: Should not do this if the graph is directed.
    m_weight_in /= 2;

    m_density = calc_density(m_weight_in, m_weight_out, 
      m_members.size(), num_vertices(g) - m_members.size());
  }

  virtual bool update(
    const Graph & g, 
    const EdgeWeightMap & edge_weight, 
    typename graph_traits<Graph>::vertex_descriptor v)
  {
    typedef typename graph_traits<Graph>::out_edge_iterator OutEdgeIt;
    float w_in = 0;
    float w_out = 0;
    OutEdgeIt ei, ei_end;
    for (tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
      if (m_members.count(target(*ei, g)) > 0) {
        w_in += get(edge_weight, *ei);
      } else {
        w_out += get(edge_weight, *ei);
      }
    }
    if (m_members.count(v) > 0) {
      float new_dens = calc_density(
        m_weight_in - w_in, 
        m_weight_out - w_out + w_in, 
        m_members.size() - 1, num_vertices(g) - (m_members.size() - 1));
      if (new_dens > m_density) {
        //cout << m_density << " < " << new_dens 
        //     << ", removing vertex " << v << endl;
        //cout << "m_weight_in=" << m_weight_in << ","
        //     << "m_weight_out=" << m_weight_out << ","
        //     << "w_in=" << w_in << ","
        //     << "w_out=" << w_in << ","
        //     << "m_members.size()=" << m_members.size() << endl;
        m_members.erase(v);
        m_weight_in -= w_in;
        m_weight_out -= (w_out - w_in);
        m_density = new_dens;
        return true;
      }
    } else {
      float new_dens = calc_density(
        m_weight_in + w_in, 
        m_weight_out - w_in + w_out, 
        m_members.size() + 1, num_vertices(g) - (m_members.size() + 1));
      if (new_dens > m_density) {
        //cout << m_density << " < " << new_dens 
        //     << ", adding vertex " << v << endl;
        //cout << "m_weight_in=" << m_weight_in << ","
        //     << "m_weight_out=" << m_weight_out << ","
        //     << "w_in=" << w_in << ","
        //     << "w_out=" << w_in << ","
        //     << "m_members.size()=" << m_members.size() << endl;
        m_members.insert(v);
        m_weight_in += w_in;
        m_weight_out += (w_out - w_in);
        m_density = new_dens;
        return true;
      }
    }
    return false;
  }

  virtual float calc_density(
    float w_in,
    float w_out,
    int vert_in,
    int vert_out) = 0;
};

template <typename Graph, typename EdgeWeightMap = unit_weight_map<Graph> >
class average_degree : public basic_density<Graph, EdgeWeightMap> {
  virtual float calc_density(
    float w_in,
    float w_out,
    int vert_in,
    int vert_out)
  {
    // cout << "density=" << (w_in / vert_in) << endl;
    return w_in / vert_in;
  }
};

template <typename Graph, typename EdgeWeightMap = unit_weight_map<Graph> >
class average_weight : public basic_density<Graph, EdgeWeightMap> {
  virtual float calc_density(
    float w_in,
    float w_out,
    int vert_in,
    int vert_out)
  {
    return w_in / vert_in;
  }
};

template <typename Graph, typename EdgeWeightMap = unit_weight_map<Graph> >
class weight_ratio : public basic_density<Graph, EdgeWeightMap> {
  virtual float calc_density(
    float w_in,
    float w_out,
    int vert_in,
    int vert_out)
  {
    return w_in / (w_in + w_out);
  }
};

template <typename Graph, typename EdgeWeightMap = unit_weight_map<Graph> >
class weight_prob : public basic_density<Graph, EdgeWeightMap> {
  virtual float calc_density(
    float w_in,
    float w_out,
    int vert_in,
    int vert_out)
  {
    return w_in / (float)(vert_in*(vert_in - 1)/2);
  }
};

#endif // _CLUSTERDENS_HPP_
