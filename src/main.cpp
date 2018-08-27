// Main file, contains the overall flow (in main() ) and some important functions

#include <iostream>
#include <set>
#include <algorithm>
#include <map>
#include <vector>
#include <cmath>
#include <stack>
#include <queue>

//Not doing much in practice
int debug = 0;
//Not the standard assert!
#define assert(a, m) {if (!(a)) {cout << m << endl; exit(0);}}


using namespace std;

inline double dist(double x, double y) { return sqrt(x*x+y*y); }
inline double dist2(double x, double y) { return (x*x+y*y); }

//Basics for 3d coordinate representation
struct point {
  double x, y, z;
  point() {}
  point(double x, double y, double z) : x(x),y(y),z(z) {}
  inline point operator-(const point&p) {
    return point(x-p.x, y-p.y, z-p.z);
  }
  inline point operator+(const point&p) {
    return point(x+p.x, y+p.y, z+p.z);
  }
  inline double operator*(const point&p) {
    return x*p.x+y*p.y+z*p.z;
  }
  inline point operator*(double f) {
    return point(x*f, y*f, z*f);
  }
};
ostream& operator<<(ostream&out, const point p) {
  out << p.x << ' ' << p.y << ' ' << p.z;
  return out;
}
inline double dist(const point&p) { return sqrt(p.x*p.x+p.y*p.y+p.z*p.z); }

//Structure for storing promising triples of hits
struct triple {
  int x, y, z; //hit ids
  triple() {}
  triple(int a, int b, int c) : x(a), y(b), z(c) {}
};
bool operator<(const triple&a, const triple&b) {
  if (a.z != b.z) return a.z < b.z;
  if (a.y != b.y) return a.y < b.y;
  return a.x < b.x;
}
bool operator==(const triple&a, const triple&b) {
  return a.x==b.x && a.y==b.y && a.z==b.z;
}

//Convert "dir" to polar (really cylindrical coordinates) coordinates, suffix p  (as in "refp") usually means polar coodinates throughout the code
point topolar(const point&dir, const point&ref, const point&refp) {
  return point(ref.x*dir.x+ref.y*dir.y, ref.x*dir.y-ref.y*dir.x, dir.z*refp.x);
}

#include "input.hpp"

//does hits a and b correspond to the same particle?
int samepart(int a, int b) {
  long long aa = truth_part[a];
  long long bb = truth_part[b];
  return aa == bb && aa;
}

#include "load.hpp"

//reused global array for storing matching hits from PolarModule->getNear function
int match[200000];
//Assignment of hits to paths
int assignment[200000];

//Number of directly adjacent hits in layers
int next_layer[48][48];
//Threshold in next_layer to consider layers adjacent
int adj_thres = 1000; //TODO: tweak


#include "polarmodule.hpp"
#include "old_pairs.hpp"

#include "io.hpp"
#include "score.hpp"
#include "density.hpp"
#include "investigate.hpp"

#include "new_pairs.hpp"

//Initialize next_layer
void loadAdjacent() {
  FILE*fp = fopen("trained/adjacency", "r");
  if (!fp) {
    cout << "Could not open adjacency" << endl;
    exit(0);
  }
  forij(48,48)
    int tmp = fscanf(fp, "%d", &next_layer[i][j]);
  fclose(fp);
}

//The next_layer adjacency matrix above is simply the output of this function summed over the 100 training events
void truthAdjacent() {
  //count_layer[48] = truth_tracks.size();
  for (auto&p : truth_tracks) {
    int prev = -1;//48;
    for (int i : p.second) {
      int m = metai[i];
      //count_layer[m]++;
      if (prev != -1) {
	next_layer[prev][m]++;
      }
      prev = m;
    }
  }
}


inline bool z_cmp(const int a, const int&b) {
  return hits[a].z < hits[b].z;
}

inline bool r_cmp(const int&a, const int&b) {
  return hits[a].x*hits[a].x+hits[a].y*hits[a].y < hits[b].x*hits[b].x+hits[b].y*hits[b].y;
}


vector<int>*readTubes() {
  vector<int>*tube = new vector<int>[48](); // List of hits in each layer
  map<long long,int> added[48];

  for (int i = 1; i < hits.size(); i++) {
    if (!assignment[i])
      tube[metai[i]].push_back(i);
  }
  for (int i = 0; i < 48; i++) {
    if (layer[i].type == Tube)
      sort(tube[i].begin(), tube[i].end(), z_cmp);
    else
      sort(tube[i].begin(), tube[i].end(), r_cmp);
  }
  return tube;
}




//Extend triples based on hits "ai" and "bi" to layer "li", add then to "triples", do this using an approximate straight line going through a and b (though taking into account some slightly curved helices by looking into an elliptic region given by "dp", "xp" and "bap". Again postfix "p" for polar coordinates)
double extendTripleLine(vector<triple>&triples, int ai, int bi, point&a, point&b, int li, PolarModule&mod, int rev = 0) {
  point d, dp, xp, bap;
  if (prepareTripleScore(ai, bi, li, d, dp, xp, bap)) return 0;

  //xp = normalize(xp)*0.999;
  //xp = point(0,0,0.5);
  const double target0 = 0.1, target = 10;//0.5;

  //double mid0 = -findDensity(dp, xp, target0, li);
  double mid = findDensity(dp, xp, target, li);

  int matches = mod.getNear(dp, xp, bap, mid, match);

  int mini;
  double best = target;
  vector<pair<double, int> > v;
  for (int i = 0; i < matches; i++) {
    int ci = match[i];
    //if (ci == ai || ci == bi) continue;
    double s = scoreTriple(ai, bi, ci);//evaluateScore(ci, dp, xp, bap);//
    v.push_back(make_pair(s, ci));
  }
  //Take only best ones
  sort(v.begin(), v.end());
  for (int i = 0; i < v.size(); i++) {
    if (i >= target) break;
    triple t(ai, bi, v[i].second);
    if (rev) swap(t.x,t.z);
    if (acceptTriple(t)) //Prune most triples
      triples.push_back(t);
  }
  return 1;
  //cout << added << endl;
}

//Extend triples using origin, similar to extendTripleLine. Except that it uses the assumption that the particle start in the origin instead of the straight line assumption. This is not used, as the triples starting from the origin are easy anyway, and therefore captured by extendTripleLine.
double extendTripleOrigin(vector<triple>&triples, int ai, int bi, point&a, point&b, int li, PolarModule&mod, double target = 0.5, int rev = 0) {
  hits[0] = polar[0] = point(0,0,0);

  point d, dp, xp, bap;
  if (prepareQuadrupleScore(0, ai, bi, li, d, dp, xp, bap)) return 0;

  //const double target0 = 0.1, target = 0.5;

  //double mid0 = findDensity(dp, xp, target0, li);
  double mid = findDensity(dp, xp, target, li);

  int matches = mod.getNear(dp, xp, bap, mid, match);

  int mini;
  double best = target;
  vector<pair<double, int> > v;
  for (int i = 0; i < matches; i++) {
    int ci = match[i];
    //double s = scoreTriple(ai, bi, ci);
    double s = evaluateScore(ci, dp, xp, bap);
    v.push_back(make_pair(s, ci));
  }
  if (!matches) return 0;
  sort(v.begin(), v.end());
  double thres = v[0].first*4;
  for (int i = 0; i < v.size(); i++) {
    if (i >= 2 || v[i].first > thres) break;// && v[i].first > mid0) break; //i >= 1 and added < 1 gives better than old
    if (rev)
      triples.push_back(triple(v[i].second, bi, ai));
    else
      triples.push_back(triple(ai, bi, v[i].second));
  }
  return 1;
  //cout << added << endl;
}



//Expand all pairs into triples (possibly many triples per pair). method 1 uses origin assumption, method 0 uses straight line assumption.
vector<triple> findTriples(vector<pair<int, int> >&pairs, PolarModule*mod, int method = 0, double target = 0.5) {
  vector<triple> triples;
  vector<double> v[48];
  for (auto p : pairs) {
    load(pairs.size(), "Find triples using origin");
    int ai = p.first, bi = p.second;
    //if (!samepart(ai, bi)) continue;
    point&a = hits[ai], &b = hits[bi];

    int added = 0;
    for (int li = metai[bi]+1; added < 100*method+1 && li < 48; li++) {
      if (next_layer[metai[bi]][li] < adj_thres) continue;
      if (method == 1 && extendTripleLine(triples, ai, bi, a, b, li, mod[li])) added++;
      if (method == 0 && extendTripleOrigin(triples, ai, bi, a, b, li, mod[li], target)) added++;
    }
    //extendTriple(triples, ai, bi, a, b, metai[ai], mod[metai[ai]]);
    //extendTriple(triples, ai, bi, a, b, metai[bi], mod[metai[bi]]);
    /*for (int li = metai[ai]-1; li >= 0; li--) {
      if (next_layer[li][metai[ai]] < adj_thres) continue;
      if (method == 1 && extendTripleLine(triples, bi, ai, a, b, li, mod[li], 1)) break;
      //if (extendTripleOrigin(triples, bi, ai, a, b, li, mod[li], 1)) break;
      }*/
  }
  /*
  for (int i = 0; i < 48; i++) {
    sort(v[i].begin(), v[i].end());
    double t = v[i].size() ? v[i][v[i].size()*99/100]+0.1 : 0.;
    cout << t << ", ";
  }
  cout << endl;
  */

  return triples;
}


//Add triples extended from a duplicate pair (hits on the same layer)
void addDuplicateTriples(vector<triple>&triples, PolarModule*mod, int ai, int bi) {
  vector<int> s[48];
  for (int li = 0; li < 48; li++) {
    if (next_layer[li][metai[ai]] < adj_thres/2 &&
	next_layer[metai[ai]][li] < adj_thres/2) continue;
    if (li == metai[ai]) continue;

    point d, dp, xp, bap;
    if (prepareTripleScore(ai, bi, li, d, dp, xp, bap) &&
	prepareTripleScore(bi, ai, li, d, dp, xp, bap)) continue;

    xp = normalize(xp)*0.99;
    //xp = point(0,0,0);
    const double target = 100;
    double mid = findDensity(dp, xp, target, li);
    int matches = mod[li].getNear(dp, xp, bap, mid, match);

    double best = target;
    vector<pair<double, int> > v;
    for (int i = 0; i < matches; i++) {
      int ci = match[i];
      double s = evaluateScore(ci, dp, xp, bap);
      v.push_back(make_pair(s, ci));
    }
    sort(v.begin(), v.end());
    for (int i = 0; i < v.size(); i++) {
      if (i >= target) break;
      s[li].push_back(v[i].second);
      //triples.push_back(triple(ai, bi, v[i].second));
    }
    //break;
  }

  for (int la = 0; la < 48; la++) {
    for (int ci : s[la]) {
      int xi = ai, yi = bi, zi = ci;
      if ((hits[zi].z-hits[yi].z)*(hits[yi].z-hits[xi].z) < 0) swap(xi, yi);
      if (metai[zi] < metai[xi]) swap(xi, zi);
      for (int li = 0; li < 48; li++) {
	if (li >= metai[xi] && li <= metai[zi]) continue;
	if (next_layer[li][metai[zi]] < adj_thres/2 &&
	    next_layer[metai[xi]][li] < adj_thres/2) continue;

	int xi_ = xi, zi_ = zi;
	if (li < metai[xi_]) swap(xi_, zi_);

	point d, dp, xp, bap;
	if (prepareQuadrupleScore(xi_, yi, zi_, li, d, dp, xp, bap)) continue;

	double target2 = 0.5;
	double tt = findDensity(dp, xp, target2, li);
	int matches = mod[li].getNear(dp, xp, bap, tt, match);
	vector<pair<double, int> > v;
	for (int i = 0; i < matches; i++) {
	  int di = match[i];
	  double s = evaluateScore(di, dp, xp, bap);
	  v.push_back(make_pair(s, di));
	}
	sort(v.begin(), v.end());
	for (int i = 0; i < v.size() && i < target2; i++) {
	  triple t(ai, ci, v[i].second);
	  if (metai[t.x] > metai[t.y]) swap(t.x,t.y);
	  if (metai[t.y] > metai[t.z]) swap(t.y,t.z);
	  if (metai[t.x] > metai[t.y]) swap(t.x,t.y);
	  assert(metai[t.x] < metai[t.y] && metai[t.y] < metai[t.z], "Wrong order!");
	  triples.push_back(t);
	}
      }
    }
  }

  /*int good = 0;
  for (int i = 0; i < 48; i++) {
    int found = 0;
    for (int j : s[i]) if (samepart(ai, j)) found++;
    if (found) good++;
  }
  if (good >= 2) triples.push_back(triple(ai, ai, ai));*/
}

//(Try to) find all duplicates in the dataset by looking for point close together given cells' data velocity direction
vector<pair<int, int> > findDuplicates(vector<int>*tube, PolarModule*mod) {
  vector<pair<int, int> > pairs;

  for (int li = 0; li < 48; li++) {
  point d, dp, xp, bap;
  double target = 0.05;

  map<pair<int, int>, double> v;
  int fails = 0;
  set<long long> found;
  for (int hit_id : tube[li]) {
    for (int k = 0; k < 2; k++) {
      d = hits[hit_id];
      dp = polar[hit_id];
      point dir = hit_dir[hit_id][k];
      xp = point(0,0,0);
      point bap = topolar(dir, d, dp);
      bap = bap*(1./(layer[li].type == Disc ? bap.z : bap.x));

      double tt = findDensity(dp, xp, target, li);
      int matches = mod[li].getNear(dp, xp, bap, tt, match);
      for (int i = 0; i < matches; i++) {
	int di = match[i];
	if (metaz[di] == metaz[hit_id]) continue;
	auto p = make_pair(hit_id, di);
	double s = getDensity3(dp, xp, evaluateScore(di, dp, xp, bap), li);
	if (!v.count(p) || s < v[p])
	  v[p] = s;
      }
    }
  }
  for (auto&i : v) {
    int a = i.first.first, b = i.first.second;
    if (a < b || !v.count(make_pair(b,a)))
      pairs.push_back(make_pair(a, b));
  }
  }
  return pairs;
}


//Try to find triples using by first finding duplicates, and then expanding them to triples in both direction
vector<triple> findTriplesDuplicates(vector<int>*tube, PolarModule*mod) {
  vector<pair<int, int> > pairs = findDuplicates(tube, mod);
  cout << pairs.size() << " pairs" << endl;
  vector<triple> triples;
  for (auto&p : pairs) {
    load(pairs.size(), "Finding triples using duplicates");
    addDuplicateTriples(triples, mod, p.first, p.second);
  }
  return triples;
}

//Cheating comparison to findTriplesDuplicates
vector<triple> findTriplesDuplicatesTruth(PolarModule*mod) {
  vector<triple> triples;
  for (auto&t : truth_tracks) {
    load(truth_tracks.size());
    if (dist(start_pos[t.first].x, start_pos[t.first].y) < 10) continue;
    map<int, vector<int> > c;
    for (int i : t.second)
      c[metai[i]].push_back(i);
    for (auto&p : c) {
      if (p.second.size() > 1) {
	for (int i = 0; i < p.second.size(); i++)
	  for (int j = 0; j < i; j++)
	    addDuplicateTriples(triples, mod, p.second[j], p.second[i]);
      }
    }
  }
  return triples;
}



//Extend the helix going through hits with ids "ai", "bi", "ci", to layer "li". Do this by looking at the intersection with the layer, and expecting around "target" continuation for each outlier triple. "li" must be after metai[ci]
int extend3(PolarModule*mod, int ai, int bi, int ci, int li, double target = 0.5) {
  point d, dp, xp, bap;
  if (prepareQuadrupleScore(ai,bi,ci,li,d,dp,xp,bap)) return -2;

  //double mins = 0.5;//d*d*1e-4;
  double tt = 400;//findDensity(dp, xp, target, li);
  double mins = 1e9;//target;//1e9;

  double fac = getDensity3(dp, xp, tt, li)/tt;
  tt = target/fac;
  int matches = mod[li].getNear(dp, xp, bap, tt, match);
  int mini = -1;
  for (int i = 0; i < matches; i++) {
    int ti = match[i];
    if (assignment[ti]) continue;
    double s = evaluateScore(ti,dp,xp,bap)*fac;
    if (s < mins) {
      mins = s;
      mini = ti;
    }
  }
  return mini;
}

//Similar to extend3, but now with "li" in between metai[bi] and metai[ci]
int extend4(PolarModule*mod, int ai, int bi, int ci, int li, double target = 0.5) {
  point d, dp, xp, bap;
  if (prepareQuadrupleScore(ai,bi,ci,li,d,dp,xp,bap,-1)) return -2;

  //double mins = 0.5;//d*d*1e-4;
  double tt = 400;//findDensity(dp, xp, target, li);
  double mins = 1e9;//target;//1e9;

  double fac = getDensity3(dp, xp, tt, li)/tt;
  tt = target/fac;
  int matches = mod[li].getNear(dp, xp, bap, tt, match);
  int mini = -1;
  for (int i = 0; i < matches; i++) {
    int ti = match[i];
    if (assignment[ti]) continue;
    double s = evaluateScore(ti,dp,xp,bap)*fac;
    if (s < mins) {
      mins = s;
      mini = ti;
    }
  }
  return mini;
}


//Expand triples into paths, do this by expanding helices through the outermost 3 point in each direction repeatedly
vector<vector<int> > findPaths(vector<triple>&triples, PolarModule*mod) {
  vector<vector<int> > paths;
  paths.reserve(triples.size());
  for (auto p : triples) {
    load(triples.size(), "Find paths", 0);
    vector<int> v;
    int ai = p.x, bi = p.y, ci = p.z, misses = 0;
    /*for (int li = metai[bi]-1; li > metai[ai]; li--) {
      if (next_layer[li][metai[bi]] < adj_thres) continue;
      int di = extend4(mod, ci, bi, ai, li, 0.1);
      if (di > 0) {
	ai = p.x = di;
	break;
      } else if (di == -1) break;
    }
    for (int li = metai[bi]+1; li < metai[ci]; li++) {
      if (next_layer[metai[bi]][li] < adj_thres) continue;
      int di = extend4(mod, ai, bi, ci, li, 0.1);
      if (di > 0) {
	ci = p.z = di;
	break;
      } else if (di == -1) break;
      }*/

    for (int li = metai[ai]-1; li >= 0; li--) {
      if (next_layer[li][metai[ai]] < adj_thres) continue;
      int di = extend3(mod, ci, bi, ai, li);
      if (di > 0) {
	v.push_back(di);
	ci = bi;
	bi = ai;
	ai = di;
	misses = 0;
      } else if (di == -1) misses++;
      if (misses == 1) break;
    }
    reverse(v.begin(), v.end());

    ai = p.x, bi = p.y, ci = p.z;

    v.push_back(ai), v.push_back(bi), v.push_back(ci);

    misses = 0;
    for (int li = metai[ci]+1; li < 48; li++) {
      if (next_layer[metai[ci]][li] < adj_thres) continue;
      int di = extend3(mod, ai, bi, ci, li);
      if (di > 0) {
	v.push_back(di);
	ai = bi;
	bi = ci;
	ci = di;
	misses = 0;
      } else if (di == -1) misses++;
      if (misses == 2) break;
    }

    v.shrink_to_fit();
    paths.push_back(v);
  }
  return paths;
}






//Attempt to add all duplicate hits (hits on same layer as an already added hit in the path) to the paths
vector<vector<int> > addDuplicates(vector<vector<int> >&paths, PolarModule*mod) {
  vector<vector<int> > extended;

  for (auto&path : paths) {
    load(paths.size(), "Add duplicates");
    if (path.size() < 3) continue;

    vector<int> ext;
    for (int i = 0; i < path.size(); i++) {
      ext.push_back(path[i]);
      int ai, bi, ci;
      if (i < path.size()-2) ai = path[i], bi = path[i+1], ci = path[i+2];
      else if (i == path.size()-1) ai = path[i-2], bi = path[i-1], ci = path[i];
      else ai = path[i-1], bi = path[i], ci = path[i+1];

      int li = metai[path[i]];
      point d, dp, xp, bap;

      //Add on average about "target" outliers to each path
      const double target = 0.1;

      if (prepareDuplicateScore(ai, bi, ci, li, d, dp, xp, bap)) continue;

      double tt = findDensity(dp, xp, target, li);

      int pi = path[i];
      int matches = mod[li].getNear(dp, xp, bap, tt, match);//target/fac

      map<int, pair<double, int> > mins;

      for (int i = 0; i < matches; i++) {
	int di = match[i];
	double s = evaluateScore(di, dp, xp, bap);//*fac;
	if (meta[di].z != meta[pi].z) {
	  int zi = metaz[di];
	  if (!mins.count(zi) || s < mins[zi].first) {
	    mins[zi] = make_pair(s, di);
	  }
	}
      }
      //vector<pair<double, int> > v;
      for (auto&p : mins)
	ext.push_back(p.second.second);
      /*sort(v.begin(), v.end());
      for (auto&p : v)
      ext.push_back(p.second);*/
    }
    ext.shrink_to_fit();
    path.clear();
    path.shrink_to_fit();
    extended.push_back(ext);
  }
  return extended;
}




//Prune away paths based on path score, also take away the start of the path if it doesn't fit well. Also remove duplicate paths using hashing
vector<vector<int> > prunePaths(vector<vector<int> >&paths) {
  vector<vector<int> > pruned;

  vector<vector<int> > hash_list;
  hash_list.resize(1<<17);
  int c = 0;
  double thres = 2.5;
  for (auto&path : paths) {
    load(paths.size(), "Pruning paths", 0);
    while (path.size() >= 4 && scoreQuadrupleDensity(path[3], path[2], path[1], path[0]) > thres) {
      path.erase(path.begin());
    }
    /*while (path.size() >= 4 && scoreQuadrupleDensity(path[path.size()-4], path[path.size()-3], path[path.size()-2], path[path.size()-1]) > 1e2) {
      path.pop_back();
    }*/
    double s = -scorepathDensity(path);
    if (s < 10 && path.size() >= 3) {
      int h = 0;
      for (int i : path) h ^= i;
      h &= (1<<17)-1;
      int found = 0;
      for (int j : hash_list[h]) {
	if (pruned[j] == path) {
	  found = 1;
	  break;
	}
      }
      if (!found) {
	hash_list[h].push_back(pruned.size());
	path.shrink_to_fit();
	pruned.push_back(path);
      }
    }
    path.clear();
    path.shrink_to_fit();
  }
  pruned.shrink_to_fit();
  return pruned;
}








//Data structure to extract index of path with lowers score and dynamically updating, supports:
// - Adding new (index, score) pair
// - Updating (index, score) pair
// - Extracting pair with lowest score
//Note that std::priority_queue is vastly faster than std::set (which was previously used, and makes for easier implementation)
struct myMap2 {
  priority_queue<pair<double, int> > pq;
  double*b;
  int realsize;
  myMap2(int n) {
    b = new double[n];
    realsize = 0;
  }
  ~myMap2() {
    delete[]b;
  }
  void add(int i, double score) {
    pq.push(make_pair(score, i));
    b[i] = score;
    realsize++;
  }
  void update(int i, double score) {
    add(i, score);
    realsize--;
    //Change *4 to a lower constant > 1 for (slightly) lower memory usage, this is currently the bottleneck for memory I think
    if (pq.size() >= realsize*4 && pq.size() > 1e6) {
      priority_queue<pair<double, int> > clean;
      pair<double, int> last(1e9, 1e9);
      while (pq.size()) {
	auto p = pq.top();
	pq.pop();
	if (b[p.second] == p.first && p != last) clean.push(p);
	last = p;
      }
      swap(pq, clean);
      //cout << pq.size() << ' ' << realsize << endl;
    }
    //static int cc = 0;
    //if (++cc%100000 == 0)
      //  cout << pq.size() * 1. / realsize << endl;
  }
  pair<int, double> pop() {
    while (pq.size() && b[pq.top().second] != pq.top().first) pq.pop();
    int i = pq.top().second;
    double score = pq.top().first;
    b[i] = -1e9;
    pq.pop();
    realsize--;
    while (pq.size() && b[pq.top().second] != pq.top().first) pq.pop();
    return make_pair(i, score);
  }
  int notempty() {
    return !pq.empty();
  }
};


//Score a path based on the highest score we could possibly get from it, not used except debugging
double scorepathTruth(vector<int>&path) {
  map<long long, double> c;
  double ma = 0;
  for (int i : path)
    if (i > 0) {
      c[truth_part[i]] += truth_weight[i];
      ma = max(ma, c[truth_part[i]]);
    }
  return ma;
}

//Insert your favourite path scoring function here :)
#define scorepath scorepathDensity

//Find assignment of hits to paths
//Do this by iteratively:
// 1. Take path with highest score
// 2. Assign all hits in that path to the path's index
// 3. Remove all hits in the path from all other paths
// 4. Repeat from step 1 until all paths are empty
vector<vector<int> > findAssignment(vector<vector<int> >&paths, PolarModule*mod, int use_trash = 1) {
  paths.insert(paths.begin(), vector<int>());
  map<int, int> map_assignment;
  myMap2 path_score(paths.size()+hits.size());
  vector<vector<pair<int, int> > > used_by;
  used_by.resize(hits.size());

  for (int i = 1; i < paths.size(); i++) {
    load(paths.size()-1, "Finding assignment (1/2)", 0);
    double score = scorepath(paths[i]);
    path_score.add(i, score);
    for (int j = 0; j < paths[i].size(); j++)
      used_by[paths[i][j]].push_back(make_pair(i, j));
  }

  int total = hits.size()-1;
  for (int i = 1; i < hits.size(); i++) {
    used_by[i].shrink_to_fit();
    if (!assignment[i] && used_by[i].empty()) total--;
  }
  /*
  int added = 0;
  for (int i = 1; i < hits.size(); i++)
    if (!used_by.count(i)) {
      lost += truth_weight[i];
      added++;
      vector<int> v = {i};
      //used_by[i].push_back(make_pair(paths.size(), 0));
      paths.push_back(v);
      double score = scorepath(paths[paths.size()-1]);
      path_score.add(paths.size()-1, score);
    }
    cout << "Added: " << added << endl;*/
  double lost = 0;
  cout << "Unused: " << lost << endl;
  lost = 0;

  /*for (int i = 1; i < hits.size(); i++)
    if (assignment[i])
    assignment[i] += 1e8;*/

  vector<vector<int> > solution_paths;
  solution_paths.push_back(vector<int>());

  int last = 0, merged = 0, thrown = 0;
  while (path_score.notempty()) {
    pair<int, double> pop = path_score.pop();
    int i = pop.first;

    if (paths[i].empty()) continue;

    int r[2] = {};
    for (int p : paths[i]) r[p<0]++;

    int trash = 0;
    if (use_trash) {
      int c = 0;
      for (int hit_id : paths[i]) if (hit_id > 0) c++;
      if (c <= 2) trash = 1;// || -pop.second >= 1e-2) trash = 1;
    }
    for (int hit_id : paths[i]) {
      if (hit_id < 0 || map_assignment.count(hit_id)) {
	if (hit_id > 0)
	  cout << "ERROR: " << hit_id << endl;
	continue;
      }
      load(total, "Finding assignment (2/2)", 0);
      map_assignment[hit_id] = trash ? 0 : solution_paths.size();
      if (trash) lost += truth_weight[hit_id];
    }
    if (!trash)
      solution_paths.push_back(paths[i]);

    for (int hit_id : paths[i]) {
      if (hit_id < 0) continue;

      vector<pair<int, int> >&used = used_by[hit_id];
      for (int k = 0; k < used.size(); k++) {
	if (used[k].first == i) continue;
	vector<int>&pi = paths[used[k].first];
	if (pi.empty()) continue;
	//cout << assignment[pi[used[k].second]] << ' '<< c << endl;
	pi[used[k].second] *= -1;
	double score = scorepath(pi);
	path_score.update(used[k].first, score);
      }
    }

    paths[i].clear();
    paths[i].shrink_to_fit();
  }
  cout << "Unused after: " << lost << endl;

  return solution_paths;
}


//Investigate some stats on the assignment, to try to see of there are any more points to gather
void investigateAssignment(vector<vector<int> >&solution_paths) {
  cout << solution_paths.size() << endl;
  int good = 0;
  map<long long, int> goodc;
  double stolen = 0, missing = 0, between = 0, duplicate = 0, wrong = 0;
  for (int i = 1; i < solution_paths.size(); i++) {
    map<long long, int> c;
    long long bestp = -1;
    int mi = solution_paths[i].size();
    int missed = 0;
    for (int j : solution_paths[i]) {
      if (j <= 0) {missed++;continue;}
      c[truth_part[j]]++;
      if (c[truth_part[j]]*2 > mi) bestp = truth_part[j];
    }
    if (bestp != -1) {
      goodc[bestp]++;
      if (goodc[bestp] == 1) {
	set<int> found, m;
	for (int j : truth_tracks[bestp]) found.insert(j);
	int mi = 1e9, ma = -1;
	for (int j : solution_paths[i]) {
	  int k = abs(j);
	  if (truth_part[k] == bestp && j <= 0) stolen += truth_weight[k];
	  if (truth_part[k] == bestp) {
	    found.erase(k);
	    mi = min(mi, metai[k]);
	    ma = max(ma, metai[k]);
	    m.insert(metai[k]);
	  }
	}
	for (int j : found) {
	  missing += truth_weight[j];
	  if (metai[j] > mi && metai[j] < ma) between += truth_weight[j];
	  if (m.count(metai[j])) duplicate += truth_weight[j];
	}
	good++;
      }
      //log(-scorepathDensity(solution_paths[i])) << endl;
    } else {
      int bad = 0;
      for (int k = 1; k < solution_paths[i].size()-1; k++) {
	auto&v = solution_paths[i];
	if (v[k-1] > 0 && v[k] < 0 && v[k+1] > 0) bad = 1;
      }
      //cerr << bad*1./solution_paths[i].size() << endl;
      for (int j : solution_paths[i])
	if (j > 0)
	  wrong += truth_weight[j];
    }
  }
  cout << good << endl;
  cout << "Stolen: " << stolen << endl;
  cout << "Missing: " << missing << endl;
  cout << "Between: " << between << endl;
  cout << "Duplicate: " << duplicate << endl;
  cout << "Wrong: " << wrong << endl;
  double outside = 0;
  for (auto&p : part_weight)
    if (!goodc.count(p.first)) outside += p.second;
  cout << "Outside: " << outside << endl;


  map<int, int> map_assignment;
  for (int i = 1; i < solution_paths.size(); i++)
    for (int j : solution_paths[i])
      if (j > 0)
	map_assignment[j] = i;

  double completely = 0;
  for (auto&t : truth_tracks) {
    map<int, int> c;
    for (int i : t.second)
      if (map_assignment[i])
	c[map_assignment[i]]++;
    int ma = 0;
    for (auto&p : c) ma = max(ma, p.second);
    if (ma <= 2) completely += part_weight[t.first];
  }
  cout << "Completely lost: " << completely << endl;
}

vector<pair<int, int> > findDuplicates(vector<int>*tube, PolarModule*mod);


//Try to gather some of the score that the main algorithm lost, gets about 0.001 bonus score if I remember correctly. So not very important
vector<vector<int> > extendPaths(vector<vector<int> >&paths) {

  /*{
    double w = 0;
    for (auto&p : truth_tracks) {
      map<int, int> a, b;
      for (int i : p.second) a[metai[i]]++, b[metai[i]] += !assignment[i];
      if (a.size() <= 2) {
	int left = 0;
	for (pair<int, int> i : b) if (a[i.first] > 1 && i.second) left++;
	if (left) {
	  cout << part_weight[p.first] << endl;
	  for (pair<int, int> i : a) cout << b[i.first] << ' ' << i.second <<endl;
	  w += part_weight[p.first];
	}
      }
    }
    cout << w << endl;
    exit(0);
    }*/

  initDensity3();

  PolarModule mod[48];
  for (int i = 0; i < 48; i++)
    mod[i] = PolarModule(i);

  vector<int>*tube = readTubes();


  for (int k = 1; k < paths.size(); k++) {
    vector<int> path;
    int last = -1;
    for (int i : paths[k])
      if (i > 0 && metai[i] != last) {
	last = metai[i];
	path.push_back(i);
      }
    if (path.size() < 3) continue;

    int misses = 0;
    int s = path.size()-3;
    int ai = path[s], bi = path[s+1], ci = path[s+2];
    for (int li = metai[ci]+1; li < 48; li++) {
      if (next_layer[metai[ci]][li] < adj_thres) continue;
      int di = extend3(mod, ai, bi, ci, li, 0.07+!(paths.size()%2)*0.2);
      if (di > 0) {
	path.push_back(di);
	assignment[di] = k;
	ai = bi;
	bi = ci;
	ci = di;
	misses = 0;
      } else if (di == -1) misses++;
      //if (misses == 2) break;
    }

    misses = 0;
    ai = path[0], bi = path[1], ci = path[2];
    for (int li = metai[ai]-1; li >= 0; li--) {
      if (next_layer[li][metai[ai]] < adj_thres) continue;
      int di = extend3(mod, ci, bi, ai, li, 0.07+!(paths.size()%2)*0.2);
      if (di > 0) {
	path.insert(path.begin(), di);
	assignment[di] = k;
	ci = bi;
	bi = ai;
	ai = di;
	misses = 0;
      } else if (di == -1) misses++;
      //if (misses == 2) break;
    }



    /*for (int i = 2; i < path.size(); i++) {
      int ai = path[i-2], bi = path[i-1], ci = path[i];
      for (int li = metai[bi]+1; li < metai[ci]; li++) {
	if (next_layer[metai[bi]][li] < adj_thres) continue;
	int di = extend4(mod, ai, bi, ci, li, 0.1);
	if (di > 0) {
	  path.insert(path.begin()+i-1, di);
	  assignment[di] = k;
	  i--;
	  break;
	}
      }
      }*/
    paths[k] = path;
  }

  {

    initDensity3();

    PolarModule mod[48];
    for (int i = 0; i < 48; i++)
      mod[i] = PolarModule(i);

    vector<int>*tube = readTubes();

    double earn = 0;
    vector<pair<int, int> > pairs = findDuplicates(tube, mod);
    //cout << pairs.size() << endl;
    for (pair<int, int> p : pairs) {
      if (assignment[p.first] || assignment[p.second]) continue;
      /*if (samepart(p.first, p.second)) {
	if (truth_tracks[truth_part[p.first]].size() <= 5) earn += part_weight[truth_part[p.first]];
	//cout << truth_tracks[truth_part[p.first]].size() << ' ' << part_weight[truth_part[p.first]] << endl;
	}*/

      int found = 0;
      {
	int ai = p.first, bi = p.second;
	vector<int> s[48];
	for (int li = 0; li < 48; li++) {
	  if (next_layer[li][metai[ai]] < adj_thres/2 &&
	      next_layer[metai[ai]][li] < adj_thres/2) continue;
	  if (li == metai[ai]) continue;

	  point d, dp, xp, bap;
	  if (prepareTripleScore(ai, bi, li, d, dp, xp, bap) &&
	      prepareTripleScore(bi, ai, li, d, dp, xp, bap)) continue;

	  xp = normalize(xp)*0.99;
	  //xp = point(0,0,0);
	  const double target = 0.5;
	  double mid = findDensity(dp, xp, target, li);
	  int matches = mod[li].getNear(dp, xp, bap, mid, match);

	  double best = target;
	  vector<pair<double, int> > v;
	  for (int i = 0; i < matches; i++) {
	    int ci = match[i];
	    double s = evaluateScore(ci, dp, xp, bap);
	    v.push_back(make_pair(s, ci));
	  }
	  sort(v.begin(), v.end());
	  for (int i = 0; i < v.size(); i++) {
	    if (i >= 1) break;
	    found = v[i].second;
	    //s[li].push_back(v[i].second);
	    //triples.push_back(triple(ai, bi, v[i].second));
	  }
	  if (v.size()) break;
	  //break;
	}
	/*{
	  if (found && samepart(found, ai) && truth_tracks[truth_part[p.first]].size() <= 5) cout << part_weight[truth_part[ai]] << endl;
	  }*/
      }

      vector<int> v;
      v.push_back(p.first);
      v.push_back(p.second);
      v.push_back(found);
      for (int i : v)
	assignment[i] = paths.size();
      paths.push_back(v);
    }
    //cout << "Max earnings: " << earn << endl;
  }
  return paths;
}


//calculate angle deviation from cell's data of the three points in a triple, return the resulting angles in a (missleading) point data structure
point scoreTripleHitDir(triple t) {
  point a = hits[t.x], b = hits[t.y], c = hits[t.z];
  //Find circle with center p, radius r, going through a, b, and c (in xy plane)
  double ax = a.x-c.x, ay = a.y-c.y, bx = b.x-c.x, by = b.y-c.y;
  double aa = ax*ax+ay*ay, bb = bx*bx+by*by;
  double idet = .5/(ax*by-ay*bx);
  point p;
  p.x = (aa*by-bb*ay)*idet;
  p.y = (ax*bb-bx*aa)*idet;
  p.z = 0;
  double r = dist(p.x, p.y), ir = 1./r;
  p.x += c.x;
  p.y += c.y;

  point ca = hits[t.z]-hits[t.x];
  double ang_ca = asin(dist(ca.x, ca.y)*.5*ir)*2;

  double ret[3];
  for (int i = 0; i < 3; i++) {
    int di = i ? i==2 ? t.z : t.y : t.x;

    double rx = hits[di].x-p.x, ry = hits[di].y-p.y;
    double cross = rx*ca.y-ry*ca.x;

    point dir;
    if (ir) {
      dir.x =-ry*ang_ca;
      dir.y = rx*ang_ca;
      dir.z = ca.z;
      if (cross < 0) dir.z *= -1;
    } else {
      dir = ca;
    }
    double angle = 1e9;
    for (int k = 0; k < 2; k++) {
      point&hd = hit_dir[di][k];
      double angle_ = acos(fabs(dir*hd)/dist(dir));
      angle = min(angle, angle_);
    }
    ret[i] = angle;
  }
  return point(ret[0], ret[1], ret[2]);
}

//Prune triples to remove duplicates, and remove more based on deviation from cell's data direction
void pruneTriples(vector<triple>&triples) {
  int start_list_len = 100;
  triple start_list[100] = {{0, 1, 2}, {11, 12, 13}, {4, 5, 6}, {0, 1, 11}, {0, 1, 4}, {0, 4, 5}, {0, 11, 12}, {18, 19, 20}, {1, 2, 3}, {5, 6, 7}, {12, 13, 14}, {13, 14, 15}, {6, 7, 8}, {2, 3, 18}, {3, 18, 19}, {19, 20, 21}, {0, 1, 3}, {0, 2, 3}, {20, 21, 22}, {0, 1, 5}, {0, 1, 12}, {1, 4, 5}, {1, 2, 18}, {11, 18, 19}, {1, 11, 12}, {7, 8, 9}, {4, 18, 19}, {14, 15, 16}, {0, 4, 6}, {18, 19, 24}, {18, 19, 36}, {2, 18, 19}, {0, 18, 19}, {1, 18, 19}, {21, 22, 23}, {11, 12, 18}, {4, 5, 18}, {0, 11, 13}, {0, 1, 18}, {24, 26, 28}, {1, 2, 11}, {1, 2, 4}, {36, 38, 40}, {15, 16, 17}, {8, 9, 10}, {5, 18, 19}, {18, 24, 26}, {18, 36, 38}, {38, 40, 42}, {12, 18, 19}, {0, 12, 13}, {40, 42, 44}, {28, 30, 32}, {18, 20, 21}, {0, 2, 18}, {26, 28, 30}, {4, 5, 7}, {19, 20, 24}, {11, 12, 14}, {19, 20, 36}, {18, 19, 21}, {6, 18, 19}, {13, 18, 19}, {2, 11, 18}, {0, 5, 6}, {0, 2, 11}, {12, 13, 18}, {19, 36, 38}, {4, 6, 7}, {19, 24, 26}, {2, 4, 18}, {0, 2, 4}, {5, 6, 18}, {3, 19, 20}, {11, 13, 14}, {11, 12, 36}, {2, 3, 19}, {1, 3, 18}, {19, 22, 23}, {2, 18, 24}, {4, 5, 24}, {7, 18, 19}, {20, 22, 23}, {24, 26, 29}, {11, 18, 36}, {19, 20, 22}, {20, 21, 37}, {15, 18, 19}, {2, 18, 36}, {14, 18, 19}, {4, 18, 24}, {8, 18, 19}, {3, 18, 20}, {36, 38, 41}, {20, 21, 25}, {9, 18, 19}, {6, 7, 26}, {13, 14, 36}, {42, 44, 46}, {30, 32, 34}};
  set<triple> start_list_set;
  for (int i = 0; i < start_list_len; i++) start_list_set.insert(start_list[i]);

  cout << "Sorting" << endl;
  sort(triples.begin(), triples.end());
  cout << "Pruning" << endl;
  for (int i = 0; i < triples.size(); i++) {
    triple&t = triples[i];
    int ok = 1;
    if (i && t.x == abs(triples[i-1].x) && t.y == triples[i-1].y && t.z == triples[i-1].z) ok = 0;
    else {
      point p = scoreTripleHitDir(t);
      double f = 4;
      if (metai[t.x] < 18) p.x *= f;
      if (metai[t.y] < 18) p.y *= f;
      if (metai[t.z] < 18) p.z *= f;
      double score = p.x+p.y+p.z;
      if (score >= 3.1) ok = 0;
    }
    if (0) {//Not used, only keep triples that are in the beginning of paths
      triple ii(metai[t.x], metai[t.y], metai[t.z]);
      if (!start_list_set.count(ii)) ok = 0;
    }
    if (!ok) {
      t.x *= -1;
    }
  }
  for (int i = 0; i < triples.size(); i++) {
    triple&t = triples[i];
    if (t.x < 0) {
      swap(t, triples[triples.size()-1]);
      triples.pop_back();
      i--;
    }
  }
}




//Deprecated
vector<pair<int, int> > findPairsDir(vector<int>*tube, PolarModule*mod, double target = 80) {
  int top_a[] = {0, 11, 4, 0, 18, 0, 1, 12, 13, 5, 6, 19, 3, 2, 20, 0, 1, 14, 7, 24, 1, 11, 4, 15, 21, 2, 5, 0, 22, 1};
  int top_b[] = {1, 12, 5, 4, 19, 11, 2, 13, 14, 6, 7, 20, 18, 3, 21, 2, 4, 15, 8, 26, 11, 18, 18, 16, 22, 18, 18, 18, 23, 18};

  vector<pair<int, int> > pairs;
  map<pair<int, int>, double> done;
  for (int kk = 0; kk < 24; kk++) {
    load(24, "Find pairs using hit_dir", 0);
    int ii = top_a[kk];
    int jj = top_b[kk];
    //for (int ii = 0; ii < 48; ii++) {
      for (int i : tube[ii]) {
	//load(hits.size()-1, "Find pairs using hit_dir", 0);
	//for (int jj = 0; jj < 48; jj++) {
	//if (ii == jj) continue;

	//double p = max(next_layer[ii][jj], next_layer[jj][ii])*1./count_layer[jj];
	//if (p < 1e-5) continue;
    	//if (max(next_layer[ii][jj], next_layer[jj][ii]) < adj_thres/100) continue;//|| ii >= 18 || jj >= 18) continue; //
	double target0 = target;//*2*sqrt(p);
	for (int k = 0; k < 2; k++) {
	  hits[0] = hits[i]+hit_dir[i][k];
	  polar[0] = point(dist(hits[0].x, hits[0].y), atan2(hits[0].y, hits[0].x), hits[0].z);

	  point d, dp, xp, bap;
	  if (prepareTripleScore(i, 0, jj, d, dp, xp, bap)) continue;
	  xp = point(0,0,0.9);
	  /*
	    double nan = dp*dp+xp*xp+bap*bap;
	    if (nan != nan) {
	    cout << hit_dir[i][k] << ' ' << (layer[jj].type==Tube) << endl;
	    cout << hits[0] << endl;
	    cout << polar[0] << endl;
	    cout << dp << endl;
	    cout << xp << endl;
	    cout << bap << endl;
	    cout << "NaN" << endl;
	    continue;
	    }
	  */

	  int target_ = target0;
	  if (ii < 18) {
	    target_ /= 3;
	    if (jj < 18) target_ *= 2;
	  }


	  double tt = findDensity(dp, xp, target_, jj);

	  int matches = mod[jj].getNear(dp, xp, bap, tt, match);

	  //double z = dist(hits[i]-d);
	  //cout << atan(sqrt(tt)/z) << ' ' << tt << ' ' << z << endl;
	  //if (atan(sqrt(tt)/z) < 0.02) continue;

	  vector<pair<double, int> > v;
	  for (int mi = 0; mi < matches; mi++) {
	    double s = evaluateScore(match[mi], dp, xp, bap);
	    v.push_back(make_pair(s, match[mi]));
	  }
	  sort(v.begin(), v.end());
	  if (v.size() > target_) v.resize(int(target_));
	  for (auto&p : v) {
	    pair<int, int> a = make_pair(i,p.second);
	    pair<int, int> b = make_pair(p.second,i);
	    if (ii >= 18 || jj >= 18) {
	      pair<int, int> add = ii < jj ? a : b;
	      pairs.push_back(add);
	    } else {
	      double s = p.first*target_/tt;
	      if (ii < jj) //Always first
		done[a] = s;
	      else if (done.count(b) && done[b]+s < target_)
		pairs.push_back(a);
	    }
	  }
	}
	//}
    }
  }
  return pairs;
}

//Deprecated
vector<pair<int, int> > findPairsDir2(vector<int>*tube, PolarModule*mod, double target = 80) {
  vector<pair<int, int> > pairs;
  map<pair<int, int>, double> done;
    for (int ii = 0; ii < 48; ii++) {
      for (int i : tube[ii]) {
	load(hits.size()-1, "Find pairs using hit_dir", 0);
	for (int jj = 0; jj < 48; jj++) {
	if (ii == jj) continue;

	//double p = max(next_layer[ii][jj], next_layer[jj][ii])*1./count_layer[jj];
	//if (p < 1e-5) continue;
    	if (max(next_layer[ii][jj], next_layer[jj][ii]) < adj_thres) continue;//|| ii >= 18 || jj >= 18) continue; //
	double target0 = target;//*2*sqrt(p);
	for (int k = 0; k < 2; k++) {
	  hits[0] = hits[i]+hit_dir[i][k];
	  polar[0] = point(dist(hits[0].x, hits[0].y), atan2(hits[0].y, hits[0].x), hits[0].z);

	  point d, dp, xp, bap;
	  if (prepareTripleScore(i, 0, jj, d, dp, xp, bap)) continue;
	  xp = point(0,0,0.9);
	  /*
	    double nan = dp*dp+xp*xp+bap*bap;
	    if (nan != nan) {
	    cout << hit_dir[i][k] << ' ' << (layer[jj].type==Tube) << endl;
	    cout << hits[0] << endl;
	    cout << polar[0] << endl;
	    cout << dp << endl;
	    cout << xp << endl;
	    cout << bap << endl;
	    cout << "NaN" << endl;
	    continue;
	    }
	  */

	  int target_ = target0;
	  if (ii < 18) {
	    target_ /= 3;
	    if (jj < 18) target_ *= 2;
	  }


	  double tt = findDensity(dp, xp, target_, jj);

	  int matches = mod[jj].getNear(dp, xp, bap, tt, match);

	  //double z = dist(hits[i]-d);
	  //cout << atan(sqrt(tt)/z) << ' ' << tt << ' ' << z << endl;
	  //if (atan(sqrt(tt)/z) < 0.02) continue;

	  vector<pair<double, int> > v;
	  for (int mi = 0; mi < matches; mi++) {
	    double s = evaluateScore(match[mi], dp, xp, bap);
	    v.push_back(make_pair(s, match[mi]));
	  }
	  sort(v.begin(), v.end());
	  if (v.size() > target_) v.resize(int(target_));
	  for (auto&p : v) {
	    pair<int, int> a = make_pair(i,p.second);
	    pair<int, int> b = make_pair(p.second,i);
	    if (ii >= 18 || jj >= 18) {
	      pair<int, int> add = ii < jj ? a : b;
	      pairs.push_back(add);
	    } else {
	      double s = p.first*target_/tt;
	      if (ii < jj) //Always first
		done[a] = s;
	      else if (done.count(b) && done[b]+s < target_)
		pairs.push_back(a);
	    }
	  }
	}
	}
    }
  }
  return pairs;
}




//Deprecated
vector<triple> findTriplesDuplicates(vector<pair<int, int> >&pairs, PolarModule*mod) {
  vector<triple> triples;

  double target = 0.1, target2 = 0.2; //target = 0.2 has about everything
  /*
  set<long long> found;
  double impossible = 0;
  for (auto&p : truth_tracks) {
    set<int> metais;
    for (int i : p.second) metais.insert(metai[i]);
    if (metais.size() == p.second.size()) impossible += part_weight[p.first];
  }
  cout << impossible << endl;*/
  int outliers = 0;
  for (auto p : pairs) {
    load(pairs.size(), "Find triples using duplicates");
    int ai = p.first, bi = p.second;
    //if (!samepart(ai, bi)) continue;

    vector<pair<double, int> > top;
    point dir = hits[bi]-hits[ai];
    for (int swap = 0; swap < 2; swap++) {
      int i = swap ? bi : ai;
      int j = swap ? ai : bi;
      point d = hits[i], dp = polar[i], xp(0,0,0);
      int li = metai[i];
      point bap = topolar(dir, d, dp);
      bap = bap*(1./(layer[li].type == Disc ? bap.z : bap.x));
      double tt = findDensity(dp, xp, target, li);
      int matches = mod[li].getNear(dp, xp, bap, tt, match);
      for (int mi = 0; mi < matches; mi++) {
	int k = match[mi];
	if (metaz[k] == metaz[i]) continue;

	/*if (samepart(k, i)) found.insert(truth_part[i]);
	else outliers++;*/

	int xi, yi, zi;
	if (dir.z*(hits[k].z-hits[ai].z) < 0) xi = k, yi = ai, zi = bi;
	else if (dir.z*(hits[k].z-hits[bi].z) > 0) xi = ai, yi = bi, zi = k;
	else xi = ai, yi = k, zi = bi;

	for (int li = metai[bi]+1; li < 48; li++) {
	  if (next_layer[metai[bi]][li] < adj_thres) continue;
	  point d, dp, xp, dirp;
	  if (prepareQuadrupleScore(xi, yi, zi, li, d,dp,xp,dirp)) continue;

	  double tt = findDensity(dp, xp, target2, li);
	  int matches = mod[li].getNear(dp, xp, bap, tt, match);
	  vector<pair<double, int> > v;
	  for (int i = 0; i < matches; i++) {
	    int mi = match[i];
	    v.push_back(make_pair(evaluateScore(mi, dp, xp, bap), mi));
	  }
	  sort(v.begin(), v.end());
	  if (matches)
	    triples.push_back(triple(ai, bi, v[0].second));
	  break;
	}

	for (int li = metai[ai]-1; li >= 0; li--) {
	  if (next_layer[li][metai[ai]] < adj_thres) continue;
	  point d, dp, xp, dirp;
	  if (prepareQuadrupleScore(zi, yi, xi, li, d,dp,xp,dirp)) continue;

	  double tt = findDensity(dp, xp, target2, li);
	  int matches = mod[li].getNear(dp, xp, bap, tt, match);
	  vector<pair<double, int> > v;
	  for (int i = 0; i < matches; i++) {
	    int mi = match[i];
	    v.push_back(make_pair(evaluateScore(mi, dp, xp, bap), mi));
	  }
	  sort(v.begin(), v.end());
	  if (matches)
	    triples.push_back(triple(v[0].second, ai, bi));
	  break;
	}
      }
    }
  }
  //cout << outliers << ' ' << pairs.size()*target*2 << endl;
  /*double score = 0;
  for(long long i : found) score += part_weight[i];
  cout << score << endl;*/
  return triples;
}


vector<pair<int, int>> allPairStarts() {
  vector<pair<int, int>> pairs;
  for (auto&t : truth_tracks) {
    //if (dist(start_pos[t.first].x, start_pos[t.first].y) < 10) continue;
    if (0) {
      int ok = 0;
      map<int, int> c;
      for (int i : t.second)
	c[metai[i]]++;
      for (auto&p : c) if (p.second > 1) ok++;
      if (ok < 1) continue;
    }
    int last = -1;
    vector<int> v;
    for (int i : t.second) {
      if (metai[i] != last) {
	v.push_back(i);
	last = metai[i];
      }
    }
    if (v.size() >= 2) {
      if (metai[v[1]] < metai[v[0]])
	swap(v[0], v[1]);
      //if (metai[v[0]] < 18) continue;
      pairs.push_back(make_pair(v[0], v[1]));
    }
  }
  return pairs;
}

vector<triple> allTripleStarts() {
  vector<triple> triples;
  for (auto&p : truth_tracks) {
    vector<int> r;
    set<int> done;
    for (int i : p.second) {
      int m = metai[i];
      if (done.lower_bound(m) == done.end()) {
	r.push_back(i);
	done.insert(m);
      }
    }
    if (r.size() >= 3) {
      triples.push_back(triple(r[0], r[1], r[2]));
    }
  }
  return triples;
}

int main(int argc, char**argv) {
  //Read which event to run from arguments, default is event # 1000
  //Supports events 0-124 for test, or 1000-1099 for validation (small dataset)
  if (argc >= 2) {
    filenum = atoi(argv[1]);
    cout << "Running on event #" << filenum << endl;
  }
  ios::sync_with_stdio(false);
  cout << fixed;

  int eval = 0;
#if defined EVAL
  eval = 1;
#endif
  if (!eval) {
    readBlacklist();
    readTruth();
    sortTracks();
    readStarts();
  }
  readHits();
  //storeAdjacent();
  loadAdjacent();
  //truthAdjacent();

  readDetectors();
  readCells();
  initHitDir();

  //investigateMissing();

  int step = 1, score = 1;
  if (eval) step = 1, score = 10;

  int iter = 0;
  while (1) {//Actually only running a single iteration
    initDensity3();

    PolarModule mod[48];
    for (int i = 0; i < 48; i++)
      mod[i] = PolarModule(i);

    vector<int>*tube = readTubes();

#if defined TRAIN
    trainNewPairs(tube);
    return 0;
#endif
    //investigateScattering();
    //findDuplicates(tube, mod);


    // Find promising pairs

    vector<pair<int, int> > origin_pairs, pairs;// = allPairStarts();
    if (step <= 1) {
      pairs = findPairsNew(tube, 3);
      //origin_pairs = findPairsNew(tube, 4);
      //pairs = findPairsOld(tube, 0, 0);//2.5+0.5, 0);
      //vector<pair<int, int> > pairs2 = findPairsDir(tube, mod, 150+50);
      //for (auto&p : pairs2) pairs.push_back(p);

      if (!eval && iter == 0) storePairs(pairs);
    } else if (step < 3) pairs = loadPairs();

    if (score <= 1) scorePairs(pairs);


    // Find promising triples
    vector<triple> triples;
    if (step <= 2) {
      {
	vector<triple> straight_triples = findTriples(pairs, mod, 1);
	for (auto&t : straight_triples) triples.push_back(t);
      }
      //scoreOutsideTriples(triples);
      //scoreOriginTriples(triples);
      /*{
	vector<triple> origin_triples = findTriples(origin_pairs, mod, 0);
	for (auto&t : origin_triples) triples.push_back(t);
      }
      scoreOutsideTriples(triples);
      scoreOriginTriples(triples);*/
      /*{
	vector<triple> origin2_triples = findTriples(pairs, mod, 0, 0.3);
	for (auto&t : origin2_triples) triples.push_back(t);
      }
      scoreOutsideTriples(triples);
      scoreOriginTriples(triples);*/
      //return 0;
      {
	vector<triple> duplicate_triples = findTriplesDuplicates(tube, mod);//findTriplesDuplicates(pairs, mod);
	//scoreTriples(triples2);
	for (auto&t : duplicate_triples) triples.push_back(t);
      }
      //scoreOutsideTriples(triples);
      //scoreOriginTriples(triples);
      //scoreTriples(triples);
      if (!eval && iter == 0) storeTriples(triples);
      pruneTriples(triples);
    } else if (step < 4) {
      triples = loadTriples();
      pruneTriples(triples);
    }


    /*{
      set<long long> missing;
      for (auto&p : part_weight) missing.insert(p.first);
      for (triple t : triples)
	if (samepart(t.x, t.y) && samepart(t.y,t.z) && missing.count(truth_part[t.x]))
	  missing.erase(truth_part[t.x]);
      FILE*fp = fopen("missing", "w");
      int n = missing.size();
      fwrite(&n, sizeof(n), 1, fp);
      for (long long t : missing)
	fwrite(&t, sizeof(t), 1, fp);
      fclose(fp);
      exit(0);
    }
    */
    //for (const triple&t : triples) {
    //  assert(metai[t.x] < metai[t.y] && metai[t.y] < metai[t.z], "Wrong order!");
    //}
    //dumpPrunedTriples(triples);

    //triples = allTripleStarts();

    if (score <= 2) scoreTriples(triples);
    pairs.clear();


    // Extend promising triples to paths

    vector<vector<int> > paths;
    if (step <= 3) {
      paths = findPaths(triples, mod);
      if (!eval && iter == 0) storePaths(paths, 1);
    } else if (step == 4) paths = loadPaths(1);

    if (score <= 3) scorePaths4(paths);

    triples.clear();


    // Prune paths

    if (step <= 4) {
      paths = prunePaths(paths);
      if (!eval && iter == 0) storePaths(paths, 2);
    } else if (step == 5) paths = loadPaths(2);


    //Add duplicates to paths

    if (step <= 5) {
      paths = addDuplicates(paths, mod);
      if (!eval && iter == 0) storePaths(paths, 3);
    } else paths = loadPaths(3);

    if (score <= 5) scorePaths(paths);


    //Find hit assignments to paths

    if (step <= 6) {
      paths = findAssignment(paths, mod);
      if (!eval && iter == 0) storePaths(paths, 4);
    } else paths = loadPaths(4);


    map<int, int> assignment_part;
    for (int i = 1; i < paths.size(); i++)
      for (int j : paths[i])
	if (j > 0)
	  assignment_part[j] = i;
	//assignment_part = loadAssignment();

    for (pair<int, int> p : assignment_part)
      assignment[p.first] = p.second;


    //Gather a little more score (~0.001) from tricks on remaining points

    //investigateAssignment(paths);
    paths = extendPaths(paths);
    //investigateAssignment(paths);


    map<int, int> map_assignment;
    for (int i = 1; i < hits.size(); i++)
      if (assignment[i]) map_assignment[i] = assignment[i];

    if (score <= 6) {
      scoreAssignment(map_assignment);
      scoreAssignment2(map_assignment);
    }

    //Write out solution to ./submissions/submission<event id>.csv
    writeSubmission(map_assignment);

    iter++;
    break;
  }
}
