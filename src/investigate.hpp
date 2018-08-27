//Functions used for investigations into the data and debugging, not important

double feature1(int a, int b) {
  double r = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      point dir = hit_dir[i ? a : b][j];
      r = max(r, fabs(dir.z));//fabs(hits[b].z-hits[a].z));
    }
  }
  return r;
}

void findDuplicatesInvestigate(vector<int>*tube, PolarModule*mod) {
  int li = 4;
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
  vector<pair<int, int>> good, bad;
  for (auto&i : v) {
    int a = i.first.first, b = i.first.second;
    if (a < b || !v.count(make_pair(b,a))) {
      //if (!(v[make_pair(b,a)]+v[make_pair(a,b)] < 0.12)) continue;
      if (samepart(a, b)) good.push_back(i.first),found.insert(truth_part[a]);
      else fails++, bad.push_back(i.first);
    }
  }
  {
    vector<double> g, b;
    for (auto&p : good) g.push_back(feature1(p.first, p.second));
    for (auto&p : bad) b.push_back(feature1(p.first, p.second));
    sort(g.begin(), g.end());
    sort(b.begin(), b.end());
    for (int i = 0; i < g.size(); i++)
      cerr << i*1./g.size() << ' ' << g[i] << endl;
    cerr << endl << endl;
    for (int i = 0; i < b.size(); i++)
      cerr << i*1./b.size() << ' ' << b[i] << endl;
  }
  exit(0);

  cout << fails << endl;

  double w = 0;
  for (long long i : found) w += part_weight[i];
  cout << w << endl;

  double full = 0;
  for (auto&p : truth_tracks) {
    int c = 0;
    for (int i : p.second) c += metai[i]==li;
    if (c > 1) full += part_weight[p.first];
  }
  cout << full << endl;
  exit(0);
}





//Validate prepareQuadrupleScore for sign = -1
void investigateScattering() {
  //int A = 0, B = 1, C = 3, D = 2;
  int A = 4, B = 5, C = 7, D = 6;

  vector<int> Di;
  for (int i = 1; i < hits.size(); i++)
    if (metai[i] == D) Di.push_back(i);

  vector<pair<double, int> > allscores;
  int sum_outliers = 0, inliers = 0;
  double expected = 0;
  for (auto&t : truth_tracks) {
    point s = start_pos[t.first];
    for (int i : t.second) {
      if (metai[i] != A) continue;
      for (int j : t.second) {
	if (metai[j] != B) continue;
	for (int k : t.second) {
	  if (metai[k] != C) continue;
	  point d, dp, xp, bap;
	  if (prepareQuadrupleScore(i, j, k, D, d, dp, xp, bap, point(layer[D].avgr, 0, layer[D].avgz), -1)) {
	    continue;
	  }

	  vector<pair<double, int> > score;
	  for (int l : Di) {
	    //if (prepareQuadrupleScore(i, j, k, D, d, dp, xp, bap, polar[l])) continue;
	    double s = evaluateScore(l, dp, xp, bap);
	    score.push_back(make_pair(s, l));
	  }
	  /*for (int l : Di) {
	    score.push_back(make_pair(scoreQuadruple2(i,j,k,l), l));
	    }*/
	  sort(score.begin(), score.end());
	  int outliers = 0, cc = 0;
	  for (auto&p : score) {
	    if (cc++ == 100) break;
	    if (!samepart(i,p.second)) {
	      outliers++;
	      allscores.push_back(make_pair(p.first, 0));
	    } else {
	      expected += p.first;//p.first*p.first*M_PI*getDensity(hits[p.second], metai[p.second]);
	      inliers++;
	      sum_outliers += outliers;
	      allscores.push_back(make_pair(p.first, 1));
	    }
	  }
	}
      }
    }
  }
  cout << inliers << ' ' << sum_outliers << endl;
  cout << expected << endl;
  cout << "Average outliers added per inlier: " << sum_outliers*1./inliers << endl;

  {
    sort(allscores.begin(), allscores.end());
    int inliers = 0, outliers = 0;
    for (auto i : allscores) {
      if (i.second) {
	inliers++;
	cerr << inliers << ' ' << outliers << endl;//log(outliers+1)/log(10) << endl;
      } else outliers++;
    }
  }
  exit(0);
}

//Investigate which particles we are missing
void investigateMissing() {
  vector<long long> missing;
  int n;
  FILE*fp = fopen("missing", "r");
  int tmp = fread(&n, sizeof(n), 1, fp);
  for (int i = 0; i < n; i++) {
    long long p;
    int tmp = fread(&p, sizeof(missing[i]), 1, fp);
    if (part_weight[p])
      missing.push_back(p);
  }
  fclose(fp);

  vector<long long> all;
  for (auto&p : truth_tracks)
    if (part_weight[p.first])
      all.push_back(p.first);

  /*map<int, int> a, b;
  for (long long p : all)
    a[truth_tracks[p].size()]++;
  for (long long p : missing)
    b[truth_tracks[p].size()]++;
  for (auto&p : a)
    cerr << p.first << ' ' << p.second << endl;
  cerr << endl << endl;
  for (auto&p : b)
  cerr << p.first << ' ' << p.second << endl;*/

  vector<double> a,b;
  for (long long p : all) {
    point s = start_mom[p];
    double v = log(dist(s));
    a.push_back(v);
  }
  for (long long p : missing) {
    point s = start_mom[p];
    double v = log(dist(s));
    b.push_back(v);
  }
  sort(a.begin(), a.end());
  sort(b.begin(), b.end());
  for (int i = 0; i < a.size(); i++)
    cerr << i*1./a.size() << ' ' << a[i] << endl;
  cerr << endl << endl;
  for (int i = 0; i < b.size(); i++)
    cerr << i*1./b.size() << ' ' << b[i] << endl;
    /*r (long long p : all) {
    point s = start_mom[p];
    double v = dist(s.x, s.y);
    if (v < 1e3)
      cerr << v << endl;
      }*/


  exit(0);
}
