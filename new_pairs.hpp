int good_pair(int a, int b) {
  if (!samepart(a, b)) return 0;
  point s = start_pos[truth_part[a]];
  if (s.x*s.x+s.y*s.y < 100) return 2;
  auto&v = truth_tracks[truth_part[a]];
  int ai = metai[ai], bi = metai[bi];
  if (ai > bi) swap(ai, bi);
  for (int i : v)
    if (metai[i] < ai || metai[i] > ai && metai[i] < bi) return 2;
  return 1;
}

double dir_miss(int ai, int bi) {
  point d = hits[ai]-hits[bi];
  return acos(max(fabs(hit_dir[ai][0]*d),
		  fabs(hit_dir[ai][1]*d))/dist(d));
}


void getFeatures3(int ai, int bi,  double*feature) {
  int x = metai[ai], y = metai[bi];
  point a = hits[ai], b = hits[bi];
  point d = a-b;

  double dr2 = dist(d.x, d.y);
  double r1 = dist(a.x, a.y);
  double r2 = dist(b.x, b.y);
  double dr = r2-r1;

  feature[0] = dir_miss(ai, bi);//bla
  feature[1] = dir_miss(bi, ai);//5
  feature[2] = wdist(a, d, 0);//4
  feature[3] = zdist2(a, b);//2
  feature[4] = wdistr(r1, dr, a.z, d.z, 1);//6
  feature[5] = wdist(a, d, 1);//0
}


void investigateNewPairs(vector<int>*tube) {
  const int n = 50;
  pair<int, int> start_list[100] = {{0, 1}, {11, 12}, {4, 5}, {0, 4}, {0, 11}, {18, 19}, {1, 2}, {5, 6}, {12, 13}, {13, 14}, {6, 7}, {2, 3}, {3, 18}, {19, 20}, {0, 2}, {20, 21}, {1, 4}, {7, 8}, {11, 18}, {1, 11}, {14, 15}, {4, 18}, {2, 18}, {21, 22}, {0, 18}, {1, 18}, {24, 26}, {36, 38}, {15, 16}, {8, 9}, {22, 23}, {9, 10}, {16, 17}, {38, 40}, {5, 18}, {18, 24}, {18, 36}, {12, 18}, {40, 42}, {28, 30}, {26, 28}, {0, 12}, {18, 20}, {6, 18}, {2, 11}, {13, 18}, {2, 4}, {0, 5}, {19, 36}, {19, 24}, {4, 6}, {19, 22}, {20, 22}, {11, 13}, {3, 19}, {7, 18}, {14, 18}, {3, 4}, {22, 25}, {1, 3}, {20, 24}, {15, 18}, {3, 11}, {22, 37}, {30, 32}, {42, 44}, {8, 18}, {9, 18}, {8, 26}, {15, 38}, {20, 36}, {14, 36}, {7, 24}, {1, 5}, {16, 18}, {22, 24}, {18, 22}, {25, 27}, {16, 40}, {10, 30}, {25, 26}, {17, 40}, {36, 39}, {1, 12}, {10, 28}, {7, 26}, {17, 42}, {24, 27}, {21, 24}, {23, 37}, {13, 36}, {15, 36}, {22, 36}, {14, 38}, {8, 28}, {19, 21}, {6, 24}, {9, 28}, {16, 38}, {0, 3}};

  const int features = 6;
  double feature[features];

  for (int i = 0; i < n; i++) {
    load(n);
    int goods = 0;
    char filename[100];
    sprintf(filename, "data2/data%d", i);
    FILE*fp = fopen(filename, "w");
    for (auto a : tube[start_list[i].first]) {
      for (auto b : tube[start_list[i].second]) {
	double dot = hits[a].x*hits[b].x+hits[a].y*hits[b].y;
	double alen = dist2(hits[a].x, hits[a].y);
	double blen = dist2(hits[b].x, hits[b].y);
	if (dot < 0 || dot*dot < alen*blen*(.7*.7)) continue;

	dot += hits[a].z*hits[b].z;
	alen += hits[a].z*hits[a].z;
	blen += hits[b].z*hits[b].z;
	if (dot < 0 || dot*dot < alen*blen*(.7*.7)) continue;

	int g = good_pair(a, b);
	if (rand()%1 == 0 && (g == 1 || g == 0 && rand()%100 == 0)) {

	  getFeatures3(a, b, feature);
	  fprintf(fp, "%d", g==1);
	  goods += g==1;
	  for (int j = 0; j < features; j++) fprintf(fp, " %lf", feature[j]);
	  fprintf(fp, "\n");
	}
      }
    }
    //if (goods < 100
    fclose(fp);
  }
  exit(0);
}

vector<pair<int, int> > findPairsNew(vector<int>*tube, int modeli = 2) {
  const int features = 6;
  double feature[features];

  const int n = 50;//TODO
  pair<int, int> start_list[100] = {{0, 1}, {11, 12}, {4, 5}, {0, 4}, {0, 11}, {18, 19}, {1, 2}, {5, 6}, {12, 13}, {13, 14}, {6, 7}, {2, 3}, {3, 18}, {19, 20}, {0, 2}, {20, 21}, {1, 4}, {7, 8}, {11, 18}, {1, 11}, {14, 15}, {4, 18}, {2, 18}, {21, 22}, {0, 18}, {1, 18}, {24, 26}, {36, 38}, {15, 16}, {8, 9}, {22, 23}, {9, 10}, {16, 17}, {38, 40}, {5, 18}, {18, 24}, {18, 36}, {12, 18}, {40, 42}, {28, 30}, {26, 28}, {0, 12}, {18, 20}, {6, 18}, {2, 11}, {13, 18}, {2, 4}, {0, 5}, {19, 36}, {19, 24}, {4, 6}, {19, 22}, {20, 22}, {11, 13}, {3, 19}, {7, 18}, {14, 18}, {3, 4}, {22, 25}, {1, 3}, {20, 24}, {15, 18}, {3, 11}, {22, 37}, {30, 32}, {42, 44}, {8, 18}, {9, 18}, {8, 26}, {15, 38}, {20, 36}, {14, 36}, {7, 24}, {1, 5}, {16, 18}, {22, 24}, {18, 22}, {25, 27}, {16, 40}, {10, 30}, {25, 26}, {17, 40}, {36, 39}, {1, 12}, {10, 28}, {7, 26}, {17, 42}, {24, 27}, {21, 24}, {23, 37}, {13, 36}, {15, 36}, {22, 36}, {14, 38}, {8, 28}, {19, 21}, {6, 24}, {9, 28}, {16, 38}, {0, 3}};

  vector<pair<int, int> > pairs;
  for (int i = 0; i < n; i++) {
    load(n, "Finding pairs");
    double weight[23];
    {
      char filename[100];
      sprintf(filename, "pair_logreg%d/model%d", modeli, i);
      FILE*fp = fopen(filename, "r");
      if (!fp) {
	//cout << "Could not open " << filename << endl;
	continue;
      }
      for (int i = 0; i < 23; i++) int tmp = fscanf(fp, "%lf", &weight[i]);
      fclose(fp);
    }

    for (auto a : tube[start_list[i].first]) {
      for (auto b : tube[start_list[i].second]) {
	double dot = hits[a].x*hits[b].x+hits[a].y*hits[b].y;
	double alen = dist2(hits[a].x, hits[a].y);
	double blen = dist2(hits[b].x, hits[b].y);
	if (dot < 0 || dot*dot < alen*blen*(.7*.7)) continue;

	dot += hits[a].z*hits[b].z;
	alen += hits[a].z*hits[a].z;
	blen += hits[b].z*hits[b].z;
	if (dot < 0 || dot*dot < alen*blen*(.7*.7)) continue;

	getFeatures3(a, b, feature);
	double&A = feature[0], &B = feature[1], &C = feature[2], &D = feature[3], &E = feature[4], &F = feature[5];
	double score = weight[0]+
		  A*(weight[1]+
		     B*weight[7]+
		     C*weight[8]+
		     D*weight[9]+
		     E*weight[10]+
		     F*weight[11])+
		  B*(weight[2]+
		     C*weight[12]+
		     D*weight[13]+
		     E*weight[14]+
		     F*weight[15])+
		  C*(weight[3]+
		     D*weight[16]+
		     E*weight[17]+
		     F*weight[18])+
		  D*(weight[4]+
		     E*weight[19]+
		     F*weight[20])+
		  E*(weight[5]+
		     F*weight[21])+
		  F*weight[6];
	if (score > weight[22])
	  pairs.push_back(make_pair(a, b));
      }
    }
  }
  return pairs;
}



double scoreTripleLogRadius_and_HitDir(int ai, int bi, int ci, double*L) {
  point a = hits[ai], b = hits[bi], c = hits[ci];
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

  for (int k = 0; k < 3; k++) {
    int di = k ? k==2 ? ci : bi : ai;
    double rx = hits[di].x-p.x, ry = hits[di].y-p.y;

    point ca = hits[ci]-hits[ai];
    double ang_ca = asin(dist(ca.x, ca.y)*.5*ir)*2;
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
    L[k] = acos(max(fabs(hit_dir[di][0]*dir),
		    fabs(hit_dir[di][1]*dir))/dist(dir));
  }
  return log(ir);
}

void investigatePruneTriples(vector<triple>&triples) {
  FILE*fp = fopen("triples", "r");
  int n;
  int tmp = fscanf(fp, "%d", &n);
  triples.resize(n);
  for (int i = 0; i < triples.size(); i++) {
    load(triples.size());
    tmp = fscanf(fp, "%d%d%d", &triples[i].x, &triples[i].y, &triples[i].z);
    triple&t = triples[i];

    int good = samepart(t.x,t.y) && samepart(t.y,t.z);
    if (rand()%20 != 0) continue;
    if (!good && rand()%100 != 0) continue;
    double A = log(scoreTriple(t.x, t.y, t.z)+1e-8);
    if (A != A) {
      cout << t.x << ' ' << t.y << ' ' << t.z << endl;
      exit(0);
    }
    double B = log(scoreTripleDensity(t.x, t.y, t.z));
    double C = log(scoreTripleDensity(t.z, t.y, t.x));
    double L[3];
    double D = scoreTripleLogRadius_and_HitDir(t.x,t.y,t.z,L);
    cerr << good << ' ' << A << ' ' << B << ' ' << C << ' ' << D << ' ' << L[0] << ' ' << L[1] << ' ' << L[2] << endl;
  }
  fclose(fp);
  exit(0);
}

int acceptTriple(const triple&t) {
  double A = log(scoreTriple(t.x, t.y, t.z)+1e-8);
  double B = log(scoreTripleDensity(t.x, t.y, t.z));
  double C = log(scoreTripleDensity(t.z, t.y, t.x));
  double L[3];
  double D = scoreTripleLogRadius_and_HitDir(t.x,t.y,t.z,L);
  double w[121] = {-13.215482291516638, -0.519368174205195, -0.6019168737814719, -0.400773825827796, -3.0689189279504614, -8.21987444638849, -1.773083608093787, -3.7271459966647913, -0.18753136282696767, -0.1700350202416788, -0.13325020734065293, -0.0712787103124509, 2.2365502305889295, -0.38264699613950004, -1.5996361946235698, 0.02602607302842127, -0.04090477074387659, -0.12800207114108786, -2.0616314706262706, 0.9350417490331662, -0.6313327964001432, 0.00830034532077729, -0.1021716887039019, 0.3719980432444666, 0.43818671158350325, 0.0338130884608543, 0.19225191422472998, -0.33226136371623366, -1.0631299954951279, -1.3290857109832128, 8.50340840417112, 4.489339976769724, -3.6866359902477703, -1.530677908510526, -0.3660375991432235, -0.2832850515900752, -0.003067393550643885, -0.06185860378584967, -0.004472073355177509, -0.034047188478205974, 0.056232208303619684, -0.09251101374546467, -0.3186456107148592, -0.011497406815609599, 0.0040898730087192275, 0.04166475101451824, 0.5313081554181062, 0.05691563704023761, 0.004054188315119864, 0.009440976068230187, 0.015452389083207108, 0.02857025202533131, -0.01788978369714811, -0.014820867725080077, -0.0032975179225221054, -0.2739810756530984, -0.209895536224461, -0.05555070596049059, -3.8393234795148445, -0.39189992715019867, 0.5302884318217037, -1.0560724732243318, 0.5808249742500916, 0.2085127159157602, -0.002879796716268462, -0.008289453513497825, -0.013327308424637882, 0.034516052559319284, 0.05612738574267425, -0.04698958101602463, 0.0007407605230615924, -0.015547995524776616, 0.06280040184070336, -0.056422842974113374, -0.02553695075115984, -0.030162351232030156, -0.216209409546151, 0.03852063554031595, -0.0693834129966963, -1.0570960495204662, 0.6811799884934827, 0.3386224510850844, -0.10244400357684635, -0.17437169642288441, 0.527447777429105, -0.0009197072806774356, -0.004512546596535816, -0.026048615023175962, -0.016165328699534447, -0.007957908851240184, -0.01677913671380496, 0.00448514125782629, -0.0164129789525374, -0.04792927265651915, 0.3459064488723725, 0.08305188504334206, -0.4177214300084773, -0.09227473501037928, 0.04508615512899353, -0.03988016215006392, 0.029600325479286028, -0.2533468783999991, -0.1438693183062194, -0.17942937900359165, -1.277174888294048, -0.12050721012197445, -1.306910361564254, -0.056617003726385146, -1.1681337555296898, 0.06259298498866638, -6.501290522349262, -10.841239719611016, 2.156866020752887, 1.3871764445557901, 4.945464722802966, -4.26463890108575, -1.510051189434741, -3.140021739172429, -5.693045331329942, 1.1610084032964856, 2.2204604560570425};
  double x[8] = {1,A,B,C,D,L[0],L[1],L[2]};
  int c = 0;
  double score = 0;
  for (int i = 0; i < 8; i++) {
    for (int j = i; j < 8; j++) {
      double a = x[i]*x[j];
      for (int k = j; k < 8; k++)
	score += a*x[k]*w[c++];
    }
  }
  //cout << score << ' ' << w[120] << endl;
  return score > w[120];
}
