//Find circle with center p, radius r, going through a, b, and c (in xy plane)
inline void circle(point&a, point&b, point&c, point&p, double&r) {
  double ax = a.x-c.x, ay = a.y-c.y, bx = b.x-c.x, by = b.y-c.y;
  double aa = ax*ax+ay*ay, bb = bx*bx+by*by;
  double idet = .5/(ax*by-ay*bx);
  p.x = (aa*by-bb*ay)*idet;
  p.y = (ax*bb-bx*aa)*idet;
  p.z = 0;
  r = dist(p.x, p.y);
  p.x += c.x;
  p.y += c.y;
}

//Different distances from origin (like how far does the line through ai-bi pass from the origin)
double wdistr(double r1, double dr, double az, double dz, double w) {
  double pp = r1*r1+az*az*w;
  double pd = r1*dr+az*dz*w;
  double dd = dr*dr+dz*dz*w;
  return sqrt(pp-pd*pd/dd);
}
double wdist(point&a, point&d, double w) {
  double pp = a.x*a.x+a.y*a.y+a.z*a.z*w;
  double pd = a.x*d.x+a.y*d.y+a.z*d.z*w;
  double dd = d.x*d.x+d.y*d.y+d.z*d.z*w;
  return sqrt(pp-pd*pd/dd);
}

double zdist(point&a, point&b) {
  static point origin(0,0,0);
  point p;
  double r;
  circle(origin, a, b, p, r);
  double ang_ab = 2*asin(dist(a.x-b.x, a.y-b.y)*.5/r);
  double ang_a = 2*asin(dist(a.x, a.y)*.5/r);
  return fabs(a.z-(b.z-a.z)*ang_a/ang_ab);
}

double zdist2(point&a, point&b) {
  static point origin(0,0,0);
  point p;
  double r;
  circle(origin, a, b, p, r);
  double ang_ab = 2*asin(dist(a.x-b.x, a.y-b.y)*.5/r);
  double ang_a = 2*asin(dist(a.x, a.y)*.5/r);

  return fabs(b.z-a.z-a.z*ang_ab/ang_a);
}

//Deprecated
void getFeatures2(int ai, int bi,  double*feature) {
  int x = metai[ai], y = metai[bi];
  point a = hits[ai], b = hits[bi];
  point d = a-b;

  double dr2 = dist(d.x, d.y);
  double r1 = dist(a.x, a.y);
  double r2 = dist(b.x, b.y);
  double dr = r2-r1;

  feature[0] = wdist(a, d, 1);
  feature[1] = wdist(a, d, 0);
  feature[2] = zdist2(a, b);
  feature[3] = wdistr(r1, dr, a.z, d.z, 1);
  //feature[4] = zdist(a, b);
}



//Deprecated
vector<pair<int, int> > findPairsOld(vector<int>*tube, double off = 0.5, double grad = 1) {
  int top_a[] = {0, 11, 4, 0, 18, 0, 1, 12, 13, 5, 6, 19, 3, 2, 20, 0, 1, 14, 7, 24, 1, 11, 4, 15, 21, 2, 5, 0, 22, 1};
  int top_b[] = {1, 12, 5, 4, 19, 11, 2, 13, 14, 6, 7, 20, 18, 3, 21, 2, 4, 15, 8, 26, 11, 18, 18, 16, 22, 18, 18, 18, 23, 18};
  //for (int i = 0; i < 30; i++) cout << top_a[i] << ' ' << top_b[i] << endl;
  const int features = 4;
  double feature[features];

  /*
  for (int i = 0; i < 30; i++) {
    char filename[100];
    sprintf(filename, "data/data%d", i);
    FILE*fp = fopen(filename, "w");
    for (auto a : tube[top_a[i]]) {
      for (auto b : tube[top_b[i]]) {

	if (hits[a].z*hits[b].z < 0) continue;
	double dot = hits[a].x*hits[b].x+hits[a].y*hits[b].y;
	double alen = dist2(hits[a].x, hits[a].y);
	double blen = dist2(hits[b].x, hits[b].y);
	if (dot < 0 || dot*dot < alen*blen*(.8*.8)) continue;

	dot += hits[a].z*hits[b].z;
	alen += hits[a].z*hits[a].z;
	blen += hits[b].z*hits[b].z;
	if (dot < 0 || dot*dot < alen*blen*(.8*.8)) continue;


	if (rand()%1 == 0 && (samepart(a, b) || rand()%10 == 0)) {

	  getFeatures2(a, b, feature);
	  fprintf(fp, "%d", samepart(a, b));
	  for (int j = 0; j < features; j++) fprintf(fp, " %lf", feature[j]);
	  fprintf(fp, "\n");
	}
      }
    }
    fclose(fp);
    cout << i << endl;
  }
  exit(0);
  */

  int tot = 0;
  //4 is bad, 6 is decent,11, 12, 13, 14, 15, (19), 24, 25, 27, 28
  vector<pair<int, int> > pairs;
  for (int i = 0; i < 24; i++) {
    double weight[17];
    {
      char filename[100];
      sprintf(filename, "pair_logreg/model%d", i);
      FILE*fp = fopen(filename, "r");
      if (!fp) cout << "Could not open " << filename << endl;
      for (int i = 0; i < 17; i++) int tmp = fscanf(fp, "%lf", &weight[i]);
      fclose(fp);
    }
    vector<double> v;
    int c = 1, adds = 0, goods = 0;
    for (auto a : tube[top_a[i]]) {
      for (auto b : tube[top_b[i]]) {
	if (hits[a].z*hits[b].z < 0) continue;
	double dot = hits[a].x*hits[b].x+hits[a].y*hits[b].y;
	double alen = dist2(hits[a].x, hits[a].y);
	double blen = dist2(hits[b].x, hits[b].y);
	if (dot < 0 || dot*dot < alen*blen*(.8*.8)) continue;

	dot += hits[a].z*hits[b].z;
	alen += hits[a].z*hits[a].z;
	blen += hits[b].z*hits[b].z;
	if (dot < 0 || dot*dot < alen*blen*(.8*.8)) continue;

	/*if (1 || samepart(a, b))
	  pairs.push_back(make_pair(a, b));
	  continue;*/

	/*if (samepart(a, b) || rand()%100 == 0) {
	  cerr << samepart(a, b) << ' ';
	  likely2(a, b);
	  }*/
	getFeatures2(a, b, feature);
	double&A = feature[0], &B = feature[1], &C = feature[2], &D = feature[3];
	double score = weight[0]+
	  A*(weight[1]+
	     B*(weight[5]+
		C*(weight[11]+
		   D*weight[15])+
		D*weight[12])+
	     C*(weight[6]+
		D*weight[13])+
	     D*weight[7])+
	  B*(weight[2]+
	     C*(weight[8]+
		D*weight[14])+
	     D*weight[9])+
	  C*(weight[3]+
	     D*weight[10])+
	  D*weight[4];

	  //if (samepart(a, b)) {+(i==4||i==11||i==14)*.5
	if (score > (i*i*.01+i*.1)*grad-off) { //weight[16]+i*.1) {//1+i*.2) {//weight[16]) {
	  pairs.push_back(make_pair(a, b));
	  adds++;
	}
	//if (c++%1000000 == 0) cout << c << " / " << tube[top_a[i]].size()*tube[top_b[i]].size() << endl;
      }
    }
    cout << i << ' ' << adds << ' ' << tube[top_a[i]].size()*tube[top_b[i]].size() << endl;
    //sort(v.rbegin(), v.rend());
    //cout << v[goods*20] << "," << endl;
    //cout << goods << ' ' << adds << endl;
    //cout << pairs.size() << endl;
  }
  cout << tot << endl;
  return pairs;
}
