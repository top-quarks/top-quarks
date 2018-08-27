//Contains functions to read the input data, and global variables to hold the data

//Where to find all the data
const char base_path[] = "/data/kaggle";
//const char base_path[] = ".";

//Which event to run, this may be overwritten by main()'s arguments
int filenum = 1000;

map<long long, vector<int> > truth_tracks; //truth hit ids in each track
point truth_pos[200000], truth_mom[200000]; //truth position and momentum
double truth_weight[200000]; //weighting of each hit
long long truth_part[200000]; //particle this hit belongs to
map<long long, double> part_weight; //weighting of each particle
map<long long, map<int, double> > metai_weight; //weighting of each particle hit, also adding duplicates

//Empirical field strengh, to scale the momentum
const double Bfield = 1673;

set<long long> blacklist;
void readBlacklist() {
  if (filenum < 1000) return;
  char file[1000];
  sprintf(file, "%s/blacklist/event%09d-blacklist_particles.csv", base_path, filenum);
  FILE*fp = fopen(file, "r");
  if (!fp) { printf("couldn't open blacklist\n"); return; }
  char tmpstr[1000];
  int tmp = fscanf(fp, "%s", tmpstr);
  long long particle_id;
  while (fscanf(fp, "%lld", &particle_id) == 1) {
    blacklist.insert(particle_id);
  }
  fclose(fp);
}

void readTruth() {
  if (filenum < 1000) return;
  char file[1000];
  sprintf(file, "%s/train_100_events/event%09d-truth.csv", base_path, filenum);
  FILE*fp = fopen(file, "r");
  if (!fp) { printf("couldn't open truth\n"); return; }
  char tmpstr[1000];
  int tmp = fscanf(fp, "%s", tmpstr);
  while (1) {
    int hit_id;
    long long particle_id;
    double tx, ty, tz, tpx, tpy, tpz, weight;
    if (fscanf(fp, "%d,%lld,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &hit_id, &particle_id, &tx, &ty, &tz, &tpx, &tpy, &tpz, &weight) != 9) break;
    if (!particle_id || blacklist.count(particle_id)) continue;

    truth_tracks[particle_id].push_back(hit_id);
    truth_pos[hit_id] = point(tx, ty, tz);
    truth_mom[hit_id] = point(tpx, tpy, tpz)*Bfield;
    truth_weight[hit_id] = weight;
    part_weight[particle_id] += weight;
    truth_part[hit_id] = particle_id;
  }
  fclose(fp);

  cout << truth_tracks.size() << " particles with truth" << endl;
}

map<long long, point> start_pos; //start position
map<long long, point> start_mom; //start momentum
map<long long, int> part_q; //start charge
map<long long, int> part_hits; // = truth_tracks[particle_id].size()

void readStarts() {
  char file[1000];
  sprintf(file, "%s/train_100_events/event%09d-particles.csv", base_path, filenum);
  FILE*fp = fopen(file, "r");
  if (!fp) { printf("couldn't open particles\n"); return; }
  char tmpstr[1000];
  int tmp = fscanf(fp, "%s", tmpstr);
  while (1) {
    long long id;
    point p, m;
    int hits;
    int q;
    if (fscanf(fp, "%lld,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d", &id, &p.x, &p.y, &p.z, &m.x, &m.y, &m.z, &q, &hits) == -1) break;
    if (blacklist.count(id)) continue;
    start_pos[id] = p;
    start_mom[id] = m*Bfield;
    part_q[id] = -q;
    part_hits[id] = hits;
  }
  fclose(fp);
  cout << start_pos.size() << " particles" << endl;
}


bool track_cmp(int a, int b) {
  point&ma = truth_mom[a];
  point&mb = truth_mom[b];
  double va = ma*ma, vb = mb*mb;
  if (fabs((va-vb)/va) > 1e-5) return va > vb;
  return (truth_pos[b]-truth_pos[a])*ma > 0;
}

//Sort the hits in each track chronologically
void sortTracks() {
  int fails = 0, goods = 0, failparts = 0;
  for (auto&p : truth_tracks) {
    vector<int> v;
    for (int hit_id : p.second) {
      v.push_back(hit_id);
    }
    sort(v.begin(), v.end(), track_cmp);
    for (int i = 0; i < v.size(); i++)
      p.second[i] = v[i];
    int bad = 0;
    for (int i = 2; i < v.size(); i++) {
      point&a = truth_pos[p.second[i-2]];
      point&b = truth_pos[p.second[i-1]];
      point&c = truth_pos[p.second[i]];
      point&ma = truth_mom[p.second[i-2]];
      point&mb = truth_mom[p.second[i-1]];
      point&mc = truth_mom[p.second[i]];
      if ((c.z-b.z)*(b.z-a.z) < 0) {
	/*cout << endl;
	cout << v[i-2].first << ' ' << v[i-1].first << ' ' << v[i].first << endl;
	cout << a.z << ' ' << b.z << ' ' << c.z << endl;
	cout << mb << endl;
	cout << ma.z << ' ' << mb.z << ' ' << mc.z << endl;*/
	fails++;
	bad++;
      }
      else goods++;
    }
    //point ma = truth_mom[v[0]], mb = truth_mom[v[1]];//v.size()-1]];
    //point m0 = start_mom[truth_part[v[0]]];
    //point dm = m0-ma*(sqrt(m0*m0)/sqrt(ma*ma));
    //if (sqrt(m0*m0) > 10)
    //  cerr << (sqrt(m0*m0)-sqrt(ma*ma)) << endl;
    //cout << sqrt(m0*m0)-sqrt(ma*ma) << ' ' << sqrt(ma*ma)-sqrt(mb*mb) << ' ' << sqrt(mb*mb) << endl;
    //cout << (sqrt(ma*ma)-sqrt(mb*mb)) << endl;///sqrt(ma*ma) << endl;
    /*if (bad) {
      static int cc = 0;
      if (cc++ == -5) {
	for (int hit_id : p.second) {
	  point p = truth_pos[hit_id];
	  //cerr << p.x << ' ' << p.y << endl;
	  cerr << p.z << ' ' << sqrt(p.x*p.x+p.y*p.y) << endl;
	}
	exit(0);
      }
      failparts++;
      cout << bad << ' ' << int(v.size())-2 << ' ' << part_q[p.first] << endl;
      }*/
  }
  /*
  cout << "fails: " << fails << endl;
  cout << "goods: " << goods << endl;
  cout << "failparts: " << failparts << endl;*/
}


//reordering of layers for approximate sorting
int topo[48], itopo[48];

void initOrder() {
  int handpicked[48] = {};
  int c = 0;
  for (int i = 0; i < 4; i++) handpicked[c++] = 7+i;
  for (int i = 0; i < 7; i++) handpicked[c++] = 7-1-i;
  for (int i = 0; i < 7; i++) handpicked[c++] = 11+i;

  for (int i = 0; i < 4; i++) handpicked[c++] = 24+i;
  for (int i = 0; i < 2; i++) handpicked[c++] = 40+i;
  for (int i = 0; i < 6; i++) {
    handpicked[c++] = 24-1-i;
    handpicked[c++] = 40-1-i;
  }

  for (int i = 0; i < 6; i++) {
    handpicked[c++] = 28+i;
    handpicked[c++] = 42+i;
  }
  for (int i = 0; i < 48; i++) {
    topo[i] = handpicked[i];
    itopo[topo[i]] = i;
  }
}



vector<point> hits; //hit position
vector<point> polar; //hit position in polar / cylindrical coordinates
vector<point> meta; //volume_id / layer_id / module_id
vector<int> metai, metaz; //ordered layer id in [0,48), and classification of z for disc layers in [0,4)
double disc_z[48][4];

const int Tube = 0, Disc = 1;

//Geometry of layer
struct Layer {
  double minr, avgr, maxr;
  double minz, avgz, maxz;
  int count;
  int type;

  double var0, var1;
};

Layer layer[48];
double z_minr[48][4], z_maxr[48][4];

//init layer geometries
void initLayers() {
  double avgz1[2][7] = {{-1500,-1300,-1100,-960,-820,-700,-600},
			{ 600, 700, 820, 960, 1100, 1300, 1500}};
  for (int k = 0; k < 2; k++)
    for (int i = 0; i < 7; i++) {
      layer[k*11+i].minr = 30;
      layer[k*11+i].maxr = 176.5;
      layer[k*11+i].avgz = avgz1[k][i];
      layer[k*11+i].type = Disc;
    }
  double avgz2[2][6] = {{-2950,-2550,-2150,-1800,-1500,-1220},
			{ 1220, 1500, 1800, 2150, 2550, 2950}};
  for (int k = 0; k < 2; k++)
    for (int i = 0; i < 6; i++) {
      layer[k*10+i+18].minr = 240;
      layer[k*10+i+18].maxr = 701;
      layer[k*10+i+18].avgz = avgz2[k][i];
      layer[k*10+i+18].type = Disc;

      layer[k*8+i+34].minr = 755;
      layer[k*8+i+34].maxr = 1018;
      layer[k*8+i+34].avgz = avgz2[k][i];
      layer[k*8+i+34].type = Disc;
    }

  double avgr1[4] = {32.3, 72.1, 116.1, 172.1};
  double avgr2[4] = {260.3, 360.2, 500.2, 660.2};
  double avgr3[2] = {820.2, 1020.2};

  for (int i = 0; i < 4; i++) {
    layer[i+7].minz =-491;
    layer[i+7].maxz = 491;
    layer[i+7].avgr = avgr1[i];
    layer[i+7].type = Tube;
  }
  for (int i = 0; i < 4; i++) {
    layer[i+24].minz =-1084;
    layer[i+24].maxz = 1084;
    layer[i+24].avgr = avgr2[i];
    layer[i+24].type = Tube;
  }
  for (int i = 0; i < 2; i++) {
    layer[i+40].minz =-1084;
    layer[i+40].maxz = 1084;
    layer[i+40].avgr = avgr3[i];
    layer[i+40].type = Tube;
  }
  Layer layer2[48];
  for (int i = 0; i < 48; i++) layer2[i] = layer[i];
  for (int i = 0; i < 48; i++) layer[i] = layer2[topo[i]];

  layer[0].var0 = 1e-3;
  layer[1].var0 = 5e-4;
  for (int i = 2; i < 18; i++) layer[i].var0 = 3e-4;
  for (int i = 18; i < 22; i++) layer[i].var0 = 5e-2;
  for (int i = 22; i < 48; i++) layer[i].var0 = i%2 || i == 22 ? 9 : 0.1;

  for (int i = 0; i < 4; i++) layer[i].var1 = 0.5;
  for (int i = 4; i < 18; i++) layer[i].var1 = 5;
  for (int i = 18; i < 24; i++) layer[i].var1 = 7;
  for (int i = 24; i < 48; i++) layer[i].var1 = i%2 ? 19 : 11;
}



void readHits() {
  initOrder();

  char file[1000];
  if (filenum >= 1000)
    sprintf(file, "%s/train_100_events/event%09d-hits.csv", base_path, filenum);
  else
    sprintf(file, "%s/test/event%09d-hits.csv", base_path, filenum);
  FILE*fp = fopen(file, "r");
  if (!fp) {
    printf("couldn't open hits\n");
    exit(1);
  }
  char tmpstr[1000];
  int tmp = fscanf(fp, "%s", tmpstr);

  // For one indexing
  hits.push_back(point(0,0,0));
  polar.push_back(point(0,0,0));
  meta.push_back(point(0,0,0));
  metai.push_back(0);

  int layers[9] = {7,4,7,6,4,6,6,2,6};
  int metai_list[9][7];
  int c = 0;
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < layers[i]; j++)
      metai_list[i][j] = c++;
  }
  cout << "Detectors: " << c << endl;

  for (int i = 0; i < 48; i++) {
    layer[i].minr = layer[i].minz = 1e9;
    layer[i].maxr = layer[i].maxz =-1e9;
  }

  while (1) {
    long long hit_id;
    double tx, ty, tz;
    int volume_id, layer_id, module_id;
    if (fscanf(fp, "%lld,%lf,%lf,%lf,%d,%d,%d", &hit_id, &tx, &ty, &tz, &volume_id, &layer_id, &module_id) == -1) break;
    if (hit_id != hits.size()) cout << "Hit id's not as expected" << endl;
    meta.push_back(point(volume_id, layer_id, module_id));
    if (volume_id <= 9)
      volume_id -= 7;
    else if (volume_id <= 14)
      volume_id -= 9;
    else
      volume_id -= 10;

    int mi = itopo[metai_list[volume_id][layer_id/2-1]];

    /*double a = rand()*(1./RAND_MAX)*M_PI*2;
    //double b = rand()*(1./RAND_MAX)*900-450;
    if (mi == 47) {
      double r = sqrt(tx*tx+ty*ty);
      tx = cos(a)*r;
      ty = sin(a)*r;
      //tz = b;
      }*/
    hits.push_back(point(tx, ty, tz));
    polar.push_back(point(sqrt(tx*tx+ty*ty), atan2(ty,tx), tz));
    metai.push_back(mi);

    double r = sqrt(tx*tx+ty*ty);
    Layer&l = layer[metai_list[volume_id][layer_id/2-1]];
    l.minr = min(l.minr, r);
    l.avgr += r;
    l.maxr = max(l.maxr, r);
    l.minz = min(l.minz, tz);
    l.avgz += tz;
    l.maxz = max(l.maxz, tz);
    l.count++;
    //cerr << tz << ' ' << r << endl;
  }
  fclose(fp);
  cout << hits.size() << " hits" << endl;

  initLayers();


  for (int hit_id = 1; hit_id < hits.size(); hit_id++) {
    metai_weight[truth_part[hit_id]][metai[hit_id]] += truth_weight[hit_id];
  }



  map<double, double> mir[48], mar[48];
  for (int i = 1; i < hits.size(); i++) {
    int mi = metai[i];
    if (layer[mi].type != Disc) continue;
    double&mir_ = mir[mi][polar[i].z];
    double&mar_ = mar[mi][polar[i].z];
    if (!mir_) mir_ = 1e9;
    mir_ = min(mir_, polar[i].x);
    mar_ = max(mar_, polar[i].x);
  }

  map<double, int> zi[48];
  for (int mi = 0; mi < 48; mi++) {
    if (layer[mi].type != Disc) continue;
    int k = 0;
    for (auto p : mir[mi]) {
      double mir_ = mir[mi][p.first]-1e-5;
      double mar_ = mar[mi][p.first]+1e-5;
      z_minr[mi][k] = mir_;
      z_maxr[mi][k] = mar_;
      disc_z[mi][k] = p.first;
      zi[mi][p.first] = k++;
    }
  }

  metaz.resize(hits.size());
  metaz[0] = 0;
  for (int i = 1; i < hits.size(); i++) {
    int mi = metai[i];
    metaz[i] = meta[i].z;
    if (layer[mi].type == Disc)
      metaz[i] = zi[mi][hits[i].z];
  }





  /*for (int i = 0; i < 48; i++) {
    layer[i].avgr /= layer[i].count;
    layer[i].avgz /= layer[i].count;
    cout << layer[i].minr << ' ' << layer[i].avgr << ' ' << layer[i].maxr << ' ' <<
      layer[i].minz << ' ' << layer[i].avgz << ' ' << layer[i].maxz << endl;
      }*/

  /*for (int i = 0; i < 48; i++) {
    for (int k = 0; k <= 100; k++)
    if (layer[i].type == Tube) {
      cerr << layer[i].minz << ' ' << layer[i].avgr << endl;
      cerr << layer[i].maxz << ' ' << layer[i].avgr << endl << endl;
    } else {
      cerr << layer[i].avgz << ' ' << layer[i].minr << endl;
      cerr << layer[i].avgz << ' ' << layer[i].maxr << endl << endl;
    }
  }
  exit(0);*/

  /*set<int> s;
  for (int i = 0; i < hits.size(); i++) {
    point p = hits[i];
    //cout << meta[i].x << endl;
    uint col = uint(meta[i].z)*0x10101*50;
    if (meta[i].x == 8 && meta[i].y == 2) {//>= 7 && meta[i].x <= 18)
      s.insert(meta[i].z);
      double r = p.z;//sqrt(p.x*p.x+p.y*p.y);
      cout << r << endl;
      cerr << p.x << ' ' << p.y << ' ' << p.z << ' ' << col << endl;
    }
    }*/
  //for (auto i : s) cout << i << endl;
  //exit(0);
}

// Volumes: 7,8,9, 12,13,14, 16,17,18
// Volume 8 is innermost, it has modules 2,4,6,8 concentric cylinders outwards


struct Detector {
  int volume_id, layer_id, module_id;
  point c;
  point rx, ry, rz;
  double d, minw, maxw, h, cell_w, cell_h;
};
map<int, Detector> detectors;

void readDetectors() {
  char file[1000];
  sprintf(file, "%s/detectors.csv", base_path);
  FILE*fp = fopen(file, "r");
  if (!fp) {
    printf("couldn't open detectors\n");
    exit(1);
  }

  char tmpstr[1000];
  int tmp = fscanf(fp, "%s", tmpstr);
  //cout << tmpstr << endl;

  Detector d;
  while (fscanf(fp, "%d,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &d.volume_id, &d.layer_id, &d.module_id, &d.c.x, &d.c.y, &d.c.z,
		// &d.rx.x, &d.rx.y, &d.rx.z,
		// &d.ry.x, &d.ry.y, &d.ry.z,
		// &d.rz.x, &d.rz.y, &d.rz.z,
		&d.rx.x, &d.ry.x, &d.rz.x,
		&d.rx.y, &d.ry.y, &d.rz.y,
		&d.rx.z, &d.ry.z, &d.rz.z,
		&d.d, &d.minw, &d.maxw, &d.h, &d.cell_w, &d.cell_h) == 21) {
    //if (d.module_id >= 10000 || d.layer_id >= 1000 || d.volume_id >= 100) cout << "What!?" << endl;
    detectors[d.volume_id*10000000+d.layer_id*10000+d.module_id] = d;
  }
  fclose(fp);
}

//pair<pair<ch0, ch1>, value>
vector<pair<pair<int, int>, double> > hit_cells[200000];
point hit_dir[200000][2]; //The two possible directions of the hit according to the cell's data for each hit

void readCells() {
  char file[1000];
  if (filenum >= 1000)
    sprintf(file, "%s/train_100_events/event%09d-cells.csv", base_path, filenum);
  else
    sprintf(file, "%s/test/event%09d-cells.csv", base_path, filenum);
  FILE*fp = fopen(file, "r");
  if (!fp) {
    printf("couldn't open cells\n");
    exit(1);
  }
  char tmpstr[1000];
  int tmp = fscanf(fp, "%s", tmpstr);
  //cout << tmpstr << endl;
  int hit_id, ch0, ch1;
  double value;
  while (fscanf(fp, "%d,%d,%d,%lf", &hit_id, &ch0, &ch1, &value) == 4) {
    hit_cells[hit_id].push_back(make_pair(make_pair(ch0, ch1), value));
    //cout << hit_id << ' ' << ch0 << ' ' << ch1 << ' ' << value << endl;
  }
  fclose(fp);
}



point normalize(point a) {
  point ret = a*(1./dist(a));
  if (ret.z < 0) ret = ret*-1;
  return ret;
}


//Calculate direction of each hit with cell's data
void initHitDir() {
  for (int hit_id = 1; hit_id < hits.size(); hit_id++) {
    point m = meta[hit_id];
    Detector&d = detectors[int(m.x)*10000000+int(m.y)*10000+int(m.z)];

    //if (!hit_cells[hit_id].size()) cout << "Hit with zero cells" << endl;
    //if (metai[hit_id] < 18) continue;

    //Use linear regression for direction
    double mx = 0, my = 0, mw = 0;
    auto&cells = hit_cells[hit_id];
    for (auto&c : cells) {
      double w = c.second;
      double x = c.first.first*d.cell_w;
      double y = c.first.second*d.cell_h;
      mw += w;
      mx += x*w;
      my += y*w;
    }
    mx /= mw;
    my /= mw;
    double mxx = 0, mxy = 0, myy = 0;
    for (auto&c : cells) {
      double w = c.second;
      double x = c.first.first*d.cell_w-mx;
      double y = c.first.second*d.cell_h-my;
      mxx += x*x*w;
      myy += y*y*w;
      mxy += x*y*w;
    }
    //Find eigenvector with minimum eigenvalue
    double a = mxx-myy, b = 2*mxy;
    double x = a+b+sqrt(a*a+b*b);
    double y =-a+b+sqrt(a*a+b*b);
    if (0) {
      double lambda = (mxx+myy+sqrt(a*a+b*b))/2;
      cout << lambda << ' ' << (mxx*x+mxy*y)/x << ' ' << (mxy*x+myy*y)/y << endl;
    }

    //Analytical formula for z
    double z = 2*d.d*(fabs(x)/d.cell_w+fabs(y)/d.cell_h+1e-8);
    x *= (cells.size()*1.-1.3);//1.3 != 1 was adjusted empirically
    y *= (cells.size()*1.-1.3);
    point d1(x,y,z), d2(x,y,-z);
    d1 = d.rx*d1.x+d.ry*d1.y+d.rz*d1.z;
    d2 = d.rx*d2.x+d.ry*d2.y+d.rz*d2.z;
    hit_dir[hit_id][0] = normalize(d1);
    hit_dir[hit_id][1] = normalize(d2);

    for (int k = 0; k < 2; k++)
      if (hit_dir[hit_id][k]*hits[hit_id] < 0)
	hit_dir[hit_id][k] = hit_dir[hit_id][k]*-1;

    //Write out missalignment to ground truth for plotting
    //if (metai[hit_id] == 0)
      //cerr << acos(max(fabs(hit_dir[hit_id][0]*truth_mom[hit_id]),
      //	       fabs(hit_dir[hit_id][1]*truth_mom[hit_id]))/dist(truth_mom[hit_id])) << endl;
  }
  //exit(0);
}
