// Contains functions relating to outlier densities and extension of helices

#define fori(a) for (int i = 0; i < a; i++)
#define forij(a, b) for (int i = 0; i < a; i++) for (int j = 0; j < b; j++)

//Decides how curved "approximately straight" helices are supposed to be
const double stretch = 0.02;

//List of coordinates for each layer, used to estimate outlier densities
vector<double> sorted_hits[48];

//Acceleration look-up table for faster lookup in sorted_hits
const int crude_steps = 1<<10;
int crudeIndex[48][crude_steps];
//Scaling to get sorted_hits in range [0,crude_steps)
pair<double, double> crudeIndex_a[48];
//Accumulated polynomial coefficients of second order polynomial for O(1) density look-up. Second dimension (20000) must be bigger than number of hits in the most populated layer
double poly[48][20000][3];


// O(1) indexing in sorted_hits
// Functionally similar to "upper_bound(sorted_hits[li].begin(), sorted_hits[li].end(), x)-sorted_hits[li].begin();"
inline int getIndex(int&li, double x) {
  int ci = x*crudeIndex_a[li].first+crudeIndex_a[li].second;
  ci = min(crude_steps-1, max(0, ci));
  int i = crudeIndex[li][ci];
  //cout << x << ' ' << sorted_hits[li][i] << ' ' << crudeIndex_a[li].first << endl;
  //if (i < 1) i = 1;
  //if (i > sorted_hits[li].size()-2) i = sorted_hits[li].size()-2;
  //static int c = 0, d = 0;
  //while (i+1 < sorted_hits[li].size() && x >= sorted_hits[li][i]) i++, c++;
  //while (i && x < sorted_hits[li][i-1]) i--, c++;

  //Might segfault sometimes :)
  while (x >= sorted_hits[li][i]) i++;
  while (x < sorted_hits[li][i-1]) i--;

  //d++;
  //if (d%1000000 == 0) cout << c*1./d << ' ' << c << ' ' << d << endl;
  /*{
    int j = upper_bound(sorted_hits[li].begin(), sorted_hits[li].end(), x)-sorted_hits[li].begin();
    if (j != i) {
      cout << i << ' ' << j << endl;
      cout << x << endl;
      cout << sorted_hits[li][i-1] << ' ' << sorted_hits[li][i] << ' ' << sorted_hits[li][i+1] << endl;
      cout << sorted_hits[li][j-1] << ' ' << sorted_hits[li][j] << ' ' << sorted_hits[li][j+1] << endl;
      cout << endl;

    }
    }*/
  return i;//max(0,min(i,int(sorted_hits[li].size())-1));
}

//Init everything needed for fast density calculations, includes most global variables above
void initDensity3() {
  vector<int>*tube = new vector<int>[48]();

  for (int i = 1; i < hits.size(); i++) {
    if (!assignment[i])
      tube[metai[i]].push_back(i);
  }

  for (int li = 0; li < 48; li++) {
    sorted_hits[li].clear();
    if (layer[li].type == Tube) {
      for (int i : tube[li])
	sorted_hits[li].push_back(hits[i].z);
    } else {
      for (int i : tube[li])
	sorted_hits[li].push_back(polar[i].x);
    }
    sorted_hits[li].push_back(-1e50);
    sorted_hits[li].push_back(1e50);
    sort(sorted_hits[li].begin(), sorted_hits[li].end());

    double minx = *next(sorted_hits[li].begin())-1e-8;
    double maxx = *next(sorted_hits[li].rbegin())+1e-8;
    //cout << maxx << ' ' << minx << endl;
    double f = crude_steps/(maxx-minx);
    crudeIndex_a[li] = make_pair(f, -minx*f);

    for (int i = 0; i < crude_steps; i++) {
      double x = (i+.5)/f+minx;
      crudeIndex[li][i] = upper_bound(sorted_hits[li].begin(), sorted_hits[li].end(), x)-sorted_hits[li].begin();
    }
    double acc[3] = {};
    for (int i = 1; i < sorted_hits[li].size(); i++) {
      for (int j = 0; j < 3; j++) poly[li][i][j] = acc[j];
      double x = sorted_hits[li][i];
      for (int j = 0; j < 3; j++)
	acc[j] += pow(x, j);
    }
  }
  delete[]tube;
}


//Get expected number of hits on layer "li" in area (in polar/cylindrical coordnates) spanned by (p-dp)^2-dot(p-dp, xp)^2 < tt
double getDensity3(point&dp, point&xp, double tt, int li) {
  Layer&l = layer[li];

  double x0, dx, dy;
  if (l.type == Tube) {
    x0 = dp.z;
    dx = xp.z;
    dy = xp.y;
  } else {
    x0 = dp.x;
    dx = xp.x;
    dy = xp.y;
  }
  double b = tt*(1-dy*dy)/(1-dx*dx-dy*dy);
  double a = sqrt(1-dx*dx-dy*dy)/((1-dy*dy)*M_PI*dp.x);
  double rx = sqrt(b);

  /*int ai = getIndex(li, x0-rx);
  int bi = getIndex(li, x0);
  int ci = getIndex(li, x0+rx);

  double ret =
    ((poly[li][bi][0]-poly[li][ai][0])*(rx-x0)+
     (poly[li][bi][1]-poly[li][ai][1])*(1)+
     (poly[li][ci][0]-poly[li][bi][0])*(rx+x0)+
     (poly[li][ci][1]-poly[li][bi][1])*(-1))*a;
  if (ret < 0) {
    cout << x0-rx << ' ' << x0 << ' ' << x0+rx << endl;
    cout << sorted_hits[li][ai] << ' ' << sorted_hits[li][bi] << ' ' << sorted_hits[li][ci] << endl;
    cout << ret << ' ' << getDensity2(dp, xp, tt, li) << endl;
    cout << ai << ' ' << bi << ' ' << ci << endl;
    cout << poly[li][ai][0] << ' ' << poly[li][ai][1] << endl;
    cout << poly[li][bi][0] << ' ' << poly[li][bi][1] << endl;
    cout << poly[li][ci][0] << ' ' << poly[li][ci][1] << endl;
    cout << endl;
    ret = 1e-8;
  }
  return ret;*/

  int ai = getIndex(li, x0-rx);
  int bi = getIndex(li, x0+rx);
  if (bi-ai > 10) {//Approximate integration by 2. order polynomial approximation to half disc
    //cout << ai << ' ' << bi << endl;
    const double A = 21*M_PI/64., B = -15*M_PI/64.;
    double ib = 1./b;
    double c0 = A+B*x0*x0*ib, c1 = -2*B*ib*x0, c2 = B*ib;
    double ret =
      ((poly[li][bi][0]-poly[li][ai][0])*c0+
       (poly[li][bi][1]-poly[li][ai][1])*c1+
       (poly[li][bi][2]-poly[li][ai][2])*c2)*a*rx;
    return max(ret,0.);
  } else { //Exact integration, uses half disc
    double density = 0;
    for(int i = ai; i < bi; i++) {
      double x = sorted_hits[li][i]-x0;
      double h = a*sqrt(b-x*x);/// *it;
      //cout << h << endl;
      density += h;
    }
    return density;
  }
}

//Find density by binary search
//This means we want to find (and return) tt such that getDensity3(dp, xp, tt, li) = target
double findDensity(point&dp, point&xp, double target, int li) {
  double Ad = 0, A = 0, B = 1, Bd;
  while (1) {
    Bd = getDensity3(dp, xp, B, li);
    //cout << B << ' ' << Bd << endl;
    if (B > 1e20) {
      cout << "No density?" << endl;
      cout << dp << ' ' << xp << ' ' << li << endl;
      exit(0);
    }
    if (Bd > target) break;
    B *= 10;
    if (target/Bd < 1e8) B = max(B, target/Bd);
  }
  double mid = B/2;
  int cc = 0;
  while (1) {
    double density = getDensity3(dp, xp, mid, li);
    if (density > target) B = mid, Bd = density;
    else A = mid, Ad = density;

    //cout << A << ' ' << mid << ' ' << B << ' ' << density << endl;
    if ((B-A) < A*1e-3 || density > target*0.9 && density < target*1.1 || cc >= 100) break;
    mid = max(A*0.9+B*0.1, min(B*0.9+A*0.1, (target-Ad)*(B-A)/(Bd-Ad)+A));
    if (++cc == 100) { //Should never happen
      cout << "Warning: Infinite loop in findDensity" << endl;
      /*cout << dp << endl;
      cout << xp << endl;
      cout << mid << endl;
      cout << target << endl;
      cout << li << endl;
      exit(0);*/
    }
  }
  return mid;
}




//Prepare ellipse equation of collision between line extrapolated through hits with id "ai" and "bi" and layer "li". Return collision coordinate "d", in polar coordinates "dp", ellipse stretching "xp", and direction of hit in polar coordnates "bap". "target" describes the layer, possibly corrected for a single point we are evaluating a helix quadruple
int prepareTripleScore(int ai, int bi, int li, point&d, point&dp, point&xp, point&bap, point target) {
  const double slack = 1.00; //No slack

  Layer&l = layer[li];
  point&a = hits[ai], &b = hits[bi];

  point ba = b-a;
  if (l.type == Tube) {
    double vv = ba.x*ba.x+ba.y*ba.y;
    double pv = ba.x*a.x+ba.y*a.y;
    double pp = a.x*a.x+a.y*a.y;
    double RR = target.x*target.x;
    double sq = pv*pv-vv*(pp-RR);
    if (sq < 0) return -1;

    double t = (-pv+sqrt(sq))/vv;
    if (t < 0 || !vv) return -1;
    d.x = a.x+ba.x*t;
    d.y = a.y+ba.y*t;
    d.z = a.z+ba.z*t;

    if (d.z < l.minz*slack || d.z > l.maxz*slack) return -1;

    dp = point(dist(d.x,d.y),atan2(d.y,d.x),d.z);

    xp = point(0, -dp.x*(ba.x*ba.x+ba.y*ba.y), ba.z);
    bap = point(ba.x*d.x+ba.y*d.y, d.x*ba.y-d.y*ba.x, ba.z*dp.x);

    bap = bap*(1./bap.x);
  } else if (l.type == Disc) {
    double t = (target.z-a.z)/ba.z;
    if (t < 0 || !ba.z) return -1;
    d.x = a.x+ba.x*t;
    d.y = a.y+ba.y*t;
    d.z = a.z+ba.z*t;

    dp = point(dist(d.x,d.y),atan2(d.y,d.x),d.z);

    if (dp.x < l.minr*(1./slack) || dp.x > l.maxr*slack) return -1;

    xp = point(ba.x*d.y-ba.y*d.x, d.x*ba.x+d.y*ba.y, 0);
    bap = point(ba.x*d.x+ba.y*d.y, d.x*ba.y-d.y*ba.x, ba.z*dp.x);

    bap = bap*(1./bap.z);
  }
  double xp2 = xp.x*xp.x+xp.y*xp.y+xp.z*xp.z;
  if (xp2)
    xp = xp*sqrt((1-stretch)/xp2);
  return 0;
}

//Default is using average position in the layer
int prepareTripleScore(int ai, int bi, int li, point&d, point&dp, point&xp, point&bap) {
  Layer&l = layer[li];
  point target(l.avgr, 0, l.avgz);
  return prepareTripleScore(ai, bi, li, d, dp, xp, bap, target);
}

//Use the prepared "dp", "xp", "bap" and return the area that is closer to the collision line (taking into account xp for elliptic behaviour) compared to the hit with id "ci"
double evaluateScore(int ci, point&dp, point&xp, point&bap) {
  point&r = polar[ci];
  point err = r-dp;
  if (err.y > M_PI) err.y -= M_PI*2;
  if (err.y <-M_PI) err.y += M_PI*2;
  err.y *= dp.x;

  err = err-bap*(layer[metai[ci]].type == Disc ? err.z : err.x);
  double r2 = err*err-pow(err*xp, 2);
  return r2;
}



//Approximate magnetic field strengh as a function of z coordinate, decays drastically near the ends
double field(double z) {
  z *= 1./2750;
  double z2 = z*z;
  return 1.002-z*3e-2-z2*(0.55-0.3*(1-z2));
}



//Similar to prepareTripleScore, but now we extend the helix that passes through hits with id "ai", "bi", "ci". Assumes li > metai[ci] if sign = 1, and metai[bi] < li < metai[ci] if sign = -1
int prepareQuadrupleScore(int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&dirp, point target, double sign = 1) {
  Layer&l = layer[li];

  point p;
  double r, ir;

  point c = hits[ci];
  point cb = hits[ci]-hits[bi];//c-b;
  double ang_cb;

  if (0) {
    point a = hits[ai], b = hits[bi], c = hits[ci];

    //TODO: test if has bieffects
    if (l.type == Disc && (c.z-b.z)*(target.z-c.z) < 0) return -1;

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

    ang_cb = asin(dist(cb.x, cb.y)*.5*ir)*2;
  }

  if (1) { //Take into account varying magnetic field strength
    point a = hits[ai], b = hits[bi], c = hits[ci];
    if (l.type == Disc && (c.z-b.z)*(target.z-c.z)*sign < 0) return -1;
    double B1 = field((a.z+b.z)*.5), B2 = field((b.z+c.z)*.5), B3;
    if (l.type == Disc) B3 = field((c.z+target.z)*.5);
    else B3 = field(c.z);
    //B1 = B2 = B3 = 1;
    double ax = b.x-a.x, ay = b.y-a.y, bx = c.x-b.x, by = c.y-b.y;
    double aa = ax*ax+ay*ay, dot = ax*bx+ay*by, cross = ax*by-ay*bx;
    double alpha = B2/(2*B3), beta = (-B1*aa-B2*dot)/(2*cross*B3);
    //alpha *= -1;
    double rx = alpha*bx-beta*by, ry = alpha*by+beta*bx;
    p = point(c.x-rx, c.y-ry, 0);
    r = dist(rx, ry);
    ir = 1./r;
    /*cout << endl;
    //cout << (rx*ax+ry*ay)/(ax*ax+ay*ay) << endl;
    //cout << (rx*bx+ry*by)/(bx*bx+by*by) << endl;

    point p(c.x-rx, c.y-ry, 0);
    cout << dist(a.x-p.x, a.y-p.y) << ' ' << dist(b.x-p.x, b.y-p.y) << ' ' << dist(c.x-p.x, c.y-p.y) << endl;
    cout << r << ' ' << dist(rx, ry) << endl;*/
    ang_cb = B3/B2*asin(dist(cb.x, cb.y)*.5*ir*B2/B3)*2;
  }

  //exit(0);

  const double slack = 1.00;

  xp = point(0,0,0); //Circle for now

  if (l.type == Tube) {
    double RR = target.x*target.x, pp = dist2(p.x, p.y);
    double s = .5+(RR-r*r)/(2*pp);
    double sq = RR/pp-s*s;

    /*
    if (truth_part[ai] == 4506417142702082LL) {
      cout << sq << endl;
      cout << metai[ai] << ' ' << metai[bi] << ' ' << metai[ci] << endl;
    }
    */

    if (sq < 0) return -1;

    double t = sqrt(sq);
    if (p.y*c.x-p.x*c.y < 0) t *= -1;
    d.x = p.x*s+p.y*t;
    d.y = p.y*s-p.x*t;

    point dc = d-c;
    double A = dist(dc.x, dc.y);
    double B = A*.5*ir;
    double ang_dc = asin(B)*2;
    if (dc.x*cb.x+dc.y*cb.y < 0) ang_cb *= -1;

    d.z = c.z+cb.z*ang_dc/ang_cb;

    if (!(d.z > l.minz*slack && d.z < l.maxz*slack)) return -1;

    point dir;
    double s_ = target.x/pp, t_ = s_*(1-s)/t;
    dir.x = p.x*s_+p.y*t_;
    dir.y = p.y*s_-p.x*t_;
    dir.z = (dc.x*dir.x+dc.y*dir.y)*ir*cb.z/(ang_cb*A*sqrt(1-B*B));

    dp = point(dist(d.x,d.y), atan2(d.y, d.x), d.z);
    dirp = point(d.x*dir.x+d.y*dir.y, d.x*dir.y-d.y*dir.x, dir.z*dp.x);
    //cout << dirp << endl; //dirp.x = l.avgr

    dirp = dirp*(1./dirp.x);

  } else if (l.type == Disc) {
    d.z = target.z;
    double fac = ang_cb/cb.z;
    double ang_dc = (d.z-c.z)*fac;

    double sa = sin(ang_dc), ca = cos(ang_dc);

    double rx = c.x-p.x, ry = c.y-p.y;
    double cross = rx*cb.y-ry*cb.x;
    if (cross < 0) sa *= -1;

    d.x = ca*rx-sa*ry+p.x;
    d.y = sa*rx+ca*ry+p.y;


    point dir;
    dir.x =-fac*(rx*sa+ry*ca);
    dir.y = fac*(rx*ca-ry*sa);
    dir.z = cross < 0 ? -1 : 1;


    dp = point(dist(d.x,d.y), atan2(d.y, d.x), d.z);

    if (!(dp.x > l.minr*(1./slack) && dp.x < l.maxr*slack)) return -1;

    dirp = point(d.x*dir.x+d.y*dir.y, d.x*dir.y-d.y*dir.x, dir.z*dp.x);
    //cout << dirp << endl; //dirp.x = l.avgr

    dirp = dirp*(1./dirp.z);
  }
  return 0;
}


//Not used
int prepareQuadrupleScore2(int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&dirp, point target) {
  Layer&l = layer[li];

  point p;
  double r, ir;

  point c = hits[ci];
  point cb = hits[ci]-hits[bi];//c-b;
  double ang_cb;

  if (0) {
    point a = hits[ai], b = hits[bi], c = hits[ci];

    //TODO: test if has bieffects
    if (l.type == Disc && (c.z-b.z)*(target.z-c.z) > 0) return -1;

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

    ang_cb = asin(dist(cb.x, cb.y)*.5*ir)*2;
  }

  if (1) {
    point a = hits[ai], b = hits[bi], c = hits[ci];
    if (l.type == Disc && (c.z-b.z)*(target.z-c.z) > 0) return -1;
    double B1 = field((a.z+b.z)*.5), B2 = field((b.z+c.z)*.5), B3;
    if (l.type == Disc) B3 = field((c.z+target.z)*.5);
    else B3 = field(c.z);
    //B1 = B2 = B3 = 1;
    double ax = b.x-a.x, ay = b.y-a.y, bx = c.x-b.x, by = c.y-b.y;
    double aa = ax*ax+ay*ay, dot = ax*bx+ay*by, cross = ax*by-ay*bx;
    double alpha = B2/(2*B3), beta = (-B1*aa-B2*dot)/(2*cross*B3);
    //alpha *= -1;
    double rx = alpha*bx-beta*by, ry = alpha*by+beta*bx;
    p = point(c.x-rx, c.y-ry, 0);
    r = dist(rx, ry);
    ir = 1./r;
    /*cout << endl;
    //cout << (rx*ax+ry*ay)/(ax*ax+ay*ay) << endl;
    //cout << (rx*bx+ry*by)/(bx*bx+by*by) << endl;

    point p(c.x-rx, c.y-ry, 0);
    cout << dist(a.x-p.x, a.y-p.y) << ' ' << dist(b.x-p.x, b.y-p.y) << ' ' << dist(c.x-p.x, c.y-p.y) << endl;
    cout << r << ' ' << dist(rx, ry) << endl;*/
    ang_cb = B3/B2*asin(dist(cb.x, cb.y)*.5*ir*B2/B3)*2;
  }

  //exit(0);

  const double slack = 1.00;

  xp = point(0,0,0); //Circle for now

  if (l.type == Tube) {
    double RR = target.x*target.x, pp = dist2(p.x, p.y);
    double s = .5+(RR-r*r)/(2*pp);
    double sq = RR/pp-s*s;

    /*
    if (truth_part[ai] == 4506417142702082LL) {
      cout << sq << endl;
      cout << metai[ai] << ' ' << metai[bi] << ' ' << metai[ci] << endl;
    }
    */

    if (sq < 0) return -1;

    double t = sqrt(sq);
    if (p.y*c.x-p.x*c.y < 0) t *= -1;
    d.x = p.x*s+p.y*t;
    d.y = p.y*s-p.x*t;

    point dc = d-c;
    double A = dist(dc.x, dc.y);
    double B = A*.5*ir;
    double ang_dc = asin(B)*2;
    if (dc.x*cb.x+dc.y*cb.y < 0) ang_cb *= -1;
    d.z = c.z+cb.z*ang_dc/ang_cb;

    if (!(d.z > l.minz*slack && d.z < l.maxz*slack)) return -1;

    point dir;
    double s_ = target.x/pp, t_ = s_*(1-s)/t;
    dir.x = p.x*s_+p.y*t_;
    dir.y = p.y*s_-p.x*t_;
    dir.z = (dc.x*dir.x+dc.y*dir.y)*ir*cb.z/(ang_cb*A*sqrt(1-B*B));

    dp = point(dist(d.x,d.y), atan2(d.y, d.x), d.z);
    dirp = point(d.x*dir.x+d.y*dir.y, d.x*dir.y-d.y*dir.x, dir.z*dp.x);
    //cout << dirp << endl; //dirp.x = l.avgr

    dirp = dirp*(1./dirp.x);

  } else if (l.type == Disc) {
    d.z = target.z;
    double fac = ang_cb/cb.z;
    double ang_dc = (d.z-c.z)*fac;

    double sa = sin(ang_dc), ca = cos(ang_dc);

    double rx = c.x-p.x, ry = c.y-p.y;
    double cross = rx*cb.y-ry*cb.x;
    if (cross < 0) sa *= -1;

    d.x = ca*rx-sa*ry+p.x;
    d.y = sa*rx+ca*ry+p.y;


    point dir;
    dir.x =-fac*(rx*sa+ry*ca);
    dir.y = fac*(rx*ca-ry*sa);
    dir.z = cross < 0 ? -1 : 1;


    dp = point(dist(d.x,d.y), atan2(d.y, d.x), d.z);

    if (!(dp.x > l.minr*(1./slack) && dp.x < l.maxr*slack)) return -1;

    dirp = point(d.x*dir.x+d.y*dir.y, d.x*dir.y-d.y*dir.x, dir.z*dp.x);
    //cout << dirp << endl; //dirp.x = l.avgr

    dirp = dirp*(1./dirp.z);
  }
  return 0;
}


//Default target is average position of layer
int prepareQuadrupleScore(int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&bap, double sign = 1) {
  Layer&l = layer[li];
  point target(l.avgr, 0, l.avgz);
  return prepareQuadrupleScore(ai, bi, ci, li, d, dp, xp, bap, target, sign);
}



//Return this if no points were found, somewhat tunable parameter
const double density_eps = 1e-6;

//map<pair<pair<int, int>, int >, double> scoreTripleDensity_mem;

//How many outliers do we expect to fit better than "ci" in the triple "ai", "bi", "ci"?
double scoreTripleDensity(int ai, int bi, int ci) {
  /*auto memi = make_pair(make_pair(ai, bi), ci);
  double&memo = scoreTripleDensity_mem[memi];
  if (!memo) {
  */
    point d, dp, xp, bap;
    if (prepareTripleScore(ai, bi, metai[ci], d, dp, xp, bap, polar[ci])) return 1e9;
    double s = evaluateScore(ci, dp, xp, bap);
    s = getDensity3(dp, xp, s, metai[ci]);
    return s+density_eps;
    /*
    memo = s+density_eps;
  }
  return memo;*/
}


//map<pair<pair<int, int>, pair<int, int> >, double> scoreQuadrupleDensity_mem;

//How many outliers do we expect to fit better than "di" in the triple "ai", "bi", "ci", "di"?
double scoreQuadrupleDensity(int ai, int bi, int ci, int di) {
  /*auto memi = make_pair(make_pair(ai, bi), make_pair(ci, di));
  double&memo = scoreQuadrupleDensity_mem[memi];
  if (!memo) {*/
    //cout << metai[ai] << ' ' << metai[bi] << ' ' << metai[ci] << ' ' << metai[di] << endl;
    point d, dp, xp, bap;
    if (prepareQuadrupleScore(ai, bi, ci, metai[di], d, dp, xp, bap, polar[di])) return 1e9;
    double s = evaluateScore(di, dp, xp, bap);
    //cout << "S0: " << s << endl;
    s = getDensity3(dp, xp, s, metai[di]);
    //if (!s) cout << "What?" << endl;
    return s+density_eps;
    /*
    memo = s+density_eps;
  }
  return memo;*/
}


//Similar to the other prepareXScore functions, but now we try to find duplicates on layer "li" to one of the input hits. This means looking for hits that are close to the helix fitted through "ai", "bi", "ci"
int prepareDuplicateScore(int ai, int bi, int ci, int li, point&d, point&dp, point&xp, point&dirp) {
  Layer&l = layer[li];
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

  int di = -1;
  if (metai[ai] == li) di = ai;
  else if (metai[bi] == li) di = bi;
  else if (metai[ci] == li) di = ci;
  else {
    cout << "prepareDuplocateScore given layeri not corresponding to ai, bi or ci" << endl;
    return -1;
  }
  d = hits[di];
  dp = polar[di];

  double rx = hits[di].x-p.x, ry = hits[di].y-p.y;

  //TODO: do with respect to nearest circle arc, not ca
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

  /*
  //dir = truth_mom[di];
  {
    point dir2 = dir*(-1./sqrt(dir*dir));
    point mom2 = truth_mom[di]*(1./dist(truth_mom[di]));
    if (dir2.z * mom2.z < 0) dir2 = dir2*-1;
    cout << dir2 << endl;
    cout << mom2 << endl << endl;
  }
  */

  xp = point(0,0,0);
  dirp = point(d.x*dir.x+d.y*dir.y, d.x*dir.y-d.y*dir.x, dir.z*dp.x);
  //cout << dirp << endl; //dirp.x = l.avgr
  if (l.type == Tube)
    dirp = dirp.x ? dirp*(1./dirp.x) : point(1,0,0);
  else
    dirp = dirp.z ? dirp*(1./dirp.z) : point(0,0,1);

  return 0;
}

//How many outliers do we expect to fit better than "di" to the triple "ai", "bi", "ci"
double scoreDuplicateDensity(int ai, int bi, int ci, int di) {
  //cout << metai[ai] << ' ' << metai[bi] << ' ' << metai[ci] << ' ' << metai[di] << endl;
  point d, dp, xp, bap;
  if (prepareDuplicateScore(ai, bi, ci, metai[di], d, dp, xp, bap)) return 1e9;
  double s = evaluateScore(di, dp, xp, bap);
  //cout << "S0: " << s << endl;
  s = getDensity3(dp, xp, s, metai[di]);
  //if (!s) cout << "What?" << endl;
  return s+density_eps;
}





// Score a full path by looking at the probability that it could happen by outliers alone
// Multiply average number of outliers gotten from
// - The first triple
// - All consecutive quadruples
// - Duplicates
// This function is very important for score, and full of tuning opportunities
double scorepathDensity(vector<int>&path) {
  static vector<int> u(48);
  u.resize(0);
  //cout << u.capacity() << endl;

  int last_metai = -1;
  for (int i : path) {
    if (i <= 0) continue;
    int mi = metai[i];
    if (mi != last_metai) {
      //if (last_metai != -1 && next_layer[last_metai][mi] < adj_thres) return -1e4;
      u.push_back(i);
      last_metai = mi;
    }
  }
  if (u.size() < 3) return -0.1;
  double prod = 1;
  double s = scoreTripleDensity(u[0], u[1], u[2]);//+1e-4;
  s = min(s, 1e4);
  prod *= s;



  //Hugely important tuning parameters
  double quad_off = 4;
  double dup_off = 10;
  double quad_max = 50;
  double quad_max1 = 10;
  double dup_max = 1;
  //if (u.size() >= 4)
  //  prod *= pow(scoreQuadrupleDensity(u[3], u[2], u[1], u[0]), 1);
  for (int i = 3; i < u.size(); i++) {
    double s = min(scoreQuadrupleDensity(u[i-3], u[i-2], u[i-1], u[i]), quad_max1);
    s *= min(scoreQuadrupleDensity(u[i], u[i-1], u[i-2], u[i-3]), quad_max1)*quad_off;
    //s = min(s, quad_max);
    prod *= s;
  }

  int j = 0;
  for (int i = 0; i < path.size(); i++) {
    if (path[i] <= 0) continue;
    while (metai[path[i]] != metai[u[j]]) j++;
    if (path[i] == u[j]) continue;
    int k = max(j, 2);
    double s = pow(scoreDuplicateDensity(u[k-2], u[k-1], u[k], path[i]), 0.92)*dup_off; // Really important magic constant of 0.92, I have no clue why
    s = min(s, dup_max);
    prod *= s;
  }

  return -prod;
}






//Score triple based on the deviation from a perfect helix, no prior that it should be straight
double scoreTriple(int ai, int bi, int ci) {
  point center;
  double radius;
  circle(hits[ai], hits[bi], hits[ci], center, radius);

  point cb = hits[ci]-hits[bi];
  point ba = hits[bi]-hits[ai];
  double ang_cb = asin(dist(cb.x, cb.y)*.5/radius)*2;
  double ang_ba = asin(dist(ba.x, ba.y)*.5/radius)*2;
  if (radius != radius || fabs(radius) > 1e50) {
    ang_cb = dist(cb.x, cb.y);
    ang_ba = dist(ba.x, ba.y);
  }
  if (ba.z*cb.z < 0) ang_ba *= -1;

  //if (dist(cb.x, cb.y)*.5/radius > M_PI/2 || dist(ba.x, ba.y)*.5/radius > M_PI/2) return 1e9;
  //-radius*2e-5+
  double x = ba.z ? (fabs(cb.z*ang_ba/ba.z-ang_cb))*radius : 1e9;
  double y = ang_cb ? (fabs(cb.z*ang_ba/ang_cb-ba.z)) : 1e9;
  double score = min(x, y);//, fabs(cb.z-ba.z*ang_cb/ang_ba)));
  /*
  cout << endl;
  cout << truth_mom[bi]*part_q[truth_part[bi]] << endl;
  point rr = hits[bi]-center;
  if (cb.x*rr.y-cb.y*rr.x > 0) ang_cb *= -1;
  cout << point(-rr.y, rr.x, cb.z/ang_cb) << endl;*/
  return score;
}
