
void storePairs(vector<pair<int, int> >&pairs) {
  FILE*fp = fopen("pairs", "w");
  fprintf(fp, "%d\n", int(pairs.size()));
  for (int i = 0; i < pairs.size(); i++)
    fprintf(fp, "%d %d\n", pairs[i].first, pairs[i].second);
  fclose(fp);
}
vector<pair<int, int> > loadPairs() {
  vector<pair<int, int> > pairs;
  FILE*fp = fopen("pairs", "r");
  int n;
  int tmp = fscanf(fp, "%d", &n);
  pairs.resize(n);
  for (int i = 0; i < pairs.size(); i++)
    tmp = fscanf(fp, "%d%d", &pairs[i].first, &pairs[i].second);
  fclose(fp);
  return pairs;
}



void storeTriples(vector<triple>&triples) {
  FILE*fp = fopen("triples", "w");
  fprintf(fp, "%d\n", int(triples.size()));
  for (int i = 0; i < triples.size(); i++)
    fprintf(fp, "%d %d %d\n", triples[i].x, triples[i].y, triples[i].z);
  fclose(fp);
}
vector<triple> loadTriples() {
  vector<triple> triples;
  FILE*fp = fopen("triples", "r");
  int n;
  int tmp = fscanf(fp, "%d", &n);
  triples.resize(n);
  for (int i = 0; i < triples.size(); i++)
    tmp = fscanf(fp, "%d%d%d", &triples[i].x, &triples[i].y, &triples[i].z);
  fclose(fp);
  return triples;
}



void storePaths(vector<vector<int> >&paths, int id) {
  char filename[100];
  sprintf(filename, "paths%d", id);
  FILE*fp = fopen(filename, "w");
  fprintf(fp, "%lu\n", paths.size());
  for (int i = 0; i < paths.size(); i++) {
    fprintf(fp, "%lu", paths[i].size());
    for (int j = 0; j < paths[i].size(); j++)
      fprintf(fp, " %d", paths[i][j]);
    fprintf(fp, "\n");
  }
  fclose(fp);
}

vector<vector<int> > loadPaths(int id) {
  vector<vector<int> > paths;
  char filename[100];
  sprintf(filename, "paths%d", id);
  FILE*fp = fopen(filename, "r");
  int n;
  int tmp;
  tmp = fscanf(fp, "%d", &n);
  paths.resize(n);
  for (int i = 0; i < paths.size(); i++) {
    int k;
    tmp = fscanf(fp, "%d", &k);
    paths[i].resize(k);
    for (int j = 0; j < k; j++)
      tmp = fscanf(fp, "%d", &paths[i][j]);
  }
  fclose(fp);
  return paths;
}


void writeSubmission(map<int, int>&assignment) {
  char filename[1000];
  sprintf(filename, "submissions/submission%d.csv", filenum);
  FILE*fp = fopen(filename, "w");
  if (!fp) {
    cout << "Could not open " << filename << " for writing submission" << endl;
    return;
  }
  fprintf(fp, "event_id,hit_id,track_id\n");

  for (int i = 1; i < hits.size(); i++) {
    int a = 0;
    if (assignment.count(i)) a = assignment[i];
    fprintf(fp, "%d,%d,%d\n", filenum, i, a);
  }
  fclose(fp);
}

map<int, int> loadAssignment() {
  char filename[1000];
  sprintf(filename, "submissions/submission%d.csv", filenum);
  FILE*fp = fopen(filename, "r");
  map<int, int> ret;
  if (!fp) {
    cout << "Could not open " << filename << " for reading assignment" << endl;
    return ret;
  }
  int tmp = fscanf(fp, "%s", filename);

  for (int i = 1; i < hits.size(); i++) {
    int a, b, c;
    int tmp = fscanf(fp, "%d,%d,%d", &a, &b, &c);
    if (c) ret[i] = c;
    //assignment[i] = c;
  }
  fclose(fp);
}
