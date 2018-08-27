//Contains functions to estimate score at different steps in the pipeline (assuming the remaining steps are done optimally)

void scorePairs(vector<pair<int, int> >&pairs) {
  set<long long> found;
  for (auto p : pairs) {
    load(pairs.size(), "Pair score:  ", 1);
    if (samepart(p.first, p.second)) {
      found.insert(truth_part[p.first]);
    }
  }
  double score = 0;
  for (long long p : found) {
    score += part_weight[p];
  }
  cout << score << " from " << pairs.size() << " pairs" << endl;
}

void scorePairs2(vector<pair<int, int> >&pairs) {
  set<long long> found;
  for (auto p : pairs) {
    load(pairs.size(), "Pair score 2:", 1);
    if (samepart(p.first, p.second)) {
      int prev1 = -1, prev2 = -1, ok = 0;
      for (int i : truth_tracks[truth_part[p.first]]) {
	int m = metai[i];
	if (m != prev1) {
	  prev2 = prev1;
	  prev1 = m;
	}
	if (i == p.second) {
	  if (prev2 == metai[p.first]) ok = 1;
	  break;
	}
      }
      if (ok)
	found.insert(truth_part[p.first]);
    }
  }
  double score = 0;
  for (long long p : found) {
    score += part_weight[p];
  }
  cout << score << " from " << pairs.size() << " pairs" << endl;
}



void scoreTriples(vector<triple>&triples) {
  set<long long> found;
  for (auto p : triples) {
    load(triples.size(), "Triple score:", 1);
    if (samepart(p.x, p.y) && samepart(p.x, p.z)) {
      found.insert(truth_part[p.x]);
    }
  }
  double score = 0;
  for (long long p : found) {
    score += part_weight[p];
  }
  cout << score << " from " << triples.size() << " triples" << endl;
}

void scoreOriginTriples(vector<triple>&triples) {
  set<long long> found;
  for (auto p : triples) {
    load(triples.size(), "Origin triple score:", 1);
    if (samepart(p.x, p.y) && samepart(p.x, p.z)) {
      found.insert(truth_part[p.x]);
    }
  }
  double score = 0;
  for (long long p : found) {
    point s = start_pos[p];
    if (s.x*s.x+s.y*s.y < 100)
      score += part_weight[p];
  }
  double total_possible = 0;
  for (auto&t : part_weight) {
    point s = start_pos[t.first];
    if (s.x*s.x+s.y*s.y < 100)
      total_possible += t.second;
  }
  cout << score/total_possible << " from " << triples.size() << " triples of " << total_possible << endl;
}

void scoreOutsideTriples(vector<triple>&triples) {
  set<long long> found;
  for (auto p : triples) {
    load(triples.size(), "Outside triple score:", 1);
    if (samepart(p.x, p.y) && samepart(p.x, p.z)) {
      found.insert(truth_part[p.x]);
    }
  }
  double score = 0;
  for (long long p : found) {
    point s = start_pos[p];
    if (s.x*s.x+s.y*s.y > 100)
      score += part_weight[p];
  }
  double total_possible = 0;
  for (auto&t : part_weight) {
    point s = start_pos[t.first];
    if (s.x*s.x+s.y*s.y > 100)
      total_possible += t.second;
  }
  cout << score/total_possible << " from " << triples.size() << " triples of " << total_possible << endl;
}

void scoreTriples2(vector<triple>&triples) {
  set<long long> found;
  for (auto p : triples) {
    load(triples.size(), "Triple score 2:", 1);
    if (samepart(p.x, p.y) && samepart(p.x, p.z)) {
      int prev1 = -1, prev2 = -1, prev3 = -1, ok = 0;
      for (int i : truth_tracks[truth_part[p.x]]) {
	int m = metai[i];
	if (m != prev1) {
	  prev3 = prev2;
	  prev2 = prev1;
	  prev1 = m;
	}
	if (i == p.z) {
	  if (prev2 == metai[p.y] && prev3 == metai[p.x]) ok = 1;
	  break;
	}
      }
      if (ok)
	found.insert(truth_part[p.x]);
    }
  }
  double score = 0;
  for (long long p : found) {
    score += part_weight[p];
  }
  cout << score << " from " << triples.size() << " triples" << endl;
}




void scorePaths(vector<vector<int> >&paths) {
  int total_length = 0;
  map<long long, double> score_part;
  for (auto&path : paths) {
    load(paths.size(), "Path score:  ", 1);

    total_length += path.size();
    map<long long, int> count;
    for (int i : path)
      count[truth_part[i]]++;
    long long part;
    int ma = 0;
    for (auto p : count) {
      if (p.second >= ma) {
	ma = p.second;
	part = p.first;
      }
    }
    double w = 0;
    for (int j : path) {
      if (truth_part[j] == part) w += truth_weight[j];
    }
    score_part[part] = max(score_part[part], w);
  }
  double score = 0;
  for (auto p : score_part) if (p.first) score += p.second;
  cout << score << " from " << paths.size() << " paths with " << total_length << " hits" << endl;
}



void scorePaths2(vector<vector<int> >&paths) {
  int total_length = 0;
  map<long long, double> score_part;
  for (auto&path : paths) {
    load(paths.size(), "Path score 2:", 1);
    total_length += path.size();
    map<long long, int> count;
    for (int i : path)
      count[truth_part[i]]++;
    long long part;
    int ma = 0;
    for (auto p : count) {
      if (p.second >= ma) {
	ma = p.second;
	part = p.first;
      }
    }
    double w = 0;
    for (int i : truth_tracks[part]) {
      for (int j : path) {
	if (truth_part[j] == part && metai[j] == metai[i]) {
	  w += truth_weight[i];
	  break;
	}
      }
    }
    score_part[part] = max(score_part[part], w);
  }
  double score = 0;
  for (auto p : score_part) if (p.first) score += p.second;
  cout << score << " from " << paths.size() << " paths with " << total_length << " hits" << endl;
}



void scorePaths3(vector<vector<int> >&paths) {
  int total_length = 0;
  map<long long, double> score_part;
  for (auto&path : paths) {
    load(paths.size(), "Path score 3:", 1);
    total_length += path.size();
    set<long long> done;
    for (int i : path) {
      long long part = truth_part[i];
      if (done.count(part)) continue;
      done.insert(part);
      double w = 0;
      for (int j : path) {
	if (truth_part[j] == part) w += truth_weight[j];
      }
      score_part[part] = max(score_part[part], w);
    }
  }
  double score = 0;
  for (auto p : score_part) if (p.first) score += p.second;
  cout << score << " from " << paths.size() << " paths with " << total_length << " hits" << endl;
}


void scorePaths4(vector<vector<int> >&paths) {
  int total_length = 0;
  map<long long, double> score_part;
  for (auto&path : paths) {
    load(paths.size(), "Path score 4:", 1);
    total_length += path.size();
    set<long long> done;
    //int minmeta = 1e9;
    //for (int i : path) minmeta = min(minmeta, metai[i]);
    for (int i : path) {
      long long part = truth_part[i];
      if (done.count(part)) continue;
      done.insert(part);
      double w = 0;
      for (int i : truth_tracks[part]) {
	for (int j : path) {
	  if (truth_part[j] == part && metai[j] == metai[i] && dist(hits[i]-hits[j]) < 10) {
	    w += truth_weight[i];
	    break;
	  }
	}
      }
      /*
      int oks = 0;
      for (int j : path) {
	if (truth_part[j] == part)
	  oks++;
	else break;
      }
      if (oks >= 3)
	for (int i : truth_tracks[part])
	  if (metai[i] < minmeta)
	    w += truth_weight[i];
      */
      score_part[part] = max(score_part[part], w);
    }
  }
  double score = 0;
  for (auto p : score_part) if (p.first) score += p.second;
  cout << score << " from " << paths.size() << " paths with " << total_length << " hits" << endl;
}



void scorePaths5(vector<vector<int> >&paths) {
  int total_length = 0;
  map<long long, double> score_part;
  for (auto&path : paths) {
    load(paths.size(), "Path score 5:", 1);
    total_length += path.size();
    set<long long> done;
    for (int i : path) {
      long long part = truth_part[i];
      if (done.count(part)) continue;
      done.insert(part);
      double bestw = 0, w = 0;
      int badstreak = 0, goodstreak = 0;
      for (int j : path) {
	if (truth_part[j] == part) {
	  badstreak = 0;
	  goodstreak++;
	  if (metai_weight[part].count(metai[j]))
	    w += metai_weight[part][metai[j]];
	} else {
	  badstreak++;
	  if (badstreak == 1) {
	    if (goodstreak >= 3)
	      bestw = max(bestw, w);
	    w = 0;
	  }
	  goodstreak = 0;
	}
      }
      if (goodstreak >= 3)
	bestw = max(bestw, w);
      score_part[part] = max(score_part[part], bestw);
    }
  }
  double score = 0;
  for (auto p : score_part) if (p.first) score += p.second;
  cout << score << " from " << paths.size() << " paths with " << total_length << " hits" << endl;
}


void scoreAssignment(map<int, int>&assignment) {
  map<int, int> track_length;
  for (auto p : assignment)
    track_length[p.second]++;

  double score = 0, score2 = 0;
  for (auto it : truth_tracks) {
    map<int, int> c;
    for (int i : it.second)
      if (assignment.count(i))
	c[assignment[i]]++;
    int pick = -1, maxlen = -1;
    for (auto p : c)
      if (p.second*2 > maxlen &&
	  p.second*2 > track_length[p.first]) pick = p.first, maxlen = p.second*2;
    if (pick == -1) continue;

    for (int i : it.second)
      if (assignment[i] == pick) {
	if (maxlen > it.second.size()) score += truth_weight[i];
	else score2 += truth_weight[i];
      }
  }
  //"Final score: " should be the same as the official score (except for blacklisted electrons)
  cout << "Final score:  " << score << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
  cout << "Short score:  " << score+score2 << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
}


void scoreAssignment2(map<int, int>&assignment) {
  map<int, int> track_length;
  for (auto p : assignment)
    track_length[p.second]++;

  double score = 0, score2 = 0, score3 = 0, score4 = 0;
  for (auto it : truth_tracks) {
    map<int, int> c;
    for (int i : it.second) {
      if (assignment.count(i))
	c[assignment[i]]++;
    }

    int pick = -1, maxlen = -1;
    for (auto p : c)
      if (p.second*2 > maxlen &&
	  p.second*2 > track_length[p.first]) pick = p.first, maxlen = p.second*2;
    if (pick == -1) continue;
    map<int, set<int> > modules;
    for (int i : it.second)
      if (assignment[i] == pick)
	modules[metai[i]].insert(i);

    int len = 0, len2 = 0;
    double s = 0, s2 = 0;
    for (int i : it.second) {
      if (modules.count(metai[i])) {
	int found = 3;
	for (int j : modules[metai[i]]) {
	  if (dist(hits[i]-hits[j]) < 10) found = min(found, 1);
	  found = min(found, 2);
	}
	if (found == 1)
	  s += truth_weight[i], len++;
	else if (found == 2)
	  s2 += truth_weight[i], len2++;
      }
    }
    len2 += len;
    if (len*2 > it.second.size()) score += s;
    else score3 += s;
    if (len2*2 > it.second.size()) score2 += s2;
    else score4 += s2;
  }
  score2 += score;
  score4 += score3;
  score3 += score;
  score4 += score2;
  cout << "No-dup score: " << score << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
  cout << "Far score:    " << score2 << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
  cout << "Near all:     " << score3 << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
  cout << "Far all:      " << score4 << " from " << assignment.size() << " hits out of " << hits.size()-1 << endl;
}
