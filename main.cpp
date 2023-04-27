#include <bits/stdc++.h>
#include <sys/time.h>
#include <atcoder/all>

using namespace std;

#pragma region prototype_declaration
/* ============================================== 
    プロトタイプ宣言はここから
   ============================================== */

/*乱数生成器*/
struct RandGenerator {
    random_device seed_gen;
    mt19937 engine;
    mt19937_64 engine64;
    static const int pshift = 1000000000;
    RandGenerator() : engine(seed_gen()), engine64(seed_gen()) {}
    /*mod以下の乱数を返す（32bit）*/
    int rand(int mod) {
        return engine() % mod;
    }
    /*mod以下の乱数を返す（64bit）*/
    long long randll(long long mod) {
        return engine64() % mod;
    } 
    /*確率pでTrueを返す*/
    bool pjudge(double p) {
        int p_int;
        if(p > 1) p_int = pshift;
        else p_int = p * pshift;
        return rand(pshift) < p_int;
    }
} ryuka;

/*タイマー*/
struct Timer {
    double global_start;
    /*現在の時刻を返す*/
    double gettime() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }
    void init() {
        global_start = gettime();
    }
    /*プログラム開始からの経過時間を返す*/
    double elapsed() {
        return gettime() - global_start;
    }
} toki;

struct Node {
    int x, y, id, is_planet;
    Node() {};
    Node(int, int, int, int);
};

struct Input {
    /*TODO: ここに入力変数を定義する*/
    int n, m;
    const int a = 5;
    vector<Node> planets;
    void read();
} input;

struct Output {
    /*TODO: ここに出力変数を定義する*/
    vector<Node> stations;
    vector<Node> route; 
    Output();
    void print();
};

/*解を管理するクラス*/
struct State {
    Output output;
    long long length;
    long long score;
    State() : score(0) {}
    static State initState();
    static State generateState(const State& input_state);
    void changeState(double, int&, int, int);
    void changeStateClimb(int&, int, int);
};

namespace Utils {
    vector<vector<int>> planetsDist;
    void initPlanetsDist();
    int calcSquareDistOnlyPlanets(const Node& a, const Node& b);
    int calcSquareDist(const Node& a, const Node& b);
    int calcWeightedSquareDist(const Node& a, const Node& b);
    bool isPlanet(const Node& node);
    bool isStart(const Node& node);
    pair<long long, long long> calcScore(const Output& output);
    long long calcScoreFromLength(long long length);
    vector<Node> initStations();
    vector<Node> solveInsertedTSP();
    vector<Node> goThroughStations(vector<Node>, const vector<Node>& stations, int);
    pair<vector<Node>, vector<Node>> optimizeStations(const vector<Node>&);
};

/*イテレーション管理クラス*/
template<class STATE>
struct IterationControl {
    int iteration_counter;
    int swap_counter;
    double start_time;
    IterationControl() : iteration_counter(0), swap_counter(0) {}
    /*山登り法*/
    STATE climb(int iter_max, STATE initial_state) {
        STATE best_state = initial_state;
        #ifdef DEBUG
        cerr << "[INFO] - IterationControl::climb - Starts climbing...\n";
        #endif
        const int rsize = best_state.output.route.size();
        for(int it = 0; it < iter_max; it++) {
            for(int i = 0; i < rsize-1; i++) {
                for(int j = i+2; j < rsize-1; j++) {
                    iteration_counter++;
                    best_state.changeStateClimb(swap_counter, i, j);
                }
            }
        }
        #ifdef DEBUG
        cerr << "[INFO] - IterationControl::climb - Iterated " << iteration_counter << " times and swapped " << swap_counter << " times.\n";
        #endif
        return best_state;
    }
    /*焼きなまし法*/
    STATE anneal(double time_limit, double temp_start, double temp_end, STATE initial_state) {
        assert(temp_start >= temp_end);
        start_time = toki.gettime();
        STATE best_state = initial_state;
        double elapsed_time = 0;
        #ifdef DEBUG
        cerr << "[INFO] - IterationControl::anneal - Starts annealing...\n";
        #endif
        const int rsize = best_state.output.route.size();
        while(elapsed_time < time_limit) {
            double normalized_time = elapsed_time / time_limit;
            double temp_current = pow(temp_start, 1.0 - normalized_time) * pow(temp_end, normalized_time);
            for(int i = 0; i < rsize-1; i++) {
                for(int j = 0; j < rsize-1; j++) {
                    iteration_counter++;
                    best_state.changeState(temp_current, swap_counter, i, j);
                }
            }
            elapsed_time = toki.gettime() - start_time;
        }
        #ifdef DEBUG
        cerr << "[INFO] - IterationControl::anneal - Iterated " << iteration_counter << " times and swapped " << swap_counter << " times.\n";
        #endif
        return best_state;
    }
};


/* ============================================== 
    プロトタイプ宣言はここまで
   ============================================== */

#pragma endregion prototype_declaration

Node::Node(int x, int y, int id, int is_planet) : 
    x(x), y(y), id(id), is_planet(is_planet) {
        ;
}

/*TODO: ここで入力を受け取る*/
void Input::read() {
    cin >> n >> m;
    planets.resize(n);
    for(int i = 0; i < n; i++) {
        cin >> planets[i].x >> planets[i].y;
        planets[i].id = i;
        planets[i].is_planet = 1;
    }
}

/*TODO：ここで出力変数を初期化する。vectorのメモリ確保など*/
Output::Output() {
}

/*TODO：ここで答えを出力する*/
void Output::print() {
    assert(stations.size() == input.m);
    for(auto e : stations) {
        cout << e.x << " " << e.y << endl;
    }
    cout << route.size() << endl;
    for(auto e : route) {
        cout << (Utils::isPlanet(e) ? 1 : 2) << " " << e.id + 1 << endl;
    }
}

/*TODO: ここで初期解を作成する*/
State State::initState() {
    State res;
    res.output.route = Utils::solveInsertedTSP();
    auto [score, length] = Utils::calcScore(res.output);
    res.score = score;
    res.length = length;
    return res;
}

void State::changeState(double temp_current, int &swap_counter, int i, int j) {
    if(i > j) swap(i, j);
    int i2 = i+1;
    int j2 = j+1;
    bool chk = (i != j) && (i != j2) && (j != i2); 
    if(!chk) return;
    long long org_score = score;
    long long new_len = length;
    new_len -= Utils::calcWeightedSquareDist(output.route[i], output.route[i2]);
    new_len -= Utils::calcWeightedSquareDist(output.route[j], output.route[j2]);
    new_len += Utils::calcWeightedSquareDist(output.route[i], output.route[j]);
    new_len += Utils::calcWeightedSquareDist(output.route[j2],output.route[i2]);
    long long new_score = Utils::calcScoreFromLength(new_len);
    long long delta = new_score - org_score;
    if(delta > 0 || ryuka.pjudge(exp(1.0 * delta / temp_current)) ) {
        swap_counter++;
        reverse(output.route.begin() + i2, output.route.begin() + j2);
        length = new_len;
        score = Utils::calcScoreFromLength(length);        
    }
}

void State::changeStateClimb(int &swap_counter, int i, int j) {
    const int i2 = i+1;
    const int j2 = j+1;
    const long long org_score = score;
    long long new_len = length;
    new_len -= Utils::calcWeightedSquareDist(output.route[i], output.route[i2]);
    new_len -= Utils::calcWeightedSquareDist(output.route[j], output.route[j2]);
    new_len += Utils::calcWeightedSquareDist(output.route[i], output.route[j]);
    new_len += Utils::calcWeightedSquareDist(output.route[j2],output.route[i2]);
    //const long long new_score = Utils::calcScoreFromLength(new_len);
    //const long long delta = new_score - org_score;
    const long long delta = new_len - length;
    if(delta < 0) {
        swap_counter++;
        reverse(output.route.begin() + i2, output.route.begin() + j2);
        length = new_len;
    }
}


/*TODO: ここでinput_stateを変化させた解を作る（局所探索）*/
State State::generateState(const State& input_state) {
    State res = input_state;
    int i = ryuka.rand(res.output.route.size() - 1);
    int j = ryuka.rand(res.output.route.size() - 1);
    if(j < i) swap(i, j);
    int i2 = i+1;
    int j2 = j+1;
    bool chk = (i != j) && (i != j2) && (j != i2); 
    if(chk) {
        res.length -= Utils::calcWeightedSquareDist(res.output.route[i], res.output.route[i2]);
        res.length -= Utils::calcWeightedSquareDist(res.output.route[j], res.output.route[j2]);
        res.length += Utils::calcWeightedSquareDist(res.output.route[i], res.output.route[j]);
        res.length += Utils::calcWeightedSquareDist(res.output.route[j2], res.output.route[i2]);
        reverse(res.output.route.begin() + i2, res.output.route.begin() + j2);
        res.score = Utils::calcScoreFromLength(res.length);
    }
    return res;
}


int Utils::calcSquareDist(const Node& a, const Node& b) {
    const int dx = a.x - b.x;
    const int dy = a.y - b.y;
    return dx * dx + dy * dy;
}

int Utils::calcSquareDistOnlyPlanets(const Node& a, const Node& b) {
    return Utils::planetsDist[a.id][b.id];
}

int Utils::calcWeightedSquareDist(const Node& a, const Node& b) {
    const int dx = a.x - b.x;
    const int dy = a.y - b.y;
    const int s = dx * dx + dy * dy;
    int res;
    if(Utils::isPlanet(a) && Utils::isPlanet(b)) res = s * input.a * input.a;
    else if(!Utils::isPlanet(a) && !Utils::isPlanet(b)) res = s;
    else res = s * input.a; 
    return res;
}

bool Utils::isStart(const Node& node) {
    return isPlanet(node) && node.id == 0;
}

bool Utils::isPlanet(const Node& node) {
    return node.is_planet;
}

/*TODO: ここでスコアを計算する*/
pair<long long, long long> Utils::calcScore(const Output& output) {
    long long sum = 0;
    const auto& route = output.route;
    for(int i = 0; i < route.size() - 1; i++) {
        const Node& cur = route[i];
        const Node& nxt = route[i+1];
        const long long d2 = Utils::calcSquareDist(cur, nxt);
        if(isPlanet(cur) && isPlanet(nxt)) {
            sum += input.a * input.a * d2; 
        } else if(!isPlanet(cur) && !isPlanet(nxt)) {
            sum += d2;
        } else {
            sum += input.a * d2;
        }
    }
    long long res = (long long)(1e9 / (1e3 + sqrt(sum)));
    return {res, sum};
}

long long Utils::calcScoreFromLength(long long length) {
    long long res = (long long)(1e9 / (1e3 + sqrt(length)));
    return res;
}

void Utils::initPlanetsDist() {
    Utils::planetsDist.resize(input.n, vector<int>(input.n));
    for(int i = 0; i < input.n; i++) {
        for(int j = i; j < input.n; j++) {
            const int dist = Utils::calcSquareDist(input.planets[i], input.planets[j]);
            Utils::planetsDist[i][j] = dist;
            Utils::planetsDist[j][i] = dist;
        }
    }
}

vector<Node> Utils::solveInsertedTSP() {
    vector<Node> route;
    route.reserve(input.n);
    vector<Node> nodes = input.planets;
    route.emplace_back(input.planets.front());
    route.emplace_back(input.planets.front());
    shuffle(nodes.begin(), nodes.end(), ryuka.engine);
    for(auto e: nodes) {
        int min_dist = 1<<30, min_id = -1;
        for(int i = 0; i < route.size() - 1; i++) {
            int dist = Utils::calcSquareDistOnlyPlanets(e, route[i]) + Utils::calcSquareDistOnlyPlanets(e, route[i+1]);
            if(min_dist > dist) {
                min_dist = dist;
                min_id = i + 1;
            }
        }
        assert(min_id != -1);
        route.insert(route.begin() + min_id, e);
    }
    return route;
}

vector<Node> Utils::initStations() {
    vector<int> cluster(input.n, 0);
    for(int i = 0; i < input.n; i++) {
        cluster[i] = ryuka.rand(input.m);
    }
    vector<int> prev;
    int count = 0;
    const int iter_max = 100;
    for(int it = 0; it < iter_max; it++) {
        vector<long long> w_x(input.m, 0);
        vector<long long> w_y(input.m, 0);
        vector<int> num(input.m, 0);
        for(int i = 0; i < input.n; i++) {
            w_x[cluster[i]] += input.planets[i].x;
            w_y[cluster[i]] += input.planets[i].y;
            num[cluster[i]]++;
        }
        for(int i = 0; i < input.m; i++) {
            if(num[i] > 0) {
                w_x[i] /= num[i];
                w_y[i] /= num[i];
            }
        }
        for(int i = 0; i < input.n; i++) {
            int min_dist = 1<<30, min_id = -1;
            for(int j = 0; j < input.m; j++) {
                int dist = Utils::calcSquareDist(input.planets[i], Node(w_x[j], w_y[j], -1, 0));
                if(dist < min_dist) {
                    min_dist = dist;
                    min_id = j;
                }
            }
            assert(min_id != -1);
            cluster[i] = min_id;
        }
        count++;
    }
    vector<Node> res(input.m);
    vector<long long> w_x(input.m, 0);
    vector<long long> w_y(input.m, 0);
    vector<int> num(input.m, 0);
    for(int i = 0; i < input.n; i++) {
        w_x[cluster[i]] += input.planets[i].x;
        w_y[cluster[i]] += input.planets[i].y;
        num[cluster[i]]++;
    }
    for(int i = 0; i < input.m; i++) {
        if(num[i] > 0) {
            w_x[i] /= num[i];
            w_y[i] /= num[i];
            w_x[i] = clamp(w_x[i] + ryuka.rand(100) - 50, 1LL, 999LL);
            w_y[i] = clamp(w_y[i] + ryuka.rand(100) - 50, 1LL, 999LL);
        } else {
            int z = ryuka.rand(input.n);
            w_x[i] = clamp(input.planets[z].x + ryuka.rand(100) - 50, 1, 999);
            w_y[i] = clamp(input.planets[z].y + ryuka.rand(100) - 50, 1, 999);
        }
    }
    for(int i = 0; i < input.m; i++) {
        res[i] = Node(w_x[i], w_y[i], i, 0);
    }   
    return res;
}

vector<Node> Utils::goThroughStations(vector<Node> route, const vector<Node>& stations, int iter_max) {
    vector<Node> nodes;
    nodes.reserve(input.n + input.m);
    for(auto e: input.planets) nodes.emplace_back(e);
    for(auto e: stations) nodes.emplace_back(e);
    for(int it = 0; it < iter_max; it++) {
        vector<Node> res;
        for(int i = 0; i < route.size() - 1; i++) {
            const Node& cur = route[i];
            const Node& nxt = route[i+1];
            res.emplace_back(cur);
            int min_dist = Utils::calcWeightedSquareDist(cur, nxt);
            int min_id = -1;
            for(int j = 0; j < nodes.size(); j++) {
                int dist = Utils::calcWeightedSquareDist(cur, nodes[j])
                            + Utils::calcWeightedSquareDist(nodes[j], nxt);
                if(min_dist > dist) {
                    min_dist = dist;
                    min_id = j;
                }
            }
            if(min_id >= 0) {
                res.emplace_back(nodes[min_id]);
            }
        }
        res.emplace_back(route.back());
        route = res;
    }
    return route;
}

pair<vector<Node>, vector<Node>> Utils::optimizeStations(const vector<Node>& route) {
    vector<int> w_x(input.m, 0);
    vector<int> w_y(input.m, 0);
    vector<int> num(input.m, 0);
    for(int i = 0; i < route.size(); i++) {
        if(Utils::isPlanet(route[i])) {
            if(i - 1 >= 0 && !Utils::isPlanet(route[i-1])) {
                w_x[route[i-1].id] += route[i].x;
                w_y[route[i-1].id] += route[i].y;
                num[route[i-1].id]++;
            } 
            if(i + 1 < route.size() && !Utils::isPlanet(route[i+1])){
                w_x[route[i+1].id] += route[i].x;
                w_y[route[i+1].id] += route[i].y;
                num[route[i+1].id]++;
            }
        }
    }
    vector<Node> stations;
    stations.reserve(input.m);
    for(int i = 0; i < input.m; i++) {
        if(num[i] > 0) {
            w_x[i] /= num[i];
            w_y[i] /= num[i];
        }
        stations.emplace_back(Node(w_x[i], w_y[i], i, 0));
    }
    vector<Node> ret_route = route;
    for(int i = 0; i < ret_route.size(); i++) {
        if(!Utils::isPlanet(ret_route[i])) {
            ret_route[i].x = stations[ret_route[i].id].x;
            ret_route[i].y = stations[ret_route[i].id].y;
        }
    }
    return {ret_route, stations};
} 

int main(int argc, char* argv[]) {
    toki.init();
    input.read();   
    Utils::initPlanetsDist();
    long long best_score = 0;
    State best;
    for(int t = 0; t < 10000000; t++) {
        if(t % 10 == 0) {
            if(toki.elapsed() > 0.9) break;
        }
        IterationControl<State> sera;
        //State ans = sera.anneal(0.01, 1e5, 1, State::initState());
        //State ans = sera.climb(0.0005, State::initState());
        //State ans = State::initState();
        State ans = sera.climb(6, State::initState());
        ans.output.stations = Utils::initStations();
        ans.output.route = Utils::goThroughStations(ans.output.route, ans.output.stations, 2);
        auto [_route, _stations] = Utils::optimizeStations(ans.output.route);
        ans.output.route = move(_route);
        ans.output.stations = move(_stations);
        //ans = sera.anneal(0.01, 1e5, 1, ans);
        //ans = sera.climb(0.0005, ans);
        ans = sera.climb(3, ans);
        ans.output.route = Utils::goThroughStations(ans.output.route, ans.output.stations, 2);
        auto [route, stations] = Utils::optimizeStations(ans.output.route);
        ans.output.route = move(route);
        ans.output.stations = move(stations);
        ans.score = Utils::calcScore(ans.output).first;
        if(ans.score > best_score) {
            best_score = ans.score;
            best = move(ans);
        }
    }
    best.output.print();
    cerr << "[INFO] - main - MyScore = " << best.score << "\n";
}