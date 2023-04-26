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
    enum Type {
        planet,
        station
    };
    int x, y, id;
    Type type;
    Node() {};
    Node(int, int, int, Type);
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
    static State initState(const vector<Node>&);
    static State generateState(const State& input_state);
    void changeState(double, int&, int, int);
};

/*イテレーション管理クラス*/
template<class STATE>
struct IterationControl {
    int iteration_counter;
    int swap_counter;
    double average_time;
    double start_time;
    IterationControl() : iteration_counter(0), swap_counter(0) {}
    /*山登り法*/
    STATE climb(double time_limit, STATE initial_state) {
        start_time = toki.gettime();
        average_time = 0;
        STATE best_state = initial_state;
        double time_stamp = start_time;
        cerr << "[INFO] - IterationControl::climb - Starts climbing...\n";
        while(time_stamp - start_time + average_time < time_limit) {
            STATE current_state = STATE::generateState(best_state);
            if(current_state.score > best_state.score) {
                swap(best_state, current_state);
                swap_counter++;
            }
            iteration_counter++;
            time_stamp = toki.gettime();
            average_time = (time_stamp - start_time) / iteration_counter;
        }
        cerr << "[INFO] - IterationControl::climb - Iterated " << iteration_counter << " times and swapped " << swap_counter << " times.\n";
        return best_state;
    }
    /*焼きなまし法*/
    STATE anneal(double time_limit, double temp_start, double temp_end, STATE initial_state) {
        assert(temp_start >= temp_end);
        start_time = toki.gettime();
        average_time = 0;
        STATE best_state = initial_state;
        double elapsed_time = 0;
        cerr << "[INFO] - IterationControl::anneal - Starts annealing...\n";
        const int rsize = best_state.output.route.size();
        while(elapsed_time + average_time < time_limit) {
            double normalized_time = elapsed_time / time_limit;
            double temp_current = pow(temp_start, 1.0 - normalized_time) * pow(temp_end, normalized_time);
            for(int i = 0; i < rsize-1; i++) {
                for(int j = 0; j < rsize-1; j++) {
                    iteration_counter++;
                    best_state.changeState(temp_current, swap_counter, i, j);
                }
            }
            elapsed_time = toki.gettime() - start_time;
            average_time = elapsed_time / iteration_counter;
        }
        cerr << "[INFO] - IterationControl::anneal - Iterated " << iteration_counter << " times and swapped " << swap_counter << " times.\n";
        return best_state;
    }
};

namespace Utils {
    int calcSquareDist(const Node& a, const Node& b);
    int calcWeightedSquareDist(const Node& a, const Node& b);
    bool isPlanet(const Node& node);
    bool isStart(const Node& node);
    pair<long long, long long> calcScore(const Output& output);
    long long calcScoreFromLength(long long length);
    vector<Node> initStations();
    vector<Node> solveInsertedTSP(const vector<Node>& stations);
    pair<vector<Node>, vector<Node>> goThroughStations(vector<Node>, int);
    vector<Node> optimizeStations(const vector<Node>&);
};

/* ============================================== 
    プロトタイプ宣言はここまで
   ============================================== */

#pragma endregion prototype_declaration

Node::Node(int x, int y, int id, Node::Type type) : 
    x(x), y(y), id(id), type(type) {
        ;
}

/*TODO: ここで入力を受け取る*/
void Input::read() {
    cin >> n >> m;
    planets.resize(n);
    for(int i = 0; i < n; i++) {
        cin >> planets[i].x >> planets[i].y;
        planets[i].id = i;
        planets[i].type = Node::Type::planet;
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
    res.output.route = Utils::solveInsertedTSP(vector<Node>());
    auto [score, length] = Utils::calcScore(res.output);
    res.score = score;
    res.length = length;
    return res;
}

State State::initState(const vector<Node>& stations) {
    State res;
    res.output.route = Utils::solveInsertedTSP(stations);
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
    return node.type == Node::Type::planet;
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

vector<Node> Utils::solveInsertedTSP(const vector<Node>& stations) {
    vector<Node> nodes, route;
    for(auto e: input.planets) if(e.id != 0) nodes.push_back(e);
    for(auto e: stations) nodes.push_back(e);
    route.push_back(input.planets.front());
    route.push_back(input.planets.front());
    shuffle(nodes.begin(), nodes.end(), ryuka.engine);
    for(auto e: nodes) {
        int min_dist = 1<<30, min_id = -1;
        for(int i = 0; i < route.size() - 1; i++) {
            int dist = Utils::calcSquareDist(e, route[i]) + Utils::calcSquareDist(e, route[i+1]);
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
        cluster[i] = i % input.m;
    }
    vector<int> prev;
    int count = 0;
    const int iter_max = 100;
    while(prev != cluster && count < iter_max) {
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
                int dist = Utils::calcSquareDist(input.planets[i], Node(w_x[j], w_y[j], -1, Node::Type::station));
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
    if(count >= iter_max) {
        cerr << "[WARNING] - Utils::initStations - Failed to k-means" << endl;
        /*
        for(int i = 0; i < input.n; i++) {
            cerr << cluster[i] << " " << input.planets[i].x << " " << input.planets[i].y << endl;
        }
        */
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
        } else {
            int z = ryuka.rand(input.n);
            w_x[i] = clamp(input.planets[z].x + ryuka.rand(100) - 50, 1, 999);
            w_y[i] = clamp(input.planets[z].y + ryuka.rand(100) - 50, 1, 999);
        }
    }
    for(int i = 0; i < input.m; i++) {
        res[i] = Node(w_x[i], w_y[i], i, Node::Type::station);
    }   
    return res;
}

pair<vector<Node>, vector<Node>> Utils::goThroughStations(vector<Node> route, int iter_max) {
    vector<Node> stations = Utils::initStations();
    vector<Node> nodes = input.planets;
    for(auto e: stations) nodes.push_back(e);
    for(int it = 0; it < iter_max; it++) {
        vector<Node> res;
        for(int i = 0; i < route.size() - 1; i++) {
            const Node& cur = route[i];
            const Node& nxt = route[i+1];
            res.push_back(cur);
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
                res.push_back(nodes[min_id]);
            }
        }
        res.push_back(route.back());
        route = res;
    }
    return {route, stations};
}

vector<Node> Utils::optimizeStations(const vector<Node>& route) {
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
    for(int i = 0; i < input.m; i++) {
        if(num[i] > 0) {
            w_x[i] /= num[i];
            w_y[i] /= num[i];
        }
        stations.push_back(Node(w_x[i], w_y[i], i, Node::Type::station));
    }
    return stations;
} 

int main(int argc, char* argv[]) {
    toki.init();
    input.read();   

    long long best_score = 0;
    State best;
    for(int t = 0; t < 40; t++) {
        IterationControl<State> sera;
        State ans = sera.anneal(0.02, 1e5, 1, State::initState());
        //State ans = sera.climb(0.8, State::initState());
        //State ans = State::initState();
        auto [route, stations] = Utils::goThroughStations(ans.output.route, 10);
        ans.output.route = route;
        ans.output.stations = Utils::optimizeStations(ans.output.route);
        ans.score = Utils::calcScore(ans.output).first;
        if(ans.score > best_score) {
            best_score = ans.score;
            best = ans;
        }
    }
    best.output.print();
    cerr << "[INFO] - main - MyScore = " << best.score << "\n";
}