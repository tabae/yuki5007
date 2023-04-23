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
    int x, y, id, to, from;
    Type type;
    Node() : to(-1), from(-1) {};
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
    static State generateState(const State& input_state);
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
        while(elapsed_time + average_time < time_limit) {
            double normalized_time = elapsed_time / time_limit;
            double temp_current = pow(temp_start, 1.0 - normalized_time) * pow(temp_end, normalized_time);
            STATE current_state = STATE::generateState(best_state);
            long long delta = current_state.score - best_state.score;
            if(delta > 0 || ryuka.pjudge(exp(1.0 * delta / temp_current)) ) {
                swap(best_state, current_state);
                swap_counter++;
            }
            iteration_counter++;
            elapsed_time = toki.gettime() - start_time;
            average_time = elapsed_time / iteration_counter;
        }
        cerr << "[INFO] - IterationControl::anneal - Iterated " << iteration_counter << " times and swapped " << swap_counter << " times.\n";
        return best_state;
    }
};

namespace Utils {
    int calcSquareDist(const Node& a, const Node& b);
    bool isPlanet(const Node& node);
    bool isStart(const Node& node);
    pair<long long, long long> calcScore(const Output& output);
    long long calcScoreFromLength(long long length);
    vector<Node> solveInsertedTSP(const vector<Node>& stations);
    vector<Node> insertStations(const vector<Node>& nodes);
    vector<Node> rearrangeRoute(const vector<Node>& nodes);
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

/*TODO: ここでinput_stateを変化させた解を作る（局所探索）*/
State State::generateState(const State& input_state) {
    State res = input_state;
    int i = ryuka.rand(res.output.route.size());
    int j = ryuka.rand(res.output.route.size());
    int i2 = res.output.route[i].to;
    int j2 = res.output.route[j].to;
    bool chk = (i != j) && (i != j2) && (j != i2); 
    if(chk) {
        res.length -= Utils::calcSquareDist(res.output.route[i], res.output.route[i2]);
        res.length -= Utils::calcSquareDist(res.output.route[j], res.output.route[j2]);
        res.length += Utils::calcSquareDist(res.output.route[i], res.output.route[j]);
        res.length += Utils::calcSquareDist(res.output.route[j2], res.output.route[i2]);
        res.output.route[i].to = j;
        res.output.route[j2].from = i2;
        int cur = j;
        while(true) {
            int nxt = res.output.route[cur].from;
            swap(res.output.route[cur].from, res.output.route[cur].to);
            if(cur == j) res.output.route[cur].from = i;
            if(cur == i2) {
                res.output.route[cur].to = j2;
                break;
            }
            cur = nxt;
        }
        res.score = Utils::calcScoreFromLength(res.length);
    }
    return res;
}

int Utils::calcSquareDist(const Node& a, const Node& b) {
    const int dx = a.x - b.x;
    const int dy = a.y - b.y;
    return dx * dx + dy * dy;
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
    Node cur = route.front();
    assert(Utils::isStart(cur));
    while(true) {
        Node nxt = route[cur.to];
        const long long d2 = Utils::calcSquareDist(cur, nxt);
        if(isPlanet(cur) && isPlanet(nxt)) {
            sum += input.a * input.a * d2; 
        } else if(!isPlanet(cur) && !isPlanet(nxt)) {
            sum += d2;
        } else {
            sum += input.a * d2;
        }
        cur = nxt;
        if(Utils::isStart(cur)) break;
    }
    long long res = (long long)(1e9 / (1e3 + sqrt(sum)));
    return {res, sum};
}

long long Utils::calcScoreFromLength(long long length) {
    long long res = (long long)(1e9 / (1e3 + sqrt(length)));
    return res;
}

vector<Node> Utils::rearrangeRoute(const vector<Node>& route) {
    vector<Node> res;
    Node cur = route.front();
    while(true) {
        Node nxt = route[cur.to];
        cur.to = route.size();
        res.push_back(cur);
        if(Utils::isStart(nxt)) {
            res.push_back(nxt);
            break;
        }
        cur = nxt;
    }
    return res;
}

vector<Node> Utils::solveInsertedTSP(const vector<Node>& stations) {
    vector<Node> nodes, route;
    for(auto e: input.planets) if(e.id != 0) nodes.push_back(e);
    for(auto e: stations) nodes.push_back(e);
    route.push_back(input.planets.front());
    route.push_back(input.planets.front());
    shuffle(route.begin(), route.end(), ryuka.engine);
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
    for(int i = 0; i < route.size() - 2; i++) route[i].to = i+1;
    for(int i = 1; i < route.size() - 1; i++) route[i].from = i-1;
    route[route.size()-2].to = 0;
    route[0].from = route.size() - 2;
    route.pop_back();
    return route;
}

vector<Node> Utils::insertStations(const vector<Node>& nodes) {
    vector<pair<int, int>> v;
    vector<Node> res = nodes;
    for(int i = 0; i < res.size(); i++) {
        v.push_back({Utils::calcSquareDist(res[i], res[res[i].to]), i});
    } 
    sort(v.begin(), v.end(), greater<pair<int,int>>());
    sort(v.begin(), v.begin() + min<int>(input.m, v.size()), [](auto l, auto r) {
        return l.second < r.second;
    });
    vector<pair<int, Node>> inserted;
    for(auto [_, i]: v) {
        int nx = min(res[i].x, res[res[i].to].x) + abs(res[i].x - res[res[i].to].x) / 2; 
        int ny = min(res[i].y, res[res[i].to].y) + abs(res[i].y - res[res[i].to].y) / 2; 
        inserted.push_back({i, Node(nx, ny, inserted.size(), Node::Type::station)});
        if(inserted.size() == input.m) break;
    }
    for(auto [i, node]: inserted) {
        int ri = res.size();
        node.to = res[i].to;
        res[i].to = ri;
        res.push_back(node);
    }
    return res;
}

int main(int argc, char* argv[]) {
    toki.init();
    input.read();   
    IterationControl<State> sera;
    State ans = sera.climb(0.9, State::initState());
    //State ans = State::initState();
    cerr << "[Debug] - main - Starts insertStations ..." << endl;
    ans.output.route = Utils::insertStations(ans.output.route);
    for(auto e: ans.output.route) {
        if(!Utils::isPlanet(e)) ans.output.stations.push_back(e);
    }
    sort(ans.output.stations.begin(), ans.output.stations.end(), [](auto l, auto r) {
        return l.id < r.id;
    });
    cerr << "[Debug] - main - Starts calcScore ..." << endl;
    ans.score = Utils::calcScore(ans.output).first;
    cerr << "[Debug] - main - Starts rearrangeRoute ..." << endl;
    ans.output.route = Utils::rearrangeRoute(ans.output.route);
    ans.output.print();
    cerr << "[INFO] - main - MyScore = " << ans.score << "\n";
}