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
    long long calcScore(const Output& output);
    vector<Node> solveInsertedTSP(const vector<Node>& stations);
    vector<Node> insertStations(const vector<Node>& nodes);
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
    res.output.route = Utils::insertStations(res.output.route);
    for(auto e: res.output.route) {
        if(!Utils::isPlanet(e)) res.output.stations.push_back(e);
    }
    sort(res.output.stations.begin(), res.output.stations.end(), [](auto l, auto r) {
        return l.id < r.id;
    });
    res.score = Utils::calcScore(res.output);
    return res;
}

/*TODO: ここでinput_stateを変化させた解を作る（局所探索）*/
State State::generateState(const State& input_state) {
    State res = input_state;
    res.score = Utils::calcScore(res.output);
    return res;
}

int Utils::calcSquareDist(const Node& a, const Node& b) {
    const int dx = a.x - b.x;
    const int dy = a.y - b.y;
    return dx * dx + dy * dy;
}

bool Utils::isPlanet(const Node& node) {
    return node.type == Node::Type::planet;
}

/*TODO: ここでスコアを計算する*/
long long Utils::calcScore(const Output& output) {
    long long sum = 0;
    const auto& route = output.route;
    for(int i = 0; i < route.size()-1; i++) {
        const long long d2 = Utils::calcSquareDist(route[i], route[i+1]);
        if(isPlanet(route[i]) && isPlanet(route[i+1])) {
            sum += input.a * input.a * d2; 
        } else if(!isPlanet(route[i]) && !isPlanet(route[i+1])) {
            sum += d2;
        } else {
            sum += input.a * d2;
        }
    }
    long long res = (long long)(1e9 / (1e3 + sqrt(sum)));
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
    return route;
}

vector<Node> Utils::insertStations(const vector<Node>& nodes) {
    vector<pair<int, int>> v;
    vector<Node> res = nodes;
    for(int i = 0; i < res.size() - 1; i++) {
        v.push_back({Utils::calcSquareDist(res[i], res[i+1]), i});
    } 
    sort(v.begin(), v.end(), greater<pair<int,int>>());
    sort(v.begin(), v.begin() + min<int>(input.m, v.size()), [](auto l, auto r) {
        return l.second < r.second;
    });
    vector<pair<int, Node>> inserted;
    for(auto [_, i]: v) {
        int nx = min(res[i].x, res[i+1].x) + abs(res[i].x - res[i+1].x) / 2; 
        int ny = min(res[i].y, res[i+1].y) + abs(res[i].y - res[i+1].y) / 2; 
        inserted.push_back({i, Node(nx, ny, inserted.size(), Node::Type::station)});
        if(inserted.size() == input.m) break;
    }
    int sl = 0;
    for(auto [i, node]: inserted) {
        res.insert(res.begin() + i + 1 + sl, node);
        sl++;
    }
    return res;
}

int main(int argc, char* argv[]) {
    toki.init();
    input.read();   
    IterationControl<State> sera;
    /*山登りの場合は、山登りする時間を第一引数で渡す*/
    // State ans = sera.climb(1.8, State::initState());
    /*焼きなましの場合は、焼きなます時間を第一引数で、初期温度を第二引数で、終了温度を第三引数で渡す*/
    // State ans = sera.anneal(1.8, temp_start, temp_end, State::initState());
    State ans = State::initState();
    ans.output.print();
    cerr << "[INFO] - main - MyScore = " << ans.score << "\n";
}