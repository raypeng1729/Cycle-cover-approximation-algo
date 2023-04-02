#include <bits/stdc++.h>
using namespace std;

#define int long long

struct DSU{
    int sz;
    vector<int> rt;
    DSU(int SZ) : sz(SZ){
        rt.assign(sz, -1);
    }
    int find(int i){
        return rt[i] == -1 ? i : rt[i] = find(rt[i]);
    }
    void join(int i, int j){
        if(find(i) == find(j)){
            return;
        }
        if(find(i) < find(j)) rt[j] = find(i);
        else rt[i] = find(j);
    }
    bool same(int i, int j){
        return find(i) == find(j);
    }
};

struct WeightGraph {
    static const int inf = INT_MAX;
    static const int maxn = 514;
    struct edge {
        int u, v, w; 
        edge(){}
        edge(int u, int v, int w): u(u), v(v), w(w) {}
    };
    int n, n_x;
    edge g[maxn * 2][maxn * 2];
    int lab[maxn * 2];
    int match[maxn * 2], slack[maxn * 2], st[maxn * 2], pa[maxn * 2];
    int flo_from[maxn * 2][maxn + 1], S[maxn * 2], vis[maxn * 2];
    vector<int> flo[maxn * 2];
    queue<int> q;
    int e_delta(const edge &e) { return lab[e.u] + lab[e.v] - g[e.u][e.v].w * 2; }
    void update_slack(int u, int x) { if (!slack[x] || e_delta(g[u][x]) < e_delta(g[slack[x]][x])) slack[x] = u; }
    void set_slack(int x) {
        slack[x] = 0;
        for (int u = 1; u <= n; ++u)
            if (g[u][x].w > 0 && st[u] != x && S[st[u]] == 0)
                update_slack(u, x);
    }
    void q_push(int x) {
        if (x <= n) q.push(x);
        else for (size_t i = 0; i < flo[x].size(); i++) q_push(flo[x][i]);
    }
    void set_st(int x, int b) {
        st[x] = b;
        if (x > n) for (size_t i = 0; i < flo[x].size(); ++i) set_st(flo[x][i], b);
    }
    int get_pr(int b, int xr) {
        int pr = find(flo[b].begin(), flo[b].end(), xr) - flo[b].begin();
        if (pr % 2 == 1) {
            reverse(flo[b].begin() + 1, flo[b].end());
            return (int)flo[b].size() - pr;
        }
        return pr;
    }
    void set_match(int u, int v) {
        match[u] = g[u][v].v;
        if (u <= n) return;
        edge e = g[u][v];
        int xr = flo_from[u][e.u], pr = get_pr(u, xr);
        for (int i = 0; i < pr; ++i) set_match(flo[u][i], flo[u][i ^ 1]);
        set_match(xr, v);
        rotate(flo[u].begin(), flo[u].begin() + pr, flo[u].end());
    }
    void augment(int u, int v) {
        for (; ; ) {
            int xnv = st[match[u]];
            set_match(u, v);
            if (!xnv) return;
            set_match(xnv, st[pa[xnv]]);
            u = st[pa[xnv]], v = xnv;
        }
    }
    int get_lca(int u, int v) {
        static int t = 0;
        for (++t; u || v; swap(u, v)) {
            if (u == 0) continue;
            if (vis[u] == t) return u;
            vis[u] = t;
            u = st[match[u]];
            if (u) u = st[pa[u]];
        }
        return 0;
    }
    void add_blossom(int u, int lca, int v) {
        int b = n + 1;
        while (b <= n_x && st[b]) ++b;
        if (b > n_x) ++n_x;
        lab[b] = 0, S[b] = 0;
        match[b] = match[lca];
        flo[b].clear();
        flo[b].push_back(lca);
        for (int x = u, y; x != lca; x = st[pa[y]])
            flo[b].push_back(x), flo[b].push_back(y = st[match[x]]), q_push(y);
        reverse(flo[b].begin() + 1, flo[b].end());
        for (int x = v, y; x != lca; x = st[pa[y]])
            flo[b].push_back(x), flo[b].push_back(y = st[match[x]]), q_push(y);
        set_st(b, b);
        for (int x = 1; x <= n_x; ++x) g[b][x].w = g[x][b].w = 0;
        for (int x = 1; x <= n; ++x) flo_from[b][x] = 0;
        for (size_t i = 0; i < flo[b].size(); ++i) {
            int xs = flo[b][i];
            for (int x = 1; x <= n_x; ++x)
                if (g[b][x].w == 0 || e_delta(g[xs][x]) < e_delta(g[b][x]))
                    g[b][x] = g[xs][x], g[x][b] = g[x][xs];
            for (int x = 1; x <= n; ++x)
                if (flo_from[xs][x]) flo_from[b][x] = xs;
        }
        set_slack(b);
    }
    void expand_blossom(int b) {
        for (size_t i = 0; i < flo[b].size(); ++i)
            set_st(flo[b][i], flo[b][i]);
        int xr = flo_from[b][g[b][pa[b]].u], pr = get_pr(b, xr);
        for (int i = 0; i < pr; i += 2) {
            int xs = flo[b][i], xns = flo[b][i + 1];
            pa[xs] = g[xns][xs].u;
            S[xs] = 1, S[xns] = 0;
            slack[xs] = 0, set_slack(xns);
            q_push(xns);
        }
        S[xr] = 1, pa[xr] = pa[b];
        for (size_t i = pr + 1; i < flo[b].size(); ++i) {
            int xs = flo[b][i];
            S[xs] = -1, set_slack(xs);
        }
        st[b] = 0;
    }
    bool on_found_edge(const edge &e) {
        int u = st[e.u], v = st[e.v];
        if (S[v] == -1) {
            pa[v] = e.u, S[v] = 1;
            int nu = st[match[v]];
            slack[v] = slack[nu] = 0;
            S[nu] = 0, q_push(nu);
        } else if (S[v] == 0) {
            int lca = get_lca(u, v);
            if (!lca) return augment(u,v), augment(v,u), true;
            else add_blossom(u, lca, v);
        }
        return false;
    }
    bool matching() {
        memset(S + 1, -1, sizeof(int) * n_x);
        memset(slack + 1, 0, sizeof(int) * n_x);
        q = queue<int>();
        for (int x = 1; x <= n_x; ++x)
            if (st[x] == x && !match[x]) pa[x] = 0, S[x] = 0, q_push(x);
        if (q.empty()) return false;
        for (; ; ) {
            while (q.size()) {
                int u = q.front(); q.pop();
             if (S[st[u]] == 1) continue;
                for (int v = 1; v <= n; ++v)
                    if (g[u][v].w > 0 && st[u] != st[v]) {
                        if (e_delta(g[u][v]) == 0) {
                            if (on_found_edge(g[u][v])) return true;
                        } else update_slack(u, st[v]);
                    }
            }
            int d = inf;
            for (int b = n + 1; b <= n_x; ++b)
                if (st[b] == b && S[b] == 1) d = min(d, lab[b] / 2);
            for (int x = 1; x <= n_x; ++x)
                if (st[x] == x && slack[x]) {
                    if (S[x] == -1) d = min(d, e_delta(g[slack[x]][x]));
                    else if (S[x] == 0) d = min(d, e_delta(g[slack[x]][x]) / 2);
                }
            for (int u = 1; u <= n; ++u) {
                if (S[st[u]] == 0) {
                    if (lab[u] <= d) return 0;
                    lab[u] -= d;
                } else if (S[st[u]] == 1) lab[u] += d;
            }
            for (int b = n + 1; b <= n_x; ++b)
                if (st[b] == b) {
                    if (S[st[b]] == 0) lab[b] += d * 2;
                    else if (S[st[b]] == 1) lab[b] -= d * 2;
                }
            q = queue<int>();
            for (int x = 1; x <= n_x; ++x)
                if (st[x] == x && slack[x] && st[slack[x]] != x && e_delta(g[slack[x]][x]) == 0)
                    if (on_found_edge(g[slack[x]][x])) return true;
            for (int b = n + 1; b <= n_x; ++b)
                if (st[b] == b && S[b] == 1 && lab[b] == 0) expand_blossom(b);
        }
        return false;
    }
    pair<long long, int> solve() {
        memset(match + 1, 0, sizeof(int) * n);
        n_x = n;
        int n_matches = 0;
        long long tot_weight = 0;
        for (int u = 0; u <= n; ++u) st[u] = u, flo[u].clear();
        int w_max = 0;
        for (int u = 1; u <= n; ++u)
            for (int v = 1; v <= n; ++v) {
                flo_from[u][v] = (u == v ? u : 0);
                w_max = max(w_max, g[u][v].w);
            }
        for (int u = 1; u <= n; ++u) lab[u] = w_max;
        while (matching()) ++n_matches;
        for (int u = 1; u <= n; ++u)
            if (match[u] && match[u] < u)
                tot_weight += g[u][match[u]].w;
        return make_pair(tot_weight, n_matches);
    }
    void add_edge(int ui, int vi, int wi) { g[ui][vi].w = g[vi][ui].w = wi; }
    void init(int _n) {
        n = _n;
        for (int u = 1; u <= n; ++u)
            for (int v = 1; v <= n; ++v)
                g[u][v] = edge(u, v, 0);
    }
};

struct Graph{
    struct edge{
        int u, v, w;
        edge(){}
        edge(int _u, int _v, int _w) : u(_u), v(_v), w(_w){}
        bool operator<(const edge &e){
            return w < e.w;
        }
    };
    struct idx{
        int v, id;
        idx(){}
        idx(int _v, int _id) : v(_v), id(_id){}
    };
    int n, B, B_square, threshold;
    vector<pair<int, int>> coord;
    vector<vector<idx>> adj_list;
    vector<vector<bool>> adj;
    vector<int> degree;
    vector<edge> edge_list;
    vector<bool> used;
    Graph(){
        cin >> n >> B;
        B_square = B * B;
        init();
        for(int i = 0; i < n; i++){
            cin >> coord[i].first >> coord[i].second;
        }
    }
    void init(){
        adj_list.resize(n);
        for(int i = 0; i < n; i++) adj_list.clear();
        coord.resize(n);
        adj.assign(n, vector<bool>(n, 0));
        degree.assign(n, 0);
        edge_list.clear();
    }
    int dist(int i, int j){
        int dx = coord[i].first - coord[j].first;
        int dy = coord[i].second - coord[j].second;
        return dx * dx + dy * dy;
    }
    void set_threshold(int th){
        threshold = th * th;
        for(int i = 0; i < n; i++){
            for(int j = i + 1; j < n; j++){
                if(dist(i, j) < threshold){
                    edge_list.push_back(edge(i, j, dist(i, j)));
                    adj[i][j] = adj[j][i] = 1;
                }
            }
        }
        sort(edge_list.begin(), edge_list.end());
    }
    void build_mst(){
        vector<vector<bool>> new_adj(n, vector<bool>(n, 0));
        DSU dsu(n);
        for(auto e : edge_list){
            if(dsu.same(e.u, e.v)){
                continue;
            }
            dsu.join(e.u, e.v);
            new_adj[e.u][e.v] = 1;
            new_adj[e.v][e.u] = 1;
            degree[e.u]++;
            degree[e.v]++;
        }
    }
    void add_matching_edge(){
        vector<int> odd_degree;
        for(int i = 0; i < n; i++){
            if(degree[i] & 1){
                odd_degree.push_back(i);
            }
        }
        int od = odd_degree.size();
        WeightGraph solver;
        solver.init(od);
        for(int i = 0; i < od; i++){
            for(int j = i + 1; j < od; j++){
                if(dist(odd_degree[i], odd_degree[j]) < threshold){
                    solver.add_edge(odd_degree[i] + 1, odd_degree[j] + 1, dist(odd_degree[i], odd_degree[j]));
                }
            }
        }
        solver.solve();
        for(int i = 1; i <= n; i++){
            if(solver.match[i] && i < solver.match[i]){
                int mat = solver.match[i];
                adj[i - 1][mat - 1] = adj[mat - 1][i - 1] = 1;
            }
        }
    }
    void dfs(int i, vector<int> &circuit){
        for(auto j : adj_list[i]){
            if(used[j.id]){
                continue;
            }
            used[j.id] = 1;
            dfs(j.v, circuit);
        }
        circuit.push_back(i);
    }
    int calculate(vector<int> v){
        int ans = v.size(), n = v.size();
        for(int i = 0; i < n; i++){
            int j = 0, k = 1, cnt = 0;
            while(j < n){
                int sum = 0;
                while(k < n && sum + dist(v[k - 1], v[k]) + (j < k - 1) * dist(v[j], v[k]) <= B * B){
                    sum += dist(j, k);
                    k++;
                }
                j = k;
                cnt++;
            }
            ans = min(ans, cnt);
            rotate(v.begin(), v.begin() + 1, v.end());
        }
        return ans;
    }
    int euler_circuit(){
        int id = 0, result = 0;
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(i == j){
                    continue;
                }
                if(adj[i][j]){
                    adj_list[i].push_back(idx(j, id++));
                }
            }
        }
        used.assign(id, 0);
        vector<bool> vis(n, 0);
        for(int i = 0; i < n; i++){
            if(!adj_list[i].size() || used[adj_list[i][0].id]){
                continue;
            }
            vector<int> circuit, shortcuted_circuit;
            dfs(i, circuit);
            for(int j : circuit){
                if(vis[j]){
                    continue;
                }
                shortcuted_circuit.push_back(j);
                vis[j] = 1;
            }
            result += calculate(shortcuted_circuit);
        }
        return result;
    }
};

void solve(){
    //format: two integers n and B and n pairs of integer representing the node's coordinate
    Graph G;
    int MCC = G.n;
    for(int i = 2; i <= G.n; i++){
        int threshold = G.B / i;
        G.init();
        G.set_threshold(threshold);
        G.build_mst();
        G.add_matching_edge();
        MCC = min(MCC, G.euler_circuit());
    }
    cout << MCC << endl;
}

signed main(){
    //freopen("input.txt", "r", stdin);
    solve();
}
