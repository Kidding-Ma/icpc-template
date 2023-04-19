# icpc 模板

## CLion

CmakeList.txt

```txt
add_definitions(-O2)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
add_executable(a a.cpp)
add_executable(b b.cpp)
```

## 生成随机数

```cpp
mt19937 rng((unsigned int) chrono::steady_clock::now().time_since_epoch().count());
```

## 文件操作

```cpp
freopen("in.txt", "r", stdin);
freopen("out.txt", "w", stdout);
```

## Exgcd

```cpp
i64 exgcd(i64 a, i64 b, i64 &x, i64 &y) {
    if (!b) {
    	x = 1;
    	y = 0;
    	return a;
    }
    i64 d = exgcd(b, a % b, x, y);
    i64 t = x;
    x = y;
    y = t - (a / b) * y;
    return d;
}
```

## Fenwick

```cpp
template <typename T>
struct Fenwick {
    int n;
    vector<T> a;
    Fenwick(const int n = 0) : n(n), a(n, T()) {}
    void modify(int i, T x) {
        for (i += 1; i <= n; i += i & -i) {
            a[i - 1] += x;
        }
    }
    T get(int i) {
        T res = T();
        for ( ; i > 0; i -= i & -i) {
            res += a[i - 1];
        }
        return res;
    }
    T sum(int l, int r) { // [l, r)
        return get(r) - get(l);
    }
    int kth(T k) {
        int x = 0;
        for (int i = 1 << __lg(n); i; i >>= 1) {
            if (x + i <= n && k >= a[x + i - 1]) {
                x += i;
                k -= a[x - 1];
            }
        }
        return x;
    }
};

constexpr int inf = 1E9;
struct Min {
    int x;
    Min(int x = inf) : x(x) {}
 
    Min &operator+=(const Min &a) {
        x = std::min(a.x, x);
        return *this;
    }
};
```

## StringHash

```cpp
constexpr int p[2] = {223333333, 773333333};
constexpr int mod[2] = {1000000033, 1000002233};
 
constexpr int N = 5E5;
 
i64 pw[2][N + 1];
 
struct StringHash {
    vector<vector<i64>> h;
    StringHash() : h(2, vector<i64>(1)) {}
    void push_back(char ch) {
        h[0].push_back((h[0].back() * p[0] + ch) % mod[0]);
        h[1].push_back((h[1].back() * p[1] + ch) % mod[1]);
    }
    pair<int, int> get(int l, int r) { // [l, r)
        return {
            (h[0][r] - h[0][l] * pw[0][r - l] % mod[0] + mod[0]) % mod[0], 
            (h[1][r] - h[1][l] * pw[1][r - l] % mod[1] + mod[1]) % mod[1]
        };
    }
};

void init() {
    pw[0][0] = pw[1][0] = 1;
    for (int i = 1; i <= N; i++) {
        pw[0][i] = pw[0][i - 1] * p[0] % mod[0];
        pw[1][i] = pw[1][i - 1] * p[1] % mod[1];
    }
}
```

## Pollard_rho

```cpp
using i128 = __int128;
 
i64 POW(i64 a, i64 b, i64 mod) {
    i64 res = 1;
    while (b) {
        if (b & 1) {
            res = (i128) res * a % mod;
        }
        b >>= 1;
        a = (i128) a * a % mod;
    }
    return res;
}
 
bool isprime(i64 p) {
    if (p < 2) {
        return 0;
    }
    i64 d = p - 1, r = 0;
    while (!(d & 1)) {
        r++;
        d >>= 1;
    }
    int prime[] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
    for (int k = 0; k < 9; k++) {
        if (p == prime[k]) {
            return true;
        }
        i64 x = POW(prime[k], d, p);
        if (x == 1 || x == p - 1) {
            continue;
        }
        for (int i = 0; i < r - 1; i++) {
            x = (i128) x * x % p;
            if (x == p - 1) {
                break;
            }
        }
        if (x != p - 1) {
            return false;
        }
    }
    return true;
}
 
mt19937 rng((unsigned int) chrono::steady_clock::now().time_since_epoch().count());
 
i64 pollard_rho(i64 x) {
    i64 s = 0, t = 0;
    i64 c = (i64) rng() % (x - 1) + 1;
    i64 val = 1;
    for (int goal = 1; ; goal <<= 1, s = t, val = 1) {
        for (int step = 1; step <= goal; step++) {
            t = ((i128) t * t + c) % x;
            val = (i128) val * abs(t - s) % x;
            if (step % 127 == 0) {
                i64 g = gcd(val, x);
                if (g > 1) {
                    return g;
                }
            }
        }
        i64 g = gcd(val, x);
        if (g > 1) {
            return g;
        }
    }
}
 
unordered_map<i64, int> primes;
 
void get(i64 x) {
    if (x < 2) {
        return;
    }
    if (isprime(x)) {
        primes[x]++;
        return;
    }
    i64 mx = pollard_rho(x);
    get(x / mx);
    get(mx);
}
 
void getprimes(i64 x) {
    primes.clear();
    get(x);
}
 
vector<i64> fac;
 
void getfac(i64 x) {
    fac.clear();
    getprimes(x);
    vector<pair<i64, int>> tmp(primes.begin(), primes.end());
    int SIZE = tmp.size();
    function<void(int, i64)> dfs = [&](int id, i64 x) {
        if (id == SIZE) {
            fac.push_back(x);
            return;
        }
        i64 p = 1;
        for (int i = 0; i <= tmp[id].second; i++) {
            if (i != 0) {
                p *= tmp[id].first;
            }
            dfs(id + 1, x * p);
        }
    };
    dfs(0, 1);
}
```

## Set

```cpp
#include "ext/pb_ds/assoc_container.hpp"

using namespace __gnu_pbds;

template <typename Key, typename Cmp_Fn = std::less<Key>>
class Set : public tree<Key, null_type, Cmp_Fn, rb_tree_tag,
                        tree_order_statistics_node_update> {};

Set<int> s;
s.insert(1);
cout << s.order_of_key(1) << '\n'; // 找 1 是第几个元素
cout << *s.find_by_order(0) << '\n'; // 括号里是几就找第几个元素，注意这里是迭代器记得加星星号
```

## SparseTable

```cpp
template <typename T>
struct SparseTable {
    int n;
    function<T(const T&, const T&)> func;
    vector<vector<T>> st;
    SparseTable(const vector<T> &a, const function<T(const T&, const T&)> &f) : n(a.size()), func(f), st(__lg(n) + 1, vector<T>(n)) {
        st[0] = a;
        int lg = __lg(n);
        for (int i = 1; i <= lg; i++) {
            for (int j = 0; j <= n - (1 << i); j++) {
                st[i][j] = func(st[i - 1][j], st[i - 1][(1 << (i - 1)) + j]);
            }
        }  	    
    }
    T get(int l, int r) {// [l, r)
    	int lg = __lg(r - l);
    	return func(st[lg][l], st[lg][r - (1 << lg)]);
    }
};
```

## UnionFind

```cpp
struct UnionFind {
    int n;
    vector<int> f;
    UnionFind(const int &n) : n(n), f(n) { 
        iota(f.begin(), f.end(), 0); 
    }
    int get(int x) {
    	while (x != f[x]) {
      	    x = f[x] = f[f[x]];
    	}
    	return x;
    }
    bool unite(int x, int y) {
    	int gx = get(x), gy = get(y);
    	if (gx != gy) {
      	    f[gx] = gy;
      	    return 1;
    	}
    	return 0;
    }
    bool united(int x, int y) {
    	return get(x) == get(y);
    }
};
```

## Trie

```cpp
template <typename T>
struct Trie {
    int tot;
    const int N = 11; //注意！注意！
    vector<int> tag;
    vector<vector<int>> t;
    Trie(const int &n = 100) : tag(n), t(n, vector<int>(N)), tot(0) {}
    void insert(const T &s) {
        int SIZE = s.size();
        int u = 0;
        for (int i = 0; i < SIZE; i++) {
            int c = s[i]; // 注意！注意！
            if (t[u][c] == 0) {
                t[u][c] = ++tot;
            }
            u = t[u][c];
            tag[u]++;
        }
    }
    int get(const T &s) {
        int u = 0;
        int res = 0;
        int SIZE = s.size();
        for (int i = 0; i < SIZE; i++) {
            int c = s[i]; // 注意！注意！
            u = t[u][c];
            if (tag[u] > 0) {
                res++;
            } else {
                return res;
            }
        }
        return res;
    }
};
```

## LCA

```cpp
vector<int> dep(n + 1, 1);
vector<vector<int>> p(n + 1, vector<int>(__lg(n + 1) + 1));
vector<i64> sum(n + 1);
function<void(int, int)> dfs = [&](int cur, int pre) {
    p[cur][0] = pre;
    for (int i = 1; i <= __lg(dep[cur]); i++) {
        p[cur][i] = p[p[cur][i - 1]][i - 1];
    }
    for (auto [nex, w] : g[cur]) {
        if (nex != pre) {
            dep[nex] = dep[cur] + 1;
            sum[nex] = sum[cur] + w;
            dfs(nex, cur);
        }
    }
};
dfs(1, 0);
auto depdis = [&](int x, int y) {
    int LCA = lca(x, y);
    return dep[x] + dep[y] - 2 * dep[LCA];
};
auto get = [&](int u, int k) {
    for (int i = __lg(dep[u]); i >= 0; i--) {
        if (k & (1 << i)) {
            u = p[u][i];
        }
    }
    return u;
};
auto lca = [&](int x, int y) {
    if (dep[x] < dep[y]) {
        swap(x, y);
    }
    for (int i = __lg(dep[x] - dep[y]); i >= 0; i--) {
        if(dep[p[x][i]] >= dep[y]) {
            x = p[x][i];
        }
    }
    if (x == y) {
        return x;
    }
    for (int i = __lg(dep[x]); i >= 0; i--) {
        if(p[x][i] != p[y][i]) {
            x = p[x][i];
            y = p[y][i];
        }
    }        
    return p[x][0];
};
auto dis = [&](int x, int y) {
    int LCA = lca(x, y);
    return sum[x] + sum[y] - 2 * sum[LCA];
};
```

## Sieve

```cpp
vector<int> isprime;
vector<int> primes;

void sieve(int N) {
    isprime.assign(N + 1, 1);
    
    for (int i = 2; i <= N; i++) {
        if (isprime[i]) {
            primes.push_back(i);
        }
        for (auto p : primes) {
            if (i * p > N) {
                break;
            }
            isprime[i * p] = 0;
            if (i % p == 0) {
                break;
            }
        }
    }
}

// 分解法
vector<int> get(int x) {
    vector<int> a;
    for (auto p : primes) {
        if (1LL * p * p > x) {
            break;
        }
        if (x % p == 0) {
            a.push_back(p);
            while (x % p == 0) {
                x /= p;
            }
        }
    }
    if (x > 1) {
        a.push_back(x);
    }
    return a;
}
```

## 树的双重心

```cpp
auto getroot = [&]() {
    int root1 = 0, root2 = 0;
    vector<i64> sz(n + 1), mx(n + 1); // 注意！注意！
    mx[0] = 2E9;
    function<void(int, int)> dfs = [&](int cur, int pre) {
        sz[cur] = 1;
        mx[cur] = 0;
        for (auto nex : g[cur]) {
            if (nex != pre) {
                dfs(nex, cur);
                sz[cur] += sz[nex];
                mx[cur] = max(mx[cur], sz[nex]);
            }
        }
        mx[cur] = max(mx[cur], n - sz[cur]);
        if (mx[cur] == mx[root1]) {
            root2 = cur;
        }
        if (mx[cur] < mx[root1]) {
            root1 = cur;
            root2 = 0;
        }
    };
    dfs(1, 0);
    return (pair<int, int>){root1, root2};
};
```

## TreeHash

```cpp
constexpr int p[2] = {223333333, 773333333};
constexpr int mod[2] = {1000000033, 1000002233};
 
using Hash = pair<i64, i64>;
 
#define x first
#define y second

auto TreeHash = [&](int root) {
    vector<i64> sz(n + 1); // 注意！注意！
    function<Hash(int, int)> dfs = [&](int cur, int pre) {
        Hash res;
        sz[cur] = 1;
        vector<Hash> s;
        for (auto nex : g[cur]) {
      	    if (nex != pre) {
                s.push_back(dfs(nex, cur));
                sz[cur] += sz[nex];
      	    }
        }
        sort(s.begin(), s.end());
        for (auto si : s) {
      	    res.x = (res.x * p[0] + si.x) % mod[0];
      	    res.y = (res.y * p[1] + si.y) % mod[1];
        }
        res.x = (res.x * p[0] + sz[cur]) % mod[0];
        res.y = (res.y * p[1] + sz[cur]) % mod[1];
        return res;
    };   
    return dfs(root, 0);
};     
```

## 组合数

```cpp
template <typename T, typename U> 
T POW(T a, U b) {
    T res = 1;
    while (b) {
        if (b & 1) {
      	    res = res * a;
    	}
        a = a * a;
    	b >>= 1;
    }
    return res;
}

template <int M>
struct Mint {
    int x;
    constexpr Mint() : x{} {}
    constexpr Mint(const i64 &x) : x{norm(x % M)} {}
    
    constexpr int norm(const int &x) const {
        if (x < 0) {
            return x + M;
        }
        if (x >= M) {
            return x - M;
        }
        return x;
    }
    constexpr int val() const {
        return x;
    }
    explicit constexpr operator int() const {
        return x;
    }
    constexpr Mint operator-() const {
        Mint res;
        res.x = norm(M - x);
        return res;
    }
    constexpr Mint inv() const {
        assert(x != 0);
        return POW(*this, M - 2);
    }
    constexpr Mint &operator*=(const Mint &rhs) {
        x = 1LL * x * rhs.x % M;
        return *this;
    }
    constexpr Mint &operator+=(const Mint &rhs) {
        x = norm(x + rhs.x);
        return *this;
    }
    constexpr Mint &operator-=(const Mint &rhs) {
        x = norm(x - rhs.x);
        return *this;
    }
    constexpr Mint &operator/=(const Mint &rhs) {
        return *this *= rhs.inv();
    }
    friend constexpr Mint operator*(const Mint &lhs, const Mint &rhs) {
        Mint res = lhs;
        res *= rhs;
        return res;
    }
    friend constexpr Mint operator+(const Mint &lhs, const Mint &rhs) {
        Mint res = lhs;
        res += rhs;
        return res;
    }
    friend constexpr Mint operator-(const Mint &lhs, const Mint &rhs) {
        Mint res = lhs;
        res -= rhs;
        return res;
    }
    friend constexpr Mint operator/(const Mint &lhs, const Mint &rhs) {
        Mint res = lhs;
        res /= rhs;
        return res;
    }
    friend constexpr std::istream &operator>>(std::istream &is, Mint &a) {
        i64 v;
        is >> v;
        a = Mint(v);
        return is;
    }
    friend constexpr std::ostream &operator<<(std::ostream &os, const Mint &a) {
        return os << a.val();
    }
    friend constexpr bool operator==(const Mint &lhs, const Mint &rhs) {
        return lhs.val() == rhs.val();
    }
    friend constexpr bool operator!=(const Mint &lhs, const Mint &rhs) {
        return lhs.val() != rhs.val();
    }
};

constexpr int P = 998244353;
constexpr int M = 1000000007;
using mint = Mint<M>;

vector<mint> fac(1, 1);
vector<mint> ifac(1, 1);

mint C(int n, int m) {
    if (m < 0 || m > n) {
        return 0;
    }
    while ((int) fac.size() < n + 1) {
    	fac.push_back(fac.back() * (int) fac.size());
    	ifac.push_back(fac.back().inv());
    }
    return fac[n] * ifac[m] * ifac[n - m];
}
```

## 莫队

```cpp
int n, m;
cin >> n >> m;
vector<int> a(n);
for (int i = 0; i < n; i++) {
    cin >> a[i];
}
struct node {
    int l, r, id;
};
vector<node> q(m);
for (int i = 0; i < m; i++) {
    int l, r;
    cin >> l >> r;
    l--, r--;
    q[i] = {l, r, id};
}
int bk = sqrt(n);
sort(q.begin(), q.end(), [&](node a, node b) {
	return ((a.l / bk) ^ (b.l / bk)) ? a.l < b.l : (((a.l / bk) & 1) ? a.r < b.r : a.r > b.r);
});
int l = 0, r = -1;
int res = 0;
vector<int> ans(m);
auto add = [&](int pos) {
    int num = a[pos];
    
};
auto del = [&](int pos) {
    int num = a[pos];
    
};
for (int i = 0; i < m; i++) {
   auto &[ql, qr, qid] = q[i];
    while (r < qr) add(++r);
    while (l > ql) add(--l);
    while (l < ql) del(l++);
    while (r > qr) del(r--);
    ans[qid] = res;
}
for (int i = 0; i < m; i++) {
    cout << ans[i] << '\n';
}
```

## 树深子树和

```cpp
vector<vector<int>> g(n);
vector<int> dep(n, 1), cnt(n);
function<void(int, int)> dfs = [&](int cur, int pre) {
    cnt[cur]++;
    for (auto &nex : g[cur]) {
    	if (x != pre) {
            dep[nex] = dep[cur] + 1;
            dfs(nex, cur);
            cnt[cur] += cnt[nex];
      	}
    }
};
dfs(0, -1);
```

## __int128

```cpp
using i128 = __int128;
 
constexpr i64 N = 1E18;
 
void print(i128 x) {
    if (x <= N) {
    	cout << i64(x);
    	return;
    }
    cout << i64(x / N) << setw(18) << setfill('0') << i64(x % N);
}
```

## Python_Input

```cpp
import sys
input = sys.stdin.buffer.readline
```

## Unordered_map

```cpp
unordered_map<int, int> mp(1024);
mp.max_load_factor(0.25);
```

## Dijkstra 

```cpp
auto dijkstra = [&](int s, int t) {
    vector<int> vis(n);
    vector<int> dis(n, 1E9);
    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> q;    
    dis[s] = 0;
    q.push({dis[s], s});
    while (!q.empty()) {
        auto [u, v] = q.top();
        q.pop();
        if (v == t) {
            cout << dis[t] << '\n';
            return;
        }
        if (!vis[v]) {
            vis[v] = 1;
            for (auto &[x, w] : g[v]) {
                if (dis[x] > dis[v] + w) {
                    dis[x] = dis[v] + w;
                    q.push({dis[x], x});
                }
            }
        }
    }
    cout << -1 << '\n';
};
```

## 树的直径

```cpp
auto bfs = [&](int s) {
    vector<int> dis(n, -1);
    std::queue<int> q;
    q.push(s);
    dis[s] = 0;

    while (!q.empty()) {
        int x = q.front();
        q.pop();

        for (auto y : g[x]) {
            if (dis[y] == -1) {
                dis[y] = dis[x] + 1;
                q.push(y);
            }
        }
    }

    return std::max_element(dis.begin(), dis.end()) - dis.begin();
};
int l = bfs(0);
int r = bfs(l);
```

## KMP

循环节

```cpp
vector<int> nex(n + 1);
int j = 0;
for (int i = 1; i < n; i++) {
    while (j && s[i] != s[j]) {
        j = nex[j];
    }
    if (s[i] == s[j]) {
        j++;
    }
    nex[i + 1] = j;
}
int ans = n;
if (n - nex[n] > 0 && (n % (n - nex[n]) == 0)) {
    ans = n - nex[n];
}
```

## 计算几何

```cpp
const double Pi = acos(-1.0);

using Real = int;
using Point = complex<Real>;

#define x real
#define y imag

constexpr Real eps = 1E-9;

int sign(const Real &x) {
    return x < -eps ? -1 : x > eps ? 1 : 0;
}

bool equal(const Real &a, const Real &b) {
    return sign(a - b) == 0;
}

const Real inf = numeric_limits<Real>::max();

struct Line {
    Point a, b;
    Line() = default;
    Line(const Point &a, const Point &b) : a(a), b(b) {}
};

struct Segment : Line {
    Segment() = default;
    using Line::Line;
    Real len() const {
        return abs(a - b);
    }
};

// a -> b
struct Ray : Line {
    Ray() = default;
    using Line::Line;
};

struct Circle {
    Point p;
    Real r;
    Circle(const Point &p, const Real &r) : p(p), r(r) {}
};

Real cross(const Point &a, const Point &b) {
    return (conj(a) * b).y();
}

Real dot(const Point &a, const Point &b) {
    return (conj(a) * b).x();
}

Real area(const Point &a, const Point &b, const Point &c) {
    return abs(cross(b - a, c - a));
}

Point rotate(const Point &p, const double &rad) {
    Point res = conj(p) * Point(sin(rad), cos(rad));
    return Point(res.y(), res.x());
}

Point rotate90(const Point &p) {
    return Point(-p.y(), p.x());
}

Real len(const Point &p) {
    return sqrt(norm(p));
}

Real dist(const Point &a, const Point &b) {
    return len(a - b);
}

Real dist(const Point &p, const Line &l) {
    return fabs(cross(p - l.a, l.b - l.a)) / dist(l.a, l.b);
}

// 1 -> 点在左侧 -1 -> 点在右侧
int onLeft(const Point &p, const Line &l) {
    return sign(cross(l.b - l.a, p - l.a));
}

bool onSegment(const Point &p, const Segment &s) {
    return sign(cross(p - s.a, s.b - s.a)) == 0 && sign(dot(p - s.a, p - s.b)) <= 0;
}

Point projection(const Point &p, const Line &l) {
    return l.a + (l.b - l.a) * dot(l.b - l.a, p - l.a) / norm(l.b - l.a);
}

// 平行
bool collinear(const Line &l1, const Line &l2) {
    return sign(cross(l1.b - l1.a, l2.b - l2.a)) == 0 && sign(cross(l1.b - l1.a, l1.b - l2.a)) == 0;
}

// 共线
bool parallel(const Line &l1, const Line &l2) {
    return sign(cross(l1.a - l1.b, l2.a - l2.b)) == 0;
}

bool intersect(const Line &l, const Segment &s) { 
    return cross(l.b - l.a, s.a - l.a) * cross(l.b - l.a, s.b - l.a) <= 0; 
}
 

bool intersect(const Line &l1, const Line &l2) {
    return collinear(l1, l2) || !parallel(l1, l2);    
}

// 两线段严格相交，不包含端点
bool intersect(const Segment &l1, const Segment &l2) {
    return onLeft(l2.a, l1) * onLeft(l2.b, l1) < 0 && onLeft(l1.a, l2) * onLeft(l1.b, l2) > 0;
}

// 非严格哪里交都算交
bool intersect(const Segment &s1, const Segment &s2) {
    return onLeft(s2.a, s1) * onLeft(s2.b, s1) <= 0 && onLeft(s1.a, s2) * onLeft(s1.b, s2) <= 0;
}

bool intersect(const Ray &l1, const Ray &l2) {
    int sgn = sign(cross(l1.b - l1.a, l2.b - l2.a));
    if (sgn == 0) {
    	return false;
    }
    return onLeft(l1.a, l2) == sgn && onLeft(l2.a, l1) == -sgn;
}

// 对称
Point reflection(const Point &p, const Line &l) {
    return p + (projection(p, l) - p) * Real(2);
}

// 三点求圆心 共圆
Point circlecenter(const Point &a, const Point &b, const Point &c) {
    double a1 = b.x() - a.x(), b1 = b.y() - a.y(), c1 = (a1 * a1 + b1 * b1) / 2.0;
    double a2 = c.x() - a.x(), b2 = c.y() - a.y(), c2 = (a2 * a2 + b2 * b2) / 2.0;
    double d = a1 * b2 - a2 * b1;
    return Point(a.x() + (c1 * b2 - c2 * b1) / d, a.y() + (a1 * c2 - a2 * c1) / d);
}

Point cpoints(const Line &l, const Segment &s) {
    return {s.a + (s.b - s.a) * cross(l.b - l.a, l.b - s.a) / cross(l.b - l.a, s.b - s.a)};
}
 
Point cpoints(const Line &f, const Line &g) {
    return {g.a + (g.b - g.a) * cross(f.b - f.a, f.b - g.a) / cross(f.b - f.a, g.b - g.a)};
}

vector<Point> cpoints(const Circle &c, const Line &l) {
    Point pr = projection(c.p, l);
    if (equal(abs(pr - c.p), c.r)) {
        return {pr};
    }
    Point e = (l.b - l.a) / abs(l.b - l.a);
    Real k = sqrt(norm(c.r) - norm(pr - c.p));
    return {pr - e * k, pr + e * k};
}
 
vector<Point> cpoints(const Circle &f, const Circle &g) {
    Real d = abs(f.p - g.p), r = f.r + g.r;
    if (sign(d - r) > 0 || sign(d + f.r - g.r) < 0 || sign(d + g.r - f.r) < 0) {
        return {};
    }
    Real a = acos((norm(f.r) - norm(g.r) + norm(d)) / (2 * f.r * d));
    Real t = arg(g.p - f.p);
    Point p = f.p + polar(f.r, t + a);
    Point q = f.p + polar(f.r, t - a);
    if (equal(p.x(), q.x()) && equal(p.y(), q.y())) {
    	return {p};
    }
    return {p, q};
}
// IN : 1
// OUT : -1
// ON : 0
int contains(const Circle &c, const Point &p) {
    int sgn = sign(abs(c.p - p) - c.r);
    return sgn <= -1 ? 1 : sgn >= 1 ? -1 : 0;
}

double angle(const Point &a, const Point &b) {
    return acos(dot(a, b) / (len(a) * len(b)));
}

// ON : -1, OUT : 0, IN : 1
int contains(const Point &p, const vector<Point> &g) {
    if (g.size() == 1) {
        return p == g[0] ? -1 : 0;
    }
    if (g.size() == 2) {
        return onSegment(p, Segment(g[0], g[1]));
    }
    if (p == g[0]) {
        return -1;
    }
    if (onLeft(p, Segment(g[0], g[1])) == -1 || onLeft(p, Segment(g[0], g.back())) == 1) {
        return 0;
    }
    const auto cmp = [&](const Point &u, const Point &v) { 
        return onLeft(v, Segment(g[0], u)) == 1; 
    };
    const size_t i = lower_bound(g.begin() + 1, g.end(), p, cmp) - g.begin();
    if (i == 1) {
        return onSegment(p, Segment(g[0], g[i])) ? -1 : 0;
    }
    if (i == g.size() - 1 && onSegment(p, Segment(g[0], g[i]))) {
        return -1;
    }
    if (onSegment(p, Segment(g[i - 1], g[i]))) {
        return -1;
    }
    return onLeft(p, Segment(g[i - 1], g[i])) > 0;
}

// -1 : on, 0 : out, 1 : in
int contains(const Point &p, const vector<Point> &g) {
    int gsz = g.size();
    bool in = false;
    for(int i = 0; i < gsz; i++) {
        Point a = g[i] - p, b = g[(i + 1) % gsz] - p;
        if (a.y() > b.y()) swap(a, b);
        if (a.y() <= 0 && 0 < b.y() && cross(a, b) < 0) in = !in;
        if (cross(a, b) == 0 && dot(a, b) <= 0) return -1;
    }
    return in ? 1 : 0;
}

vector<Point> get(const vector<Point> &a) {
    int n = a.size();
 
    sort(a.begin(), a.end(), [&](auto a, auto b) {
        if (a.x() != b.x()) {
            return a.x() < b.x();
        }
        return a.y() < b.y();
    });
 
    vector<Point> b;
    for (int i = 0; i < n; i++) {
        while (b.size() > 1 && cross(b[(int) b.size() - 1] - b[(int) b.size() - 2], a[i] - b[(int) b.size() - 2]) <= 0) {
            b.pop_back();
        }
        b.push_back(a[i]);
    }
    int k = b.size();
    for (int i = n - 2; i >= 0; i--) {
        while (b.size() > k && cross(b[(int) b.size() - 1] - b[(int) b.size() - 2], a[i] - b[(int) b.size() - 2]) <= 0) {
            b.pop_back();
        }
        b.push_back(a[i]);        
    }
    b.pop_back();
    return b;
}
```

## Bitset

```cpp
bitset() // 每一位都是 false。
bitset(unsigned long val) // 设为 val 的二进制形式。
bitset(const string& str) // 设为 01 串 str。
count() // 返回 true 的数量。
size() // 返回 bitset 的大小。
test(pos) // 它和 vector 中的 at() 的作用是一样的，和 [] 运算符的区别就是越界检查。
any() // 若存在某一位是 true 则返回 true，否则返回 false。
none() // 若所有位都是 false 则返回 true，否则返回 false。
all() // C++11，若所有位都是 true 则返回 true，否则返回 false。
set() // 将整个 bitset 设置成 true。
set(pos, val = true) // 将某一位设置成 true/false。
reset() // 将整个 bitset 设置成 false。
reset(pos) // 将某一位设置成 false。相当于 set(pos, false)。
flip() // 翻转每一位。（0\leftrightarrow1，相当于异或一个全是 1 的 bitset）
flip(pos) // 翻转某一位。
to_string() // 返回转换成的字符串表达。
to_ulong() // 返回转换成的 unsigned long 表达 (long 在 NT 及 32 位 POSIX 系统下与 int 一样，在 64 位 POSIX 下与 long long 一样）。
to_ullong() // C++11，返回转换成的 unsigned long long 表达。
一些文档中没有的成员函数：
_Find_first() // 返回 bitset 第一个 true 的下标，若没有 true 则返回 bitset 的大小。
_Find_next(pos) // 返回 pos 后面（下标严格大于 pos 的位置）第一个 true 的下标，若 pos 后面没有 true 则返回 bitset 的大小。
```

优化埃氏筛

```cpp
bitset<N> vis;

void Prime(int n) {
    vis.set();
    vis[0] = vis[1] = 0;
    for (int i = 2; i * i <= n; i++) {
        if (vis[i]) {
            for (int j = i << 1; j <= n; j += i) vis[j] = 0;
        }
    }
}
```

## 最大字段和

```cpp
template<class Info,
        class Merge = std::plus<Info>>
struct SegmentTree {
    const int n;
    const Merge merge;
    std::vector<Info> info;
    SegmentTree(int n) : n(n), merge(Merge()), info(4 << std::__lg(n)) {}
    SegmentTree(std::vector<Info> init) : SegmentTree(init.size()) {
        std::function<void(int, int, int)> build = [&](int p, int l, int r) {
            if (r - l == 1) {
                info[p] = init[l];
                return;
            }
            int m = (l + r) / 2;
            build(2 * p, l, m);
            build(2 * p + 1, m, r);
            pull(p);
        };
        build(1, 0, n);
    }
    void pull(int p) {
        info[p] = merge(info[2 * p], info[2 * p + 1]);
    }
    void modify(int p, int l, int r, int x, const Info &v) {
        if (r - l == 1) {
            info[p] = v;
            return;
        }
        int m = (l + r) / 2;
        if (x < m) {
            modify(2 * p, l, m, x, v);
        } else {
            modify(2 * p + 1, m, r, x, v);
        }
        pull(p);
    }
    void modify(int p, const Info &v) {
        modify(1, 0, n, p, v);
    }
    Info rangeQuery(int p, int l, int r, int x, int y) {
        if (l >= y || r <= x) {
            return Info();
        }
        if (l >= x && r <= y) {
            return info[p];
        }
        int m = (l + r) / 2;
        return merge(rangeQuery(2 * p, l, m, x, y), rangeQuery(2 * p + 1, m, r, x, y));
    }
    Info rangeQuery(int l, int r) {
        return rangeQuery(1, 0, n, l, r);
    }
};
 
struct Info {
    i64 sum;
    i64 lmx;
    i64 rmx;
    i64 mx;
    Info(const i64 &x = 0) : sum(x), lmx(x), rmx(x), mx(x) {}
};
 
Info operator+(const Info &a, const Info &b) {
    Info res;
    res.sum = a.sum + b.sum;
    res.lmx = max(a.sum + b.lmx, a.lmx);
    res.rmx = max(b.sum + a.rmx, b.rmx);
    res.mx = max({a.mx, b.mx, a.rmx + b.lmx});
    return res;
}
```

## LazySegmentTree

区间修改，区间求和

```cpp
template<class Info, class Tag>
struct LazySegmentTree {
    const int n;
    std::vector<Info> info;
    std::vector<Tag> tag;
    LazySegmentTree(int n) : n(n), info(4 << std::__lg(n)), tag(4 << std::__lg(n)) {}
    LazySegmentTree(std::vector<Info> init) : LazySegmentTree(init.size()) {
        std::function<void(int, int, int)> build = [&](int p, int l, int r) {
            if (r - l == 1) {
                info[p] = init[l];
                return;
            }
            int m = (l + r) / 2;
            build(2 * p, l, m);
            build(2 * p + 1, m, r);
            pull(p);
        };
        build(1, 0, n);
    }
    void pull(int p) {
        info[p] = info[2 * p] + info[2 * p + 1];
    }
    void apply(int p, const Tag &v) {
        info[p].apply(v);
        tag[p].apply(v);
    }
    void push(int p) {
        apply(2 * p, tag[p]);
        apply(2 * p + 1, tag[p]);
        tag[p] = Tag();
    }
    void modify(int p, int l, int r, int x, const Info &v) {
        if (r - l == 1) {
            info[p] = v;
            return;
        }
        int m = (l + r) / 2;
        push(p);
        if (x < m) {
            modify(2 * p, l, m, x, v);
        } else {
            modify(2 * p + 1, m, r, x, v);
        }
        pull(p);
    }
    void modify(int p, const Info &v) {
        modify(1, 0, n, p, v);
    }
    Info rangeQuery(int p, int l, int r, int x, int y) {
        if (l >= y || r <= x) {
            return Info();
        }
        if (l >= x && r <= y) {
            return info[p];
        }
        int m = (l + r) / 2;
        push(p);
        return rangeQuery(2 * p, l, m, x, y) + rangeQuery(2 * p + 1, m, r, x, y);
    }
    Info rangeQuery(int l, int r) {
        return rangeQuery(1, 0, n, l, r);
    }
    void rangeApply(int p, int l, int r, int x, int y, const Tag &v) {
        if (l >= y || r <= x) {
            return;
        }
        if (l >= x && r <= y) {
            apply(p, v);
            return;
        }
        int m = (l + r) / 2;
        push(p);
        rangeApply(2 * p, l, m, x, y, v);
        rangeApply(2 * p + 1, m, r, x, y, v);
        pull(p);
    }
    void rangeApply(int l, int r, const Tag &v) {
        return rangeApply(1, 0, n, l, r, v);
    }
};
 
struct Tag {
    int x;
    Tag(int x = 0) : x(x) {};
    void apply(const Tag &t) {
        x += t.x;
    }
};
 
struct Max {
    int x;
    Max(int x = -1E9) : x(x) {}
    void apply(const Tag &t) {
        x += t.x;
    }
};
 
Max operator+(const Max &a, const Max &b) {
    return std::max(a.x, b.x);
}
```

## 区间求和

```cpp
template<class Info, class Tag>
struct LazySegmentTree {
    const int n;
    std::vector<Info> info;
    std::vector<Tag> tag;
    LazySegmentTree(int n) : n(n), info(4 << std::__lg(n)), tag(4 << std::__lg(n)) {}
    LazySegmentTree(std::vector<Info> init) : LazySegmentTree(init.size()) {
        std::function<void(int, int, int)> build = [&](int p, int l, int r) {
            if (r - l == 1) {
                info[p] = init[l];
                return;
            }
            int m = (l + r) / 2;
            build(2 * p, l, m);
            build(2 * p + 1, m, r);
            pull(p);
        };
        build(1, 0, n);
    }
    void pull(int p) {
        info[p] = info[2 * p] + info[2 * p + 1];
    }
    void apply(int p, int l, int r, const Tag &v) {
        info[p].apply(l, r, v);
        tag[p].apply(v);
    }
    void push(int p, int l, int r) {
        int m = (l + r) / 2;
        apply(2 * p, l, m, tag[p]);
        apply(2 * p + 1, m, r, tag[p]);
        tag[p] = Tag();
    }
    void modify(int p, int l, int r, int x, const Info &v) {
        if (r - l == 1) {
            info[p] = v;
            return;
        }
        int m = (l + r) / 2;
        push(p, l, r);
        if (x < m) {
            modify(2 * p, l, m, x, v);
        } else {
            modify(2 * p + 1, m, r, x, v);
        }
        pull(p);
    }
    void modify(int p, const Info &v) {
        modify(1, 0, n, p, v);
    }
    Info rangeQuery(int p, int l, int r, int x, int y) {
        if (l >= y || r <= x) {
            return Info();
        }
        if (l >= x && r <= y) {
            return info[p];
        }
        int m = (l + r) / 2;
        push(p, l, r);
        return rangeQuery(2 * p, l, m, x, y) + rangeQuery(2 * p + 1, m, r, x, y);
    }
    Info rangeQuery(int l, int r) {
        return rangeQuery(1, 0, n, l, r);
    }
    void rangeApply(int p, int l, int r, int x, int y, const Tag &v) {
        if (l >= y || r <= x) {
            return;
        }
        if (l >= x && r <= y) {
            apply(p, l, r, v);
            return;
        }
        int m = (l + r) / 2;
        push(p, l, r);
        rangeApply(2 * p, l, m, x, y, v);
        rangeApply(2 * p + 1, m, r, x, y, v);
        pull(p);
    }
    void rangeApply(int l, int r, const Tag &v) {
        return rangeApply(1, 0, n, l, r, v);
    }
};
 
struct Tag {
    i64 x;
    Tag(i64 x = 0) : x(x) {};
    void apply(const Tag &t) {
        x += t.x;
    }
};
 
struct Info {
    i64 x;
    Info(i64 x = 0) : x(x) {}
    void apply(const int &l, const int &r, const Tag &t) {
        x += (r - l) * t.x;
    }
};
 
Info operator+(const Info &a, const Info &b) {
    return a.x + b.x;
}
```

## ChthollyTree

```cpp
template <typename T> 
struct ChthollyTree {
    map<int, T> mp;
    ChthollyTree() { mp[0] = 0; }
    void split(int x) {
        auto it = prev(mp.upper_bound(x));
        mp[x] = it->second;
    }
    void assign(int l, int r, T x) {
        split(l);
        split(r + 1);
        auto it = mp.find(l);
        while (it->first != r + 1) {
        it = mp.erase(it);
        }
        mp[l] = x;
    }
    void add(int l, int r, T x) {
        split(l);
        split(r + 1);
        auto it = mp.find(l);
        while (it->first != r + 1) {
            it->second += x;
            it = next(it);
        }
    }
    T get(int l, int r, int x, int y) {
        i64 res = 0;
        split(l);
        split(r + 1);
        auto it = mp.find(l);
        while (it->first != r + 1) {
            res = (res + POW(it->second, x, y) * (next(it)->first - it->first) % y) % y;
            it = next(it);
        }
        return res;
    }
    T nth_element(int l, int r, int k) {
        split(l);
        split(r + 1);
        auto it = mp.find(l);
        vector<pair<T, T>> a;
        while (it->first != r + 1) {
            a.push_back({it->second, next(it)->first - it->first});
            it = next(it);
        }
        sort(a.begin(), a.end());
        for (int i = 0, cnt = 0;; i++) {    
            if (k <= (cnt += a[i].second)) {
                return a[i].first;
            }
        }
    }
};
```

## Crt

$n\equiv a_i(\mod r_i)$

```cpp
i64 crt(const vector<i64> &a, const vector<int> &r) {
    i64 n = 1, ans = 0;
    int k = r.size();
    for (int i = 0; i < k; i++) {
    	n *= r[i];
    }
    for (int i = 0; i < k; i++) {
    	i64 m = n / r[i], b, y;
    	exgcd(m, r[i], b, y);
    	ans = (ans + a[i] * m * b % n) % n;
    }
    return (ans % n + n) % n;
}
```

## 主席树
```cpp
struct node {
    node *l, *r;
    i64 sum;
    node() : l(nullptr), r(nullptr), sum(0) {}
};

node *add(node *t, int l, int r, int x) {
    node *rt = new node();
    if (t) {
        *rt = *t;
    }
    rt->sum++;
    if (r > l) {
        int mid = (l + r) >> 1;
        if (x < mid) {
            rt->l = add(rt->l, l, mid, x);
        } else {
            rt->r = add(rt->r, mid + 1, r, x);
        }    
    }
    return rt;
}

int get(node *t1, node *t2, int l, int r, int k) {
    if (r == l) {
        return l - 1;
    }
    int mid = (l + r) >> 1;
    i64 res1 = 0;
    if (t2) {
        if (t2->l) {
           res1 = t2->l->sum;
        }
    }
    i64 res2 = 0;
    if (t1) {
        if (t1->l) {
            res2 = t1->l->sum;
        }
    }
    int num = res1 - res2;
    if (num >= k) {
        return get(t1 ? t1->l : nullptr, t2 ? t2->l : nullptr, l, mid, k);
    }
    return get(t1 ? t1->r : nullptr, t2 ? t2->r : nullptr, mid + 1, r, k - num);
}
vector<node *> tree(n + 1);
for (int i = 0; i < n; i++) {
    int x = lower_bound(v.begin(), v.end(), a[i]) - v.begin();
    tree[i + 1] = add(tree[i], 1, SIZE, x);
}
get(tree[l - 1], tree[r], 1, SIZE, k)
```


## 树上启发式合并

$\text{totcnt}$ 表示子树中出现了多少种不同的颜色， $\text{res}$ 表示子树中出现次数等于出现最多颜色出现次数的颜色数
```cpp
int ans = 0;
i64 res = 0, mxc = 0, totcnt = 0;
int tot = -1;

auto add = [&](int cur) {
    if (cnt[c[cur]] == 0) {
        totcnt++;
    }
    cnt[c[cur]]++;
    if (cnt[c[cur]] == mxc) {
        res++;
    } else if (cnt[c[cur]] > mxc) {
        mxc = cnt[c[cur]];
        res = 1;
    }
};

auto del = [&](int cur) {
    cnt[c[cur]]--;
    if (cnt[c[cur]] == 0) {
        totcnt--;
    }
    res = mxc = 0;
};

function<void(int, int)> dfs = [&](int cur, int pre) {
    l[cur] = ++tot;
    node[tot] = cur;
    siz[cur] = 1;
    for (auto &nex : g[cur]) {
        if (nex != pre) {
            dfs(nex, cur);
            siz[cur] += siz[nex];
            if (son[cur] == -1 || siz[son[cur]] < siz[nex]) {
                son[cur] = nex;
            }
        }
    }
    r[cur] = tot;
};

function<void(int, int, bool)> dfs_again = [&](int cur, int pre, bool ok) {
    for (auto &nex : g[cur]) {
        if (nex != pre && nex != son[cur]) {
            dfs_again(nex, cur, 0);
        }
    }
    if (son[cur] != -1) {
        dfs_again(son[cur], cur, 1);
    }
    for (auto &nex : g[cur]) {
        if (nex != pre && nex != son[cur]) {
            for (int i = l[nex]; i <= r[nex]; i++) {
                add(node[i]);
            }
        }
    }
    add(cur);
    ans += (res == totcnt);
    if (!ok) {
        for (int i = l[cur]; i <= r[cur]; i++) {
            del(node[i]);
        }
    }
};

dfs(0, -1);
dfs_again(0, -1, 0);
```
## zkwSegmentTree
```cpp
template<class Info,
        class Merge = std::plus<Info>>
struct SegmentTree {
    const int n;
    const Merge merge;
    std::vector<Info> info;
    int N;
    SegmentTree(int n) : n(n), merge(Merge()) {
        N = 1;
        while (N < n) {
            N <<= 1;
        }
        info.assign(2 * N, Info());
    }
    SegmentTree(std::vector<Info> &init) : SegmentTree(init.size()) {
        for (int i = 0; i < n; i++) {
            info[N + i] = init[i];
        }
        for (int i = N - 1; i; i--) {
            pull(i);
        }
    }
    void pull(int p) {
        info[p] = merge(info[2 * p], info[2 * p + 1]);
    }
    void modify(int p, const Info &v) {
        p += N;
        info[p] = v;
        for (int i = 1; (1 << i) <= N; i++) {
            pull(p >> i);
        }
    }
    Info rangeQuery(int l, int r) {
        Info res = Info();
        l += N;
        r += N;
        while (l < r) {
            if (l & 1) res = merge(res, info[l++]);
            if (r & 1) res = merge(info[--r], res);
            l >>= 1;
            r >>= 1;
        }
        return res;
    }
};
```
## 表达式板子
```cpp
auto lv = [&](const char &ch) {
    if (ch == '+') {
        return 1;
    }
    return 2;
};
    
int N = s.size();
vector<i64> NUM;
vector<char> OP;
    
auto work = [&]() {
    char op = OP.back();
    OP.pop_back();
    auto r = NUM.back();
    NUM.pop_back();
    auto l = NUM.back();
    NUM.pop_back();
    if (op == '*') {
        NUM.push_back(l * r);
    }
    if (op == '+') {
        NUM.push_back(l + r);
    }
};
    
auto run = [&](string &s) {
    for (int i = 0; i < N; ) {
        if (s[i] == '&' || s[i] == '$') {
            i++;
            continue;
        }
        if (isdigit(s[i])) {
            i64 res = 0;
            for ( ; i < N && isdigit(s[i]); i++) {
                res *= 10;
                res += s[i] - '0';
            }
            NUM.push_back(res);
        } else {
            if (s[i] == '(') {
                OP.push_back('(');
            } else if (s[i] == ')') {
                while (OP.back() != '(') {
                    work();
                }
                OP.pop_back();
            } else {
                while (!OP.empty() && OP.back() != '(' && lv(OP.back()) >= lv(s[i])) {
                    work();
                }
                OP.push_back(s[i]);
            }
            i++;
        }
    }
    while (!OP.empty()) {
        work();
    }
    i64 res = NUM.back();
    NUM.pop_back();
    return res;
};
```
