# $\text{icpc-模板}$

## $\text{CLion}$

$\text{CmakeList.txt}$

```txt
add_definitions(-O2)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
add_executable(a a.cpp)
add_executable(b b.cpp)
```

## $\text{random}$

```cpp
mt19937 rng((unsigned int) chrono::steady_clock::now().time_since_epoch().count());
```

## $\text{exgcd}$

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

## $\text{fenwick}$

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

## $\text{StringHash}$

```cpp
constexpr int B[2] = {223333333, 773333333};
constexpr int P[2] = {1000000033, 1000002233};
 
vector<i64> p[2];

void init(int N) {
    p[0].assign(N, 0);
    p[1].assign(N, 0);
    p[0][0] = p[1][0] = 1;
    for (int i = 1; i <= N; i++) {
        p[0][i] = p[0][i - 1] * B[0] % P[0];
        p[1][i] = p[1][i - 1] * B[1] % P[1];
    }
}

struct StringHash {
    vector<vector<i64>> h;
    StringHash() : h(2, vector<i64>(1)) {}
    void push_back(char ch) {
        h[0].push_back((h[0].back() * B[0] + ch) % P[0]);
        h[1].push_back((h[1].back() * B[1] + ch) % P[1]);
    }
    pair<int, int> to_pair(int l, int r) { // [l, r)
        return {
            (h[0][r] - h[0][l] * p[0][r - l] % P[0] + P[0]) % P[0], 
            (h[1][r] - h[1][l] * p[1][r - l] % P[1] + P[1]) % P[1]
        };
    }
    i64 to_i64(int l, int r) {
        return 
        (((h[0][r] - h[0][l] * p[0][r - l] % P[0] + P[0]) % P[0]) << 30) +  
        (h[1][r] - h[1][l] * p[1][r - l] % P[1] + P[1]) % P[1];        
    }
};
```

## $\text{pollard rho}$

```cpp
using i128 = __int128;
 
i64 power(i64 a, i64 b, i64 m) {
    i64 res = 1;
    for (; b; b >>= 1, a = i128(a) * a % m) {
        if (b & 1) {
            res = i128(res) * a % m;
        }
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
    for (auto a : prime) {
        if (p == a) {
            return true;
        }
        i64 x = power(a, d, p);
        if (x == 1 || x == p - 1) {
            continue;
        }
        for (int i = 0; i < r - 1; i++) {
            x = i128(x) * x % p;
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
    i64 c = i64(rng()) % (x - 1) + 1;
    i64 val = 1;
    for (int goal = 1; ; goal <<= 1, s = t, val = 1) {
        for (int step = 1; step <= goal; step++) {
            t = (i128(t) * t + c) % x;
            val = i128(val) * abs(t - s) % x;
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

unordered_map<i64, int> getprimes(i64 x) {
    unordered_map<i64, int> p;
    function<void(i64)> get = [&](i64 x) {
        if (x < 2) {
            return;
        }
        if (isprime(x)) {
            p[x]++;
            return;
        }
        i64 mx = pollard_rho(x);
        get(x / mx);
        get(mx);
    };
    get(x);
    return p;
}
 
vector<i64> getfac(i64 x) {
    vector<i64> fac;
    auto primes = getprimes(x);
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
    return fac;
}
```

## $\text{Set}$

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

## $\text{SparseTable}$

```cpp
template <typename T>
struct SparseTable {
    int n;
    function<T(const T&, const T&)> func;
    vector<vector<T>> a;
    SparseTable(const vector<T> &init, const function<T(const T&, const T&)> &f) : n(init.size()), func(f) {
        int lg = __lg(n);
        a.assign(lg + 1, vector<T>(n));
        a[0] = init;
        for (int i = 1; i <= lg; i++) {
            for (int j = 0; j <= n - (1 << i); j++) {
                a[i][j] = func(a[i - 1][j], a[i - 1][(1 << (i - 1)) + j]);
            }
        }  	    
    }
    T get(int l, int r) {// [l, r)
    	int lg = __lg(r - l);
    	return func(a[lg][l], a[lg][r - (1 << lg)]);
    }
};
```

## $\text{UnionFind}$

小的并入大的。

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
    	x = get(x), y = get(y);
    	if (x != y) {
      	    f[y] = x;
      	    return 1;
    	}
    	return 0;
    }
    bool united(int x, int y) {
    	return get(x) == get(y);
    }
};
```

## $\text{Trie}$

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

## $\text{lca}$

```cpp
int lg = __lg(n + 1);
vector<vector<int>> p(lg + 1, vector<int>(n + 1));
vector<int> dep(n + 1);
function<void(int, int)> dfs = [&](int cur, int pre) {
    dep[cur] = dep[pre] + 1;
    p[0][cur] = pre;
    for (int i = 1; i <= lg; i++) {
        p[i][cur] = p[i - 1][p[i - 1][cur]];
    }
    for (auto &nex : g[cur]) {
        if (nex != pre) {
            dfs(nex, cur);
        }
    }
};
dfs(s, 0);
auto lca = [&](int x, int y) {
    if (dep[x] < dep[y]) {
        swap(x, y);
    }
    for (int i = lg; i >= 0; i--) {
        if(dep[x] - dep[y] >= (1 << i)) {
            x = p[i][x];
        }
    }
    if (x == y) {
        return x;
    }
    for (int i = lg; i >= 0; i--) {
        if(p[i][x] != p[i][y]) {
            x = p[i][x];
            y = p[i][y];
        }
    }        
    return p[0][x];
};
```

```cpp
int lg = __lg(n + 1);
vector<vector<int>> g(n + 1);
vector<vector<i64>> sum(lg + 1, vector<i64>(n + 1)); 
vector<vector<int>> p(lg + 1, vector<int>(n + 1));

vector<int> dep(n + 1);
function<void(int, int)> dfs = [&](int cur, int pre) {
    dep[cur] = dep[pre] + 1;
    p[0][cur] = pre;
    sum[0][cur] = g[cur].size();
    for (int i = 1; i <= lg; i++) {
        p[i][cur] = p[i - 1][p[i - 1][cur]];
        sum[i][cur] = sum[i - 1][cur] + sum[i - 1][p[i - 1][cur]];
    }
    for (auto &nex : g[cur]) {
        if (nex != pre) {
            dfs(nex, cur);
        }
    }
};
dfs(1, 0);
auto dis = [&](int x, int y) {
    if (dep[x] < dep[y]) {
        swap(x, y);
    }
    i64 l = 0, r = 0;
    for (int i = lg; i >= 0; i--) {
        if(dep[x] - dep[y] >= (1 << i)) {
            l = l + sum[i][x];
            x = p[i][x];
        }
    }
    if (x == y) {
        return l + sum[0][x];
    }
    for (int i = lg; i >= 0; i--) {
        if(p[i][x] != p[i][y]) {
            l = l + sum[i][x];
            r = r + sum[i][y];
            x = p[i][x];
            y = p[i][y];
        }
    }        
    return l + sum[0][p[0][x]] + sum[0][y] + r + sum[0][x];
}; // 点权
```


## $\text{sieve}$

```cpp
vector<int> minp, primes;
 
void sieve(int n) {
    minp.assign(n + 1, 0);
    
    for (int i = 2; i <= n; i++) {
        if (minp[i] == 0) {
            minp[i] = i;
            primes.push_back(i);
        }
        
        for (auto p : primes) {
            if (i * p > n) {
                break;
            }
            minp[i * p] = p;
            if (p == minp[i]) {
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

vector<int> get1(int x) {
    vector<int> a;
    while (x > 1) {
        int p = minp[x];
        x /= p;
        a.push_back(p);
    }
	return a;
}
```

## $\text{Graphs}$

树的重心

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

树深子树和
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

树的直径

```cpp
auto bfs = [&](int s) {
    vector<int> dis(n, -1);
    queue<int> q;
    q.push(s);
    dis[s] = 0;

    while (!q.empty()) {
        int cur = q.front();
        q.pop();

        for (auto &nex : g[cur]) {
            if (dis[nex] == -1) {
                dis[nex] = dis[cur] + 1;
                q.push(y);
            }
        }
    }

    return max_element(dis.begin(), dis.end()) - dis.begin();
};
int l = bfs(0);
int r = bfs(l);
```

## $\text{TreeHash}$

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

## $\text{Comb}$

组合数

```cpp
template <typename T, typename U> 
T power(T a, U b) {
    T res = T(1);
    for (; b; b >>= 1, a *= a) {
        if (b & 1) {
      	    res = res * a;
    	}
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
        return power(*this, M - 2);
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
using mint = Mint<P>;

vector<mint> fac;
vector<mint> ifac;

mint C(int n, int m) {
    if (m < 0 || m > n) {
        return 0;
    }
    return fac[n] * ifac[m] * ifac[n - m];
}

void init(int n) {
    fac.assign(n + 1, mint());
    ifac.assign(n + 1, mint());
    
    fac[0] = 1;
    for (int i = 1; i <= n; i++) {
        fac[i] = fac[i - 1] * i;
    }

    ifac[n] = fac[n].inv();
    for (int i = n - 1; i >= 0; i--) {
        ifac[i] = ifac[i + 1] * (i + 1);
    }
}
```

## $\text{mo's algorithm}$

莫队

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

## $\text{__int128}$

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

## $\text{Python}$
$\text{str}$ 别用，最后多两个字符。

```cpp
import sys
input = sys.stdin.buffer.readline
```

## $\text{unordered_map}$

```cpp
unordered_map<int, int> mp(1024);
mp.max_load_factor(0.25);
```

## $\text{dijkstra}$

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


## $\text{kmp}$

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

## $\text{Geometry}$

计算几何

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

bool operator==(const Point &a, const Point &b) {
    return abs(a.x() - b.x()) <= eps && abs(a.y() - b.y()) <= eps;
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

double angle(const Point &a, const Point &b) {
  	return acos(dot(a, b) / (len(a) * len(b)));
}

struct Circle {
    Point p;
    Real r;
    Circle(const Point &p, const Real &r) : p(p), r(r) {}
};

Real cross(const Point &a, const Point &b) {
    return (conj(a) * b).imag();
}

Real dot(const Point &a, const Point &b) {
    return (conj(a) * b).real();
}

Real area(const Point &a, const Point &b, const Point &c) {
    return abs(cross(b - a, c - a));
}

// 逆
Point rotate(const Point &p, const double &rad) {
    return p * Point(cos(rad), sin(rad));
}

// 顺
Point rotate(const Point &p, const double &rad) {
    return p * Point(cos(rad), -sin(rad));
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
    long double a1 = b.x() - a.x(), b1 = b.y() - a.y(), c1 = (a1 * a1 + b1 * b1) / 2.0;
    long double a2 = c.x() - a.x(), b2 = c.y() - a.y(), c2 = (a2 * a2 + b2 * b2) / 2.0;
    long double d = a1 * b2 - a2 * b1;
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

## $\text{bitset}$

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

$\text{优化埃氏筛}$

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
## $\text{SegmentTree}$
最大字段和

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

## $\text{LazySegmentTree}$

区间修改，区间求和

```cpp
template<class Info, class Tag>
struct LazySegmentTree {
    int n;
    std::vector<Info> info;
    std::vector<Tag> tag;
    LazySegmentTree() : n(0) {}
    LazySegmentTree(int n_, Info v_ = Info()) {
        init(n_, v_);
    }
    template<class T>
    LazySegmentTree(std::vector<T> init_) {
        init(init_);
    }
    void init(int n_, Info v_ = Info()) {
        init(std::vector(n_, v_));
    }
    template<class T>
    void init(std::vector<T> init_) {
        n = init_.size();
        info.assign(4 << std::__lg(n), Info());
        tag.assign(4 << std::__lg(n), Tag());
        std::function<void(int, int, int)> build = [&](int p, int l, int r) {
            if (r - l == 1) {
                info[p] = init_[l];
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
    template<class F>
    int findFirst(int p, int l, int r, int x, int y, F pred) {
        if (l >= y || r <= x || !pred(info[p])) {
            return -1;
        }
        if (r - l == 1) {
            return l;
        }
        int m = (l + r) / 2;
        push(p);
        int res = findFirst(2 * p, l, m, x, y, pred);
        if (res == -1) {
            res = findFirst(2 * p + 1, m, r, x, y, pred);
        }
        return res;
    }
    template<class F>
    int findFirst(int l, int r, F pred) {
        return findFirst(1, 0, n, l, r, pred);
    }
    template<class F>
    int findLast(int p, int l, int r, int x, int y, F pred) {
        if (l >= y || r <= x || !pred(info[p])) {
            return -1;
        }
        if (r - l == 1) {
            return l;
        }
        int m = (l + r) / 2;
        push(p);
        int res = findLast(2 * p + 1, m, r, x, y, pred);
        if (res == -1) {
            res = findLast(2 * p, l, m, x, y, pred);
        }
        return res;
    }
    template<class F>
    int findLast(int l, int r, F pred) {
        return findLast(1, 0, n, l, r, pred);
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

区间求和

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

## $\text{ChthollyTree}$

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

## $\text{crt}$

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

## $\text{PersistentSegmentTree}$

主席树

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

## $\text{dsu on tree}$

树上启发式合并

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
## $\text{zkwSegmentTree}$
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
## $\text{Expression}$

表达式板子

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

## $\text{rope}$

可持续化数组

```cpp
//rope
#include "ext/rope"
using namespace __gnu_cxx;
rope<int> test;
test.push_back(x);//在末尾添加x
test.insert(pos,x);//在pos插入x　　
test.erase(pos,x);//从pos开始删除x个
test.copy(pos,len,x);//从pos开始到pos+len为止用x代替
test.replace(pos,x);//从pos开始换成x
test.substr(pos,x);//提取pos开始x个
test.at(x)/[x];//访问第x个元素
vector<rope<int>> a(n);
a[1] = a[0];
```

## $\text{Flow}$

```cpp
template <typename T>
struct Flow {
    struct _Edge {
        int to;
        T cap;
        _Edge(int to, T cap) : to(to), cap(cap) {}
    };
    
    int n;
    vector<_Edge> e;
    vector<vector<int>> g;
    vector<int> cur, h;
    
    Flow(const int &n = 0) : n(n), g(n) {}
    
    bool bfs(int s, int t) {
        h.assign(n, -1);
        queue<int> que;
        h[s] = 0;
        que.push(s);
        while (!que.empty()) {
            const int u = que.front();
            que.pop();
            for (int i : g[u]) {
                auto [v, c] = e[i];
                if (c > 0 && h[v] == -1) {
                    h[v] = h[u] + 1;
                    if (v == t) {
                        return true;
                    }
                    que.push(v);
                }
            }
        }
        return false;
    }
    
    T dfs(int u, int t, T f) {
        if (u == t) {
            return f;
        }
        auto r = f;
        for (int &i = cur[u]; i < int(g[u].size()); ++i) {
            const int j = g[u][i];
            auto [v, c] = e[j];
            if (c > 0 && h[v] == h[u] + 1) {
                auto a = dfs(v, t, min(r, c));
                e[j].cap -= a;
                e[j ^ 1].cap += a;
                r -= a;
                if (r == 0) {
                    return f;
                }
            }
        }
        return f - r;
    }
    void addEdge(int u, int v, T c) {
        g[u].push_back(e.size());
        e.emplace_back(v, c);
        g[v].push_back(e.size());
        e.emplace_back(u, 0);
    }
    T maxflow(int s, int t) {
        T ans = 0;
        while (bfs(s, t)) {
            cur.assign(n, 0);
            ans += dfs(s, t, numeric_limits<T>::max());
        }
        return ans;
    }
    
    vector<bool> minCut() {
        vector<bool> c(n);
        for (int i = 0; i < n; i++) {
            c[i] = (h[i] != -1);
        }
        return c;
    }
    
    struct Edge {
        int from;
        int to;
        T cap;
        T flow;
    };
    vector<Edge> edges() {
        vector<Edge> a;
        for (int i = 0; i < e.size(); i += 2) {
            Edge x;
            x.from = e[i + 1].to;
            x.to = e[i].to;
            x.cap = e[i].cap + e[i + 1].cap;
            x.flow = e[i + 1].cap;
            a.push_back(x);
        }
        return a;
    }
};
```

## $\text{HLD}$

```cpp
struct HLD {
    int n;
    std::vector<int> siz, top, dep, parent, in, out, seq;
    std::vector<std::vector<int>> adj;
    int cur;
    
    HLD() {}
    HLD(int n) {
        init(n);
    }
    void init(int n) {
        this->n = n;
        siz.resize(n);
        top.resize(n);
        dep.resize(n);
        parent.resize(n);
        in.resize(n);
        out.resize(n);
        seq.resize(n);
        cur = 0;
        adj.assign(n, {});
    }
    void addEdge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    void work(int root = 0) {
        top[root] = root;
        dep[root] = 0;
        parent[root] = -1;
        dfs1(root);
        dfs2(root);
    }
    void dfs1(int u) {
        if (parent[u] != -1) {
            adj[u].erase(std::find(adj[u].begin(), adj[u].end(), parent[u]));
        }
        
        siz[u] = 1;
        for (auto &v : adj[u]) {
            parent[v] = u;
            dep[v] = dep[u] + 1;
            dfs1(v);
            siz[u] += siz[v];
            if (siz[v] > siz[adj[u][0]]) {
                std::swap(v, adj[u][0]);
            }
        }
    }
    void dfs2(int u) {
        in[u] = cur++;
        seq[in[u]] = u;
        for (auto v : adj[u]) {
            top[v] = v == adj[u][0] ? top[u] : v;
            dfs2(v);
        }
        out[u] = cur;
    }
    int lca(int u, int v) {
        while (top[u] != top[v]) {
            if (dep[top[u]] > dep[top[v]]) {
                u = parent[top[u]];
            } else {
                v = parent[top[v]];
            }
        }
        return dep[u] < dep[v] ? u : v;
    }
    
    int dist(int u, int v) {
        return dep[u] + dep[v] - 2 * dep[lca(u, v)];
    }
    
    int jump(int u, int k) {
        if (dep[u] < k) {
            return -1;
        }
        
        int d = dep[u] - k;
        
        while (dep[top[u]] > d) {
            u = parent[top[u]];
        }
        
        return seq[in[u] - dep[u] + d];
    }
    
    bool isAncester(int u, int v) {
        return in[u] <= in[v] && in[v] < out[u];
    }
    
    int rootedParent(int u, int v) {
        std::swap(u, v);
        if (u == v) {
            return u;
        }
        if (!isAncester(u, v)) {
            return parent[u];
        }
        auto it = std::upper_bound(adj[u].begin(), adj[u].end(), v, [&](int x, int y) {
            return in[x] < in[y];
        }) - 1;
        return *it;
    }
    
    int rootedSize(int u, int v) {
        if (u == v) {
            return n;
        }
        if (!isAncester(v, u)) {
            return siz[v];
        }
        return n - siz[rootedParent(u, v)];
    }
    
    int rootedLca(int a, int b, int c) {
        return lca(a, b) ^ lca(b, c) ^ lca(c, a);
    }
};
```

## $\text{MinCostFlow}$

```cpp
template <class T>
struct MinCostFlow {
    struct Edge {
        int v, c;
        T f;
        Edge(int v, int c, T f) : v(v), c(c), f(f) {}
    };
    int n;
    vector<Edge> e;
    vector<vector<int>> g;
    vector<T> h, dis;
    vector<int> pre;
    vector<T> p;
    bool dijkstra(int s, int t) {
        dis.assign(n, numeric_limits<T>::max());
        pre.assign(n, -1);
        priority_queue<pair<T, int>, vector<pair<T, int>>, greater<>> que;
        dis[s] = 0;
        que.emplace(0, s);
        while (!que.empty()) {
            auto [d, u] = que.top();
            que.pop();
            if (dis[u] != d) continue;
            for (int &i : g[u]) {
                auto [v, c, f] = e[i];
                if (c > 0 && dis[v] > d + h[u] - h[v] + f) {
                    dis[v] = d + h[u] - h[v] + f;
                    pre[v] = i;
                    que.emplace(dis[v], v);
                }
            }
        }
        return dis[t] != numeric_limits<T>::max();
    }
    MinCostFlow(const int &n = 0) : n(n), g(n) {}
    void addEdge(int u, int v, int c, T f) {
        g[u].push_back(e.size());
        e.emplace_back(v, c, f);
        g[v].push_back(e.size());
        e.emplace_back(u, 0, -f);
    }
    pair<int, T> maxflow(int s, int t) {
        int flow = 0;
        T cost = 0;
        h.assign(n, 0);
        
        for (int i = 0; i < n; i++) {
            for (auto &j : g[i]) {
                if (e[j].v > i) {
                    h[e[j].v] = min(h[e[j].v], h[i] + e[j].f);
                }
            }
        }
        
        while (dijkstra(s, t)) {
            for (int i = 0; i < n; ++i) {
                h[i] += dis[i];
            }
            int aug = numeric_limits<int>::max();
            for (int i = t; i != s; i = e[pre[i] ^ 1].v) {
                aug = min(aug, e[pre[i]].c);
            }
            for (int i = t; i != s; i = e[pre[i] ^ 1].v) {
                e[pre[i]].c -= aug;
                e[pre[i] ^ 1].c += aug;
            }
            flow += aug;
            p.push_back(aug * h[t]);
            cost += aug * h[t];
        }
        return {flow, cost};
    }
};
```

## $\text{Inv}$

```cpp
void exgcd(ll a,ll b,ll& d,ll& x,ll& y)
{
    if(!b) { d = a; x = 1; y = 0; }
    else{ exgcd(b, a%b, d, y, x); y -= x*(a/b); }
}
 
ll inv(ll a, ll p)
{
    ll d, x, y;
    exgcd(a, p, d, x, y);
    return d == 1 ? (x+p)%p : -1;
}
```

## $\text{partial_sum}$

```cpp
vector<i64> sum(a.begin(), a.end());
partial_sum(sum.begin(), sum.end(), sum.begin());
```

## $\text{pr}$

```cpp
i64 mul(i64 a, i64 b, i64 m) {
    return static_cast<__int128>(a) * b % m;
}
i64 power(i64 a, i64 b, i64 m) {
    i64 res = 1 % m;
    for (; b; b >>= 1, a = mul(a, a, m))
        if (b & 1)
            res = mul(res, a, m);
    return res;
}
bool isprime(i64 n) {
    if (n < 2)
        return false;
    static constexpr int A[] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
    int s = __builtin_ctzll(n - 1);
    i64 d = (n - 1) >> s;
    for (auto a : A) {
        if (a == n)
            return true;
        i64 x = power(a, d, n);
        if (x == 1 || x == n - 1)
            continue;
        bool ok = false;
        for (int i = 0; i < s - 1; ++i) {
            x = mul(x, x, n);
            if (x == n - 1) {
                ok = true;
                break;
            }
        }
        if (!ok)
            return false;
    }
    return true;
}
std::vector<i64> factorize(i64 n) {
    std::vector<i64> p;
    std::function<void(i64)> f = [&](i64 n) {
        if (n <= 10000) {
            for (int i = 2; i * i <= n; ++i)
                for (; n % i == 0; n /= i)
                    p.push_back(i);
            if (n > 1)
                p.push_back(n);
            return;
        }
        if (isprime(n)) {
            p.push_back(n);
            return;
        }
        auto g = [&](i64 x) {
            return (mul(x, x, n) + 1) % n;
        };
        i64 x0 = 2;
        while (true) {
            i64 x = x0;
            i64 y = x0;
            i64 d = 1;
            i64 power = 1, lam = 0;
            i64 v = 1;
            while (d == 1) {
                y = g(y);
                ++lam;
                v = mul(v, std::abs(x - y), n);
                if (lam % 127 == 0) {
                    d = std::gcd(v, n);
                    v = 1;
                }
                if (power == lam) {
                    x = y;
                    power *= 2;
                    lam = 0;
                    d = std::gcd(v, n);
                    v = 1;
                }
            }
            if (d != n) {
                f(d);
                f(n / d);
                return;
            }
            +0;
        }
    };
    f(n);
    std::sort(p.begin(), p.end());
    return p;
}
```

## $\text{manacher}$

```cpp
// ---
string t = "#";
for (int i = 0; i < n; i++) {
    t += s[i];
    t += '#';
}

int m = t.size();
vector<int> r(m);
for (int i = 0, j = 0; i < m; i++) {
    if (2 * j - i >= 0) {
        r[i] = max(0, min(j + r[j] - i, r[2 * j - i]));
    }
    while (i - r[i] >= 0 && i + r[i] < m && t[i - r[i]] == t[i + r[i]]) {
        r[i]++;
    }
    if (i + r[i] > j + r[j]) {
        j = i;
    }
}
```

## $\text{match}$

匈牙利二分图匹配

```cpp
constexpr int N = 1 << 16;
 
void solve() {
  int n;
  cin >> n;
  int lg = __lg(n);
  vector<bool> used(n);
  vector<int> ans(n);
  for (int i = 0; i < n; i += 2) {
    cin >> ans[i];
    used[ans[i]] = 1;
  }
  vector<vector<int>> g(n);
  for (int i = 1; i < n; i += 2) {
    int x = ans[i - 1] ^ ans[i + 1];
    for (int j = 0; j < lg; j++) {
      if ((x >> j & 1) && !used[ans[i - 1] ^ (1 << j)]) {
        g[i].push_back(ans[i - 1] ^ (1 << j));
      }
    }
  }
  for (int i = 0; i < lg; i++) {
    if (!used[ans[n - 2] ^ (1 << i)]) {
      g[n - 1].push_back(ans[n - 2] ^ (1 << i));
    }
  }
  vector<int> match(n, -1);
  bitset<N> vis;
  function<bool(int)> dfs = [&](int pos) {
    for (auto num : g[pos]) {
      if (!vis[num]) {
        vis[num] = 1;
        if (match[num] == -1 || dfs(match[num])) {
          match[num] = pos;
          return true;
        }
      }
    }
    return false;
  };
  for (int i = 1; i < n; i += 2) {
    vis.reset();
    dfs(i);
  }
  for (int i = 0; i < n; i++) {
    if (!used[i]) {
      ans[match[i]] = i;
    }
  }
  for (int i = 0; i < n; i++) {
    cout << ans[i] << " \n"[i == n - 1];
  }
}
```

## $\text{java}$

```java
import java.io.*;
import java.util.*;
 
public class Main {
    static BufferedReader Input = new BufferedReader(new InputStreamReader(System.in));
    static PrintWriter Print = new PrintWriter(new OutputStreamWriter(System.out));
 
    static String input() throws IOException {
        return Input.readLine();
    }
 
    static void printf(String format, Object... args) {
        Print.printf(format, args);
    }
 
    static void print(Object x) {
        Print.println(x);
    }
 
    static void solve() throws IOException {
    }
 
 
    public static void main(String[] args) throws IOException {
        int t = Integer.parseInt(input());
 
        while (t-- > 0) {
            solve();
        }
 
        Print.flush();
    }
}
```

## $\text{one mod hash}$

```cpp
constexpr int B = 114514;
constexpr i64 P = 1000000000039;

i64 *p;

void init(int N) {
    p = new i64 [N + 1];
    for (int i = 0; i <= N; i++) {
        p[i] = 0;
    } 
    p[0] = 1;
    for (int i = 1; i <= N; i++) {
        p[i] = p[i - 1] * B % P;
    }
}

struct StringHash {
    vector<i64> h;
    StringHash() : h(1) {}
    void push_back(char ch) {
        h.push_back((h.back() * B + ch) % P);
    }
    i64 toi64(int l, int r) { // [l, r)
        return (h[r] + __int128(h[l]) * (P - p[r - l])) % P;
    }
};
/*
    const ll P1=1e14+31,P2=1e14+67;
    const int BASE1=777,BASE2=4396;
*/
```

## $\text{sg}$

打表大法

```cpp
constexpr int N = 1E4;
 
int f[N + 1];
 
void solve() {
  int k;
  cin >> k;
  vector<int> a(k);
  for (int i = 0; i < k; i++) {
    cin >> a[i];
  }
  sort(a.begin(), a.end());
  auto get = [&]() {
    for (int i = 0; i <= N; i++) {
      bool vis[N + 1];
      for (int j = 0; j <= N; j++) {
        vis[j] = 0;
      }
      for (int j = 0; j < k; j++) {
        if (i - a[j] >= 0) {
          vis[f[i - a[j]]] = 1;
        }
      }
      for (int j = 0; j <= N; j++) {
        if (!vis[j]) {
          f[i] = j;
          break;
        }
      }
    }
  };
  get();
  int m;
  cin >> m;
  for (int i = 0; i < m; i++) {
    int l;
    cin >> l;
    int ans = 0;
    for (int j = 0; j < l; j++) {
      int x;
      cin >> x;
      ans ^= f[x];
    }
    cout << (ans != 0 ? "W" : "L");
  }
}
```

## $\text{Math}$

![002](https://cdn.luogu.com.cn/upload/image_hosting/8e4osaum.png)

## $\text{phi}$

```cpp
vector<i64> phi;
vector<int> minp, primes;
 
void sieve(int n) {
    minp.assign(n + 1, 0);
    phi.assign(n + 1, 0);

    phi[1] = 1;
    for (int i = 2; i <= n; i++) {
        if (minp[i] == 0) {
            minp[i] = i;
            phi[i] = i - 1;
            primes.push_back(i);
        }
        
        for (auto p : primes) {
            if (i * p > n) {
                break;
            }
            minp[i * p] = p;
            if (i % p == 0) {
                phi[i * p] = phi[i] * p;
                break;                
            } else {
                phi[i * p] = phi[i] * phi[p];
            }
        }
    }
} 
```

## $\text{LCP}$

```cpp
void solve() {
    int n;
    std::cin >> n;
    
    std::vector<std::string> s(n);
    for (int i = 0; i < n; i++) {
        std::cin >> s[i];
    }
    
    int ans = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            int l1 = s[i].size(), l2 = s[j].size();
            std::vector f(l1 + 1, std::vector(l2 + 1, 0));
            for (int x = l1 - 1; x >= 0; x--) {
                for (int y = l2 - 1; y >= 0; y--) {
                    f[x][y] = s[i][x] == s[j][y] ? 1 + f[x + 1][y + 1] : 0;
                    ans = std::max(ans, f[x][y]);
                }
            }
        }
    }
    std::cout << ans << "\n";
}
```

