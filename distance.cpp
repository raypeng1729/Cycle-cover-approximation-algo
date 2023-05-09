#include <bits/stdc++.h>
using namespace std;

signed main(){
    freopen("test_data.txt", "r", stdin);
    int n; cin >> n;
    double x[n], y[n];
    for(int i = 0; i < n; i++) cin >> x[i] >> y[i];
    long double sum = 0;
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            sum += sqrtl((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));
        }
    }
    sum /= n * (n - 1) / 2;
    cout << sum;
}