#include <bits/stdc++.h>
using namespace std;

void solve()
{
    string s;
    int i = 0;
    cin >> s;
    for (auto c : s)
    {
        if (c == '+')
            i++;
        else
            i--;
    }
    cout << i << "\n";
    return;
}
int32_t main()
{
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    // int t;
    // cin >> t;
    // while (t--)
    solve();
    return 0;
}