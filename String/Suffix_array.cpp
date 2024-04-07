#include <bits/stdc++.h>
using namespace std;
int main()
{
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    string s;
    cin >> s;
    s += '$';
    vector<pair<char, int>> a(s.size());
    vector<int> c(s.size()), p(s.size());
    for (int i = 0; i < s.size(); i++)
        a[i] = {s[i], i};
    sort(a.begin(), a.end());
    for (int i = 0; i < s.size(); i++)
        p[i] = a[i].second;
    c[p[0]] = 0;
    for (int i = 1; i < s.size(); i++)
    {
        if (a[i].first == a[i - 1].first)
            c[p[i]] = c[p[i - 1]];
        else
            c[p[i]] = c[p[i - 1]] + 1;
    }
    int k = 0;
    while ((1 << k) <= s.size())
    {
        vector<pair<pair<int, int>, int>> v(s.size());
        for (int i = 0; i < s.size(); i++)
            v[i] = {{c[i], c[(i + (1 << k)) % s.size()]}, i};
        sort(v.begin(), v.end());
        for (int i = 0; i < s.size(); i++)
            p[i] = v[i].second;
        c[p[0]] = 0;
        for (int i = 1; i < s.size(); i++)
        {
            if (v[i].first == v[i - 1].first)
                c[p[i]] = c[p[i - 1]];
            else
                c[p[i]] = c[p[i - 1]] + 1;
        }
        k++;
    }
    for (auto i : p)
        cout << i << " ";
    cout << "\n";
    return 0;
}