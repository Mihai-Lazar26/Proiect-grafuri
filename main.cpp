#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <queue>
#include <stack>
#include <algorithm>
#include <climits>
#include <unordered_map>

using namespace std;

ifstream fin("date.in");
ofstream fout("date.out");

#define NMAX 100001

class Graf{
    vector<int> nod[NMAX];
    vector<pair<int, int>> nodCosturi[NMAX];
    int nrNoduri, nrMuchi;
    bool orientat = false, neorientat = false;

    vector<int> left[NMAX];
    vector<int> right[NMAX];
    int leftNoduri, rightNoduri;


    typedef pair<int, pair<int, int>> costMuchi;

    void actBFS(int S, vector<int> &viz);
    void actDFS(int S, vector<bool> &viz);
    void biconexeDFS(int S, int Tata, vector<int> &viz, vector<int> &nma, stack<int> &stiva,
                       vector<vector<int>> &componente);

    void ctcDFS(int S, vector<int> &viz, vector<int> &low, vector<bool> &inStruct, stack<int> &stiva,
                  vector<vector<int>> &componente, int &id);

    void dfsSt(int S, vector<bool> &viz, vector<int> &res);

    void dfsMC(int S, int Tata, vector<int> &viz, vector<int> &nma, vector<vector<int>> &componente);

    void actApmPrim(int S, vector<int> &viz, priority_queue<costMuchi, vector<costMuchi>, greater<costMuchi>> &pqMuchi,
                    vector<costMuchi> &results);

    void actDijkstra(int S, vector<int> &dist, vector<bool> &viz,
                     priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> &pqMuchi);

    int bfsMaxflow(int s, int d, vector<int>& parent, vector<vector<int>>& capacitate, vector<vector<int>>& flux,
                   vector<bool>& viz);

    bool dfsCuplaj(int nod, vector<int> &pairL, vector<int> &pairR, vector<bool> &visited);

public:

    void setOrientat(bool orientat = true){
        if(orientat) this->orientat = true, this->neorientat = false;
        else this->orientat = false, this->neorientat = true;
    }

    void setNeorientat(bool neorientat = true){
        if(neorientat) this->neorientat = true, this->orientat = false;
        else this->neorientat = false, this->orientat = true;
    }

    void citire(int nrNoduri, int nrMuchi, vector<pair<int, int>> muchi, bool orientat = true){
        if(orientat) setOrientat();
        else setNeorientat();
        this->nrNoduri = nrNoduri;
        this->nrMuchi = nrMuchi;
        for(auto i : muchi){
            this->nod[i.first].push_back(i.second);
            if(this->neorientat) this->nod[i.second].push_back(i.first);

        }
    }

    void citireCosturi(int nrNoduri, int nrMuchi, vector<pair<int, pair<int, int>>> muchi, bool orientat = true){
        if(orientat) setOrientat();
        else setNeorientat();
        this->nrNoduri = nrNoduri;
        this->nrMuchi = nrMuchi;
        for(auto i : muchi){
            int nodMe = i.second.first;
            int nodT = i.second.second;
            int cost = i.first;

            this->nodCosturi[nodMe].push_back(make_pair(nodT, cost));
            if(this->neorientat) this->nodCosturi[nodT].push_back(make_pair(nodMe, cost));
        }
    }

    void setNrNoduri(int n){
        this->nrNoduri = n;
    }

    void setNrMuchi(int n){
        this->nrMuchi = n;
    }

    void addMuchie(int x, int y){
        nod[x].push_back(y);
    }

    void addMuchieCost(int x, int y, int c){
        nodCosturi[x].push_back(make_pair(y, c));
    }

    void setLeftNoduri(int n){
        this->leftNoduri = n;
    }

    void setRightNoduri(int n){
        this->rightNoduri = n;
    }

    void addLRMuchie(int x, int y){
        left[x].push_back(y);
        right[y].push_back(x);
    }


    vector<int> BFS(int S);

    int DFS();

    vector<vector<int>> compBiconexe();

    vector<vector<int>> CTC();

    vector<int> ST();

    bool HavelHakimi(vector<int> grade);

    vector<vector<int>> muchiCritice();

    vector<pair<int, pair<int, int>>> apmPrim();

    int findParinte(int x, vector<int> &paduri);

    void disjoint1(int x, int y, vector<int> &paduri);

    bool disjoint2(int x, int y, vector<int> &paduri);

    vector<int> dijkstra();

    pair<vector<int>, bool> bellmanFord(int start);

    int maxflow(int s, int t, vector<vector<int>> &capacitate);

    vector<vector<int>> royfloyd(int n, vector<vector<int>> matAd);

    int darb();

    vector<int> cicluEuler();

    int hamilton(vector<vector<int>> &costuriMuchi);

    vector<pair<int, int>> hopcroftKarp();


};
Graf graf;


void Graf::actBFS(int S, vector<int> &viz){

    queue<int> coada;

    coada.push(S);
    viz[S] = 0;
    while(!coada.empty()){
        int current = coada.front();
        coada.pop();
        int dist = viz[current] + 1;
        for(auto vecin : nod[current]){
            if(viz[vecin] == -1){
                viz[vecin] = dist;
                coada.push(vecin);
            }
        }
    }
}

vector<int> Graf::BFS(int S){
    vector<int> viz(nrNoduri+1, -1);

    actBFS(S, viz);
    return viz;
}

void Graf::actDFS(int S, vector<bool> &viz){
    viz[S] = 1;
    for(auto vecini : nod[S]){
        if(!viz[vecini]) actDFS(vecini, viz);
    }
}

int Graf::DFS(){
    int componente = 0;

    vector<bool> viz(nrNoduri+1, 0);
    for(int i = 1; i <= nrNoduri; i++){
        if(!viz[i]) actDFS(i, viz), ++componente;
    }
    return componente;
}

void Graf::biconexeDFS(int S, int Tata, vector<int> &viz, vector<int> &nma, stack<int> &stiva,
                       vector<vector<int>> &componente){
    nma[S] = viz[S] = viz[Tata] + 1;
    stiva.push(S);
    for(auto vecin : nod[S]){
        if(vecin != Tata){
            if(viz[vecin]){
                nma[S] = min(nma[S], viz[vecin]);
            }
            else{
                biconexeDFS(vecin, S, viz, nma, stiva, componente);
                nma[S] = min(nma[S], nma[vecin]);

                if(viz[S] <= nma[vecin]){
                    vector<int> comp;
                    while(stiva.top() != vecin){
                        comp.push_back(stiva.top());
                        stiva.pop();
                    }
                    comp.push_back(vecin);
                    stiva.pop();
                    comp.push_back(S);
                    componente.push_back(comp);
                }
            }
        }
    }

}

vector<vector<int>> Graf::compBiconexe(){

    vector<int> viz(nrNoduri+1, 0);
    vector<int> nma(nrNoduri+1, 0);
    stack<int> stiva;
    vector<vector<int>> componente;
    biconexeDFS(1, 0, viz, nma, stiva, componente);
    return componente;
}

void Graf::ctcDFS(int S, vector<int> &viz, vector<int> &low, vector<bool> &inStruct, stack<int> &stiva,
                  vector<vector<int>> &componente, int &id){
    stiva.push(S);
    inStruct[S] = true;
    viz[S] = low[S] = id++;

    for(auto vecin : nod[S]){
        if(viz[vecin] == -1) ctcDFS(vecin, viz, low, inStruct, stiva, componente, id);
        if(inStruct[vecin]) low[S] = min(low[S], low[vecin]);
    }

    if(viz[S] == low[S]){
        vector<int> componenta;
        int node = stiva.top();
        while(node != S){
            inStruct[node] = false;
            componenta.push_back(node);
            stiva.pop();
            node = stiva.top();
        }
        componenta.push_back(node);
        inStruct[node] = false;
        stiva.pop();
        componente.push_back(componenta);
    }
}

vector<vector<int>> Graf::CTC(){

    vector<int> viz(nrNoduri+1, -1);
    vector<int> low(nrNoduri+1, 0);
    vector<bool> inStruct(nrNoduri+1, false);
    stack<int> stiva;
    vector<vector<int>> componente;
    int id = 1;

    for(int i = 1; i <= nrNoduri; i++){
        if(viz[i] == -1) ctcDFS(i, viz, low, inStruct, stiva, componente, id);
    }
    return componente;
}

void Graf::dfsSt(int S, vector<bool> &viz, vector<int> &res){
    viz[S] = 1;
    for(auto vecini : nod[S]){
        if(!viz[vecini]) dfsSt(vecini, viz, res);
    }
    res.push_back(S);
}

vector<int> Graf::ST(){
    vector<bool> viz(nrNoduri + 1, 0);
    vector<int> res;

    for(int i = 1; i <= nrNoduri; i++){
        if(!viz[i]) dfsSt(i, viz, res);
    }
    reverse(res.begin(), res.end());
    return res;
}

bool Graf::HavelHakimi(vector<int> grade){
    sort(grade.begin(), grade.end(), greater<int>());


    while(grade.size() && grade[0]){
        int prim = grade[0];
        grade.erase(grade.begin());

        if(prim > grade.size()) return false;

        for(auto& gr : grade){
            if(gr - 1 < 0)
                return false;
            else{
                --gr;
                --prim;
            }
            if(prim == 0) break;
        }
        if(prim) return false;

        sort(grade.begin(), grade.end(), greater<int>());

    }

    if(grade.size() && grade[0] == 0) return true;
    if(grade.size() == 0) return true;

    return false;
}

void Graf::dfsMC(int S, int Tata, vector<int> &viz, vector<int> &nma, vector<vector<int>> &componente){
    nma[S] = viz[S] = viz[Tata] + 1;
    for(auto vecin : nod[S]){
        if(vecin != Tata){
            if(viz[vecin]){
                nma[S] = min(nma[S], viz[vecin]);
            }
            else{
                dfsMC(vecin, S, viz, nma, componente);
                nma[S] = min(nma[S], nma[vecin]);

                if(viz[S] < nma[vecin])
                {
                    vector<int> temp;
                    temp.push_back(S);
                    temp.push_back(vecin);
                    componente.push_back(temp);
                }
            }
        }
    }
}

vector<vector<int>> Graf::muchiCritice(){

    vector<int> viz(nrNoduri+1, 0);
    vector<int> nma(nrNoduri+1, 0);
    vector<vector<int>> componente;
    dfsMC(1,0, viz, nma, componente);
    return componente;
}


void Graf::actApmPrim(int S, vector<int> &viz, priority_queue<costMuchi, vector<costMuchi>, greater<costMuchi>> &pqMuchi,
                      vector<costMuchi> &results){

    viz[S] = 1;
    for(auto vecini : nodCosturi[S]){
        int nod = vecini.first;
        int cost = vecini.second;
        if(!viz[nod])
            pqMuchi.push(make_pair(cost, make_pair(S, nod)));
    }

    while(!pqMuchi.empty()){
        auto top = pqMuchi.top();
        pqMuchi.pop();

        int target = top.second.second;

        if(!viz[target]){
            results.push_back(top);
            actApmPrim(target, viz, pqMuchi, results);
        }

    }

}


vector<pair<int, pair<int, int>>> Graf::apmPrim(){
    priority_queue<costMuchi, vector<costMuchi>, greater<costMuchi>> pqMuchi;
    vector<int> viz;
    viz.assign(nrNoduri+1, 0);
    vector<costMuchi> results;
    actApmPrim(1, viz, pqMuchi, results);



    return results;

}


int Graf::findParinte(int x, vector<int> &paduri){
    if(paduri[x] == x) return x;

    paduri[x] = findParinte(paduri[x], paduri);

    return paduri[x];

}

void Graf::disjoint1(int x, int y, vector<int> &paduri){
    int parx = findParinte(x, paduri);
    int pary = findParinte(y, paduri);

    if(parx != pary)
        paduri[pary] = parx;
}

bool Graf::disjoint2(int x, int y, vector<int> &paduri){

    int parx = findParinte(x, paduri);
    int pary = findParinte(y, paduri);

    return (parx == pary);
}

void Graf::actDijkstra(int S, vector<int> &dist, vector<bool> &viz,
                     priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> &pqMuchi){

    dist[S] = 0;
    pqMuchi.push(make_pair(0, S));

    while(!pqMuchi.empty()){
        auto top = pqMuchi.top();
        pqMuchi.pop();

        int target = top.second;
        if(!viz[target]){
            viz[target] = true;

            for(auto vecini : nodCosturi[target]){
                int nod = vecini.first;
                int cost = vecini.second;

                if(dist[target] + cost < dist[nod]){
                    dist[nod] = dist[target] + cost;
                    pqMuchi.push(make_pair(dist[nod], nod));
                }
            }
        }


    }

}


vector<int> Graf::dijkstra(){

    vector<int> dist;
    dist.assign(nrNoduri+1, INT_MAX);

    vector<bool> viz;
    viz.assign(nrNoduri+1, false);

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pqMuchi;

    actDijkstra(1, dist, viz, pqMuchi);

    for(int i = 2; i<=nrNoduri; i++)
        if(dist[i] == INT_MAX) dist[i] = 0;

    return dist;
}

pair<vector<int>, bool> Graf::bellmanFord(int start){

    vector<int> dist;
    dist.assign(nrNoduri + 1, INT_MAX);

    vector<bool> inCoada;
    inCoada.assign(nrNoduri+1, false);

    vector<int> cnt;
    cnt.assign(nrNoduri+1, 0);

    dist[start] = 0;

    queue<int> coada;

    coada.push(start);
    inCoada[start] = true;

    while(!coada.empty()){
        int nod = coada.front();
        coada.pop();
        inCoada[nod] = false;
        cnt[nod]++;
        if(cnt[nod] >= nrNoduri){
            return make_pair(dist, true);
        }

        for(auto vecini : nodCosturi[nod]){
            int vec = vecini.first;
            int cost = vecini.second;

            if(dist[nod]+cost < dist[vec]){
                dist[vec] = dist[nod] + cost;
                if(!inCoada[vec]){
                    coada.push(vec);
                    inCoada[vec] = true;
                }
            }
        }
    }


    return make_pair(dist, false);

}



int Graf::bfsMaxflow(int s, int d, vector<int>& parent, vector<vector<int>>& capacitate, vector<vector<int>>& flux,
                     vector<bool>& viz){


    for(int i = 1; i<=nrNoduri; i++){
        viz[i] = false;
    }
    viz[s] = true;
    queue<int> q;
    q.push(s);

    while(!q.empty()){
        int act = q.front();
        q.pop();
        if(act != d)
        for(int next : nod[act]){
            if(capacitate[act][next] == flux[act][next] || viz[next])
                continue;

            viz[next] = true;
            q.push(next);
            parent[next] = act;
        }
    }

    return viz[d];
}


int Graf::maxflow(int s, int d, vector<vector<int>> &capacitate){
    vector<int> parent(nrNoduri+1, 0);
    vector<vector<int>> flux(nrNoduri+1, vector<int>(nrNoduri+1, 0));
    vector<bool> viz(nrNoduri+1, false);
    int flow;

    for(flow = 0; bfsMaxflow(s, d, parent, capacitate, flux, viz);){
        for(int act : nod[d]){
            if(flux[act][d] == capacitate[act][d] || viz[act] == false)
                continue;

            parent[d] = act;
            int fmin = INT_MAX;
            for(int nod = d; nod != s; nod = parent[nod]){
                fmin = min(fmin, capacitate[parent[nod]][nod] - flux[parent[nod]][nod]);
            }

            if(fmin == 0) continue;

            for(int nod = d; nod != s; nod = parent[nod]){
                flux[parent[nod]][nod] += fmin;
                flux[nod][parent[nod]] -= fmin;
            }
            flow += fmin;
        }
    }

    return flow;
}

vector<vector<int>> Graf::royfloyd(int n, vector<vector<int>> dist){

    for(int k = 1; k <= n; k++){
        for(int i =1; i <=n; i++){
            for(int j = 1; j<= n; j++){
                if(dist[i][k] + dist[k][j] < dist[i][j] && dist[i][k] != INT_MAX && dist[k][j] != INT_MAX){
                    dist[i][j] = dist[i][k] + dist[k][j];
                }
            }
        }
    }

    return dist;
}

int Graf::darb(){
    queue<int> q;
    vector<int> dist(nrNoduri+1, 0);
    vector<bool> viz(nrNoduri+1, false);

    viz[1]= true;
    q.push(1);
    int last;

    while(!q.empty()){
        int cur = q.front();
        q.pop();

        last = cur;

        for(auto vec : nod[cur]){
            if(!viz[vec]){
                viz[vec] = true;
                q.push(vec);
            }

        }
    }

    replace(viz.begin(), viz.end(), true, false);

    q.push(last);
    viz[last] = true;
    dist[last] = 1;

    int maxi = 1;

    while(!q.empty()){
        int cur = q.front();
        q.pop();

        last = cur;

        for(auto vec : nod[cur]){
            if(!viz[vec]){
                dist[vec] = dist[cur] + 1;
                maxi = max(maxi, dist[vec]);
                viz[vec] = true;
                q.push(vec);
            }

        }
    }

    return maxi;
}


vector<int> Graf::cicluEuler(){
    vector<int> path;
    vector<bool> vizitat(nrMuchi, false);
    vector<vector<int>> muchiRamase(nrNoduri+1);
    vector<pair<int, int>> muchi;

    int nr = 0;

    for(int i = 1; i <= nrNoduri; i++)
        for(int j = 0; j < nod[i].size(); j++){
            muchiRamase[i].push_back(nr);
            muchiRamase[nod[i][j]].push_back(nr);
            muchi.push_back(make_pair(i, nod[i][j]));
            nr++;
        }

    for(int i = 1; i <= nrNoduri; i++)
        if(muchiRamase[i].size()%2 != 0)
            return path;

    stack<int> st;
    st.push(1);

    while(!st.empty()){
        int at = st.top();
        if(!muchiRamase[at].empty()){
            int muchieCurenta = muchiRamase[at].back();
            muchiRamase[at].pop_back();

            if(!vizitat[muchieCurenta]){
                vizitat[muchieCurenta] = true;
                int to;
                if(at == muchi[muchieCurenta].first)
                    to = muchi[muchieCurenta].second;
                else
                    to = muchi[muchieCurenta].first;
                st.push(to);
            }
        }
        else{
            st.pop();
            path.push_back(at);
        }
    }
    path.pop_back();
    return path;
}

int Graf::hamilton(vector<vector<int>> &costuriMuchi){
    int inf = 100000000;
    vector<vector<int>> best(nrNoduri, vector<int>(1 << nrNoduri, inf));

    best[0][1] = 0;
    for(int bstate = 1; bstate < (1 << nrNoduri) - 1; ++bstate){
        for(int n = 0; n < nrNoduri; ++n){
            if(best[n][bstate] == inf) continue;
            for(auto vec : nodCosturi[n]){
                int vecin = vec.first;
                int cost = vec.second;

                if((1 << vecin) & bstate) continue;

                int nstate = (1 << vecin) | bstate;

                if(best[vecin][nstate] > best[n][bstate] + cost){
                    best[vecin][nstate] = best[n][bstate] + cost;
                }
            }
        }
    }

    int mini = inf, sf = (1 << nrNoduri) - 1;
    for(int i = 1; i < nrNoduri; i++){
        if(costuriMuchi[i][0]!=inf && best[i][sf] != inf){
            mini = min(mini, best[i][sf] + costuriMuchi[i][0]);
        }
    }

    return mini;
}


bool Graf::dfsCuplaj(int nod, vector<int> &pairL, vector<int> &pairR, vector<bool> &visited){

    if(visited[nod]){
        return false;
    }

    visited[nod] = true;

    for(auto vec : left[nod]){
        if(pairR[vec] == -1){
            pairL[nod] = vec;
            pairR[vec] = nod;

            return true;
        }
    }

    for(auto vec : left[nod]){
        if(dfsCuplaj(pairR[vec], pairL, pairR, visited)){
            pairR[vec] = nod;
            pairL[nod] = vec;

            return true;
        }
    }

    return false;
}

vector<pair<int, int>> Graf::hopcroftKarp(){
    vector<pair<int, int>> M;

    vector<int> pairL(leftNoduri+1, -1);
    vector<int> pairR(rightNoduri+1, -1);
    vector<int> distance(leftNoduri+1, 0);
    vector<bool> visited;

    int ok = 1;

    while(ok){
        ok = 0;
        visited.assign(leftNoduri+1, false);

        for(int i = 1; i <= leftNoduri; ++i){
            if(pairL[i] == -1){
                if(dfsCuplaj(i, pairL, pairR, visited)){
                    ok = 1;
                }
            }

        }
    }

    for(int i = 1; i <= leftNoduri; ++i){
        if(pairL[i] != -1){
            M.push_back(make_pair(i, pairL[i]));
        }
    }

    return M;
}




void rezolvareBFS(){
    vector<pair<int, int>> muchi;
    int n, m, s;
    fin>>n>>m>>s;
    for(int i = 1; i <= m; i++){
        int x, y;
        fin>>x>>y;
        muchi.push_back(make_pair(x, y));
    }
    graf.citire(n, m, muchi);
    vector<int> sol = graf.BFS(s);
    for(int i = 1; i <= n; i++){
        fout<<sol[i]<<" ";
    }

}

void rezolvareDFS(){
    vector<pair<int, int>> muchi;
    int n, m;
    fin>>n>>m;
    for(int i = 1; i <= m; i++){
        int x, y;
        fin>>x>>y;
        muchi.push_back(make_pair(x, y));
    }
    graf.citire(n, m, muchi, false);

    int res = graf.DFS();
    fout<<res;
}

void rezolvareBiconexe(){
    vector<pair<int, int>> muchi;
    int n, m;
    fin>>n>>m;
    for(int i = 1; i <= m; i++){
        int x, y;
        fin>>x>>y;
        muchi.push_back(make_pair(x, y));
    }
    graf.citire(n, m, muchi, false);

    vector<vector<int>> sol = graf.compBiconexe();
    int nrComp = sol.size();
    fout<<nrComp<<'\n';
    for(int i = 0; i < nrComp; i++)
    {
        for(auto nod : sol[i])
            fout<<nod<<' ';
        fout<<'\n';
    }
}

void rezolvareCTC(){
    vector<pair<int, int>> muchi;
    int n, m;
    fin>>n>>m;
    for(int i = 1; i <= m; i++){
        int x, y;
        fin>>x>>y;
        muchi.push_back(make_pair(x, y));
    }
    graf.citire(n, m, muchi);
    vector<vector<int>> sol = graf.CTC();
    fout<<sol.size()<<'\n';
    for(auto comp : sol){
        for(auto nod : comp)
            fout<<nod<<' ';
        fout<<'\n';
    }
}

void rezolvareST(){
    vector<pair<int, int>> muchi;
    int n, m;
    fin>>n>>m;
    for(int i = 1; i <= m; i++){
        int x, y;
        fin>>x>>y;
        muchi.push_back(make_pair(x, y));
    }
    graf.citire(n, m, muchi);

    vector<int> sol = graf.ST();
    for(auto nod : sol)
        fout<<nod<<' ';
}

void rezolvareHH(){
    int n;
    vector<int> grade;
    fin>>n;
    for(int i = 1; i <= n; i++)
    {
        int x;
        fin>>x;
        grade.push_back(x);
    }

    bool sol = graf.HavelHakimi(grade);

    if(sol) fout<<"Da";
    else fout<<"Nu";
}

void rezolvareMC(){
    vector<pair<int, int>> muchi;
    int n, m;
    fin>>n>>m;
    for(int i = 1; i <= m; i++){
        int x, y;
        fin>>x>>y;
        muchi.push_back(make_pair(x, y));
    }
    graf.citire(n, m, muchi, false);
    vector<vector<int>> sol = graf.muchiCritice();
    for(auto muchie : sol){
        for(auto nod : muchie)
            fout<<nod<<" ";
        fout<<"\n";
    }
}

void rezolvareAPM(){
    vector<pair<int, pair<int, int>>> muchi;
    int n, m;
    fin>>n>>m;
    for(int i = 1; i<= m; i++){
        int x, y, cost;
        fin>>x>>y>>cost;

        muchi.push_back(make_pair(cost, make_pair(x, y)));
    }

    graf.citireCosturi(n, m, muchi, false);
    vector<pair<int, pair<int, int>>> results = graf.apmPrim();

    int sum = 0;
    for(auto res : results){
        sum+= res.first;
    }
    fout<<sum<<"\n"<<results.size()<<"\n";

    for(auto res : results){
        fout<<res.second.first<<" "<<res.second.second<<"\n";
    }

}

void rezolvareDisjoint(){
    int n, m;
    fin>>n>>m;
    vector<int> paduri;
    paduri.assign(n+1, 0);
    for(int i = 1; i<=n; i++)
        paduri[i] = i;

    for(int i = 1; i <= m; i++){
        int cod, x, y;
        fin>>cod>>x>>y;

        if(cod == 1)
            graf.disjoint1(x, y, paduri);

        if(cod == 2){
            bool res = graf.disjoint2(x, y, paduri);
            if(res) fout<<"DA\n";
            else fout<<"NU\n";
        }

    }
}


void rezolvareDijkstra(){
    int n, m;
    vector<pair<int, pair<int, int>>> muchi;

    fin>>n>>m;

    for(int i = 1; i <= m; i++){
        int a,b,c;
        fin>>a>>b>>c;

        muchi.push_back(make_pair(c, make_pair(a,b)));
    }

    graf.citireCosturi(n, m, muchi);

    muchi.clear();

    auto res = graf.dijkstra();

    for(int i = 2; i<=n; i++)
        fout<<res[i]<<" ";
}

void rezolvareBellmanFord(){
    int n, m;
    vector<pair<int, pair<int, int>>> muchi;

    fin>>n>>m;

    for(int i = 1; i <= m; i++){
        int a,b,c;
        fin>>a>>b>>c;

        muchi.push_back(make_pair(c, make_pair(a,b)));
    }

    graf.citireCosturi(n, m, muchi);

    muchi.clear();

    auto results = graf.bellmanFord(1);

    if(results.second) fout<<"Ciclu negativ!\n";
    else{
        for(int i = 2; i<=n; i++)
            fout<<results.first[i]<<" ";
    }


}

void rezolvareMaxFlow(){

    int n, m;
    fin>>n>>m;

    graf.setNrNoduri(n);
    graf.setNrMuchi(m);

    vector<vector<int>> capacitati(n+1, vector<int>(n+1, 0));
    for(int i = 1; i <= m; ++i){
        int x, y, z;
        fin>>x>>y>>z;
        graf.addMuchie(x, y);
        graf.addMuchie(y, x);
        capacitati[x][y] = z;
    }

    int res = graf.maxflow(1, n, capacitati);
    fout<<res;
}


void rezolvareRoyFloyd(){
    int n;
    vector<vector<int>> matAd(NMAX, vector<int>(NMAX, 0));
    fin>>n;
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            fin>>matAd[i][j];
            if(i != j && matAd[i][j] == 0)
                matAd[i][j] = INT_MAX;
        }
    }

    auto dist = graf.royfloyd(n, matAd);

    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= n; j++){
            if(dist[i][j] == INT_MAX) fout<<0<<" ";
            else fout<<dist[i][j]<<" ";
        }
        fout<<"\n";
    }

}

void rezolvareDarb(){
    int n;
    vector<pair<int, int>> muchi;
    fin>>n;
    for(int i = 1; i < n; i++){
        int x, y;
        fin>>x>>y;
        muchi.push_back(make_pair(x, y));
    }

    graf.citire(n, n-1, muchi, false);

    int sol = graf.darb();
    fout<<sol;
}

void rezolvareCicluEuler(){
    int n, m;
    fin>>n>>m;
    graf.setNrNoduri(n);
    graf.setNrMuchi(m);
    graf.setNeorientat();

    for(int i = 1; i <= m; i++){
        int x, y;
        fin>>x>>y;
        graf.addMuchie(x, y);
    }

    auto res = graf.cicluEuler();
    if(!res.empty()){
        for(auto it : res){
            fout<<it<<" ";
        }
    }
    else fout<<-1;
}

void rezolvareHamilton(){
    int inf = 100000000;
    int n, m;
    fin>>n>>m;
    vector<vector<int>> costuriMuchi(n, vector<int>(n, inf));
    graf.setNrMuchi(m);
    graf.setNrNoduri(n);
    graf.setOrientat();

    for(int i = 1; i <= m; i++){
        int x, y, c;
        fin>>x>>y>>c;
        graf.addMuchieCost(x, y, c);
        costuriMuchi[x][y] = c;
    }

    int res = graf.hamilton(costuriMuchi);
    if(res != inf) fout<<res;
    else fout<<"Nu exista solutie";
}

void rezolvareCuplaj(){
    int n, m, e;
    fin>>n>>m>>e;
    graf.setLeftNoduri(n);
    graf.setRightNoduri(m);
    graf.setNrMuchi(e);

    for(int i = 1; i <= e; ++i){
        int l, r;
        fin>>l>>r;
        graf.addLRMuchie(l, r);
    }

    auto res = graf.hopcroftKarp();
    fout<<res.size()<<"\n";
    for(auto val : res){
        fout<<val.first<<" "<<val.second<<"\n";
    }
}


int main(){
    
}
