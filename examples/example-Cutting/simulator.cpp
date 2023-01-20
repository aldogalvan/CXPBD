#include "simulator.h"
#include <iostream>
#include <queue>
#include <set>
#include <map>



// Function to count the minimum
// number of color required
void minimumColors(int N, int E,
                   int U[], int V[])
{

    // Create array of vectors
    // to make adjacency list
    vector<int> adj[N];

    // Initialise colors array to 1
    // and count array to 0
    vector<int> count(N, 0);
    vector<int> colors(N, 1);

    // Create adjacency list of
    // a graph
    for (int i = 0; i < N; i++) {
        adj[V[i] - 1].push_back(U[i] - 1);
        count[U[i] - 1]++;
    }

    // Declare queue Q
    queue<int> Q;

    // Traverse count[] and insert
    // in Q if count[i] = 0;
    for (int i = 0; i < N; i++) {
        if (count[i] == 0) {
            Q.push(i);
        }
    }

    // Traverse queue and update
    // the count of colors
    // adjacent node
    while (!Q.empty()) {
        int u = Q.front();
        Q.pop();

        // Traverse node u
        for (auto x : adj[u]) {
            count[x]--;

            // If count[x] = 0
            // insert in Q
            if (count[x] == 0) {
                Q.push(x);
            }

            // If colors of child
            // node is less than
            // parent node, update
            // the count by 1
            if (colors[x] <= colors[u]) {
                colors[x] = 1 + colors[u];
            }
        }
    }

    // Stores the minimumColors
    // requires to color the graph.
    int minColor = -1;

    // Find the maximum of colors[]
    for (int i = 0; i < N; i++) {
        minColor = max(minColor, colors[i]);
    }

    // Print the minimum no. of
    // colors required.
    std::cout << minColor << std::endl;
}

void meshObject::createClothMesh(float scale)
{

}



void simulator::initializeBendingConstraint(float alpha, float beta)
{

    std::map<std::pair<int, int>, std::vector<int>> edgemap;
    int nfaces = object->nfaces;
    for (int i = 0; i < nfaces; i++) {
        for (int j = 0; j < 3; j++) {
            int nextj = (j + 1) % 3;
            int v1 = object->F(i, j);
            int v2 = object->F(i, nextj);
            if (v1 > v2)
                std::swap(v1, v2);
            edgemap[std::pair<int, int>(v1, v2)].push_back(i);
        }
    }

    int nhinges = 0;
    for (auto it : edgemap) {
        if (it.second.size() == 2)
            nhinges++;
    }
    H_.resize(nhinges, 4);
    int idx = 0;
    for (auto it : edgemap) {
        if (it.second.size() != 2)
            continue;
        std::set<int> hingeverts;
        for (int j = 0; j < 3; j++) {
            hingeverts.insert(object->F(it.second[0], j));
            hingeverts.insert(object->F(it.second[1], j));
        }
        int colidx = 0;
        for (auto v : hingeverts) {
            H_(idx, colidx) = v;
            colidx++;
        }
        idx++;
    }

}

void simulator::initializeEdgeConstraint(float alpha, float beta)
{
    ec = new EdgeConstraint;

    // set values
    ec->h_alpha = alpha;
    ec->h_beta = beta;

    // defines number of constraints
    ec->h_nconstraints = object->nedges;

    // allocates space
    ec->h_e0 = (float*)malloc(object->nedges*sizeof (float));

    for ( int i = 0 ; i < object->nedges; i++)
    {
        int v1 = object->E(i,0);
        int v2 = object->E(i,1);

        auto const p1 = object->x.row(v1);
        auto const p2 = object->x.row(v2);

        ec->h_e0[i] = (p2 - p1).norm();

    }

    ec->computeGraph(object->E);
    ec->transferToGPU();
}

void simulator::initializeFixedPointConstraint(void)
{
    fpc = new FixedPointConstraint;

    float zmin = numeric_limits<float>::max();
    float ep = 1e-3;

    // find minimum
    for (int i = 0; i < object->x.rows(); i++)
    {
        if (object->x(i,2) < zmin)
        {
            zmin = object->x(i,2);
        }
    }

    fpc->h_nconstraints = 0;
    std::vector<int> id;

    // constrain to minimum
    for (int i = 0; i < object->x.rows(); i++)
    {
        if (object->x(i,2) <= zmin + ep)
        {
            fpc->h_nconstraints += 1;
            id.emplace_back(i);
        }
    }

    fpc->h_i = id.data();
    fpc->transferToGPU();

}

void BendingConstraint::computeGraph(MatrixXi E)
{

    // graph algorithm
    vector<set<int>> nlist(E.rows());

    for (int c = 0 ; c < E.rows() ; c++)
    {
        int c1 = E(c,0);
        int c2 = E(c,1);
        int c3 = E(c,2);
        int c4 = E(c,3);

        set<int> adj;

        for (int a = 0; a < E.rows(); a++)
        {
            if (a == c)
                continue;

            int a1 = E(a,0);
            int a2 = E(a,1);
            int a3 = E(a,2);
            int a4 = E(a,3);

            if (c1 == a1 || c1 == a2 || c1 == a3 || c1 == a4 ||
                c2 == a1 || c2 == a2 || c2 == a3 || c2 == a4 ||
                c3 == a1 || c3 == a2 || c3 == a3 || c3 == a4 ||
                c4 == a1 || c4 == a2 || c4 == a3 || c4 == a4 )
            {
                adj.insert(a);
            }
        }

        nlist[c] = adj;
    }

    // graph coloring
    vector<int> nodecolor(nlist.size());

    vector<set<int>> colorpool;
    colorpool.emplace_back(set<int>());

    std::fill(nodecolor.begin(), nodecolor.end(),-1);

    for (int v = nlist.size() - 1; v >= 0; v--)
    {

        if (nodecolor[v] == -1)
        {
            set<int> neighbor_colors;

            // check neighbors
            for (auto it : nlist[v])
            {

                if (nodecolor[it] != -1)
                    neighbor_colors.insert(nodecolor[it]);

            }

            if (neighbor_colors.size() == 0)
            {
                nodecolor[v] = 0;
            }
            else
            {
                int coloridx = 0;
                for (;;)
                {
                    if (!neighbor_colors.count(coloridx))
                    {
                        nodecolor[v] = coloridx;

                        if (colorpool.size() <= coloridx)
                            colorpool.emplace_back(set<int>());

                        colorpool[coloridx].insert(v);

                        break;
                    }
                    coloridx++;
                }
            }
        }
    }

    // find minimum number of colors
    int N = nlist.size();
    set<pair<int,int>> E_;
    for (int u = 0; u < nlist.size(); u++)
    {
        for (auto v : nlist[u])
        {
            E_.insert(pair<int,int>(u,v));
        }
    }

    int NE = E_.size();
    int U[NE]; int V[NE];
    int idx = 0;

    for (auto it : E_)
    {
        U[idx] = it.first; V[idx] = it.second;
        idx++;
    }

    //minimumColors(N,NE,U,V);

    // TODO: BETTER GRAPH COLORING ALGORITHM / VERY IMPORTANT

    // try another graph coloring algorithm
    /*
    std::fill(nodecolor.begin(), nodecolor.end(),-1);
    int nodeidx = 0;
    int coloridx = 0;


    for (;;)
    {
        // check if all nodes have been colored
        bool flag = 1;

        for (int u = 0 ; u < nlist.size() ; u++)
        {
            if (nodecolor[u] == -1)
            {
                flag = 0;
                break;
            }
        }

        if (flag == 1)
        {
            break;
        }

        if (nodecolor[nodeidx] == -1)
        {
            nodecolor[nodeidx] = coloridx;
        }

        flag = 0;
        for (auto v : nlist[nodeidx])
        {
            if (nodecolor[v] == nodeidx)
            {

            }
        }
    }
     */

    // get the max color
    int maxcolor = -1;
    for (auto it : nodecolor)
    {
        assert(it != -1);

        if (it > maxcolor)
        {
            maxcolor = it;
        }
    }

    // now find the maximum number of nodes in each color
    int maxpool = 0;
    for (auto it :colorpool)
    {
        if (it.size() > maxpool)
            maxpool = it.size();
    }

    Matrix<int,Dynamic,Dynamic,RowMajor> G(maxcolor,maxpool);
    G.fill(-1);

    for (int i = 0; i < maxcolor ; i++)
    {
        int ctr = 0;

        for (const auto it : colorpool[i])
        {
            G(i,ctr) = it;
            ctr++;
        }
    }

    this->h_graph = G.data();

}

void EdgeConstraint::computeGraph(MatrixXi E)
{

    // graph algorithm
    vector<set<int>> nlist(E.rows());

    for (int c = 0 ; c < E.rows() ; c++)
    {
        int c1 = E(c,0);
        int c2 = E(c,1);

        set<int> adj;

        for (int a = 0; a < E.rows(); a++)
        {
            if (a == c)
                continue;

            int a1 = E(a,0);
            int a2 = E(a,1);

            if (c1 == a1 || c1 == a2 || c2 == a1 || c2 == a2)
            {
                adj.insert(a);
            }
        }

        nlist[c] = adj;
    }

    // graph coloring
    vector<int> nodecolor(nlist.size());

    vector<set<int>> colorpool;
    colorpool.emplace_back(set<int>());

    std::fill(nodecolor.begin(), nodecolor.end(),-1);

    for (int v = 0; v < nlist.size(); v++)
    {

        if (nodecolor[v] == -1)
        {
            set<int> neighbor_colors;

            // check neighbors
            for (auto it : nlist[v])
            {

                if (nodecolor[it] != -1)
                    neighbor_colors.insert(nodecolor[it]);

            }

            if (neighbor_colors.size() == 0)
            {
                nodecolor[v] = 0;
            }
            else
            {
                int coloridx = 0;
                for (;;)
                {
                    if (!neighbor_colors.count(coloridx))
                    {
                        nodecolor[v] = coloridx;

                        if (colorpool.size() <= coloridx)
                            colorpool.emplace_back(set<int>());

                        colorpool[coloridx].insert(v);

                        break;
                    }
                    coloridx++;
                }
            }
        }
    }

    // find minimum number of colors
    int N = nlist.size();
    set<pair<int,int>> E_;
    for (int u = 0; u < nlist.size(); u++)
    {
        for (auto v : nlist[u])
        {
            E_.insert(pair<int,int>(u,v));
        }
    }

    int NE = E_.size();
    int U[NE]; int V[NE];
    int idx = 0;

    for (auto it : E_)
    {
        U[idx] = it.first; V[idx] = it.second;
        idx++;
    }

    //minimumColors(N,NE,U,V);

    // get the max color
    int maxcolor = -1;
    for (auto it : nodecolor)
    {
        assert(it != -1);

        if (it > maxcolor)
        {
            maxcolor = it;
        }
    }

    // now find the maximum number of nodes in each color
    int maxpool = 0;
    for (auto it :colorpool)
    {
        if (maxpool < it.size())
            maxpool = it.size();
    }

    Matrix<int,Dynamic,Dynamic,RowMajor> G(maxcolor,maxpool);
    G.fill(-1);

    for (int i = 0; i < maxcolor ; i++)
    {
        int ctr = 0;

        for (const auto it : colorpool[i])
        {
            // TODO: REMEMBER
            G(i,ctr) = it;
            ctr++;
        }

    }

    this->h_graph = G.data();
    this->maxcolor = maxcolor;
    std::cout << maxcolor << std::endl;
    this->maxelem = maxpool;

}