#include "simulator.h"
#include "tetgen/tetgen.h"
#include <iostream>
#include <queue>
#include <set>


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

void meshObject::createTetrahedralMesh(float scale)
{

    tetgenio input;

    // TetGen switches
    char TETGEN_SWITCHES[] = "pq1.414a0.0002";

    if (input.load_off("/home/agalvan-admin/CXPBD/resources/ducky/cube.off")) {
        // use TetGen to tetrahedralize our mesh
        tetgenio output;
        tetrahedralize(TETGEN_SWITCHES, &input, &output);

        this->nvertices = output.numberofpoints;
        Matrix<float,Dynamic,Dynamic,RowMajor> points(output.numberofpoints, 3);

        // create a vertex in the object for each point of the result
        for (int p = 0, pi = 0; p < output.numberofpoints; ++p, pi += 3)
        {
            points.row(p) = RowVector3f(output.pointlist[pi + 0],
                                               output.pointlist[pi + 1],
                                               output.pointlist[pi + 2]);
        }

        // sets the vertices of the mesh
        this->x = scale*points;

        Matrix<int,Dynamic,Dynamic,RowMajor> faces(output.numberoftrifaces, 3);
        this->nfaces = output.numberoftrifaces;

        // create a triangle for each face on the surface
        for (int t = 0, ti = 0; t < output.numberoftrifaces; ++t, ti += 3) {
            unsigned int vi[3];

            for (int i = 0; i < 3; ++i)
            {
                int tc = output.trifacelist[ti + i];
                vi[i] = tc;
                int pi = tc * 3;
            }
            //unsigned int index = a_object->newTriangle(p[1], p[0], p[2]);
            //a_chai3dMesh->newTriangle(vi[0], vi[1], vi[2]);
            faces.row(t) = RowVector3i(vi[0], vi[1], vi[2]);
        }

        // sets the faces of the mesh
        this->F = faces;

        this->nelements = output.numberoftetrahedra;
        Matrix<int,Dynamic,Dynamic,RowMajor> tetrahedra(output.numberoftetrahedra, 4);
        MatrixXd tetrahedra_centroids(output.numberoftetrahedra, 3);

        for (int t = 0, ti = 0; t < output.numberoftetrahedra; ++t, ti += 4) {

            int v0 = output.tetrahedronlist[ti + 0];
            int v1 = output.tetrahedronlist[ti + 1];
            int v2 = output.tetrahedronlist[ti + 2];
            int v3 = output.tetrahedronlist[ti + 3];

            RowVector4i tetrahedron;
            tetrahedron[0] = v0;
            tetrahedron[1] = v1;
            tetrahedron[2] = v2;
            tetrahedron[3] = v3;


            tetrahedra.row(t) = (tetrahedron);

        }

        // compute the neighbors for each tetrahedron (share faces)
        vector<set<int>> tetrahedra_neighbors[output.numberoftetrahedra];

        for (int i = 0; i < output.numberoftetrahedra; i++) {
            set<int> temp;
            temp.insert(tetrahedra(i, 0));
            temp.insert(tetrahedra(i, 1));
            temp.insert(tetrahedra(i, 2));
            temp.insert(tetrahedra(i, 3));

            for (int j = 0; j < output.numberoftetrahedra; j++) {
                if (i != j) {
                    temp.insert(tetrahedra(i, 0));
                    temp.insert(tetrahedra(i, 1));
                    temp.insert(tetrahedra(i, 2));
                    temp.insert(tetrahedra(i, 3));

                    if (temp.size() == 5) {
                        tetrahedra_neighbors->at(i).insert(j);
                    }
                }
            }

        }

        // get all the edges of our tetrahedra
        set<pair<int, int>> springs;

        for (int t = 0, ti = 0; t < output.numberoftetrahedra; ++t, ti += 4)
        {
            // store each edge of the tetrahedron as a pair of indices
            for (int i = 0; i < 4; ++i) {
                int v0 = output.tetrahedronlist[ti + i];
                for (int j = i + 1; j < 4; ++j)
                {
                    int v1 = output.tetrahedronlist[ti + j];
                    springs.insert(pair<int, int>(min(v0, v1), max(v0, v1)));
                }
            }
        }

        Matrix<int,Dynamic,Dynamic,RowMajor> edges(springs.size(),2);

        int i = 0;

        for (auto sp : springs)
        {
            edges.row(i) = Eigen::RowVector2i(sp.first, sp.second);
            i++;
        }

        this->E = edges;
        this->h_ei = this->E.data();
        this->nedges = edges.rows();

        // set the tetrahedra position
        this->T = tetrahedra;

        // set initial velocity to zero
        this->xdot = this->x;
        this->xdot.setZero();

        // create pointer arrays
        this->h_ti = this->T.data();
        this->h_x = this->x.data();
        this->h_xlast = this->x.data();
        this->h_x0 = this->x.data();
        this->h_xdot = this->xdot.data();
        this->h_fi = this->F.data();

        // create mass array
        this->h_m =  (float*)malloc(this->nvertices*sizeof (float));
        for (int i = 0; i < this->nvertices; i++)
        {
            this->h_m[i] = 0.001;
        }

        set<int> interior_set;
        // find the indices that are in the middle
        for (int i = 0; i < this->nvertices; i++)
        {
            if (abs(this->h_x[i,1]) < 0.1)
            {
                interior_set.insert(i);
            }
        }

        auto x_copy = x;
        x_copy.col(1) = -x_copy.col(1);

    }


}


void simulator::initializeVolumeConstraint(float alpha,float beta)
{
    // creates neohookean constraint
    vc = new VolumeConstraint;

    // set values
    vc->h_alpha = alpha;
    vc->h_beta = beta;

    // defines number of constraints
    vc->h_nconstraints = object->nelements;

    // allocates space
    vc->h_v0 = (float*)malloc(object->nelements*sizeof (float));

    for ( int i = 0 ; i < object->nelements; i++)
    {

        int v1 = object->T(i,0);
        int v2 = object->T(i,1);
        int v3 = object->T(i,2);
        int v4 = object->T(i,3);

        auto const p1 = object->x.row(v1);
        auto const p2 = object->x.row(v2);
        auto const p3 = object->x.row(v3);
        auto const p4 = object->x.row(v4);

        Matrix3f Dm;
        Dm.col(0) = (p1 - p4).transpose();
        Dm.col(1) = (p2 - p4).transpose();
        Dm.col(2) = (p3 - p4).transpose();

        vc->h_v0[i]     = abs((1. / 6.) * Dm.determinant());

    }
    vc->computeGraph(object->T);
    vc->transferToGPU();

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

void NeohookeanConstraint::computeGraph(MatrixXi E)
{

}

void VolumeConstraint::computeGraph(MatrixXi E)
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
    this->maxcolor = maxcolor;
    std::cout << maxcolor << std::endl;
    this->maxelem = maxpool;

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