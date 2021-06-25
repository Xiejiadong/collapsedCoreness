#include <bits/stdc++.h>
#include "glist.h"

#define PII pair<int,int>
#define inf 0x3fffffff
using namespace std;
const int maxn=10000,maxm=100000;
int n,m,ppp[maxn],pre[maxn],followers[maxn],collapsed[maxn],update[maxn],critical[maxn];
struct nodes{
    int coreness,degree,binpos,valid;
}node[maxn],tmp[maxn],save[maxn],tmp2[maxn];
struct Sort{
    int l,r;
}a[maxn];
vector<int> edge[maxn];
struct Edge{
    vector<int> More,Equal,Less;
}edges[maxn];
int budget;
vector<int> bin[maxn],hier,tid_;
int bb[maxn],cur_k;
int lowerbound[maxn],upperbound[maxn];
struct LCPSTreeVertex{
	vector<int> delta;
	vector<LCPSTreeVertex*> son;
	LCPSTreeVertex *p;
	int id;
	short k;
};
vector<LCPSTreeVertex*> t_;

void ReadEdgesS(int n, int m)
{
    srand(time(NULL));
    int v1,v2,k=0;
    for (int i=1;i<=m;i++)
    {
        scanf("%d%d", &v1, &v2);
        if (v1 > v2) swap(v1, v2);
        const auto e = make_pair(v1, v2);
        edge[v1].push_back(v2);
        edge[v2].push_back(v1);
        node[v1].degree++,node[v2].degree++;
    }
    for (int i=0;i<=n;i++)
        save[i]=node[i];
}

void __adjust_k(int target_k)
{
	if (cur_k < target_k) {
		int times = target_k - cur_k;
		hier.push_back(times * -2);
	}
	else if (cur_k > target_k) {
		int times = cur_k - target_k;
		hier.push_back(times * -2 + 1);
	}
	cur_k = target_k;
}

pair<char, int> decode_adjust(int mask)
{
	mask = -mask;
	if (mask & 1) {
		return make_pair(']', mask / 2 + 1);
	}
	else {
		return make_pair('[', mask / 2);
	}
}

void get_v(LCPSTreeVertex * rt, vector<int> & res)
{
	for (int v : rt->delta) {
		res.push_back(v);
	}
	for (auto p : rt->son) {
		get_v(p, res);
	}
}

vector<int> get_v(int v)
{
	vector<int> res;
	get_v(t_[tid_[v]], res);
	return res;
}

void __build_hier(vector<vector<int>> graph,vector<int> core,int n,int max_coreness)
{
	int kmax = max_coreness;
	vector<vector<int>*> R;
	R.resize(kmax + 1);
	for (int i = 0; i <= kmax; i++) {
		R[i] = new vector<int>;
	}
	vector<bool> vis(n, false);
	cur_k = -1;

	for (int i = 0; i < n; i++) {
		int cur_v = i;
		if (vis[cur_v]) continue;

		int cur_rmax = 0;
		R[0]->push_back(cur_v);
		while (true) {
			while (true) {
				if (cur_rmax < 0) break;
				if (R[cur_rmax]->empty()) {
					cur_rmax--;
					continue;
				}
				cur_v = R[cur_rmax]->back();
				R[cur_rmax]->pop_back();
				if (!vis[cur_v]) break;
			}
			if (cur_rmax < 0) break;
			
			if (cur_k > cur_rmax) __adjust_k(cur_rmax);
			int cv = core[cur_v];
			if (cv > cur_rmax) __adjust_k(cv);
			hier.push_back(cur_v);
			vis[cur_v] = true;
            for (auto nbr : graph[cur_v])
            {
				if (vis[nbr]) continue;
				int mn = min(core[cur_v], core[nbr]);
				cur_rmax = max(mn, cur_rmax);
				R[mn]->push_back(nbr);
			}
		}
		__adjust_k(-1);
	}
	for (int i = 0; i <= kmax; i++) {
		delete R[i];
	}
}

void __build_tree()
{
	LCPSTreeVertex *root = new LCPSTreeVertex;
	t_.push_back(root);
	LCPSTreeVertex *cur_tv = root;
	int cur_k = -1;
	for (int v : hier) {
		if (v < 0) {
			char op;
			int times;
			tie(op, times) = decode_adjust(v);
			if (op == '[') {
				while (times--) {
					LCPSTreeVertex *nxt_tv = new LCPSTreeVertex;
					nxt_tv->p = cur_tv;
					cur_tv->son.push_back(nxt_tv);
					cur_tv = nxt_tv;
					t_.push_back(nxt_tv);
				}
			}
			else if (op == ']') {
				while (times--) cur_tv = cur_tv->p;
			}
		}
		else {
			cur_tv->delta.push_back(v);
		}
	}
	hier.clear();
}

LCPSTreeVertex* dfs(LCPSTreeVertex* cur, vector<LCPSTreeVertex*> &store, LCPSTreeVertex* p)
{
	if (cur->delta.empty()) {
		auto pt = cur->son[0];
		delete cur;
		return dfs(pt, store, p);
	}
	else {
		cur->p = p;
		store.push_back(cur);
		for (int i = 0; i < cur->son.size(); i++) {
			cur->son[i] = dfs(cur->son[i], store, cur);
		}
		return cur;
	}
}

void __compress_tree(int n)
{
	vector<LCPSTreeVertex*> new_tv;
	LCPSTreeVertex *root = t_[0];
	for (int i = 0; i < root->son.size(); i++) {
		root->son[i] = dfs(root->son[i], new_tv, root);
		root->son[i]->p = NULL;
	}
	tid_ = vector<int>(n, 0);
	for (int i = 0; i < new_tv.size(); i++) {
		new_tv[i]->id = i;
		for (auto v : new_tv[i]->delta) {
			tid_[v] = i;
		}
	}
	t_ = new_tv;
	delete root;
}

void pre_LCPS(vector<vector<int>> graph,vector<int> core,int n,int max_coreness)
{
    hier.clear();cur_k=0;
    t_.clear();tid_.clear();
    __build_hier(graph,core,n,max_coreness);
	__build_tree();
	__compress_tree(n);
}


struct CCT{
    int id,layer;
}tree[maxn];
int nodeCCT,CCT_size[maxn],CCT_ance_size[maxn],CCT_child_size[maxn],CCT_fa[maxn];
vector<int> CCTchild[maxn];

void BuildCCT(int n,vector<int> G,vector<int> core,int fa,int Min_core)
{
    nodeCCT++;
    //cout<<nodeCCT<<endl;
    CCT_fa[nodeCCT]=fa;
    CCTchild[nodeCCT].clear();
    CCT_size[nodeCCT]=0;
    CCT_child_size[nodeCCT]=0;
    CCT_ance_size[nodeCCT]=CCT_ance_size[fa]+CCT_size[fa];
    for (auto x:G)
        if (Min_core==core[x])
            tree[x].id=nodeCCT,CCT_size[nodeCCT]++;
        else break;
    for (auto x:G)
        if (tree[x].id==0)
        {
            vector<vector<int> > bin(n+1);
            vector<int> nodes;
            nodes=get_v(x);
            int Min_coreness=core[x];
            for (auto y:nodes)
                bin[core[y]].push_back(y),
                Min_coreness=min(Min_coreness,core[y]);
            nodes.clear();
            for (int i=Min_coreness;i<=n;i++)
                    for (auto y:bin[i])
                        nodes.push_back(y);
            if (nodes.size()>0)
            {
                int id=nodeCCT,child=nodeCCT+1;
                CCTchild[id].push_back(child);
                BuildCCT(n,nodes,core,nodeCCT,Min_coreness);
                CCT_child_size[id]+=CCT_child_size[child]+CCT_size[child];
            }
        }
}
int ub(int x)
{
    x=tree[x].id;
    return CCT_ance_size[x]+CCT_size[x];
}
int upper_bound(int n,int x,vector<int> core)
{
    unordered_set<int> s;
    int sum=0;
    for (auto y:edges[x].Equal)
        if (!collapsed[y]&&critical[y])
        {
            if (s.find(tree[y].id)==s.end()) sum+=CCT_size[tree[y].id],s.insert(tree[y].id);
        }
    for (auto y:edges[x].Less)
        if (!collapsed[y]&&critical[y])
        {
            if (s.find(tree[y].id)==s.end()) sum+=CCT_size[tree[y].id],s.insert(tree[y].id);
        }
    return sum;
}
int lower_bound(int x)
{
    int sum=0;
    for (auto y:edges[x].Equal)
        if (!collapsed[y]&&critical[y]) sum++;
    for (auto y:edges[x].Less)
        if (!collapsed[y]&&critical[y]) sum++;
    return sum;
}

int main(int argc, char** argv)
{
    if (argc<3)
    {
        printf("Usage: ./core data budge\n");
        return 0;
    }
    int n, m, m2;
    char input[100],output[100];
    budget=atoi(argv[2]);
    sprintf(input,"%s.txt",argv[1]);
    sprintf(output,"%s-%d-GCC.txt",argv[1],budget);
    freopen(input,"r",stdin);
    freopen(output,"w",stdout);
    scanf("%d%d",&n,&m);
    // read the graph
    ReadEdgesS(n, m);
    // initialize the core component
    core::CoreMaintenance* cm = nullptr;
    cm = new core::GLIST(n);
    // create the adjacent list representation
    vector<vector<int>> graph(n);
    for (int i=0;i<n;i++)
        graph[i]=edge[i];
    vector<int> core(n);
    cm->ComputeCore(graph, true, core);
    for (int i=0;i<n;i++)
        ppp[i]=core[i],pre[i]=core[i];
    auto Starttime=(double)clock();
    int Max_follower=0,Max_follower_pos=0;
    for (int x=0;x<n;x++)
    {
        edges[x].More.clear();
        edges[x].Less.clear();
        edges[x].Equal.clear();
        for  (auto y:graph[x])
            if (core[y]>core[x]) edges[x].More.push_back(y);
            else if (core[y]<core[x]) edges[x].Less.push_back(y);
            else edges[x].Equal.push_back(y);
    }
    for (int i=0;i<n;i++)
    {
        int x=i;
        for (auto y:edges[x].Less)
        {
            if (x<y) cm->Remove(x, y, graph, core);
            else cm->Remove(y, x, graph, core);
        }
        for (auto y:edges[x].Equal)
        {
            if (x<y) cm->Remove(x, y, graph, core);
            else cm->Remove(y, x, graph, core);
        }
        followers[i]=0;
        for (int j=0;j<n;j++)
            if (i!=j&&core[j]<ppp[j]) followers[i]++;
        for (auto y:edges[x].Less)
        {
            if (x<y) cm->Insert(x, y, graph, core);
            else cm->Insert(y, x, graph, core);
        }
        for (auto y:edges[x].Equal)
        {
            if (x<y) cm->Insert(x, y, graph, core);
            else cm->Insert(y, x, graph, core);
        }
        collapsed[i]=0;
        update[i]=0;
        if (followers[i]>Max_follower) Max_follower=followers[i],Max_follower_pos=i;
    }
    auto Endtime=(double)clock();
    int ans=0;
    for (int i=1;i<=budget;i++)
    {
        printf("%d\t%d\t",i,ans+Max_follower);
        ans+=Max_follower;
        Endtime=(double)clock();
        printf("%0.5f\n",(Endtime-Starttime)/CLOCKS_PER_SEC);
        //collapse node
        int x=Max_follower_pos;
        for (auto y:edges[x].More)
            if(!collapsed[x]&&!collapsed[y])
            {
                if (x<y) cm->Remove(x, y, graph, core);
                else cm->Remove(y, x, graph, core);
            }
        for (auto y:edges[x].Less)
            if(!collapsed[x]&&!collapsed[y])
            {
                if (x<y) cm->Remove(x, y, graph, core);
                else cm->Remove(y, x, graph, core);
            }
        for (auto y:edges[x].Equal)
            if(!collapsed[x]&&!collapsed[y])
            {
                if (x<y) cm->Remove(x, y, graph, core);
                else cm->Remove(y, x, graph, core);
            }
        collapsed[Max_follower_pos]=1;
        //update neighbor
        for (int x=0;x<n;x++)
        {
            edges[x].More.clear();
            edges[x].Less.clear();
            edges[x].Equal.clear();
            for  (auto y:graph[x])
                if (core[y]>core[x]) edges[x].More.push_back(y);
                else if (core[y]<core[x]) edges[x].Less.push_back(y);
                else edges[x].Equal.push_back(y);
        }
        //get lowest coreness
        int Min_coreness=pre[Max_follower_pos];
        for (int i=0;i<n;i++)
        {
            if (!collapsed[i]&&core[i]<pre[i]) Min_coreness=min(Min_coreness,core[i]);
            else update[i]=1;
            bin[i].clear();
            tree[i].id=0;
            critical[i]=0;
        }
        bin[n].clear();
        //get critical vertex
        int Min_coreness_all=n+100,Max_coreness=0;
        for (int i=0;i<n;i++)
            if(!collapsed[i])
            {
                int More_degree=0;
                for (auto x:edges[i].More)
                    if (!collapsed[x]) More_degree++;
                for (auto x:edges[i].Equal)
                    if (!collapsed[x]) More_degree++;
                if (More_degree==core[i])
                {
                    critical[i]=1;
                    for (auto x:edges[i].Equal)
                        if (!collapsed[x]) update[x]=1;
                    for (auto x:edges[i].More)
                        if (!collapsed[x]) update[x]=1;
                }
                Min_coreness_all=min(Min_coreness_all,core[i]);
                Max_coreness=max(Max_coreness,core[i]);
                bin[core[i]].push_back(i);
            }
        //build core component tree
        vector<int> G;
        for (int i=1;i<=n;i++)
            for (auto x:bin[i])
                G.push_back(x);
        pre_LCPS(graph,core,n,Max_coreness);
        nodeCCT=0;
        BuildCCT(n,G,core,0,Min_coreness_all);
        for (int i=0;i<n;i++)
          pre[i]=core[i];
        Max_follower=0,Max_follower_pos=0;
        int sum=0,change_sum=0,ub_num=0,update_num=0,min_coreness_num=0,num=0,ub2_num=0,lb_num=0;
        vector<vector<int> > bound(n+1);
        for (int i=0;i<n;i++)
            if (!collapsed[i])
            {
                upperbound[i]=upper_bound(n,i,core);
                lowerbound[i]=lower_bound(i);
                bound[upperbound[i]].push_back(i);
            }
        for (int upb=n;upb>=0;upb--)
            for (auto i:bound[upb])
        //for (int i=0;i<n;i++)
            {
                int UB,LB;
                if (!collapsed[i])
                {
                    num++;
                    UB=upperbound[i];
                    LB=lowerbound[i];
                    if (UB<=Max_follower) ub_num++;
                    if (!update[i]) update_num++;//follow=0
                    if (core[i]<Min_coreness) min_coreness_num++;//<min change coreness
                    if (LB==UB) lb_num++;
                }
                if (update[i] //follow=0
                    &&!collapsed[i]
                    &&core[i]>=Min_coreness //< Min changing nodes' coreness must not change
                    &&core[i]>0//&&ub(i)>Max_follower
                    &&UB>Max_follower
                    ) //upper_bound
                {
                    if (LB==UB) followers[i]=UB; //lower_bound
                    else
                    {
                        int x=i;
                        for (auto y:edges[x].Less)
                        {
                            if (x<y) cm->Remove(x, y, graph, core);
                            else cm->Remove(y, x, graph, core);
                        }
                        for (auto y:edges[x].Equal)
                        {
                            if (x<y) cm->Remove(x, y, graph, core);
                            else cm->Remove(y, x, graph, core);
                        }
                        int pre_followers=followers[i];
                        followers[i]=0;
                        for (int j=0;j<n;j++)
                            if (i!=j&&core[j]<pre[j]&&!collapsed[j]) followers[i]++;
                        sum++;
                        if (pre_followers==followers[i]) change_sum++;
                        for (auto y:edges[x].Less)
                        {
                            if (x<y) cm->Insert(x, y, graph, core);
                            else cm->Insert(y, x, graph, core);
                        }
                        for (auto y:edges[x].Equal)
                        {
                            if (x<y) cm->Insert(x, y, graph, core);
                            else cm->Insert(y, x, graph, core);
                        }
                    }
                }
                if (!collapsed[i]&&UB>Max_follower)
                    if (followers[i]>Max_follower) Max_follower=followers[i],Max_follower_pos=i;
            }
    }
    return 0;
}
