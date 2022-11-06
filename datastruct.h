#include<vector>
#include<array>
#include<memory.h>
#include<math.h>
#include<algorithm>
using namespace std;
const int maxn = 1e2+10;
struct vert{
    double x, y, z;
};

template<typename index>
struct KD_tree{
    vector<vert> v, P, nodes;
    vector<index> indexes;
    vector<vector<int>> edges;
    vector<vert> dn;
    vector<int> mark1;
    vector<vector<int>> mark2;
    int lc[maxn],rc[maxn],tab[maxn],tab1[maxn];
    double L[maxn],R[maxn],D[maxn],U[maxn],F[maxn],
            B[maxn],lc1[maxn],rc1[maxn],L1[maxn],
            R1[maxn],D1[maxn],U1[maxn],F1[maxn],B1[maxn];
    double eta=1,lambda=0.1,minn=100000;
    int w=-1; 

    KD_tree(vector<vert> &_P, vector<index> &_indexes, vector<vert> &_nodes, vector<vector<int>> &_edges) : P(_P), indexes(_indexes), nodes(_nodes), edges(_edges){
        // memset(lc, 0, sizeof(lc));
        // memset(rc, 0, sizeof(rc));
        // memset(R, 0, sizeof(R));
        // memset(L, 0, sizeof(L));
        // memset(B, 0, sizeof(B));
        // memset(F, 0, sizeof(F));
        // memset(U, 0, sizeof(U));
        // memset(D, 0, sizeof(D));
        // memset(tab, 0, sizeof(tab));
        // memset(tab1, 0, sizeof(tab1));
        // memset(lc1, 0, sizeof(lc1));
        // memset(rc1, 0, sizeof(rc1));
        // memset(R1, 0, sizeof(R1));
        // memset(L1, 0, sizeof(L1));
        // memset(B1, 0, sizeof(B1));
        // memset(F1, 0, sizeof(F1));
        // memset(U1, 0, sizeof(U1));
        // memset(D1, 0, sizeof(D1));
        // vector<vert> tmp1;
        // dn = tmp1;
        // vector<vert> tmp4;
        // v = tmp4;
        // vector<int> tmp2;
        // mark1 = tmp2;
        // vector<vector<int>> tmp3;
        // mark2 = tmp3;
        for(size_t i=0;i<P.size();i++){
            v.push_back(P[i]);
            tab[i]=i;
        }
    }
    void set_tmp_array(){
        for(size_t i=0;i<P.size();i++){
            lc1[i]=lc[i];
            rc1[i]=rc[i];
            U1[i]=U[i];
            D1[i]=D[i];
            R1[i]=R[i];
            L1[i]=L[i];
            B1[i]=B[i];
            F1[i]=F[i];
            tab1[i]=tab[i];
            v.push_back(P[i]);
        }
    }
    
    void clear(){
        memset(lc, 0, sizeof(lc));
        memset(rc, 0, sizeof(rc));
        memset(R, 0, sizeof(R));
        memset(L, 0, sizeof(L));
        memset(B, 0, sizeof(B));
        memset(F, 0, sizeof(F));
        memset(U, 0, sizeof(U));
        memset(D, 0, sizeof(D));
    }
    double norm(){
        double ans=0;
        for(size_t i=1;i<dn.size();i++) ans=ans+(dn[i].x)*(dn[i].x)+(dn[i].y)*(dn[i].y)+(dn[i].z)*(dn[i].z);
        ans=sqrt(ans);
        return ans;
    }
    void maintain(int x) {
        L[x]=R[x]=v[tab[x]].x;
        D[x]=U[x]=v[tab[x]].y;
        F[x]=B[x]=v[tab[x]].z;
        if(lc[x]) L[x]=min(L[x],L[lc[x]]),R[x]=max(R[x],R[lc[x]]),D[x]=min(D[x],D[lc[x]]),U[x]=max(U[x], U[lc[x]]),F[x]=min(F[x],F[lc[x]]),B[x]=max(B[x],B[lc[x]]);
        if(rc[x]) L[x]=min(L[x],L[rc[x]]),R[x]=max(R[x],R[rc[x]]),D[x]=min(D[x],D[rc[x]]),U[x]=max(U[x], U[rc[x]]),F[x]=min(F[x],F[rc[x]]),B[x]=max(B[x],B[rc[x]]);
    }
    int build(int l,int r){
        if(l>r) return 0;
        if(l==r){
            maintain(l);
            return l;
        }
        int mid=(l+r)>>1;
        double av1=0,av2=0,av3=0,va1=0,va2=0,va3=0;  // average variance
        for(int i=l;i<=r;i++) av1+=v[tab[i]].x,av2+=v[tab[i]].y,av3+=v[tab[i]].z;
        av1/=(r-l+1);
        av2/=(r-l+1);
        av3/=(r-l+1);
        for (int i=l;i<=r;i++){
            va1+=(av1-v[tab[i]].x)*(av1-v[tab[i]].x);
            va2+=(av2-v[tab[i]].y)*(av2-v[tab[i]].y);
            va3+=(av3-v[tab[i]].z)*(av3-v[tab[i]].z);
        }
        // if((va1>va2)&&(va1>va3)) nth_element(tab+l,tab+mid,tab+r+1,cmp1);
        if((va1>va2)&&(va1>va3)) nth_element(tab+l,tab+mid,tab+r+1,[&](int i, int j){return v[i].x<v[j].x;});
        else if((va2>va1)&&(va2>va3)) nth_element(tab+l,tab+mid,tab+r+1,[&](int i, int j){return v[i].y<v[j].y;});
        else if((va3>va2)&&(va3>va1)) nth_element(tab+l,tab+mid,tab+r+1,[&](int i, int j){return v[i].z<v[j].z;});
        lc[mid]=build(l,mid-1);
        rc[mid]=build(mid+1,r);
        maintain(mid);
        return mid;
    }
    double dist(vert V, int b){
        double ans=0;
        if(L[b]>V.x) ans+=(L[b]-V.x)*(L[b]-V.x);
        if(R[b]<V.x) ans+=(R[b]-V.x)*(R[b]-V.x);
        if(D[b]>V.y) ans+=(D[b]-V.y)*(D[b]-V.y);
        if(U[b]<V.y) ans+=(U[b]-V.y)*(U[b]-V.y);
        if(F[b]>V.z) ans+=(F[b]-V.z)*(F[b]-V.z);
        if(B[b]<V.z) ans+=(B[b]-V.z)*(B[b]-V.z);
        return ans;
    }
    void query(int l, int r,vert V){
        if(l>r) return;
        int mid=(l+r)>>1;
        double t=(v[tab[mid]].x-V.x)*(v[tab[mid]].x-V.x)+(v[tab[mid]].y-V.y)*(v[tab[mid]].y-V.y)+(v[tab[mid]].z-V.z)*(v[tab[mid]].z-V.z);
        if(t<minn){
            w=tab[mid];
            minn=t;
        }
        if(l==r) return;
        double distl=1000000, distr=1000000;
        if(lc[mid]) distl=dist(V,lc[mid]);
        if(rc[mid]) distr=dist(V,rc[mid]);
        //cout<<l<<" "<<" "<<mid<<" "<<r<<" "<<minn<<" "<<distl<<" "<<distr<<endl;
        if ((distl<minn)&&(distr<minn)){
            if(distl<distr) {
            query(l,mid-1,V);
            if(distr<minn) query(mid+1,r,V);
            } 
            else{
            query(mid+1,r,V);
            if(distl<minn) query(l,mid-1,V);
            }
        } 
        else{
            if(distl<minn) query(l,mid-1,V);
            if(distr<minn) query(mid+1,r,V);
        }
    }

    void iterative_update(){
        // set_tmp_array(); 
        for(int k=1;k<10000;k++){
            vector<int> tmp;
            mark1.push_back(-1);
            mark2.push_back(tmp);
            v.clear();
            set_tmp_array();
            for(size_t i=1;i<nodes.size();i++){
                minn=10000000;
                w=-1;
                query(1,(int)P.size()-1, nodes[i]); 
                if(w>0) mark1.push_back(w);
                mark2.push_back(tmp);
            }
            v.clear();
            for(size_t i=0;i<nodes.size();i++){
                v.push_back(nodes[i]);
                tab[i]=i;
            } 
            clear();
            build(1,(int)nodes.size()-1);
            for(size_t i=0;i<mark2.size();i++) mark2[i].push_back(-1);
            for(size_t i=1;i<P.size();i++){
                minn=10000000;
                w=-1;
                query(1,(int)nodes.size()-1,P[i]); 
                if(w>0) mark2[w].push_back(i);
            }
            vert tmp4{0, 0, 0};
            for(size_t i=0;i<nodes.size();i++){
                dn.push_back(tmp4);
            }  
            for(size_t i=1;i<nodes.size();i++){
                dn[i].x=dn[i].x+2*(nodes[i].x-P[mark1[i]].x)/(double)(nodes.size()-1);
                dn[i].y=dn[i].y+2*(nodes[i].y-P[mark1[i]].y)/(double)(nodes.size()-1);
                dn[i].z=dn[i].z+2*(nodes[i].z-P[mark1[i]].z)/(double)(nodes.size()-1);
                for(size_t j=1;j<mark2[i].size();j++){
                    dn[i].x=dn[i].x+2*(nodes[i].x-P[mark2[i][j]].x)/(double)(P.size()-1);
                    dn[i].y=dn[i].y+2*(nodes[i].y-P[mark2[i][j]].y)/(double)(P.size()-1);
                    dn[i].z=dn[i].z+2*(nodes[i].z-P[mark2[i][j]].z)/(double)(P.size()-1);
                }
                tmp4.x = 0;
                tmp4.y = 0;
                tmp4.z = 0;
                for(size_t j=1;j<edges[i].size();j++){
                    tmp4.x = tmp4.x+2*lambda*nodes[edges[i][j]].x/(double)(edges[i].size()-1);
                    tmp4.y = tmp4.y+2*lambda*nodes[edges[i][j]].y/(double)(edges[i].size()-1);
                    tmp4.z = tmp4.z+2*lambda*nodes[edges[i][j]].z/(double)(edges[i].size()-1);
                }
                dn[i].x=dn[i].x+2*lambda*nodes[i].x-tmp4.x;
                dn[i].y=dn[i].y+2*lambda*nodes[i].y-tmp4.y;
                dn[i].z=dn[i].z+2*lambda*nodes[i].z-tmp4.z;
                for(size_t j=1;j<edges[i].size();j++){
                    dn[edges[i][j]].x=dn[edges[i][j]].x-2*lambda*nodes[i].x/(double)(edges[i].size()-1)+tmp4.x/(double)(edges[i].size()-1);
                    dn[edges[i][j]].y=dn[edges[i][j]].y-2*lambda*nodes[i].y/(double)(edges[i].size()-1)+tmp4.y/(double)(edges[i].size()-1);
                    dn[edges[i][j]].z=dn[edges[i][j]].z-2*lambda*nodes[i].z/(double)(edges[i].size()-1)+tmp4.z/(double)(edges[i].size()-1);
                }  
            }
            for(size_t i=1;i<nodes.size();i++){
                nodes[i].x=nodes[i].x-eta*dn[i].x;
                nodes[i].y=nodes[i].y-eta*dn[i].y;
                nodes[i].z=nodes[i].z-eta*dn[i].z;
            }
            double N=norm();
            cout<<k<<' '<<N<<endl;
            if(N<0.01) break;
            mark1.clear();
            mark2.clear();
            dn.clear(); 
        }
    }
};
