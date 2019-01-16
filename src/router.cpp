/**************************************************************************
 * File       [ main.cpp ]
 * Author     [ wlkb83 ]
 * Synopsis   [ demonstration for the usage of parser.h ]
 * Usage      [ ./parser [inputfileName] ]
 * Date       [ 2014/12/28 created ]
**************************************************************************/

#include "parser.h"
#include <iostream>
#include <vector>
#include <string.h>
#include <deque>
#include <fstream>
#include <algorithm>
#include <climits>
using namespace std;

//record the edge weight
typedef struct Weight_Matrix{
    pair<int, int> u;
    int weight;
    int d;
}weight_matrix;

//record the routing path
typedef struct Route_Matrix{
    pair<int, int> v;
}route_matrix;

void dijkstra(vector< vector < pair < int, int> > >,int, int, int,int);

void check_adjacency(pair<int,int>&,weight_matrix**,int**,int**,int,int,int);

void relaxation_y(pair<int,int>&,pair<int,int>&,weight_matrix**,int**,int,int,int);
void relaxation_x(pair<int,int>&,pair<int,int>&,weight_matrix**,int**,int,int,int);

bool compare_weight(pair<int,int>&,pair<int, int>&,weight_matrix**,int**,int);

int calculate_weight(int**,pair<int,int>&,int);

void addline(weight_matrix**,route_matrix**,vector< vector<pair<int, int> >  >,int,int**,int**);

void write_to_file(vector< vector< pair< int, int> > >&, weight_matrix** ,route_matrix**, int);

//queue use for BFS
deque<pair< int,int> > queue;
ofstream outfile;
int main(int argc, char **argv)
{
    if( argc < 3 ){ cout << "Usage: ./parser [input_file_name] [output_file_name]" << endl; return 1; }
    outfile.open(argv[2],ios::trunc);
    AlgParser parser;
    
    // read the file in the first argument
    if( ! parser.read( argv[1] ) ) { return 1; }
    //build grid
    
    int h = parser.gNumHTiles(), v = parser.gNumVTiles();
    
    //capacity
    int capacity = parser.gCapacity();
    //num net
    int num_net = parser.gNumNets();
    //use pair to store the points of the nets
    pair<int, int> tmp(0,0);
    vector < pair < int, int > > tmp2(2,tmp);
    vector< vector< pair < int, int > > >points(num_net,tmp2);
    for (int idNet = 0; idNet < parser.gNumNets(); ++idNet){
        pair<int, int> posS = parser.gNetStart( idNet );
        pair<int, int> posE = parser.gNetEnd( idNet );
        points[idNet][0] = posS;
        points[idNet][1] = posE;
//        cout << idNet << " " << posS.first << " " << posS.second << " "
//                             << posE.first << " " << posE.second << endl;
    }
    
    dijkstra(points, h, v, num_net, capacity);
    outfile.close();
    return 0;
}
void dijkstra(vector< vector< pair< int, int> > > start_and_end, int horizontal, int vertical,int net_num, int capacity){
    int **trace_x = new int*[horizontal];
    int **trace_y = new int*[horizontal];
    for(int i = 0 ; i < horizontal;i++){
        trace_x[i] = new int[vertical]();
        trace_y[i] = new int[vertical]();
    }
    
    //for every net num
    for(int i = 0 ; i < net_num;i++){
        //open up weight matrix
        weight_matrix** weight = new weight_matrix*[horizontal];
        route_matrix ** route = new route_matrix*[horizontal];
        for(int j = 0; j < horizontal;j++){
            weight[j] = new weight_matrix[vertical];
            route[j] = new route_matrix[vertical];
        }
        //initial weight matrix
        for(int j = 0 ; j < horizontal;j++){
            for(int k = 0; k < vertical;k++){
                weight[j][k].weight = INT_MAX;
                weight[j][k].d = 0;
            }
        }
        //set starting point for weight matrix
        weight[start_and_end[i][0].first][start_and_end[i][0].second].weight = 0;
        weight[start_and_end[i][0].first][start_and_end[i][0].second].u.first = start_and_end[i][0].first;
        weight[start_and_end[i][0].first][start_and_end[i][0].second].u.second = start_and_end[i][0].second;
        //push starting point into queue
        queue.push_back(start_and_end[i][0]);

        while(!queue.empty()){
            pair<int,int> pre_front = queue.front();
            queue.pop_front();
            check_adjacency(pre_front, weight, trace_x, trace_y, capacity, horizontal, vertical);
        }
        
        addline(weight, route, start_and_end, i, trace_x, trace_y);
        write_to_file(start_and_end, weight, route, i);
    }
}
void check_adjacency(pair<int ,int>&check_point,weight_matrix** weight, int **trace_x,int **trace_y,int capacity,int horizontal,int vertical){
    pair<int,int> up(check_point.first,check_point.second + 1);
    pair<int,int> down(check_point.first,check_point.second - 1);
    pair<int,int> left(check_point.first-1,check_point.second);
    pair<int,int> right(check_point.first + 1,check_point.second);
    
    //check up point out of range or not
    relaxation_y(check_point, up, weight, trace_y, capacity, horizontal, vertical);
    relaxation_y(check_point, down, weight, trace_y, capacity, horizontal, vertical);
    relaxation_x(check_point, right, weight, trace_x, capacity, horizontal, vertical);
    relaxation_x(check_point, left, weight, trace_x, capacity, horizontal, vertical);
}

void relaxation_y(pair<int ,int> &check_point,pair<int,int> &point ,weight_matrix** weight,int **trace_y,int capacity,int horizontal,int vertical){
    
    if(point.first >= 0 && point.first < horizontal && point.second >= 0 && point.second < vertical){
        if(compare_weight(check_point, point, weight, trace_y, capacity)){
            //v.d = u.d + w(u,v)
            weight[point.first][point.second].weight = weight[check_point.first][check_point.second].weight + calculate_weight(trace_y, point, capacity);
            weight[point.first][point.second].d = weight[check_point.first][check_point.second].d + 1;
            weight[point.first][point.second].u.first = check_point.first;
            weight[point.first][point.second].u.second = check_point.second;
            
            if(find(queue.begin(), queue.end(), point) == queue.end()){
                queue.push_back(point);
            }
        }
        
    }
}

void relaxation_x(pair<int ,int> &check_point,pair<int,int> &point ,weight_matrix** weight, int **trace_x,int capacity,int horizontal,int vertical){
    
    if(point.first >= 0 && point.first < horizontal && point.second >= 0 && point.second < vertical){
        if(compare_weight(check_point, point, weight, trace_x, capacity)){
            //v.d = u.d + w(u,v)
            weight[point.first][point.second].weight = weight[check_point.first][check_point.second].weight + calculate_weight(trace_x, point, capacity);
            //distance + 1
            weight[point.first][point.second].d = weight[check_point.first][check_point.second].d + 1;
            //adjust the previous node
            weight[point.first][point.second].u.first = check_point.first;
            weight[point.first][point.second].u.second = check_point.second;
            
            if(find(queue.begin(), queue.end(), point) == queue.end()){
                queue.push_back(point);
            }
        }
    }
}

bool compare_weight(pair<int,int>& curr_point,pair<int, int>&next_point,weight_matrix **weight,int**trace_route,int capacity){
    bool result = false;
    // if v.d > u.d + w(u,v)
    if(weight[next_point.first][next_point.second].weight > weight[curr_point.first][curr_point.second].weight + calculate_weight(trace_route, next_point, capacity)){
        result = true;
    }
    return result;
}
int calculate_weight(int **trace_route,pair<int,int> &point,int capacity){
    int weight = 1;
    int demand = trace_route[point.first][point.second];
    
    if(demand!=0){
        weight = 2^demand;
        if(demand >= capacity)
            weight = 999;
    }
    return weight;
}

void addline(weight_matrix** w,route_matrix** r,vector< vector<pair<int, int> > > start_and_end,int i,int** trace_x,int** trace_y){
    pair<int,int> previous;
    pair<int,int> diff;
    int end_x = start_and_end[i][1].first;
    int end_y = start_and_end[i][1].second;
    r[end_x][end_y].v.first = end_x;
    r[end_x][end_y].v.second = end_y;
    //traverse from the end point to the start point
    while(end_x != start_and_end[i][0].first || end_y!= start_and_end[i][0].second){
        previous.first = w[end_x][end_y].u.first;
        previous.second = w[end_x][end_y].u.second;
        diff.first = previous.first - end_x;
        diff.second = previous.second - end_y;
        //last line is vertical
        if(diff.first == 0){
            trace_y[end_x][end_y]++;
        }
        //last line is horizontal
        else if(diff.second == 0){
            trace_x[end_x][end_y]++;
        }
        //set the routing table
        r[previous.first][previous.second].v.first = end_x;
        r[previous.first][previous.second].v.second = end_y;
        end_x = previous.first;
        end_y = previous.second;
    }
}

void write_to_file(vector< vector< pair< int, int> > >&start_and_end,weight_matrix** w, route_matrix** r, int i){
    outfile << i << " " << w[start_and_end[i][1].first][start_and_end[i][1].second].d << endl;
    
    // c:checking point  n:next point
    pair<int, int> check_point = start_and_end[i][0];
    pair<int, int> next_point;
    
    while( check_point.first != start_and_end[i][1].first|| check_point.second != start_and_end[i][1].second)
    {
        outfile << check_point.first << " " << check_point.second << " ";
        next_point.first= r[check_point.first][check_point.second].v.first;
        next_point.second = r[check_point.first][check_point.second].v.second;
        outfile << next_point.first << " " << next_point.second << endl;
        check_point.first = next_point.first;
        check_point.second= next_point.second;
    }
}
    

