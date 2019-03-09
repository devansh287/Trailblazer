#include "trailblazer.h"
#include "queue.h"
#include "pqueue.h"
#include <climits>

using namespace std;
bool depthFirstSearchHelper (BasicGraph& graph, Vertex* start, Vertex* end, Vector <Vertex*>& path, Set <Vertex*> visited);
Vector <Vertex*> HashMap_to_vector (Vertex* key, HashMap <Vertex*, Vertex*>& previous);
Vector <Edge*> vector_to_edges (BasicGraph& graph, Vector <Vertex*> v);
Vector <Vertex*> min_without_edge (BasicGraph& graph, Vertex* start, Vertex* end, Edge* avoid, double& weight);
double difference_between_vectors (Vector <Vertex*> original, Vector <Vertex*> alternate);

//Depth first search
Vector<Vertex*> depthFirstSearch (BasicGraph& graph, Vertex* start, Vertex* end)
{
    Vector<Vertex*> path;

    //set of vertexes that have been visited
    Set <Vertex*> visited;

    //helper function that recursively executes the DFS algorithm
    depthFirstSearchHelper (graph, start, end, path, visited);
    return path;
}

//aforementioned helper function
bool depthFirstSearchHelper (BasicGraph& graph, Vertex* start, Vertex* end, Vector <Vertex*>& path, Set <Vertex*> visited)
{
    //return value
    bool b = false;

    //if the end has not been reached
    if (start != end)
    {
        //start has been visited
        visited.add(start);

        //choose by adding the vertex to the path
        path.add(start);
        start->setColor(GREEN);

        //running the algorithm on all neighbours of vertex
        for (Vertex* vertex : graph.getNeighbors(start))
        {
            if (!visited.contains(vertex)) b = depthFirstSearchHelper (graph, vertex, end, path, visited);
            if (b) break;
        }

        // in case the path has not been found from the neighbours

        if (!b)
        {
            start->setColor(GRAY);

            //UNCHOOSE
            path.removeValue(start);
        }
    }

    //when we find the end
    else
    {
        path.add (start);
        visited.add (start);
        start->setColor(GREEN);
        b = true;
    }
    return b;
}

//breadth first seach
Vector<Vertex*> breadthFirstSearch(BasicGraph& graph, Vertex* start, Vertex* end)
{
    Vector<Vertex*> path;

    //queue for the vertexes to be visited
    Queue <Vertex*> toVisit;

    //set of visited vertexes
    Set <Vertex*> visited;

    //a map that keeps a record of each vertex's previous  vertex
    HashMap <Vertex*, Vertex*> previous;

    toVisit.enqueue(start);
    path.add(start);
    visited.add(start);

    //indicator variable that detects when we reach the end
    bool flag = true;

    while (!toVisit.isEmpty() && flag)
    {
        //vertex that is dequeued
        Vertex* toExplore = toVisit.dequeue();
        toExplore->setColor(GREEN);

        //checking all neighbours
        for (Vertex* neighbor : graph.getNeighbors(toExplore))
        {
            if (!visited.contains(neighbor))
            {
                //setting neighbor's previous vertex to toExplore
                previous.add(neighbor, toExplore);
                toVisit.enqueue (neighbor);
                visited.add(neighbor);
                neighbor->setColor(YELLOW);

                //ending the loop when the end if found
                if(neighbor == end)
                {
                    flag = false;
                }
            }
        }
    }

    //helper function that converts a hash map key into a vector using its previous
    path = HashMap_to_vector (end, previous);

    return path;
}

//aforementioned helper function
Vector <Vertex*> HashMap_to_vector (Vertex* key, HashMap <Vertex*, Vertex*>& previous)
{
    Vector <Vertex*> path;
    path.add(key);
    while (previous.containsKey (key))
    {
        path.insert(0,previous [key]);
        key = previous [key];
    }

    return path;
}

Vector<Vertex*> dijkstrasAlgorithm(BasicGraph& graph, Vertex* start, Vertex* end)
{
    Vector<Vertex*> path;

    //priority queue
    PriorityQueue <Vertex*> pq;

    //set of visited vertexes
    Set <Vertex*> visited;

    //map that tracks the previous vertex
    HashMap <Vertex*, Vertex*> previous;

    //map that keeps track of each vertex's cost
    HashMap <Vertex*, double> costs;

    //initializing the above data structures with start vertex
    pq.enqueue(start, 0);
    visited.add(start);
    costs [start] = 0;

    bool flag = true;

    while (!pq.isEmpty() && flag)
    {
        Vertex* explore = pq.dequeue();
        explore->setColor(GREEN);
        visited.add(explore);

        //if the end has been found
        if (explore == end) break;

        for (Vertex* neighbor : graph.getNeighbors(explore))
        {
            if(!visited.contains(neighbor))
            {
                //finding the cost of visiting that neighbour
                double cost = costs [explore] + graph.getEdge(explore, neighbor)->cost;
                neighbor->setColor(YELLOW);

                //updating the data structures if the neighbour is unexplored or if it is cheaper to reach

                if (!costs.containsKey(neighbor))
                {
                    pq.enqueue(neighbor, cost);
                    previous [neighbor] =  explore;
                    costs [neighbor] = cost;
                }

                else if (costs [neighbor] > cost)
                {
                    pq.changePriority (neighbor, cost);
                    previous [neighbor] =  explore;
                    costs [neighbor] = cost;
                }
            }

        }
    }

    //using the helper function
    path = HashMap_to_vector(end, previous);
    return path;
}

//everything is the same here, except the heuristic function
Vector<Vertex*> aStar(BasicGraph& graph, Vertex* start, Vertex* end)
{
    double weight = 0;
    //using a helper function that I used for alternate path
    return min_without_edge (graph, start, end, nullptr, weight);
}

//alternate path algorithm
Vector<Vertex*> alternatePath(BasicGraph& graph, Vertex* start, Vertex* end, double difference)
{
    Vector<Vertex*> final_path;

    //finding the shortest path using the a star algorithm
    Vector<Vertex*> min_path = aStar(graph,start,end);

    //helper function that converts the vectors in the min path to edges of those edges
    Vector <Edge*> min_path_edges = vector_to_edges (graph, min_path);

    //map of candidate paths along with their weights
    HashMap <Vector <Vertex*>, double> candidates;

    for (int i = 0; i < min_path_edges.size(); i++)
    {
        double weight = 0.0;

        //finding the alternate path with its weight using a helper function that uses
        //a* algorithm

        Vector <Vertex*> probable = min_without_edge(graph, start, end, min_path_edges [i], weight);

        //if the difference between the paths is greater than the required differnce, add it to the candidates map
        if (difference_between_vectors(min_path, probable) > difference)
        {
            candidates [probable] = weight;
        }
    }

    //finding the path with the minimum weight by looping through the map
    double min_weight = POSITIVE_INFINITY;

    for (Vector <Vertex*> probable: candidates)
    {
        if (min_weight > candidates [probable])
        {
            final_path = probable;
            min_weight = candidates [probable];
        }
    }

    return final_path;
}

//helper function that converts vector of vertices to a vector of edges
Vector <Edge*> vector_to_edges (BasicGraph& graph, Vector <Vertex*> v)
{
    Vector <Edge*> s;
    int i;
    int n = v.size();
    for (i = 0; i < n-2; i++)
    {
        Edge* e = graph.getEdge(v[i] , v[i+1]);
        s.add(e);
    }

    return s;
}

//helper function that uses a* algorithm to find the shortest alternate path
Vector <Vertex*> min_without_edge (BasicGraph& graph, Vertex* start, Vertex* end, Edge* avoid, double& weight)
{
    Vector<Vertex*> path;
    PriorityQueue <Vertex*> pq;

    //set of visited vertexes
    Set <Vertex*> visited;

    //map that maps a vertex's previous vertex
    HashMap <Vertex*, Vertex*> previous;

    //map that maps a vertex's cost
    HashMap <Vertex*, double> costs;

    //basically, everything is the same as dijkstra's algorithm except that
    //we enqueue vertexes with their heuristic values added

    pq.enqueue(start, heuristicFunction(start, end));
    visited.add(start);
    costs [start] = 0;

    bool flag = true;

    while (!pq.isEmpty() && flag)
    {
        Vertex* explore = pq.dequeue();
        visited.add(explore);
        explore->setColor(GREEN);
        if (explore == end) break;

        for (Vertex* neighbor : graph.getNeighbors(explore))
        {
            //this is important: here we check if the edge we're exploring is not the one
            //we wanna skip
            if(!visited.contains(neighbor) && (graph.getEdge(explore, neighbor) != avoid))
            {
                double cost = costs [explore] + graph.getEdge(explore, neighbor)->cost;
                neighbor->setColor(YELLOW);

                if (!costs.containsKey(neighbor))
                {
                    pq.enqueue(neighbor, cost + heuristicFunction(neighbor,end));
                    previous [neighbor] =  explore;
                    costs [neighbor] = cost;
                }

                else if (costs [neighbor] > cost)
                {
                    pq.changePriority (neighbor, cost + heuristicFunction(neighbor,end));
                    previous [neighbor] =  explore;
                    costs [neighbor] = cost;
                }
            }

        }
    }

    //setting the weight to the cost of the destination vertex
    weight = costs [end];

    path = HashMap_to_vector(end, previous);
    return path;
}

//helper function that finds the difference between two vectors
double difference_between_vectors (Vector <Vertex*> original, Vector <Vertex*> alternate)
{
    double not_in_original = 0.0;
    for (int i = 0; i < alternate.size(); i++)
    {
        if(!original.contains(alternate[i])) not_in_original++;
    }

    double original_total = original.size() * 1.0;
    double difference = not_in_original/original_total;
    return difference;

}

//Kruskal's algorithm
Set<Edge*> kruskal(BasicGraph& graph)
{
    Set<Edge*> mst;
    PriorityQueue <Edge*> pq;

    //a vector, which contains all clusters
    //clusters are hash sets of vertexes connected to each other
    Vector <HashSet<Vertex*>> clusters;

    //inserting an edge into the priority queue as per the cost
    for (Edge* e: graph.getEdgeSet())
    {
        pq.enqueue(e, e->cost);
    }

    //initializing each vertex to be a separate cluster in its own right
    Set <Vertex*> vertices = graph.getVertexSet();
    for (Vertex* current: vertices)
    {
        HashSet <Vertex*> cluster;
        cluster.add(current);
        clusters.add(cluster);
    }

    //creating the maze now
    while (!pq.isEmpty())
    {
        Edge* e = pq.dequeue();
        Vertex* v1 = e->start;
        Vertex* v2 = e->finish;

        //indicators that store the location of each cluster v1 and v2 are in
        int cluster1 = 0;
        int cluster2 = 0;

        //indicator that checks if they're in the same cluster or not
        bool indicator = false;
        for (int i = 0; i < clusters.size(); i++)
        {
            //if the vertexes are in the same cluster, then break and do nothing
           if (clusters[i].contains(v1) && clusters[i].contains(v2))
           {
               indicator = true;
               break;
           }
           else
           {
               //find out which cluster these vertexes are in
               if (clusters[i].contains(v1)) cluster1 = i;
               else if (clusters[i].contains(v2)) cluster2 = i;
           }

         }

        //if the vertex's are not in the same cluster
        if (!indicator)
        {
            //merge the clusters, add it to the vector, and delete the original two clusters
            HashSet <Vertex*> merged = clusters[cluster1] + clusters[cluster2];
            clusters.add(merged);
            if (cluster1 > cluster2)
            {
                clusters.remove(cluster1);
                clusters.remove(cluster2);
            }
            else
            {
                clusters.remove(cluster2);
                clusters.remove(cluster1);
            }
            //add the edge to the graph
            mst.add(e);
        }

    }



    return mst;
}
