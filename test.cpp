#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "main.h"

TEST_CASE("Test findShortestPath") {
    MultiGraph graph(6);
    graph.addEdge(0, 1, 2);
    graph.addEdge(0, 2, 4);
    graph.addEdge(1, 2, 1);
    graph.addEdge(1, 3, 7);
    graph.addEdge(2, 3, 1);
    graph.addEdge(2, 4, 5);
    graph.addEdge(3, 4, 3);
    graph.addEdge(3, 5, 2);
    graph.addEdge(4, 5, 4);

    SUBCASE("Shortest path exists") {
        vector<int> expected = {0, 1, 2, 3, 5};
        vector<int> path = graph.findShortestPath(0, 5);
        CHECK(path == expected);
    }
}

TEST_CASE("Test findMaxMatching") {
    MultiGraph graph(6);
    graph.addEdge(0, 1, 1);
    graph.addEdge(0, 2, 1);
    graph.addEdge(1, 2, 1);
    graph.addEdge(1, 3, 1);
    graph.addEdge(2, 3, 1);
    graph.addEdge(2, 4, 1);
    graph.addEdge(3, 4, 1);
    graph.addEdge(3, 5, 1);
    graph.addEdge(4, 5, 1);

    vector<Edge> expected = {Edge(0, 2, 1), Edge(1, 3, 1), Edge(2, 4, 1)};
    vector<Edge> matchingEdges;
    int maxMatching = graph.findMaxMatching(matchingEdges);
    CHECK(maxMatching == 3);
}

TEST_CASE("Test findMaxFlow") {
    MultiGraph graph(6);
    graph.addEdge(0, 1, 16);
    graph.addEdge(0, 2, 13);
    graph.addEdge(1, 2, 10);
    graph.addEdge(1, 3, 12);
    graph.addEdge(2, 1, 4);
    graph.addEdge(2, 4, 14);
    graph.addEdge(3, 2, 9);
    graph.addEdge(3, 5, 20);
    graph.addEdge(4, 3, 7);
    graph.addEdge(4, 5, 4);

    int source = 0;
    int sink = 5;
    vector<int> path;
    int expectedMaxFlow = 23;
    int actualMaxFlow = graph.findMaxFlow(source, sink, path);

    CHECK(expectedMaxFlow == actualMaxFlow);
    CHECK(path == vector<int>({0, 2, 4, 3, 5}));
}