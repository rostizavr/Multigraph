#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>
#include <limits>

using namespace std;

// Структура для представления ребра в мультиграфе
struct Edge
{
    int source;
    int destination;
    int capacity;

    Edge(int source, int destination, int capacity)
    {
        this->source = source;
        this->destination = destination;
        this->capacity = capacity;
    }
};

// Класс для представления мультиграфа
class MultiGraph
{
private:
    int numVertices;                    // Количество вершин
    vector<vector<Edge>> adjacencyList; // Список смежности

    // Обход в глубину (DFS) для поиска компонент и эйлерова цикла
    void dfs(int vertex, vector<bool>& visited, stack<int>& verticesStack)
    {
        visited[vertex] = true;

        for (const Edge& edge : adjacencyList[vertex]) {
            int nextVertex = edge.destination;
            if (!visited[nextVertex]) {
                dfs(nextVertex, visited, verticesStack);
            }
        }

        verticesStack.push(vertex);
    }

    // Поиск в ширину (BFS) для поиска увеличивающего пути в резидуальной сети
    bool bfs(const vector<vector<int>>& residualNetwork, int source, int sink,
             vector<int>& parent)
    {
        vector<bool> visited(numVertices, false);
        queue<int> q;
        q.push(source);
        visited[source] = true;
        parent[source] = -1;

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (int v = 0; v < numVertices; v++) {
                if (!visited[v] && residualNetwork[u][v] > 0) {
                    q.push(v);
                    parent[v] = u;
                    visited[v] = true;
                }
            }
        }

        return visited[sink];
    }

    // BFS для поиска увеличивающей цепи в паросочетании
    bool bfsForMatching(vector<int>& matching)
    {
        vector<int> dist(numVertices, -1);
        queue<int> q;

        for (int u = 0; u < numVertices; u++) {
            if (matching[u] == -1) {
                dist[u] = 0;
                q.push(u);
            }
        }

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (const Edge& edge : adjacencyList[u]) {
                int v = edge.destination;
                if (matching[v] == -1) {
                    return true; // Найдена увеличивающая цепь
                }
                if (dist[matching[v]] == -1) {
                    dist[matching[v]] = dist[u] + 1;
                    q.push(matching[v]);
                }
            }
        }

        return false; // Увеличивающая цепь не найдена
    }

    // DFS для поиска увеличивающей цепи и установки новых пар в паросочетание
    bool dfsForMatching(int u, vector<int>& matching)
    {
        for (const Edge& edge : adjacencyList[u]) {
            int v = edge.destination;
            if (matching[v] == -1) {
                matching[u] = v;
                matching[v] = u;
                return true; // Найдена увеличивающая цепь
            }
        }

        for (const Edge& edge : adjacencyList[u]) {
            int v = edge.destination;
            if (matching[v] != -1) {
                int w = matching[v];
                if (dfsForMatching(w, matching)) {
                    matching[u] = v;
                    matching[v] = u;
                    return true; // Найдена увеличивающая цепь
                }
            }
        }

        return false; // Увеличивающая цепь не найдена
    }

    // Вспомогательная функция для рекурсивного поиска максимальной клики
    void findMaximumCliqueUtil(vector<int>& currentClique, vector<int>& maxClique,
                               vector<bool>& visited)
    {
        // Базовый случай: все вершины уже посещены
        if (currentClique.size() == numVertices) {
            maxClique = currentClique;
            return;
        }

        // Найти первую непосещенную вершину
        int startVertex = -1;
        for (int i = 0; i < numVertices; i++) {
            if (!visited[i]) {
                startVertex = i;
                break;
            }
        }

        // Перебираем соседей непосещенной вершины
        for (const Edge& edge : adjacencyList[startVertex]) {
            int neighbor = edge.destination;
            if (!visited[neighbor]) {
                // Проверяем, является ли текущая клика допустимой
                bool isClique = true;
                for (int vertex : currentClique) {
                    if (!isAdjacent(vertex, neighbor)) {
                        isClique = false;
                        break;
                    }
                }

                // Если текущая клика является допустимой, добавляем вершину в клику и
                // продолжаем рекурсивно
                if (isClique) {
                    currentClique.push_back(neighbor);
                    visited[neighbor] = true;
                    findMaximumCliqueUtil(currentClique, maxClique, visited);
                    visited[neighbor] = false;
                    currentClique.pop_back();
                }
            }
        }
    }

    // Проверка, являются ли две вершины смежными
    bool isAdjacent(int vertex1, int vertex2)
    {
        for (const Edge& edge : adjacencyList[vertex1]) {
            if (edge.destination == vertex2) {
                return true;
            }
        }
        return false;
    }

public:
    MultiGraph(int numVertices)
    {
        this->numVertices = numVertices;
        adjacencyList.resize(numVertices);
    }

    // Добавление ребра в мультиграф
    void addEdge(int source, int destination, int capacity)
    {
        Edge edge(source, destination, capacity);
        adjacencyList[source].push_back(edge);
    }

    // Алгоритм поиска максимальной клики
    vector<int> findMaximumClique()
    {
        vector<int> currentClique; // Текущая клика
        vector<int> maxClique;     // Максимальная клика

        vector<bool> visited(numVertices,
                             false); // Массив для отслеживания посещенных вершин

        findMaximumCliqueUtil(currentClique, maxClique, visited);

        return maxClique;
    }

    // Поиск максимального потока в сети с использованием алгоритма
    // Форда-Фалкерсона
    int findMaxFlow(int source, int sink, vector<int>& path)
    {
        // Создаем резидуальную сеть с начальным потоком 0
        vector<vector<int>> residualNetwork(numVertices,
                                            vector<int>(numVertices, 0));
        for (int u = 0; u < numVertices; u++) {
            for (const Edge& edge : adjacencyList[u]) {
                residualNetwork[edge.source][edge.destination] += edge.capacity;
            }
        }

        vector<int> parent(numVertices); // Массив для хранения родительских вершин

        int maxFlow = 0;

        // Пока существует увеличивающий путь в резидуальной сети
        while (bfs(residualNetwork, source, sink, parent)) {
            int pathFlow = numeric_limits<int>::max();

            // Находим минимальную пропускную способность в пути
            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                pathFlow = min(pathFlow, residualNetwork[u][v]);
            }

            // Обновляем значения ребер в резидуальной сети
            for (int v = sink; v != source; v = parent[v]) {
                int u = parent[v];
                residualNetwork[u][v] -= pathFlow;
                residualNetwork[v][u] += pathFlow;
            }

            maxFlow += pathFlow;
        }

        // Сохраняем путь потока
        path.clear();
        int currentVertex = sink;
        while (currentVertex != source) {
            path.push_back(currentVertex);
            currentVertex = parent[currentVertex];
        }
        path.push_back(source);
        reverse(path.begin(), path.end());

        return maxFlow;
    }

    // Поиск эйлерова цикла в мультиграфе
    vector<int> findEulerianCycle()
    {
        vector<int> eulerianCycle;

        // Проверяем, есть ли эйлеров цикл в мультиграфе
        if (!isEulerian()) {
            return eulerianCycle;
        }

        int startVertex = 0;
        while (!adjacencyList[startVertex].size()) {
            startVertex++;
        }

        stack<int> currentPath;
        vector<int> circuit;

        currentPath.push(startVertex);
        int currentVertex = startVertex;

        while (!currentPath.empty()) {
            if (adjacencyList[currentVertex].size()) {
                currentPath.push(currentVertex);

                int nextVertex = adjacencyList[currentVertex].back().destination;
                adjacencyList[currentVertex].pop_back();

                currentVertex = nextVertex;
            } else {
                circuit.push_back(currentVertex);

                currentVertex = currentPath.top();
                currentPath.pop();
            }
        }

        // Переворачиваем цикл, чтобы получить эйлеров цикл
        for (int i = circuit.size() - 1; i >= 0; i--) {
            eulerianCycle.push_back(circuit[i]);
        }

        return eulerianCycle;
    }

    // Проверка наличия эйлерова цикла в мультиграфе
    bool isEulerian()
    {
        for (int i = 0; i < numVertices; i++) {
            if (adjacencyList[i].size() % 2 != 0) {
                return false;
            }
        }

        return true;
    }

    // Поиск сильных компонент в мультиграфе с использованием алгоритма Косарайю
    vector<vector<int>> findStronglyConnectedComponents()
    {
        vector<vector<int>> components;
        vector<bool> visited(numVertices, false);
        stack<int> verticesStack;

        // Первый проход DFS
        for (int i = 0; i < numVertices; i++) {
            if (!visited[i]) {
                dfs(i, visited, verticesStack);
            }
        }

        // Транспонирование мультиграфа
        MultiGraph transposedGraph = getTransposedGraph();

        visited.assign(numVertices, false);

        // Второй проход DFS
        while (!verticesStack.empty()) {
            int vertex = verticesStack.top();
            verticesStack.pop();

            if (!visited[vertex]) {
                vector<int> component;
                transposedGraph.dfs(vertex, visited, verticesStack);
                components.push_back(component);
            }
        }

        return components;
    }

    // Получение транспонированного мультиграфа
    MultiGraph getTransposedGraph()
    {
        MultiGraph transposedGraph(numVertices);

        for (int u = 0; u < numVertices; u++) {
            for (const Edge& edge : adjacencyList[u]) {
                transposedGraph.addEdge(edge.destination, edge.source, edge.capacity);
            }
        }

        return transposedGraph;
    }

    // Поиск максимального паросочетания в мультиграфе с использованием алгоритма
    // Хопкрофта-Карпа
    int findMaxMatching(vector<Edge>& matchingEdges)
    {
        vector<int> matching(numVertices, -1); // Массив, хранящий пары вершин

        int maxMatching = 0;

        while (bfsForMatching(matching)) {
            for (int u = 0; u < numVertices; u++) {
                if (matching[u] == -1 && dfsForMatching(u, matching)) {
                    maxMatching++;
                }
            }
        }

        // Сохраняем найденные ребра максимального паросочетания
        matchingEdges.clear();
        for (int u = 0; u < numVertices; u++) {
            if (matching[u] != -1) {
                int v = matching[u];
                matchingEdges.push_back(Edge(u, v, 1));
            }
        }

        return maxMatching;
    }

    // Поиск кратчайшего пути между двумя вершинами с использованием алгоритма
    // Дейкстры
    vector<int> findShortestPath(int source, int destination)
    {
        vector<int> distance(numVertices, numeric_limits<int>::max());
        distance[source] = 0;

        priority_queue<pair<int, int>, vector<pair<int, int>>,
                greater<pair<int, int>>> pq;
        pq.push(make_pair(0, source));

        vector<int> parent(numVertices, -1);

        while (!pq.empty()) {
            int u = pq.top().second;
            pq.pop();

            for (const Edge& edge : adjacencyList[u]) {
                int v = edge.destination;
                int weight = edge.capacity;

                if (distance[u] != numeric_limits<int>::max() &&
                    distance[u] + weight < distance[v]) {
                    distance[v] = distance[u] + weight;
                    parent[v] = u;
                    pq.push(make_pair(distance[v], v));
                }
            }
        }

        // Восстановление пути из родительских вершин
        vector<int> path;
        int current = destination;
        while (current != -1) {
            path.push_back(current);
            current = parent[current];
        }

        reverse(path.begin(), path.end());

        return path;
    }
};