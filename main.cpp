#include <algorithm>
#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <numeric>
#include <functional>
#include <iterator>

#include "channel.h"
#include "MCEngine.h"
#include "Metropolis-Hasting.h"
#include "utils.h"

struct Edge {
    int source;
    int destination;
};

struct EdgeProperties {
    int capacity;
    int flow;
};

class Graph {
public:
    explicit Graph(size_t vertices_num) {
        adjacency_lists_.resize(vertices_num);
    }
    void AddEdge(int source, int dest) {
        adjacency_lists_[source].push_back(edges_.size());
        edges_.push_back({source, dest});
    }
    IteratorRange<std::vector<int>::const_iterator> GetOutgoingEdges(int source) const {
        return {adjacency_lists_[source].begin(), adjacency_lists_[source].end()};
    }
    int GetTarget(int edge_id) const {
        return edges_[edge_id].destination;
    }
    int GetSource(int edge_id) const {
        return edges_[edge_id].source;
    }
    size_t VertexCount() const {
        return adjacency_lists_.size();
    }

private:
    std::vector<std::vector<int>> adjacency_lists_;
    std::vector<Edge> edges_;
};

template <class Predicate>
class FilteredGraph {
public:
    FilteredGraph(Graph graph, Predicate predicate)
        : graph_(graph), predicate_(predicate) {}
    IteratorRange<FilterIterator<std::vector<int>::const_iterator, Predicate>>
        GetOutgoingEdges(int source) const {
        auto edges = graph_.GetOutgoingEdges(source);
        auto begin = FilterIterator<std::vector<int>::const_iterator, Predicate>
            (edges.begin(), edges.end(), predicate_);
        auto end = FilterIterator<std::vector<int>::const_iterator, Predicate>
            (edges.end(), edges.end(), predicate_);
        return {begin, end};
    }
    int GetTarget(int edge_id) const {
        return graph_.GetTarget(edge_id);
    }
    int GetSource(int edge_id) const {
        return graph_.GetSource(edge_id);
    }
    size_t VertexCount() const {
        return graph_.VertexCount();
    }

private:
    Graph graph_;
    Predicate predicate_;
};

template <class Graph>
class BFSVisitor {
public:
    explicit BFSVisitor(std::vector<int>* dist)
        : distance_(*dist) {}
    void ExamineEdge(int edge_id, const Graph& graph) {
        distance_[graph.GetTarget(edge_id)] = distance_[graph.GetSource(edge_id)] + 1;
    }
    bool IsVisited(int vertex) {
        return distance_[vertex] != -1;
    }
    void DiscoverVertex(int vertex) {
        distance_[vertex] = 0;
    }
    void Update() {
        distance_.assign(distance_.size(), -1);
    }

private:
    std::vector<int>& distance_;
};

template <class Graph, class Visitor>
void BreadthFirstSearch(int origin_vertex, const Graph& graph,
                        Visitor visitor) {
    visitor.Update();
    std::queue<int> order({origin_vertex});
    visitor.DiscoverVertex(origin_vertex);
    while (!order.empty()) {
        auto source = order.front();
        order.pop();
        for (auto& edge : graph.GetOutgoingEdges(source)) {
            auto target = graph.GetTarget(edge);
            if (!visitor.IsVisited(target)) {
                visitor.DiscoverVertex(target);
                visitor.ExamineEdge(edge, graph);
                order.push(target);
            }
        }
    }
}

class FlowNetwork {
public:
    FlowNetwork(Graph& graph, int source, int sink, std::vector<EdgeProperties>& edges_properties)
        : graph_(graph), source_(source), sink_(sink), edges_properties_(edges_properties) {}
    FilteredGraph<std::function<bool(int)>> ResidualNetworkView() const {
        return FilteredGraph<std::function<bool(int)>>(graph_, [&](int edge_id) {
            return edges_properties_[edge_id].capacity > edges_properties_[edge_id].flow;
        });
    }
    void SendFlow(int edge_id, int pushed) {
        edges_properties_[edge_id].flow += pushed;
        edges_properties_[edge_id ^ 1].flow -= pushed;
    }
    int PushFlow(std::vector<std::vector<int>::const_iterator>& ptrs, const std::vector<int>& distance,
        int start, int flow) {
        if (start == sink_ || flow == 0) {
            return flow;
        }
        auto edges = graph_.GetOutgoingEdges(start);
        while (ptrs[start] != graph_.GetOutgoingEdges(start).end()) {
            int id = *ptrs[start], target = graph_.GetTarget(id);
            if (distance[target] != distance[start] + 1) {
                ++ptrs[start];
                continue;
            }
            int pushed = PushFlow(ptrs, distance, target,
                std::min(edges_properties_[id].capacity - edges_properties_[id].flow, flow));
            if (pushed > 0) {
                SendFlow(id, pushed);
                return pushed;
            }
            ++ptrs[start];
        }
        return 0;
    }
    void InitializePointers(std::vector<std::vector<int>::const_iterator>& pointers) {
        for (int i = 0; i < static_cast<int>(pointers.size()); ++i) {
            pointers[i] = graph_.GetOutgoingEdges(i).begin();
        }
    }
    int Dinic() {
        int answer = 0;
        std::vector<int> distance(graph_.VertexCount());
        std::vector<std::vector<int>::const_iterator> pointers(graph_.VertexCount());
        do {
            InitializePointers(pointers);
            BFSVisitor<FilteredGraph<std::function<bool(int)>>> visitor(&distance);
            BreadthFirstSearch(source_, ResidualNetworkView(), visitor);
            while (int pushed = PushFlow(pointers, distance, source_, inf_)) {
                answer += pushed;
            }
        } while (distance[sink_] != -1);
        return answer;
    }

    friend class FlowNetworkBuilder;
private:
    Graph graph_;
    int source_, sink_;
    std::vector<EdgeProperties> edges_properties_;
    const int inf_ = 100000000;
};

class FlowNetworkBuilder {
public:
    explicit FlowNetworkBuilder(size_t vertices_num) : graph_(vertices_num) {
        source_ = 0;
        sink_ = vertices_num - 1;
    }
    void AddEdge(int source, int destination, int capacity) {
        graph_.AddEdge(source, destination);
        edges_properties_.push_back({capacity, 0});
        graph_.AddEdge(destination, source);
        edges_properties_.push_back({0, 0});
    }
    FlowNetwork Build() {
        return {graph_, source_, sink_, edges_properties_};
    }
    int GetSource() const {
        return source_;
    }
    int GetSink() const {
        return sink_;
    }

private:
    Graph graph_;
    int source_, sink_;
    std::vector<EdgeProperties> edges_properties_;
};

double BidirectionalSinglePath(double lambda_lhs, double lambda_rhs,
    uint64_t l_lhs, uint64_t l_rhs, double interest_rate, double on_chain_cost,
    double transaction_fee) {
    
    std::random_device rd;
    std::mt19937 rng(rd());

    std::exponential_distribution<double> exp_lhs(lambda_lhs), exp_rhs(lambda_rhs);

    double next_transaction_time_lhs = exp_lhs(rng);
    double next_transaction_time_rhs = exp_rhs(rng);

    double total_cost = 0.0;

    Market market(interest_rate);
    Channel channel(l_lhs, l_rhs, &market);

    constexpr double discount_factor_threshold = 1e-10;
    constexpr double additional_cost_threshold = 1e-5;

    while (true) {
        if (market.Discount(1.0) < discount_factor_threshold) {
            return MAXFLOAT;
        }

        if (next_transaction_time_lhs < next_transaction_time_rhs && channel.ProcessPayment(1)) {
            market.Update(next_transaction_time_lhs);
            total_cost += market.Discount(transaction_fee);
            next_transaction_time_rhs -= next_transaction_time_lhs;
            next_transaction_time_lhs = exp_lhs(rng);
            continue;
        }

        if (next_transaction_time_lhs > next_transaction_time_rhs && channel.ProcessPayment(-1)) {
            market.Update(next_transaction_time_rhs);
            total_cost += market.Discount(transaction_fee);
            next_transaction_time_lhs -= next_transaction_time_rhs;
            next_transaction_time_rhs = exp_rhs(rng);
            continue;
        }


        double additional_cost = market.Discount(channel.GetOpportunityCost())
            + market.Discount(2 * on_chain_cost);

        if (additional_cost < additional_cost_threshold) {
            return total_cost;
        }

        total_cost += additional_cost;
        channel.Reset(l_lhs, l_rhs);
    }

}

double UnidirectionalSinglePath(double lambda, uint64_t m,
    double interest_rate, double on_chain_cost, double transaction_fee) {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::exponential_distribution<double> exp(lambda);

    Market market(interest_rate);
    Channel channel(m, 0u, &market);

    double total_cost = 0.0;
    constexpr double additional_cost_threshold = 1e-5;

    while (true) {

        if (channel.ProcessPayment(1)) {
            market.Update(exp(rng));
            total_cost += market.Discount(transaction_fee);
            continue;
        }

        double additional_cost = market.Discount(channel.GetOpportunityCost())
            + market.Discount(2 * on_chain_cost);

        if (additional_cost < additional_cost_threshold) {
            return total_cost;
        }

        total_cost += additional_cost;

        channel.Reset(m, 0);
    }
}

double OnChainSinglePath(double lambda, double on_chain_cost,
    double interest_rate) {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::exponential_distribution<double> exp(lambda);

    double total_cost = 0.0;
    constexpr double additional_cost_threshold = 1e-5;

    Market market(interest_rate);

    while (true) {
        market.Update(exp(rng));
        double additional_cost = market.Discount(on_chain_cost);

        if (additional_cost < additional_cost_threshold) {
            return total_cost;
        }

        total_cost += additional_cost;
    }
}

void TestUnidirectional(double lambda, uint64_t m,
    double interest_rate, double on_chain_cost, double transaction_fee) {
    MCEngine engine;

    auto func = [lambda, m, interest_rate, on_chain_cost, transaction_fee]() {
        return UnidirectionalSinglePath(lambda, m, interest_rate, on_chain_cost, transaction_fee);
    };

    auto data = engine.Simulate<decltype(func), double>(func);

    std::cout << std::accumulate(data.begin(), data.end(), 0.0) / data.size() << std::endl;
}

void Test(double lambda, double on_chain_cost,
    double interest_rate) {
    MCEngine engine;

    auto func = [lambda, on_chain_cost, interest_rate]() {
        return OnChainSinglePath(lambda, on_chain_cost, interest_rate);
    };

    auto data = engine.Simulate<decltype(func), double>(func);

    std::cout << std::accumulate(data.begin(), data.end(), 0.0) / data.size() << std::endl;
}

int64_t Generator(int64_t balance) {
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<> dist(-100, 100);

    return std::max<int64_t>(balance + dist(rng), 1);
}

int64_t FindBidirectionalSymmetricOptimalBalance(double lambda_lhs = 1.0,
    double lambda_rhs = 1.0, double interest_rate = 1e-3, double on_chain_cost = 0.25,
    double transaction_fee = 1e-5) {
    MetropolisHasting mh(0.1);

    auto func = [lambda_lhs, lambda_rhs, interest_rate, on_chain_cost, transaction_fee](int64_t balance) {
        MCEngine engine;

        auto cost = [
            balance, lambda_lhs, lambda_rhs,
            interest_rate, on_chain_cost, transaction_fee
        ]() {
            return BidirectionalSinglePath(lambda_lhs, lambda_rhs, balance, balance,
                interest_rate, on_chain_cost, transaction_fee);
        };

        auto data = engine.Simulate<decltype(cost), double>(cost, 10u);

        return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    };

    return mh.Solver<decltype(func), decltype(Generator), int64_t, true>(func, Generator);
}

int main() {

    std::cout << FindBidirectionalSymmetricOptimalBalance() << std::endl;

    // TestUnidirectional(1.0, 10, 1e-2, 0.5);

    // std::cout << BidirectionalSinglePath(1.0, 1.0, 10, 10, 1e-3, 0.25) << std::endl;;

    return 0;
}
