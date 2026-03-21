#include <cstddef>
#include <iostream>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <functional>
#include <cmath>

template<typename coeffType>
class PolyTrie
{
private:
    struct TrieNode 
    {
        std::map<int, std::unique_ptr<TrieNode>> childs;
        std::optional<coeffType> coeff;
    };

private:
    std::unique_ptr<TrieNode> m_root;
    std::vector<std::string> m_var_names;

public:
    PolyTrie(std::vector<std::string> const &init_var_names)
    {
        if (init_var_names.empty())
        {
            throw std::invalid_argument("Poly can't be initialized without variables");
        }
        m_root = std::make_unique<TrieNode>();
        m_var_names = init_var_names;
    }

public:
    void add_term(std::vector<int> degrees, coeffType coeff)
    {
        if (degrees.size() != m_var_names.size())
        {
            throw std::invalid_argument("Number of degrees doesn't match number of variables");
        }
        if (coeff == 0) {return;}
        if (m_root == nullptr)
        {
            m_root = std::make_unique<TrieNode>();
        }
        TrieNode *current = m_root.get();

        for (auto deg : degrees)
        {
            auto iter = current->childs.find(deg);
            if (iter == current->childs.end())
            {
                current->childs[deg] = std::make_unique<TrieNode>();
            }
            current = current->childs[deg].get();
        }
        if (current->coeff.has_value())
        {
            current->coeff = current->coeff.value() + coeff;
            if (current->coeff.value() == coeffType(0)) 
            {
                current->coeff.reset();
            }
        }
        else 
        {
            current->coeff = coeff;
        }
    }

    std::map<std::vector<int>, coeffType> get_supp() const
    {
        std::map<std::vector<int>, coeffType> supp;
        std::vector<int> current_degrees(m_var_names.size(), 0);
        
        std::function<void(const TrieNode*, size_t)> dfs_collect = 
            [&](const TrieNode* node, size_t depth) 
            {
                if (depth == m_var_names.size()) 
                {
                    if (node->coeff.has_value() && node->coeff.value() != coeffType(0)) 
                    {
                        supp[current_degrees] = node->coeff.value();
                    }
                    return;
                }
                
                for (auto const &[degree, child] : node->childs)
                {
                    current_degrees[depth] = degree;
                    dfs_collect(child.get(), depth + 1);
                }
            };

        dfs_collect(m_root.get(), 0);
        return supp;
    }

public:
    void print() const
    {
        auto supp = get_supp();
        bool is_first = true;
        
        for (auto const &[degrees, coeff] : supp)
        {
            if (coeff == coeffType(0)) continue;
            
            if (is_first)
            {
                if (coeff < coeffType(0))
                    std::cout << '-';
                is_first = false;
            }
            else
            {
                std::cout << (coeff < coeffType(0) ? " - " : " + ");
            }
            
            coeffType abs_coeff = coeff < coeffType(0) ? -coeff : coeff;
            bool has_variables = false;
            
            for (size_t i = 0; i < m_var_names.size(); ++i)
            {
                if (degrees[i] != 0)
                {
                    has_variables = true;
                    break;
                }
            }
            
            if (abs_coeff != coeffType(1) || !has_variables)
            {
                std::cout << abs_coeff;
            }
            
            for (size_t i = 0; i < m_var_names.size(); ++i)
            {
                if (degrees[i] != 0)
                {
                    std::cout << m_var_names[i];
                    if (degrees[i] != 1)
                    {
                        std::cout << '^' << degrees[i];
                    }
                    if (i < m_var_names.size() - 1) {std::cout << " * ";};

                }
            }
        }
        
        if (is_first)
        {
            std::cout << '0';
        }
        std::cout << std::endl;
    }


    // arithmetic operations
public:
    PolyTrie &operator+=(PolyTrie const &other)
    {
        std::function<void(TrieNode*, const TrieNode*, 
                       std::vector<int>&, size_t)> dfs_add = 
            [&](TrieNode* node1, const TrieNode* node2, 
                std::vector<int>& degrees, size_t depth)
            {
                if (depth == m_var_names.size())
                {
                    if (node1->coeff.has_value() && node2->coeff.has_value())
                    {
                        node1->coeff = node1->coeff.value() + node2->coeff.value();
                    }
                    else if(!(node1->coeff.has_value()) && node2->coeff.has_value())
                    {
                        node1->coeff = node2->coeff.value();
                    }
                    return;
                }

                std::set<int> all_degrees;
                for (auto const &[deg, child] : node1->childs) all_degrees.insert(deg);
                for (auto const &[deg, child] : node2->childs) all_degrees.insert(deg);

                for( auto deg : all_degrees)
                {
                    degrees[depth] = deg;
                    auto it1 = node1->childs.find(deg);
                    auto it2 = node2->childs.find(deg);
                    if (it1 == node1->childs.end())
                    {
                        node1->childs[deg] = std::make_unique<TrieNode>();
                        it1 = node1->childs.find(deg);
                    }

                    TrieNode *child2 = nullptr;
                    if (it2 != node2->childs.end())
                    {
                        child2 = it2->second.get();
                    }

                    dfs_add(it1->second.get(), child2, degrees, depth + 1);
                }
            };

        std::vector<int> degrees(m_var_names.size(), 0);
        dfs_add(m_root.get(), other.m_root.get(), degrees, 0);

        return *this;
    }

    PolyTrie operator+(PolyTrie const &other) const
    {
        PolyTrie obj = *this;
        return obj += other;
    }

    PolyTrie &operator-=(PolyTrie const &other)
    {
        std::function<void(TrieNode*, const TrieNode*, 
                       std::vector<int>&, size_t)> dfs_add = 
            [&](TrieNode* node1, const TrieNode* node2, 
                std::vector<int>& degrees, size_t depth)
            {
                if (depth == m_var_names.size())
                {
                    if (node1->coeff.has_value() && node2->coeff.has_value())
                    {
                        node1->coeff = node1->coeff.value() - node2->coeff.value();
                    }
                    else if(!(node1->coeff.has_value()) && node2->coeff.has_value())
                    {
                        node1->coeff = node2->coeff.value();
                    }
                    return;
                }

                std::set<int> all_degrees;
                for (auto const &[deg, child] : node1->childs) all_degrees.insert(deg);
                for (auto const &[deg, child] : node2->childs) all_degrees.insert(deg);

                for( auto deg : all_degrees)
                {
                    degrees[depth] = deg;
                    auto it1 = node1->childs.find(deg);
                    auto it2 = node2->childs.find(deg);
                    if (it1 == node1->childs.end())
                    {
                        node1->childs[deg] = std::make_unique<TrieNode>();
                        it1 = node1->childs.find(deg);
                    }

                    TrieNode *child2 = nullptr;
                    if (it2 != node2->childs.end())
                    {
                        child2 = it2->second.get();
                    }

                    dfs_add(it1->second.get(), child2, degrees, depth + 1);
                }
            };

        std::vector<int> degrees(m_var_names.size(), 0);
        dfs_add(m_root.get(), other.m_root.get(), degrees, 0);

        return *this;
    }

    PolyTrie operator-(PolyTrie const &other) const
    {
        PolyTrie obj = *this;
        return obj -= other;
    }

    PolyTrie &operator*=(PolyTrie const &other)
    {
        PolyTrie result(m_var_names);
        auto terms1 = get_supp();
        auto terms2 = other.get_supp();
        for (auto const &[deg1, coeff1] : terms1)
        {
            for(auto const &[deg2, coeff2] : terms2)
            {
                std::vector<int> new_degs(m_var_names.size());
                for (auto i = 0; i < m_var_names.size(); i++)
                {
                    new_degs[i] = deg1[i] + deg2[i];
                }
                result.add_term(new_degs, coeff1 * coeff2);
            }
        }
        *this = std::move(result);
        return *this;
    }

PolyTrie operator*(PolyTrie const &other) const
    {
        PolyTrie obj = *this;
        return obj *= other;
    }


    // equality operations

private:
    int lex_order_compare(std::vector<int> const &term1_degs, std::vector<int> const &term2_degs) const
    {
        for (auto i = 0; i < m_var_names.size(); i++)
        {
            if (term1_degs[i] > term2_degs[i]){return 1;}

            else if (term1_degs[i] < term2_degs[i]) {return -1;}
        }
        return 0;
    }

public:
    bool operator==(PolyTrie const &other) const
    {
        auto terms1 = get_supp();
        auto terms2 = other.get_supp();

        if (terms1.size() != terms2.size()) {return false;}

        for (auto const &[degs, coeff] : terms1)
        {
            auto it2 = terms2.find(degs);
            if (it2 == terms2.end()) {return false;}
            if (coeff != it2.value()) {return false;}
        }
        return true;
    }

    bool operator!=(PolyTrie const &other)
    {
        return !(*this == other);
    }


public:
    /*coeffType evaluate(std::vector<coeffType> const &point) const
    {
        std::function<coeffType(const TrieNode*, size_t )> eval_node = 
        [&](const TrieNode *node, size_t depth)
        {
            if (depth == m_var_names.size())
            {
                return (node->coeff.has_value() ? node->coeff.value() : coeffType(0));
            }

            coeffType x_i = point[depth];
            coeffType result = coeffType(0);
            for (auto rev_it = node->childs.rbegin(); rev_it != node->childs.rend(); rev_it++)
            {
                coeffType node_value = eval_node(rev_it->second.get(), depth + 1);
                result = result * x_i + node_value;
            }
            return result;
        };
        return eval_node(m_root.get(), 0);
    }*/

    coeffType evaluate(std::vector<coeffType> const &point) const
    {
        if (point.size() != m_var_names.size())
        {
            throw std::invalid_argument("Размерность точки не совпадает с числом переменных");
        }
        
        coeffType result = coeffType(0);
        auto terms = get_supp();
        
        for (auto const &[degrees, coeff] : terms)
        {
            coeffType term_value = coeff;
            for (size_t i = 0; i < m_var_names.size(); ++i)
            {
                if (degrees[i] != 0)
                {
                    term_value *= std::pow(point[i], degrees[i]);
                }
            }
            result += term_value;
        }
        return result;
    }

    int deg_of_uniformity() const
    {
        auto terms = get_supp();
        if (terms.empty())
        {
            return 0;
        }

        std::optional<int> max_deg;
        for (auto const &[degs, coeff] : terms)
        {
            int curent_total_deg = 0;
            for (auto deg : degs)
            {
                curent_total_deg += deg;
            }
            
            if (!max_deg.has_value()) 
            {
                max_deg = curent_total_deg;
            }
            else if (max_deg != curent_total_deg) 
            {
                return {};
            }
        }
        return max_deg.value();
    }

    PolyTrie homogeneous_part(int d) const
    {
        PolyTrie h(m_var_names);
        auto terms = get_supp();
        for (auto const &[degs, coeff] : terms)
        {
            int curent_f_term_deg = 0;
            for (auto deg : degs)
            {
                curent_f_term_deg += deg;
            }
            if (d == curent_f_term_deg)
            {
                h.add_term(degs, coeff);
            }
        }
        return h;
    }
};