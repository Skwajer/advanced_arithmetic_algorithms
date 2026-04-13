#include <algorithm>
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
#include <type_traits>
#include <complex>

namespace order 
{
    struct Lex
    {
        bool operator()(const std::vector<int>& a,
                        const std::vector<int>& b) const
        {
            for (size_t i = 0; i < a.size(); i++)
            {
                if (a[i] != b[i])
                    return a[i] < b[i];
            }
            return false;
        }
    };

    struct InvLex
    {
        bool operator()(const std::vector<int>& a,
                        const std::vector<int>& b) const
        {
            for (size_t i = a.size() - 1; i > 0; i--)
            {
                if (a[i] != b[i])
                    return a[i] < b[i];
            }
            return false;
        }
    };

    struct GrLex
    {
        bool operator()(const std::vector<int>& a,
                        const std::vector<int>& b) const
        {
            int sum_deg_a = 0, sum_deg_b = 0;

            for (int x : a) sum_deg_a += x;
            for (int x : b) sum_deg_b += x;

            if (sum_deg_a != sum_deg_b)
                return sum_deg_a < sum_deg_b;

            for (size_t i = 0; i < a.size(); i++)
            {
                if (a[i] != b[i])
                    return a[i] < b[i];
            }
            return false;
        }
    };

    struct GrevLex
    {
        bool operator()(const std::vector<int>& a,
                        const std::vector<int>& b) const
        {
            int deg_a = 0, deg_b = 0;

            for (int x : a) deg_a += x;
            for (int x : b) deg_b += x;

            if (deg_a != deg_b)
                return deg_a < deg_b;

            for (int i = (int)a.size() - 1; i >= 0; i--)
            {
                if (a[i] != b[i])
                    return a[i] > b[i];
            }
            return false;
        }
    };
}

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
    PolyTrie(const PolyTrie& other)
    {
        m_var_names = other.m_var_names;
        m_root = clone(other.m_root.get());
    }

    PolyTrie& operator=(const PolyTrie& other)
    {
        if (this != &other)
        {
            m_var_names = other.m_var_names;
            m_root = clone(other.m_root.get());
        }
        return *this;
    }

    std::unique_ptr<TrieNode> clone(const TrieNode* node)
    {
        if (!node) return nullptr;

        auto new_node = std::make_unique<TrieNode>();

        if (node->coeff.has_value())
            new_node->coeff = node->coeff;

        for (const auto& [deg, child] : node->childs)
        {
            new_node->childs[deg] = clone(child.get());
        }

        return new_node;
    }

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
        if (coeff == coeffType(0)) {return;}
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
        
        bool is_negative = false;
        if constexpr (std::is_same_v<coeffType, std::complex<double>>) {
            is_negative = coeff.real() < 0 || (coeff.real() == 0 && coeff.imag() < 0);
        } else {
            is_negative = coeff < coeffType(0);
        }
        
        if (is_first)
        {
            if (is_negative) std::cout << '-';
            is_first = false;
        }
        else
        {
            std::cout << (is_negative ? " - " : " + ");
        }
        
        coeffType abs_coeff;
        if constexpr (std::is_same_v<coeffType, std::complex<double>>) 
        {
            abs_coeff = std::abs(coeff);
        } else {
            abs_coeff = is_negative ? -coeff : coeff;
        }
        
        bool has_variables = false;
        for (size_t i = 0; i < m_var_names.size(); ++i) {
            if (degrees[i] != 0) {
                has_variables = true;
                break;
            }
        }
        
        bool print_coeff = true;
        if constexpr (std::is_same_v<coeffType, std::complex<double>>) 
        {
            if (abs_coeff == coeffType(1) && has_variables) {
                print_coeff = false;
            }
        } else {
            if (abs_coeff == coeffType(1) && has_variables) {
                print_coeff = false;
            }
        }
        
        if (print_coeff) 
        {
            std::cout << abs_coeff;
        }
        
        bool first_var = true;
        for (size_t i = 0; i < m_var_names.size(); ++i) {
            if (degrees[i] != 0) {
                if (!first_var && !print_coeff) 
                {
                } else if (!first_var) {
                    std::cout << " * ";
                }
                std::cout << m_var_names[i];
                if (degrees[i] != 1) {
                    std::cout << '^' << degrees[i];
                }
                first_var = false;
            }
        }
    }
    
    if (is_first) {
        std::cout << '0';
    }
    std::cout << std::endl;
}


    // arithmetic operations
public:

    bool isZero() const 
    {
        return get_supp().empty();
    }

    PolyTrie operator*(const coeffType& scalar) const 
    {
        PolyTrie result(m_var_names);
        auto terms = get_supp();
        for (const auto& [degs, coeff] : terms) {
            result.add_term(degs, coeff * scalar);
        }
        return result;
    }

    PolyTrie &operator+=(PolyTrie const &other)
    {
        std::function<void(TrieNode*, const TrieNode*, size_t)> dfs_add =
            [&](TrieNode* node1, const TrieNode* node2, size_t depth)
        {
            if (node2 == nullptr)
                return;

            if (depth == m_var_names.size())
            {
                if (node2->coeff.has_value())
                {
                    if (node1->coeff.has_value())
                    {
                        node1->coeff = node1->coeff.value() + node2->coeff.value();
                    }
                    else
                    {
                        node1->coeff = node2->coeff.value();
                    }

                    if (node1->coeff.value() == coeffType(0))
                    {
                        node1->coeff.reset();
                    }
                }
                return;
            }

            std::set<int> degrees;

            for (auto const &[deg, _] : node2->childs)
                degrees.insert(deg);

            for (auto const &[deg, _] : node1->childs)
                degrees.insert(deg);

            for (int deg : degrees)
            {
                TrieNode* child1;
                auto it1 = node1->childs.find(deg);

                if (it1 == node1->childs.end())
                {
                    node1->childs[deg] = std::make_unique<TrieNode>();
                    child1 = node1->childs[deg].get();
                }
                else
                {
                    child1 = it1->second.get();
                }

                const TrieNode* child2 = nullptr;
                auto it2 = node2->childs.find(deg);
                if (it2 != node2->childs.end())
                {
                    child2 = it2->second.get();
                }

                dfs_add(child1, child2, depth + 1);
            }
        };

        dfs_add(m_root.get(), other.m_root.get(), 0);
        return *this;
    }

    PolyTrie operator+(PolyTrie const &other) const
    {
        PolyTrie obj = *this;
        return obj += other;
    }

    PolyTrie &operator-=(PolyTrie const &other)
{
    std::function<void(TrieNode*, const TrieNode*, size_t)> dfs_sub =
        [&](TrieNode* node1, const TrieNode* node2, size_t depth)
    {
        if (node2 == nullptr)
            return;

        if (depth == m_var_names.size())
        {
            if (node2->coeff.has_value())
            {
                if (node1->coeff.has_value())
                {
                    node1->coeff = node1->coeff.value() - node2->coeff.value();
                }
                else
                {
                    node1->coeff = -node2->coeff.value();
                }

                if (node1->coeff.value() == coeffType(0))
                {
                    node1->coeff.reset();
                }
            }
            return;
        }

        std::set<int> degrees;

        for (auto const &[deg, _] : node2->childs)
            degrees.insert(deg);

        for (auto const &[deg, _] : node1->childs)
            degrees.insert(deg);

        for (int deg : degrees)
        {
            TrieNode* child1;
            auto it1 = node1->childs.find(deg);

            if (it1 == node1->childs.end())
            {
                node1->childs[deg] = std::make_unique<TrieNode>();
                child1 = node1->childs[deg].get();
            }
            else
            {
                child1 = it1->second.get();
            }

            const TrieNode* child2 = nullptr;
            auto it2 = node2->childs.find(deg);
            if (it2 != node2->childs.end())
            {
                child2 = it2->second.get();
            }

            dfs_sub(child1, child2, depth + 1);
        }
    };

    dfs_sub(m_root.get(), other.m_root.get(), 0);
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

private:
bool is_devides(std::vector<int> const &a, std::vector<int> const &b)
{
    bool devides = true;
    for (int i = 0; i < a.size(); i++)
    {
        if (a[i] > b[i]) return false;
    }
    return true;
}

std::vector<int> devide_monoms(std::vector<int> divisible, std::vector<int> const &divisor)
{
    for (int i = 0; i < divisible.size(); i++)
    {
        divisible[i] -= divisor[i];
    }
    return divisible;
}

public:
template <typename Comparator>
std::pair<std::vector<PolyTrie<coeffType>>, PolyTrie<coeffType>> 
divide(
    std::vector<PolyTrie <coeffType> > &divisors,
    Comparator comp)
{
    std::vector<PolyTrie> qs;
    qs.reserve(divisors.size());

    for (size_t i = 0; i < divisors.size(); i++)
    {
        qs.emplace_back(m_var_names);
    }
    PolyTrie r(m_var_names);
    PolyTrie p(m_var_names);
    p = *this;

    while(!p.get_supp().empty())
    {
        bool devided = false;
        auto lead_monom_p = p.leading_monomial(comp);
        auto lead_coeff_p = p.leading_coeff(comp);

        for (int i = 0; i < divisors.size(); i++)
        {
            auto p_lead_monom = p.leading_monomial(comp);
            auto f_i_lead_monom = (divisors[i]).leading_monomial(comp);
            auto p_lead_coeff = p.leading_coeff(comp);
            auto f_i_lead_coeff = (divisors[i]).leading_coeff(comp);
            if (is_devides(f_i_lead_monom, p_lead_monom))
            {
                auto new_degs = devide_monoms(p_lead_monom, f_i_lead_monom);
                auto new_coeff = p_lead_coeff / f_i_lead_coeff;
                PolyTrie t(m_var_names);
                t.add_term(new_degs, new_coeff);
                qs[i] += t;
                p -= t * divisors[i];
                devided = true;
                break;
            }
        }

        if (!devided)
        {
            PolyTrie t(m_var_names);
            t.add_term(lead_monom_p, lead_coeff_p);
            r += t;
            p -= t;
        }
    }

    return std::make_pair(std::move(qs), std::move(r));;
}



public:
template<typename Comparator>
std::vector<int> find_multideg(Comparator comp) const
{
    auto terms = get_supp();
    std::vector<int> multideg(m_var_names.size());
    auto it = std::max_element(
            terms.begin(),
            terms.end(),
            [&](auto const &a, auto const &b)
            {
                return comp(a.first, b.first);
            });
    return it->first;
}

template<typename Comparator>
coeffType leading_coeff(Comparator comp) const
{
    auto terms = get_supp();
    auto deg = find_multideg(comp);
    return terms.at(deg);
}

template<typename Comparator>
std::vector<int> leading_monomial(Comparator comp) const
{
    return find_multideg(comp);
}

template<typename Comparator>
std::pair<std::vector<int>, coeffType> leading_term(Comparator comp) const
{
    return {find_multideg(comp), leading_coeff(comp)};
}


    // equality operations
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
            throw std::invalid_argument("point dimension doesn't match the number of variables");
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
    
    std::vector<int> lcm_monoms(const std::vector<int>& a, const std::vector<int>& b) const 
    {
        std::vector<int> res(a.size());
        for (size_t i = 0; i < a.size(); ++i) 
        {
            res[i] = std::max(a[i], b[i]);
        }
        return res;
    }

    template<typename Comparator>
    PolyTrie S_poly(const PolyTrie& g, Comparator comp) const 
    {

        auto lf = this->leading_monomial(comp);
        auto lg = g.leading_monomial(comp);
        auto cf = this->leading_coeff(comp);
        auto cg = g.leading_coeff(comp);

        auto lcm = lcm_monoms(lf, lg);


        std::vector<int> factor_f(lcm.size());
        std::vector<int> factor_g(lcm.size());
        for (size_t i = 0; i < lcm.size(); ++i) {
            factor_f[i] = lcm[i] - lf[i];
            factor_g[i] = lcm[i] - lg[i];
        }

        auto term1_coeff = coeffType(1) / cf;
        PolyTrie term1(this->m_var_names);
        term1.add_term(factor_f, term1_coeff);
        term1 *= *this;

        auto term2_coeff = coeffType(1) / cg;
        PolyTrie term2(this->m_var_names);
        term2.add_term(factor_g, term2_coeff);
        term2 *= g;

        return term1 - term2;
    }

public:
    template<typename Comparator>
    bool S_poly_reduces_to_zero(PolyTrie<coeffType>& gi, 
                                PolyTrie<coeffType>& gj,
                                std::vector<PolyTrie<coeffType>>& basis,
                                Comparator comp) 
    {
        auto S = gi.S_poly(gj, comp);
        auto [quotients, remainder] = S.divide(basis, comp);
        return remainder.isZero();
    }

    template<typename Comparator>
    bool isGroebnerBasis(std::vector<PolyTrie<coeffType>>& basis, Comparator comp) 
    {
        size_t s = basis.size();
        for (size_t i = 0; i < s; ++i) 
        {
            for (size_t j = i + 1; j < s; ++j) 
            {
                if (!S_poly_reduces_to_zero(basis[i], basis[j], basis, comp)) 
                {
                    return false;
                }
            }
        }
        return true;
    }
};