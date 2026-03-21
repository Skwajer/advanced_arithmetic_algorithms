#include <exception>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <map>
#include <functional>

template<typename coeffType>
class PolyTrie
{
private:
    struct TrieNode 
    {
        std::unordered_map<int, std::unique_ptr<TrieNode>> childs;
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
                    if (i != m_var_names.size() - 1) std::cout << " * ";

                }
            }
        }
        
        if (is_first)
        {
            std::cout << '0';
        }
        std::cout << std::endl;
    }

    PolyTrie &operator+=(PolyTrie const &other)
    {

        return *this;
    }

    PolyTrie &operator-=(PolyTrie const &other)
    {

        return *this;
    }
    PolyTrie &operator*=(PolyTrie const &other)
    {
        return *this;
    }
};