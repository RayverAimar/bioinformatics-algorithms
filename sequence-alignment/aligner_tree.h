#include <iostream>
#include <fstream>
#include <queue>
#include <string>
#include <vector>

#define GAP '-'
#define INF ((1 << 16) - 1)

#define EXPORT // Define it if want to > file.txt

class AlignerNode
{
private:
public:
    AlignerNode();
    AlignerNode(const unsigned short int, const unsigned short int, const char = '$', const char = '$');
    ~AlignerNode();

    void go_for_diag(const char &, const char &);
    void go_for_top(const char &);
    void go_for_left(const char &);

    void branch_diag_top(const char &, const char &);
    void branch_diag_left(const char &, const char &);
    void branch_top_left(const char &, const char &);
    void branch_diag_top_left(const char &, const char &);

    void free_memory();

    bool is_leaf();

    std::string A, B;
    unsigned short int *i, *j;
    AlignerNode *children[3]; // 0 -> Diag, 1 -> Top, 2 -> Left
};

AlignerNode::AlignerNode()
{
    *i = INF;
    *j = INF;
    children[0] = children[1] = children[2] = nullptr;
}

AlignerNode::AlignerNode(const unsigned short int _i, const unsigned short int _j, const char _A, const char _B)
{
    i = new unsigned short int;
    j = new unsigned short int;
    
    *i = _i;
    *j = _j;
    A.push_back(_A);
    B.push_back(_B);
    children[0] = children[1] = children[2] = nullptr;
}

AlignerNode::~AlignerNode()
{
    for (int i = 0; i < 3; i++)
    {
        delete children[i];
    }
}

void AlignerNode::free_memory()
{
    i = j = nullptr;
}

void AlignerNode::go_for_diag(const char &_A, const char &_B)
{
    A.push_back(_A);
    B.push_back(_B);
    (*i)--, (*j)--;
}

void AlignerNode::go_for_top(const char &_A)
{
    A.push_back(_A);
    B.push_back(GAP);
    (*i)--;
}

void AlignerNode::go_for_left(const char &_B)
{
    A.push_back(GAP);
    B.push_back(_B);
    (*j)--;
}

void AlignerNode::branch_diag_top_left(const char &_A, const char &_B)
{
    children[0] = new AlignerNode(*i - 1, *j - 1, _A, _B); // Diagonal Branch
    children[1] = new AlignerNode(*i - 1, *j, _A, GAP);    // Top Branch
    children[2] = new AlignerNode(*i, *j - 1, GAP, _B);    // Left Branch
}

void AlignerNode::branch_diag_top(const char &_A, const char &_B)
{
    children[0] = new AlignerNode(*i - 1, *j - 1, _A, _B);
    children[1] = new AlignerNode(*i - 1, *j, _A, GAP);
}

void AlignerNode::branch_diag_left(const char &_A, const char &_B)
{
    children[0] = new AlignerNode(*i - 1, *j - 1, _A, _B);
    children[2] = new AlignerNode(*i, *j - 1, GAP, _B);
}

void AlignerNode::branch_top_left(const char &_A, const char &_B)
{
    children[1] = new AlignerNode(*i - 1, *j, _A, GAP);
    children[2] = new AlignerNode(*i, *j - 1, GAP, _B);
}

bool AlignerNode::is_leaf()
{
    return (children[0] == nullptr) && (children[1] == nullptr) && (children[2] == nullptr);
}

class AlignerTree
{
private:
    
public:
    AlignerTree(const unsigned short int, const unsigned short int);
    AlignerTree();
    ~AlignerTree();

    unsigned int get_solutions();
    bool leaves_are_aligned();
    void dfs(AlignerNode *&, std::string = "", std::string = "");
    void save(std::string);

    unsigned int solutions;
    unsigned int iterator = 0;
    std::queue<AlignerNode*> leaves;
    std::vector<std::string> possible_alignments;
    AlignerNode *root;
};

AlignerTree::AlignerTree()
{
    root = nullptr;
    solutions = 0;
    leaves.push(root);
}

AlignerTree::AlignerTree(const unsigned short int _i, const unsigned short int _j)
{
    root = new AlignerNode(_i, _j);
    solutions = 0;
    leaves.push(root);
}

AlignerTree::~AlignerTree()
{
    // delete root; // Not deleting whole tree :|
}

bool AlignerTree::leaves_are_aligned()
{
    return leaves.empty();
}

void AlignerTree::dfs(AlignerNode *&cur_node,
                      std::string Alignment_A,
                      std::string Alignment_B)
{
    if(cur_node->is_leaf())
    {
        Alignment_A  += cur_node->A;
        Alignment_B  += cur_node->B;

        #ifndef EXPORT
            /*for(int i = Alignment_A.size(); i > 0; i--)
            {
                std::cout << Alignment_A[i];
            }
            std::cout << ",";
            for(int i = Alignment_B.size(); i > 0; i--)
            {
                std::cout << Alignment_B[i];
            }
            std::cout << "\n";*/
        #else
            possible_alignments.push_back("");
            for(int i = Alignment_A.size(); i > 0; i--)
            {
                possible_alignments.back().push_back(Alignment_A[i]);
            }
            possible_alignments.push_back("");
            for(int i = Alignment_B.size(); i > 0; i--)
            {
                possible_alignments.back().push_back(Alignment_B[i]);
            }
        #endif

        return;
    }

    for (int i = 0; i < 3; i++)
    {
        if(cur_node->children[i])
            dfs(cur_node->children[i],
                Alignment_A + cur_node->A,
                Alignment_B + cur_node->B);
    }
}

unsigned int AlignerTree::get_solutions()
{
    return solutions;
}

void AlignerTree::save(std::string file_name)
{
    std::ofstream file(file_name);
    if(!file.is_open())
    {
        std::cerr << "Failed to open the file!" << std::endl;
    }
    for(int i = 0; i < possible_alignments.size(); i++)
    {
        for(int j = 1; j < possible_alignments[i].size(); j++)
        {
            file << possible_alignments[i][j];
        }
        file << "\n";
    }
    file.close();
    std::cout << "Data saved succesfully at  " << file_name << "!" << std::endl;
}