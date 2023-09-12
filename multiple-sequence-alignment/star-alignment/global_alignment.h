#include <chrono>
#include <iomanip>
#include <tuple>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define N_PROTEINS 3
#define LINES_PER_PROTEIN 19
#define GAP '-'

using namespace std;

typedef std::pair<std::string, std::string> SequenceType;

enum Direction
{
    UNKNOWN,
    DIAG = 1 << 0,
    TOP = 1 << 1,
    LEFT = 1 << 2,
    TOP_DIAG = TOP + DIAG,
    LEFT_DIAG = LEFT + DIAG,
    TOP_LEFT = TOP + LEFT,
    TOP_LEFT_DIAG = DIAG + TOP + LEFT,
    DONE = 99,
};

class GlobalSequenceAlignerManager
{
public:
    GlobalSequenceAlignerManager();
    GlobalSequenceAlignerManager(const std::string &,
                                 const std::string &,
                                 const int = 1,
                                 const int = -1,
                                 const int = -2,
                                 const int = -5);

    void align();
    void print_scoring_matrix(const int & = 4);
    void print_trail_matrix(const int & = 4);

    int get_score();
    SequenceType get_alignment();

private:
    SequenceType get_alignment(const int &,
                               const int &,
                               std::string = "",
                               std::string = "");
    int penalty_per_gap(unsigned int, int, int);
    void init_matrices();
    unsigned int get_trail(const int &,
                           const int &,
                           const int &,
                           const int &);

    int **M;
    unsigned int **M_trail;
    unsigned int columns, rows;
    std::string A, B;
    int gap_penalty, match_reward, mismatch_penalty, first_gap_penalty;
};

GlobalSequenceAlignerManager::GlobalSequenceAlignerManager()
{
    M = nullptr;
}

GlobalSequenceAlignerManager::GlobalSequenceAlignerManager(const std::string &_A,
                                                           const std::string &_B,
                                                           const int _match_reward,
                                                           const int _mismatch_penalty,
                                                           const int _gap_penalty,
                                                           const int _f_gap_penalty) : A(_A),
                                                                                       B(_B),
                                                                                       gap_penalty(_gap_penalty),
                                                                                       match_reward(_match_reward),
                                                                                       mismatch_penalty(_mismatch_penalty),
                                                                                       first_gap_penalty(_f_gap_penalty),
                                                                                       columns(A.size() + 1),
                                                                                       rows(B.size() + 1)
{
    columns = A.size() + 1;
    rows = B.size() + 1;
    this->init_matrices();
}

void GlobalSequenceAlignerManager::init_matrices()
{
    /* Allocating memory for matrices */

    M = new int *[columns];
    M_trail = new unsigned int *[columns];
    for (int i = 0; i < columns; i++)
    {
        M[i] = new int[rows];
        M_trail[i] = new unsigned int[rows];
    }

    /* Initializing matrices */

    for (int i = 0; i < columns; i++)
    {
        for (int j = 0; j < rows; j++)
        {
            M_trail[i][j] = UNKNOWN;
        }
    }

    for (int i = 0; i < columns; i++)
    {
        M[i][0] = i * gap_penalty;
        M_trail[i][0] = TOP;
    }

    for (int i = 0; i < rows; i++)
    {
        M[0][i] = i * gap_penalty;
        M_trail[0][i] = LEFT;
    }

    M_trail[0][0] = DONE;
}

void GlobalSequenceAlignerManager::print_scoring_matrix(const int &padding)
{
    std::cout << "\n*** Displaying Scoring Matrix ***\n"
              << std::endl;
    std::cout << std::setw(padding << 1) << "-";
    for (int i = 0; i < B.size(); i++)
        std::cout << std::setw(padding) << B[i];
    std::cout << std::endl;

    for (int i = 0; i < columns; i++)
    {
        if (i == 0)
            std::cout << std::setw(padding) << "-";
        else
            std::cout << std::setw(padding) << A[i - 1];
        for (int j = 0; j < rows; j++)
        {
            std::cout << std::setw(padding) << M[i][j];
        }
        std::cout << std::endl;
    }
}

void GlobalSequenceAlignerManager::print_trail_matrix(const int &padding)
{
    std::cout << "\n*** Displaying Trail Matrix ***\n"
              << std::endl;

    std::cout << std::setw(padding << 1) << "-";
    for (int i = 0; i < B.size(); i++)
        std::cout << std::setw(padding) << B[i];
    std::cout << std::endl;

    for (int i = 0; i < columns; i++)
    {
        if (i == 0)
            std::cout << std::setw(padding) << "-";
        else
            std::cout << std::setw(padding) << A[i - 1];
        for (int j = 0; j < rows; j++)
        {
            std::cout << std::setw(padding) << M_trail[i][j];
        }
        std::cout << std::endl;
    }
}

unsigned int GlobalSequenceAlignerManager::get_trail(const int &max_value,
                                                     const int &match,
                                                     const int &top,
                                                     const int &left)
{
    return (match == max_value ? 1 : 0) + // << 0
           ((top == max_value ? 1 : 0) << 1) +
           ((left == max_value ? 1 : 0) << 2);
}

int GlobalSequenceAlignerManager::penalty_per_gap(unsigned int direction, int i, int j)
{
    if (M_trail[i][j] == direction)
        return gap_penalty;
    else
        return first_gap_penalty;
}

void GlobalSequenceAlignerManager::align()
{
    for (int i = 1; i < columns; i++)
    {
        for (int j = 1; j < rows; j++)
        {
            int match = M[i - 1][j - 1] +
                        ((A[i - 1] == B[j - 1]) ? match_reward : mismatch_penalty);
            int top = M[i - 1][j] + penalty_per_gap(TOP, i - 1, j);
            int left = M[i][j - 1] + penalty_per_gap(LEFT, i, j - 1);

            M[i][j] = std::max(match, std::max(top, left));
            M_trail[i][j] = this->get_trail(M[i][j], match, top, left);
        }
    }
}

int GlobalSequenceAlignerManager::get_score()
{
    return M[A.size()][B.size()];
}

SequenceType GlobalSequenceAlignerManager::get_alignment()
{
    return get_alignment(A.size(), B.size());
}

SequenceType GlobalSequenceAlignerManager::get_alignment(const int &i,
                                                         const int &j,
                                                         std::string a,
                                                         std::string b)
{
    if (M_trail[i][j] == DONE)
    {
        std::reverse(a.begin(), a.end());
        std::reverse(b.begin(), b.end());
        return std::make_pair(a, b);
    }
    if (M_trail[i][j] & DIAG)
    {
        if (A[i - 1] == B[j - 1])
        {
            return get_alignment(i - 1, j - 1, a + A[i - 1], b + B[j - 1]);
        }
        else
        {
            return get_alignment(i - 1, j - 1, a + A[i - 1], b + B[j - 1]);
        }
    }
    if (M_trail[i][j] & TOP)
    {
        return get_alignment(i - 1, j, a + A[i - 1], b + GAP);
    }
    //if (M_trail[i][j] & LEFT)
    //{
    //    return get_alignment(i, j - 1, a + GAP, b + B[j - 1]);
    //}
    return get_alignment(i, j - 1, a + GAP, b + B[j - 1]);
}

std::vector<SequenceType> process_txt_sequences()
{
    std::fstream file;

    std::vector<SequenceType> sequences;

    file.open("sequences.txt");
    if (!file.is_open())
    {
        std::cout << "Error while opening file" << std::endl;
    }

    for (int i = 0; i < N_PROTEINS; i++)
    {
        SequenceType current_sequence;
        std::string trimmed_string, temp;
        for (int j = 0; j < LINES_PER_PROTEIN; j++)
        {
            std::getline(file, temp);
            if (j == 0)
            {
                current_sequence.first = temp;
                continue;
            }
            std::string cleaned_str = temp.substr(10, temp.size() - 10);
            for (int k = 0; k < (int(cleaned_str.size()) / 11) + 1; k++)
                trimmed_string.append(cleaned_str.substr(k * 11, 10));
        }
        current_sequence.second = trimmed_string;
        sequences.push_back(current_sequence);
    }
    file.close();
    return sequences;
}