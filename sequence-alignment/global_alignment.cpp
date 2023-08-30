#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <iomanip>
#include <stack>

using namespace std;

typedef std::tuple<int, int, int> IntergerTuple;

enum Direction
{
    UNKNOWN,
    DONE,
    TOP,
    LEFT,
    DIAG,
    TOP_LEFT,
    TOP_DIAG,
    LEFT_DIAG,
    TOP_LEFT_DIAG = TOP + LEFT + DIAG,
};

class State
{
public:
    State();
    State(int, int);

    int value, came_from;

private:
};

State::State() : value(0), came_from(UNKNOWN)
{
}

State::State(int _value, int _came_from) : value(_value), came_from(_came_from)
{
}

class GlobalSequenceAligner
{
public:
    GlobalSequenceAligner();
    GlobalSequenceAligner(const std::string &, const std::string &, const int = 1, const int = -1, const int = -2);

    void align();
    void print_matrix(const int & = 4);
    void get_best_alignments(const char = '-');

private:
    void init_matrix();
    void copy_cur_state(std::vector<std::stack<char>> &,
                        std::vector<std::stack<char>> &,
                        std::vector<std::size_t> &,
                        std::vector<std::size_t> &,
                        const int &
                        );
    void push_nucleotids(std::stack<char> &,
                        std::stack<char> &,
                        const char,
                        const char
                        );
    State get_max_state(const int &, const int &, const int &);

    std::vector<std::vector<State>> M;
    std::string A, B;
    int gap_penalty, match_reward, mismatch_penalty;
};

GlobalSequenceAligner::GlobalSequenceAligner()
{
}

GlobalSequenceAligner::GlobalSequenceAligner(const std::string &_A, const std::string &_B,
                                             const int _match_reward, const int _mismatch_penalty,
                                             const int _gap_penalty) : A(_A),
                                                                       B(_B),
                                                                       gap_penalty(_gap_penalty),
                                                                       match_reward(_match_reward),
                                                                       mismatch_penalty(_mismatch_penalty)
{
    M = std::vector<std::vector<State>>(A.size() + 1, std::vector<State>(B.size() + 1));
    this->init_matrix();
}

void GlobalSequenceAligner::init_matrix()
{
    for (int i = 0; i < M.size(); i++)
    {
        M[i][0].value = gap_penalty * i;
        M[i][0].came_from = TOP;
    }
    for (int i = 0; i < M[0].size(); i++)
    {
        M[0][i].value = gap_penalty * i;
        M[0][i].came_from = LEFT;
    }

    M[0][0].came_from = DONE;
}

void GlobalSequenceAligner::print_matrix(const int &padding)
{
    std::cout << "\n*** Displaying Matrix ***\n" << std::endl;
    
    std::cout << setw(padding << 1) << "-";
    for (int i = 0; i < B.size(); i++)
        std::cout << std::setw(padding) << B[i];
    std::cout << std::endl;

    for (int i = 0; i < M.size(); i++)
    {
        if (i == 0)
            std::cout << std::setw(padding) << "-";
        else
            std::cout << std::setw(padding) << A[i - 1];
        for (int j = 0; j < M[i].size(); j++)
        {
            std::cout << std::setw(padding) << M[i][j].value;
        }
        std::cout << std::endl;
    }
}

State GlobalSequenceAligner::get_max_state(const int &match, const int &top_gap, const int &left_gap)
{
    int max_value = std::max(match, std::max(top_gap, left_gap));
    IntergerTuple direction_tuple = IntergerTuple((match - max_value) == 0 ? 1 : 0,
                                                  (top_gap - max_value) == 0 ? 1 : 0,
                                                  (left_gap - max_value) == 0 ? 1 : 0); // int(DIAG, TOP, LEFT);
    int came_from = DIAG * get<0>(direction_tuple) + TOP * get<1>(direction_tuple) + LEFT * get<2>(direction_tuple);
    return State(max_value, came_from);
}

void GlobalSequenceAligner::align()
{
    for (int i = 1; i < M.size(); i++)
    {
        for (int j = 1; j < M[i].size(); j++)
        {
            int match = M[i - 1][j - 1].value + ((A[i - 1] == B[j - 1]) ? 1 : -1);
            int top_gap = M[i - 1][j].value + gap_penalty;
            int left_gap = M[i][j - 1].value + gap_penalty;

            M[i][j] = this->get_max_state(match, top_gap, left_gap);
        }
    }
}

void GlobalSequenceAligner::copy_cur_state(std::vector<std::stack<char>> &AlignmentAs,
                                           std::vector<std::stack<char>> &AlignmentBs,
                                           std::vector<std::size_t> &indexAs,
                                           std::vector<std::size_t> &indexBs,
                                           const int &current_state_id)
{
    AlignmentAs.push_back(AlignmentAs[current_state_id]);
    AlignmentBs.push_back(AlignmentBs[current_state_id]);
    indexAs.push_back(indexAs[current_state_id]);
    indexBs.push_back(indexBs[current_state_id]);
}

void GlobalSequenceAligner::push_nucleotids(  std::stack<char> &A,
                            std::stack<char> &B,
                            const char to_push_A,
                            const char to_push_B
                        )
{
    A.push(to_push_A);
    B.push(to_push_B);
}

void GlobalSequenceAligner::get_best_alignments(const char gap)
{
    /*
        Solution
            -TCGAA
            ATCGTA

        How to Read Solution
            DIAG : Both Nucleotids are OK
            TOP  : Top Nucleotid gets a GAP (string B)
            LEFT : BOT Nucleotid gets a GAP (string A)

        Penalty for first GAP if there is more than one solution
            Penalty = 5

        Return
            Vector of alignments
    */
    std::vector<std::stack<char>> Alignments_A(1), Alignments_B(1);
    std::vector<size_t> indexAs = {A.size()};
    std::vector<size_t> indexBs = {B.size()};

    while (M[indexAs[indexAs.size() - 1]][indexBs[indexBs.size() - 1]].came_from != DONE)
    {
        const int current_alignments = Alignments_A.size();
        for (int i = 0; i < current_alignments; i++)
        {
            if (M[indexAs[i]][indexBs[i]].came_from == DIAG)
            {
                this->push_nucleotids(Alignments_A[i], Alignments_B[i], A[--indexAs[i]], B[--indexBs[i]]);
            }
            else if (M[indexAs[i]][indexBs[i]].came_from == TOP)
            {
                this->push_nucleotids(Alignments_A[i], Alignments_B[i], A[--indexAs[i]], gap);
            }
            else if (M[indexAs[i]][indexBs[i]].came_from == LEFT)
            {
                this->push_nucleotids(Alignments_A[i], Alignments_B[i], gap, B[--indexBs[i]]);
            }
            else if (M[indexAs[i]][indexBs[i]].came_from == TOP_LEFT)
            {
                this->copy_cur_state(Alignments_A, Alignments_B, indexAs, indexBs, i);
                int last_idx = Alignments_A.size() - 1;
                
                // TOP BRANCH
                this->push_nucleotids(Alignments_A[i], Alignments_B[i], A[--indexAs[i]], gap);

                // LEFT BRANCH
                this->push_nucleotids(Alignments_A[last_idx], Alignments_B[last_idx], gap, B[--indexBs[last_idx]]);
            }
            else if (M[indexAs[i]][indexBs[i]].came_from == TOP_DIAG)
            {
                this->copy_cur_state(Alignments_A, Alignments_B, indexAs, indexBs, i);
                int last_idx = Alignments_A.size() - 1;

                // TOP BRANCH
                this->push_nucleotids(Alignments_A[i], Alignments_B[i], A[--indexAs[i]], gap);

                // DIAG BRANCH
                this->push_nucleotids(Alignments_A[last_idx], Alignments_B[last_idx], A[--indexAs[last_idx]], B[--indexBs[last_idx]]);
            }
            else if (M[indexAs[i]][indexBs[i]].came_from == LEFT_DIAG)
            {
                this->copy_cur_state(Alignments_A, Alignments_B, indexAs, indexBs, i);
                int last_idx = Alignments_A.size() - 1;

                // LEFT BRANCH
                this->push_nucleotids(Alignments_A[i], Alignments_B[i], gap, B[--indexBs[i]]);

                // DIAG BRANCH
                this->push_nucleotids(Alignments_A[last_idx], Alignments_B[last_idx], A[--indexAs[last_idx]], B[--indexBs[last_idx]]);
            }
            else if (M[indexAs[i]][indexBs[i]].came_from == TOP_LEFT_DIAG)
            {
                this->copy_cur_state(Alignments_A, Alignments_B, indexAs, indexBs, i);
                this->copy_cur_state(Alignments_A, Alignments_B, indexAs, indexBs, i);
                int last_idx = Alignments_A.size() - 1;

                // TOP BRANCH
                this->push_nucleotids(Alignments_A[i], Alignments_B[i], A[--indexAs[i]], gap);

                // LEFT BRANCH
                this->push_nucleotids(Alignments_A[last_idx], Alignments_B[last_idx], gap, B[--indexBs[last_idx]]);

                // DIAG BRANCH
                this->push_nucleotids(Alignments_A[last_idx - 1], Alignments_B[last_idx - 1], A[--indexAs[last_idx - 1]], B[--indexBs[last_idx - 1]]);
            }
        }
    }

    std::cout << "\n*** Displaying Results ***\n" << std::endl;

    for (int i = 0; i < Alignments_A.size(); i++)
    {
        std::cout << "[" << i + 1 << "] ";
        do
        {
            std::cout << Alignments_B[i].top();
            Alignments_B[i].pop();
        } while (!Alignments_B[i].empty());
        
        std::cout << std::endl << setw(4) << " ";
        do{
            std::cout << Alignments_A[i].top();
            Alignments_A[i].pop();
        } while (!Alignments_A[i].empty());
        
        std::cout << std::endl;
    }
}

int main()
{
    std::string A = "AAAC";
    std::string B = "AGC";

    GlobalSequenceAligner Aligner(A, B);
    Aligner.align();
    Aligner.print_matrix();
    Aligner.get_best_alignments();
}