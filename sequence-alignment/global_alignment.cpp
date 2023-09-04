#include <chrono>
#include <iomanip>
#include <stack>
#include <tuple>
#include <iomanip>

#include "aligner_tree.h"

#define N_PROTEINS 3
#define LINES_PER_PROTEIN 19
#define MASK_LEFT 1 << 0
#define MASK_DIAG 1 << 1
#define MASK_TOP 1 << 2

using namespace std;

typedef std::tuple<int, int, int> IntergerTuple;
typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::duration<float> fsec;

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

class GlobalSequenceAligner
{
public:
    GlobalSequenceAligner();
    GlobalSequenceAligner(const std::string &,
                          const std::string &,
                          const int = 1,
                          const int = -1,
                          const int = -2);

    void align();
    void print_scoring_matrix(const int & = 4);
    void print_trail_matrix(const int & = 4);
    void save_best_alignments(std::string, std::string, const char = '-');

private:
    void init_matrices();
    unsigned int get_trail(const int &,
                           const int &,
                           const int &,
                           const int &);

    int **M;
    unsigned int **M_trail;
    std::string A, B;
    int gap_penalty, match_reward, mismatch_penalty;
    unsigned int columns, rows;
};

GlobalSequenceAligner::GlobalSequenceAligner()
{
    M = nullptr;
}

GlobalSequenceAligner::GlobalSequenceAligner(const std::string &_A,
                                             const std::string &_B,
                                             const int _match_reward,
                                             const int _mismatch_penalty,
                                             const int _gap_penalty) : 
    A(_A),
    B(_B),
    gap_penalty(_gap_penalty),
    match_reward(_match_reward),
    mismatch_penalty(_mismatch_penalty),
    columns(A.size() + 1),
    rows(B.size() + 1)
{
    this->init_matrices();
}

void GlobalSequenceAligner::init_matrices()
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

void GlobalSequenceAligner::print_scoring_matrix(const int &padding)
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

void GlobalSequenceAligner::print_trail_matrix(const int &padding)
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

unsigned int GlobalSequenceAligner::get_trail(const int &max_value,
                                              const int &match,
                                              const int &top,
                                              const int &left)
{
    return (match == max_value ? 1 : 0) + // << 0
           ((top == max_value ? 1 : 0) << 1) +
           ((left == max_value ? 1 : 0) << 2);
}

void GlobalSequenceAligner::align()
{
    auto start = Time::now();
    for (int i = 1; i < columns; i++)
    {
        for (int j = 1; j < rows; j++)
        {
            int match = M[i - 1][j - 1] +
                        ((A[i - 1] == B[j - 1]) ? match_reward : mismatch_penalty);
            int top = M[i - 1][j] + gap_penalty;
            int left = M[i][j - 1] + gap_penalty;

            M[i][j] = std::max(match, std::max(top, left));
            M_trail[i][j] = this->get_trail(M[i][j], match, top, left);
        }
    }
    auto end = Time::now();
    fsec duration = end - start;
    #ifndef EXPORT
        std::cout << "\n*** Alignment Finished ***\n" << std::endl;
        std::cout << "Elapsed time for alignment: " <<
            std::fixed << std::setprecision(10) <<
            duration.count() << " seconds" << std::endl;
    #endif
}

void GlobalSequenceAligner::save_best_alignments(std::string file_path, std::string log_path ,const char gap)
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
    AlignerTree *tree = new AlignerTree(A.size(), B.size());

    #ifndef EXPORT
        std::cout << "\n*** Getting Best Alignments ***" << std::endl;
    #endif
    auto start = Time::now();

    while (!tree->leaves_are_aligned())
    {
        int current_i = *(tree->leaves.front()->i);
        int current_j = *(tree->leaves.front()->j);

        if (M_trail[current_i][current_j] == DIAG)
        {
            tree->leaves.front()->go_for_diag(A[--current_i], B[--current_j]);
        }
        else if (M_trail[current_i][current_j] == TOP)
        {
            tree->leaves.front()->go_for_top(A[--current_i]);
        }
        else if (M_trail[current_i][current_j] == LEFT)
        {
            tree->leaves.front()->go_for_left(B[--current_j]);
        }
        else if (M_trail[current_i][current_j] == TOP_DIAG)
        {
            tree->leaves.front()->branch_diag_top(A[--current_i], B[--current_j]);
            tree->leaves.push(tree->leaves.front()->children[0]);
            tree->leaves.push(tree->leaves.front()->children[1]);
            tree->leaves.front()->free_memory();
            tree->leaves.pop();
        }
        else if (M_trail[current_i][current_j] == LEFT_DIAG)
        {
            tree->leaves.front()->branch_diag_left(A[--current_i], B[--current_j]);
            tree->leaves.push(tree->leaves.front()->children[0]);
            tree->leaves.push(tree->leaves.front()->children[2]);
            tree->leaves.front()->free_memory();
            tree->leaves.pop();
        }
        else if (M_trail[current_i][current_j] == TOP_LEFT)
        {
            tree->leaves.front()->branch_diag_top_left(A[--current_i], B[--current_j]);
            tree->leaves.push(tree->leaves.front()->children[1]);
            tree->leaves.push(tree->leaves.front()->children[2]);
            tree->leaves.front()->free_memory();
            tree->leaves.pop();
        }
        else if (M_trail[current_i][current_j] == TOP_LEFT_DIAG)
        {
            tree->leaves.front()->branch_diag_top_left(A[--current_i], B[--current_j]);
            tree->leaves.push(tree->leaves.front()->children[0]);
            tree->leaves.push(tree->leaves.front()->children[1]);
            tree->leaves.push(tree->leaves.front()->children[2]);
            tree->leaves.front()->free_memory();
            tree->leaves.pop();
        }
        else if (M_trail[current_i][current_j] == DONE)
        {
            tree->solutions++;
            tree->leaves.front()->free_memory();
            tree->leaves.pop();
        }
    }

    #ifndef EXPORT
        auto end = Time::now();
        fsec duration = end - start;
        std::cout << "\n*** Searching Finished ***\n" << std::endl;
        std::cout << "Elapsed time for searching best alignments: "
                << std::setprecision(10) << duration.count()
                << " seconds." << std::endl;
        std::cout << "There were found "
                << tree->get_solutions()
                << " possible solutions."
                << std::endl;

        std::cout << "\n*** Displaying possible Alignments ***\n" << std::endl;
    #else
        //std::cout << "A,B\n";
    #endif
    
    tree->dfs(tree->root);
    tree->save(file_path);

    std::ofstream file(log_path, std::ios::app);
    file << to_string(M[A.size()][B.size()]) << "," << std::to_string(tree->possible_alignments.size() / 2) << "\n";
    file.close();
    #ifndef EXPORT
        std::cout << "Maximum score for best alignments is: "
                << M[A.size()][B.size()] << std::endl;
    #endif
}

void test()
{
    std::fstream file;
    file.open("sequences.txt");

    if (!file.is_open())
    {
        std::cout << "Error while opening file" << std::endl;
    }

    std::vector<std::string> adn_strings;

    for (int i = 0; i < N_PROTEINS; i++)
    {
        std::string trimmed_string, temp;
        for (int j = 0; j < LINES_PER_PROTEIN; j++)
        {
            std::getline(file, temp);
            if (j == 0) continue;
            std::string cleaned_str = temp.substr(10, temp.size() - 10);
            for (int k = 0; k < (int(cleaned_str.size()) / 11) + 1; k++)
                trimmed_string.append(cleaned_str.substr(k * 11, 10));
        }
        adn_strings.push_back(trimmed_string);
    }
    file.close();

    int sequence_length = 60;

    /*
    adn_strings[0] -> Bacteria
    adn_strings[1] -> Sars-Cov
    adn_strings[2] -> Influenza
    */

    std::string output_path = "./data/sarscov_influenza_";
    std::string output_path_log = "./data/sarscov_influenza_logs.txt";

    std::ofstream header(output_path_log);
    header << "chunk,file_name,score,total_alignments\n";
    header.close();

    for(int i = 0; i < adn_strings[1].size(); i+=sequence_length)
    {
        std::string output_name = output_path + std::to_string(i) + "_" + std::to_string(i + sequence_length) + ".txt";
        std::ofstream log_file(output_path_log, std::ios::app);

        log_file << std::to_string(i) << "-" << std::to_string(i+sequence_length) << ",";
        log_file << output_name << ",";
        log_file.close();
        
        GlobalSequenceAligner Aligner(adn_strings[1].substr(i, sequence_length), adn_strings[2].substr(i, sequence_length));
        Aligner.align();
        Aligner.save_best_alignments(output_name, output_path_log);
    }

    //GlobalSequenceAligner Aligner(adn_strings[1].substr(0, sequence_length), adn_strings[2].substr(0, sequence_length));

    //Aligner.align();
    //Aligner.save_best_alignments();
}

void basic_test()
{
    std::string A = "TAGCATGTCTAGC";
    std::string B = "TAGCAGC";

    GlobalSequenceAligner Aligner(A, B);
    Aligner.align();
    Aligner.save_best_alignments("./data/output.txt", "data/log.txt");
    //Aligner.print_scoring_matrix();
    //Aligner.print_trail_matrix();
}

int main()
{
    //basic_test();
    test();
}