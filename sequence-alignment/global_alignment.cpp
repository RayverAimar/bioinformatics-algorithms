#include <chrono>
#include <iomanip>
#include <stack>
#include <tuple>
#include <iomanip>
#include <algorithm>

#include "aligner_tree.h"

#ifdef __LINUX__
    #include <sys/types.h>
    #include <sys/stat.h>
#else
    #include <Windows.h>
#endif

#undef max
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

class GlobalSequenceAlignerManager
{
public:
    GlobalSequenceAlignerManager();
    GlobalSequenceAlignerManager(const std::string &,
                          const std::string &,
                          const std::string &,
                          const std::string &,
                          const int = 1,
                          const int = -1,
                          const int = -2,
                          const int = -5);

    void align();
    void print_scoring_matrix(const int & = 4);
    void print_trail_matrix(const int & = 4);
    void save_best_alignments(std::string, std::string, const char = '-');
    void traceback_alignments(const int &, const int &, std::string = "", std::string = "");
    void save_alignments();

private:
    int penalty_per_gap(unsigned int, int, int);
    void init_matrices();
    unsigned int get_trail(const int &,
                           const int &,
                           const int &,
                           const int &);

    int **M;
    unsigned int **M_trail;
    std::string A, B, alignment_name, logs_path;
    std::ofstream file;
    int gap_penalty, match_reward, mismatch_penalty, first_gap_penalty;
    int total_gaps = 0, total_matchs = 0, total_mismatchs = 0;
    int total_alignments = 0;
    unsigned int columns, rows;
    fsec alignment_duration, reconstruction_duration;
};

GlobalSequenceAlignerManager::GlobalSequenceAlignerManager()
{
    M = nullptr;
}

GlobalSequenceAlignerManager::GlobalSequenceAlignerManager(const std::string &_A,
                                             const std::string &_B,
                                             const std::string &_alignment_name,
                                             const std::string &_logs_path,
                                             const int _match_reward,
                                             const int _mismatch_penalty,
                                             const int _gap_penalty,
                                             const int _f_gap_penalty) : 
    A(_A),
    B(_B),
    alignment_name(_alignment_name),
    gap_penalty(_gap_penalty),
    match_reward(_match_reward),
    mismatch_penalty(_mismatch_penalty),
    first_gap_penalty(_f_gap_penalty),
    logs_path(_logs_path),
    columns(A.size() + 1),
    rows(B.size() + 1)
{
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
    else return -5;
}

void GlobalSequenceAlignerManager::align()
{
    auto start = Time::now();
    for (int i = 1; i < columns; i++)
    {
        for (int j = 1; j < rows; j++)
        {
            int match = M[i - 1][j - 1] +
                        ((A[i - 1] == B[j - 1]) ? match_reward : mismatch_penalty);
            int top = M[i - 1][j] + penalty_per_gap(TOP, i-1, j);
            int left = M[i][j - 1] + penalty_per_gap(LEFT, i, j-1);

            M[i][j] = std::max(match, std::max(top, left));
            M_trail[i][j] = this->get_trail(M[i][j], match, top, left);
        }
    }
    auto end = Time::now();
    alignment_duration = end - start;
}

void GlobalSequenceAlignerManager::traceback_alignments(const int &i, const int &j, std::string a, std::string b)
{
    if(M_trail[i][j] == DONE)
    {
        total_alignments++;
        std::reverse(a.begin(), a.end());
        std::reverse(b.begin(), b.end());
        /*
        for(int idx = a.size()-1; idx >= 0; idx--)
            file << a[idx];
        file << "\n";
        for(int idx = b.size()-1; idx >= 0; idx--)
            file << b[idx];
        file << "\n";
        */
        file << a << "\n";
        file << b << "\n";
        return;
    }
    if(M_trail[i][j] & DIAG)
    {
        if(A[i-1] == B[j-1]) total_matchs++;
        else total_mismatchs++;
        traceback_alignments(i-1,j-1, a+A[i-1], b+B[j-1]);
    }
    if (M_trail[i][j] & TOP)
    {
        total_gaps++;
        traceback_alignments(i-1,j, a+A[i-1], b+GAP);
    }
    if(M_trail[i][j] & LEFT)
    {
        total_gaps++;
        traceback_alignments(i,j-1, a+GAP, b+B[j-1]);
    }
}

void GlobalSequenceAlignerManager::save_alignments()
{
    file.open(alignment_name, std::ios::app);
    auto start = Time::now();
    traceback_alignments(A.size(), B.size());
    auto end = Time::now();
    file.close();
    std::cout << "Data saved succesfully at  " << alignment_name << "!" << std::endl;
    reconstruction_duration = end - start;
    
    file.open(logs_path, std::ios::app);
    file << M[A.size()][B.size()] << "," 
         << total_alignments << ","
         << std::fixed << std::setprecision(10)
         << alignment_duration.count() << ","
         << reconstruction_duration.count() << ","
         << total_gaps/total_alignments<< ","
         << total_matchs/total_alignments << ","
         << total_mismatchs/total_alignments <<"\n";
    file.close();
}

void GlobalSequenceAlignerManager::save_best_alignments(std::string file_path,
                                                std::string log_path,
                                                const char gap)
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
    auto end = Time::now();
    reconstruction_duration = end - start;
    tree->dfs(tree->root);
    tree->save(file_path);

    std::ofstream file(log_path, std::ios::app);
    file << to_string(M[A.size()][B.size()]) << "," 
         << std::to_string(tree->possible_alignments.size() / 2) << ","
        << std::fixed << std::setprecision(10) << alignment_duration.count()
        << "," << reconstruction_duration.count() <<"\n";
    file.close();
}

class GlobalSequenceAligner
{
private:
    void make_folder(const std::string& folder_name);
public:
    GlobalSequenceAligner(const std::string&, const std::string&, const std::string&, const std::string&);
    void fit(int = 540);
    std::string A, B;
    std::string A_name, B_name;
    std::string directory;
};

GlobalSequenceAligner::GlobalSequenceAligner(const std::string &_A,
                                            const std::string &_B,
                                            const std::string &_A_name,
                                            const std::string &_B_name) :
    A(_A),
    B(_B),
    A_name(_A_name),
    B_name(_B_name)
{
    directory = A_name + "_x_" + B_name;
    make_folder(directory);
}

void GlobalSequenceAligner::make_folder(const std::string& folder_name)
{
    #ifdef __LINUX__
        if (mkdir(folderName.c_str(), 0777) == 0)
        {
            std::cout << "\"" << folder_name << "\"" <<" folder created successfully." << std::endl;
        }
        else
        {
            std::cerr << "Error creating folder." << std::endl;
        }
    #else
        int bufferSize = MultiByteToWideChar(CP_UTF8, 0, folder_name.c_str(), -1, nullptr, 0);
        if (bufferSize == 0)
        {
            std::cerr << "Error converting string to wide string." << std::endl;
        }
        wchar_t* wideFolderName = new wchar_t[bufferSize];
        MultiByteToWideChar(CP_UTF8, 0, folder_name.c_str(), -1, wideFolderName, bufferSize);

        if (CreateDirectoryW(wideFolderName, nullptr) || ERROR_ALREADY_EXISTS == GetLastError())
        {
            std::cout << "\"" << folder_name << "\"" <<" folder created successfully." << std::endl;
        }
        else
        {
            std::cerr << "Failed to create the folder. Error code: " << GetLastError() << std::endl;
        }

        delete[] wideFolderName;
    #endif
}

void GlobalSequenceAligner::fit(int chunk_size)
{
    /*
    if ((A.size() / chunk_size) != (B.size() / chunk_size))
    {
        std::cout << "\n[ERROR] Sequences sizes are to far apart. Exiting...\n" << std::endl;
        return;
        //Create a logic try to align all possble sequences (even if they are of different size)
    }
    */
    int chunks = A.size() / chunk_size;
    std::string logs_path = directory + "/logs.csv";
    std::ofstream logs(logs_path);
    logs << "chunk,file_name,score,total_alignments,alignment_time,reconstruction_time,mean_gaps,mean_match,mean_mismatch\n";
    logs.close();
    for(int i = 0; i < chunks; i++)
    {
        
        std::string chunk_file_name = directory + 
                                      "/" +
                                      to_string(i * chunk_size) +
                                      "_" +
                                      to_string((i * chunk_size) + chunk_size) +
                                      ".txt";
        std::ofstream log_file(logs_path, std::ios::app);
        log_file << i * chunk_size << "-" << (i * chunk_size) + chunk_size << ",";
        log_file << chunk_file_name << ",";
        log_file.close();

        std::string temporal_aligner_name = directory + 
                                            "/" +
                                            to_string(i * chunk_size) +
                                            "_" +
                                            to_string((i * chunk_size) + chunk_size) +
                                            ".txt";

        GlobalSequenceAlignerManager temporalAligner(
                A.substr(i * chunk_size, chunk_size),
                B.substr(i * chunk_size, chunk_size),
                temporal_aligner_name,
                logs_path
        );
        temporalAligner.align();
        temporalAligner.save_alignments();
        //temporalAligner.print_scoring_matrix();
        //temporalAligner.print_trail_matrix();
    }

    if((A.size() % chunk_size) == 0 && (B.size() % chunk_size) == 0) 
        return; // No remainder 

    // Aligning remaning of string
    std::string A_remainder = A.substr(chunks * chunk_size);
    std::string B_remainder = B.substr(chunks * chunk_size);
    std::string last_chunk_file_name = directory + 
                                            "/" +
                                            to_string(chunks * chunk_size) +
                                            "_" +
                                            to_string(A.size()) +
                                            ".txt";
    logs.open(logs_path, std::ios::app);
    logs << chunks * chunk_size << "," << A.size() << ",";
    logs << last_chunk_file_name << ",";
    logs.close();
    GlobalSequenceAlignerManager ResidualAligner(A_remainder,
                                                B_remainder,
                                                last_chunk_file_name,
                                                logs_path
    );
    ResidualAligner.align();
    ResidualAligner.save_alignments();
}



void basic_test()
{
    std::string A = "TAGCATGTCTAGC";
    std::string B = "TAGCAGC";

    GlobalSequenceAlignerManager Aligner(A, B, "Random_alignment", "logs.csv");
    Aligner.align();
    Aligner.save_best_alignments("./data/output.txt", "data/log.txt");
    Aligner.print_scoring_matrix();
    Aligner.print_trail_matrix();
}

void basic_test_2()
{
    std::string A = "caggtaacaa";
    std::string B = "gactatacta";
    std::string A_name = "random_string_1";
    std::string B_name = "random_string_2"; 
    GlobalSequenceAligner Aligner(A, B, A_name, B_name);
    Aligner.fit();
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

    /*
    adn_strings[0] -> Bacteria
    adn_strings[1] -> Sars-Cov
    adn_strings[2] -> Influenza
    */
    std::string A_name = "Sars_Cov";
    std::string B_name = "Bacteria"; 
    //Aligning Sars-Cov with Bacteria sequences of nucleotids
    GlobalSequenceAligner Aligner(adn_strings[1], adn_strings[2], A_name, B_name);
    Aligner.fit(1080);
}

int main()
{
    test();
}