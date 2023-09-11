#include <chrono>
#include <iomanip>
#include <tuple>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#ifdef __LINUX__
#include <sys/types.h>
#include <sys/stat.h>
#else
#include <Windows.h>
#endif

#undef max
#define N_PROTEINS 3
#define LINES_PER_PROTEIN 19
#define GAP '-'

using namespace std;

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
    void save_alignments();
    void save_scoring_matrix();
    void save_trail_matrix();

private:
    void traceback_alignments(const int &,
                              const int &,
                              std::string = "",
                              std::string = "",
                              const int & = 0,
                              const int & = 0,
                              const int & = 0);
    std::vector<int> get_maximums();
    int penalty_per_gap(unsigned int, int, int);
    void init_matrices();
    unsigned int get_trail(const int &,
                           const int &,
                           const int &,
                           const int &);

    int **M;
    unsigned int **M_trail;
    unsigned int columns, rows;
    std::string A, B, project_name, folder, logs_path;
    std::string sequences_path, matrix_score_path, matrix_trail_path;
    std::ofstream file;
    int gap_penalty, match_reward, mismatch_penalty, first_gap_penalty;
    size_t total_gaps = 0, total_matchs = 0;
    size_t total_mismatchs = 0, total_alignments = 0;
    fsec alignment_duration, reconstruction_duration;
};

GlobalSequenceAlignerManager::GlobalSequenceAlignerManager()
{
    M = nullptr;
}

GlobalSequenceAlignerManager::GlobalSequenceAlignerManager(const std::string &_A,
                                                           const std::string &_B,
                                                           const std::string &_folder,
                                                           const std::string &_project_name,
                                                           const int _match_reward,
                                                           const int _mismatch_penalty,
                                                           const int _gap_penalty,
                                                           const int _f_gap_penalty) : 
    A(_A),
    B(_B),
    folder(_folder),
    project_name(_project_name),
    gap_penalty(_gap_penalty),
    match_reward(_match_reward),
    mismatch_penalty(_mismatch_penalty),
    first_gap_penalty(_f_gap_penalty),
    columns(A.size() + 1),
    rows(B.size() + 1)
{
    columns = A.size() + 1;
    rows = B.size() + 1;
    sequences_path = folder + project_name + "_sequences.txt";
    matrix_score_path = folder + project_name + "_matrix_score.csv";
    matrix_trail_path = folder + project_name + "_matrix_trail.csv";
    logs_path = folder + "logs.csv";
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
        M[i][0] = 0;
        M_trail[i][0] = DONE;
    }

    for (int i = 0; i < rows; i++)
    {
        M[0][i] = 0;
        M_trail[0][i] = DONE;
    }
}


void GlobalSequenceAlignerManager::save_scoring_matrix()
{
    file.open(matrix_score_path);
    file << "\\,-,";
    for(int i = 0; i < B.size(); i++)
    {
        file << B[i];
        if(i != B.size() - 1)
            file << ",";
    }
    file << "\n";
    for (int i = 0; i < columns; i++)
    {
        if (i == 0)
            file << "-,";
        else
            file << A[i - 1] << ",";
        for (int j = 0; j < rows; j++)
        {
            file << M[i][j];
            if (j != rows - 1)
                file << ",";
        }
        file << "\n";
    }
    file.close();
}

void GlobalSequenceAlignerManager::save_trail_matrix()
{
    file.open(matrix_trail_path);
    file << "\\,-,";
    for(int i = 0; i < B.size(); i++)
    {
        file << B[i];
        if(i != B.size() - 1)
            file << ",";
    }
    file << "\n";
    for (int i = 0; i < columns; i++)
    {
        if (i == 0)
            file << "-,";
        else
            file << A[i - 1] << ",";
        for (int j = 0; j < rows; j++)
        {
            file << M_trail[i][j];
            if (j != rows - 1)
                file << ",";
        }
        file << "\n";
    }
    file.close();
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
    if(max_value == 0) return DONE;
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
    auto start = Time::now();
    for (int i = 1; i < columns; i++)
    {
        for (int j = 1; j < rows; j++)
        {
            int match = M[i - 1][j - 1] +
                        ((A[i - 1] == B[j - 1]) ? match_reward : mismatch_penalty);
            int top = M[i - 1][j] + penalty_per_gap(TOP, i - 1, j);
            int left = M[i][j - 1] + penalty_per_gap(LEFT, i, j - 1);

            M[i][j] = std::max(match, std::max(top, std::max(left, 0)));
            M_trail[i][j] = this->get_trail(M[i][j], match, top, left);
        }
    }
    auto end = Time::now();
    alignment_duration = (end - start);
}

void GlobalSequenceAlignerManager::traceback_alignments(const int &i,
                                                        const int &j,
                                                        std::string a,
                                                        std::string b,
                                                        const int &matches,
                                                        const int &mismatches,
                                                        const int &gaps)
{
    if (M_trail[i][j] == DONE)
    {
        total_alignments++;
        total_gaps += gaps;
        total_mismatchs += mismatches;
        total_matchs += matches;
        std::reverse(a.begin(), a.end());
        std::reverse(b.begin(), b.end());
        file << a << "\n";
        file << b << "\n";
        return;
    }
    if (M_trail[i][j] & DIAG)
    {
        if (A[i - 1] == B[j - 1])
        {
            traceback_alignments(i - 1, j - 1, a + A[i - 1], b + B[j - 1],
                                 matches + 1, mismatches, gaps);
        }
        else
        {
            traceback_alignments(i - 1, j - 1, a + A[i - 1], b + B[j - 1],
                                 matches, mismatches + 1, gaps);
        }
    }
    if (M_trail[i][j] & TOP)
    {
        traceback_alignments(i - 1, j, a + A[i - 1], b + GAP,
                             matches, mismatches, gaps + 1);
    }
    if (M_trail[i][j] & LEFT)
    {
        traceback_alignments(i, j - 1, a + GAP, b + B[j - 1],
                             matches, mismatches, gaps + 1);
    }
}

std::vector<int> GlobalSequenceAlignerManager::get_maximums()
{
    int max_value = -1e9;

    for(int i = 0; i < columns; i++)
    {
        for(int j = 0; j < rows; j++)
        {
            if(M[i][j] > max_value)
            {
                max_value = M[i][j];
            }
        }
    }

    std::vector<int> maximums;

    for(int i = 0; i < columns; i++)
    {
        for (int j = 0; j < rows; j++)
        {
            if(M[i][j] == max_value)
            {
                maximums.push_back(i);
                maximums.push_back(j);
            }
        }
    }
    return maximums;
}

void GlobalSequenceAlignerManager::save_alignments()
{
    std::vector<int> maximums = get_maximums();
    file.open(sequences_path, std::ios::app);
    auto start = Time::now();
    for(int i = 0; i < maximums.size(); i+=2)
    {
        traceback_alignments(maximums[i], maximums[i+1]);
    }
    auto end = Time::now();
    file.close();
    std::cout << "Data saved succesfully at " << sequences_path << "!" << std::endl;
    reconstruction_duration = (end - start);

    file.open(logs_path, std::ios::app);
    file << project_name << ","
         << sequences_path << ","
         << M[maximums[0]][maximums[1]] << ","
         << total_alignments << ","
         //<< std::fixed << std::setprecision(20)
         << (float)alignment_duration.count() << ","
         << (float)reconstruction_duration.count() << ","
         << total_gaps / total_alignments << ","
         << total_matchs / total_alignments << ","
         << total_mismatchs / total_alignments << "\n";
    file.close();
}

class GlobalSequenceAligner
{
private:
    void make_folder(const std::string &folder_name);

public:
    GlobalSequenceAligner(const std::string &,
                          const std::string &,
                          const std::string &,
                          const std::string &,
                          const int = 1,
                          const int = -1,
                          const int = -2,
                          const int = -5
                        );
    void fit(int = 540);
    std::string A, B, directory;
    std::string A_name, B_name;
    int match_score, mismatch_score;
    int gap_penalty, first_gap_penalty;
};

GlobalSequenceAligner::GlobalSequenceAligner(const std::string &_A,
                                             const std::string &_B,
                                             const std::string &_A_name,
                                             const std::string &_B_name,
                                             const int _match_score,
                                             const int _mismatch_score,
                                             const int _gap_penalty,
                                             const int _first_gap_penalty
                                             ): 
    A(_A),
    B(_B),
    A_name(_A_name),
    B_name(_B_name),
    match_score(_match_score),
    mismatch_score(_mismatch_score),
    gap_penalty(_gap_penalty),
    first_gap_penalty(_first_gap_penalty)
{
    directory = A_name + "_x_" + B_name;
    make_folder(directory);
}

void GlobalSequenceAligner::make_folder(const std::string &folder_name)
{
#ifdef __LINUX__
    if (mkdir(folderName.c_str(), 0777) == 0)
    {
        std::cout << "\"" << folder_name << "\""
                  << " folder created successfully." << std::endl;
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
    wchar_t *wideFolderName = new wchar_t[bufferSize];
    MultiByteToWideChar(CP_UTF8, 0, folder_name.c_str(), -1, wideFolderName, bufferSize);

    if (CreateDirectoryW(wideFolderName, nullptr) || ERROR_ALREADY_EXISTS == GetLastError())
    {
        std::cout << "\"" << folder_name << "\""
                  << " folder created successfully." << std::endl;
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

    for (int i = 0; i < chunks; i++)
    {
        std::string project_name = to_string(i*chunk_size) +
                                   "_" +
                                   to_string((i * chunk_size) + chunk_size);
        GlobalSequenceAlignerManager temporalAligner(
            A.substr(i * chunk_size, chunk_size),
            B.substr(i * chunk_size, chunk_size),
            directory + "/",
            project_name,
            match_score,
            mismatch_score,
            gap_penalty,
            first_gap_penalty
        );
        temporalAligner.align();
        temporalAligner.save_alignments();
        temporalAligner.save_scoring_matrix();
        temporalAligner.save_trail_matrix();
    }

    if ((A.size() % chunk_size) == 0 && (B.size() % chunk_size) == 0)
        return; // No remainder

    // Aligning remaning of string
    std::string A_remainder = A.substr(chunks * chunk_size);
    std::string B_remainder = B.substr(chunks * chunk_size);
    std::string last_project_name = to_string(chunks * chunk_size) + "_END";
    GlobalSequenceAlignerManager ResidualAligner(A_remainder,
                                                 B_remainder,
                                                 directory + "/",
                                                 last_project_name,
                                                 match_score,
                                                 mismatch_score,
                                                 gap_penalty,
                                                 first_gap_penalty);
    ResidualAligner.align();
    ResidualAligner.save_alignments();
    ResidualAligner.save_scoring_matrix();
    ResidualAligner.save_trail_matrix();
}

typedef std::pair<std::string, std::string> SequenceType;

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

void align_sequences(const std::vector<SequenceType> &sequences)
{
    for (int i = 0; i < sequences.size(); i++)
    {
        for (int j = i; j < sequences.size(); j++)
        {
            if (i == j)
                continue;
            GlobalSequenceAligner Aligner(sequences[i].second,
                                          sequences[j].second,
                                          sequences[i].first,
                                          sequences[j].first,
                                          1, -1,-2, -5); // Current Alignment Configuration
            Aligner.fit(1080);
        }
    }
}

void export_sequences(const std::vector<SequenceType> &sequences)
{
    std::ofstream file("processed_sequences.txt");
    for(int i = 0; i < sequences.size(); i++)
    {
        file << sequences[i].first << "\n";
        file << sequences[i].second << "\n";
    }
}

int main()
{
    GlobalSequenceAligner Aligner("GGTTGACTA","TGTTACGG", "random_string_A", "random_string_B", 3,-3,-2,-2);
    Aligner.fit();
}