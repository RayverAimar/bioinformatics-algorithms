#include <iostream>
#include <utility>
#include <vector>

#include "./global_alignment.h"

using namespace std;

typedef std::pair<std::string, std::string> SequenceType;

class StarAligner
{
private:
public:
    StarAligner();
    StarAligner(const std::vector<SequenceType> &);

    int get_center_of_sequences();
    void fill_scores_matrix();
    void print_scores_matrix();
    void align();

    std::vector<SequenceType> Sequences;
    std::vector<std::vector<int>> scores_matrix;
};

StarAligner::StarAligner()
{
}

StarAligner::StarAligner(const std::vector<SequenceType> &_Sequences) : Sequences(_Sequences)
{
    scores_matrix = std::vector<std::vector<int>>(Sequences.size() + 1, std::vector<int>(Sequences.size() + 1, 0));
}

int StarAligner::get_center_of_sequences()
{
    fill_scores_matrix();
    int max_score = -1e9, idx;
    for (int i = 0; i < scores_matrix.size() - 1; i++)
    {
        if (max_score < scores_matrix[i][Sequences.size()])
        {
            max_score = scores_matrix[i][Sequences.size()];
            idx = i;
        }
    }
    return idx;
}

void StarAligner::print_scores_matrix()
{
    const int padding = 4;
    std::cout << "\n*** Displaying Scores Matrix ***\n"
              << std::endl;

    std::cout << std::setw(padding + 1) << " ";
    for (int i = 0; i < scores_matrix.size() - 1; i++)
        std::cout << std::setw(padding) << "S" + to_string(i + 1);
    std::cout << std::endl;
    for (int i = 0; i < scores_matrix.size(); i++)
    {
        for (int j = 0; j < scores_matrix[i].size(); j++)
        {
            if (j == 0 && i != scores_matrix.size() - 1)
                std::cout << std::setw(padding) << "S" << i + 1;
            if (i == scores_matrix.size() - 1 && j == 0)
                std::cout << std::setw(padding + 1) << " ";
            if (i == j)
                std::cout << std::setw(padding) << "\\";
            else
                std::cout << std::setw(padding) << scores_matrix[i][j];
        }
        std::cout << std::endl;
    }
}

void StarAligner::fill_scores_matrix()
{
    for (int i = 0; i < scores_matrix.size() - 1; i++)
    {
        for (int j = i + 1; j < scores_matrix[i].size() - 1; j++)
        {
            GlobalSequenceAlignerManager Aligner(Sequences[i].second, Sequences[j].second, 1, -1, -2, -2);
            Aligner.align();
            scores_matrix[i][j] = scores_matrix[j][i] = Aligner.get_score();
        }
    }

    for (int i = 0; i < scores_matrix.size() - 1; i++)
    {
        int score = 0;
        for (int j = 0; j < scores_matrix[i].size(); j++)
        {
            if (j == scores_matrix[i].size() - 1)
                scores_matrix[i][j] = scores_matrix[j][i] = score;
            else
                score += scores_matrix[i][j];
        }
    }
}

void StarAligner::align()
{
    int idx_center = get_center_of_sequences();
    for (int i = 0; i < Sequences.size(); i++)
    {
        if (i == idx_center)
            continue;
        GlobalSequenceAlignerManager Aligner(Sequences[idx_center].second, Sequences[i].second, 1, -1, -2, -2);
        Aligner.align();
        SequenceType alignment = Aligner.get_alignment();
        std::cout << alignment.first << "\n" << alignment.second << "\n\n";
    }
}

int main()
{
    std::vector<SequenceType> sequences = {std::make_pair("S1", "ATTGCCATT"), std::make_pair("S2", "ATGGCCATT"), std::make_pair("S3", "ATCCAATTTT"), std::make_pair("S4", "ATCTTCTT"), std::make_pair("S5", "ACTGACC")};
    StarAligner MultipleAligner(sequences);
    MultipleAligner.align();
    MultipleAligner.print_scores_matrix();

}