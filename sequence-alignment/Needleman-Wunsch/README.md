<h1 align="center">Needleman-Wunsch Algorithm </h1>

The current repository has an scalable implementation of the [Needleman-Wunsch algorithm](https://en.wikipedia.org/wiki/Needleman–Wunsch_algorithm), is well used it may give statistics of different alignments configurations, like reconstruction time, mean gaps per alignemtns, mean mismatches or mean matches, it can be later analyzed to see which configuration is the best for your DNA sequences.

## Installation
Make sure you've got a modern g++ version in your environment PATH. Also make sure you've got the python libraries used in this project, if not run the following command:

```bash
$ pip install -r requirements.txt
```

## Usage
To use the default alignment configuration just run any of the created scripts either in cmd or powershell consoles.

### CMD console
```bash
$ ./trigger.bat
```

### Powershell console
```bash
$ .\trigger.ps1
```
It will compile and run the cpp file and the current py dot_plot_maker. Once done, folders with the names of the differences sequences allocated in the [sequences.txt](https://github.com/RayverAimar/bioinformatics-algorithms/blob/master/sequence-alignment/Needleman-Wunsch/sequences.txt) file will appear with their proper trail_matrix, score_matrix and optimal sequences.

Also in the [/dot_plots](https://github.com/RayverAimar/bioinformatics-algorithms/blob/master/sequence-alignment/Needleman-Wunsch/dot_plots) folder will appear the possible dot plots per each sequence.

If you run the project with different configurations you will likely see how the possible alignments decrease or increase due to the "initial rupture penalty", you can configure it when calling the GlobalSequenceAligner, here is an example:

```cpp
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
                                          1, -1,-2, -5); // Current Alignment Configuration [Change configuration]
            Aligner.fit(1080);
        }
    }
}
```
Where the first highlighter parameter indicates the score when two nucleotids are aligned, the second when are not aligned, the third indicates the penalty per rupture and the fourth the penalty per first rupture.
The current configuration will give optimal answers but a lot of them will have few gaps and a lot of missmatched nucleotids.


The following configuration will give better results but it will take more time as there will be a lot of similar results.
```cpp
...
            GlobalSequenceAligner Aligner(sequences[i].second,
                                          sequences[j].second,
                                          sequences[i].first,
                                          sequences[j].first,
                                          1, -1,-2, -2); // Current Alignment Configuration [Change configuration]
            Aligner.fit(40);
...
```

## Analysis

Each time a folder is created, there will be a logs.csv file where some statitics per alignment are store, the last configuration (chunksize of 40 with (1,-1,-2,-2)) gave the following partial results.

| score | total_alignments | alignment_time | reconstruction_time | mean_gaps | mean_match | mean_mismatch |
|-------|------------------|----------------|---------------------|-----------|------------|---------------|
| -11   | 32               | 4.01e-05       | 0.0001074           | 6         | 19         | 18            |
| -15   | 432              | 3.79e-05       | 0.0012144           | 9         | 19         | 15            |
| -10   | 6                | 3.87e-05       | 5.89e-05            | 4         | 18         | 20            |
| -12   | 864              | 4.34e-05       | 0.0036288           | 8         | 20         | 16            |
| -2    | 18               | 5.4e-05        | 0.0001302           | 4         | 22         | 16            |

Meanwhile if the configuration changes to (1,-1,-2,-5) the results will be:

| score | total_alignments | alignment_time | reconstruction_time | mean_gaps | mean_match | mean_mismatch |
|-------|------------------|----------------|---------------------|-----------|------------|---------------|
| -16   | 1                | 4.77e-05       | 1.42e-05            | 0         | 12         | 28            |
| -20   | 1                | 6.69e-05       | 2.15e-05            | 0         | 10         | 30            |
| -16   | 2                | 4.87e-05       | 3.04e-05            | 1         | 13         | 26            |
| -20   | 1                | 4.59e-05       | 1.98e-05            | 0         | 10         | 30            |
| -18   | 1                | 5.19e-05       | 2.26e-05            | 4         | 17         | 21            |

If results are analyzed the mean missmatch in the second configuration is greater, mean gaps are zero in different solutions. The second configuration also give alignments with worse score but alignments remain to be optimal. Second configuration tends to be quicker than the first configuration as it only finds 1 or 2 optimal alignments meanwhile first configuration finds 864 alignments in sequences of length of 40 nucleotids, which is not scalable if trying to align sequences of greater lenth. So you can analyze and decide which configuration will likely fit to your sequences.