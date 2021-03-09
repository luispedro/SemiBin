# S<sup>3</sup>N<sup>2</sup>Bin

S<sup>3</sup>N<sup>2</sup>Bin is a command line tool for metagenomic binning with semi-supervised siamese neural network using additional information from reference genomes. It will output the reconstructed bins in single sample/co-assembly/multi-samples binning mode.

### Single sample binning

- Single sample assembly.
- Mapping reads to contig fasta file.
- Binning.

### Co-assembly binning

- Co-assembly.
- Mapping reads from several samples to the co-assembly fasta file.
- Binning.

### Multi-samples binning

- Concatenate all contigs from all samples together. 
- Mapping reads from several samples to the concatenated fasta file.
- Binning for every sample.

## Commands

### Basic commands

* `-i/--input-fasta` : Path to the input contig fasta file(.fasta/.gz/.bz2).

* `-b/--input-bam`: Path to the input bam files. If there are several samples, you can input several bam files.

* `-c/--cannot-link`:  Path to the input cannot-link file generated from other additional biological information. One row for each cannot-link constraint.The file format: contig_1,contig_2.

* `-o/--output`: Output directory (will be created if non-existent).

* `-s/--separator`: Used when multi-samples binning to separete sample name and contig name.(None means using single sample and co-assemble binning)

### Optional commands

* `-p/--processes`: Number of subprocess used in processing bam files(default:whole).

* `--epoches`: Number of epoches used in the training process(default:20).

* `--batch-size`: Batch size used in the training process.(default:2048).

* `--max-edges`: The maximun number of edges that can be connected to one contig(default:200).

* `--max-node`: Percentage of contigs that considered to be binned(default:1).

	
