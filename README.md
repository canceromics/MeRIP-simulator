# MeRIP-simulator

A tool to simulate fasta sequence with m6A peak


## Table of Contents
* [Requirements](#Requirements)
* [Installation](#Installation)
* [Usage](#Usage)
* [Options](#Options)
* [License](#License)

## Requirements

* JDK 8

## Installation

* This tool can be installed by instructions as follows:

```
git clone https://github.com/canceromics/MeRIP-simulator.git
cd circm6a/m6asimulate/src
javac -d ./ -classpath ./lib/* ./sim/*.java ./genome/*.java ./note/*.java
jar -cvmf META-INF/MANIFEST.MF ../../m6asim.jar *
```
The tool is generated as m6asim.jar in this directory.

## Usage

* Start from genome sequence file and annotation file.

```
java -Xmx8g -jar m6asim.jar -g <genome.fa> -r <annotation.gtf> -p <peak.bed> -m <mutation.bed> -ip <ip.fastq> -in <input.fastq> [options]
```
Running this instruction will result in getting files named `peak.bed, mutation.bed, ip.fastq, input.fastq`. You can provide `peak.bed, mutation.bed` as well.

## Options

* Here are definitions of headers in output file named `(output_dir/file_prefix)_circRNAs.txt`

| Option       | Description                           |
| ---------- | ------------------------------------ |
| -bn | Background number (> 0) for simulation.  (10 default)|
| -en | Peak enrichment number (>= 0) for simulation. (15 default) |
| -sl | Read length in pair_end simulation, single-end mode enable when no greater than read length. (200 ~ 450 recommend) (350 default) |
| -rl | Alignment mate length (> 0) in simulation. (150 default) |
| -mrl | Minimum alignment mate length (> 0) in simulation. (5 default) |
| -rs | Random segment size in simulation. Therefore, simulated sequence length between alignment length ¡À random segment size. |
| -mp | Proportion (0.0 ~ 1.0) of mutation reads provided. (0.5 default)|
| -q | Qualities of bases simulated from this |


## License
Licensed GPLv3 for open source use or contact zuoLab (zuozhx@sysucc.org.cn) for commercial use.
