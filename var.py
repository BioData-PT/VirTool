import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from jinja2 import Template
import numpy as np
import argparse
import pysam
import gzip
from Bio import SeqIO

def read_size(query_name):
    folder = args.run.split('/')
    reads_folder = ''.join(folder[:len(folder)-1])
    input_r1 = reads_folder + "/Results/reads/" + query_name + "_R1.fastq.gz" 
    input_r2 = reads_folder + "/Results/reads/" + query_name + "_R2.fastq.gz"
    count = 0
    size = 0
    with gzip.open(input_r1, 'rt') as f:
        for record in SeqIO.parse(f, 'fastq'):
            size += len(record.seq)
            count += 1

    with gzip.open(input_r2, 'rt') as f:
        for record in SeqIO.parse(f, 'fastq'):
            size += len(record.seq)
            count += 1

    read_length = round(float(size) / float(count), 2)    
    return read_length


def merge_check(query_name):
    input_file = args.run + "/" + query_name + "_merge.hist"
    with open(input_file, 'r') as f:
        size = [int(i.split('\t')[1].strip()) for i in f.readlines()]
    data = float(sum(size)) / float(len(size))
    if size:
        data = float(sum(size)) / float(len(size))
    else:
        data = 0

    return data


def genome_stats(ref, consensus):
    with open(ref, 'r') as f:
        f.readline()
        ref_len = len(f.readline())

    with open(consensus, 'r') as g:
        g.readline()
        seq = g.readline()
        ind = 0
        for i in seq:
            if i == 'N' or i == 'n' or i == '-':
                ind += 1

        seq_n = len(seq) - ind

    real_seq = [ref_len - seq_n, (float(seq_n)/float(ref_len))*100]
    return real_seq


def amplicons(amp, consensus):
    with open(amp, 'r') as g:
        amplicons = {}
        for line in g:
            line = line.split('\t')
            amplicons[line[0]] = [line[1], line[2].strip()]

    with open(consensus, 'r') as f:
        f.readline()
        f.readline()
        f.readline()
        seq = f.readline().strip('\n')
        insertion = []
        i = 1
        for base in seq:
            if base == 'n' or base == '-':
                insertion.append(i)
            i += 1
        last = 0
        ranges = []
        ranger = []
        for item in insertion:
            if item == last + 1:
                ranger.append(item)
                last = item
            else:
                ranges.append(ranger)
                ranger = []
                ranger.append(item)
                last = item
        ranges.append(ranger)

    return ranges


def parse_bam(args, positions):
    """Parse all positions from Bam pileup"""

    with open(args.ref, 'r') as f:
        f.readline()
        ref_seq = f.readline().strip()
    samfile = pysam.AlignmentFile(args.bam, 'rb')
    data = {}
    coverage = []

    for column in samfile.pileup(max_depth=60000, stepper='nofilter'):
        if (column.pos) in positions:
            position = {'a': 0,
                        'c': 0,
                        'g': 0,
                        't': 0
                        }

            for read in column.pileups:
                if not read.is_del and not read.is_refskip:
                    base = read.alignment.query_sequence[read.query_position]
                    position[base.lower()] += 1

            reads_num = sum([i for _, i in position.items()])
            coverage.append(reads_num)
            sort = sorted(position.items(), key=lambda x: x[1])

            if reads_num:
                allele_freq = (sort[-1][1] * 100.0) / reads_num
                allele_nuc = sort[-1][0]
                allele2_freq = (sort[-2][1] * 100.0) / reads_num
                allele2_nuc = sort[-2][0]
            else:
                allele_freq = 0
                allele_nuc = '-'
                allele2_freq = 0
                allele2_nuc = '-'

            data[column.reference_pos] = [reads_num,
                                          position['a'],
                                          position['c'],
                                          position['g'],
                                          position['t'],
                                          ref_seq[column.pos],
                                          allele_nuc,
                                          round(allele_freq, 2),
                                          allele2_nuc,
                                          round(allele2_freq, 2)]
    samfile.close()
    return data


def check_primers(position):
    with open(args.bed, 'r') as f:
        for line in f:
            line = line.split('\t')
            if int(position) > (int(line[1]) - 1) and int(position) < int(line[2]):
                data = line[3]
                return data


def get_ct(query_name):
    with open(args.info, 'r') as f:
        f.readline()
        for line in f:
            line = line.split('\t')
            if query_name == line[1]:
                value = line[2].strip()
                return value


def variation(ref, query):
    variation = {}
    for i in range(len(ref)):
        if ref[i] != query[i] and query[i] != '-' and query[i] != 'n':
            primer = check_primers(i)
            variation[i] = primer

    return variation


def process(args):
    list_var = []
    n_ranges = []

    # Get sample name, consensus....
    with open(args.query, 'r') as f:
        f.readline()
        ref = f.readline()
        query_name = f.readline().split('_')[1].split('.')[0]
        query = f.readline()

    # Get Ct Value
    ct_value = get_ct(query_name)

    # Get variation between reference and consensus
    var = variation(ref, query)

    # Check read depth and nucleotide composition in Variation sites
    report = parse_bam(args, var.keys())

    # Prepare data for Table1 html
    for key, value in report.items():
        anys = []
        anys.append(key + 1)
        for i in value:
            anys.append(i)
        anys.append(var[key])
        list_var.append(anys)

    

    # Finds NNN in consensus sequence
    ranges = amplicons(args.amp, args.query)

    # Prepare data to Table2 html
    for lista in ranges:
        if len(lista):
            if len(lista) == 1:
                n_ranges.append([lista[0], lista[0], len(lista)])
            else:
                n_ranges.append([lista[0], lista[-1], len(lista)])

    # Coverage by position
    samfile = pysam.AlignmentFile(args.bam, 'rb')
    pos_coverage = []
    for column in samfile.pileup(max_depth=60000, stepper='nofilter'):
        if column.n:
            pos_coverage.append([column.pos, column.n])
        else:
            pos_coverage.append([column.pos, '0'])
    coverage_total = 0
    for position in pos_coverage:
        coverage_total += position[1]
    coverage_mean = round(float(coverage_total) / float(29903), 2)

    # Set Table1
    header = ['Pos',
              'Reads',
              'A',
              'C',
              'G',
              'T',
              'Ref',
              'Base1',
              'Freq1',
              'Base2',
              'Freq2',
              'Primer']
    if list_var:
        table1 = pd.DataFrame(list_var)
        table1.columns = header
        html_table1 = table1.to_html(index=False, justify='center')
    else:
        html_table1 = "<p> NO VARIATION </p>"

    # Set Table2
    table2_header = ['Beg', 'End', 'Len']
    if n_ranges:
        table2 = pd.DataFrame(n_ranges)
        table2.columns = table2_header
        html_table2 = table2.to_html(index=False, justify='center')

    # Get Genome stats
    percent = genome_stats(args.ref, args.cons)
    missing_bases = percent[0]
    genome_percent = round(percent[1], 2)

    # Save Coverage Plot to file
    if pos_coverage:
        df = pd.DataFrame(pos_coverage)
        df.columns = ['Position', 'Num']
        plt.figure(figsize=(10, 4))
        plt.title(query_name)
        plt.bar(df['Position'], df['Num'], width=1.2, linewidth=0)
        plt.ylim(ymax=300)
        plt.xlim(xmin=0, xmax=29903)
        plt.xticks(np.arange(0, 29903, 1000))
        plt.xticks(fontsize=5, rotation='vertical')
        plt.yticks(fontsize=5)
        # plt.tight_layout()
        plt.savefig(args.run + '/' + query_name + '.png')  # , quality=100)
        # plt.show()
    else:
        with open(args.run + '/' + query_name + '.png', 'w') as w:
            w.write("  ")
    
    # Get mean read length
    read_length = read_size(query_name)

    #Get merge info
    merge_info = round(merge_check(query_name), 2)

    
    # Print Report txt
    print('{}\t{}\t{}\t{}\t{}\t{}'.format(query_name, coverage_mean, genome_percent, missing_bases, read_length, merge_info))



    # Output template
    output = Template("""<!DOCTYPE html>
    <html lang="en">
    <head>
    <h1> {{ query_name }} </h1>
    <h3> Ct Value - {{ ct_value }} </h3>
    <h2> Stats </h2>
    <h3> Missing Bases- {{ missing_bases }} </h3>
    <h3> Genome - {{ genome_percent }}% </h3>
    <h3> Coverage Mean - {{ coverage_mean }} </h3>
    </head>
    <body>
    <center>
    <h2> Reads/Positions </h2>
    <img src={{ image_file }}>
    <h2> Variation </h2>
    <table align='center'>{{ html_table1 }}</table>
    <h2> Gaps Analysis </h2>
    <table align='center'> {{ html_table2 }} </table>
    </body>
    </center>
    </html>
    <br>
    <br>
    <p style='page-break-before: always'>\n""")

    image_file = query_name + '.png'
    msg = output.render(query_name=query_name,
                        ct_value=ct_value,
                        missing_bases=missing_bases,
                        genome_percent=genome_percent,
                        coverage_mean=coverage_mean,
                        image_file=image_file,
                        html_table1=html_table1,
                        html_table2=html_table2
                        )
    with open(args.run + '/' + query_name + 'OUTPUT.html', 'w') as f:
        f.write(msg)


def arg_parse():
    parser = argparse.ArgumentParser(description="Format Variation output")
    parser.add_argument('run', help="Run Folder")
    parser.add_argument('ref', help="Reference File")
    parser.add_argument('query', help='Alignment File')
    parser.add_argument('info', help="Sample Info")
    parser.add_argument('bam', help="Bam for nucleotide count")
    parser.add_argument('bed', help="Bed file with primers")
    parser.add_argument('amp', help="Bed file with amplicons")
    parser.add_argument('cons', help="Consensus from iVar")

    return parser.parse_args()


if __name__ == '__main__':
    args = arg_parse()
    process(args)
