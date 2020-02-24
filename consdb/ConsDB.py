import argparse
import bz2
import ftplib
try:
    import manhole
except ModuleNotFoundError:
    pass
import multiprocessing as mp
import numpy as np
import os
import re
import RSEntry as rse
import SlimRSCollection as srsc
from string import ascii_lowercase
import sys
import time
import traceback
import urllib.error
import urllib.request

# List of all chromosomes
CHRS = [str(i) for i in range(1,23)] + ['X', 'Y', 'MT']

def check_consdb_sorted(fn):
    """
    Helper function to check if a ConsDB file is sorted by position.

    Parameters:
    fn: File to check
    """

    chrom = 0
    cur = 0
    for line in rse.RSCollection.open(fn):
        try:
            line = line.decode()
        except AttributeError:
            pass

        c, p = re.split('[:,]', line)[1:3]
        c = rse.RSCollection.chrom_to_int(c)
        p = int(p)

        if c > chrom:
            chrom = c
            cur = 0
        if p < cur or c < chrom:
            return(False)
        cur = p

    return(True)

def check_downloads(db, chrs, fp='.', force=False):
    """
    Function to download all necessary database files. Valid databases are:
    * dbSNP
    * 1000gp
    * gnomAD
    Database names are not case sensitive.

    Parameters:
    db: Database to use
    chrs: List of chromosomes to use
    fp: File path to store files
    force: Whether or not to re-download files that are already downloaded
    """

    db = db.lower()
    if db == 'dbsnp':
        file_base = ['refsnp-chr{}.json.bz2']
        url_base = 'ftp.ncbi.nlm.nih.gov/snp/latest_release/JSON/{}'
        ext = re.escape('.json.bz2')
    elif db == '1000gp':
        # Also need to download index file and pop file
        fns = ['1000genomes.sequence.index', '20131219.populations.tsv']
        urls = [('ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/'
            '1000_genomes_project/{}'),
            'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/{}']
        for i in range(len(fns)):
            fn = fns[i]
            fn_full = f'{fp}/{fn}'
            print(f'Checking for {fn}...', flush=True)
            if os.path.isfile(fn_full) and not force:
                print(f'{fn} already present.', flush=True)
            else:
                print(f'Downloading {fn}...', flush=True)
                urllib.request.urlretrieve(urls[i].format(fn), fn_full)

        file_base = [('ALL.chr{}.shapeit2_integrated_snvindels_v2a_27022019.'
            'GRCh38.phased.vcf.gz')]
        url_base = ('ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/'
            'data_collections/1000_genomes_project/release/'
            '20190312_biallelic_SNV_and_INDEL/{}')
        ext = re.escape('.vcf.gz')
    elif db == 'gnomad':
        file_base = ['gnomad.genomes.r3.0.sites.chr{}.vcf.bgz']
        url_base = ('https://storage.googleapis.com/gnomad-public/release/'
            '3.0/vcf/genomes/{}')
        ext = re.escape('.vcf.bgz')
    
    # Attempt to download using FTP, otherwise use HTTPS
    try:
        if url_base[:3] == 'ftp':
            serv, *directory = url_base[6:].split('/')
        else:
            raise ftplib.socket.gaierror
        ftp = ftplib.FTP(serv)
        ftp.login()
        ftp.cwd('/'.join(directory[:-1]))
        all_fns = ftp.nlst()
        use_ftp = True
    except ftplib.socket.gaierror:
        use_ftp = False

    for c in chrs:
        if use_ftp:
            print('Using FTP', flush=True)
            fn_match = f'\.(?:chr)?{c}\..*{ext}$'
            print(f'Looking for remote files matching {fn_match}...',
                flush=True)
            chr_fns = [fn for fn in all_fns if re.search(fn_match, fn)]
            if len(chr_fns) == 0:
                print(f'No file found for chr {c}.', flush=True)
                continue
            for fn in chr_fns:
                fn_full = f'{fp}/{fn}'
                
                print(f'Checking for {fn}...', flush=True)
                if os.path.isfile(fn_full) and not force:
                    print(f'{fn} already present.', flush=True)
                    continue

                with open(fn_full, 'wb') as ftp_fp:
                    print(f'Downloading {fn}...', flush=True)
                    ftp.retrbinary(f'RETR {fn}', ftp_fp.write)
        else:
            print('Using HTTPS', flush=True)
            for i in range(len(file_base)):
                fn = file_base[i].format(c)
                fn_full = f'{fp}/{fn}'

                print(f'Checking for {fn}...', flush=True)
                if os.path.isfile(fn_full) and not force:
                    print(f'{fn} already present.', flush=True)
                    continue

                try:
                    print(f'Downloading {fn}...', flush=True)
                    urllib.request.urlretrieve(url_base.format(fn), fn_full)
                except (urllib.error.HTTPError, ftplib.error_perm,
                    urllib.error.URLError) as e:
                    print(f'Error downloading {fn}.', flush=True)
                    continue
    if use_ftp:
        ftp.quit()

def col_from_db(db, chrs, fp='.', chr_path=None, chr_maj=False, chr_all=False,
    quiet=False):
    """
    Create a RSCollection object for all variants on the given chromosome in the
    given database.
    Valid databases are:
    * dbSNP
    * 1000gp
    * gnomAD
    Database names are not case sensitive.

    Parameters:
    db: Name of database to use
    chrs: List of chromosomes to use
    fp: File path where the database files are stored
    chr_path: File path to save individual chromosome ConsDB files
    chr_maj: Indicates whether to save major allele variants
    chr_all: Indicates whether to save all variants
    quiet: Do not print log/progress messages
    """

    db = db.lower()
    if db == 'dbsnp':
        file_base = ['{}/refsnp-chr{}.json.bz2']
        db_parse = rse.RSCollection.from_dbsnp
    elif db == '1000gp':
        file_base = [('{}/ALL.chr{}.shapeit2_integrated_snvindels_v2a_27022019.'
            'GRCh38.phased.vcf.gz')]
        db_parse = lambda fn, quiet: rse.RSCollection.from_1000gp(fn,
            f'{fp}/1000genomes.sequence.index',
            f'{fp}/20131219.populations.tsv', quiet)
    elif db == 'gnomad':
        file_base = ['{}/gnomad.genomes.r2.1.1.sites.{}.vcf.bgz',
            '{}/gnomad.exomes.r2.1.1.sites.{}.vcf.bgz']
        db_parse = rse.RSCollection.from_gnomad

    rsc = rse.RSCollection()
    for c in chrs:
        for fn in file_base:
            fn = fn.format(fp, c)
            if not os.path.isfile(fn): continue
            print(f'Loading {fn}...', flush=True)
            s = time.time()
            rsc += db_parse(fn, quiet=quiet)
            e = time.time()
            print(f'Finished {fn}. Took {e-s} s.', flush=True)
        if chr_path:
            if chr_all:
                fn = f'{chr_path}/chr{c}_rscol.gz'
                print(f'Saving chromosome {c}...', flush=True)
                rsc.dump(fn, c)
                print(f'Chromosome {c} saved.', flush=True)
            if chr_maj:
                fn = f'{chr_path}/chr{c}_maj.gz'
                print(f'Saving chromosome {c} major alleles...', flush=True)
                rsc.get_major().dump(fn, c)
                print(f'Chromosome {c} major alleles saved.', flush=True)

    return(rsc)

def consdb_to_fasta(consdb_fn, in_fa, out_fa):
    if check_consdb_sorted(consdb_fn):
        lines = rse.RSCollection.open(consdb_fn)
    else:
        print((f'WARNING: File {consdb_fn} is not sorted, this operation may '
            'take a lot of time/RAM. Consider pre-sorting.'))

        lines = sorted([line.decode() if type(line) == bytes else line \
            for line in rse.RSCollection.open(consdb_fn)],
            key=rse.RSCollection.sort_rsidx_line)

    ## Load original genome
    orig_gen, chr_header = load_genome(in_fa)

    ## Go through ConsDB file and create new genome
    new_gen = {}
    gen_ctr = {c: 0 for c in orig_gen.keys()}
    for line in lines:
        try:
            line = line.decode()
        except AttributeError:
            pass

        if line[0] == '#':
            continue

        line = re.split('[:,]', line.strip())
        chrom = line[1]
        pos = int(line[2]) - 1
        rec_ref = line[3]
        rec_alt = line[4]

        # Make sure that we have a reference genome for this chromosome
        try:
            g = orig_gen[chrom]
        except KeyError:
            print(f'No reference found for chr{chrom}, skipping variants.')
            orig_gen[chrom] = None
            continue
        if g is None:
            continue

        if pos >= len(g):
            continue

        if gen_ctr[chrom] > pos:
            print((f'Skipping variant '
                f'{line[1]}\t{line[2]}\t{line[3]}\t{line[4]}, overlapped by '
                'prior variant.'))
            continue

        try:
            n = new_gen[chrom]
        except KeyError:
            n = new_gen[chrom] = []

        if g[pos:pos+len(rec_ref)] != rec_ref:
            raise AssertionError(('Reference genome and VCF file disagree at '
                f'{line[1]} {line[2]} (ref: {g[pos:len(rec_ref)]}, '
                f'vcf: {rec_ref}).'))

        n.append(g[gen_ctr[chrom]:pos])
        n.append(rec_alt)
        gen_ctr[chrom] = pos + len(rec_ref)

    ## Print new genome
    with open(out_fa, 'w') as fp:
        for c in CHRS:
            try:
                g = new_gen[c]
                if gen_ctr[c] < len(orig_gen[c]):
                    g.append(orig_gen[c][gen_ctr[c]:])
            except KeyError:
                if c in orig_gen and orig_gen[c]:
                    g = orig_gen[c]
                else:
                    continue


            fp.write(chr_header[c])
            g = ''.join(g)
            for i in range(0, len(g), 70):
                fp.write(f'{g[i:i+70]}\n')

    try:
        lines.close()
    except AttributeError:
        pass

def filter_vcf_pers(fn_in, fn_out, pers_id, het=None):
    """
    Filter a VCF file for a given individual.
    Possible choices for het are:
    * 'rand': Randomly pick a haplotype at each position
    * [0,1]: Choose either left or right haplotype at each position

    Parameters:
    fn_in: VCF file to filter
    fn_out: Output filename
    pers_id: ID of the individual to filter
    het: How to handle heterozygous variants
    """

    try:
        het = int(het)
    except ValueError:
        het = het.lower()

    if het not in {'rand', 0, 1}:
        raise ValueError(f'Bad option for het: {het}')

    
    samp_idx = None
    pos_set = set()

    with open(fn_out, 'w') as fp:
        for line in rse.RSCollection.open(fn_in):
            try:
                line = line.decode()
            except AttributeError:
                pass

            if line[:2] == '##':
                fp.write(line)
                continue

            line = np.asarray(line.strip().split('\t'))
            if line[0] == '#CHROM' and samp_idx is None:
                samp_idx = np.asarray(line) == pers_id
                samp_idx[:9] = True
                line = '\t'.join(line[samp_idx])
                fp.write(f'{line}\n')
                continue

            if line[0][0] == '#':
                continue

            rec_alts = line[4].split(',')
            rec_info = {}
            for inf in line[7].split(';'):
                inf = inf.split('=')
                k = inf[0]
                try:
                    v = inf[1].split(',')
                    if len(v) == 1:
                        v = v[0]
                except IndexError:
                    v = 0
                rec_info[k] = v

            line = line[samp_idx]

            # Only want one first variant at any given position, just take the
            #  first one seen
            if tuple(line[:2]) in pos_set:
                continue

            gt = re.split(r'[|/]', line[-1])
            if (len(np.unique(gt)) == 1):
                if gt[0] == '0':
                    continue

                pos_set.add(tuple(line[:2]))

                # Need to adjust alt index (0 is ref in VCF files)
                line[4] = rec_alts[int(gt[0])-1]
                line[7] = ';'.join([f'{k}={v[int(gt[0])-1]}' if type(v) is list \
                    else (f'{k}={v}' if v != 0 else k) \
                    for k,v in rec_info.items()])

                line = '\t'.join(line)
                fp.write(f'{line}\n')
            elif het == 'rand':
                gt = np.random.choice([0,1])
                if gt == 0:
                    continue

                pos_set.add(tuple(line[:2]))

                # Need to adjust alt index (0 is ref in VCF files)
                line[4] = rec_alts[gt-1]
                line[7] = ';'.join([f'{k}={v[gt-1]}' if type(v) is list \
                    else (f'{k}={v}' if v != 0 else k) \
                    for k,v in rec_info.items()])

                line = '\t'.join(line)
                fp.write(f'{line}\n')
            elif type(het) == int:
                try:
                    gt = int(gt[het])
                except IndexError:
                    raise ValueError(f'Bad argument for het: {het}.')
                if gt == 0:
                    continue

                pos_set.add(tuple(line[:2]))

                # Need to adjust alt index (0 is ref in VCF files)
                line[4] = rec_alts[gt-1]
                line[7] = ';'.join([f'{k}={v[gt-1]}' if type(v) is list \
                    else (f'{k}={v}' if v != 0 else k) \
                    for k,v in rec_info.items()])

                line = '\t'.join(line)
                fp.write(f'{line}\n')

def find_rsid_frags(c, fp_in, fns, vcf=True):
    """
    Build dictionary mapping from (chr, pos) -> list of files containing that
    position. This is used in the load_and_save method so that chunks can be
    saved safely without worrying about having two entries at the same position
    in the final ConsDB file.

    Parameters:
    c: Chromosome
    fp_in: Input file path
    fns: List of fil fragment names to look for variants in
    vcf: Indicates if the passed file fragments came from a VCF file
    """

    frag_dict = {}
    n_frags = len(fns)

    if vcf:
        rs_match = ('^(?:chr)?([0-9XYMT]{1,2})\s+([0-9]+)\s+([.;a-z0-9]+)\s+'
            '([A-Z]+)\s+([A-Z,]+)')
    else:
        rs_match = '"refsnp_id":"([0-9]+)"'
    count = 1
    for fn in fns:
        if count % 15 == 0:
            print((f'chr {c}: Finding rsid fragments '
                f'(fragment {count}/{n_frags})'), flush=True)
        count += 1

        with open(f'{fp_in}/{fn}', 'r') as fp:
            for line in fp:
                m = re.search(rs_match, line)
                try:
                    if vcf:
                        # Use (chr, pos) to avoid entries that may need to be
                        #  grouped together but one is missing rsid
                        # Need to account for cases with multiple alts
                        rsid = []
                        for alt in m.group(5).split(','):
                            rsid.append((m.group(1), m.group(2)))
                    else:
                        rsid = [m.group(1)]
                except AttributeError as e:
                    if line[0] != '#' and '<' not in line.split('\t')[4]:
                        print(f'chr {c}: rsid not found in {line}', flush=True)
                        raise e
                    continue

                for rs in rsid:
                    try:
                        frag_dict[rs].add(fn)
                    except KeyError:
                        frag_dict[rs] = {fn}

    return(frag_dict)

def fragment_files(c, db, fp_in='.', fp_out='./frags/', nlines=10000):
    """
    Function to create file fragments of nlines for the given chromosome.
    Designed to approximate the split GNU program.
    Valid databases are:
    * dbSNP
    * 1000gp
    * gnomAD
    Database names are not case sensitive.

    Parameters:
    c: Chromosome to use
    db: Database being used
    fp_in: Input file path
    fp_out: Output file path to store file fragments
    nlines: Number of lines each fragment should contain
    """

    c_fp = f'{fp_out}/chr{c}_frag_'

    db = db.lower()
    if db == 'dbsnp':
        ext = re.escape('.json.bz2')
    elif db == 'gnomad':
        ext = re.escape('vcf.bgz')
    elif db == '1000gp':
        ext = re.escape('.vcf.gz')
    else:
        print(f'db: {db}')
        raise AssertionError

    fns = os.listdir(fp_in)
    fn_match = f'\.(?:chr)?{c}\..*{ext}$'
    chr_fns = [fn for fn in fns if re.search(fn_match, fn)]
    frag_fns = []

    for fn in chr_fns:
        fn = f'{fp_in}/{fn}'

        fn_lines = 0
        for line in rse.RSCollection.open(fn):
            try:
                line = line.decode()
            except AttributeError:
                pass

            if line[0] == '#':
                continue
            else:
                fn_lines += 1

        # Calculate number of fragment files needed and max number of letters
        #  needed
        n_files = np.ceil(fn_lines / nlines)
        n_frag_letter = int(np.ceil(np.log(n_files)/np.log(26)))

        header = []
        frag = np.zeros(n_frag_letter, dtype=int)
        counter = 0
        fp_frag = open(f'{c_fp}{idx_to_frag(frag)}', 'w')
        frag_fns.append(f'{c_fp}{idx_to_frag(frag)}')
        for line in rse.RSCollection.open(fn):
            try:
                line = line.decode()
            except AttributeError:
                pass

            if line[0] == '#':
                header.append(line)
                continue

            # If we've written nlines to the current file, open up the next file
            if counter == nlines:
                counter = 0
                frag[-1] += 1
                i = -1
                # Increment any indices that need to be incremented
                while frag[i] == 26:
                    frag[i] = 0
                    frag[i-1] += 1
                    i -= 1

                fp_frag.close()
                fp_frag = open(f'{c_fp}{idx_to_frag(frag)}', 'w')
                frag_fns.append(f'{c_fp}{idx_to_frag(frag)}')

            # Write the VCF header at the beginning of any new file
            if counter == 0:
                fp_frag.write(''.join(header))

            fp_frag.write(line)

            counter += 1

    return(frag_fns)

def idx_to_frag(idx):
    """
    Convert a list of integers in the range [0,25] to an alphabetic string of
    the same length. Helper method for fragment_files().
    0 -> a
    1 -> b
    ...

    Parameters:
    idx: List of alphabetic indeces
    """

    return(''.join([ascii_lowercase[i] for i in idx]))

def join_fns(fns, fn_out):
    """
    Join multiple files together by simple concatenation. Intended to be used to
    combine multiple ConsDB files together after parallelised parsing.

    Parameters:
    fns: List of files
    fn_out: File to output to
    """

    if os.path.isfile(fn_out):
        os.remove(fn_out)

    with open(fn_out, 'wb') as fp_out:
        for fn in fns:
            for line in open(fn, 'rb'):
                fp_out.write(line)

def join_vcfs(fns, fn_out):
    """
    Join multiple VCF files together. Files will be added in order, so they
    should already be sorted by chromosome.

    Parameters:
    fns: List of files
    fn_out: File to output to
    """

    # Get meta information
    header_lines = []
    filedate = False
    for fn in fns:
        for line in open(fn, 'r'):
            line = line.strip()
            if line[:2] != '##':
                break

            if line not in header_lines:
                if 'fileDate' in line:
                    if filedate:
                        continue
                    filedate = True
                header_lines.append(line)

    # Get header line
    for line in open(fns[0], 'r'):
        line = line.strip()
        if line.split('\t')[0] == '#CHROM':
            header_lines.append(line)
            break
    
    with open(fn_out, 'w') as fp_out:
        fp_out.write('\n'.join(header_lines) + '\n')
        for fn in fns:
            for line in open(fn, 'r'):
                if line[0] == '#':
                    continue

                fp_out.write(line)

def load_genome(gen_fn):
    """
    Load a genome from a FASTA file. Returns a dictionary of chrom -> DNA seq
    and a dictionary of chromosome headers.

    Parameters:
    gen_fn: Genome file
    """

    orig_gen = {}
    chr_header = {}
    chrom = ''
    for line in open(gen_fn, 'r'):
        if line[0] == '>':
            if chrom != '':
                orig_gen[chrom] = ''.join(orig_gen[chrom])

            chrom = line.split()[0].strip('>chr')
            chr_header[chrom] = line
            orig_gen[chrom] = []
            continue

        orig_gen[chrom].append(line.strip())
    orig_gen[chrom] = ''.join(orig_gen[chrom])

    return(orig_gen, chr_header)

def make_cons_vcf(fn_in, fn_out, pop=None, quiet=True):
    """
    Make a consensus VCF file using the BitRSCollection class.

    Parameters:
    fn_in: ConsDB file to load
    fn_out: VCF file to save to
    pop: Which population to use
    quiet: Load quietly
    """

    if pop:
        col = srsc.BitRSCollection.load_from_file_pop(fn_in, pop, quiet)
    else:
        col = srsc.BitRSCollection.load_from_file_full(fn_in, quiet)

    col.dump_vcf(fn_out, cons=True, is_maj=('maj' in fn_in))

def merge_files(c, fps, fp_out, fmt='chr{}_{}.gz', merge_all=False,
    merge_maj=False):
    """
    Function to merge multiple ConsDB files. Used as a wrapper for the
    corresponding function in the RSEntry module.

    Parameters:
    c: Chromsome to use
    fps: File paths to look for input files
    fp_out: Output file path
    fmt: Format to use to look for files
    merge_all: Indicates whether to merge files with all variants
    merge_maj: Indicates whether to merge files with major allele variants
    """

    fns = []
    # Check how many placeholders in format to avoid loading/saving the same
    #  file multiple times
    if len(re.findall('\{\}', fmt)) == 1:
        fns += [('.', [f'{fp}/{fmt.format(c)}' for fp in fps \
            if os.path.isfile(f'{fp}/{fmt.format(c)}')])]
    else:
        if merge_all:
            fns += [('rscol', [f'{fp}/{fmt.format(c, "rscol")}' for fp in fps \
                if os.path.isfile(f'{fp}/{fmt.format(c, "rscol")}')])]
        if merge_maj:
            fns += [('maj', [f'{fp}/{fmt.format(c, "maj")}' for fp in fps \
                if os.path.isfile(f'{fp}/{fmt.format(c, "maj")}')])]

    if all([len(i[1]) == 0 for i in fns]):
        print(f'No files found for chr {c}.', flush=True)
        return

    for (ext, fn_list) in fns:
        # Know each file is only present once, and extra arguments to format are
        #  ignored, so don't need to check for multiple fmt placeholders
        fn_out = f'{fp_out}/{fmt.format(c, ext)}'
        rse.RSCollection.merge_files(fn_list, fn_out, c)

def mp_load_and_save(c, db, fp_in='.', fp_out=None, store_all=False,
    store_maj=False, quiet=False):
    """
    Function to load an entire database file for the given chromosome and save
    it as a ConsDB file.
    Valid databases are:
    * dbSNP
    * 1000gp
    * gnomAD
    Database names are not case sensitive.

    Parameters:
    c: Chromosome to use
    db: Database being used
    fp_in: Input file path of the database files
    fp_out: Output file path for ConsDB files
    store_all: Indicates whether to store ConsDB files with all variants
    store_maj: Indicates whether to store ConsDB files with major allele
        variants
    quiet: Disable log/progress messages
    """

    db = db.lower()
    if db == 'dbsnp':
        ext = re.escape('.json.bz2')
        db_parse = lambda fn: rse.RSCollection.from_dbsnp(fn, quiet)
    elif db == '1000gp':
        ext = re.escape('.vcf.gz')
        db_parse = lambda fn: rse.RSCollection.from_1000gp(
            fn, f'{fp_in}/1000genomes.sequence.index',
            f'{fp_in}/20131219.populations.tsv', quiet)
    elif db == 'gnomad':
        ext = re.escape('.vcf.bgz')
        db_parse = lambda fn: rse.RSCollection.from_gnomad(fn, quiet)

    fns = os.listdir(fp_in)
    fn_match = f'\.(?:chr)?{c}\..*{ext}$'
    chr_fns = [fn for fn in fns if re.search(fn_match, fn)]

    rsc = rse.RSCollection()
    for fn in chr_fns:
        # Copying the code from col_from_db here to avoid having to return
        #  large variables
        fn = f'{fp_in}/{fn}'
        if not os.path.isfile(fn):
            print(f'{fn} not found.', flush=True)
            continue
        print(f'Loading {fn}...', flush=True)
        s = time.time()
        rsc += db_parse(fn)
        e = time.time()
        print(f'Finished {fn}. Took {e-s} s.', flush=True)

    if fp_out:
        if store_all:
            print(f'Saving chromosome {c}...', flush=True)
            rsc.dump(f'{fp_out}/chr{c}_rscol.gz', c)
            print(f'Chromosome {c} saved.', flush=True)
        if store_maj:
            print(f'Saving chromosome {c} major alleles...', flush=True)
            rsc.get_major().dump(f'{fp_out}/chr{c}_maj.gz', c)
            print(f'Chromosome {c} major alleles saved.', flush=True)



def mp_load_and_save_chr(c, db, fp_db, fp_in, fp_out, store_all=False,
    store_maj=False, quiet=False):
    """
    Function to load database fragment files for the given chromosome and save
    as a ConsDB file.
    Valid databases are:
    * dbSNP
    * 1000gp
    * gnomAD
    Database names are not case sensitive.

    Parameters:
    c: Chromosome to use
    db: Database being used
    fp_db: File path containing the original database files
    fp_in: File path containing the database file fragments
    fp_out: Output file path for ConsDB files
    store_all: Indicates whether to store ConsDB files with all variants
    store_maj: Indicates whether to store ConsDB files with major allele
        variants
    quiet: Whether to disable log messages
    """

    db = db.lower()
    if db == 'dbsnp':
        db_parse = lambda fn: rse.RSCollection.from_dbsnp(fn, quiet)
        vcf = False
    elif db == '1000gp':
        db_parse = lambda fn, quiet: rse.RSCollection.from_1000gp(fn,
            f'{fp_db}/1000genomes.sequence.index',
            f'{fp_db}/20131219.populations.tsv', quiet)
        vcf = True
    elif db == 'gnomad':
        db_parse = lambda fn: rse.RSCollection.from_gnomad(fn, quiet)
        vcf = True
    else:
        print(f'db: {db}', flush=True)
        raise AssertionError

    if db == 'gnomad':
        match = [f'chr{c}_frag_exomes_[a-z]+', f'chr{c}_frag_genomes_[a-z]+']
    else:
        match = [f'chr{c}_frag_[a-z]+']
    fns = [fn for fn in os.listdir(fp_in) for m in match if re.search(m, fn)]
    n_frags = len(fns)
    # Create dictionary of the last fragment a (chrom, pos) was found in to know
    #  when it's ok to save an entry
    rsid_frag_dict = find_rsid_frags(c, fp_in, fns, vcf)
    if not quiet:
        print(f'chr {c}:', len(rsid_frag_dict), 'rsids total', flush=True)

    # Remove target files if they exist so we can open in append mode
    save_fn = f'{fp_out}/chr{c}_rscol.gz'
    idx_file = f'{save_fn}.idx'
    maj_fn = f'{fp_out}/chr{c}_maj.gz'
    maj_idx = f'{maj_fn}.idx'
    if os.path.exists(save_fn):
        os.remove(save_fn)
    if os.path.exists(idx_file):
        os.remove(idx_file)
    if os.path.exists(maj_fn):
        os.remove(maj_fn)
    if os.path.exists(maj_idx):
        os.remove(maj_idx)

    if not quiet:
        print (f'Loading chromosome {c}...', flush=True)
    s = time.time()
    all_size = 0
    maj_size = 0
    count = 1
    rsc = rse.RSCollection()
    seen_frags = set()
    write_rsids = set()
    for fn in fns:
        if count % 15 == 0 and not quiet:
            print(f'chr {c}: Loading fragment {count}/{n_frags}', flush=True)

        rsc += db_parse(f'{fp_in}/{fn}', quiet=quiet)
        count += 1
        
        frag = fn.split('_')[-1]
        if db == 'dbsnp':
            write_rsids.update([k for k in rsc.entries.keys() \
                if rsid_frag_dict[k[0]].issubset(seen_frags)])
        else:
            write_rsids.update([k for k in rsc.entries.keys() \
                if rsid_frag_dict[(str(k[1]), str(k[2]))].issubset(seen_frags)])

        if not quiet:
            print(fn)
            print(len(write_rsids), 'rsids in write_rsids', flush=True)
            print(len(rsc), 'rsids in rsc', flush=True)
        if (len(write_rsids) >= 50000 or fn == fns[-1]) and fp_out:

            if store_all:
                if not quiet:
                    print(f'Writing {len(write_rsids)} entries to chromosome '
                        f'{c}...', flush=True)
                all_size = rsc.dump(save_fn, c, rsids=write_rsids,
                    old_size=all_size, append=True)
                if not quiet:
                    print(f'Wrote to chromosome {c}.', flush=True)
            if store_maj:
                maj_size = rsc.get_major().dump(maj_fn, c, rsids=write_rsids,
                    old_size=maj_size, append=True)

            rsc = rse.RSCollection()
            write_rsids = set()

    e = time.time()
    if not quiet:
        print(f'Finished chromosome {c}. Took {e-s} s.', flush=True)

def parse_chr_file(fn):
    """
    Function to extract chromosomes to use from the given file.

    Parameters:
    fn: File name to use
    """

    chrs = np.loadtxt(fn, dtype=object)
    chr_idx = np.asarray([c in CHRS for c in chrs])
    for c in chrs[~chr_idx]:
        print('Unknown chromosome {}.'.format(c))
        
    chrs = chrs[chr_idx]
    print('Chromosomes to use: {}'.format(', '.join(chrs)))

    return(chrs)

def print_traceback(c, err, tb):
    """
    Function to print a traceback from given error and traceback objects.

    Parameters:
    c: Chromosome to use
    err: Exception object
    tb: Traceback object
    """

    tb_str = '##########################\n'

    tb_str += (f'Error "{c.__name__}: {err}" occured on line {tb.tb_lineno} '
        f'in {tb.tb_frame.f_code.co_filename}.\n')
    tb_str += 'Traceback:\n'

    while tb.tb_next:
        tb_str += (f'line {tb.tb_lineno} in {tb.tb_frame.f_code.co_filename} '
            '->\n')
        tb = tb.tb_next
    tb_str += f'line {tb.tb_lineno} in {tb.tb_frame.f_code.co_filename}\n'

    tb_str += '##########################'
    print(tb_str, flush=True)

def queue_wrap(q, func, *args, **kwargs):
    """
    Wrapper function around whatever function is being run in parallel. This
    allows for catching errors thrown from a process in a multiprocessing Pool.
    Normally these errors aren't thrown until all processes finish, which can
    result in a lot of wasted time.

    q: Queue object used for sending error messages to the parent process
    func: Function to be called
    args: Positional arguments to be sent to func
    kwargs: Keyword arguments to be sent to func
    """

    try:
        res = func(*args, **kwargs)
    except Exception as e:
        exc = sys.exc_info()
        print(f'Sending message to queue.', flush=True)
        print_traceback(*exc)
        q.put((None,e))
    else:
        q.put(('END',None))

    return(res)

def vcf_to_fasta(vcf_fn, in_fa, out_fa):
    ## Load original genome
    orig_gen, chr_header = load_genome(in_fa)

    ## Go through VCF file and create new genome
    new_gen = {}
    gen_ctr = {c: 0 for c in orig_gen.keys()}
    for line in open(vcf_fn, 'r'):
        if line[0] == '#':
            continue

        line = line.strip().split('\t')
        chrom = line[0].strip('chr')
        pos = int(line[1]) - 1
        rec_ref = line[3]
        # If multiple alts, just take the first one
        rec_alt = line[4].split(',')[0]

        # Make sure that we have a reference genome for this chromosome
        try:
            g = orig_gen[chrom]
        except KeyError:
            print(f'No reference found for chr{chrom}, skipping variants.')
            orig_gen[chrom] = None
            continue
        if g is None:
            continue

        if pos >= len(g):
            continue

        if gen_ctr[chrom] > pos:
            print((f'Skipping variant '
                f'{line[0]}\t{line[1]}\t{line[3]}\t{line[4]}, overlapped by '
                'prior variant.'))
            continue

        try:
            n = new_gen[chrom]
        except KeyError:
            n = new_gen[chrom] = []

        if g[pos:pos+len(rec_ref)] != rec_ref:
            raise AssertionError(('Reference genome and VCF file disagree at '
                f'{line[0]} {line[1]} (ref: {g[pos:len(rec_ref)]}, '
                f'vcf: {rec_ref}).'))

        n.append(g[gen_ctr[chrom]:pos])
        n.append(rec_alt)
        gen_ctr[chrom] = pos + len(rec_ref)

    ## Print new genome
    with open(out_fa, 'w') as fp:
        for c in CHRS:
            try:
                g = new_gen[c]
                if gen_ctr[c] < len(orig_gen[c]):
                    g.append(orig_gen[c][gen_ctr[c]:])
            except KeyError:
                if c in orig_gen and orig_gen[c]:
                    g = orig_gen[c]
                else:
                    continue


            fp.write(chr_header[c])
            g = ''.join(g)
            for i in range(0, len(g), 70):
                fp.write(f'{g[i:i+70]}\n')



################################################################################
def get_args():
    """
    Parse command line arguments.

    Valid run modes are:
    Parse
    Filter
    Merge
    Cons
    """

    run_mode = sys.argv[1].lower()

    if run_mode == 'parse':
        parser = argparse.ArgumentParser(prog='ConsDB Parse',
            description='Parse a variant database and output ConsDB files.')

        ## General arguments
        parser.add_argument('-o', help='File to store output.')
        parser.add_argument('-maj', help='File to store major alleles.')
        parser.add_argument('-chr', default=['all'], nargs='+',
            help=('Which chromosome(s) to use. Can either be a list of '
                'numbers/letters or a file containing each chromosome to use '
                'on a new line. Use "all" for all chromosomes.'))
        parser.add_argument('-quiet', action='store_true',
            help='Parse quietly.')

        ## Database arguments
        parser.add_argument('-db', default='1000GP', help= ('Database to use '
            '[dbSNP 1000GP gnomAD]. Uses 1000GP if no argument is supplied.'))
        parser.add_argument('-db_path', default='.', help='Where to store '
            'database downloads. Defaults to run path.')
        parser.add_argument('-db_force', action='store_true', help='Whether or '
            'not to force database downloading (overwrite existing files).')
        parser.add_argument('-db_mp', action='store_true', 
            help=('Whether or not to use multiprocessing for downloading.'))
        parser.add_argument('-db_proc', type=int, default=12,
            help='Number of downloads to run concurrently.')

        ## Multiprocessing arguments
        parser.add_argument('-mp', action='store_true',
            help=('Whether or not to process files using multiple cores.'))
        parser.add_argument('-mp_path', default='.',
            help=('Where to store output files. Defaults to run path.'))
        parser.add_argument('-mp_o', action='store_true',
            help='Whether or not to store files containing all allele.')
        parser.add_argument('-mp_maj', action='store_true',
            help='Whether or not to store files containing only major allele.')
        parser.add_argument('-mp_proc', type=int, default=12,
            help='Number of processes to run concurrently.')
        parser.add_argument('-preprocess', action='store_true',
            help='Whether or not to fragment database files before processing.')
        parser.add_argument('-pp_path', default='./frags/',
            help='Directory to store fragment files. Defaults to run path.')
        parser.add_argument('-pp_nlines', type=int, default=100000,
            help='Number of lines per fragment file.')
        parser.add_argument('-pp_clean', action='store_true',
            help='Whether or not to remove fragment files.')
        parser.add_argument('-pp_proc', type=int, default=12,
            help='Number of processes to run concurrently.')
        parser.add_argument('-no_frag', action='store_true',
            help='Skip fragmenting database files.')

        args = parser.parse_args(sys.argv[2:])

        if not args.o and not args.maj and not \
            (args.mp and (args.mp_o or args.mp_maj)):
            raise RuntimeError(
                'Must supply an output argument (-o, -maj, -mp_o, -mp_maj).')
        if args.mp and not args.mp_o and not args.mp_maj:
            raise RuntimeError(('Must select a multiprocessing output '
                '(-mp_o and/or -mp_maj).'))
        if args.o and args.mp and not args.mp_o:
            raise RuntimeError(('Must use -mp_o when specifying -o in '
                'multiprocessing mode.'))
        if args.maj and args.mp and not args.mp_maj:
            raise RuntimeError(('Must use -mp_maj when specifying -maj in '
                'multiprocessing mode.'))

        return(args)

    ## Filter VCF file
    if run_mode == 'filter':
        parser = argparse.ArgumentParser(prog='ConsDB Filter',
            description='Filter a VCF file based on ConsDB files.')

        parser.add_argument('-i', required=True,
            help='Input VCF for filtering.')
        parser.add_argument('-o', required=True, help='VCF filtering output.')
        parser.add_argument('-consdb_path',
            help='Directory containing ConsDB files.')
        parser.add_argument('-pop', help='Population to use for filtering.')
        parser.add_argument('-samp', help='Sample ID to use for filtering.')
        parser.add_argument('-het',
            help='Option for handling heterozygous variants.')
        parser.add_argument('-log', help='Log file to use for filtering.')
        parser.add_argument('-cons', action='store_true', help=('Making a '
            'consenus so only keep first major allele at each position.'))
        parser.add_argument('-keep_samps', action='store_true',
            help='Keep sample information in the filtered VCF.')
        parser.add_argument('-quiet', action='store_true',
            help='Filter quietly.')

        args = parser.parse_args(sys.argv[2:])

        if args.pop and args.samp:
            raise RuntimeError('Can only specify one of (-pop, -samp).')
        if (args.pop or args.samp is None) and args.consdb_path is None:
            raise RuntimeError(('Must specify -consdb_path when not filtering '
                'by sample.'))

        return(args)

    ## Merge files
    if run_mode == 'merge':
        parser = argparse.ArgumentParser(prog='ConsDB Merge',
            description='Merge multiple ConsDB files.')

        parser.add_argument('-i', required=True, nargs='+',
            help='Input files to be merged.')
        parser.add_argument('-o', required=True, help='Merged output file.')
        parser.add_argument('-chr', default=['all'], nargs='+',
            help=('Which chromosome(s) to use. Can either be a list of '
                'numbers/letters or a file containing each chromosome to use '
                'on a new line. Use "all" for all chromosomes.'))
        parser.add_argument('-fmt', default='chr{}_{}.gz',
            help=('Format used to find files to merge. Use {} as a wildcard, '
                'with one for chromosome, and the second for rscol/maj '
                '(if desired).'))
        
        ## Multiprocessing arguments       
        parser.add_argument('-mp', action='store_true',
            help='Use multiprocessing.')
        parser.add_argument('-mp_proc', type=int, default=12,
            help='Number of processes to run concurrently.')
        parser.add_argument('-merge_all', action='store_true',
            help='Merge ConsDB files with all alleles.')
        parser.add_argument('-merge_maj', action='store_true',
            help='Merge ConsDB files with major alleles.')

        args = parser.parse_args(sys.argv[2:])

        if all([os.path.isfile(fn) for fn in args.i]):
            args.inp_type = 'file'
        elif all([os.path.isdir(fn) for fn in args.i]):
            args.inp_type = 'dir'
        else:
            raise RuntimeError('Input arguments must be all files or all '
                'directories.')

        if args.inp_type == 'file' and \
            (os.path.exists(args.o) != os.path.isfile(args.o)):
            raise RuntimeError(('Input arguments and output argument must be '
                'the same type (all files or all directories).'))
        if args.inp_type == 'dir':
            if os.path.exists(args.o) != os.path.isdir(args.o):
                raise RuntimeError(('Input arguments and output argument must '
                    'be the same type (all files or all directories).'))
            if not args.merge_all and not args.merge_maj:
                raise RuntimeError(('Must select -merge_all and/or -merge_maj '
                    'when merging from directories.'))

        if args.inp_type == 'file' and args.mp:
            raise RuntimeError(('Multiprocessing not currently supported when '
                'using a list of files as the input.'), flush=True)

        return(args)

    ## Make consensus VCF file
    if run_mode == 'cons':
        parser = argparse.ArgumentParser(prog='ConsDB Cons',
            description='Create consensus VCF file.')

        parser.add_argument('-i', required=True, nargs='+',
            help='Input ConsDB files or directory path to use.')
        parser.add_argument('-o', required=True, help='Output directory.')
        parser.add_argument('-fmt', default='chr{}_{}.gz',
            help=('Format used to find files. Use {} as a wildcard, with one '
                'for chromosome, and the second for rscol/maj (if desired).'))
        parser.add_argument('-pop', help='Population to make consensus for.')
        parser.add_argument('-chr', default=['all'], nargs='+',
            help=('Which chromosome(s) to use. Can either be a list of '
                'numbers/letters or a file containing each chromosome to use '
                'on a new line. Use "all" for all chromosomes.'))
        parser.add_argument('-join', help=('If present, file to concatenate '
            'individual chromosome VCF files to.'))
        parser.add_argument('-clean', action='store_true',
            help='Delete individual chromosome VCF files when done.')
        parser.add_argument('-v', action='store_true', help='Load verbosely.')

        ## Multiprocessing arguments       
        parser.add_argument('-mp', action='store_true',
            help='Use multiprocessing.')
        parser.add_argument('-mp_proc', type=int, default=12,
            help='Number of processes to run concurrently.')

        args = parser.parse_args(sys.argv[2:])

        if all([os.path.isfile(fn) for fn in args.i]):
            args.inp_type = 'file'
        elif all([os.path.isdir(fn) for fn in args.i]):
            args.inp_type = 'dir'
        else:
            raise RuntimeError('Input arguments must be all files or all '
                'directories.')

        if args.inp_type == 'dir' and len(args.i) > 1:
                raise RuntimeError('Only one input directory is supported.')

        return(args)

    ## Make a consensus FASTA file, either from a VCF file or from a ConsDB file
    if run_mode == 'fa':
        parser = argparse.ArgumentParser(prog='ConsDB FA',
            description='Create consensus FASTA file.')

        parser.add_argument('-ref', required=True,
            help='Reference genome FASTA file.')
        parser.add_argument('-t', required=True, help='File to transform with.')
        parser.add_argument('-o', required=True, help='Output file.')

        args = parser.parse_args(sys.argv[2:])

        # Check for VCF file suffix, otherwise assume it's a ConsDB file
        if re.search('\.vcf(?:\.gz)?$', args.t):
            args.ft = 'vcf'
        else:
            args.ft = 'consdb'

        return(args)

def main():
    try:
        run_mode = sys.argv[1].lower()
    except IndexError:
        raise RuntimeError('No run mode given.')

    if run_mode not in {'parse', 'filter', 'merge', 'cons', 'fa'}:
        raise RuntimeError(f'Unknown run mode {sys.argv[1]}.')
    
    args = get_args()

    if run_mode not in {'filter', 'fa'}:
        # Convert passed chromosomes argument to a list of chroms to use
        if args.chr == ['all']:
            chrs = CHRS
        elif os.path.isfile(args.chr[0]):
            chrs = parse_chr_file(args.chr[0])
        else:
            chrs = args.chr

    if run_mode == 'parse':
        # Open a manhole for process monitoring
        #  (https://github.com/ionelmc/python-manhole)
        try:
            manhole.install()
        except NameError:
            pass

        # Ensure all needed database files are downloaded (using multiprocessing
        #  if the appropriate command-line argument was passed)
        if args.db_mp:
            nproc = min(mp.cpu_count(), args.db_proc, len(chrs))
            with mp.Pool(processes=nproc) as pool:
                pool.starmap(check_downloads, [(args.db, [c], args.db_path,
                    args.db_force) for c in chrs])
        else:
            check_downloads(args.db, chrs, args.db_path, args.db_force)

        ## If multiprocessing is enabled
        if args.mp:
            rsc_queue = None

            m = mp.Manager()
            msg_q = m.Queue()
            nproc = min(mp.cpu_count(), args.mp_proc, len(chrs))

            # If using preprocessing (file fragmentation), fragment files as
            #  needed
            pp_proc = min(mp.cpu_count(), args.pp_proc, len(chrs))
            if args.preprocess and not args.no_frag:
                frag_fns = {}
                with mp.Pool(processes=pp_proc) as pool:
                    if not os.path.exists(args.pp_path):
                        os.mkdir(args.pp_path)

                    fn_lists = pool.starmap(fragment_files, [(c, args.db,
                        args.db_path, args.pp_path, args.pp_nlines) for c in chrs])
                    # Keep track of fragment files in order to delete them later
                    #  if desired
                    frag_fns = set(np.concatenate(fn_lists))

            with mp.Pool(processes=nproc) as pool:
                if args.preprocess:
                    mp_args = [(msg_q, mp_load_and_save_chr, c, args.db,
                        args.db_path, args.pp_path, args.mp_path,
                        args.mp_o, args.mp_maj, args.quiet) for c in chrs]
                    # Use async so we can check for error messages
                    rsc_list = pool.starmap_async(queue_wrap, mp_args)

                    # Check for messages from all processes until all chrs have
                    #  finished. If we get an error message, raise the error.
                    num_finished = 0
                    while num_finished < len(chrs):
                        print('Waiting for message from queue.', flush=True)
                        # Main process waits here until something is received
                        #  from the queue (either an error message or a finished
                        #  message)
                        m, err = msg_q.get()
                        if issubclass(err.__class__, Exception):
                            print('Error message received from queue.',
                                flush=True)
                            raise err
                        if m == 'END':
                            num_finished += 1
                            print((f'Finished {num_finished}/{len(chrs)} '
                                'chromosomes.'), flush=True)
                            continue
                        print('Queue received message:', m, err, flush=True)

                    # Remove fragment directory if requested
                    if args.pp_clean:
                        for fn in frag_fns:
                            os.remove(fn)

                else:
                    rsc_list = pool.starmap(mp_load_and_save,
                        [(c, args.db, args.db_path, args.mp_path,
                            args.mp_o, args.mp_maj, args.quiet) for c in chrs])

            # Join the separate chromosome ConsDB files into one large
            #  file for all alleles and one for major alleles (as requested)
            if args.o:
                fns = [f'{args.mp_path}/chr{c}_rscol.gz' for c in chrs]
                join_fns(fns, args.o)
            if args.maj:
                fns = [f'{args.mp_path}/chr{c}_maj.gz' for c in chrs]
                join_fns(fns, args.maj)

        else:
            rsc = col_from_db(args.db, chrs, args.db_path, args.mp_path,
                args.mp_maj, args.mp_o, args.quiet)

            if args.o:
                rsc.dump_full(args.o)
            if args.maj:
                rsc.get_major().dump_full(args.maj)

    ## Merge files
    if run_mode == 'merge':
        if args.inp_type == 'file':
            rse.RSCollection.merge_files(args.i, args.o)
        else:
            args.i = args.i[0]
            if args.mp:
                nproc = min(mp.cpu_count(), args.mp_proc, len(chrs))
                with mp.Pool(nproc) as pool:
                    fun_args = [(c, args.i, args.o, args.fmt,
                        args.merge_all, args.merge_maj) for c in chrs]
                    pool.starmap(merge_files, fun_args)
            else:
                for c in chrs:
                    merge_files(c, args.i, args.o, args.fmt, args.merge_all,
                        args.merge_maj)

    ## Filter VCF file
    if run_mode == 'filter':
        if args.samp:
            filter_vcf_pers(args.i, args.o, args.samp, args.het)
        else:
            srsc.BitRSCollection.filter_vcf(args.consdb_path, args.i, args.o,
                args.pop, args.log, args.cons, args.keep_samps, args.quiet)

    ## Make consensus VCF file
    if run_mode == 'cons':
        if args.inp_type == 'file':
            fns = args.i
        else:
            args.i = args.i[0]
            os.listdir(args.i)

            if args.pop is None:
                fns = []
                out_fns = []
                for c in chrs:
                    if os.path.isfile(f'{args.i}/{args.fmt.format(c, "maj")}'):
                        fns.append(f'{args.i}/{args.fmt.format(c, "maj")}')
                    elif os.path.isfile(f'{args.i}/{args.fmt.format(c, "rscol")}'):
                        fns.append(f'{args.i}/{args.fmt.format(c, "rscol")}')
                    else:
                        print(f'WARNING: No file found for chr {c}')
                        continue

                    out_fns.append((f'{args.o}/chr{c}_'
                        f'{args.pop if args.pop else "pan"}_cons.vcf'))
            else:
                fns = []
                out_fns = []
                for c in chrs:
                    if os.path.isfile(f'{args.i}/{args.fmt.format(c, "rscol")}'):
                        fns.append(f'{args.i}/{args.fmt.format(c, "rscol")}')
                    else:
                        print(f'WARNING: No file found for chr {c}')
                        continue

                    out_fns.append((f'{args.o}/chr{c}_'
                        f'{args.pop if args.pop else "pan"}_cons.vcf'))

        cmd_args = ((fns[i], out_fns[i], args.pop, not args.v) \
            for i in range(len(fns)))
        if args.mp:
            nproc = min(mp.cpu_count(), args.mp_proc, len(fns))
            with mp.Pool(nproc) as pool:
                pool.starmap(make_cons_vcf, cmd_args)
        else:
            for args_list in cmd_args:
                make_cons_vcf(*args_list)

        if args.join:
            join_vcfs(out_fns, args.join)
        if args.clean:
            for fn in out_fns:
                os.remove(fn)

    ## Make consensus FASTA file
    if run_mode == 'fa':
        if args.ft == 'vcf':
            vcf_to_fasta(args.t, args.ref, args.o)
        else:
            consdb_to_fasta(args.t, args.ref, args.o)

if __name__ == '__main__':
    main()
