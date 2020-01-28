import bz2
import gzip
import itertools
import numpy as np
import operator
import os
import re
import time

class BitRSCollection():
    BASE_ENC = {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3
    }

    """
        A sub-class to represent a variant using bit-packing. Provides a much
        more compact object which makes large-scale operations faster.

        Attributes
        --------------------
        entries : dict
            Dictionary that maps from (rsid, chrom, pos) to RSEntry object
        chr_pos_table : dict
            Dictionary that maps from (chrom, pos) to (rsid, chrom, pos)
        rsid_table : dict
            Dictionary that maps from rsid to (rsid, chrom, pos)

        Methods
        --------------------
        add_entry(rsidx, e)
            Add entry from a variant list
        add_entry_line(e)
            Add entry from a line in a ConsDB file
        add_entry_from_args(chrom, rsid, pos, quiet=False)
            Add an entry based on the given arguments
        dump_vcf(fn, cons=False, is_maj=False, append=False)
            Create a VCF file containing the variants in the BitRSCollection
                object
        get_by_chr_pos(chrom, pos)
            Return a list of (rsid, chrom, pos) at given chrom, pos and a list of
                corresponding entries
        get_by_rsid(rsid)
            Return a list of (rsid, chrom, pos) with given rsid and a list of
                corresponding entries
        get_major(mut=True)
            Return a RSCollection object containing all major alleles

        Static Methods
        --------------------
        chrom_to_int(c)
            Function to use to sort chromosomes (numeric < X < Y < M)
        decode_bit(bit_code)
            Convert from bit-packed int to list of variants
        encode_bit(var_list)
            Convert from list of variants to bit-packed int
        filter_vcf(rsc_dir, fn_in, fn_out, pop=None, log_fn=None, cons=False,
            keep_samps=False, quiet=False)
            Filter a given VCF file using ConsDB files stored in the given dir
        get_chrom_from_filename(fn)
            Helper method to parse a filename and find which chromosome it is
        load_from_file_full(fn, quiet=False)
            Load an entire ConsDB file
        load_from_file_pop(fn, pop, quiet=False)
            Load an entire ConsDB file, setting the allele frequency as the AF
            of the given population
        open(fn)
            Helper method to appropriately open the given file
        sort_rsidx(rsidx)
            Helper method to sort a list of rsidx
        var_list_to_vcf(rsidx, var_list, cons=False)
            Helper method to convert a list of variants to vcf entry(ies)
    """

    def __init__(self):
        # Keys are (rsid, chrom, pos)
        # All values in key are strings
        self.entries = {}

        # Dict of (chr, pos) -> rsidx or rsid -> rsidx
        # These make finding entries by position/rsid much faster.
        self.chr_pos_table = {}
        self.rsid_table = {}

    def __add__(self, rsc):
        """Implement addition for two RSCollection objects."""

        new_col = BitRSCollection()

        for r, e in self.entries.items():
            new_col.add_entry(r, e)

        for r, e in rsc.entries.items():
            new_col.add_entry(r, e)

        return(new_col)

    def __eq__(self, rsc):
        """Implement equality checking for two RSCollection objects."""

        if set(self.entries.keys()) != set(rsc.entries.keys()):
            return(False)

        for rsidx, e in self.entries.items():
            if e != rsc.entries[rsidx]:
                return(False)

        return(True)

    def __getitem__(self, key):
        """Allow access to entries via the [] operator."""

        return(self.entries[key])

    def __iadd__(self, rsc):
        """Implement incremental addition."""

        for rsidx, e in rsc.entries.items():
            self.add_entry(rsidx, e)

        return(self)

    def __len__(self):
        """Implement len operator."""

        return(len(self.entries))

    def add_entry(self, rsidx, e):
        """
        Add entry from a variant list.

        Parameters:
        rsidx: Tuple of (rsid, chrom, pos)
        e: List of variants (output from decode_bit)
        """

        if rsidx not in self.entries:
            var_list = e
            added = True
        else:
            var_list = BitRSCollection.decode_bit(self.entries[rsidx])
            added = False
            for var in e:
                if any([var[0] == v[0] and var[1] == v[1] \
                    for v in var_list]):
                    continue

                var_list += var
                added = True

        if added:
            self.entries[rsidx] = BitRSCollection.encode_bit(var_list)

            try:
                self.chr_pos_table[(rsidx[1], rsidx[2])].append(rsidx)
                self.rsid_table[rsidx[0]].append(rsidx)
            except KeyError:
                self.chr_pos_table[(rsidx[1], rsidx[2])] = [rsidx]
                self.rsid_table[rsidx[0]] = [rsidx]

    def add_entry_line(self, e):
        """
        Add entry from a line in a ConsDB file.

        Parameters:
        e: Line from a ConsDB file containing the entry to add
        """
        e = e.split(':')
        rsidx = tuple(e[0].split(','))

        if rsidx in self.entries:
            raise KeyError(f'{rsidx} already present.')

        e_split = e[1].split(';')
        
        maj_fn = lambda maj,minor: 0 if int(maj) > int(minor) \
            else (1 if int(minor) > int(maj) \
                else (2 if int(minor) == int(maj) == 0 else 3))
        e_split = [v.split(',')[:2] + [maj_fn(*v.split(',')[2:4])] \
            + v.split(',')[5:] for v in e_split \
            if len(v.split(',')[0]) < 50 and len(v.split(',')[1]) < 50 \
            and 'N' not in v.split(',')[0] and 'N' not in v.split(',')[1]]

        self.entries[rsidx] = BitRSCollection.encode_bit(e_split)
        
        try:
            self.chr_pos_table[(rsidx[1], rsidx[2])].append(rsidx)
            self.rsid_table[rsidx[0]].append(rsidx)
        except KeyError:
            self.chr_pos_table[(rsidx[1], rsidx[2])] = [rsidx]
            self.rsid_table[rsidx[0]] = [rsidx]

    def add_entry_from_args(self, chrom, rsid, pos, quiet=False):
        """
        Add an entry based on the given arguments.

        Parameters:
        chrom: Chromosome of the entry being added
        rsid: RefSNP ID of the entry being added
        pos: Position of the entry being added
        quiet: Suppress progress information being printed
        """

        rsid = str(rsid)
        chrom = str(chrom)
        pos = str(pos)

        rsidx = (rsid, chrom, pos)
        if rsidx in self.entries:
            if not quiet:
                print('Entry with chrom {} and RefSNP ID {} '
                      'already present.'.format(chrom, rsid))
            return

        self.entries[rsidx] = 0

        try:
            self.chr_pos_table[(chrom, pos)].append(rsidx)
            self.rsid_table[rsid].append(rsidx)
        except KeyError:
            self.chr_pos_table[(chrom, pos)] = [rsidx]
            self.rsid_table[rsid] = [rsidx]

    def dump_vcf(self, fn, cons=False, is_maj=False, append=False):
        """
        Create a VCF file containing the variants in the BitRSCollection object.
        If saving a consensus VCF file, makes sure that there is only one
        alternate allele per line. The append argument allows support for
        loading/saving in a piecewise fashion.

        Parameters:
        fn: File to save VCF to
        cons: If this VCF file is a consensus VCF
        is_maj: If this RSCollection object already contains only major alleles
        append: If appending to an existing VCF file
        """

        if not append:
            header=['##fileformat=VCFv4.3',
                f'##fileDate={time.strftime("%d%m%Y_%Hh%Mm%Ss")}',
                '##source=ConsDB',
                ('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele '
                    'frequency in the given population, in the range (0,1).">'),
                ('##INFO=<ID=VT,Number=A,Type=String,Description="Variant type '
                    'of the given alternate allele.">')]
            chroms = (sorted({k[0] for k in self.chr_pos_table.keys()},
                key=BitRSCollection.chrom_to_int))
            header.extend([f'##contig=<ID={c}>' for c in chroms])
            header.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')
            header = '\n'.join(header)

            read_mode = 'w'
        else:
            read_mode = 'a'

        if cons and not is_maj:
            maj = self.get_major()
            rsidx_list = sorted(maj.entries.keys(),
                key=BitRSCollection.sort_rsidx)
            entries = maj.entries
        else:
            rsidx_list = sorted(self.entries.keys(),
                key=BitRSCollection.sort_rsidx)
            entries = self.entries


        with open(fn, read_mode) as fp_out:
            if not append:
                fp_out.write(f'{header}\n')
            for rsidx in rsidx_list:
                vcf = BitRSCollection.var_list_to_vcf(rsidx,
                    BitRSCollection.decode_bit(entries[rsidx]), cons)
                fp_out.write(f'{vcf}\n')

    def get_by_chr_pos(self, chrom, pos):
        """
        Return a list of all (rsid, chrom, pos) combinations that match the
        given chrom and position, and a list of corresponding entries.

        Parameters:
        chrom: Chromosome to get
        pos: Position to get
        """

        if type(chrom) == int:
            chrom = str(chrom)
        if type(pos) == int:
            pos = str(pos)

        try:
            rsidx = self.chr_pos_table[(chrom, pos)]
        except KeyError:
            return([], [])

        return(rsidx, [self.entries[r] for r in rsidx])

    def get_by_rsid(self, rsid):
        """
        Return a list of all (rsid, chrom, pos) combinations that match the
        given RefSNP ID, and a list of corresponding entries.

        Parameters:
        rsid: RefSNP ID number to get
        """

        if type(rsid) == int:
            rsid = str(rsid)

        try:
            rsidx = self.rsid_table[rsid]
        except KeyError:
            return([], [])

        return(rsidx, [self.entries[r] for r in rsidx])


    def get_major(self, mut=True):
        """
        Return a BitRSCollection object containing all major alleles.

        Parameters:
        mut: Whether to only include variants that are different between
            reference and alternate allele
        """

        major_col = BitRSCollection()
        for rsidx, b in self.entries.items():
            var_list = [v for v in BitRSCollection.decode_bit(b) \
                if (v[4] or not mut) and \
                (v[2] == 0 or (v[2] == 2 and v[3] >= 0.5))]

            if len(var_list) > 0:
                major_col.add_entry(rsidx, var_list)

        return(major_col)

    @staticmethod
    def chrom_to_int(c):
        """
        Function to use to sort chromosomes (numeric < X < Y < M)

        Parameters:
        c: Chromosome to convert
        """

        try:
            return(int(c))
        except ValueError:
            c = c.upper()
            if c == 'X':
                return(23)
            if c == 'Y':
                return(24)
            if c == 'M' or c == 'MT':
                return(25)
            else:
                raise ValueError(f'Unknown chromosome {c}.')

    @staticmethod
    def decode_bit(bit_code):
        """
        Convert from bit-packed int to list of variants.

        Parameters:
        bit_code: int that can be unpacked into a list of variants
        """

        var_mask = (1 << 5) - 1 # 32-1 = 0b11111 (5 bits)
        n_base_mask = (1 << 6) - 1 # 64-1 = 0b111111 (6 bits)
        base_mask = (1 << 2) - 1 # 4-1 = 0b11 (2 bits)
        maj_mask = (1 << 2) - 1 # 4-1 = 0b11 (2 bits)
        af_mask = (1 << 10) - 1 # 0b1111111111 (10 bits)
        var_type_mask = (1 << 2) - 1 # 4-1 = 0b11 (2 bits)

        # Get number of variants encoded in list and shift bit_code to remove
        #  the used bits
        n_vars = bit_code & var_mask
        bit_code >>= 5

        var_list = []
        for _ in range(n_vars):
            # Read bit telling whether variant is a SNP (1) or an indel (0)
            snp = bit_code & 1
            bit_code >>= 1

            if snp:
                # If variant is a SNP, only have one base to read
                n_ref = 1
                n_var = 1
            else:
                # If variant is an indel, get indel reference length
                n_ref = bit_code & n_base_mask
                bit_code >>= 6

            # Build reference string by appending the nucleotide that matches
            #  the base-2 encoding for each encoded nucleotide
            ref = ''
            for _ in range(n_ref):
                ref += next(b for b, b_bin in BitRSCollection.BASE_ENC.items()\
                    if b_bin == (bit_code & base_mask))
                bit_code >>= 2

            if not snp:
                # Get indel variant length
                n_var = bit_code & n_base_mask
                bit_code >>= 6

            # Build variant string
            var = ''
            for _ in range(n_var):
                var += next(b for b, b_bin in BitRSCollection.BASE_ENC.items()\
                if b_bin == (bit_code & base_mask))
                bit_code >>= 2

            # Get major allele encoding (based on number of studies)
            # Key:
            #  0: Number of major studies > number of minor studies
            #  1: Number of minor studies > number of major studies
            #  2: No studies recorded for minor or major
            #  3: Same number of studies for major and minor (non-zero)
            maj = bit_code & maj_mask
            bit_code >>= 2

            # Get allele frequency (encoded to 3 decimal places)
            af = (bit_code & af_mask) / 1000
            bit_code >>= 10

            # Get variant type
            # Key:
            #  0: No variant (ref == var)
            #  1: SNP
            #  2: indel
            var_type = bit_code & var_type_mask
            bit_code >>= 2

            var_list.append((ref, var, maj, af, var_type))

        return(var_list)

    @staticmethod
    def encode_bit(var_list):
        """
        Create a bit-packed int from a list of variants.

        Parameters:
        var_list: List of variants
        """

        n_vars = len(var_list)
        bit_code = 0
        for var in var_list:
            ref_n_bases = len(var[0])
            ref_code = 0
            for b in ''.join(reversed(var[0])):
                ref_code <<= 2

                b_bin = BitRSCollection.BASE_ENC.get(b)
                if b_bin is None:
                    continue

                ref_code |= b_bin

            var_n_bases = len(var[1])
            var_code = 0
            for b in ''.join(reversed(var[1])):
                var_code <<= 2

                b_bin = BitRSCollection.BASE_ENC.get(b)
                if b_bin is None:
                    continue

                var_code |= b_bin

            snp = ref_n_bases == 1 and var_n_bases == 1

            maj = var[2]

            try:
                # This will only run/be necessary if the list passed in was
                #  parsed directly from a line in a ConsDB file. If encoding
                #  a previously encoded list, the af entry in the list should
                #  already be an int and won't need to be parsed.

                # Calculate all allele frequencies, then pick the appropriate
                #  one to be representative of the variant
                afs = np.asarray([int(a.split('/')[0]) / int(a.split('/')[1]) \
                    for a in var[3].split('.')])
                # Pick highest allele frequency if the allele is a major allele
                #  (based on studies). Otherwise, pick the largest af < 0.5
                if maj != 1:
                    af = int(max(afs) * 1000)
                else:
                    af = int(max(afs[afs < 0.5]) * 1000)
            except AttributeError:
                # If the af entry is already an int
                af = int(var[3] * 1000)

            if type(var[4]) == int:
                var_type = var[4]
            else:
                if var[4] == '':
                    var_type = 0
                elif var[4] == 'SNP':
                    var_type = 1
                else:
                    var_type = 2

            # Variant type
            bit_code <<= 2
            bit_code |= var_type

            # Allele frequency
            bit_code <<= 10
            bit_code |= af

            # Major/minor
            bit_code <<= 2
            bit_code |= maj

            # Store var bases
            bit_code <<= 2 * var_n_bases
            bit_code |= var_code

            if not snp:
                # Number of bases in var
                bit_code <<= 6
                bit_code |= var_n_bases

            # Store ref bases
            bit_code <<= 2 * ref_n_bases
            bit_code |= ref_code

            if not snp:
                # Number of bases in ref
                bit_code <<= 6
                bit_code |= ref_n_bases

            # SNP bit
            bit_code <<= 1
            bit_code |= snp

        # Number of vars
        bit_code <<= 5
        bit_code |= n_vars

        return(bit_code)

    @staticmethod
    def filter_vcf(rsc_dir, fn_in, fn_out, pop=None, log_fn=None, cons=False,
        keep_samps=False, quiet=False):
        """
        Filter a given VCF file using ConsDB files stored in the given dir.

        Remove all records that call a variant that is a major allele and write
        all other records.

        Parameters:
        rsc_dir: Directory containing ConsDB files
        fn_in: Input VCF file to filter
        fn_out: Output filename
        pop: Population to use
        log_fn: Optional file to store progress/log output
        cons: Making a consensus VCF (keep major alleles instead of minor)
        keep_samps: Write sample information
        quiet: Disable progress/log output to stdout
        """

        if log_fn:
            if os.path.exists(log_fn):
                os.remove(log_fn)

        # If making a VCF for a consensus genome, only want to keep one allele
        #  at each position, so keep track of positions we've added an allele at
        if cons:
            seen_pos = set()

        kept = 0
        removed = 0
        loaded = {}
        with open(fn_out, 'w') as fp:
            for line in BitRSCollection.open(fn_in):
                try:
                    line = line.decode()
                except AttributeError:
                    pass

                # Metadata line so just print it
                if line[:2] == '##':
                    fp.write(line)
                    continue

                # Get rid of individual sample information
                line = line.strip().split('\t')
                if not keep_samps:
                    line = line[:8]

                if line[0] == '#CHROM':
                    line = '\t'.join(line)
                    fp.write(f'{line}\n')
                    continue

                if line[0][0] == '#':
                    continue

                c = line[0].strip('chr')
                pos = line[1]
                # If we already have an allele at this position, skip to the
                #  next line
                if cons and (c, pos) in seen_pos:
                    continue

                rec_ref = line[3]
                rec_alts = line[4].split(',')
                rec_info = {}
                for inf in line[7].split(';'):
                    inf = inf.split('=')
                    k = inf[0]
                    try:
                        v = inf[1].split(',')
                    except IndexError:
                        v = 0
                    rec_info[k] = v
                line = '\t'.join(line)
                
                try:
                    rsc = loaded[c]
                    # Previously tried to load this chromosome but couldn't find
                    #  a file for it
                    if rsc is None:
                        fp.write(f'{line}\n')
                        continue
                except KeyError:
                    # If doing a pan consensus, first attempt to just load the
                    #  major alleles file, otherwise need to load whole file
                    if pop is None:
                        fn = f'{rsc_dir}/chr{c}_maj.gz'
                        maj_file = True
                        if not os.path.exists(fn):
                            fn = f'{rsc_dir}/chr{c}_rscol.gz'
                            maj_file = False
                    else:
                        fn = f'{rsc_dir}/chr{c}_rscol.gz'

                    try:
                        if pop is None:
                            rsc = BitRSCollection.load_from_file_full(fn,
                                quiet)
                        else:
                            rsc = BitRSCollection.load_from_file_pop(fn,
                                pop, quiet)
                    except FileNotFoundError:
                        print((f'No file found for chr {c}, keeping '
                            'all alleles.'))
                        loaded[c] = None
                        fp.write(f'{line}\n')
                        continue

                    if pop or not maj_file:
                        rsc = rsc.get_major()

                    loaded[c] = rsc

                _, e_list = rsc.get_by_chr_pos(c, pos)
                if cons and len(e_list) == 0:
                    # Want to keep major alleles but this rec wasn't in ConsDB
                    #  file, so not a known major allele
                    if log_fn:
                        with open(log_fn, 'a') as fp_log:
                            fp_log.write(f'Removing rec {line}\n')
                    if not quiet:
                        print(f'Removing rec {line}', flush=True)
                    
                    removed += len(rec_alts)
                    continue

                # Set of (ref, var) for every major allele at this position
                maj_alleles = {tuple(v[:2]) for e in e_list \
                    for v in BitRSCollection.decode_bit(e)}

                if cons:
                    # Index of alts that are either not major alleles or are
                    #  equal to the reference sequence
                    rec_idx = np.asarray([(rec_ref, a) not in maj_alleles \
                        or a == rec_ref for a in rec_alts])
                    # Only want to take the first major allele
                    rec_idx[np.argmin(rec_idx)+1:] = True
                else:
                    # Index of alts that are either major alleles or are equal
                    #  to the reference sequence
                    rec_idx = [(rec_ref, a) in maj_alleles \
                        or a == rec_ref for a in rec_alts]

                # Either this (chr, pos) combination is/isn't a major allele
                #  or none of the alternate alleles are in the ConsDB file
                if sum(rec_idx) == 0:
                    # Write the record to the output file and inform the user
                    if log_fn:
                        with open(log_fn, 'a') as fp_log:
                            fp_log.write(f'Keeping rec {line}\n')
                    if not quiet:
                        print(f'Keeping rec {line}', flush=True)

                    if cons:
                        seen_pos.add((c, pos))
                    fp.write(f'{line}\n')
                    kept += len(rec_alts)
                    continue

                # All alternate alleles in this record are/aren't major alleles,
                #  so get rid of the whole record
                if sum(rec_idx) == len(rec_idx):
                    if log_fn:
                        with open(log_fn, 'a') as fp_log:
                            fp_log.write(f'Removing rec {line}\n')
                    if not quiet:
                        print(f'Removing rec {line}', flush=True)
                    
                    removed += len(rec_alts)
                    continue

                # Only keep the alts that are/aren't major alleles
                rec_alts = [rec_alts[i] for i in range(len(rec_alts)) \
                    if not rec_idx[i]]
                # Need to also edit the INFO dict to remove info corresponding
                #  to removed entries
                for k, v in rec_info:
                    if not type(v) is list:
                        continue
                    rec_info[k] = [v[i] for i in range(len(v)) \
                        if not rec_idx[i]]

                line = line.split('\t')
                line[4] = rec_alts
                line[7] = ';'.join([f'{k}={",".join(v)}' if type(v) is list \
                    else k for k,v in rec_info.items()])

                if cons:
                    seen_pos.add((c, pos))
                fp.write(f'{line}\n')
                kept += len(rec_alts) - sum(rec_idx)
                removed += sum(rec_idx)

                if log_fn or not quiet:
                    msg = 'Filtered the following alleles:\n'
                    for a in [rec_alts[i] for i in range(len(rec_alts)) \
                        if rec_idx[i]]:
                        msg += f'\tREF: {rec_ref}, ALT: {a}\n'
                    msg += f'Kept rec {line}\n'
                    if log_fn:
                        with open(log_fn, 'a') as fp_log:
                            fp_log.write(msg)
                    if not quiet:
                        print(msg[:-1], flush=True)

            if log_fn:
                with open(log_fn, 'a') as fp_log:
                    fp_log.write(f'Kept {kept} alleles.\n')
                    fp_log.write(f'Removed {removed} alleles.\n')
            print(f'Kept {kept} alleles.', flush=True)
            print(f'Removed {removed} alleles.', flush=True)

    @staticmethod
    def get_chrom_from_filename(fn):
        """
        Helper method to parse a filename and find which chromosome it is.

        Parameters:
        fn: Filename to use
        """

        match = 'chr([0-9XY]+)'
        try:
            chrom = re.search(match, fn).group(1)
        except AttributeError:
            print('Unable to interpret chromosome from file name ({}), '
                           'using -1 as a placeholder.'.format(fn))
            chrom = -1
        
        return(chrom)

    @staticmethod
    def load_from_file_full(fn, quiet=False):
        """
        Load an entire ConsDB file.

        Parameters:
        fn: Filename to load
        quiet: Disable progress/log output to stdout
        """

        rsc = BitRSCollection()
        add_entry_line = rsc.add_entry_line
        if not os.path.isfile(fn):
            print(f'{fn} not found.')
            return(None)
        if 'maj' in fn:
            print_ln = 10000
        else:
            print_ln = 2500000
        c = BitRSCollection.get_chrom_from_filename(fn)
        with BitRSCollection.open(fn) as fp:
            for n, line in enumerate(fp):
                if n % print_ln == 0 and not quiet:
                    print(f'chr {c}: Reading line {n}.', flush=True)
                # Need to decode line from bytes to string if the file is
                #  compressed
                try:
                    line = line.decode('utf8').strip()
                except AttributeError:
                    pass

                add_entry_line(line)

        return(rsc)

    @staticmethod
    def load_from_file_pop(fn, pop, quiet=False):
        """
        Load an entire ConsDB file, setting the allele frequency as the AF of
        the given population.
        Skip variants that are not present in the given population.

        Parameters:
        fn: Filename to use
        pop: Population to load
        quiet: Disable progress/log output to stdout
        """

        rsc = BitRSCollection()
        add_entry_line = rsc.add_entry_line
        if not os.path.isfile(fn):
            print(f'{fn} not found.')
            return(None)
        c = BitRSCollection.get_chrom_from_filename(fn)
        for n, line in enumerate(BitRSCollection.open(fn)):
            if n % 2500000 == 0 and not quiet:
                print(f'chr {c}: Reading line {n}.', flush=True)
            # Need to decode line from bytes to string if the file is
            #  compressed
            try:
                line = line.decode('utf8').strip()
            except AttributeError:
                pass

            line = line.split(':')
            var_list = np.asarray(line[1].split(';'))
            # Check for given population in the pop afs section of the alleles
            afs = np.asarray([re.search(f'{pop}~([0-9/.]+)', v) \
                for v in var_list])
            
            # Get rid of alleles that don't contain an af for the pop
            idx = [m is not None for m in afs]
            var_list = var_list[idx]

            afs = [m.group(1) for m in afs[idx]]
            for i in range(len(var_list)):
                temp_var = var_list[i].split(',')
                # Replace AF with pop AF
                temp_var[5] = afs[i]
                var_list[i] = ','.join(temp_var)

            line = (f'{line[0]}:{";".join(var_list)}')
            if line[-1] != ':':
                add_entry_line(line)

        if len(rsc.entries) == 0:
            rsc = None
        return(rsc)


    @staticmethod
    def open(fn):
        """
        Helper method to appropriately open the given file.

        Parameters:
        fn: Filename to open
        """

        if fn.split('.')[-1] == 'bz2':
            return(bz2.open(fn, 'rb'))
        elif fn.split('.')[-1] == 'gz':
            return(gzip.open(fn, 'rb'))

        return(open(fn, 'r'))

    @staticmethod
    def sort_rsidx(rsidx):
        """
        Helper method to sort a list of rsidx. Meant to be used as a
        key for sorting. Returns the chromosome and position of the rsidx, both
        in int form.

        Parameters:
        rsidx: List/tuple of (rsid, chrom, pos)
        """

        _, chrom, pos = rsidx
        chrom = BitRSCollection.chrom_to_int(chrom)
        pos = int(pos)

        return(chrom, pos)

    @staticmethod
    def var_list_to_vcf(rsidx, var_list, cons=False):
        """
        Helper method to convert a list of variants to vcf entry(ies).

        Parameters:
        rsidx: List/tuple of (rsid, chrom, pos)
        var_list: List of variants (as returned from BitRSCollection.decode_bit)
        cons: Making a consensus VCF file
        """

        vcf_list = []
        for ref, vl in itertools.groupby(sorted(var_list),
            key=operator.itemgetter(0)):
            vl = list(vl)
            if cons:
                vl = [vl[0]]
            chrom = rsidx[1]
            pos = rsidx[2]
            rsid = f'rs{rsidx[0]}' if rsidx[0] != '-1' else '.'
            alts = ','.join([v[1] for v in vl])
            qual = '.'
            filt = 'PASS'
            afs = ','.join([f'{v[3]:0.3f}' for v in vl])
            vts = ','.join(['SNP' if len(ref) == len(v[1]) else \
                ('ins' if len(ref) < len(v[1]) else 'del') for v in vl])
            info = f'AF={afs};VT={vts}'

            vcf_list.append('\t'.join([chrom, pos, rsid, ref, alts, qual, filt,
                info]))
            # Only want one line for each position
            if cons:
                break

        return('\n'.join(vcf_list))