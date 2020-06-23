import bz2
import gzip
import itertools
import json
import numpy as np
import os
import re
import time

class RSEntry:
    """
    A class to hold all variants at a given position.

    Attributes
    --------------------
    chrom : str
        Which chromosome the entry is on (1-22, X, Y)
    rsid : int
        The RefSNP ID number for the entry
    pos : int
        Base position of the entry (1-indexed)
    vars : dict
        Dictionary that maps variant code strings to RSVariant objects

    Methods
    --------------------
    add_var(var)
        Add a RSVariant object to the RSEntry
    add_var_from_args(pos, ref, var, major=0, minor=0, clin=[], afs=[],
        var_type='', pop_afs={})
        Add a variant to the RSEntry based on the given arguments
    all_vars_empty()
        Check if all variants in the RSEntry are empty
    get_major_alleles()
        Return a list of major RSVariants
    get_major_alleles_pop(pop)
        Return a list of major RSVariants for a given population
    pick_major_allele(var_list)
        Reproducibly choose a representative variant
    to_vcf(pop=None, cons=False, is_maj=False)
        Create a VCF record from variants
    """
    
    def __init__(self, chrom, rsid, pos=None):
        """
        Initialize an instance of the RSEntry class. Requires a chromosome and
        a RSID number. Base position is optional, and if not passed will be
        inferred from the first added variant.

        Parameters:
        chrom: Chromosome of the entry
        rsid: RefSNP ID number of the entry
        pos: Position of the entry
        """

        self.chrom = chrom
        self.rsid = int(rsid)
        try:
            self.pos = int(pos)
        except TypeError:
            self.pos = pos
        self.vars = {}

    def __add__(self, rse):
        """Implement addition for two RSEntry objects."""

        # Check that chromosome and rsid are the same
        if self.chrom != rse.chrom:
            raise ValueError('Given operands have different chromosomes.')
        if self.rsid != rse.rsid:
            raise ValueError('Given operands have different RefSNP ID numbers.')

        # Create new RSEntry object and add all variants from both operands
        new_rse = RSEntry(self.chrom, self.rsid)
        for v in self.vars.values():
            new_rse.add_var(v)
        for v in rse.vars.values():
            new_rse.add_var(v)

        return(new_rse)

    def __eq__(self, rse):
        """Implement equality checking for two RSEntry objects."""

        # Check that chromosome and rsid are the same, and that both operands
        #  contain all of the same variants
        if self.chrom != rse.chrom or self.rsid != rse.rsid or \
            set(self.vars.keys()) != set(rse.vars.keys()):
            return(False)

        for vc, v in self.vars.items():
            if v != rse.vars[vc]:
                return(False)
        return(True)

    def __getitem__(self, key):
        """Allow access to variants via the [] operator."""

        return(self.vars[key])
    
    def __iadd__(self, rse):
        """Implement incremental addition."""

        # Check that chromosome and rsid are the same
        if self.chrom != rse.chrom:
            raise ValueError('Given operands have different chromosomes.')
        if self.rsid != rse.rsid:
            raise ValueError('Given operands have different RefSNP ID numbers.')

        for v in rse.vars.values():
            self.add_var(v)

        return(self)

    def __len__(self):
        """Implement len operator."""

        return(len(self.vars))

    def __repr__(self):
        """Implement repr operator."""

        base_string = f'{self.rsid},{self.chrom},{self.pos}:'
        string = base_string
        string += ';'.join([repr(v) for v in self.vars.values()])

        return(string)            

    def __str__(self):
        """Implement str operator."""

        string = f'rs{self.rsid}:\n'
        for v in self.vars.values():
            string += str(v)
            string += '\n'
            
        return(string)
            
    def add_var(self, var):
        """
        Add a variant by passing an RSVariant object. Note that this creates a
        new RSVar object and doesn't use the passed object.
    
        Parameters:
        var: RSVar object to add
        """

        self.add_var_from_args(var.pos, var.ref, var.var, var.major, var.minor,
                               var.clin, var.afs, var.var_type, var.pop_afs)

    def add_var_from_args(self, pos, ref, var, major=0, minor=0, clin=[],
        afs=[], var_type='', pop_afs={}):
        """
        Add a variant by passing all necessary information.

        Parameters:
        pos: Variant position
        ref: Reference sequence at the variant position
        var: Variant sequence
        major: Number of studies in dbSNP listing this variant as major
        minor: Number of studies in dbSNP listing this variant as minor
        clin: Clinical significance of this variant
        afs: Allele frequency of this allele
        var_type: Variant type (SNP or indel) of this allele
        pop_afs: Population allele frequencies
        """

        # Set position value for the class if one hasn't already been set.
        # If the position has been set and it doesn't match the given position,
        #  raise an error.
        if self.pos is None:
            self.pos = pos
        elif self.pos != pos:
            raise ValueError(('Given position does not match the position of '
                'the RSEntry object.'))

        vc = f'{pos}{ref}->{var}'

        # If we don't have a variant matching the given var info, make a new
        #  RSVar object and add it to the dict. Otherwise, append all the new
        #  information to the existing var.
        v = self.RSVar(pos, ref, var, major, minor, clin, afs, var_type,
            pop_afs)
        try:
            self.vars[vc] += v
        except KeyError:
            self.vars[vc] = v

    def all_vars_empty(self):
        """Check if all variants are empty."""

        return(all([v.is_empty() for v in self.vars.values()]))

    def get_major_alleles(self):
        """Get a list of variants representing a major allele (AF >= 0.5)."""

        if len(self.vars) == 0:
            return([])

        major_vars = [v for v in self.vars.values() \
            if v.is_major()]

        # If no variants have af >= 0.5, take the one with the highest allele
        #  frequency as the major allele.
        if len(major_vars) == 0:
            afs = np.asarray([max(v.calc_afs()) for v in self.vars.values()])
            if np.all(afs == 0.5) or np.all(afs < 0.5):
                return([])

            major_vars = [v for v in self.vars.values() \
                if max(v.calc_afs()) == max(afs)]

        return(major_vars)

    def get_major_alleles_pop(self, pop):
        """
        Get a list of variants representing a major allele (AF >= 0.5) for the
        given population.

        Parameters:
        pop: Population to use
        """

        if len(self.vars) == 0:
            return([])

        major_vars = [v for v in self.vars.values() if v.pop_af(pop) >= 0.5]
        return(major_vars)

    def pick_major_allele(self, var_list):
        """
        Reproducibly choose a representative variant.

        Parameters:
        var_list: List of RSVar objects to pick from
        """

        # Sort by af
        # If there are ties, prefer SNPs over indels and pur<->pur over
        #  pur<->pyr
        afs = np.asarray([v.af() for v in var_list])
        af_idx = np.flip(np.argsort(afs))
        afs = afs[af_idx]
        sorted_vars = np.asarray(var_list)[af_idx]

        max_af = afs[0]
        if np.sum(afs == max_af) == 1:
            return(sorted_vars[0])

        filt_vars = [v for v in sorted_vars if v.var_type == 'SNP']
        if len(filt_vars) == 0:
            filt_vars = sorted_vars
        elif len(filt_vars) == 1:
            return(filt_vars[0])
        else:
            filt_vars = sorted(filt_vars,
                key=lambda v: (v.ref in {'A','T'} == v.var in {'A','T'}),
                reverse=True)
            return(filt_vars[0])

        # If only indels, choose by smallest number if inserted/deleted bases
        filt_vars = sorted(filt_vars, key=lambda v: abs(len(v.ref)-len(v.var)))
        return(filt_vars[0])

    def to_vcf(self, pop=None, cons=False, is_maj=False):
        """
        Create a VCF record from variants. If preparing a consensus VCF
        file, makes sure that there is only one alternate allele per line.

        Parameters:
        pop: Population to use
        cons: If this VCF file is a consensus VCF
        is_maj: If this RSCollection object already contains only major alleles
        """

        if cons and not is_maj:
            if pop:
                var_list = self.get_major_alleles_pop(pop)
            else:
                var_list = self.get_major_alleles()
        else:
            var_list = list(self.vars.values())

        vcf_list = []
        chrom = self.chrom
        pos = str(self.pos)
        rsid = f'rs{self.rsid}' if self.rsid != -1 else '.'
        for ref, vl in itertools.groupby(sorted(var_list,
            key=lambda v: (v.ref, v.var)), key=lambda v: v.ref):
            vl = list(vl)
            if cons:
                vl = [vl[0]]
            alts = ','.join([v.var for v in vl])
            qual = '.'
            filt = 'PASS'
            afs = np.asarray([v.pop_af(pop) if pop else v.af() for v in vl])
            afs = np.trunc(afs * 1000) / 1000
            afs = ','.join([f'{a:0.3f}' for a in afs])
            vts = ','.join(['SNP' if len(ref) == len(v.var) else \
                ('ins' if len(ref) < len(v.var) else 'del') for v in vl])
            info = f'AF={afs};VT={vts}'

            vcf_list.append('\t'.join([chrom, pos, rsid, ref, alts, qual, filt,
                info]))
            # Only want one line for each position
            if cons:
                break

        return('\n'.join(vcf_list))
    
    class RSVar:
        """
        A sub-class to represent a variant.

        Attributes
        --------------------
        pos : int
            Base position of the variant
        ref : str
            Reference base(s)
        var : str
            Variant base(s)
        major : int
            Number of studies that show this variant with AF >= 0.5
        minor : int
            Number of studies that show this variant with AF < 0.5
        clin : list
            List of clinical relevances
        afs : list
            List of allele frequencies, stores as tuples of 
            (# patients showing var, # patients in study)
        var_type : str
            What type of variant (SNP, indel, or blank)
        pop_afs : dict
            Dict containing population allele frequencies

        Methods
        --------------------
        af()
            Get a representative allele frequency for the variant
        calc_afs()
            Calculate float allele frequencies from afs list
        calc_pop_afs(pop)
            Calculate float allele frequencies for the given population
        is_major()
            Return if the variant represents a major allele
        is_empty()
            Return if the variant has no studies supporting it
        pop_af(pop)
            Get a representative population allele frequency for the variant
            and given population.
        var_code()
            Get the representative variant code for the variant
        """

        def __init__(self, pos, ref, var, major=0, minor=0, clin=[], afs=[],
            var_type='', pop_afs={}):
            """
            Initialize an instance of the RSVar class. Requires a position,
            reference allele, and variant allele. All other information is 
            oprtional, and will either be left blank or inferred (in the case
            of var_type).

            Parameters:
            pos: Variant position
            ref: Reference sequence at the variant position
            var: Variant sequence
            major: Number of studies in dbSNP listing this variant as major
            minor: Number of studies in dbSNP listing this variant as minor
            clin: Clinical significance of this variant
            afs: Allele frequency of this allele
            var_type: Variant type (SNP or indel) of this allele
            pop_afs: Population allele frequencies
            """

            if pos is None or ref is None or var is None:
                raise TypeError('Values required for pos, ref, and var.')

            # pos is 1-indexed
            self.pos = pos
            self.ref = ref
            self.var = var
            self.major = major
            self.minor = minor
            if clin == ['']:
                clin = []
            self.clin = clin
            self.afs = afs
            self.var_type = var_type
            if self.var_type is None:
                if len(self.ref) == len(self.var):
                    var_type = 'SNP'
                else:
                    var_type = 'indel'
            # Dict of population -> af
            self.pop_afs = pop_afs

        def __add__(self, v):
            """Implement addition for two RSVar objects."""

            if self.var_code() != v.var_code():
                raise ValueError('Given variants are incompatible.')

            major = self.major + v.major
            minor = self.minor + v.minor
            clin = list(set(self.clin).union(v.clin))
            afs = list(set(self.afs).union(v.afs))
            pop_afs = {pop: afs for pop,afs in self.pop_afs.items()}
            for pop, p_afs in v.pop_afs.items():
                try:
                    pop_afs[pop] += [a for a in p_afs if a not in pop_afs[pop]]
                except KeyError:
                    pop_afs[pop] = p_afs

            return(RSEntry.RSVar(self.pos, self.ref, self.var, major, minor,
                clin, afs, self.var_type, pop_afs))

        def __eq__(self, v):
            """Implement equality checking for two RSVar objects."""

            # Check that all parameters are the same
            if self.var_code() != v.var_code() or self.major != v.major or \
                self.minor != v.minor or set(self.clin) != set(v.clin) or \
                set(self.afs) != set(v.afs) or self.var_type != v.var_type or \
                set(self.pop_afs) != set(v.pop_afs):
                return(False)

            # Check that the af list for each pop is the same
            for pop, afs in self.pop_afs.items():
                if set(afs) != set(v.pop_afs[pop]):
                    return(False)

            return(True)

        def __iadd__(self, v):
            """Implement incremental addition."""

            if self.var_code() != v.var_code():
                raise ValueError('Given variants are incompatible.')

            self.major += v.major
            self.minor += v.minor
            self.clin += [c for c in v.clin if c not in self.clin]
            self.afs += [a for a in v.afs if a not in self.afs]
            for pop, afs in v.pop_afs.items():
                try:
                    self.pop_afs[pop] += [a for a in afs \
                        if a not in self.pop_afs[pop]]
                except KeyError:
                    self.pop_afs[pop] = afs

            return(self)

        def __repr__(self):
            """Implement repr operator."""

            string = ','.join([self.ref, self.var,
                str(self.major), str(self.minor),
                '.'.join(self.clin),
                '.'.join([f'{af[0]}/{af[1]}' for af in self.afs]),
                self.var_type,
                '|'.join([f'{pop}~{".".join([f"{a[0]}/{a[1]}" for a in af])}' \
                    for pop, af in self.pop_afs.items()])])
            return(string)

        def __str__(self):
            """Implement str operator."""

            string = '\t'.join([str(self.pos), f'{self.ref}->{self.var}',
                str(self.major), str(self.minor), str(self.af()),
                ', '.join(self.clin), self.var_type])
            return(string)

        def af(self):
            """Get a representative allele frequency for the variant."""

            afs = self.calc_afs()
            if len(afs) == 0:
                return(None)
            if self.is_major() or self.major == self.minor:
                return(max(afs))

            return(max(afs[afs < 0.5]))

        def calc_afs(self):
            """Calculate float allele frequencies from afs list."""

            return(np.asarray([a[0]/a[1] for a in self.afs]))

        def calc_pop_afs(self, pop):
            """
            Calculate float allele frequencies from afs list.

            Parameters:
            pop: Population to use
            """

            return(np.asarray([a[0]/a[1] for a in self.pop_afs[pop]]))

        def is_major(self):
            """Return if the variant represents a major allele."""

            maj_afs = [c for c in self.calc_afs() if c >= 0.5]
            return(self.major > self.minor or \
                (self.major == self.minor and len(maj_afs) > 0) or \
                len(maj_afs) > len(self.afs)/2)
            
        def is_empty(self):
            """Return if the variant has no studies supporting it."""

            return(self.major == self.minor == 0)

        def pop_af(self, pop):
            """
            Get a representative population allele frequency for the variant
            and given population.

            Parameters:
            pop: Population to use
            """

            afs = self.calc_pop_afs(pop)
            return(max(afs))

        def var_code(self):
            """Get the representative variant code for the variant."""

            return(f'{self.pos}{self.ref}->{self.var}')

class RSCollection:
    """
        A class to represent a collection of entries.

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
        add_entry(e)
            Add entry from an RSEntry object
        add_entry_from_args(chrom, rsid, pos, quiet=False)
            Add an entry based on the given arguments
        dump(fn, idx_file, c, rsids = ['all'], old_size=0, append=False,
            chunksize=10000)
            Save entries for the given rsids to the given file
        dump_full(fn)
            Quickly save the whole RSCollection object to the given file
        dump_chrs(chrs, fp_out, store_all=False, store_maj=False)
            Save entries for the given chromosomes
        dump_vcf(fn, pop=None, cons=False, is_maj=False)
            Save entries as a VCF file
        get_by_chr(chrom)
            Return a RSCollection object containing all entries from given chrom
        get_by_chr_pos(chrom, pos)
            Return a RSCollection object containing all entries at given chrom
            and given position
        get_by_rsid(rsid)
            Return a list of RSEntry objects that match the given rsid
        get_major(mut=True)
            Return a RSCollection object containing all major alleles

        Static Methods
        --------------------
        chrom_to_int(c)
            Function to use to sort chromosomes (numeric < X < Y < M)
        from_1000gp(fn, index_fn, superpop_fn, quiet=False)
            Create a new RSCollection object from a 1000 Genomes Project VCF file
        from_gnomad(fn, quiet=False)
            Create a new RSCollection object from a gnomAD VCF file
        from_dbsnp(fn, quiet=False)
            Create a new RSCollection object from a dbSNP JSON file
        get_chrom_from_filename(fn)
            Helper method to parse a filename and find which chromosome it is
        load_from_file_by_chr_pos(fn, idx_file, c, pos, chunk_idx_dict=None,
            ret_chunk=False)
            Find an entry in the given file using given chromosome and position
        load_from_file_by_rsid(fn, idx_file, rsid, chunk_idx_dict=None,
            ret_chunk=False)
            Find an entry in the given file using given rsid
        load_from_file_full(fn)
            Load an entire ConsDB file
        load_from_file_pops(fn, pops, cons=False)
            Only load variants that are major variants for pop
        make_chunk_idx_dict(idx_fn, key_fields=[1,2])
            Make a dict representation of the idx file
        merge_files(fn_list, out_fn, c=None)
            Merge multiple ConsDB files into one
        open(fn)
            Helper method to appropriately open the given file
        parse_file_line(d)
            Create a RSEntry object from given line from ConsDB file
        sort_rsidx(rsidx)
            Helper method to sort a list of rsidx
        sort_rsidx_line(line_split)
            Helper method to sort lines from a ConsDB file
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

        new_col = RSCollection()

        for e in self.entries.values():
            new_col.add_entry(e)

        for e in rsc.entries.values():
            new_col.add_entry(e)

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

        for e in rsc.entries.values():
            self.add_entry(e)

        return(self)

    def __len__(self):
        """Implement len operator."""

        return(sum([len(e) for e in self.entries.values()]))

    def __repr__(self):
        """Implement repr operator."""

        string = ''
        for rse in self.entries.values():
            string += repr(rse)
            string += '\n'
        return(string)

    def __str__(self):
        """Implement str operator."""

        string = 'rsid:\n'
        string += '\t'.join(['Position', 'Variant', '#_Major', '#_Minor',
                             'Clin_Rel', 'Var_Type'])
        string += '\n'
        for e in self.entries.values():
            string += str(e)

        return(string)

    def add_entry(self, e):
        """
        Add entry from an RSEntry object. Note that this does not create a new
        RSEntry object.

        Parameters:
        e: RSEntry object
        """

        rsid = str(e.rsid)
        pos = str(e.pos)
        entries_to_add = [e]
        # If either the entry to add or at least one of the entries in self
        #  has an rsid of -1 and the other entry has a non -1 rsid, then
        #  attempt to correct the -1 rsid
        if (e.chrom, pos) in self.chr_pos_table and (rsid == '-1' != \
            any(r[0] == '-1' for r in self.chr_pos_table[(e.chrom, pos)])):
            chr_pos_list = self.chr_pos_table[(e.chrom, pos)]
            for i in range(len(chr_pos_list)):
                rsidx = chr_pos_list[i]
                test_e = self.entries[rsidx]
                # Can only fix rsid if the two entries have the same (ref, var),
                #  so find all such combinations
                shared = set(test_e.vars.keys()).intersection(
                    set(e.vars.keys()))
                if len(shared) > 0:
                    if rsidx[0] == '-1' and rsid == '-1':
                        continue
                    if rsidx[0] == '-1':
                        test_e.rsid = int(rsid)
                        # Remove the correscted RSEntry from internal data
                        #  structures and add it to the list to be added later
                        self.entries.pop(rsidx)
                        self.chr_pos_table[(e.chrom, pos)].pop(i)
                        if len(self.chr_pos_table[(e.chrom, pos)]) == 0:
                            self.chr_pos_table.pop((e.chrom, pos))
                        entries_to_add.append(test_e)

                    elif rsid == '-1':
                        rsid = test_e.rsid
                        e.rsid = test_e.rsid
                    else:
                        continue
                    print(f'rsid corrected at {e.chrom}, {pos}.',
                        flush=True)
                    break

        rsidx = (rsid, e.chrom, pos)
        for rse in entries_to_add:
            try:
                self.entries[rsidx] += rse
                print(f'{rsidx} appended to.', flush=True)
            except KeyError:
                self.entries[rsidx] = rse

        try:
            if rsidx not in self.chr_pos_table[(e.chrom, pos)]:
                self.chr_pos_table[(e.chrom, pos)].append(rsidx)
                self.rsid_table[rsid].append(rsidx)
        except KeyError:
            self.chr_pos_table[(e.chrom, pos)] = [rsidx]
            self.rsid_table[rsid] = [rsidx]

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
        pos = str(pos)
        rsidx = (rsid, chrom, pos)
        # If an entry with the given args already exists, no need to create a
        #  new empty one
        if rsidx in self.entries:
            if not silent:
                print(f'Entry with chrom {chrom} and RefSNP ID {rsid} '
                      'already present.')
            return

        # Add the new empty RSEntry and populate the rsidx into the
        #  lookup tables
        self.entries[rsidx] = RSEntry(chrom, rsid, pos)
        try:
            self.chr_pos_table[(e.chrom, e.pos)].append(rsidx)
            self.rsid_table[e.rsid].append(rsidx)
        except KeyError:
            self.chr_pos_table[(e.chrom, e.pos)] = [rsidx]
            self.rsid_table[e.rsid] = [rsidx]

    def dump(self, fn, idx_file, c, rsids=['all'], old_size=0, append=False,
        chunksize=10000):
        """
        Save entries for the given rsids to the given file.

        Split the entries to write into multiple chunks, then write each chunk
        separately to the gzip file. This allows for random access when loading.

        Parameters:
        fn: Filename to save to
        idx_file: Filename to save index information to
        c: Which chromosome is being saved
        rsids: Which RefSNP IDs/(chr, pos) to save
        old_size: Prior size of the output file (used for random access)
        append: Whether to append to the given output file
        chunksize: Number of entries to write at a time
        """

        # Use GZIP compression algorithm for reduced compression ratio

        if np.array_equal(rsids, ['all']):
            keys = list(self.entries.keys())
        else:
            # If the entries in the rsids list are lists/tuples of length 2,
            #  consider them to be (chr, pos) combinations
            if len(next(iter(rsids))) == 2:
                keys = [r for chr_pos, rsidx in self.chr_pos_table.items() \
                    for r in rsidx if chr_pos in rsids]
            else:
                keys = [k for k in self.entries.keys() if k in rsids]

        # If there are no entries to write, just exit
        if len(keys) == 0:
            return(old_size)

        # Group keys into lists of chunksize
        keys = [keys[i:i+chunksize] for i in range(0, len(keys), chunksize)]

        # File to store compressed chunk indexes
        if idx_file is None:
            idx_file = f'{fn}.idx'

        # Remove target files if they exist so we can open in append mode
        if not append:
            if os.path.exists(fn):
                os.remove(fn)
            if os.path.exists(idx_file):
                os.remove(idx_file)

        new_size = 0
        print(f'chr {c}: {len(keys)} chunks to write of chunksize {chunksize}.',
            flush=True)
        for key_chunk in keys:
            # Write all entries in this chunk
            with gzip.open(fn, 'ab') as fp:
                write_str = '\n'.join([repr(self.entries[k]) \
                    for k in key_chunk])
                write_str += '\n'
                fp.write(write_str.encode())

            # Write line in index file
            # All rsid,chrom,pos in the chunk, separated by ";", the number of
            #  bytes in this chunk, and the start position of this chunk
            new_size = os.stat(fn).st_size
            with open(idx_file, 'a') as fp:
                k_str = [",".join([str(i) for i in k]) for k in key_chunk]
                fp.write((f'{";".join(k_str)} '
                    f'{new_size-old_size} {old_size}\n'))

            old_size = new_size

        return(new_size)

    def dump_full(self, fn):
        """
        Quickly save the whole RSCollection object to the given file.

        Parameters:
        fn: Filename to save to
        """

        with gzip.open(fn, 'wb') as fp:
            fp.write(repr(self).encode())

    def dump_chrs(self, chrs, fp_out, store_all=False, store_maj=False):
        """
        Save entries for the given chromosomes.

        Parameters:
        chrs: Chromosomes to save
        fp_out: Directory to store the files
        store_all: Save file with all entries
        store_maj: Save file with only major alleles
        """

        for c in chrs:
            rsc = self.get_by_chr(c)
            
            if store_all:
                rsc.dump(f'{fp_out}/chr{c}_rscol.gz')
            if store_maj:
                rsc.get_major().dump(f'{fp_out}/chr{c}_maj.gz')


    def dump_vcf(self, fn, pop=None, cons=False, is_maj=False):
        """
        Save entries as a VCF file.

        Parameters:
        fn: Filename to save to
        pop: Population to use
        con: Storing a consensus VCF
        is_maj: RSCollection already contains only major alleles
        """

        header=['##fileformat=VCFv4.3',
            f'##fileDate={time.strftime("%d%m%Y_%Hh%Mm%Ss")}',
            '##source=ConsDB',
            ('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency '
                'in the given population, in the range (0,1).">'),
            ('##INFO=<ID=VT,Number=A,Type=String,Description="Variant type of '
                'the given alternate allele.">')]
        chroms = (sorted({k[0] for k in self.chr_pos_table.keys()},
            key=RSCollection.chrom_to_int))
        header.extend([f'##contig=<ID={c}>' for c in chroms])
        header.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')

        # Get major alleles if necessary
        if cons and not is_maj:
            maj = self.get_major()
            rsidx_list = sorted(maj.entries.keys(),
                key=RSCollection.sort_rsidx)
            entries = maj.entries
            is_maj = True
        else:
            rsidx_list = sorted(self.entries.keys(),
                key=RSCollection.sort_rsidx)
            entries = self.entries

        with open(fn, 'w') as fp_out:
            header = '\n'.join(header)
            fp_out.write(f'{header}\n')
            for rsidx in rsidx_list:
                fp_out.write(f'{entries[rsidx].to_vcf(pop, cons, is_maj)}\n')

    def get_by_chr(self, chrom):
        """
        Return a RSCollection object containing all entries from given chrom.

        Parameters:
        chrom: Chromosome to get
        """

        if type(chrom) == int:
            chrom = str(chrom)

        # Loop through all entries and add them to the new RSCollection if they
        #  have the right chromosome
        rsc = RSCollection()
        for (_, c, _), e in self.entries.items():
            if c == chrom:
                rsc.add_entry(e)

        return(rsc)

    def get_by_chr_pos(self, chrom, pos):
        """
        Return a list of all entries with given chromosome and position.

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
            return([])

        return([self.entries[r] for r in rsidx])

    def get_by_rsid(self, rsid):
        """
        Return a list of all entries with given RefSNP ID number.

        Parameters:
        rsid: RefSNP ID number to get
        """

        if type(rsid) == int:
            rsid = str(rsid)

        try:
            rsidx = self.rsid_table[rsid]
        except KeyError:
            return([])

        return([self.entries[r] for r in rsidx])


    def get_major(self, mut=True):
        """
        Return a RSCollection object containing all major alleles.

        Parameters:
        mut: Only include variants that are different between reference and
            alternate allele
        """

        major_col = RSCollection()
        for (rsid, chrom, pos), rse in self.entries.items():
            rse_maj = RSEntry(chrom, rsid, pos)
            for v in rse.get_major_alleles():
                if v.var_type != '' or not mut:
                    rse_maj.add_var(v)

            if len(rse_maj) > 0:
                major_col.add_entry(rse_maj)
                
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
    def from_1000gp(fn, index_fn, superpop_fn, quiet=False):
        """
        Create a new RSCollection object from a 1000 Genomes Project VCF file.

        Parameters:
        fn: Filename to load from
        index_fn: Filename of the 1000 Genomes .index file to use
        superpop_fn: Filename containing superpopulation information
        quiet: Disable log/progress messages
        """

        # Construct sample -> pop dict
        d = np.loadtxt(index_fn, usecols=(9,10), dtype=str, delimiter='\t')
        pop_dict = {p[0]: p[1] for p in d}

        # Construct pop -> super-pop dict
        d = np.loadtxt(superpop_fn, skiprows=1, delimiter='\t', dtype=str,
            usecols=(1,2))
        superpop_dict = {sp[0]: sp[1] for sp in d \
            if sp[0] != '' and sp[1] != ''}

        fn_base = os.path.basename(fn)

        chrom = RSCollection.get_chrom_from_filename(fn)
        rscol = RSCollection()
        add_entry = rscol.add_entry

        samp_dict = {}
        for n, line in enumerate(RSCollection.open(fn)):
            if n % 10000 == 0 and not quiet:
                print(f'{fn_base}: Processing line {n}...',
                    flush=True)
            try:
                line = line.decode().strip()
            except AttributeError:
                line = line.strip()

            # Skip metadata lines
            if line[:2] == '##':
                continue
            elif line[0] == '#':
                # Build dictionary from sample to population information
                # Sample IDs are given as all entries in the header line after
                #  the standard VCF headers
                rec_samples = line.split('\t')[9:] 
                for s in rec_samples:
                    try:
                        p = pop_dict[s]
                    except KeyError:
                        if not quiet:
                            print(f'{s} not found in index file.', flush=True)
                        continue
                    try:
                        super_p = superpop_dict[p]
                    except KeyError:
                        if not quiet:
                            print(f'{p} not found in superpop file.',
                                flush=True)
                        continue
                    p = f'{super_p}_{p}'

                    samp_dict[s] = (p, super_p)
                continue

            line_split = line.split('\t')
            c = line_split[0]
            rec_pos = int(line_split[1])
            try:
                # Get rid of the 'rs' at the beginning of the id
                rsids = [int(rs[2:]) if rs != '.' else -1 \
                    for rs in line_split[2].split(';')]
            except ValueError:
                # Use -1 as a placeholder
                rsids = [-1]

            rec_ref = line_split[3]
            if len(rec_ref) > 50 or 'N' in rec_ref:
                continue
            # Skip SVs in alt alleles
            rec_alts = [a for a in line_split[4].split(',') if '<' not in a]
            if len(rec_alts) == 0: continue
            if len(rsids) == 1:
                rsids = rsids*len(rec_alts)
            rec_info = line_split[7]
            rec_gts = line_split[9:]
            # Find allele counts in the info
            ACs = re.search('(?:^|;)AC=([0-9,]+)(?:;|$)', rec_info)
            if ACs is None:
                if not quiet:
                    print((f'{fn_base}: Skipped ({rsids}, {c}, {rec_pos}): '
                        'no ACs'), flush=True)
                continue
            ACs = [int(a) for a in ACs.group(1).split(',')]
            if all([a == 0 for a in ACs]):
                if not quiet:
                    print((f'{fn_base}: Skipped ({rsids}, {c}, {rec_pos}): '
                        'all ACs 0'), flush=True)
                continue

            # Find total number of samples
            AN = re.search('(?:^|;)AN=([0-9]+);', rec_info)
            try:
                AN = int(AN.group(1))
            except AttributeError:
                AN = sum([s.count('|') + s.count('/') + 1 for s in rec_gts]) \
                    - sum([s.count('.') for s in rec_gts])
            if AN == 0:
                if not quiet:
                    print(f'{fn_base}: Skipped ({rsids}, {c}, {rec_pos}): AN 0',
                        flush=True)
                continue

            # Populate population allele frequencies by counting all alleles in
            #  the sample information
            pop_allele_dict = {}
            for i in range(len(rec_samples)):
                alleles = [int(a) for a in re.split('[|/]', rec_gts[i]) \
                    if a != '.']

                try:
                    p, super_p = samp_dict[rec_samples[i]]
                except KeyError:
                    continue

                if p not in pop_allele_dict:
                    pop_allele_dict[p] = {}
                if super_p not in pop_allele_dict:
                    pop_allele_dict[super_p] = {}

                for a in alleles:
                    if a == '.': continue
                    try:
                        pop_allele_dict[p][a] += 1
                    except KeyError:
                        pop_allele_dict[p][a] = 1

                    try:
                        pop_allele_dict[super_p][a] += 1
                    except KeyError:
                        pop_allele_dict[super_p][a] = 1

            for i in range(len(rec_alts)):
                # Nobody has the variant or it's too long or contains N base
                if ACs[i] == 0 or 'N' in rec_alts[i] or len(rec_alts[i]) > 50:
                    continue

                try:
                    rsid = rsids[i]
                except IndexError:
                    rsid = -1

                try:
                    rse = rscol[(str(rsid), str(c), str(rec_pos))]
                except KeyError:
                    rse = RSEntry(c, rsid, rec_pos)
                add_var_from_args = rse.add_var_from_args

                afs = [(ACs[i], AN)]
                var_type = 'SNP' if len(rec_ref) == len(rec_alts[i]) \
                    else 'indel'

                # Allele count info isn't stored for each pop in the VCF file
                # Need to use the index file and build a dict to check each
                #  population.
                pop_afs = {p: [(af_list[i+1], sum(af_list.values()))] \
                    for p, af_list in pop_allele_dict.items() if i+1 in af_list}

                add_var_from_args(rec_pos, rec_ref, rec_alts[i], afs=afs,
                    var_type=var_type, pop_afs=pop_afs)

                if (str(rsid), str(c), str(rec_pos)) not in rscol.entries:
                    add_entry(rse)

        return(rscol)

    @staticmethod
    def from_gnomad(fn, quiet=False):
        """
        Create a new RSCollection object from a gnomAD VCF file.

        Parameters:
        fn: Filename to load
        quiet: Whether to disable log messages
        """

        chrom = RSCollection.get_chrom_from_filename(fn)
        rscol = RSCollection()
        add_entry = rscol.add_entry

        # Get all of the relevant IDs for constructing sub-pops
        re_search = '(?:for|in) samples(?: |$)of (.+) ancestry">'
        id_search = 'ID=AF_([a-z_]+),'
        ids = set()
        for n, line in enumerate(RSCollection.open(fn)):
            if n % 10000 == 0 and not quiet:
                print(f'Processing line {n}...', flush=True)
            try:
                line = line.decode().strip()
            except AttributeError:
                line = line.strip()

                if line[:2] == '##':
                    if re.search(re_search, line):
                        try:
                            p = re.search(id_search, line).group(1)
                        except AttributeError:
                            continue

                        ids.add(p)
                    continue
            if type(ids) == set:
                ids = list(ids)

            line_split = line.split('\t')
            c = line_split[0].strip('chr')
            rec_pos = int(line_split[1])
            try:
                # Get rid of the 'rs' at the beginning of the id
                rsids = [int(rs[2:]) if rs != '.' else -1 \
                    for rs in line_split[2].split(';')]
            except ValueError:
                # Use -1 as a placeholder
                rsids = [-1]

            if len(rsids) == 1:
                rsids = rsids*len(rec_alts)
                
            if n % 50000 == 0 and not quiet:
                print(f'chr {c}: Processed {n} records.', flush=True)

            rec_ref = line_split[3]
            if len(rec_ref) > 50 or 'N' in rec_ref:
                continue
            rec_alts = line_split[4].split(',')
            rec_info = line_split[7]
            ACs = re.search('(?:^|;)AC=([0-9,]+);', rec_info)
            if ACs is None:
                continue
            ACs = [int(a) for a in ACs.group(1).split(',')]

            AN = re.search('(?:^|;)AN=([0-9]+);', rec_info)
            if AN is None:
                continue
            AN = int(AN.group(1))
            if AN == 0:
                continue

            # Build pop_afs dict (only need to do it once for each rec)
            rec_pop_afs = {}
            for p in ids:
                m = dict(re.findall(f';(A[CN])_{p}=([0-9]+)', rec_info))
                if len(m) == 0:
                    continue
                m['AC'] = [int(c) for c in m['AC'].split(',')]
                m['AN'] = int(m['AN'])
                rec_pop_afs[p] = m

            for i in range(len(rec_alts)):
                # Nobody has the variant or it's too long or contains N base
                if ACs[i] == 0 or 'N' in rec_alts[i] or len(rec_alts[i]) > 50:
                    continue

                rsid = rsids[i]
                try:
                    rse = rscol[(str(rsid), str(c), str(rec_pos))]
                except KeyError:
                    rse = RSEntry(c, rsid, rec_pos)
                add_var_from_args = rse.add_var_from_args
                
                afs = [(ACs[i], AN)]
                var_type = 'SNP' if len(rec_ref) == len(rec_alts[i]) \
                    else 'indel'

                pop_afs = {p: [(pop['AC'][i], pop['AN'])] \
                    for p, pop in rec_pop_afs.items() \
                    if pop['AC'][i] > 0 and pop['AN'] > 0}

                add_var_from_args(rec_pos, rec_ref, rec_alts[i], afs=afs,
                    var_type=var_type, pop_afs=pop_afs)

                if (str(rsid), str(c), str(rec_pos)) not in rscol.entries:
                    add_entry(rse)

        return(rscol)

    @staticmethod
    def from_dbsnp(fn, quiet=False):
        """
        Create a new RSCollection object from a dbSNP JSON file.
        
        Convert each line in the file into a JSON object and perform various
        checks/calculations. Each line becomes one RSEntry object.

        Parameters:
        fn: Filename to load
        quiet: Whether to disable log messages
        """

        with RSCollection.open(fn) as fin:
            counter = 0
            chrom = RSCollection.get_chrom_from_filename(fn)
            rscol = RSCollection()
            for line in fin:
                if counter % 100000 == 0 and not quiet:
                    print(f'chr {chrom}: Processed {counter} lines.',
                        flush=True)
                counter += 1
                
                try:
                    line = line.decode()
                except AttributeError:
                    pass

                rs_obj = json.loads(line)
                rsid = rs_obj['refsnp_id']
                rse = RSEntry(chrom, rsid)

                # Check for valid entries
                if not 'primary_snapshot_data' in rs_obj or \
                   len(rs_obj['primary_snapshot_data']\
                       ['placements_with_allele']) == 0 or \
                   len(rs_obj['primary_snapshot_data']\
                       ['allele_annotations']) == 0:
                    continue

                if all([len(x['frequency']) == 0 for x in
                           rs_obj['primary_snapshot_data']\
                           ['allele_annotations']]):
                    continue

                for a in rs_obj['primary_snapshot_data']\
                    ['placements_with_allele']:

                    # Check for valid entries
                    if not a['is_ptlp'] or \
                       len(a['placement_annot']\
                           ['seq_id_traits_by_assembly']) == 0 or \
                        not 'GRCh38' in a['placement_annot']\
                        ['seq_id_traits_by_assembly'][0]['assembly_name']:
                        continue

                    for i in a['alleles']:
                        spdi = i['allele']['spdi']
                        if len(spdi['deleted_sequence']) > 50 or \
                            len(spdi['inserted_sequence']) > 50 or \
                            'N' in spdi['deleted_sequence'] or \
                            'N' in spdi['inserted_sequence']:
                            continue
                        # Need to adjust position from 0-indexed to 1-indexed
                        rse.add_var_from_args(pos=int(spdi['position'])+1,
                                    ref=spdi['deleted_sequence'],
                                    var=spdi['inserted_sequence'])
                
                temp_vars = np.asarray(list(rse.vars.values()))
                blank_vars = np.asarray([False]*len(rse.vars))
                for i in range(len(rse.vars)):
                    v = temp_vars[i]

                    major = 0
                    minor = 0
                    clin = []
                    afs = []
                    var_type = ''

                    try:
                        for j in rs_obj['primary_snapshot_data']\
                            ['allele_annotations'][i]['frequency']:
                            if j['observation']['deleted_sequence'] != v.ref or\
                               j['observation']['inserted_sequence'] != v.var:
                                continue
                        
                            if j['allele_count'] >= 0.5 * j['total_count']:
                                major += 1
                            else:
                                minor += 1
                            afs += [(j['allele_count'], j['total_count'])]
                        
                        for j in rs_obj['primary_snapshot_data']\
                            ['allele_annotations'][i]['clinical']:
                            clin += j['clinical_significances']

                        if major == 0 and minor == 0:
                            blank_vars[i] = True

                        v.major = major
                        v.minor = minor

                        if v.ref != v.var:
                            if len(v.ref) == len(v.var):
                                var_type = 'SNP'
                            else:
                                var_type = 'indel'
                            
                        v.clin = clin
                        v.afs = afs
                        v.var_type = var_type
                    except IndexError as e:
                        print(rsid)
                        print(str(e))

                if rse.all_vars_empty():
                    continue

                pop_afs = RSCollection.dbSNP_pops(rsid,
                    fp='/'.join(fn.split('/')[:-1]))
                temp_vars = temp_vars[~blank_vars]
                for v in temp_vars:
                    try:
                        v.pop_afs = pop_afs[(v.ref, v.var)]
                    except KeyError:
                        continue
                rse.vars = {v.var_code(): v for v in temp_vars}

                rscol.add_entry(rse)

        return(rscol)

    @staticmethod
    def get_chrom_from_filename(fn):
        """
        Helper method to parse a filename and find which chromosome it is.

        Parameters:
        fn: Filename to use
        """

        match = '(?:chr|sites\.)([0-9XYMT]+)'

        try:
            chrom = re.search(match, fn).group(1)
        except AttributeError:
            print(f'Unable to interpret chromosome from file name ({fn}), '
                'using -1 as a placeholder.')
            chrom = '-1'
        
        return(chrom)

    @staticmethod
    def load_from_file_by_chr_pos(fn, idx_file, c, pos, chunk_idx_dict=None,
        ret_chunk=False):
        """
        Find an entry in the given file using given chromosome and position.

        If ret_chunk is True, return both the (rsid, chrom, pos) of the matching
        entry and the RSCollection chunk containing the entry. Otherwise, return
        just the matching entry.

        Parameters:
        fn: Filename to use
        idx_file: Filename of the index file (will default to {fn}.idx)
        c: Chromosome to load
        pos: Position to load
        chunk_idx_dict: Dictionary mapping (rsid,chrom,pos) combinations to the
            chunk in the ConsDB file containing that entry
        ret_chunk: Whether to return the entire chunk loaded (as opposed to only
            returning the RSEntry object)
        """

        if idx_file is None:
            idx_file = f'{fn}.idx'

        if not os.path.exists(idx_file):
            raise RuntimeError(f'Index file {idx_file} not found.')

        idx_d = np.loadtxt(idx_file, dtype=object)

        if chunk_idx_dict:
            try:
                chunk_idx = chunk_idx_dict[(c, pos)][0]
            except KeyError:
                print(f'Unable to load chr {c} pos {pos} from {fn}.')
                if ret_chunk:
                    return(None, None)
                else:
                    return(None)
        else:
            # If no index dict, manually search for the given rsidx
            idx_match = f'(?:^|;)([-0-9]+,{c},{pos})[; ]'
            idx_matches = [re.search(idx_match, l) for l in idx_d[:,0]]
            idx_matches = [m.group(1) if m is not None else None \
                for m in idx_matches]
            try:
                chunk_idx = np.nonzero(idx_matches)[0][0]
            except IndexError:
                print(f'Unable to load chr {c} pos {pos} from {fn}.')
                if ret_chunk:
                    return(None, None)
                else:
                    return(None)

            rsidx = tuple(idx_matches[chunk_idx].split(','))

        print(f'Loading chunk {chunk_idx} of chr {c}')
        try:
            chunk_info = idx_d[chunk_idx,1:].astype(int)
            idx_match = f'(?:^|;)([-0-9]+),({c}),({pos})[; ]'
            m = re.search(idx_match, idx_d[chunk_idx,0])
            rsidx = (m.group(1), m.group(2), m.group(3))
        except IndexError as e:
            print(c, pos, flush=True)
            raise e


        rsc = RSCollection()
        add_entry = rsc.add_entry
        with open(fn, 'rb') as fp:
            # Seek to the beginning of the chunk, then read the appropriate
            #  amount of bytes
            fp.seek(chunk_info[1])
            r = gzip.decompress(fp.read(chunk_info[0]))
            for e in r.decode().strip().split('\n'):
                try:
                    add_entry(RSCollection.parse_file_line(
                        e.strip().split(':')))
                except AttributeError as er:
                    print(e)
                    raise er
                except ValueError as er:
                    print(e)
                    raise er

        if ret_chunk:
            return(rsidx, rsc)
        else:
            return(rsc[rsidx])

    @staticmethod
    def load_from_file_by_rsid(fn, idx_file, rsid, chunk_idx_dict=None,
        ret_chunk=False):
        """
        Find an entry in the given file using given rsid.

        If ret_chunk is True, return both the (rsid, chrom, pos) of the matching
        entry and the RSCollection chunk containing the entry. Otherwise, return
        just the matching entry.

        Use temporary files for compatibility with existing functions that take
        filename arguments.

        Parameters:
        fn: Filename to load
        idx_file: Filename of the index file (will default to {fn}.idx)
        rsid: RefSNP ID to load
        chunk_idx_dict: Dictionary mapping (rsid,chrom,pos) combinations to the
            chunk in the ConsDB file containing that entry
        ret_chunk: Whether to return the entire chunk loaded (as opposed to only
            returning the RSEntry object)
        """

        if idx_file is None:
            idx_file = f'{fn}.idx'

        if not os.path.exists(idx_file):
            raise RuntimeError(f'Index file {idx_file} not found.')

        idx_d = np.loadtxt(idx_file, dtype=object)

        if chunk_idx_dict:
            try:
                chunk_idx = chunk_idx_dict[(rsid)][0]
            except KeyError:
                print(f'Unable to load RSID {rsid} from {fn}.')
                if ret_chunk:
                    return(None, None)
                else:
                    return(None)
        else:
            idx_match = f'(?:^|;)({rsid},[0-9XYM]{1,2},[0-9]+)[; ]'
            idx_matches = [re.search(idx_match, l) for l in idx_d[:,0]]
            idx_matches = [m.group(1) if m is not None else None \
                for m in idx_matches]
            try:
                chunk_idx = np.nonzero(idx_matches)[0][0]
            except IndexError:
                print(f'Unable to load RSID {rsid} from {fn}.')
                if ret_chunk:
                    return(None, None)
                else:
                    return(None)

            rsidx = tuple(idx_matches[chunk_idx].split(','))

        print(f'Loading chunk {chunk_idx} of chr {c}')
        try:
            chunk_info = idx_d[chunk_idx,1:].astype(int)
            idx_match = f'(?:^|;)({rsid}),([0-9XYM]{1,2}),([0-9]+)[; ]'
            m = re.search(idx_match, idx_d[chunk_idx,0])
            rsidx = (m.group(1), m.group(2), m.group(3))
        except IndexError as e:
            print(rsid, flush=True)
            raise e

        rsc = RSCollection()
        add_entry = rsc.add_entry
        with open(fn, 'rb') as fp:
            fp.seek(chunk_info[1])
            r = gzip.decompress(fp.read(chunk_info[0]))
            for e in r.decode().strip().split('\n'):
                try:
                    add_entry(RSCollection.parse_file_line(
                        e.strip().split(':')))
                except AttributeError as er:
                    print(e)
                    raise er
                except ValueError as er:
                    print(e)
                    raise er

        if ret_chunk:
            return(rsidx, rsc)
        else:
            return(rsc[rsidx])

    @staticmethod
    def load_from_file_full(fn):
        """
        Load an entire ConsDB file.

        Parameters:
        fn: Filename to load
        """

        rsc = RSCollection()
        add_entry = rsc.add_entry
        with RSCollection.open(fn) as fp:
            for line in fp:
                # Need to decode line from bytes to string if the file is
                #  compressed
                try:
                    line = line.decode()
                except AttributeError:
                    pass

                add_entry(RSCollection.parse_file_line(line.strip().split(':')))

        return(rsc)

    @staticmethod
    def load_from_file_pops(fn, pops, cons=False):
        """
        Only load variants that are major variants for pop.

        Parameters:
        fn: Filename to load
        pops: Populations to load
        cons: RSCollection object is being used for a consensus
            (will only load one allele for each position in this case)
        """

        rsc = [RSCollection() for _ in pops]
        if not os.path.isfile(fn):
            print(f'{fn} not found.')
            return([None]*len(pops))
        c = RSCollection.get_chrom_from_filename(fn)
        for n, line in enumerate(RSCollection.open(fn)):
            if n % 2500000 == 0:
                print(f'chr {c}: Reading line {n}.', flush=True)
            # Need to decode line from bytes to string if the file is
            #  compressed
            try:
                line = line.decode('utf8').strip()
            except AttributeError:
                line = line.strip()

            e = RSCollection.parse_file_line(line.split(':'))
            for i in range(len(pops)):
                p = pops[i]
                add_e = RSEntry(e.chrom, e.rsid, e.pos)
                add_vars = []
                for v in e.vars.values():
                    try:
                        pop_af = max(v.calc_pop_afs(p))
                    except KeyError:
                        continue

                    if pop_af >= 0.5:
                        add_vars.append(v)

                if len(add_vars) == 0:
                    continue

                if cons:
                    add_e.add_var(max(add_vars, key=lambda v: v.pop_af(p)))
                else:
                    for v in add_vars:
                        add_e.add_var(v)
                if len(add_e.vars) > 0:        
                    rsc[i].add_entry(add_e)

        return(rsc)


    @staticmethod
    def make_chunk_idx_dict(idx_fn, key_fields=[1,2]):
        """
        Make a dict representation of the idx file.
        Useful when indexing into a file multiple times in order to avoid
        reconstructing this same dict every time.

        Use fields [1,2] if loading by chrom, pos. Use field 0 if loading by
        RSID.

        Parameters:
        idx_fn: Filename of the ConsDB index file
        key_fields: Which fields in the (rsid,chrom,pos) combination are to be
            used as the keys in the dict
        """

        idx_d = np.loadtxt(idx_fn, dtype=object)

        # Make a list of the desired key fields for each rsid,chrom,pos
        #  combination in each line
        keys = [tuple(np.asarray(k.split(','))[key_fields]) \
            for l in idx_d[:,0] for k in l.split(';')]

        # Make a list of lists of line numbers, one for each rsid,chrom,pos on
        #  the line
        line_nums = [j for i in range(idx_d.shape[0]) \
            for j in [i]*(idx_d[i,0].count(';')+1)]
        chunk_idx_dict = {coord: [] for coord in keys}
        for i in range(len(keys)):
            k = keys[i]
            chunk_idx_dict[k] += [line_nums[i]]

        return(chunk_idx_dict)

    @staticmethod
    def merge_files(fn_list, out_fn, c=None):
        """
        Merge multiple ConsDB files into one.

        Parameters:
        fn_list: List of filenames to merge
        out_fn: Filename to use for the resulting merged file
        c: Chromosome of the input files
        """

        idx_file = f'{out_fn}.idx'
        if os.path.exists(out_fn):
            os.remove(out_fn)
        if os.path.exists(idx_file):
            os.remove(idx_file)

        try:
            d_list = [[line.strip().split(':') \
                for _, line in enumerate(RSCollection.open(fn))] \
                for fn in fn_list]
        except TypeError:
            d_list = [[line.decode().strip().split(':') \
                for _, line in enumerate(RSCollection.open(fn))] \
                for fn in fn_list]
        d_list = [sorted(d, key=RSCollection.sort_rsidx_line) for d in d_list]
        ## Make sets of rsids in each file, and then also sets of rsids that
        ##  we've seen for each files. Use a cursor and keep going down the list
        ##  and adding to an RSCollection object, and also add to each seen set.
        ##  When an rsid is in all 3 seen sets or one/some of the files don't
        ##  have the rsid, add it to the list to be written. Call dump after
        ##  some number has been reached and then start over.
        rsid_sets = [set([l[0] for l in d]) for d in d_list]
        chr_pos_sets = [set([tuple(l[0].split(',')[1:]) for l in d]) \
            for d in d_list]
        if len(d_list) > 1:
            all_shared = chr_pos_sets[0].intersection(chr_pos_sets[1])
            for s in chr_pos_sets[2:]:
                all_shared = all_shared.intersection(s)
            print(f'{len(all_shared)} shared entries.')

        seen_sets = [set() for _ in d_list]
        write_rsids = set()

        col = RSCollection()
        old_size = 0
        nlines = max([len(d) for d in d_list])
        for i in range(nlines):
            if i % 10000 == 0:
                print(f'Processed {i}/{nlines} lines.', flush=True)
            cur = [d[i] if i < len(d) else None for d in d_list]
            for j in range(len(d_list)):
                try:
                    rsid = cur[j][0].split(',')
                except TypeError:
                    continue
                chr_pos = tuple(rsid[1:])

                col.add_entry(RSCollection.parse_file_line(cur[j]))
                seen_sets[j].add(chr_pos)

                if all([chr_pos not in chr_pos_sets[k] or \
                    chr_pos in seen_sets[k] for k in range(len(d_list))]):
                    write_rsids.add(chr_pos)

            if len(write_rsids) > 50000 or i == (nlines - 1):
                old_size = col.dump(out_fn, c, rsids=write_rsids,
                    old_size=old_size, append=True)
                col = RSCollection()
                write_rsids = set()


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
        elif fn.split('.')[-1] == 'bgz':
            return(gzip.open(fn, 'rb'))

        return(open(fn, 'r'))

    @staticmethod
    def parse_file_line(d):
        """
        Create a RSEntry object from given line from rsc file.

        Parameters:
        d: List containing the (rsid,chrom,pos) combination as the first entry
            and the rest of the entry string as the second entry (can be created
            by splitting the RSEntry string by ':')
        """

        rsid, chrom, pos = d[0].split(',')
        e = RSEntry(chrom, rsid, pos)
        for v in d[1].split(';'):
            ref, var, major, minor, clin, afs, var_type, pop_afs = v.split(',')
            if len(ref) > 50 or len(var) > 50 or 'N' in ref or 'N' in var:
                continue

            afs = [(int(af.split('/')[0]),
                    int(af.split('/')[1])) for af in afs.split('.')]
            if pop_afs == '':
                pop_afs = {}
            else:
                pop_afs = {af.split('~')[0]: [(int(a.split('/')[0]),
                    int(a.split('/')[1])) for a in af.split('~')[1].split('.')] \
                    for af in pop_afs.split('|')}
            e.add_var_from_args(int(pos), ref, var, int(major), int(minor),
                clin.split('.'), afs, var_type, pop_afs)
        return(e)

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
        chrom = RSCollection.chrom_to_int(chrom)
        pos = int(pos)

        return(chrom, pos)

    @staticmethod
    def sort_rsidx_line(line_split):
        """
        Helper method to sort lines from a ConsDB file. Meant to be used as a
        key for sorting. Returns the chromosome and position of the line, both
        in int form.

        Parameters:
        line_split: List consisting of ['rsid,chrom,pos', 'rest of line']
        """

        try:
            line_split = line_split.split(':')
        except AttributeError:
            pass

        rsid, chrom, pos = line_split[0].split(',')
        chrom = RSCollection.chrom_to_int(chrom)
        pos = int(pos)

        return(chrom, pos)