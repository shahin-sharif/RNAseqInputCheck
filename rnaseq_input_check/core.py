"""Read-only RNA-seq input validation. Python standard library only."""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import html
import itertools
import json
import math
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

VERSION = '0.1.0'


def open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8') if str(path).endswith('.gz') else open(path, encoding='utf-8')


def read_json(path):
    def invalid_constant(value):
        raise ValueError(f'Non-finite JSON constant: {value}')
    return json.loads(path.read_text(), parse_constant=invalid_constant)


def resolve(value, base):
    p = Path(value).expanduser()
    return p.resolve() if p.is_absolute() else (base / p).resolve()


class Checker:
    def __init__(self, config, base):
        self.cfg, self.base = config, base
        self.checks, self.samples, self.metrics, self.artifacts = [], [], {}, {}
        self.tx, self.sequences, self.gene_sequences = {}, Counter(), defaultdict(set)
        self.inputs = set()
        self.annotation_ok = False

    def add(self, status, code, scope, message):
        self.checks.append(dict(status=status, check=code, scope=str(scope), message=str(message)))

    def record(self, path):
        self.inputs.add(path)
        if not path.is_file():
            raise ValueError(f'File does not exist: {path}')
        return path

    def stage(self, name, fn):
        try:
            fn()
        except (OSError, ValueError, KeyError, TypeError, csv.Error, EOFError, subprocess.SubprocessError) as e:
            self.add('ERROR', name + '.exception', name, str(e))

    def sample_sheet(self):
        path = self.record(resolve(self.cfg['samples'], self.base))
        with open_text(path) as f:
            reader = csv.DictReader(f, delimiter='\t')
            cols = reader.fieldnames or []
            if len(cols) != len(set(cols)):
                raise ValueError('Duplicate sample-sheet column names')
            if not {'sample_id', 'condition'} <= set(cols):
                raise ValueError('Sample sheet requires sample_id and condition columns')
            if not {'bam', 'quant'} & set(cols):
                raise ValueError('Sample sheet requires a bam or quant column (or both)')
            rows = list(reader)
        if not rows:
            raise ValueError('Empty sample sheet')
        seen = set()
        paths = defaultdict(dict)
        for line, row in enumerate(rows, 2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f'Incorrect number of fields on line {line}')
            row = {k: v.strip() for k, v in row.items()}
            sid = row['sample_id']
            if not sid or not row['condition'] or sid in seen:
                raise ValueError(f'Empty or duplicate sample ID / empty condition on line {line}')
            seen.add(sid)
            for kind in ('bam', 'quant'):
                if kind not in cols:
                    continue
                if not row[kind]:
                    raise ValueError(f'{sid}: empty {kind} path; omit the entire column for an unused input type')
                p = resolve(row[kind], path.parent)
                if kind == 'quant' and p.is_dir():
                    p = p / 'quant.sf'
                self.record(p)
                if p in paths[kind]:
                    raise ValueError(f'{sid} and {paths[kind][p]} refer to the same {kind} file: {p}')
                paths[kind][p] = sid
                row[kind] = str(p)
            self.samples.append(row)
        groups = Counter(r['condition'] for r in rows)
        self.metrics['samples_per_condition'] = dict(groups)
        self.add('PASS', 'samples.structure', 'samples', f'{len(rows)} unique sample IDs; paths exist and are distinct within input type')
        if len(groups) < 2:
            self.add('WARN', 'samples.groups', 'samples', 'Only one condition: input QC is possible, differential comparisons are not')
        for group, n in groups.items():
            if n < 2:
                self.add('WARN', 'samples.replication', group, 'Only one sample; biological replication cannot be established')
        self.add('INFO', 'samples.independence', 'samples', 'Distinct paths and labels do not establish independent biological replicates. Review clone, culture, donor and batch relationships.')

    def annotation(self):
        path = self.record(resolve(self.cfg['gtf'], self.base))
        genes = defaultdict(set)
        symbols = defaultdict(set)
        rows, exons = 0, 0
        bad = []
        tx = {}
        seqmax = Counter()
        for line_no, line in enumerate(open_text(path), 1):
            if not line.strip() or line.startswith('#'):
                continue
            rows += 1
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                if len(bad) < 10: bad.append(f'line {line_no}: expected nine fields')
                continue
            seq, _, feature, start, end, _, strand, _, attrs = fields
            try:
                start, end = int(start), int(end)
                if start < 1 or end < start or not seq:
                    raise ValueError()
            except ValueError:
                if len(bad) < 10: bad.append(f'line {line_no}: invalid coordinates')
                continue
            if feature != 'exon':
                continue
            exons += 1
            pairs = re.findall(r'(\S+)\s+"([^"]*)"\s*;', attrs)
            a = dict(pairs)
            critical = [key for key, _ in pairs if key in ('gene_id', 'transcript_id')]
            if len(critical) != len(set(critical)) or not a.get('gene_id') or not a.get('transcript_id') or strand not in ('+', '-'):
                if len(bad) < 10: bad.append(f'line {line_no}: exon requires unique quoted gene_id, transcript_id and +/- strand')
                continue
            tid, gid = a['transcript_id'], a['gene_id']
            value = (gid, seq, strand)
            if tid in tx and tx[tid] != value:
                if len(bad) < 10: bad.append(f'line {line_no}: transcript {tid} maps to multiple genes, sequences or strands')
            tx[tid] = value
            seqmax[seq] = max(seqmax[seq], end)
            self.sequences[seq] += 1
            genes[gid].add(seq)
            if a.get('gene_name'):
                symbols[a['gene_name']].add(gid)
        self.tx, self.seqmax, self.gene_sequences = tx, seqmax, genes
        self.metrics['annotation'] = dict(rows=rows, exons=exons, transcripts=len(tx), genes=len(genes), sequences=len(seqmax))
        if not tx:
            bad.append('No usable exon transcript IDs found')
        multi = {g: sorted(s) for g, s in genes.items() if len(s) > 1}
        if multi:
            bad.append(f'{len(multi)} stable gene IDs occur on multiple sequences; inspect annotation before gene grouping')
            self.artifacts['genes_multiple_sequences.json'] = multi
        collision = {s: sorted(ids) for s, ids in symbols.items() if len(ids) > 1}
        if collision:
            self.artifacts['shared_gene_symbols.json'] = collision
            self.add('WARN', 'gtf.shared_symbols', 'gtf', f'{len(collision)} gene symbols refer to multiple stable gene IDs. Group by gene_id, never by gene_name; see shared_gene_symbols.json.')
        if bad:
            raise ValueError('; '.join(bad))
        self.annotation_ok = True
        self.add('PASS', 'gtf.structure', 'gtf', f'{len(tx):,} transcripts with exons; stable gene IDs and transcript IDs preserved including versions')

    def salmon(self):
        if not self.samples or 'quant' not in self.samples[0]:
            self.add('SKIP', 'salmon.not_requested', 'salmon', 'No quant column')
            return
        baseline, lengths, hashes = None, None, {}
        for sample in self.samples:
            sid = sample['sample_id']
            def one():
                nonlocal baseline, lengths
                path = Path(sample['quant'])
                ids, lens, total_tpm, total_count, expressed = [], [], 0., 0., 0
                seen = set()
                with open_text(path) as f:
                    r = csv.DictReader(f, delimiter='\t')
                    required = {'Name', 'Length', 'EffectiveLength', 'TPM', 'NumReads'}
                    if not required <= set(r.fieldnames or []):
                        raise ValueError('quant.sf requires Name, Length, EffectiveLength, TPM, NumReads')
                    if len(r.fieldnames) != len(set(r.fieldnames)):
                        raise ValueError('Duplicate quant.sf column names')
                    for number, row in enumerate(r, 2):
                        tid = row['Name']
                        if self.cfg.get('strip_pipe', False):
                            tid = tid.split('|', 1)[0]
                        if not tid or tid in seen:
                            raise ValueError(f'line {number}: empty or duplicate transcript ID after configured normalization: {tid}')
                        seen.add(tid)
                        vals = [float(row[k]) for k in ('Length', 'EffectiveLength', 'TPM', 'NumReads')]
                        if any(not math.isfinite(v) or v < 0 for v in vals) or vals[0] <= 0:
                            raise ValueError(f'line {number}: lengths, TPM and counts must be finite and nonnegative; Length must be positive')
                        if vals[1] == 0 and (vals[2] > 0 or vals[3] > 0):
                            raise ValueError(f'line {number}: expressed transcript has zero effective length')
                        ids.append(tid); lens.append(vals[0])
                        total_tpm += vals[2]; total_count += vals[3]; expressed += vals[2] > 0
                if not math.isfinite(total_count) or not math.isfinite(total_tpm):
                    raise ValueError('Quantification totals overflow finite numeric range')
                if not ids or total_count <= 0:
                    raise ValueError('Empty quantification or no assigned fragments')
                self.add('PASS', 'salmon.structure', sid, f'{len(ids):,} unique transcript IDs; finite numeric values')
                self.add('PASS' if abs(total_tpm - 1e6) <= 100 else 'WARN', 'salmon.tpm_sum', sid, f'TPM sum {total_tpm:.6f}; expected approximately 1,000,000 (tolerance 100)')
                if baseline is None:
                    baseline, lengths = ids, dict(zip(ids, lens))
                else:
                    same_set = set(ids) == set(baseline)
                    self.add('PASS' if same_set else 'ERROR', 'salmon.id_set', sid, 'Same transcript ID set as first sample' if same_set else 'Transcript ID set differs from first sample')
                    if same_set:
                        self.add('PASS' if ids == baseline else 'WARN', 'salmon.id_order', sid, 'Transcript order matches' if ids == baseline else 'Order differs: downstream code must join by ID, not row position')
                        changed = sum(lengths[t] != length for t, length in zip(ids, lens))
                        self.add('ERROR' if changed else 'PASS', 'salmon.length_consistency', sid, f'{changed} transcript lengths differ from first sample')
                if self.annotation_ok:
                    missing = sorted(set(ids) - self.tx.keys())
                    self.add('ERROR' if missing else 'PASS', 'salmon.gtf_coverage', sid, f'{len(missing):,} quantified transcript IDs lack GTF exons; versions are matched exactly')
                    if missing:
                        self.artifacts.setdefault('missing_transcripts.json', {})[sid] = missing
                else:
                    self.add('SKIP', 'salmon.gtf_coverage', sid, 'Annotation did not pass structural validation')
                meta_path = path.parent / 'aux_info' / 'meta_info.json'
                metric = dict(transcripts=len(ids), expressed_transcripts=expressed, tpm_sum=total_tpm, estimated_fragments=total_count)
                if meta_path.is_file():
                    self.record(meta_path)
                    meta = read_json(meta_path)
                    if not isinstance(meta, dict):
                        raise ValueError('Salmon metadata must be a JSON object')
                    types = meta.get('library_types')
                    metric['library_types'] = types
                    metric['mapping_percentage'] = meta.get('percent_mapped')
                    expected = self.cfg.get('expected_library_type')
                    if expected:
                        self.add('PASS' if types == [expected] else 'ERROR', 'salmon.library_type', sid, f'Expected {[expected]}, metadata reports {types}; metadata is not an independent strandedness experiment')
                    for key in ('index_seq_hash', 'index_name_hash'):
                        if meta.get(key):
                            if key in hashes and hashes[key] != meta[key]:
                                self.add('ERROR', 'salmon.index_hash', sid, f'{key} differs between samples')
                            else:
                                hashes[key] = meta[key]
                else:
                    self.add('WARN', 'salmon.metadata', sid, 'Missing aux_info/meta_info.json; library type and index metadata cannot be checked')
                self.metrics.setdefault('salmon', {})[sid] = metric
            self.stage('salmon.' + sid, one)

    def bam(self):
        if not self.samples or 'bam' not in self.samples[0]:
            self.add('SKIP', 'bam.not_requested', 'bam', 'No bam column')
            return
        exe = shutil.which(self.cfg.get('samtools', 'samtools'))
        if not exe:
            raise ValueError('BAM inputs requested but samtools is not available on PATH; install it or set samtools to its executable path')
        timeout = self.cfg.get('command_timeout', 120)
        def run(*args):
            p = subprocess.run([exe, *args], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)
            if p.returncode:
                raise ValueError(f'samtools {args[0]} failed: {p.stderr.strip()[:2000]}')
            return p.stdout
        self.metrics['samtools'] = run('--version-only').strip()
        dictionary = None
        fai = None
        if self.cfg.get('fai'):
            fp = self.record(resolve(self.cfg['fai'], self.base))
            fai = {}
            for line in fp.read_text().splitlines():
                a = line.split('\t')
                if len(a) < 2 or a[0] in fai or int(a[1]) <= 0:
                    raise ValueError('Invalid or duplicate FASTA index entries')
                fai[a[0]] = int(a[1])
        for s in self.samples:
            sid, p = s['sample_id'], Path(s['bam'])
            def one():
                nonlocal dictionary
                run('quickcheck', '-v', str(p))
                self.add('PASS', 'bam.quickcheck', sid, 'Header and end-of-file check passed; this does not decode all alignments')
                hdr = run('view', '-H', str(p))
                sq, order = {}, None
                for line in hdr.splitlines():
                    a = line.split('\t')
                    if a[0] not in ('@SQ', '@HD'): continue
                    d = dict(x.split(':', 1) for x in a[1:] if ':' in x)
                    if a[0] == '@HD': order = d.get('SO')
                    if a[0] == '@SQ':
                        if not d.get('SN') or d['SN'] in sq or int(d.get('LN', '0')) <= 0:
                            raise ValueError('Invalid or duplicate BAM sequence dictionary entry')
                        sq[d['SN']] = int(d['LN'])
                if not sq: raise ValueError('No reference dictionary in BAM header')
                self.add('PASS' if order == 'coordinate' else 'ERROR', 'bam.sort_order', sid, f'Declared sort order: {order}; coordinate is required. Header alone does not prove sorting.')
                if dictionary is None: dictionary = sq
                else:
                    self.add('PASS' if list(sq.items()) == list(dictionary.items()) else 'ERROR', 'bam.dictionary', sid, 'Reference dictionary comparison with first BAM (including order)')
                indexes = [Path(str(p) + '.bai'), p.with_suffix('.bai'), Path(str(p) + '.csi'), p.with_suffix('.csi')]
                found = [i for i in indexes if i.is_file()]
                if not found:
                    self.add('ERROR', 'bam.index', sid, 'No .bai or .csi index found; coordinate-sort and index before analysis')
                else:
                    for i in found: self.record(i)
                    if any(i.stat().st_mtime < p.stat().st_mtime for i in found):
                        self.add('WARN', 'bam.index_age', sid, 'An index is older than the BAM; regenerate it if the BAM changed')
                    run('view', '-c', str(p), next(iter(sq)) + ':1-1')
                    self.add('PASS', 'bam.index_query', sid, 'An indexed region query succeeded; not proof that the index belongs to this exact BAM')
                metric = dict(sequences=len(sq), sort_order=order)
                if self.annotation_ok:
                    shared = set(sq) & self.seqmax.keys()
                    missing = set(self.seqmax) - set(sq)
                    overruns = [k for k in shared if self.seqmax[k] > sq[k]]
                    primary = lambda n: bool(re.fullmatch(r'(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)', n))
                    missing_primary = sorted(n for n in missing if primary(n))
                    self.add('ERROR' if not shared or missing_primary else 'PASS', 'bam.annotation_names', sid,
                             f'{len(shared)} annotated sequences present; {len(missing)} absent; missing human-style primary names: {missing_primary}. Exact names are required; aliases are not silently applied.')
                    if missing:
                        self.artifacts.setdefault('gtf_sequences_absent_from_bam.json', {})[sid] = sorted(missing)
                        self.add('WARN', 'bam.annotation_scope', sid, 'Some annotation sequences are absent from BAM. Their features cannot be quantified from this alignment reference; inspect the sequence audit, even if they are alternate loci.')
                    self.add('ERROR' if overruns else 'PASS', 'bam.annotation_coordinates', sid, f'{len(overruns)} shared sequences have exon ends beyond BAM reference lengths')
                    metric.update(shared_gtf_sequences=len(shared), absent_gtf_sequences=len(missing), coordinate_overruns=overruns)
                else:
                    self.add('SKIP', 'bam.annotation_names', sid, 'Annotation did not pass structural validation')
                if fai is not None:
                    absent = set(sq) - set(fai)
                    mismatch = [k for k in set(sq) & set(fai) if sq[k] != fai[k]]
                    self.add('ERROR' if absent or mismatch else 'PASS', 'bam.fai', sid, f'FASTA index: {len(absent)} BAM sequence names absent, {len(mismatch)} different lengths. Matching lengths do not establish sequence identity.')
                self.metrics.setdefault('bam', {})[sid] = metric
            self.stage('bam.' + sid, one)
        self.add('INFO', 'bam.scope', 'bam', 'No full BAM decoding, read-pair integrity, NH-tag, empirical strandedness or sequence identity test is performed by version 0.1.0.')

    def design(self):
        if not self.samples:
            raise ValueError('Sample sheet must pass first')
        formula = self.cfg.get('design', '~ condition')
        if not isinstance(formula, str) or not formula.startswith('~'):
            raise ValueError('Design must begin with ~')
        # Deliberately restricted grammar: additive, interaction and factorial terms.
        rhs = formula[1:].strip()
        if not re.fullmatch(r'[A-Za-z0-9_\s+:*]+', rhs):
            raise ValueError('Supported design syntax: ~ factor + numeric + factor:numeric or factor*numeric; intercept included. R functions, parentheses and subtraction are not supported.')
        terms = []
        for term in rhs.split('+'):
            term = term.strip()
            if term == '1': continue
            if '*' in term:
                if ':' in term: raise ValueError('Do not mix * and : in one term')
                parts = [p.strip() for p in term.split('*')]
                if len(parts) > 4: raise ValueError('At most four variables per factorial term')
                terms.extend(':'.join(c) for n in range(1, len(parts) + 1) for c in itertools.combinations(parts, n))
            else: terms.append(term)
        terms = list(dict.fromkeys(terms))
        variables = list(dict.fromkeys(p.strip() for term in terms for p in term.split(':')))
        numeric = self.cfg.get('numeric_covariates', [])
        refs = self.cfg.get('reference_levels', {})
        if set(numeric) - set(variables) or set(refs) - set(variables):
            raise ValueError('numeric_covariates/reference_levels names must occur in the design')
        basis = {}
        for var in variables:
            if not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_]*', var) or any(var not in s or not s[var] for s in self.samples):
                raise ValueError(f'Missing or invalid design variable: {var}')
            values = [s[var] for s in self.samples]
            if var in numeric:
                vals = list(map(float, values))
                if not all(math.isfinite(v) for v in vals): raise ValueError(f'Nonfinite covariate {var}')
                basis[var] = [(var, vals)]
            else:
                levels = sorted(set(values))
                ref = refs.get(var, levels[0])
                if ref not in levels: raise ValueError(f'Reference {ref} not found in {var}')
                if len(levels) < 2: raise ValueError(f'{var} has only one level; remove this constant factor from the design')
                basis[var] = [(var + '[' + level + ']', [float(v == level) for v in values]) for level in levels if level != ref]
        columns = [('Intercept', [1.] * len(self.samples))]
        for term in terms:
            for combo in itertools.product(*(basis[v.strip()] for v in term.split(':'))):
                columns.append((':'.join(c[0] for c in combo), [math.prod(c[1][i] for c in combo) for i in range(len(self.samples))]))
        # Scale columns before elimination to reduce unit-related rank errors.
        a = [[col[1][i] / max(max(abs(v) for v in col[1]), 1e-300) for col in columns] for i in range(len(self.samples))]
        rank = 0
        for j in range(len(columns)):
            if rank == len(a): break
            pivot = max(range(rank, len(a)), key=lambda i: abs(a[i][j]))
            if abs(a[pivot][j]) < 1e-10: continue
            a[rank], a[pivot] = a[pivot], a[rank]
            value = a[rank][j]
            a[rank] = [v / value for v in a[rank]]
            for i in range(rank + 1, len(a)):
                fac = a[i][j]
                a[i] = [x - fac * y for x, y in zip(a[i], a[rank])]
            rank += 1
        self.metrics['design'] = dict(formula=formula, columns=[c[0] for c in columns], rank=rank, samples=len(a), residual_df=len(a) - rank)
        self.artifacts['design_matrix.tsv'] = [['sample_id'] + [c[0] for c in columns]] + [[s['sample_id']] + [c[1][i] for c in columns] for i, s in enumerate(self.samples)]
        self.add('PASS' if rank == len(columns) and len(a) > rank else 'ERROR', 'design.rank', 'design', f'{len(a)} samples, {len(columns)} coefficients, rank {rank}, residual degrees of freedom {len(a) - rank}. Rank deficiency indicates confounding; zero residual degrees of freedom cannot estimate within-group variability.')
        self.add('INFO', 'design.scope', 'design', 'This is a restricted formula/rank check, not a fitted statistical model or confirmation of biological independence. Coefficient names here use factor[level] notation.')

    def r_packages(self):
        packages = self.cfg.get('r_packages', [])
        if not packages:
            self.add('SKIP', 'r.dependencies', 'R', 'Optional R dependency probe not requested')
            return
        if any(not isinstance(p, str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9.]*', p) for p in packages):
            raise ValueError('Invalid R package name')
        exe = shutil.which(self.cfg.get('rscript', 'Rscript'))
        if not exe: raise ValueError('Rscript executable not found')
        # Arguments are passed as argv, never evaluated as package names or shell commands.
        code = '''for (p in commandArgs(TRUE)) { ok <- requireNamespace(p, quietly=TRUE); cat(p, if(ok) as.character(packageVersion(p)) else "MISSING", sep="\\t"); cat("\\n") }; if ("IsoformSwitchAnalyzeR" %in% commandArgs(TRUE) && requireNamespace("IsoformSwitchAnalyzeR",quietly=TRUE)) { cat("API_fixStringTieAnnotationProblem\\t", "fixStringTieAnnotationProblem" %in% names(formals(IsoformSwitchAnalyzeR::importRdata)), "\\n",sep="") }'''
        p = subprocess.run([exe, '--vanilla', '-', *packages], input=code + '\n', capture_output=True, text=True, timeout=self.cfg.get('command_timeout', 120))
        if p.returncode: raise ValueError(p.stderr[-2000:])
        rows = dict(line.split('\t', 1) for line in p.stdout.splitlines() if '\t' in line)
        self.metrics['r_packages'] = rows
        for name in packages:
            value = rows.get(name, 'MISSING')
            self.add('ERROR' if value == 'MISSING' else 'PASS', 'r.package', name, value)
        if 'IsoformSwitchAnalyzeR' in packages:
            self.add('PASS' if rows.get('API_fixStringTieAnnotationProblem') == 'TRUE' else 'ERROR', 'r.isoform_api', 'R', 'importRdata must support fixStringTieAnnotationProblem; with a curated GENCODE GTF downstream pipelines should explicitly set it FALSE to preserve gene_id. This probe does not run an import.')
        self.add('INFO', 'r.scope', 'R', 'Package loading and one isoform import argument are checked; end-to-end package compatibility is not established.')

    def run(self):
        for name, fn in [('samples', self.sample_sheet), ('gtf', self.annotation), ('salmon', self.salmon), ('bam', self.bam), ('design', self.design), ('r', self.r_packages)]:
            print(f'Checking {name}...', file=sys.stderr, flush=True)
            self.stage(name, fn)
        return self

    def write(self, out):
        summary = dict(Counter(c['status'] for c in self.checks))
        result = dict(tool='RNAseqInputCheck', version=VERSION, created_utc=datetime.now(timezone.utc).isoformat(), summary=summary, checks=self.checks, metrics=self.metrics,
                      limitations='Structural preflight only. A passing report does not establish biological validity, reference sequence identity, statistical correctness or successful execution of downstream pipelines.')
        (out / 'report.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
        (out / 'resolved_config.json').write_text(json.dumps(self.cfg, indent=2) + '\n')
        with (out / 'checks.tsv').open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=['status', 'check', 'scope', 'message'], delimiter='\t'); w.writeheader(); w.writerows(self.checks)
        with (out / 'input_manifest.tsv').open('w', newline='') as f:
            w = csv.writer(f, delimiter='\t'); w.writerow(['path', 'bytes', 'mtime_ns', 'sha256'])
            for p in sorted(self.inputs):
                if not p.is_file(): continue
                stat = p.stat(); digest = ''
                if self.cfg.get('hash_inputs', False):
                    h = hashlib.sha256()
                    with p.open('rb') as stream:
                        for block in iter(lambda: stream.read(1024 * 1024), b''): h.update(block)
                    digest = h.hexdigest()
                w.writerow([str(p), stat.st_size, stat.st_mtime_ns, digest])
        for name, obj in self.artifacts.items():
            if name.endswith('.json'): (out / name).write_text(json.dumps(obj, indent=2) + '\n')
            else:
                with (out / name).open('w', newline='') as f: csv.writer(f, delimiter='\t').writerows(obj)
        escape = lambda x: html.escape(str(x))
        rows = ''.join('<tr class="' + c['status'] + '">' + ''.join('<td>' + escape(c[k]) + '</td>' for k in ('status', 'check', 'scope', 'message')) + '</tr>' for c in self.checks)
        state = 'ACTION REQUIRED' if summary.get('ERROR') else ('REVIEW WARNINGS' if summary.get('WARN') else 'REQUESTED CHECKS PASSED')
        page = f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>RNAseqInputCheck report</title><style>body{{font:16px system-ui;max-width:1250px;margin:40px auto;padding:0 20px;color:#182332}}h1{{color:#126b73}}table{{border-collapse:collapse;width:100%}}td,th{{padding:10px;text-align:left;border-bottom:1px solid #ddd;vertical-align:top;overflow-wrap:anywhere}}.ERROR td:first-child{{color:#b31515;font-weight:bold}}.WARN td:first-child{{color:#8a5000;font-weight:bold}}.PASS td:first-child{{color:#126339}}pre{{white-space:pre-wrap;background:#eef4f6;padding:18px}}</style><h1>RNAseqInputCheck</h1><h2>{state}</h2><p>Version {VERSION} · {escape(result['created_utc'])}</p><p>{escape(summary)}</p><p>{escape(result['limitations'])}</p><p>ERROR: resolve before analysis. WARN: review and document a decision. SKIP: check not performed. INFO: interpretation or limitations. Inputs were not modified. Report paths and sample labels may be sensitive; review before sharing.</p><table><thead><tr><th>Status</th><th>Check</th><th>Scope</th><th>Finding</th></tr></thead><tbody>{rows}</tbody></table><h2>Metrics</h2><pre>{escape(json.dumps(self.metrics, indent=2))}</pre></html>'''
        (out / 'report.html').write_text(page)
        return summary


ALLOWED = {'samples', 'gtf', 'fai', 'design', 'numeric_covariates', 'reference_levels', 'expected_library_type', 'strip_pipe', 'samtools', 'rscript', 'r_packages', 'command_timeout', 'hash_inputs'}


def validate_config(cfg):
    if not isinstance(cfg, dict): raise ValueError('Configuration must be a JSON object')
    unknown = set(cfg) - ALLOWED
    if unknown: raise ValueError(f'Unknown configuration keys: {sorted(unknown)}')
    for key in ('samples', 'gtf'):
        if not isinstance(cfg.get(key), str) or not cfg[key]: raise ValueError(f'Required nonempty string: {key}')
    for key in ('fai', 'design', 'expected_library_type', 'samtools', 'rscript'):
        if key in cfg and (not isinstance(cfg[key], str) or not cfg[key]): raise ValueError(f'{key} must be a nonempty string')
    for key in ('strip_pipe', 'hash_inputs'):
        if key in cfg and not isinstance(cfg[key], bool): raise ValueError(f'{key} must be true or false')
    for key in ('numeric_covariates', 'r_packages'):
        if key in cfg and (not isinstance(cfg[key], list) or not all(isinstance(x, str) for x in cfg[key])): raise ValueError(f'{key} must be a list of strings')
    if 'reference_levels' in cfg and (not isinstance(cfg['reference_levels'], dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in cfg['reference_levels'].items())): raise ValueError('reference_levels must map variable names to level names')
    if 'command_timeout' in cfg and (type(cfg['command_timeout']) not in (int, float) or not math.isfinite(cfg['command_timeout']) or cfg['command_timeout'] <= 0): raise ValueError('command_timeout must be a positive finite number of seconds')
    if 'expected_library_type' in cfg and cfg['expected_library_type'] not in ('ISR', 'ISF', 'IU', 'OSR', 'OSF', 'OU', 'MSR', 'MSF', 'MU', 'SR', 'SF', 'U'):
        raise ValueError('Invalid expected_library_type; provide an explicit Salmon library type, not A')


def main(argv=None):
    p = argparse.ArgumentParser(description='Read-only RNA-seq input preflight: sample sheet, GTF, Salmon, BAM, design and optional R dependencies.')
    p.add_argument('--version', action='version', version=VERSION)
    p.add_argument('--config', type=Path, help='JSON configuration (paths relative to this file)')
    p.add_argument('--samples', help='TSV sample sheet (paths inside it relative to its directory)')
    p.add_argument('--gtf', help='GTF annotation, optionally gzip compressed')
    p.add_argument('--out', required=True, type=Path, help='New output directory; existing directories are never overwritten')
    p.add_argument('--fail-on-warning', action='store_true', help='Exit 1 for warnings as well as errors')
    args = p.parse_args(argv)
    try:
        if args.config:
            cfg = read_json(args.config); base = args.config.resolve().parent
            if args.samples or args.gtf: raise ValueError('Use --config OR --samples with --gtf, not both')
        else:
            cfg = dict(samples=args.samples, gtf=args.gtf); base = Path.cwd()
        validate_config(cfg)
        for key in ('samples', 'gtf', 'fai'):
            if key in cfg: cfg[key] = str(resolve(cfg[key], base))
        args.out.mkdir(parents=True, exist_ok=False)
    except (ValueError, OSError) as e:
        print(f'Configuration/output error: {e}', file=sys.stderr); return 2
    checker = Checker(cfg, base).run()
    try:
        summary = checker.write(args.out)
    except (OSError, ValueError) as e:
        print(f'Cannot finish report: {e}', file=sys.stderr); return 2
    print(f'Report: {(args.out / "report.html").resolve()}')
    print(json.dumps(summary, sort_keys=True))
    return 1 if summary.get('ERROR') or (args.fail_on_warning and summary.get('WARN')) else 0
