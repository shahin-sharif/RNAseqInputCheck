import contextlib
import gzip
import io
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from rnaseq_input_check.core import Checker, main


class InputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.gtf = self.root / 'genes.gtf'
        self.gtf.write_text('chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id "G1.1"; transcript_id "T1.1"; gene_name "GENE";\nchr1\ttest\texon\t200\t299\t.\t+\t.\tgene_id "G2.1"; transcript_id "T2.1"; gene_name "GENE";\n')
        for name in ['wt1', 'wt2', 'ko1', 'ko2']:
            d = self.root / name
            (d / 'aux_info').mkdir(parents=True)
            (d / 'quant.sf').write_text('Name\tLength\tEffectiveLength\tTPM\tNumReads\nT1.1\t100\t50\t600000\t60\nT2.1\t100\t50\t400000\t40\n')
            (d / 'aux_info' / 'meta_info.json').write_text(json.dumps(dict(library_types=['ISR'], percent_mapped=90, index_seq_hash='same')))
        self.sheet = self.root / 'samples.tsv'
        self.sheet.write_text('sample_id\tcondition\tquant\tbatch\nwt1\tWT\twt1\tA\nwt2\tWT\twt2\tA\nko1\tKO\tko1\tB\nko2\tKO\tko2\tB\n')
        self.cfg = dict(samples='samples.tsv', gtf='genes.gtf', reference_levels={'condition': 'WT'}, expected_library_type='ISR')

    def tearDown(self):
        self.temp.cleanup()

    def check(self):
        with contextlib.redirect_stderr(io.StringIO()):
            return Checker(self.cfg, self.root).run()

    def errors(self, checker):
        return [c for c in checker.checks if c['status'] == 'ERROR']

    def test_valid_and_shared_symbols_remain_separate(self):
        c = self.check()
        self.assertEqual(self.errors(c), [])
        self.assertEqual(c.tx['T1.1'][0], 'G1.1')
        self.assertEqual(c.artifacts['shared_gene_symbols.json']['GENE'], ['G1.1', 'G2.1'])

    def test_exact_versions(self):
        f = self.root / 'wt1' / 'quant.sf'
        f.write_text(f.read_text().replace('T1.1', 'T1.2'))
        c = self.check()
        self.assertTrue(any(x['check'] == 'salmon.gtf_coverage' for x in self.errors(c)))

    def test_duplicate_paths(self):
        self.sheet.write_text(self.sheet.read_text().replace('wt2\tA', 'wt1\tA'))
        self.assertTrue(any('same quant' in c['message'] for c in self.errors(self.check())))

    def test_confounding(self):
        self.cfg['design'] = '~ condition + batch'
        self.assertTrue(any(c['check'] == 'design.rank' for c in self.errors(self.check())))

    def test_numeric_interaction(self):
        self.sheet.write_text('sample_id\tcondition\tquant\ttime\nwt1\tWT\twt1\t0\nwt2\tWT\twt2\t1\nko1\tKO\tko1\t0\nko2\tKO\tko2\t1\n')
        self.cfg.update(design='~ condition*time', numeric_covariates=['time'])
        c = self.check()
        self.assertEqual(c.metrics['design']['rank'], 4)
        self.assertTrue(any(x['check'] == 'design.rank' for x in self.errors(c)))  # saturated

    def test_gzip_gtf(self):
        with gzip.open(self.root / 'genes.gtf.gz', 'wt') as f: f.write(self.gtf.read_text())
        self.cfg['gtf'] = 'genes.gtf.gz'
        self.assertEqual(self.errors(self.check()), [])

    def test_transcript_maps_to_multiple_genes(self):
        self.gtf.write_text(self.gtf.read_text().replace('T2.1', 'T1.1'))
        self.assertTrue(any(c['check'] == 'gtf.exception' for c in self.errors(self.check())))

    def test_nonfinite_counts(self):
        f = self.root / 'wt1' / 'quant.sf'
        f.write_text(f.read_text().replace('600000', 'nan'))
        self.assertTrue(any(c['check'] == 'salmon.wt1.exception' for c in self.errors(self.check())))

    def test_index_metadata(self):
        p = self.root / 'ko1' / 'aux_info' / 'meta_info.json'
        p.write_text(json.dumps(dict(library_types=['ISF'], index_seq_hash='different')))
        codes = {c['check'] for c in self.errors(self.check())}
        self.assertIn('salmon.library_type', codes)
        self.assertIn('salmon.index_hash', codes)

    def test_order_difference_is_warning(self):
        p = self.root / 'ko1' / 'quant.sf'
        a = p.read_text().splitlines()
        p.write_text('\n'.join([a[0], a[2], a[1]]) + '\n')
        c = self.check()
        self.assertEqual(self.errors(c), [])
        self.assertTrue(any(x['check'] == 'salmon.id_order' and x['status'] == 'WARN' for x in c.checks))

    def test_report_escapes_html(self):
        c = self.check(); c.add('INFO', 'test', '<script>', '<img src=x>')
        out = self.root / 'out'; out.mkdir(); c.write(out)
        page = (out / 'report.html').read_text()
        self.assertNotIn('<script>', page)
        self.assertIn('&lt;script&gt;', page)
        self.assertIn('&lt;img src=x&gt;', page)

    def test_cli_existing_output_and_unknown_key(self):
        path = self.root / 'config.json'; path.write_text(json.dumps(self.cfg))
        with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(['--config', str(path), '--out', str(self.root / 'out')]), 0)
            self.assertEqual(main(['--config', str(path), '--out', str(self.root / 'out')]), 2)
            self.cfg['typo'] = True; path.write_text(json.dumps(self.cfg))
            self.assertEqual(main(['--config', str(path), '--out', str(self.root / 'new')]), 2)

    @unittest.skipUnless(shutil.which('Rscript'), 'Rscript required')
    def test_real_r_package_probe(self):
        self.cfg['r_packages'] = ['stats', 'RNAseqInputCheckNonexistentPackage']
        c = self.check()
        self.assertTrue(any(x['scope'] == 'stats' and x['status'] == 'PASS' for x in c.checks))
        self.assertTrue(any(x['scope'] == 'RNAseqInputCheckNonexistentPackage' and x['status'] == 'ERROR' for x in c.checks))

    def test_invalid_metadata_object(self):
        p = self.root / 'wt1' / 'aux_info' / 'meta_info.json'
        p.write_text('[]')
        self.assertTrue(any(x['check'] == 'salmon.wt1.exception' for x in self.errors(self.check())))

    def test_fail_on_warning(self):
        path = self.root / 'config.json'; path.write_text(json.dumps(self.cfg))
        with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main(['--config', str(path), '--out', str(self.root / 'out'), '--fail-on-warning']), 1)

    @unittest.skipUnless(shutil.which('samtools'), 'real samtools required; CI installs it')
    def test_real_bam_alias_mismatch_and_index(self):
        sam = self.root / 'tiny.sam'
        # One shared alternate sequence must not hide a missing primary chromosome.
        sam.write_text('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:1\tLN:1000\n@SQ\tSN:ALT\tLN:1000\nr1\t0\t1\t10\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n')
        self.gtf.write_text(self.gtf.read_text() + 'ALT\ttest\texon\t1\t10\t.\t+\t.\tgene_id "G3"; transcript_id "T3";\n')
        for name in ['wt1', 'wt2', 'ko1', 'ko2']:
            bam = self.root / (name + '.bam')
            subprocess.run(['samtools', 'view', '-b', '-o', str(bam), str(sam)], check=True)
            subprocess.run(['samtools', 'index', str(bam)], check=True)
        lines = self.sheet.read_text().splitlines()
        self.sheet.write_text(lines[0] + '\tbam\n' + '\n'.join(x + '\t' + x.split('\t')[0] + '.bam' for x in lines[1:]) + '\n')
        c = self.check()
        self.assertTrue(any(x['check'] == 'bam.annotation_names' for x in self.errors(c)), self.errors(c))
        self.gtf.write_text(self.gtf.read_text().replace('chr1\t', '1\t'))
        self.assertEqual(self.errors(self.check()), [])
        (self.root / 'wt1.bam.bai').unlink()
        self.assertTrue(any(x['check'] == 'bam.index' for x in self.errors(self.check())))


if __name__ == '__main__':
    unittest.main()
